# %%
import MDAnalysis as mda
import numpy as np
from dataclasses import dataclass
from tqdm import tqdm
from multiprocessing import Pool, cpu_count
from MDAnalysis.lib.distances import capped_distance, calc_angles

@dataclass
class HBonds:
    """ Record each hydrogen bond global atom index (0-based)

    Args:
        d (int): donor index (0-based)
        a (int): acceptor index (0-based)
        h (int): hydrogen index (0-based)
        dist (float): distance of d-a
        angle (float): angle between h-d-a
    """
    d: int = None 
    a: int = None
    h: int = None 
    dist: float = None 
    angle: float = None 

    def __repr__(self):
        return f"d={self.d}, a={self.a}, h={self.h}, dist={self.dist:.2f}, angle={self.angle:.2f}"
    
    def __hash__(self):
        return hash((self.d, self.a, self.h))
    
    def __eq__(self, other):
        if not isinstance(other, HBonds):
            return False
        return (self.d, self.a, self.h) == (other.d, other.a, other.h)

class GMXHBonds:
    """ @brief Calculate hydrogen bonds between two selections
    """
    def __init__(self, tpr:str, xtc:str, 
                 selA: str, selB: str, 
                 rcut: float=3.5, 
                 acut: float=30.0):
        """ @brief The method is same as `gmx hbond` analysis and obtain identical results
        Args:
            tpr (str): topology file
            xtc (str): trajectory file
            selA (str): selection for donors
            selB (str): selection for acceptors
            rcut (float): hydrogen bond distance cutoff, unit is angstrom
            acut (float): hydrogen bond angle cutoff, unit is degree
        Note:
            OH, NH groups are regarded as donors, O is an acceptor always, N is an acceptor by default
        """
        self.tpr = tpr
        self.xtc = xtc
        self.selA = selA
        self.selB = selB
        self.rcut = rcut
        self.acut = acut
        self.is_donor = {}
        self.is_acceptor = {}
        self.neighH = {}
        self.s1 : mda.AtomGroup = None
        self.s2 : mda.AtomGroup = None
        self.u = None
        self.results = None
        self.times = []
        self._prepare()

    @staticmethod
    def _get_donors(sel: mda.AtomGroup, neighH: dict, 
                    donorlist: list=['O','N']):
        """ 
        Return global index for donors
        """
        donors = []
        for atom in sel.atoms:
            if atom.name[0] in donorlist and \
                atom.index in neighH and \
                len(neighH[atom.index]) > 0:
                donors.append(atom.index)
        return donors
    
    @staticmethod
    def _get_acceptors(sel: mda.AtomGroup, acceptorlist: list=['O','N']):
        """ 
        Return global index for acceptors
        """
        acceptors = []
        for atom in sel.atoms:
            if atom.name[0] in acceptorlist:
                acceptors.append(atom.index)
        return acceptors
    
    def _prepare(self):
        # load trajectory
        self.u = mda.Universe(self.tpr, self.xtc)

        # use static selection
        self.s1 = self.u.select_atoms(self.selA)
        self.s2 = self.u.select_atoms(self.selB)
        if self.s1.n_atoms == 0 or self.s2.n_atoms == 0:
            raise ValueError("No atoms found in selection")
        #print(f'{self.s1.n_atoms} atoms found in selection 1')
        #print(f'{self.s2.n_atoms} atoms found in selection 2')

        # build hydrogen atom neighbor list for each atom
        names = self.u.atoms.names
        for [a, b] in self.u.bonds.to_indices(): # global index
            if a not in self.neighH:
                self.neighH[a] = []
            if names[b][0] == 'H':
                self.neighH[a].append(b)
            if b not in self.neighH:
                self.neighH[b] = []
            if names[a][0] == 'H':
                self.neighH[b].append(a)

        # get donors and acceptors for two selections
        donors = self._get_donors(self.s1, self.neighH) + self._get_donors(self.s2, self.neighH)
        acceptors = self._get_acceptors(self.s2) + self._get_acceptors(self.s1)
        self.is_donor = {idx: True for idx in donors}
        self.is_acceptor = {idx: True for idx in acceptors}

        self.times = [ts.time for ts in self.u.trajectory]

    def _single_frame_run(self, fr: int):
        # goto fr
        self.u.trajectory[fr]
        ts = self.u.trajectory.ts

        # do run
        pairs, distances = capped_distance(self.s1.positions, self.s2.positions,
                                    max_cutoff=self.rcut, 
                                    box=ts.dimensions,
                                    method='pkdtree',
                                    return_distances=True)
        positions = self.u.atoms.positions
        hbdists = []
        # hydroge-donor-acceptor pairs (0-based index)
        h, d, a = [], [], []
        for k, [i, j] in enumerate(pairs):
            idx1, idx2 = self.s1[i].index, self.s2[j].index
            if idx1 == idx2: continue
            # idx1=donor, idx2=acceptor
            if self.is_donor.get(idx1, False) and self.is_acceptor.get(idx2, False):
                # pre add
                for hidx in self.neighH[idx1]:
                    d.append(idx1)
                    a.append(idx2)
                    h.append(hidx)
                    hbdists.append(distances[k])
            # idx1=acceptor, idx2=donor
            if self.is_acceptor.get(idx1, False) and self.is_donor.get(idx2, False):
                # pre add
                for hidx in self.neighH[idx2]:
                    d.append(idx2)
                    a.append(idx1)
                    h.append(hidx)
                    hbdists.append(distances[k])

        # calculate angles in radians
        ph, pd, pa = positions[h], positions[d], positions[a]
        angles = calc_angles(np.asarray(ph), np.asarray(pd), np.asarray(pa), 
                             box=ts.dimensions)
        cut = np.deg2rad(self.acut)
        indexs = np.where(angles < cut)[0]
        # remove duplicates
        results = list(set([HBonds(d[i], a[i], h[i], hbdists[i], np.rad2deg(angles[i])) for i in indexs]))
        return results

    def run(self, nproc = None):
        # all frame analysis
        frames = [x for x in range(self.u.trajectory.n_frames)]
        with Pool(processes=cpu_count() if nproc is None else nproc) as work:
            # show progress
            self.results = list(tqdm(work.imap(self._single_frame_run, frames), total=len(frames)))
        return self.results

    @property
    def counts(self):
        """ @brief Return times (ps) & numbers dict for each frame
        """
        return {
            'times': np.array(self.times, dtype=np.float32),
            'numbers' : np.array([len(i) for i in self.results], dtype=np.int64)
        }
    
    def details(self, fr: int):
        """ @brief Return hydrogen bond details of fr-th frame

        Example:
        >>> get the first frame (fr=0) data
        >>> frdata = hb.details(0)
        >>> for detail in frdata:
        >>>     print(detail['donor'])
        >>>     print(detail['acceptor'])
        >>>     print(detail['hydrogen'])
        >>>     print(detail['distance'])
        >>>     print(detail['angle'])
        """
        if fr >= len(self.results) or fr < 0:
            raise ValueError(f"fr= {fr} out of range, should be in [0, {len(self.results)-1}]")
        data = []
        resids, resnames, atomnames = self.u.atoms.resids, self.u.atoms.resnames, self.u.atoms.names
        for hb in self.results[fr]:
            hbdetail = {
                'donor':    (resids[hb.d], resnames[hb.d], atomnames[hb.d]),
                'acceptor': (resids[hb.a], resnames[hb.a], atomnames[hb.a]),
                'hydrogen': (resids[hb.h], resnames[hb.h], atomnames[hb.h]),
                'distance': hb.dist,
                'angle':    hb.angle
            }
            data.append(hbdetail)
        return data

    def _test(self):
        from MDAnalysis.analysis.hydrogenbonds import HydrogenBondAnalysis
        u = mda.Universe(self.tpr, self.xtc)
        hb = HydrogenBondAnalysis(u, between=[self.selA, self.selB],
                                 d_h_cutoff=0.108,
                                 d_a_cutoff=self.rcut,
                                 d_h_a_angle_cutoff=138.1)
        hb.run()
        print(hb.count_by_time())


# %%
if __name__ == '__main__':
    hb = GMXHBonds("test/1EBZ.tpr", "test/1EBZ.xtc", 
                "protein", "resname BEC")
    hb.run()

    import matplotlib.pyplot as plt
    plt.plot(hb.counts['times'], hb.counts['numbers'])
    print(hb.counts['times'])
    print(hb.counts['numbers'])

    frdata = hb.details(5)
    for detail in frdata:
        # each hbond detail
        print(detail['donor'])
        print(detail['acceptor'])
        print(detail['hydrogen'])
        print(detail['distance'])
        print(detail['angle']) 



