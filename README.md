# Hydrogen bonds Analysis script

This script is used to analyze hydrogen bonds. The method is same as `gmx hbond`

# Purpose
Resolving the discrepancies between [MDAnalysis](https://docs.mdanalysis.org/2.10.0/documentation_pages/analysis/hydrogenbonds.html#MDAnalysis.analysis.hydrogenbonds.hbond_analysis.HydrogenBondAnalysis.results.hbonds) and [gmx hbond](https://manual.gromacs.org/current/onlinehelp/gmx-hbond.html) results, which stem from different hydrogen-bond criteria, with those used by gmx hbond being more widely accepted.

# Usage
```python
hb = GMXHBonds("wat/1EBZ.tpr", "wat/1EBZ.xtc", "protein", "resname BEC")
hb.run()

# get times and hb numbers 
print(hb.counts)

# get hbonds details for frame = 5
frdata = hb.details(5)
for detail in frdata:
    # each hbond detail
    print(detail['donor'])
    print(detail['acceptor'])
    print(detail['hydrogen'])
    print(detail['distance'])
    print(detail['angle'])
```


# Basic steps
* Load gromacs tpr topology and xtc/trr trajectory by MDAnalysis
* Find all donors and acceptors by given two atom selections
* Calculate Donor-Acceptor distance and Hydrogen-Donor-Acceptor angle for each frame
* Obtain all hydrogen bonds within given dinstance cutoff and angle cutoff

