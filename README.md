# HbondFinder
Modification of Mdtraj to customize the search of hydrogen bonds between ligand and protein in PDB trajectories.
## Input example
  python HbondFinder/main.py configuration.yaml
## YAML input configuration example (more examples in "tests" folder)
  path: "OUT/trajectory_*.pdb"

  ligname: "GRW"

  specific_hbonds: {92: ["O"],
                  94: ["O", "H"]}
                  
  angle: 2.094

  distance: 0.25

  pseudo: False

  outpath: "hbond_analysis"

  cpus: 8

  report_out: "report_hb_"

  report_in: "report_"

## Documentation

- **path**: Pattern to PDB files to analyze (glob formatt).
- **ligname**: Residue name of the LIGAND (it must be unique).
- **specific_hbonds**: Dictionary of {RESIDUE SEQUENCE NUMBER: [ "ATOM NAME 1", ATOM NAME 2", ...],...}. The counting of hydrogen bonds will be focused only on this atoms.
- **angle**: Angle X-H-X cutoff (in rad) to accept an H-bond.
- **distance**: H-X distance cutoff (in nm) to accept an H-bond.
- **pseudo**: True/False flag to also accept H-C hydrogen bonds (pseudo hydrogen bonds).
- **outpath**: Output folder to save modified report files.
- **cpus**: Number of processors to paralelize and speed up the H-bond counting.
- **top**: Topology file path. (Currently unavailable function)
- **report_out**: Report's output prefix. 
- **report_in**: Report's input prefix. 


