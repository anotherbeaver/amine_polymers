# amine_polymers
Sketching and simulating generic polymer chains with functional "amines" as sites for transient dynamic crosslinking.
## Generate data file
From the repo directory,
```
python generate_data.py   
```
## Run simulation
Depends on location of LAMMPS binary
```
../lammps/build/lmp -in in.amine_polymers 
```