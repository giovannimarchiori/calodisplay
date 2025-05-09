# calodisplay

## Description
Allegro event display (with an emphasis on the calorimeter system).
Note: development of this tool started in early 2023 as a way of helping debugging the ECAL digitisation with the new segmentation classes under development.
I became aware only later of existing event display tools. It is likely that those tools (e.g. DDEve) can - after proper customisation for the non projective ECAL cell geometry for instance - provide similar functionality as well as many more features with a more robust code.

## Installation
Make a local clone of the gitlab repository.

Setup ROOT and the environment with:
```
source setup.sh
```

Compile:
```
make
```

Download detector model file and an example event file:
```
source getdata.sh
```


## Usage
Setup ROOT and the environment with:
```
source setup.sh
```

Execute with e.g.
```
calodisplay -g data/allegro_o1_v03.root  -e data/<events.root> -c config.json --fulldet --dohcal --doendcap
```
Execute with `calodisplay -h` to see available options

To create the root file with the detector geometry (if that on the website is outdated):
- create the GDML file from Geant4 with ddsim (see `run_all_chain.sh`)
- convert the GDML to ROOT with `util/gdmltoroot.C`:
```
cd util
root
.L gdmltoroot.C+
gdmltoroot("../data/ALLEGRO_o1_v03_noDCHcells.gdml","../data/ALLEGRO_o1_v03_noDCHcells.root", "world")
```
- for allegro_o1_v03, I had to remove from the gdml all DCH cells since they use a volume (G4TwistedTube) not available in ROOT

## Code structure
The package is organised in the following subdirectories:
- `src`: contains the source code of the event display
- `util`: contains the script to convert a GDML file to ROOT
- `data`: contains the ROOT files with the detector models
The main directory contains the Makefile, the setup shell script, and the json files containing configuration parameters for the event display.
Compiling the code creates:
- `build`: for intermediate object files
- `dist`: for the shared library, the ROOT dictionary file and the executable.


## Author
Giovanni Marchiori (APC Paris), giovanni.marchiori@cern.ch