# calodisplay

## Description
Allegro event display (with an emphasis on the calorimeter system).
Note: development of this tool started in early 2023 as a way of helping debugging the ECAL digitisation with the new segmentation classes under development.
I became aware only later of existing event display tools. It is likely that those tools (e.g. DDEve) can - after proper customisation for the non projective ECAL cell geometry for instance - provide similar functionality as well as many more features with a more robust code.

Works with ALLEGRO models o1_v03 and o2_v01 as of July 2026.

## Installation
Make a local clone of the gitlab repository.

Setup ROOT and the environment with:
```
source setup.sh
```
(if ROOT is already set up, i.e. `ROOTSYS` is set, it is used as is; otherwise, on lxplus and fcc-ironic ROOT is set up from the nightly key4hep stack).
`setup.sh` also checks whether calodisplay is compiled, up to date with the sources, and built with a ROOT version compatible with the one set up, and tells you whether to run `make -j8` or `make clean` followed by `make -j8`.

Compile:
```
make -j8
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
calodisplay -g data/ALLEGRO_o1_v03_noDCHcells.root  -e data/allegro_o1_v03_evts_10_pdg_11_MomentumMinMax_10.0_10.0_GeV_ThetaMinMax_90.0_90.0_PhiMinMax_0_360_HCal_ON_digi_reco.root -c config.json --fulldet --dohcal --doendcap
```
Execute with `calodisplay -h` to see available options:
- `-g|--geometry <geomFile.root>`: ROOT file with the detector geometry
- `-e|--events <eventFile.root>`: EDM4hep file with the events to display
- `-c|--config <configFile.json>`: configuration file (default: `config.json`)
- `-s|--skipEvents <nevents>`: number of the first event to display
- `--nog4`: use a simplified geometry rather than the G4 geometry in the geometry file
- `--fulldet`: show the full detector (including tracker, muon detector, interaction region)
- `--dohcal`: show HCAL volume, hits and cells
- `--doendcap`: show ECAL and (optionally) HCAL endcaps
- `--sw`: show sliding-window clusters
- `--notopo`: do not show topoclusters
- `--debug`: display extended information on the console during execution

`calodisplay` can be run from any directory: relative paths of the configuration, geometry and field map files, and the icons, are looked up first in the current directory and then in the installation directory (`$CALODISPLAY_DIR`, set by `setup.sh`, or otherwise the parent of the directory containing the executable).

The ECAL barrel layout is deduced from the name of the geometry file: names containing `o1_v03` or `o2_v01` correspond to the ECAL barrel with 11 layers, the others to the older version with 12 layers. A warning is printed if the number of layers in the geometry does not match.

### Navigation
The navigation bar at the bottom contains:
- a button to take a screenshot of the selected view
- the event number box: type an event number and press Enter to go to that event, or use the up/down arrows to go to the next/previous event. The range of the available events is shown next to the box. An error pop-up is shown for events out of range.
- a button to hide/show the calorimeter readout grid

### Configuration
The configuration is read from a JSON file (see `config.json`), containing a map of key-value strings. Keys not specified in the file keep their default values (see `src/DisplayConfig.cpp`). The main keys are:
- `drawXXX` (`true`/`false`): whether to draw a given type of objects, e.g. `drawGenParticles`, `drawTracks`, `drawVertexHits`, `drawMainTrackerDigis`, `drawECalBarrelHits`, `drawECalBarrelCells`, `drawHCalEndcapCells`, `drawMuonHits`, `drawTopoClusters`, `drawCaloClusters`, `drawClusterDirection`, ...
- collection names, e.g. `genParticles`, `tracks`, `vertexBarrelHits`, `mainTrackerDigis`, `ecalBarrelHits`, `ecalBarrelCells`, `hcalEndcapCells`, `muonBarrelCells`, `topoClusters`, `topoClusterCells`, `caloClusters`, `caloClusterCells`, ...
- energy thresholds in GeV: `energyThresholdParticles`, `energyThresholdHits`, `energyThresholdCells`, `energyThresholdClusters`
- `bField`: `GlobalSolenoid` for an ideal solenoid field, anything else to use the field map in `bFieldMap` (default: `data/ALLEGRO_fieldmap_2T_20260424.root`). If the field map cannot be read, the ideal solenoid field is used.
- `bkgColor3D`: background colour of the 3D view

### Input event files
The event file is a podio/EDM4hep ROOT file, e.g. produced by the ALLEGRO simulation, digitisation and reconstruction chain.
- The cellID encodings of the calorimeter cell collections are read from the metadata of the file (`<collection>__CellIDEncoding`). If they are not found, default encodings are used.
- The cells of the clusters are taken from the dedicated cluster cell collections (`topoClusterCells`, `caloClusterCells`) if present in the file. Otherwise, the links stored in the clusters to the cells of the cell collections in the file (e.g. `ECalBarrelModuleThetaMergedPositioned`) are followed.

### Detector geometry
To create the root file with the detector geometry (if that on the website is outdated):
- create the GDML file from Geant4 with ddsim (do e.g. `ddsim --numberOfEvents 0 --enableGun --compactFile $K4GEO/FCCee/ALLEGRO/compact/ALLEGRO_o1_v03/ALLEGRO_o1_v03.xml  --geometry.dumpGDML ALLEGRO_o1_v03.gdml --runType run`)
- convert the GDML to ROOT with `util/gdmltoroot.C`:
```
cd util
root
.L gdmltoroot.C+
gdmltoroot("../data/ALLEGRO_o1_v03_noDCHcells.gdml","../data/ALLEGRO_o1_v03_noDCHcells.root", "world")
```
- for allegro_o1_v03, I had to remove from the gdml all DCH cells since they use a volume (G4TwistedTube) not available in ROOT. Note that the latest version (July 2026) uses a very detailed model of the vertex detector that is making the rendering pretty slow, so consider hiding the elements/show the envelope only/..


## Code structure
The package is organised in the following subdirectories:
- `src`: contains the source code of the event display
- `util`: contains the script to convert a GDML file to ROOT
- `data`: contains the ROOT files with the detector models and the magnetic field map
- `icons`: contains the icons of the navigation bar
The main directory contains the Makefile, the setup shell script, and the json files containing configuration parameters for the event display.
Compiling the code creates:
- `build`: for intermediate object files and the ROOT dictionary source
- `dist`: for the shared library, the ROOT dictionary file and the executable.


## Troubleshooting
In case of errors, check that eve is installed
```
root-config --has-opengl
root-config --evelibs
```

Check that you can run the eve tutorials:
```
root $ROOTSYS/tutorials/eve/calorimeters.C
```

In case of certificate validation problems with davix, you can add `Davix.GSI.CACheck: n` to your .rootrc.

If failure persists, contact the author ;)


## Author
Giovanni Marchiori (APC Paris), giovanni.marchiori@cern.ch

## Acknowledgements
Many thanks to the Eve author, Matevz Tadel, for many useful suggestions.
