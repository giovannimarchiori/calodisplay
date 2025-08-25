#!/bin/bash

url=https://gmarchio.web.cern.ch/fcc/allegro/calodisplay/data
file1=ALLEGRO_o1_v03_noDCHcells.root
file2=allegro_v03_evts_100_pdg_11_MomentumMinMax_20.0_20.0_GeV_ThetaMinMax_20.0_160.0_PhiMinMax_0_360_HCal_ON_digi_reco.root

download_file() {
    local url="$1"
    local file="$2"
    curl "$url/$file" -o data/$file || { echo "Download failed"; exit 1; }
}

mkdir -p data
for file in $file1 $file2
do
    if ! test -f ./data/$file; then
	download_file $url $file
    fi
done
