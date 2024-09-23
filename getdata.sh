#!/bin/bash

url=https://gmarchio.web.cern.ch/fcc/allegro/calodisplay/data
file=ALLEGRO_o1_v03_noDCHcells.root

download_file() {
    local url="$1"
    local file="$2"
    #wget -P data "$url/$file" || { echo "Download failed"; exit 1; }
    curl "$url/$file" -o data/$file || { echo "Download failed"; exit 1; }
}

mkdir -p data
if ! test -f ./data/$file; then
    download_file $url $file
fi
