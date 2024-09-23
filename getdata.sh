#!/bin/bash

url=https://gmarchio.web.cern.ch/fcc/allegro/calodisplay/data/
file=ALLEGRO_o1_v03_noDCHcells.root

download_file() {
  local url="$1"
  wget -P data "$url" || { echo "Download failed"; exit 1; }
}

if ! test -f ./data/$file; then
  download_file $url$file
fi
