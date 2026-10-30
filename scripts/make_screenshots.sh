#!/bin/bash
#
# Event display for the ALLEGRO detector for FCC-ee
# author: Giovanni Marchiori (giovanni.marchiori@cern.ch)
#
# make_screenshots.sh: produce single-particle events with run.py and save
# screenshots of their event displays with calodisplay
#
# usage: make_screenshots.sh [gen|shots|all] [-n <events>] [-o <outdir>] [-r <rundir>] [--only <name>] [--width <pixels>]
#   gen   : only produce the events (needs the key4hep environment of run.py)
#   shots : only save the screenshots of the events already produced
#           (needs the environment of calodisplay, see setup.sh)
#   all   : both (default), if the two environments are compatible
#
# The events and the options of run.py and calodisplay are set in the list
# SAMPLES below. The file produced for each sample is recorded in
# <outdir>/<name>.file, which is read by the "shots" step.
#

# name | options of run.py | options of calodisplay
# (the framing of the projections is set with --focus cluster|track|track+cluster[:<scale>] or,
# view by view, with --camera <view>=<values>, printed by the Camera button of
# the display and, for each screenshot, in <outdir>/<name>.log)
SAMPLES=(
  "electron_endcap | -p e- --pmin 20 --pmax 20 --thmin 15 --thmax 30 | --fulldet --focus track+cluster"
  "photon_barrel   | -p gamma --pmin 20 --pmax 20 --thmin 70 --thmax 110 | --fulldet --focus cluster"
  "piminus_endcap  | -p pi- --pmin 20 --pmax 20 --thmin 15 --thmax 30 --hcal | --fulldet --focus track+cluster --set drawTopoClusters=true --set drawCaloClusters=false"
  "kaon0L_barrel   | -p kaon0L --pmin 20 --pmax 20 --thmin 70 --thmax 110 --hcal | --fulldet --focus cluster --set drawTopoClusters=true --set drawCaloClusters=false --set topoClusters=AugmentedCaloTopoClusters"
  "muminus_barrel  | -p mu- --pmin 20 --pmax 20 --thmin 70 --thmax 110 --hcal --muon | --fulldet --focus track+cluster"
)

# options of calodisplay for all the samples: no simulated hits (only the
# digis and the cells), and no cluster barycenters
COMMON_DISPOPTS="--set drawVertexHits=false --set drawMainTrackerHits=false --set drawSiWrapperHits=false
  --set drawECalBarrelHits=false --set drawECalEndcapHits=false --set drawHCalBarrelHits=false
  --set drawHCalEndcapHits=false --set drawMuonHits=false --set drawClusterBarycenters=false"

STEP=all
NEVENTS=5
CALODISPLAY_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
OUTDIR=$PWD/screenshots
RUNDIR=$CALODISPLAY_DIR/../nightly/run
ONLY=""
WIDTH=1920

while [ $# -gt 0 ]; do
  case $1 in
    gen|shots|all) STEP=$1 ;;
    -n) NEVENTS=$2; shift ;;
    -o) OUTDIR=$2; shift ;;
    -r) RUNDIR=$2; shift ;;
    --only) ONLY=$2; shift ;;
    --width) WIDTH=$2; shift ;;
    -h|--help) sed -n '2,20p' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
    *) echo "Unknown option $1"; exit 1 ;;
  esac
  shift
done
mkdir -p "$OUTDIR"
OUTDIR=$(cd "$OUTDIR" && pwd)

trim() { local s="$1"; s="${s#"${s%%[![:space:]]*}"}"; echo "${s%"${s##*[![:space:]]}"}"; }

# produce the events of a sample with run.py, and record the output file
generate() {
  local name=$1 runopts=$2
  # extra options
  runopts="${runopts} --savehits --savecells --conformaltracks"
  echo "=== $name: producing $NEVENTS events with run.py $runopts"
  local stamp
  stamp=$(mktemp)
  (cd "$RUNDIR" && python run.py -n "$NEVENTS" $runopts) || { echo "ERROR: run.py failed for $name"; rm -f "$stamp"; return 1; }
  # the reconstructed file produced by this job: the newest one in root/
  local file
  file=$(find "$RUNDIR/root" -maxdepth 1 -name "*_reco.root" -newer "$stamp" -printf "%T@ %p\n" | sort -n | tail -1 | cut -d' ' -f2-)
  rm -f "$stamp"
  if [ -z "$file" ]; then
    echo "ERROR: no reconstructed file found in $RUNDIR/root for $name"
    return 1
  fi
  echo "$file" > "$OUTDIR/$name.file"
  echo "=== $name: events in $file"
}

# save the screenshots of the events of a sample
screenshots() {
  local name=$1 dispopts=$2
  if [ ! -f "$OUTDIR/$name.file" ]; then
    echo "ERROR: no events for $name (run the gen step first)"
    return 1
  fi
  local file
  file=$(cat "$OUTDIR/$name.file")
  echo "=== $name: screenshots of $file"
  (cd "$CALODISPLAY_DIR" && $VGLRUN dist/calodisplay $COMMON_DISPOPTS $dispopts -e "$file" -c config.json \
     --screenshots "$OUTDIR/$name" --nevents "$NEVENTS" --width "$WIDTH" < /dev/null > "$OUTDIR/$name.log" 2>&1) \
    || { echo "ERROR: calodisplay failed for $name, see $OUTDIR/$name.log"; return 1; }
  grep "Screenshot saved" "$OUTDIR/$name.log" | sed 's/^/    /'
}

# calodisplay needs an X display: without one, start a private virtual one
XVNC_PID=""
if [ "$STEP" != gen ]; then
  if [ -z "$DISPLAY" ]; then
    if ! command -v Xvnc > /dev/null; then
      echo "ERROR: no X display (DISPLAY not set) and no Xvnc to start one"
      exit 1
    fi
    XDISP=""
    for d in $(seq 90 99); do
      [ -e /tmp/.X11-unix/X$d ] || [ -e /tmp/.X$d-lock ] || { XDISP=$d; break; }
    done
    if [ -z "$XDISP" ]; then
      echo "ERROR: no free display (:90 to :99) for a virtual X server"
      exit 1
    fi
    Xvnc :$XDISP -SecurityTypes None -localhost -geometry 1920x1200 > /dev/null 2>&1 &
    XVNC_PID=$!
    trap 'kill $XVNC_PID 2> /dev/null' EXIT
    export DISPLAY=:$XDISP
    sleep 2
  elif command -v vglrun > /dev/null && command -v xdpyinfo > /dev/null &&
       xdpyinfo 2> /dev/null | grep -qE "^ +(TIGERVNC|VNC-EXTENSION)$"; then
    # hardware OpenGL in a VNC session (detected as in setup.sh)
    VGLRUN=vglrun
  fi
fi

for sample in "${SAMPLES[@]}"; do
  IFS='|' read -r name runopts dispopts <<< "$sample"
  name=$(trim "$name"); runopts=$(trim "$runopts"); dispopts=$(trim "$dispopts")
  [ -n "$ONLY" ] && [ "$name" != "$ONLY" ] && continue
  if [ "$STEP" != shots ]; then
    generate "$name" "$runopts" || continue
  fi
  if [ "$STEP" != gen ]; then
    screenshots "$name" "$dispopts"
  fi
done
