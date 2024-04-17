# setup ROOT
echo $HOSTNAME
if [[ $HOSTNAME = apcatlas01.in2p3.fr ]]; then
    source /usr/local/root-6.30.06/install/bin/thisroot.sh
elif [[ $HOSTNAME == lxplus*.cern.ch ]]; then
    # expect root to be setup via key4hep
    echo "ROOT should be setup via key4hep or LCG"
else
    # on my laptop
    setroot
fi
    
# root should be setup via key4hep
export DYLD_LIBRARY_PATH=${DYLD_LIBRARY_PATH}:${PWD}/dist
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${PWD}/dist
export PATH=${PATH}:${PWD}/dist

