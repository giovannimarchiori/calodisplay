# setup ROOT
HOSTNAME=`hostname`
echo $HOSTNAME
if [[ $HOSTNAME = apcatlas01.in2p3.fr ]]; then
    source /usr/local/root-6.38.00/install/bin/thisroot.sh
elif [[ $HOSTNAME == lxplus*.cern.ch || $HOSTNAME == fcc-ironic*.cern.ch ]]; then
    # expect root to be setup via key4hep
    echo "ROOT should be setup via key4hep or LCG"
else
    # on my laptop
    setroot
fi

# directory containing this script, so that it can be sourced from anywhere
export CALODISPLAY_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

# append a directory to a path variable, unless it is already there
# (and without an empty entry if the variable was not set)
append_path() {
    local var=$1 dir=$2
    local current=${!var}
    case ":${current}:" in
        *":${dir}:"*) ;;
        *) export "${var}=${current:+${current}:}${dir}" ;;
    esac
}
append_path DYLD_LIBRARY_PATH "${CALODISPLAY_DIR}/dist"
append_path LD_LIBRARY_PATH "${CALODISPLAY_DIR}/dist"
append_path PATH "${CALODISPLAY_DIR}/dist"
unset -f append_path

