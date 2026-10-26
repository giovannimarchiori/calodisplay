echo ""

HOSTNAME=`hostname`
echo "Host name = ${HOSTNAME}"

# setup ROOT, unless it is already set up
if [[ -n $ROOTSYS ]]; then
    echo "ROOT already set up"
elif [[ $HOSTNAME = apcatlas01.in2p3.fr ]]; then
    source /usr/local/root-6.38.00/install/bin/thisroot.sh
elif [[ $HOSTNAME == lxplus*.cern.ch || $HOSTNAME == fcc-ironic*.cern.ch ]]; then
    # setup ROOT via the nightly key4hep stack
    source /cvmfs/sw-nightlies.hsf.org/key4hep/setup.sh
else
    # on my laptop
    setroot
fi
if command -v root-config >/dev/null 2>&1; then
    echo "Using ROOT $(root-config --version) from $(root-config --prefix)"
else
    echo "WARNING: ROOT not found (root-config is not in the PATH)"
    # return, not exit: the script is sourced, exit would close the shell
    return 1
fi

echo ""

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

# check that the calodisplay binaries exist, are up to date, and were
# built with a ROOT compatible with the one set up above
check_build() {
    local dist="${CALODISPLAY_DIR}/dist"
    local exe="${dist}/calodisplay" lib="${dist}/calodisplay-library.so"
    local runExample="calodisplay -g <geometry.root> -e <events.root> (see 'calodisplay -h' for all options)"
    if [[ ! -f $exe || ! -f $lib ]]; then
        echo "calodisplay is not compiled yet: run 'make -j8' in ${CALODISPLAY_DIR}"
        echo "Then run it with e.g.: ${runExample}"
        return
    fi
    local rootVersion=$(root-config --version) rootPrefix=$(root-config --prefix)

    # ROOT used for the build: written by the Makefile, or (for older builds)
    # deduced from the ROOT library directory recorded in the binary
    local buildVersion="" buildPrefix=""
    if [[ -f ${dist}/build-info.txt ]]; then
        buildVersion=$(sed -n 's/^ROOT_VERSION=//p' "${dist}/build-info.txt")
        buildPrefix=$(sed -n 's/^ROOT_PREFIX=//p' "${dist}/build-info.txt")
    elif command -v readelf >/dev/null 2>&1; then
        local rpath=$(readelf -d "$lib" | sed -n 's/.*R[UN]*PATH.*\[\(.*\)\]/\1/p' | cut -d: -f1)
        if [[ -n $rpath && $rpath -ef $(root-config --libdir) ]]; then
            buildVersion=$rootVersion
            buildPrefix=$rootPrefix
        elif [[ -n $rpath ]]; then
            buildPrefix=$(dirname "$rpath")
        fi
    fi

    if [[ -z $buildVersion ]]; then
        if [[ -n $buildPrefix ]]; then
            echo "WARNING: calodisplay was built with ROOT from ${buildPrefix},"
        else
            echo "WARNING: calodisplay was built with an unknown version of ROOT,"
        fi
        echo "         but ROOT ${rootVersion} from ${rootPrefix} is set up:"
        echo "         in ${CALODISPLAY_DIR} run 'make clean' and then 'make -j8'"
        echo "Then run it with e.g.: ${runExample}"
        return
    fi
    # different major.minor versions are not binary compatible
    if [[ ${buildVersion%.*} != ${rootVersion%.*} ]]; then
        echo "WARNING: calodisplay was built with ROOT ${buildVersion} but ROOT ${rootVersion} is set up:"
        echo "         in ${CALODISPLAY_DIR} run 'make clean' and then 'make -j8'"
        echo "Then run it with e.g.: ${runExample}"
        return
    fi
    if [[ $buildVersion != "$rootVersion" || ! $buildPrefix -ef $rootPrefix ]]; then
        echo "WARNING: calodisplay was built with ROOT ${buildVersion} from ${buildPrefix},"
        echo "         a different installation than ROOT ${rootVersion} from ${rootPrefix}."
        echo "         In case of problems, in ${CALODISPLAY_DIR} run 'make clean' and then 'make -j8'"
    fi

    # compiler used for the build (written by the Makefile) and the one that
    # make would use now: an incremental build with a different compiler mixes
    # them, and the library can then fail to link (e.g. the GCC of the system
    # and the one of key4hep)
    local buildCxx="" cxx=$($(root-config --cxx) --version 2>/dev/null | head -1)
    [[ -f ${dist}/build-info.txt ]] && buildCxx=$(sed -n 's/^CXX=//p' "${dist}/build-info.txt")
    if [[ -n $buildCxx && $buildCxx != "$cxx" ]]; then
        echo "WARNING: calodisplay was built with the compiler"
        echo "           ${buildCxx}"
        echo "         but the one set up now is"
        echo "           ${cxx}"
        echo "         to recompile it with this one, in ${CALODISPLAY_DIR} run 'make clean' and then 'make -j8'"
    fi
    # objects already compiled with different compilers (GCC)
    if command -v strings >/dev/null 2>&1; then
        local compilers=$(for o in "${CALODISPLAY_DIR}"/build/*.o; do
                              [[ -f $o ]] && strings -a "$o" | grep -m1 '^GCC:'
                          done | sort -u)
        if (( $(grep -c . <<< "$compilers") > 1 )); then
            echo "WARNING: the objects of calodisplay were compiled with different compilers:"
            sed 's/^/           /' <<< "$compilers"
            echo "         in ${CALODISPLAY_DIR} run 'make clean' and then 'make -j8'"
        fi
    fi

    # sources modified after the build
    local older=$exe
    [[ $lib -ot $exe ]] && older=$lib
    local modified=$(cd "${CALODISPLAY_DIR}" && find src -newer "$older" -type f ! -name '*~' | sort)
    if [[ -n $modified ]]; then
        echo "WARNING: calodisplay is not up to date, these files were modified after the build:"
        echo "$modified" | sed 's/^/           /'
        echo "         in ${CALODISPLAY_DIR} run 'make -j8'"
        echo "Then run it with e.g.: ${runExample}"
        return
    fi

    echo "calodisplay is up to date (built on $(date -r "$exe" '+%Y-%m-%d %H:%M') with ROOT ${buildVersion})"
    echo "Run it with e.g.: ${runExample}"
}
check_build
unset -f check_build

# in a VNC session OpenGL is rendered in software (slow): remind to use
# vglrun (VirtualGL) to render on the GPU, if available
check_vnc() {
    [[ -z $DISPLAY ]] && return
    # on macOS DISPLAY points to XQuartz, which is started by any X client
    # (e.g. xdpyinfo): skip the check, VirtualGL is for Linux servers anyway
    [[ $(uname) == Darwin ]] && return
    local vnc=0
    if command -v xdpyinfo >/dev/null 2>&1; then
        # VNC servers (TigerVNC, TurboVNC, ...) provide a VNC X extension
        xdpyinfo 2>/dev/null | grep -qE "^ +(TIGERVNC|VNC-EXTENSION)$" && vnc=1
    else
        # no xdpyinfo: look for an Xvnc server for this display
        local display=${DISPLAY#*:}
        display=${display%%.*}
        pgrep -f "Xvnc :${display}( |$)" >/dev/null 2>&1 && vnc=1
    fi
    [[ $vnc = 0 ]] && return
    echo ""
    if command -v vglrun >/dev/null 2>&1; then
        echo "You are in a VNC session: to render on the GPU instead of in software, run with vglrun, e.g.:"
        echo "    vglrun calodisplay -g <geometry.root> -e <events.root>"
    else
        echo "You are in a VNC session: OpenGL is rendered in software and can be slow"
        echo "(VirtualGL is not installed, see the README)"
    fi
}
check_vnc
unset -f check_vnc

echo ""
