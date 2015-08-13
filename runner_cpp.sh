#!/bin/bash
#
# Flush disks if nobody is on the computer
#
# Ken O. Burtch
# CVS: $Header$

shopt -s -o nounset

# Global Declarations

declare -rx SCRIPT=${0##*/}
declare -rx ifort=”/opt/intel/Compiler/11.1/059/bin/ia32/ifort”

# Sanity checks
# SCRIPT is the name of this script
# the who command - man 1 who
# the sync command - man 1 sync
# the wc command - man 1 wc

if test -z “$BASH” ; then
printf “$SCRIPT:$LINENO: please run this script with the BASH shell\n” >&2
exit 192
fi

#if test ! -x “$ifort” ; then
#printf “$SCRIPT:$LINENO: the command $ifort is not available – aborting\n” >&2
#exit 192
#fi

# Arguments processing

way=/home/mk/science/Cpp_programms/gain/src/
cd /home/mk/science/Cpp_programms/gain/src/

# Compilation and linking

g++ -Wall -o gain.out ${way}main_SBE.cpp ${way}gain.cpp ${way}polarization1.cpp ${way}int_matrix.cpp ${way}exchange.cpp ${way}Fermi_level.cpp ${way}E_field.cpp

# Executing

./gain.out
rm ./gain.out

# Vizualization

# gnuplot ./gp1.gp
# evince ./gain_pic.pdf

# Cleanup
exit 0
# all is well
