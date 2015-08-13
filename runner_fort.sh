#!/bin/bash
#
# Flush disks if nobody is on the computer
#
# Ken O. Burtch
# CVS: $Header$

shopt -s -o nounset

# Global Declarations

declare -rx SCRIPT=${0##*/}
declare -rx ifort=”/opt/intel/fc/10.1.015/bin/ifort”

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

cd /home/mk/science/Fortran_programm/gain/src/

way=/home/mk/science/Fortran_programm/gain/src/
MKLPATH=/opt/intel/Compiler/11.1/059/mkl/lib/64
MKLINCLUDE=/opt/intel/Compiler/11.1/059/mkl/include

# Compilation and linking

export PATH=/usr/local/dislin/:$PATH
export DISLIN=/usr/local/dislin
export LD_LIBRARY_PATH=/usr/local/dislin/:$LD_LIBRARY_PATH

#cp /usr/local/dislin/g95/dislin.f90 /home/mk
#ifort -c /home/mk/dislin.f90 -L/usr/local/dislin -ldislin
#ifort -I/home/mk/ -g -traceback -heap-arrays -o gain.out ${way}gain.f90 ${way}Fermi_level.f90 ${way}polarization1.f90 ${way}int_matrix.f90 ${way}exchange.f90 ${way}dielectrical_const.f90 ${way}E_field.f90 ${way}concentration.f90 -L/usr/local/dislin/lib -ldislin -ldislnc -L$MKLPATH -I$MKLINCLUDE -lmkl_intel -lmkl_intel_thread -lmkl_core -liomp5 -lpthread

ifort -L$MKLPATH -I$MKLINCLUDE -g -traceback -heap-arrays -o gain.out ${way}gain.f90 ${way}Fermi_level.f90 ${way}polarization1.f90 ${way}int_matrix.f90 ${way}exchange.f90 ${way}dielectrical_const.f90 ${way}E_field.f90 ${way}concentration.f90 ${way}polyint.f90 ${way}f_f.f90 -lmkl_intel -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -L/usr/local/dislin/lib -ldislin -ldislnc

# Executing

./gain.out
rm ./gain.out

# Vizualization

# gnuplot ./gp1.gp
# evince ./gain_pic.pdf

# Cleanup
exit 0
# all is well
