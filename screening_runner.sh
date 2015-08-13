#!/bin/bash
#
# Flush disks if nobody is on the computer
#
# Ken O. Burtch
# CVS: $Header$


STARTM=$(date +%s.%N)

#-----------------------------------------------------------------

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

# Pathes

way=/home/mk/science/Fortran_programm/gain/src/	#path to fortran sources
way1=/home/mk/matlab/work/Новая\ папка\ \(2\)/			#path to matlab sources

MKLPATH=/opt/intel/Compiler/11.1/059/mkl/lib/64			#path to MKL library
MKLINCLUDE=/opt/intel/Compiler/11.1/059/mkl/include		#path to MKL includes

cd ${way}

# Band structure computing

#/home/mk/science/Fortran_programm/gain_with_adopt/band_structure.sh "${way1}" qw_fd3_1 "'seg=0','piezo=1','do_not_comp_optics','save',2.0e-9"

# Compilation and linking

export PATH=/usr/local/dislin/:$PATH
export DISLIN=/usr/local/dislin
export LD_LIBRARY_PATH=/usr/local/dislin/:$LD_LIBRARY_PATH

#cp /usr/local/dislin/g95/dislin.f90 /home/mk
#ifort -c /home/mk/dislin.f90 -L/usr/local/dislin -ldislin
#ifort -I/home/mk/ -g -traceback -heap-arrays -o gain.out ${way}gain.f90 ${way}Fermi_level.f90 ${way}polarization1.f90 ${way}int_matrix.f90 ${way}exchange.f90 ${way}dielectrical_const.f90 ${way}E_field.f90 ${way}concentration.f90 -L/usr/local/dislin/lib -ldislin -ldislnc -L$MKLPATH -I$MKLINCLUDE -lmkl_intel -lmkl_intel_thread -lmkl_core -liomp5 -lpthread

ifort  -L$MKLPATH -I$MKLINCLUDE -module ${way} -c ${way}fitting.f90  -lmkl_intel -lmkl_intel_thread -lmkl_core -lmkl_lapack -liomp5 -lpthread -lguide
ifort -L$MKLPATH -I$MKLINCLUDE -g -traceback -heap-arrays -o screening.out ${way}fitting.f90 ${way}screening.f90 ${way}Fermi_level.f90  -lmkl_intel -lmkl_intel_thread -lguide -lmkl_core -lmkl_lapack -liomp5 -lpthread -L/usr/local/dislin/lib -ldislin -ldislnc

# Executing

./screening.out
rm ./screening.out

#-----------------------------------------------------------------

STOPM=$(date +%s.%N)

printf "Elapsed:    %.3F\n"  $(echo "$STOPM - $STARTM"|bc )

# Cleanup
exit 0
# all is well
