#!/bin/bash
#
# Flush disks if nobody is on the computer
#
# Ken O. Burtch
# CVS: $Header$

# Arguments processing

cd "${1}"

matlab_exec=matlab
X="${2}(${3})"
echo ${X} > matlab_command_${3}.m
cat matlab_command_${3}.m

# Executing

${matlab_exec} -nojvm -nodisplay -nosplash < matlab_command_${3}.m

# Cleanup

rm matlab_command_${3}.m

exit 0
# all is well

