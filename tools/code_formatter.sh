#!/bin/bash

# Fortran code formatter
# Both fprettify and findent are used
# Both can be installed through conda or pip, but the packaged fprettify is rather old
# pip install findent fprettify or conda install -c conda-forge findent fprettify
# It is recommanded to install the development version of fprettify
# git clone https://github.com/pseewald/fprettify.git && cd fprettify && pip install .
# If you install fprettify from package managers, two patchs are needed due to bugs https://github.com/pseewald/fprettify/issues/112 and https://github.com/pseewald/fprettify/issues/58
# Download the patch here https://github.com/pseewald/fprettify/commit/bf41457c782997e378144297e764092d09a3fa5f.patch/ and https://github.com/pseewald/fprettify/commit/9e0108c34f651b9444582715d281672e59fa0d9c.patch
# and patch file $your_python_library_path/site-packages/fprettify/__init__.py
# A simple way to determine the path
# Run fprettify -i 4 --case 1 0 0 0 src/sqm/constants.F90 and it will stuck forver
# Press ctrl-c to force quit and you will see the path of the source file in the error message
# If you use the development version, only the first patch is needed

# if any part in pipeline fails, return error status
set -o pipefail

# if the line length is logner than default 132, fprettify will not format the line
# 300 should be very safe
fprettify_cmd="fprettify -i 4 -l 300 --case 1 0 0 0"
findent_cmd="findent -i4 -Rr"
real_to_double='s/_REAL_ function/double precision function/g'
double_to_real='s/double precision function/_REAL_ function/g'
files_with_errors=""

# if CLI arguments are not provided, format all source files in src
if [[ -z $@ ]]; then
    files=$(find src -name '*.F90')
else
    files=$@
fi
for i in $files
do 
    echo $i
    if grep -q "_REAL_ function" $i; then
        echo "Function defined with _REAL_ preprecessor found."
        # findent doesn't work with function like "_REAL_ function"
        processed=$(cat $i | sed "$real_to_double" | $fprettify_cmd | $findent_cmd | sed "$double_to_real")
        status=$?
    else
        processed=$(cat $i | $fprettify_cmd | $findent_cmd)
        status=$?
    fi
    if [[ $status -eq 0 ]]; then
        echo "$processed" > $i
    else
        files_with_errors="$files_with_errors$i "
        echo "Error encountered when processing file $i"
    fi
done
if [[ ! -z $files_with_errors ]]; then
    echo "The folowwing files were not succcessfully formatted: $files_with_errors"
fi
