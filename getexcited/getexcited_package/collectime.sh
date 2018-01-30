#!/bin/bash

if [ -f timing.out ]; then
    rm timing.out
fi
## for new code (nexmd) using intel and gnu compilers ##
tail -14 md.out | awk 'NR%2==0' | head -n-1 >> timing.out

## for new code (nexmd) using pgi compiler ##
#tail -17 md.out | awk 'NR%2==0' | head -n-2 >> timing.out

## for old code (naesmd) ##
#tail -5 md.out | head -n-4 >> timing.out
