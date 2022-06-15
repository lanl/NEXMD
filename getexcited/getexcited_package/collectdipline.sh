#!/bin/bash

if [ -f dipline.out ]; then
    rm dipline.out
fi

freq=$1

grep -n ' Begin classical propagation step #' md.out | awk '{print(substr($1, 1, length($1)-1) "   " $7}' | awk -v myvar="$freq" 'NR == 1 || NR % myvar == 0' >> dipline.out)
