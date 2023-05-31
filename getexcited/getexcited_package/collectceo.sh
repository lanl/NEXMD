#!/bin/bash

if [ -f ceo.out ]; then
    rm ceo.out
fi
grep -A $1 'Frequencies (eV) and Oscillator strengths (unitless)' md.out | tail -n+3 >> ceo.out
