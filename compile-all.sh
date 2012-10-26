#!/bin/bash

CC="gcc"
CCFLAGS="-lm"
F77="gfortran"
F77FLAGS="-w"


cd elast_source
if [ "$1" == clean ]; then
rm elast
else
$CC $CCFLAGS -o elast elast.c
fi
cd ..


cd gosia_source
if [ "$1" == clean ]; then
rm gosia_20110524.2
else
$F77 $F77FLAGS -o gosia_20110524.2 gosia_20110524.2.f
fi
cd ..

./setup.sh

