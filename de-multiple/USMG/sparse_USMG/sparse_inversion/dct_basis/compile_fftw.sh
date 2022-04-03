#!/bin/bash

set -x

# RSG 90
# Using Intel(R) Parallel Studio XE 2013 Update 3 for Linux*
#source /home/fb/apps/intel2013/bin/compilervars.sh intel64

#PWI-200
#FFTW3=/home/app/fftw322intel

#RSG 90
# FFTW 3.2.2
FFTW3=/mnt/data5/apps/fftw322

FFTW3_INC=$FFTW3/include
LFFTW3_LIB="-L${FFTW3}/lib -lfftw3 -L${FFTW3}/lib -lfftw3f"


# FB INC.
FBINC=/home/fb/rtm2d_modelingonly/include

INC=" -I$FFTW3_INC -I$FBINC -I. -I../"

src=./$1.c
out=./$1.exe

CC=icc

$CC -std=c99 ${src} ${INC} ${LFFTW3_LIB} -lm -o ${out}
