#!/bin/bash

set -x

# RSG 90
# Using Intel(R) Parallel Studio XE 2013 Update 3 for Linux*
#source /home/fb/apps/intel2013/bin/compilervars.sh intel64

Platform=90
#Platform=200

if [ $Platform = 90 ]
then
	# MPICH 2.1.4
	MPICH214=/home/fb/apps/mpich214intel2013
	MPICH=$MPICH214
	FFTW3=/mnt/data5/apps/fftw322
fi

if [ $Platform = 200 ]
then
	# MPICH 3.0.3
	MPICH303=/usr/local
	MPICH=$MPICH303
	FFTW3=/home/app/fftw322intel
fi

# MPICH
export PATH=$MPICH/bin:$PATH
MPICH_INC=$MPICH/include
MPICH_LIB=$MPICH/lib

# FFTW
FFTW3_INC=$FFTW3/include
LFFTW3_LIB="-L${FFTW3}/lib -lfftw3 -L${FFTW3}/lib -lfftw3f"

# Numerical Recipe 2
NR2_INC=/home/fb/numerical_recipe_code/ANSI_C

# Package dir.
RTMROOT=/home/fb/rtm2d_modelingonly
DoubleBeam_INC=${RTMROOT}/double_beam_forming_2d/newsrc
Intp_INC=${RTMROOT}/main/interpolation_B_spline
TTINV_INC=${RTMROOT}/traveltimeinversion/main
#DoubleBeam_INC_OLD=${RTMROOT}/double_beam_forming_2d/src

INC=" -I$MPICH_INC -I$FFTW3_INC -I$NR2_INC -I$RTMROOT/header -I$RTMROOT/include -I$DoubleBeam_INC -I$Intp_INC -I$TTINV_INC -I. "

src=./$1.c
out=./$1.exe

#CC=icc
CC=${MPICH}/bin/mpicc

$CC -std=c99 -openmp -static ${src} ${INC} ${LFFTW3_LIB} -lm -o ${out}
#$CC -std=c99 -openmp ${src} ${INC} ${LFFTW3_LIB} -lm -o ${out}
#$CC -std=c99 -openmp ${src} ${INC} -lm -o ${out}
