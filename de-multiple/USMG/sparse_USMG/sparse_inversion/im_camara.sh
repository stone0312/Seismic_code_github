#!/bin/bash

dir=/home/fb/rtm2d_modelingonly/sparse_inversion/camaraman
file=$dir/$1

n1=256
n2=256
wbox=800
hbox=800

ximage	\
	n1=$n1 n2=$n2\
	xbox=50 ybox=50 legend=1 \
	hbox=$hbox wbox=$wbox\
	<$file title="$file"&
