#!/bin/bash

dir=/mnt/data1/fb/datatmp/fivelayermodel/gradient.OffMax6000m.Niteration30
file=$dir/$1

n1=601
n2=1401
wbox=1400
hbox=600

ximage	\
	n1=$n1 n2=$n2\
	xbox=50 ybox=50 legend=1 \
	hbox=$hbox wbox=$wbox\
	<$file title="$file"&
