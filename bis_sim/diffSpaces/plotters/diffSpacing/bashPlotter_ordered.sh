#!/bin/bash

declare -a arr=("1.0" "2.0" "3.0" "4.0" "5.0")
for f in *.xml; do
	for i in "${arr[@]}"; do
		python finalPlotter_diffSpacing_ordered.py "${f}" "${i}"  
	done;
done
