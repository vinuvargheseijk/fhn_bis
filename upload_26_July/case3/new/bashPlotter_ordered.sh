#!/bin/bash

declare -a arr=("1.0" "2.0" "3.0" "4.0" "5.0")
for i in "${arr[@]}"; do
	python finalPlotter_model2_ordered.py trial.xml "${i}"  
done
