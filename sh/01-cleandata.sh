#!/bin/bash

direxpr="Data/Expression"

for t in Healthy Basal Her2 LumA LumB
do 
	python py/01-clean-data.py ${direxpr}/${t}-All.txt
done