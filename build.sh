#!/bin/bash

if [ "$1" = "all" ]; then
	all_folders=("semi-honest_linear_logistic" "semi-honest_NN" "malicious_linear_logistic" "malicious_NN" "benchmarks/SecureML" "benchmarks/client-aided_SecureML" "benchmarks/client-aided_SecureML_NN")
	ORIGINAL_DIR=$(pwd)
	for one_folder in "${all_folders[@]}"; do
		cp -r data "$one_folder"
		cd "$one_folder"
		sudo mkdir -p build
		cd build
		sudo cmake ..
		sudo make
		cd "$ORIGINAL_DIR"
	done
else
	cp -r data $1
	cd $1
	sudo mkdir build
	cd build
	sudo cmake ..
	sudo make
fi


