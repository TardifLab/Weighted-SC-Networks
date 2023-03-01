#!/bin/bash
#
# Manually changes individual bvals to desired value
#
# INPUTS:
# 	$1 : bvals file
# 	$2 : Target b values (array)
#
# OUTPUTS:
#	Modified bvals file saved in same dir as original
#
#
# Mark C Nelson (2021), McConnell Brain Imaging Centre, MNI, McGill University
#------------------------------------------------------------------------------------------

# Inputs (shift chops off input 1 and allows remaining inputs to be assigned to array)
 file_orig=$1
 shift
 values_target=("$@")

echo "Modifying bvals file at: $file_orig"
echo "Target bvals: ${values_target[@]}"

# Setup
BASE_DIR=$(dirname $file_orig)
file_mod=${BASE_DIR}/bvals_mod
bvals_string=$(cat $file_orig)					# load bvals to string
bvals_array=($bvals_string)					# convert to array
thr=50								# threshold for modifying bvalue

# Change values
for Io in "${!bvals_array[@]}"; do 							# indices in original array as outter loop var
	for Vt in "${values_target[@]}"; do 						# target values as inner loop var
		let Vo=bvals_array[$Io] 						# original value
		let diff=$Vo-$Vt 							# compute difference
		if [ ${diff##*[+-]} -le $thr ]; then
#			echo "* Pass : $Io.    Changing $Vo --> $Vt"
			bvals_array[$Io]=$Vt;
			break;
		fi
	done
done

# Write file
echo ${bvals_array[@]} > ${file_mod}



