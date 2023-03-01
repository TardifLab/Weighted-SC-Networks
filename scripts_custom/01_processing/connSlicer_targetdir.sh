#!/bin/bash
#
# Slices connectomes according to a give look up table
# Should be happening in ${MICAPIPE}/functions/03_SC.sh, but sometimes needs to be repeated...
#
# Inputs:
# 	$1 : target dir (needs to be adjusted for dir structure
#
# 2021 Mark C Nelson, McConnell Brain Imaging Centre, MNI, McGill
#------------------------------------------------------------------------------------------------------------------------------------

# Paths
 source /PATH/TO/YOUR/init.sh

# shared files
 lut_sc="${MICAPIPE}/parcellations/lut/lut_subcortical-cerebellum_mics.csv"

 targetdir=$1
 echo "$targetdir"

 for parc in $(ls ${targetdir}) ; do
 	echo "Parcellation: $parc"
 	lut="${MICAPIPE}/parcellations/lut/lut_${parc}_mics.csv"

	for weight in $(ls ${targetdir}/${parc}) ; do
		echo "... Weight: $weight"

 		for i in $(ls ${targetdir}/${parc}/${weight}) ; do
			echo "... ... Connectome : ${i##*/}"
			conn=${targetdir}/${parc}/${weight}/${i}
			Rscript ${MICAPIPE}/functions/connectome_slicer.R --conn=$conn --lut1=${lut_sc} --lut2=${lut} --mica=${MICAPIPE}

		done 	# over conns
 	done		# over weights
 done 			# over dirs
