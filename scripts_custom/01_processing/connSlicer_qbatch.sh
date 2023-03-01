#!/bin/bash
#
# Slices connectomes according to a give look up table (using the BIC batch)
# Should be happening in ${MICAPIPE}/functions/03_SC.sh, but sometimes needs to be repeated...
#
#	Inputs:
# 		$1 : sub ID
# 		$2 : FUNC_ID (DEPRECATED!)
#
# 2021 Mark C Nelson, McConnell Brain Imaging Centre, MNI, McGill
#------------------------------------------------------------------------------------------------------------------------------------

# Paths
 source /PATH/TO/YOUR/init.sh

# shared files
 lut_sc="${MICAPIPE}/parcellations/lut/lut_subcortical-cerebellum_mics.csv"

# subject-specific paths
 id=HC0"$1"
 DWI_DIR=${OUT_DIR}/sub-${id}/ses-01/proc_dwi

 for i in ${DWI_DIR}/sc/connectomes/* ; do

	echo "Connectome : ${i##*/}"
	PARC=$(echo ${i##*/} | cut -d '_' -f 3)
	WEIGHT=$(echo ${i##*/} | cut -d '_' -f 5 | cut -d '.' -f 1)

	echo "Parcellation: ${PARC}"
	echo "Weight: ${WEIGHT}"
	lut="${MICAPIPE}/parcellations/lut/lut_${PARC}_mics.csv"
	Rscript ${MICAPIPE}/functions/connectome_slicer.R --conn=$i --lut1=${lut_sc} --lut2=${lut} --mica=${MICAPIPE}
 done	# over connectomes

