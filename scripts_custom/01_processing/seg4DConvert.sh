#!/bin/bash
#
# This script is a quick workaround for the consistent failure of first_mult_bcorr in 01_proc-struc_volumetric.sh
# see first_run_all (https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FIRST/UserGuide)
#
#
# 2021 Mark C Nelson, McConnell Brain Imaging Centre, MNI, McGill
#------------------------------------------------------------------------------------------------------------------------------------

source /PATH/TO/YOUR/init.sh

for id in 17 30
do

# Define subject specific variables
subject_dir=${OUT_DIR}/sub-HC0${id}/ses-01
proc_struct=$subject_dir/proc_struct
dir_volum=$proc_struct/volumetric
T1nativepro_brain=${proc_struct}/HC0${id}_t1w_*mm_nativepro_brain.nii.gz
T1fast_seg=$proc_struct/first/HC0${id}_t1w_*mm_nativepro_all_fast_firstseg.nii.gz
T1fast_org=$proc_struct/first/HC0${id}_t1w_*mm_nativepro_all_fast_origsegs.nii.gz
T1_seg_subcortex=${dir_volum}/HC0${id}_t1w_0.8mm_nativepro_subcortical.nii.gz

Info "Converting subject-HC0${id}"
# Convert 4D segmentation file to 3D using fsl
 first_mult_bcorr -i ${T1nativepro_brain} -u ${T1fast_org}  -c ${T1fast_seg} -o ${T1_seg_subcortex}

 if [[ -f ${T1_seg_subcortex} ]]; then
     echo "${T1_seg_subcortex} successfully created"
 else
     echo "ISSUE with subject-HC0${id}"
 fi

done
