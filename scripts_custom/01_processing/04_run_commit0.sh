#!/bin/bash
#
# This script uses tools from MRtrix3 to prepare subject tractograms for COMMIT by removing all streamlines which do not connect 2 gray matter regions.
#
# 1. Computes a single connectome using the parcellation given by $parc_name
#
# 2. Reconstructs the tractogram using the node assignments from this single parcellated connectome excluding any streamlines that fail
#    to connect to 2 gray matter nodes.
#
#
# Inputs:
#       $1 : subject (i.e. 01, 02, ... , 50)
#       $2 : function id (DEPRECATED)
#
# 2022 Mark C Nelson, McConnell Brain Imaging Centre, MNI, McGill
#------------------------------------------------------------------------------------------------------------------------------------

# options & setup
  source /YOUR/PATH/TO/init.sh
  id="HC0$1"
  subject=sub-${id}
  tracts=5M
  threads=6
  SES=ses-01
  parc_name=schaefer-100

# Subject dirs
  subject_dir=${OUT_DIR}/${subject}/${SES}
  proc_dwi=${subject_dir}/proc_dwi
  proc_sc=${proc_dwi}/sc
  tmp=${proc_sc}/tmp
  dir_connectome_tmp=${proc_sc}/connectomes/tmp
  dir_volum=${subject_dir}/proc_struct/volumetric
  util_lut=${MICAPIPE}/parcellations/lut
  dir_warp=${subject_dir}/xfms

# Files
  tractogram=${proc_sc}/DWI_tractogram_${tracts}.tck
  filtertck=${proc_sc}/DWI_tractogram_${tracts}_filt-gm.tck
  assignments=${dir_connectome_tmp}/${id}_${tracts}_full_${parc_name}_nodeassignments.txt
  lut="${util_lut}/lut_${parc_name}_mics.csv"
  lut_sc="${util_lut}/lut_subcortical-cerebellum_mics.csv"
  fod=${tmp}/${id}_space-dwi_model-CSD_map-FOD_desc-wmNorm.nii.gz
  dwi_cere_NL=${proc_dwi}/${id}_dwi_cerebellum_NL.nii.gz
  dwi_subc_NL=${proc_dwi}/${id}_dwi_subcortical_NL.nii.gz
  dwi_all=$tmp/${id}_${parc_name}-full_dwi.nii.gz                                                                         	# Segmentation in dwi space
  seg=$dir_volum/${id}_t1w_0.8mm_nativepro_${parc_name}.nii.gz 									# Segmentation in native space

# transforms
  mat_dwi_affine=${dir_warp}/${id}_dwi_to_nativepro_0GenericAffine.mat
  dwi_SyN_str=${dir_warp}/${id}_space-dwi_from-dwi_to-dwi_mode-image_desc-SyN_
  dwi_SyN_warp=${dwi_SyN_str}1Warp.nii.gz
  dwi_SyN_affine=${dwi_SyN_str}0GenericAffine.mat

# Check inputs
  if [ ! -f $lut_sc ]; then echo "Can't find the subcortical-cerebellar LUT!"; exit; fi
  if [ ! -f $fod ]; then echo "Can't find the wm fod for registration of parcellations!"; exit; fi
  if [ ! -f $tractogram ]; then echo "Subject $id doesn't have a tractogram!"; exit; fi
  if [ ! -f $mat_dwi_affine ]; then echo "Subject $id doesn't have a dwi to native affine transform!"; exit; fi
  if [ ! -f $dwi_SyN_warp ]; then echo "Subject $id doesn't have a non-linear warp in dwi space!"; exit; fi
  if [ ! -f $dwi_SyN_affine ]; then echo "Subject $id doesn't have a non-linear affine to dwi space!"; exit; fi
  if [ ! -f $dwi_cere_NL ]; then echo "Subject $id doesn't have cerebellar segmentation registered to dwi space!"; exit; fi
  if [ ! -f $dwi_subc_NL ]; then echo "Subject $id doesn't have a subcortical segmentation registered to dwi space!"; exit; fi

# Begin processing
  echo "*------------------- Begin tractogram RECON for subject: $id ---------------------*"

# Timer
  aloita=$(date +%s)

# -------------------------------------------------------------------------------------------------------------------------------------------------
# 		FILTERING ALL STREAMLINES WHICH DO NOT CONNECT TWO GRAY MATTER NODES

## 1. Get node assignments for all streamlines

  echo "Generating Connectome in Parcellation: ** $parc_name **"

  if [[ ! -d ${dir_connectome_tmp} ]]; then mkdir -p ${dir_connectome_tmp} ; fi


# Take cortical parcellation into DWI space & add subcortical rois (non-linear)
  if [ ! -f $dwi_all ]; then echo "Running non-linear transformation for full segmentation to dwi space (incl. subcortex & cerebellum)"
  	dwi_cortex=$tmp/${id}_${parc_name}-cor_dwi.nii.gz
        dwi_cortexSub=$tmp/${id}_${parc_name}-sub_dwi.nii.gz
      # Registeration
        antsApplyTransforms -d 3 -e 3 -i $seg -r $fod -n GenericLabel -t $dwi_SyN_warp -t $dwi_SyN_affine -t [$mat_dwi_affine,1] -o $dwi_cortex -v -u int
      # Remove the medial wall
        for i in 1000 2000; do fslmaths $dwi_cortex -thr $i -uthr $i -binv -mul $dwi_cortex  $dwi_cortex ; done
      # add subcortical rois
        fslmaths $dwi_cortex -binv -mul $dwi_subc_NL -add $dwi_cortex $dwi_cortexSub -odt int
      # add cerebellar rois
        fslmaths $dwi_cortex -binv -mul $dwi_cere_NL -add $dwi_cortexSub $dwi_all -odt int
  else echo "File already exists: $dwi_all"; fi


# Compute connectome
  nom=$dir_connectome_tmp/${id}_${tracts}_unfiltered_full_${parc_name}
  conn="${nom}_nos.txt"
  if [ ! -f $conn ]; then echo "Connectome: $conn"
  tck2connectome -nthreads $threads $tractogram $dwi_all $conn \
                 -scale_invnodevol -out_assignments $assignments -quiet
  Rscript ${MICAPIPE}/functions/connectome_slicer.R --conn=$conn --lut1=${lut_sc} --lut2=${lut} --mica=${MICAPIPE}
  else echo "$conn ALREADY EXISTS!" ; fi


# -------------------------------------------------------------------------------------------------------------------------------------------------
## 2. Reconstruct filtered tractogram

  echo "Beginning Reconstruction of Tractogram: $filtertck"

  # This super annoying step required all node indices be giving manually for all subjects. Nodes not indicated would have all of their streamlines removed!
  # Be careful... Come up with a smarter solution
  connectome2tck $tractogram $assignments $filtertck -keep_self -files single \
    -nodes 10,11,12,13,16,17,18,26,49,50,51,52,53,54,58,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,139,141,143,149,150,1001,1002,1003,1004,1005,1006,1007,1008,1009,1010,1011,1012,1013,1014,1015,1016,1017,1018,1019,1020,1021,1022,1023,1024,1025,1026,1027,1028,1029,1030,1031,1032,1033,1034,1035,1036,1037,1038,1039,1040,1041,1042,1043,1044,1045,1046,1047,1048,1049,1050,2001,2002,2003,2004,2005,2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019,2020,2021,2022,2023,2024,2025,2026,2027,2028,2029,2030,2031,2032,2033,2034,2035,2036,2037,2038,2039,2040,2041,2042,2043,2044,2045,2046,2047,2048,2049,2050 \
    -exclusive

# -------------------------------------------------------------------------------------------------------------------------------------------------
## 3. Running SIFT2 on filtered tractogram

# Files
  fod_wmN=${proc_dwi}/${id}_wm_fod_norm.mif
  weights_sift2=${proc_sc}/SIFT2_${tracts}_filt-gm.txt

# Run SIFT2
  if [[ ! -f $weights_sift2 ]]; then echo "Running SIFT2 on filtered tractogram"
        tcksift2 -nthreads $threads $filtertck $fod_wmN $weights_sift2
  else echo "Subject ${id} has SIFT2 weights"; fi




# Notification of completion
  eri=$(echo "$lopuu - $aloita" | bc)
  eri=$(echo print $eri/60 | perl)

  echo "*------------ Connectome processing for $id complete -------------*"
  echo "  Proc time				: $eri minutes"

# Remove temporary Connectome & NodeAssignments files
  echo "  Deleting connectome & node assignment file "
  rm -rf $dir_connectome_tmp
