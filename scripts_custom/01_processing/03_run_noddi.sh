#!/bin/bash
#
# Runs NODDI on preprocessed diffusion data (micapipe preprocessing modules are run before this!)
# Requires Matlab, AMICO & MRtrix! 
#	See: AMICO_process.m 
#
# INPUTS:
# 	$1 : SUBJECT (i.e. 01, 02, ... , 50)
# 	$2 : FUNCTION ID (DEPRECATED)
#
# OUTPUTS:
#	See ${NODDI_DIR}, ${NODDI_DIR}/AMICO & ${NODDI_DIR}/kernals
#
#
# 2021 Mark C Nelson, McConnell Brain Imaging Centre, MNI, McGill
#------------------------------------------------------------------------------------------------------------------------------------

# load dependencies & variables
  source /YOUR/PATH/TO/init.sh
  bvals_target=("0" "300" "700" "2000") 				# bvals in your data (real bvals will be adjusted to these idealised values)

# define subject files
  ID="HC0$1"
  DWI_DIR=${OUT_DIR}/sub-${ID}/ses-01/proc_dwi
  NODDI_DIR=${DWI_DIR}/noddi
  dwi_corr_mif=${DWI_DIR}/${ID}_dwi_corr.mif
  dwi_corr_nii=${NODDI_DIR}/${ID}_dwi_corr.nii
  bvecs_eddy=${DWI_DIR}/eddy/dwi_post_eddy.eddy_rotated_bvecs
  bvals_orig=${NODDI_DIR}/bvals_orig
  bvals_mod=${NODDI_DIR}/bvals_mod
  bvecs=${NODDI_DIR}/bvecs
  scheme=${NODDI_DIR}/NODDI.scheme
  dwi_mask_orig=${DWI_DIR}/${ID}_dwi_mask.nii.gz
  dwi_mask_new=${NODDI_DIR}/${ID}_dwi_mask.nii.gz
  dwi_mask_unzip=${NODDI_DIR}/${ID}_dwi_mask.nii

# mkdir if necessary
  if [ ! -d ${NODDI_DIR} ]; then 	mkdir ${NODDI_DIR} ; fi

# extract bvals from corrected DWI.mif file
  if [[ ! -f ${bvals_orig} ]] ; then 	mrinfo -export_grad_fsl ${bvecs} ${bvals_orig} ${dwi_corr_mif}
  else 			      		echo "${bvals_orig} EXISTS!!" ; fi

# ensure B0s are EXACTLY 0 (required for NODDI code)
  if [[ ! -f ${bvals_mod} ]] ; then       bvalsModify.sh ${bvals_orig} "${bvals_target[@]}"
  else                                    echo "${bvals_mod} EXISTS!!" ; fi

# Combine with eddy rotated bvecs file to create scheme for each subject
  if [[ ! -f ${scheme} ]]; then 		bvecs_bvals2camino.pl -vec ${bvecs_eddy} -val ${bvals_mod} -o ${scheme}
  else 			      		echo "${scheme} EXISTS!!" ; fi

# convert .mif to .nii for corr diff data
  if [[ ! -f ${dwi_corr_nii} ]] ; then 	mrconvert ${dwi_corr_mif} ${dwi_corr_nii}
  else                            	echo "${dwi_corr_nii} EXISTS!!" ; fi

# gunzip diff brain mask (should QC these masks first?)
  if [[ ! -f ${dwi_mask_unzip} ]] ; then cp ${dwi_mask_orig} ${dwi_mask_new} ; gunzip ${dwi_mask_new}
  else                                	echo "${dwi_mask_new} EXISTS!!" ; fi

#--------------------------------------------------------------------------------------------------------------------------------------------------#
##----- Run NODDI (requires Matlab) ---------------------- ---------------------------------------------------------------------------------------##

cd ${NODDI_DIR} 								# easier to be in dir, AMICO_process sets up out paths unfavorably

# call Ilanas dope wrapper (would be smarter not to recompute kernels for each subj...)
 matlab -nodisplay -nojvm -r "try; AMICO_process( './' , './' , '$(basename ${dwi_corr_nii})' , '$(basename ${dwi_mask_unzip})' , '$(basename ${scheme})' ); catch; end; quit"

# Double check output. Clear trash if completed
 if [[ ! -f "${NODDI_DIR}/AMICO/NODDI/FIT_ICVF.nii" ]] ; then 	echo "ERROR: NO ICVF FOUND FOR ${ID}"
 else 							     	rm -r ${dwi_mask_unzip} ${dwi_corr_nii} ${NODDI_DIR}/kernels ; fi
