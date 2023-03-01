#!/bin/bash
#
#	Shell wrapper for running the original COMMIT implementation (regularization: commit.solvers.non_negative)
# 	IC weights obtained for streamlines will be used to compute connectome edge weights
#
# Inputs:
#       $1 : subject (i.e. 01, 02, ... , 50)
#       $2 : function id (DEPRECATED)
#
#  * * Requires Python & virtual environment with necessary dependencies to be setup already!
#
# 2021 Mark C Nelson, McConnell Brain Imaging Centre, MNI, McGill
#------------------------------------------------------------------------------------------------------------------------------------

 source /YOUR/PATH/TO/init.sh
 ID="HC0$1"

# Subject dirs
 subProc_dir=${OUT_DIR}/sub-${ID}/ses-01/proc_dwi
 commit_dir=${subProc_dir}/commit
 tmp=${commit_dir}/tmp
 if [ ! -d ${tmp} ]; then mkdir -p ${tmp} ;  fi

# Files
 tractogram=${subProc_dir}/sc/DWI_tractogram_5M_filt-gm.tck
 wm_fod_mif=${subProc_dir}/${ID}_wm_fod_norm.mif
 wm_fod_nii=${tmp}/${ID}_wm_fod_norm.nii.gz
 f_5tt=${subProc_dir}/${ID}_dwi_5tt.nii.gz
 wm_mask=${tmp}/${ID}_dwi_wm_mask.nii.gz
 dwi_b0=${subProc_dir}/${ID}_dwi_b0.nii.gz
 dwi_corr_mif=${subProc_dir}/${ID}_dwi_corr.mif
 dwi_corr_nii=${tmp}/${ID}_dwi_corr.nii.gz

 bvecs_eddy=${subProc_dir}/eddy/dwi_post_eddy.eddy_rotated_bvecs
 bvals_orig=${tmp}/bvals_orig
 bvals_mod=${tmp}/bvals_mod
 bvecs=${tmp}/bvecs
 scheme=${commit_dir}/commit.scheme

# options
 bvals_target=("0" "300" "700" "2000")

#----------- FILE CHECK ------------#
if [ ! -f ${tractogram} ];                                      then echo -e "Error: Input tractogram doesn't exist! ...\n EXITING \n"; exit 1; fi
if [ ! -f ${wm_fod_nii} ] && [ ! -f ${wm_fod_mif} ];            then echo -e "Error: Subject $ID doesn't have a WM fod...\n EXITING \n"; exit 1; fi
if [ ! -f ${f_5tt} ] && [ ! -f ${wm_mask} ];                    then echo -e "Error: Subject $ID doesn't have a WM mask or a 5TT file for extraction...\n EXITING \n"; exit 1; fi
if [ ! -f ${dwi_b0} ];                                          then echo -e "Error: Subject $ID doesn't have a B0 DWI image (use alternate reference)...\n EXITING \n"; exit 1; fi
if [ ! -f ${dwi_corr_mif} ] && [ ! -f ${dwi_corr_nii} ];        then echo -e "Error: Subject $ID doesn't have corrected DWI...\n EXITING \n"; exit 1; fi
if [ ! -f ${bvecs_eddy} ];                                      then echo -e "Error: Subject $ID doesn't have an eddy corrected b-vectors file...\n EXITING \n"; exit 1; fi

#-------- FILE CONVERIONS ----------#

# DWI: convert .mif to .nii for corr diff data
 if [[ ! -f ${dwi_corr_nii} ]] ; then   mrconvert ${dwi_corr_mif} ${dwi_corr_nii}
 else                                   echo -e "${dwi_corr_nii} EXISTS...\n" ; fi

# wm fod: convert .mif to .nii
 if [[ ! -f ${wm_fod_nii} ]] ; then   mrconvert ${wm_fod_mif} ${wm_fod_nii}
 else                                   echo -e "${wm_fod_nii} EXISTS...\n" ; fi

# wm mask: extract from 5tt & convert to 3D image
 if [[ ! -f ${wm_mask} ]] ; then   mrconvert -coord 3 2 -axes 0,1,2 ${f_5tt} ${wm_mask}
 else                                   echo -e "${wm_mask} EXISTS...\n" ; fi

#-------- SCHEME PREP ----------#
 if [[ ! -f ${scheme} ]] ; then echo " Creating scheme file  "

	# extract bvals from corrected DWI.mif file
   	 mrinfo -export_grad_fsl ${bvecs} ${bvals_orig} ${dwi_corr_mif}

	# ensure B values are in correct format for commit
	 bvalsModify.sh ${bvals_orig} "${bvals_target[@]}"

	# Combine with eddy rotated bvecs file to create scheme for each subject
  	 bvecs_bvals2camino.pl -vec ${bvecs_eddy} -val ${bvals_mod} -o ${scheme}
 else
 	 echo "${scheme} EXISTS!!"
 fi

#---------- COMMIT CALL -------------#

 echo " *.*.*----------------------- Running COMMIT-PREP on $ID... ----------------------*.*.*"
 echo "... TRACTOGRAM            : $tractogram "
 echo "... OUTPUT                : $commit_dir "

 echo "Activating python venv at : $pyvenv_commit  "
 source ${pyvenv_commit}/bin/activate

 echo "Calling python script to run COMMIT1"
 python3 ${scripts}/04_run_commit2.py $ID $tractogram $commit_dir 			# call to custom vanilla COMMIT script
