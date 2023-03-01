#!/bin/bash
#
# DWI POST structural TRACTOGRAPHY processing with bash:
#
# POST processing workflow for diffusion MRI TRACTOGRAPHY.
#
# This workflow makes use of MRtrix3
#
# Atlas an templates are avaliable from:
#
# https://github.com/MICA-MNI/micaopen/templates
#
#   ARGUMENTS order:
#   $1 : BIDS directory
#   $2 : participant
#   $3 : Out parcDirectory
#
# source code from micalab: https://github.com/MICA-MNI/micapipe
#
# edited 2021 Mark C Nelson, McConnell Brain Imaging Centre, MNI, McGill
#-------------------------------------------------------------------------------

BIDS=$1
id=$2
out=$3
SES=$4
PROC=$5
nocleanup=$6
tracts=$7
autoTract=$8
threads=$9
here=$(pwd)

#------------------------------------------------------------------------------#
# qsub configuration
  if [ "$PROC" = "qsub-MICA" ] || [ "$PROC" = "qsub-all.q" ];then
	export MICAPIPE=/data_/mica1/01_programs/micapipe
    	source ${MICAPIPE}/functions/init.sh;
  fi

# source utilities
source $MICAPIPE/functions/utilities.sh

# Assigns variables names
bids_variables $BIDS $id $out $SES

# Check inputs: DWI post TRACTOGRAPHY
# structural data
  T1str_nat=${id}_t1w_${res}mm_nativepro
  T1_seg_cerebellum=${dir_volum}/${T1str_nat}_cerebellum.nii.gz
  T1_seg_subcortex=${dir_volum}/${T1str_nat}_subcortical.nii.gz

# transforms
  mat_dwi_affine=${dir_warp}/${id}_dwi_to_nativepro_0GenericAffine.mat
  dwi_SyN_str=${dir_warp}/${id}_space-dwi_from-dwi_to-dwi_mode-image_desc-SyN_
  dwi_SyN_warp=${dwi_SyN_str}1Warp.nii.gz
  dwi_SyN_Invwarp=${dwi_SyN_str}1InverseWarp.nii.gz
  dwi_SyN_affine=${dwi_SyN_str}0GenericAffine.mat

# dwi
  fod_wmN=${proc_dwi}/${id}_wm_fod_norm.mif
  dwi_5tt=${proc_dwi}/${id}_dwi_5tt.nii.gz
  dwi_b0=${proc_dwi}/${id}_dwi_b0.nii.gz
  dwi_mask=${proc_dwi}/${id}_dwi_mask.nii.gz
  dwi_gmwmi=${proc_dwi}/${id}_space-dwi_desc-gmwmi-mask.mif

# Check inputs
  if [ ! -f $fod_wmN ]; then Error "Subject $id doesn't have WM FOD:\n\t\tRUN -proc_dwi"; exit; fi
  if [ ! -f $dwi_gmwmi ]; then Error "Subject $id doesn't have gray matter white matter mask:\n\t\tRUN -proc_dwi"; exit; fi
  if [ ! -f $dwi_b0 ]; then Error "Subject $id doesn't have dwi_b0:\n\t\tRUN -proc_dwi"; exit; fi
  if [ ! -f $mat_dwi_affine ]; then Error "Subject $id doesn't have an affine mat from T1nativepro to DWI space:\n\t\tRUN -proc_dwi"; exit; fi
  if [ ! -f $dwi_5tt ]; then Error "Subject $id doesn't have 5tt in dwi space:\n\t\tRUN -proc_dwi"; exit; fi
  if [ ! -f $T1_seg_cerebellum ]; then Error "Subject $id doesn't have cerebellar segmentation:\n\t\tRUN -post_structural"; exit; fi
  if [ ! -f $T1_seg_subcortex ]; then Error "Subject $id doesn't have subcortical segmentation:\n\t\tRUN -post_structural"; exit; fi
  if [ ! -f $dwi_mask ]; then Error "Subject $id doesn't have DWI binary mask:\n\t\tRUN -proc_dwi"; exit; fi
  if [ ! -f $dwi_SyN_warp ]; then Error "Subject $id doesn't have a non-linear transform to dwi space:\n\t\tRUN -proc_dwi"; exit; fi

#------------------------------------------------------------------------------#
Title "Running MICA POST-DWI processing (Tractography)"
  micapipe_software
  Info "Number of streamlines: $tracts"
  Info "Auto-tractograms: $autoTract"
  Info "Not erasing temporal dir: $nocleanup"
  Info "MRtrix will use $threads threads"

#	Timer
  aloita=$(date +%s)
  here=$(pwd)
  Nstep=0

# Hardcode temp dir
  tmp=${proc_sc}/tmp
  if [[ ! -d ${tmp} ]] ; then Do_cmd mkdir -p ${tmp} ; fi

# TRAP in case the script fails
  trap cleanup INT TERM

# Create Connectomes directory for the outputs
  [[ ! -d $dir_connectome ]] && Do_cmd mkdir -p $dir_connectome
  [[ ! -d $dir_QC_png ]] && Do_cmd mkdir -p $dir_QC_png
  Do_cmd cd ${tmp}

# -----------------------------------------------------------------------------------------------
# Compute fod map for non linear registration
  fod=${tmp}/${id}_space-dwi_model-CSD_map-FOD_desc-wmNorm.nii.gz
  Do_cmd mrconvert -coord 3 0 $fod_wmN $fod

# -----------------------------------------------------------------------------------------------
# Prepare the segmentatons (linear)
  parcellations=$(find ${dir_volum} -name "*.nii.gz" ! -name "*cerebellum*" ! -name "*subcortical*")
  T1_seg_cerebellum=${dir_volum}/${T1str_nat}_cerebellum.nii.gz
  T1_seg_subcortex=${dir_volum}/${T1str_nat}_subcortical.nii.gz
  dwi_cere=${proc_dwi}/${id}_dwi_cerebellum.nii.gz
  dwi_subc=${proc_dwi}/${id}_dwi_subcortical.nii.gz

  if [[ ! -f $dwi_cere ]]; then Info "Registering Cerebellar parcellation to DWI-b0 space"
      	Do_cmd antsApplyTransforms -d 3 -e 3 -i $T1_seg_cerebellum -r $dwi_b0 -n GenericLabel -t [$mat_dwi_affine,1] -o $dwi_cere -v -u int

      # Threshold cerebellar nuclei (29,30,31,32,33,34) and add 100
#        Do_cmd fslmaths $dwi_cere -uthr 28 $dwi_cere
      	Do_cmd fslmaths $dwi_cere -bin -mul 100 -add $dwi_cere $dwi_cere
  else Info "Subject ${id} has a Cerebellar segmentation in DWI space"; fi

  if [[ ! -f $dwi_subc ]]; then Info "Registering Subcortical parcellation to DWI-b0 space"
      	Do_cmd antsApplyTransforms -d 3 -e 3 -i $T1_seg_subcortex -r $dwi_b0 -n GenericLabel -t [$mat_dwi_affine,1] -o $dwi_subc -v -u int
      # Remove brain-stem (label 16)
#        Do_cmd fslmaths $dwi_subc -thr 16 -uthr 16 -binv -mul $dwi_subc $dwi_subc
  else Info "Subject ${id} has a Subcortical segmentation in DWI space"; fi

# -----------------------------------------------------------------------------------------------

# Prepare the segmentatons (non-linear registration)
dwi_cere_NL=${proc_dwi}/${id}_dwi_cerebellum_NL.nii.gz
dwi_subc_NL=${proc_dwi}/${id}_dwi_subcortical_NL.nii.gz

if [[ ! -f $dwi_cere_NL ]]; then Info "Registering Cerebellar parcellation to DWI-b0 space using NON-LINEAR transformation"
      	Do_cmd antsApplyTransforms -d 3 -r $fod -i $T1_seg_cerebellum -n GenericLabel -t $dwi_SyN_warp -t $dwi_SyN_affine -t [$mat_dwi_affine,1] -o $dwi_cere_NL -v -u int
      # Threshold cerebellar nuclei (29,30,31,32,33,34) and add 100
#        Do_cmd fslmaths $dwi_cere -uthr 28 $dwi_cere
      	Do_cmd fslmaths $dwi_cere_NL -bin -mul 100 -add $dwi_cere_NL $dwi_cere_NL
else Info "Subject ${id} has a Cerebellar segmentation registered to DWI space using a NON-LINEAR transformation"; fi

if [[ ! -f $dwi_subc_NL ]]; then Info "Registering Subcortical parcellation to DWI-b0 space using a NON-LINEAR transformation"
      	Do_cmd antsApplyTransforms -d 3 -r $fod -i $T1_seg_subcortex -n GenericLabel -t $dwi_SyN_warp -t $dwi_SyN_affine -t [$mat_dwi_affine,1] -o $dwi_subc_NL -v -u int
      # Remove brain-stem (label 16)
#        Do_cmd fslmaths $dwi_subc_NL -thr 16 -uthr 16 -binv -mul $dwi_subc_NL $dwi_subc_NL 											# COMMIT REQUIRES FULL BRAIN
else Info "Subject ${id} has a Subcortical segmentation registered to DWI space using a NON-LINEAR transformation"; fi

  if [[ -f $dwi_cere_NL ]]; then ((Nstep++)); fi
  if [[ -f $dwi_subc_NL ]]; then ((Nstep++)); fi

# -----------------------------------------------------------------------------------------------

# Generate probabilistic tracts
  Info "Building the ${tracts} streamlines connectome!!!"
  tck=${proc_sc}/DWI_tractogram_${tracts}.tck
  weights_sift2=${proc_sc}/SIFT2_${tracts}.txt
  tdi=${proc_dwi}/${id}_tdi_iFOD2-${tracts}.mif

  if [[ ! -f $tck ]]; then Info "Creating new tractogram"
	Do_cmd tckgen -nthreads $threads \
    		$fod_wmN \
    		$tck \
    		-act $dwi_5tt \
	        -crop_at_gmwmi \
		-seed_gmwmi $dwi_gmwmi \
    		-maxlength 400 \
    		-minlength 10 \
    		-angle 22.5 \
    		-backtrack \
    		-select ${tracts} \
    		-step .5 \
    		-cutoff 0.06 \
    		-algorithm iFOD2
  else Info "Subject ${id} has existing tractogram!!!"; fi

# Exit if tractography fails
  if [ ! -f $tck ]; then Error "Tractogram failed, check the logs: $(ls -Art ${dir_logs}/post-dwi_*.txt | tail -1)"; exit;
  else ((Nstep++)); fi

# json file of tractogram
#   tck_json iFOD2 0.5 22.5 0.06 400 10 seed_gmwmi $tck

# SIFT2
  if [[ ! -f $weights_sift2 ]]; then Info "Running SIFT2"
	Do_cmd tcksift2 -nthreads $threads $tck $fod_wmN $weights_sift2
  else Info "Subject ${id} has SIFT2 weights"; fi
  if [[ -f $weights_sift2 ]]; then ((Nstep++)); fi

# TDI for QC
  if [[ ! -f $tdi ]]; then
	Info "Creating a Track Density Image (tdi) of the $tracts connectome for QC"
	Do_cmd tckmap -vox 1,1,1 -dec -nthreads $threads $tck $tdi -force
  else Info "EXISTING FILE: ${tdi}" ; fi

# -----------------------------------------------------------------------------------------------
# Compute Auto-Tractography
  if [ $autoTract == "TRUE" ]; then
    	Info "Running Auto-tract"
    	autoTract_dir=$proc_dwi/auto_tract
    	[[ ! -d $autoTract_dir ]] && Do_cmd mkdir -p $autoTract_dir
    	fa_niigz=$tmp/${id}_dti_FA.nii.gz
    	Do_cmd mrconvert $fa $fa_niigz
    	echo -e "\033[38;5;118m\nCOMMAND -->  \033[38;5;122m03_auto_tracts.sh -tck $tck -outbase $autoTract_dir/${id} -mask $dwi_mask -fa $fa_niigz -tmpDir $tmp -keep_tmp  \033[0m"
    	${MICAPIPE}/functions/03_auto_tracts.sh -tck $tck -outbase $autoTract_dir/${id}_${tracts} -mask $dwi_mask -fa $fa_niigz -weights $weights -tmpDir $tmp -keep_tmp
  fi

# -----------------------------------------------------------------------------------------------
# Clean temporal directory
  if [[ $nocleanup == "FALSE" ]]; then Do_cmd rm -rf $tmp; else Info "Mica-pipe tmp directory was not erased: \n\t\t\t${tmp}"; fi
  Do_cmd cd $here

# QC notification of completition
  lopuu=$(date +%s)
  eri=$(echo "$lopuu - $aloita" | bc)
  eri=$(echo print $eri/60 | perl)

# Notification of completition
  if [ "$Nstep" -eq 4 ]; then status="COMPLETED"; else status="ERROR! Double check that tractogram & sift weights were generated!"; fi

  Title "DWI-post TRACTOGRAPHY processing ended in \033[38;5;220m $(printf "%0.3f\n" "${eri}") minutes \033[38;5;141m:
  \t\tSteps completed: $(printf "%01d" "$Nstep")/4
  \tStatus          : ${status}
  \tCheck logs      : $(ls "${dir_logs}"/proc-dwi_*.txt)"

# Print QC stamp
  echo "${id}, post_dwi, $status N=$(printf "%02d" $Nstep)/2, $(whoami), $(uname -n), $(date), $(printf "%0.3f\n" ${eri}), $PROC" >> ${out}/brain-proc.csv
  bids_variables_unset
