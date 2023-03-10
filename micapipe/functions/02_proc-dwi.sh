#!/bin/bash
#
# DWI structural processing with bash:
#
# Preprocessing workflow for diffusion MRI.
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
threads=$7
unring=$8
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

# Check inputs: DWI
  if [ "${#bids_dwis[@]}" -lt 1 ]; then Error "Subject $id doesn't have DWIs:\n\t\t TRY <ls -l ${subject_bids}/dwi/>"; exit; fi
  if [ ! -f ${T1_MNI152_InvWarp} ]; then Error "Subject $id doesn't have T1_nativepro warp to MNI152.\n\t\tRun -proc_structural"; exit; fi
  if [ ! -f ${T1nativepro} ]; then Error "Subject $id doesn't have T1_nativepro.\n\t\tRun -proc_structural"; exit; fi
  if [ ! -f ${bids_inv1} ]; then Error "Subject $id doesn't have inversion 1 T1map.\n\t\tRun -proc_structural"; exit; fi
  if [ ! -f ${bids_T1map} ]; then Error "Subject $id doesn't have MP2RAGE T1map.\n\t\tRun -proc_structural"; exit; fi
  if [ ! -f ${T15ttgen} ]; then Error "Subject $id doesn't have a 5tt volume in nativepro space.\n\t\tRun -proc_structural"; exit; fi

# CHECK if PhaseEncodingDirection and TotalReadoutTime exist
  for i in ${bids_dwis[*]}; do
  	json=$(echo "${i}" | awk -F ".nii" '{print $1 ".json"}')
  	ped=$(grep PhaseEncodingDirection\": "${json}" | awk -F "\"" '{print $4}')
  	trt=$(grep TotalReadoutTime "${json}" | awk -F " " '{print $2}')
  	if [[ -z "$ped" ]]; then Error "PhaseEncodingDirection is missing in $json"; exit; fi
  	if [[ -z "$trt" ]]; then Error "TotalReadoutTime is missing in $json"; exit; fi
  done

#------------------------------------------------------------------------------#
  Title "Running MICA Diffusion Weighted Imaging processing"
  micapipe_software
  bids_print.variables-dwi
  Info "Not erasing temporal dir: $nocleanup"
  Info "ANTs and MRtrix will use $threads threads"
  Info "Perform unringing: $unring"

#	Timer
  aloita=$(date +%s)
  Nsteps=0

# Hardcode temp dir
  tmp=${proc_dwi}/tmp
  if [ ! -d ${tmp} ]; then Do_cmd mkdir -p ${tmp} ;  fi

# TRAP in case the script fails
  trap cleanup INT TERM

  Do_cmd cd $tmp

#------------------------------------------------------------------------------#
# DWI processing
# Image denoising must be performed as the first step of the image-processing pipeline.
# Interpolation or smoothing in other processing steps, such as motion and distortion correction,
# may alter the noise characteristics and thus violate the assumptions upon which MP-PCA is based.

  dwi_cat=${tmp}/dwi_concatenate.mif
  dwi_dns=${tmp}/dwi_concatenate_denoised.mif
  dwi_corr_biased=${proc_dwi}/${id}_dwi_corr_biased.mif
  dwi_corr=${proc_dwi}/${id}_dwi_corr.mif
  dwi_res=${proc_dwi}/${id}_dwi_residuals.mif
#  b0_refacq=$(echo ${bids_dwis[0]} | awk -F 'acq-' '{print $2}'| sed 's:_dwi.nii.gz::g')

  if [[ ! -f $dwi_dns ]]; then
	Info "DWI denoise and concatenation"
      # Concatenate shells -if only one shell then just convert to mif and rename.
      	for dwi in ${bids_dwis[@]}; do
	    	dwi_nom=$(echo "${dwi##*/}" | awk -F ".nii" '{print $1}') 									   # strips path & ext
	    	bids_dwi_str=$(echo "$dwi" | awk -F . '{print $1}') 										   # strips ext only
            	Do_cmd mrconvert $dwi -json_import ${bids_dwi_str}.json -fslgrad ${bids_dwi_str}.bvec ${bids_dwi_str}.bval ${tmp}/${dwi_nom}.mif   # convert to mif with added metadata
#	    	Do_cmd dwiextract ${tmp}/${dwi_nom}.mif "${tmp}/${dwi_nom}_b0.mif" -bzero 							   # no need to extract b0s if not running registration
#	    	Do_cmd mrmath "${tmp}/${dwi_nom}_b0.mif" mean "${tmp}/${dwi_nom}_b0.nii.gz" -axis 3 -nthreads $threads
        done

       # Concatenate shells and convert to mif.
        Info "Concatenatenating shells"
        if [ "${#bids_dwis[@]}" -eq 1 ]; then
                cp $tmp/*.mif "$dwi_cat"
        else
		Do_cmd mrcat $tmp/*.mif "$dwi_cat" -nthreads $threads
        fi

      # Denoise DWI and calculate residuals
      	Do_cmd dwidenoise $dwi_cat $dwi_dns -nthreads $threads
      	Do_cmd mrcalc $dwi_cat $dwi_dns -subtract $dwi_res -nthreads $threads

      # Unringing (NOT TESTED)
	if [ $unring == "TRUE" ]; then
		Info "Performing unringing using mrtrix3s mrdegibbs..."
      		Do_cmd mrdegibbs $dwi_dns $dwi_dns -nthreads $threads
	fi

      	if [[ -f ${dwi_dns} ]]; then ((Nsteps++)); fi
  else
      	Info "Skipping denoising & concatenation because subject ${id} has a denoised DWI file"; ((Nsteps++))
  fi

#------------------------------------------------------------------------------#
# dwifslpreproc and TOPUP preparations
  if [[ ! -f $dwi_corr ]]; then
      	Info "--------------- DWI distortion & bias field correction ---------------"

      	ReadoutTime=$(mrinfo "$dwi_dns" -property TotalReadoutTime)
      	pe_dir=$(mrinfo "$dwi_dns" -property PhaseEncodingDirection)
      	shells=($(mrinfo "$dwi_dns" -shell_bvalues))

      # Exclude shells with 0 value (with threshold)
      	for i in "${!shells[@]}"; do if [ ${shells[i]%.*} -le 15 ]; then unset 'shells[i]'; fi; done

      # Remove slices to make an even number of slices in all directions (requisite for dwi_preproc-TOPUP).
      	dwi_4proc=${tmp}/dwi_dns_even.mif
      	dim=$(mrinfo "$dwi_dns" -size)
      	dimNew=($(echo $dim | awk '{for(i=1;i<=NF;i++){$i=$i-($i%2);print $i-1}}'))
      	mrconvert $dwi_dns $dwi_4proc -coord 0 0:${dimNew[0]} -coord 1 0:${dimNew[1]} -coord 2 0:${dimNew[2]} -coord 3 0:end -force

      # Get the mean b-zero (un-corrected)
      	dwiextract -nthreads $threads $dwi_dns - -bzero | mrmath - mean ${tmp}/b0_meanMainPhase.mif -axis 3

      # Processing the reverse encoding b0
      	if [ -f $dwi_reverse ]; then
            	b0_pair_tmp=${tmp}/b0_pair_tmp.mif
            	b0_pair=${tmp}/b0_pair.mif
            	dwi_reverse_str=$(echo $dwi_reverse | awk -F . '{print $1}')
	    	rpe_dim=$(mrinfo "$dwi_reverse" -ndim)

	      # Mean reverse phase b0
            	if [[ -f "${dwi_reverse_str}.bvec" ]] && [[ -f "${dwi_reverse_str}.bval" ]]; then
       			Do_cmd mrconvert $dwi_reverse -json_import ${dwi_reverse_str}.json -fslgrad ${dwi_reverse_str}.bvec ${dwi_reverse_str}.bval ${tmp}/b0_ReversePhase.mif
	    	else
			Do_cmd mrconvert "$dwi_reverse" -json_import "$dwi_reverse_str.json" "${tmp}/b0_ReversePhase.mif"
            	fi

	      # Concatenate the pe and rpe b0s
		dwiextract ${tmp}/b0_ReversePhase.mif - -bzero | mrmath - mean ${tmp}/b0_meanReversePhase.mif -axis 3 -nthreads $threads
            	Do_cmd mrcat ${tmp}/b0_meanMainPhase.mif ${tmp}/b0_meanReversePhase.mif $b0_pair_tmp -nthreads $threads

              # Remove slices to make an even number of slices in all directions (requisite for dwi_preproc-TOPUP).
            	dim=$(mrinfo "$b0_pair_tmp" -size)
            	dimNew=($(echo $dim | awk '{for(i=1;i<=NF;i++){$i=$i-($i%2);print $i-1}}'))
            	mrconvert $b0_pair_tmp $b0_pair -coord 0 0:${dimNew[0]} -coord 1 0:${dimNew[1]} -coord 2 0:${dimNew[2]} -coord 3 0:end -force
            	opt="-rpe_pair -align_seepi -se_epi ${b0_pair}"
      	else
	    	Info "Reverse phase encoding image was not found it will be omitted"
            	opt='-rpe_none'
      	fi

      	Info "dwifslpreproc parameters:"
      	Note "Shell values        :" " ${shells[*]%.*} "
      	Note "DWI main dimensions :" " $(mrinfo $dwi_4proc -size) "
      	if [ -f $dwi_reverse ]; then Note "DWI rpe dimensions  :" " $(mrinfo $b0_pair -size) " ; fi
      	Note "pe_dir              :" " $pe_dir "
      	Note "Readout Time        :" " $ReadoutTime "

      # Preprocess each shell
      # DWIs all acquired with a single fixed phase encoding; but additionally a
      # pair of b=0 images with reversed phase encoding to estimate the inhomogeneity field:
      	echo -e "COMMAND --> dwifslpreproc $dwi_4proc $dwi_corr $opt -pe_dir $pe_dir -readout_time $ReadoutTime -eddy_options \" --data_is_shelled --slm=linear\" -nthreads $threads -nocleanup -scratch $tmp -force"
      	dwifslpreproc $dwi_4proc $dwi_corr_biased $opt -pe_dir $pe_dir -readout_time $ReadoutTime -eddy_options " --data_is_shelled --slm=linear" -nthreads $threads -nocleanup -scratch $tmp -force
      # Step QC
      	if [[ ! -f ${dwi_corr_biased} ]]; then Error "dwifslpreproc failed, check the logs"; exit;
      	else
              # eddy_quad Quality Check
	  	Do_cmd cd $tmp/dwifslpreproc*
          	Do_cmd eddy_quad dwi_post_eddy -idx eddy_indices.txt -par eddy_config.txt -m eddy_mask.nii -b bvals -o ${dir_QC}/eddy_QC
	  	Do_cmd cd $tmp
              # Copy eddy parameters
          	eddy_DIR=${proc_dwi}/eddy
          	if [ ! -d ${eddy_DIR} ]; then Do_cmd mkdir ${eddy_DIR}; fi
          	Do_cmd cp -rf ${tmp}/dwifslpreproc*/*eddy.eddy* ${eddy_DIR}
	  	Do_cmd chmod 770 -R ${eddy_DIR}/*
      	fi

      # Bias field correction DWI
	# note: mtnormalise (later) also corrects for bias fields.
	# 	Keeping dwibiascorrect here may improve registration & mask estimation (if strong bias fields are present)
	# 	Keeping/removing dwibiascorrect at this stage does not impact mtnormalise later on
      	Do_cmd dwibiascorrect ants $dwi_corr_biased $dwi_corr -force -nthreads $threads -scratch $tmp -nocleanup
      # Step QC
      	if [[ ! -f ${dwi_corr} ]]; then Error "bias correction failed, check the logs"; exit;
      	else
		Do_cmd rm $dwi_corr_biased; ((Nsteps++))
      	fi
  else
      	Info "Skipping distortion & bias field correction because subject ${id} has a fully processed DWI..."; ((Nsteps++))
  fi

#------------------------------------------------------------------------------#
# Registration of corrected DWI-b0 to T1nativepro
  dwi_mask=${proc_dwi}/${id}_dwi_mask.nii.gz
  dwi_b0=${proc_dwi}/${id}_dwi_b0.nii.gz 					# This should be a NIFTI for compatibility with ANTS
  T1nativepro_in_dwi=${proc_dwi}/${id}_t1w_in_dwi.nii.gz
  str_dwi_affine=${dir_warp}/${id}_dwi_to_nativepro_
  mat_dwi_affine=${str_dwi_affine}0GenericAffine.mat

  if [[ ! -f $T1nativepro_in_dwi ]]; then
      	Info "Linear registration between DWI-b0 and T1nativepro"
      # Corrected DWI-b0s mean for registration
      	dwiextract -force -nthreads $threads $dwi_corr - -bzero | mrmath - mean $dwi_b0 -axis 3 -force

      # Register DWI-b0 mean corrected to T1nativepro & qT1
      	Do_cmd antsRegistrationSyN.sh -d 3 -f $T1nativepro_brain -m $dwi_b0 -o $str_dwi_affine -t a -n $threads -p d
      # Apply inverse transformation T1nativepro to DWI-b0 space
      	Do_cmd antsApplyTransforms -d 3 -i $T1nativepro -r $dwi_b0 -t [$mat_dwi_affine,1] -o $T1nativepro_in_dwi -v -u int
      	if [[ -f ${T1nativepro_in_dwi} ]] ; then ((Nsteps++)); fi

       #------------------------------------------------------------------------------#

      	if [[ ! -f ${dwi_mask} ]]; then
      		Info "Creating DWI binary mask of processed volumes"
      	      # Create a binary mask of the DWI
      		Do_cmd antsApplyTransforms -d 3 -i $MNI152_mask \
              		-r ${dwi_b0} \
              		-n GenericLabel -t [$mat_dwi_affine,1] -t [${T1_MNI152_affine},1] -t ${T1_MNI152_InvWarp} \
              		-o ${tmp}/dwi_mask.nii.gz -v
      		Do_cmd maskfilter ${tmp}/dwi_mask.nii.gz erode -npass 1 $dwi_mask
      		if [[ -f ${dwi_mask} ]]; then ((Nsteps++)); fi
      	else
		Info "Subject ${id} has a DWI mask"; ((Nsteps++))
	fi
  else
      	Info "Subject ${id} has an affine transformation from T1w to DWI-b0 space"; Nsteps=$((Nsteps + 2))
  fi

#------------------------------------------------------------------------------#
# Registration of qT1 to T1w nativepro
  qT1_in_t1w=${proc_struct}/${id}_qt1_in_nativepro.nii.gz
  str_qt1_to_t1w_affine=${dir_warp}/${id}_qt1_to_t1w_
  mat_qt1_to_t1w_affine=${str_qt1_to_t1w_affine}0GenericAffine.mat

  if [[ ! -f $qT1_in_t1w ]]; then
      	Info "RIGID registration qT1 to T1w nativepro space"
      # Compute transform for rigid registration of qT1 to T1w nativepro
      	Do_cmd antsRegistrationSyN.sh -d 3 -f $T1nativepro -m $bids_inv1 -o $str_qt1_to_t1w_affine -t r -n $threads -p d

      # Apply transformation qT1 to T1w nativepro
      	Do_cmd antsApplyTransforms -d 3 -i $bids_T1map -r $T1nativepro -t $mat_qt1_to_t1w_affine -o $qT1_in_t1w -v -u int
      	if [[ -f ${qT1_in_t1w} ]]  ; then ((Nsteps++)); fi
  else
      	Info "Subject ${id} has an affine transformation from qT1 to nativepro"; Nsteps=$((Nsteps + 2))
  fi

#------------------------------------------------------------------------------#
# Get some basic metrics.
  dwi_dti=$proc_dwi/${id}_dti.mif
  dwi_FA=$proc_dwi/${id}_dti_FA.mif
  dwi_ADC=$proc_dwi/${id}_dti_ADC.mif
  dwi_AD=$proc_dwi/${id}_dti_AD.mif
  dwi_RD=$proc_dwi/${id}_dti_RD.mif
  if [[ ! -f $dwi_RD ]]; then
      	Info "Calculating basic DTI metrics"
      	dwi2tensor -mask $dwi_mask -nthreads $threads $dwi_corr $dwi_dti
      	tensor2metric -nthreads $threads -fa ${dwi_FA} -adc ${dwi_ADC} -ad ${dwi_AD} -rd ${dwi_RD} $dwi_dti
      # Step QC
      	if [[ -f ${dwi_RD} ]]; then ((Nsteps++)); fi
  else
      	Info "Subject ${id} has DWI tensor metrics"; ((Nsteps++))
  fi

#------------------------------------------------------------------------------#
# Response function and Fiber Orientation Distribution
  fod_wmN=$proc_dwi/${id}_wm_fod_norm.mif
  fod_gmN=$proc_dwi/${id}_gm_fod_norm.mif
  fod_csfN=$proc_dwi/${id}_csf_fod_norm.mif
  if [[ ! -f $fod_wmN ]]; then
      	Info "Calculating Multi-Shell Multi-Tissue, Response function and Fiber Orientation Distribution"
       	rf=dhollander
      # Response function
       	rf_wm=${tmp}/${id}_response_wm_${rf}.txt
        rf_gm=${tmp}/${id}_response_gm_${rf}.txt
        rf_csf=${tmp}/${id}_response_csf_${rf}.txt
      # Fiber Orientation Distriution
        fod_wm=${tmp}/${id}_wm_fod.mif
        fod_gm=${tmp}/${id}_gm_fod.mif
        fod_csf=${tmp}/${id}_csf_fod.mif

        Do_cmd dwi2response $rf -nthreads $threads $dwi_corr $rf_wm $rf_gm $rf_csf -mask $dwi_mask
        Do_cmd dwi2fod -nthreads $threads msmt_csd $dwi_corr \
                $rf_wm $fod_wm \
                $rf_gm $fod_gm \
                $rf_csf $fod_csf \
                -mask $dwi_mask
      	if [ "${#shells[@]}" -ge 2 ]; then
            	Do_cmd mtnormalise $fod_wm $fod_wmN $fod_gm $fod_gmN $fod_csf $fod_csfN -nthreads $threads -mask $dwi_mask
	    	Do_cmd mrinfo "$fod_wmN" -json_all "${fod_wmN/mif/json}"
            	Do_cmd mrinfo "$fod_gmN" -json_all "${fod_gmN/mif/json}"
            	Do_cmd mrinfo "$fod_csfN" -json_all "${fod_csfN/mif/json}"
      	else
            	Do_cmd mtnormalise -nthreads $threads -mask $dwi_mask $fod_wm $fod_wmN
	    	Do_cmd mrinfo "$fod_wmN" -json_all "${fod_wmN/mif/json}"
      	fi
      	if [[ -f ${fod_wmN} ]]; then ((Nsteps++)); fi
  else
      	Info "Subject ${id} has Fiber Orientation Distribution files"; ((Nsteps++))
  fi

#------------------------------------------------------------------------------#
# Non-linear registration between DWI space and T1w nativepro & qT1 map (Solution to poor linear registrations to dwi space)
  dwi_SyN_str=${dir_warp}/${id}_space-dwi_from-dwi_to-dwi_mode-image_desc-SyN_
  dwi_SyN_warp=${dwi_SyN_str}1Warp.nii.gz
  dwi_SyN_Invwarp=${dwi_SyN_str}1InverseWarp.nii.gz
  dwi_SyN_affine=${dwi_SyN_str}0GenericAffine.mat
  dwi_5tt=${proc_dwi}/${id}_dwi_5tt.nii.gz
  dwi_in_T1nativepro=${proc_struct}/${id}_dwi_in_nativepro.nii.gz                                # Only for QC
  T1nativepro_in_dwi_brain=${proc_dwi}/${id}_t1w_in_dwi_brain.nii.gz
  T1nativepro_in_dwi_NL=${proc_dwi}/${id}_t1w_in_dwi_NL.nii.gz
  qT1_in_dwi_NL=${proc_dwi}/${id}_qt1_in_dwi_NL.nii.gz
  fod=${tmp}/${id}_space-dwi_model-CSD_map-FOD_desc-wmNorm.nii.gz

  if [[ ! -f $T1nativepro_in_dwi_NL ]] && [[ ! -f $qT1_in_dwi_NL ]] && [[ ! -f $dwi_5tt ]]  ; then
    	Info "Non-linear registration from Native space to DWI"
    	Do_cmd fslmaths $T1nativepro_in_dwi -mul $dwi_mask $T1nativepro_in_dwi_brain
    	Do_cmd mrconvert -coord 3 0 $fod_wmN $fod

    	if [[ ! -f $dwi_SyN_warp ]] ; then
   	      # compute native to dwi non-linear transform
    		Do_cmd antsRegistrationSyN.sh -d 3 -m $T1nativepro_in_dwi_brain -f $fod -o $dwi_SyN_str -t s -n $threads
  		if [[ -f $dwi_SyN_warp ]]; then ((Nsteps++)); fi
    	else
    		Info "Subject ${id} has a non linear transform from native to diffusion space"; ((Nsteps++))
    	fi

    	Info "Registering the T1map, T1w-nativepro and 5TT to DWI-b0 space, and DWI-b0 to native space"
      # Apply transformation DWI-b0 space to T1nativepro
    	Do_cmd antsApplyTransforms -d 3 -r $T1nativepro_brain -i $dwi_b0 -t $mat_dwi_affine -t [$dwi_SyN_affine,1] -t $dwi_SyN_Invwarp -o $dwi_in_T1nativepro -v -u int
      # Apply transformation T1nativepro to DWI space
    	Do_cmd antsApplyTransforms -d 3 -r $fod -i $T1nativepro -t $dwi_SyN_warp -t $dwi_SyN_affine -t [$mat_dwi_affine,1] -o $T1nativepro_in_dwi_NL -v -u int
      # Apply transformation T1map to DWI space
    	Do_cmd antsApplyTransforms -d 3 -r $fod -i $bids_T1map -t $dwi_SyN_warp -t $dwi_SyN_affine -t [$mat_dwi_affine,1] -t $mat_qt1_to_t1w_affine  -o $qT1_in_dwi_NL -v -u int
      # Apply transformation 5TT to DWI space
    	Do_cmd antsApplyTransforms -d 3 -r $fod -i $T15ttgen -t $dwi_SyN_warp -t $dwi_SyN_affine -t [$mat_dwi_affine,1] -o $dwi_5tt -v -e 3 -n linear
    	if [[ -f $dwi_5tt ]]; then ((Nsteps++)); fi
  else
    	Info "Subject ${id} has a non-linear registration from T1w_dwi-space to DWI"; Nsteps=$((Nsteps + 2))
  fi

#------------------------------------------------------------------------------#
# Gray matter White matter interface mask
  dwi_gmwmi=${proc_dwi}/${id}_space-dwi_desc-gmwmi-mask.mif
  if [[ ! -f $dwi_gmwmi ]]; then
      	Info "Calculating Gray matter White matter interface mask"
      	5tt2gmwmi $dwi_5tt $dwi_gmwmi; ((Nsteps++))
  else
      	Info "Subject ${id} has Gray matter White matter interface mask"; ((Nsteps++))
  fi

#------------------------------------------------------------------------------#
# QC of the tractography
  tracts=1M
  tdi_1M=${proc_dwi}/${id}_tdi_iFOD1-1M.mif
  tckjson=${proc_dwi}/${id}_tdi_iFOD1-1M.json

  if [[ ! -f $tdi_1M ]]; then
      	Info "Creating a track density image for quality check"
      	tck_1M=${tmp}/${id}_tdi_iFOD1-1M.tck
      	Do_cmd tckgen -nthreads $threads \
          	$fod_wmN \
          	$tck_1M \
          	-act $dwi_5tt \
          	-crop_at_gmwmi \
          	-backtrack \
          	-seed_gmwmi $dwi_gmwmi \
          	-maxlength 400 \
	  	-minlength 10 \
          	-angle 22.5 \
          	-power 1.0 \
          	-select ${tracts} \
          	-step .5 \
          	-cutoff 0.06 \
          	-algorithm iFOD1

      	Do_cmd tckmap -vox 1,1,1 -dec -nthreads $threads $tck_1M $tdi_1M
      # Step QC
      	if [[ -f ${tdi_1M} ]]; then ((Nsteps++)); fi
      tck_json iFOD1 0.5 22.5 0.06 400 10 seed_gmwmi $tck_1M
  else
      	Info "Subject ${id} has a Tract Density Image for QC 1M streamlines"; ((Nsteps++))
  fi

# -----------------------------------------------------------------------------------------------
# Clean temporal directory
  Do_cmd cd $here
  if [[ $nocleanup == "FALSE" ]]; then Do_cmd rm -rf $tmp; else Info "Mica-pipe tmp directory was not erased: \n\t\t\t${tmp}"; fi

# QC notification of completition
  lopuu=$(date +%s)
  eri=$(echo "$lopuu - $aloita" | bc)
  eri=$(echo print "$eri"/60 | perl)

# Notification of completition
  if [ "$Nsteps" -eq 11 ]; then status="COMPLETED"; else status="ERROR DWI is missing a processing step: "; fi

  Title "DWI processing ended in \033[38;5;220m $(printf  "%0.3f\n"  "${eri}") minutes \033[38;5;141m:
  \tSteps completed : $(printf "%02d" "$Nsteps")/11
  \tStatus          : ${status}
  \tCheck logs      : $(ls "${dir_logs}"/proc-dwi_*.txt)"

# Print QC stamp
  echo "${id}, proc_dwi, "$status", $(printf "%02d" $Nsteps)/10, $(whoami), $(uname -n), $(date), $(printf "%0.3f\n" ${eri}), $PROC" >> ${out}/brain-proc.csv
  bids_variables_unset


## scratch
#------------------------------------------------------------------------------#
# Non-linear registration between masked T1w and b0-dwi
# str_dwiT1_2_b0=${dir_warp}/${id}_dwiT1_to_b0_
# dwiT1_2_b0_warp=${str_dwiT1_2_b0}1Warp.nii.gz
#
# if [[ ! -f $dwiT1_2_b0_warp ]]; then
#     dwi_T1_masked=$tmp/dwi_T1_masked.nii.gz
#     dwi_b0_masked=$tmp/dwi_b0_masked.nii.gz
#     tmp_mask=${tmp}/dwi_mask_eroded.nii.gz
#
#     Do_cmd maskfilter $dwi_mask erode -npass 5 $tmp_mask
#     Do_cmd ImageMath 3 dwi_b0_scaled.nii.gz RescaleImage $dwi_b0 0 100
#     Do_cmd fslmaths dwi_b0_scaled.nii.gz -mul -1 -add 99 -mul $tmp_mask dwi_b0_inv.nii.gz
#
#     Do_cmd fslmaths $T1nativepro_in_dwi -mul $tmp_mask $dwi_T1_masked
#     Do_cmd ImageMath 3 dwi_b0_matched.nii.gz HistogramMatch dwi_b0_inv.nii.gz $dwi_T1_masked
#
#     Do_cmd antsRegistrationSyN.sh -d 3 -x $tmp_mask -m $dwi_T1_masked -f dwi_b0_matched.nii.gz -o $str_dwiT1_2_b0 -t bo -n $threads -p d
#     Do_cmd antsApplyTransforms -d 3 -e 3 -i $T1nativepro_in_dwi -r $dwi_b0 -n linear -t $dwiT1_2_b0_warp -o $T1nativepro_in_dwi -v
#     Do_cmd antsApplyTransforms -d 3 -e 3 -i $dwi_5tt -r $dwi_b0 -n linear -t $dwiT1_2_b0_warp -o $dwi_5tt -v
# else
#       Info "Subject ${id} has a non-linear registration from dwi-T1w to dwi-b0"; ((Nsteps++))
# fi

#------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
      # Rigid registration between shells (soluton to deal with misaligned B0 in HC030)
#       n=$((${#bids_dwis[*]} - 1))
#       if [[ ${#bids_dwis[*]} -gt 1 ]]; then
#               b0_ref=${tmp}/$(echo "${bids_dwis[0]}" | awk -F "dwi/" '{print $2}' | awk -F ".nii" '{print $1}')_b0.nii.gz         # selects 1st b0 image as ref
#               for ((i=1; i<=$n; i++)); do                                                                                         # counts up to n from 1 (skips reference image)
#                       dwi_nom=$(echo ${bids_dwis[i]} | awk -F "dwi/" '{print $2}' | awk -F ".nii" '{print $1}')                   # strips path & ext 
#                       bids_dwi_str=$(echo ${bids_dwis[i]} | awk -F . '{print $1}')                                                # full path without extension
#                       b0_acq=$(echo ${bids_dwis[i]} | awk -F 'acq-' '{print $2}'| sed 's:_dwi.nii.gz::g')                         # gets b value & dirs tag
#                       b0_nom="${tmp}/$(echo ${bids_dwis[i]} | awk -F "dwi/" '{print $2}' | awk -F ".nii" '{print $1}')_b0.nii.gz" # replaces path & ext
#                       b0_run="acq-${b0_acq}"                                                                                      # builds new acquisition tag
#                       b0mat_str="${dir_warp}/${id}_from-${b0_acq}_to-${b0_refacq}_mode-image_desc-rigid_"
#                       b0mat="${b0mat_str}0GenericAffine.mat"
#                       b0run_2_b0ref="${tmp}/${id}_from-${b0_acq}_to-${b0_refacq}.nii.gz"

#                       Info "Registering ${b0_acq} to ${b0_refacq}"
#                       Do_cmd antsRegistrationSyN.sh -d 3 -m "$b0_nom" -f "$b0_ref"  -o "$b0mat_str" -t r -n "$threads" -p d       # compute transforms on b0 image
#                       mrconvert ${tmp}/${dwi_nom}.mif ${tmp}/${dwi_nom}.nii.gz
#                       Do_cmd antsApplyTransforms -d 3 -e 3 -i "${tmp}/${dwi_nom}.nii.gz" -r "$b0_ref" -t "$b0mat" -o "${tmp}/${dwi_nom}_in-${b0_refacq}.nii.gz" -v -u int   # regist$
#                       Do_cmd mrconvert ${tmp}/${dwi_nom}_in-${b0_refacq}.nii.gz \
#                                        -json_import ${bids_dwi_str}.json \
#                                        -fslgrad ${bids_dwi_str}.bvec ${bids_dwi_str}.bval \
#                                         ${tmp}/${dwi_nom}_Ralign.mif -force -quiet                                                # convert registered images back to mif for concat$
#               done
#       fi

      # Concatenate shells and convert to mif.
#       Info "Concatenatenating shells"
#       dwi_0=$(echo ${bids_dwis[0]} | awk -F "dwi/" '{print $2}' | awk -F ".nii" '{print $1}')         # 1st dwi file, stripped of path & ext
#       if [ "${#bids_dwis[@]}" -eq 1 ]; then
#               cp ${tmp}/${dwi_0}.mif $dwi_cat
#       else
#               Do_cmd mrcat ${tmp}/${dwi_0}.mif "${tmp}/*_Ralign.mif" $dwi_cat -nthreads $threads
#       fi
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
