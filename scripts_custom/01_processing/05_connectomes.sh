#!/bin/bash
#
# Generates structural connectivity matrices using the following parameters:
#
# 	(1) Parcellations: ("aparc-a2009s" "aparc" "economo" "glasser-360" "vosdewael-100" "vosdewael-200" "vosdewael-300" "vosdewael-400" "schaefer-100" ...
# 			    "schaefer-200" "schaefer-300" "schaefer-400" "schaefer-500" "schaefer-600" "schaefer-700" "schaefer-800" "schaefer-900" "schaefer-1000")
#
# 	(2) Edge weights: ("nos" "los" "fa" "qt1" "rd" "ad" "adc" "icvf-noddi" "od-noddi" "icvf-commit"  "COMMIT"  "SIFT2" )
#
#
# Inputs:
#       $1 : subject (i.e. 01, 02, ... , 50)
#       $2 : function id (DEPRECATED)
#
#  NOTE: This function was cleaned up before release but not debugged. If problems arise, check variable names related to the mutiple filtering techniques.
#
# 2021 Mark C Nelson, McConnell Brain Imaging Centre, MNI, McGill
#------------------------------------------------------------------------------------------------------------------------------------

# options & setup
  source /PATH/TO/YOUR/init.sh
  id="HC0$1"
  subject=sub-${id}
  tracts=5M
  threads=6
  SES=ses-01
  model=StickZeppelinBall
  STAT=median
  Nparc=0

# Subject dirs
  subject_dir=${OUT_DIR}/${subject}/${SES}
  proc_dwi=${subject_dir}/proc_dwi
  proc_commit=${proc_dwi}/commit
  proc_noddi=${proc_dwi}/noddi
  proc_sc=${proc_dwi}/sc
  tmp=${proc_sc}/tmp
  dir_volum=${subject_dir}/proc_struct/volumetric
  dir_warp=${subject_dir}/xfms
  util_lut=${MICAPIPE}/parcellations/lut
  dir_connectome=${proc_sc}/connectomes
  dir_assignments=${proc_sc}/nodeassingments


# Details for the base tractogram which is used to compute NoS and SIFT2
  filtername=filt-gm                                                                    # filtered by removing streamlines not connecting to gray matter rois
  tractogram=${proc_sc}/DWI_tractogram_${tracts}_${filtername}.tck
#  dir_weights=${proc_sc}/weights/${filtername} 					# Options if comparison of filtering methods is desired
#  dir_connectome=${proc_sc}/connectomes/${filtername}
#  dir_assignments=${proc_sc}/nodeassingments/${filtername}
  weights_sift2=${proc_sc}/SIFT2_${tracts}_${filtername}.txt


# transforms
  mat_dwi_affine=${dir_warp}/${id}_dwi_to_nativepro_0GenericAffine.mat
  dwi_SyN_str=${dir_warp}/${id}_space-dwi_from-dwi_to-dwi_mode-image_desc-SyN_
  dwi_SyN_warp=${dwi_SyN_str}1Warp.nii.gz
  dwi_SyN_affine=${dwi_SyN_str}0GenericAffine.mat

# Files
  lut_sc="${util_lut}/lut_subcortical-cerebellum_mics.csv"
  fod=${tmp}/${id}_space-dwi_model-CSD_map-FOD_desc-wmNorm.nii.gz
#  weights_sift2=${proc_sc}/SIFT2_${tracts}.txt 					# sift2 weights computed on unfiltered tractogram (not used here)
  weights_commit=${proc_commit}/Results_${model}/streamline_weights.txt
  dti_fa=${proc_dwi}/${id}_dti_FA.mif
  dti_ad=${proc_dwi}/${id}_dti_AD.mif
  dti_rd=${proc_dwi}/${id}_dti_RD.mif
  dti_adc=${proc_dwi}/${id}_dti_ADC.mif
  qT1_in_dwi_NL=${proc_dwi}/${id}_qt1_in_dwi_NL.nii.gz
  noddi_icvfMap=${proc_noddi}/AMICO/NODDI/FIT_ICVF.nii
  noddi_odMap=${proc_noddi}/AMICO/NODDI/FIT_OD.nii
  map_icvf_commit=${proc_commit}/Results_${model}/compartment_IC.nii.gz
  dwi_cere_NL=${proc_dwi}/${id}_dwi_cerebellum_NL.nii.gz
  dwi_subc_NL=${proc_dwi}/${id}_dwi_subcortical_NL.nii.gz


# Check inputs
  if [ ! -f $lut_sc ]; then echo "Can't find the subcortical-cerebellar LUT!"; exit; fi
  if [ ! -f $fod ]; then echo "Can't find the wm fod for registration of parcellations!"; exit; fi
  if [ ! -f $tractogram ]; then echo "Subject $id doesn't have a tractogram!"; exit; fi
  if [ ! -f $weights_commit ]; then echo "Subject $id doesn't have streamline weights from commit!"; exit; fi
  if [ ! -f $weights_sift2 ]; then echo "Subject $id doesn't have sift2 weights!"; exit; fi
  if [ ! -f $map_icvf_commit ]; then echo "Subject $id doesn't have an icvf map from commit!"; exit; fi
  if [ ! -f $noddi_icvfMap ]; then echo "Subject $id doesn't have an icvf map from noddi!"; exit; fi
  if [ ! -f $noddi_odMap ]; then echo "Subject $id doesn't have an od map from noddi!"; exit; fi
  if [ ! -f $qT1_in_dwi_NL ]; then echo "Subject $id doesn't have a qt1 map in dwi space!"; exit; fi
  if [ ! -f $dti_adc ]; then echo "Subject $id doesn't have an adc map!"; exit; fi
  if [ ! -f $dti_rd ]; then echo "Subject $id doesn't have an rd map!"; exit; fi
  if [ ! -f $dti_ad ]; then echo "Subject $id doesn't have an ad map!"; exit; fi
  if [ ! -f $dti_fa ]; then echo "Subject $id doesn't have an fa map!"; exit; fi
  if [ ! -f $mat_dwi_affine ]; then echo "Subject $id doesn't have a dwi to native affine transform!"; exit; fi
  if [ ! -f $dwi_SyN_warp ]; then echo "Subject $id doesn't have a non-linear warp in dwi space!"; exit; fi
  if [ ! -f $dwi_SyN_affine ]; then echo "Subject $id doesn't have a non-linear affine to dwi space!"; exit; fi
  if [ ! -f $dwi_cere_NL ]; then echo "Subject $id doesn't have cerebellar segmentation registered to dwi space!"; exit; fi
  if [ ! -f $dwi_subc_NL ]; then echo "Subject $id doesn't have a subcortical segmentation registered to dwi space!"; exit; fi



# ------------------------- ADDITIONAL TRACTOGRAM FILTERING WITH COMMIT WEIGHTS ------------------------------ #
# Removes all streamlines with a COMMIT weight of 0 within machine precision
# This filtered tractogram will be used to compute COMMIT and all tractometry networks


  filtername_commit=filt-commit
  threshold=0.000000000001 													# chosen to eliminate streamlines = 0 at machine precision
  tckOUT=${proc_sc}/DWI_tractogram_${tracts}_${filtername_commit}.tck
  weights_commitOUT=${proc_commit}/Results_${model}/streamline_weights_${filtername_commit}.txt

# Use COMMIT weights to filter tractogram
  if [ ! -f $tckOUT ]; then echo "Removing Streamlines with COMMIT weight < $threshold"
    tckedit -minweight $threshold -tck_weights_in $weights_commit -tck_weights_out $weights_commitOUT $tractogram $tckOUT
  else echo " A COMMIT-Filtered Tractogram already exists for this subject  "; fi

# Define directories & necessary files
  tractogram_commit=${proc_sc}/DWI_tractogram_${tracts}_${filtername_commit}.tck
  weights_commit=${proc_commit}/Results_${model}/streamline_weights_${filtername_commit}.txt 						# Pre-filtering weights are not used again
  dir_weights=${proc_sc}/weights/${filtername_commit}
#  dir_connectome_commit=${proc_sc}/connectomes/${filtername_commit} 									# Options if comparison of filtering methods is desired
#  dir_assignments_commit=${proc_sc}/nodeassingments/${filtername_commit}
#  weights_sift2_commit=${proc_sc}/SIFT2_${tracts}_${filtername_commit}.txt 								# Not used here

# -------------------------------------------------------------------------------------------------------------------- #


# Ensure SIFT2 has been run for this version of the tractogram
  if [[ ! -f $weights_sift2 ]]; then echo "Running SIFT2 on  $filtername  tractogram"
        tcksift2 -nthreads $threads $tractogram $fod $weights_sift2
  else echo "Subject ${id} has SIFT2 weights for  $filtername  tractogram "; fi



##------------------------------- Begin Main Processing Script-----------------------------##

# Timer
 aloita=$(date +%s)


  echo "*------------------- Begin structural connectome generation for subject: $id ---------------------*"
  echo "  Tractograms: $filtername & $filtername_commit"

# -------------------------------------------------------------------------------------------------------------------------------------------------

# File names
  weights_fa=${dir_weights}/${id}_${STAT}-fa.csv
  weights_qt1=${dir_weights}/${id}_${STAT}-qt1.csv
  weights_rd=${dir_weights}/${id}_${STAT}-rd.csv
  weights_ad=${dir_weights}/${id}_${STAT}-ad.csv
  weights_adc=${dir_weights}/${id}_${STAT}-adc.csv
  weights_icvf_noddi=${dir_weights}/${id}_${STAT}-icvf-noddi.csv
  weights_od_noddi=${dir_weights}/${id}_${STAT}-od-noddi.csv
  weights_icvf_commit=${dir_weights}/${id}_${STAT}-icvf-commit.csv 				# sampling of volumetric ICVF image output by COMMIT

  if [[ ! -d ${dir_weights} ]]; then mkdir -p ${dir_weights} ; fi

# Sample streamlines
  if [[ ! -f ${weights_fa} ]]; then echo "Sampling fractional anisotropy"
        tcksample $tractogram_commit $dti_fa $weights_fa -stat_tck $STAT
  fi

  if [[ ! -f $weights_qt1 ]]; then echo "Sampling qT1"
        tcksample $tractogram_commit $qT1_in_dwi_NL $weights_qt1 -stat_tck $STAT
  fi

  if [[ ! -f $weights_rd ]]; then echo "Sampling radial diffusivity"
        tcksample $tractogram_commit $dti_rd $weights_rd -stat_tck $STAT
  fi

  if [[ ! -f $weights_ad ]]; then echo "Sampling axial diffusivity"
        tcksample $tractogram_commit $dti_ad $weights_ad -stat_tck $STAT
  fi

  if [[ ! -f $weights_adc ]]; then echo "Sampling the apparent diffusion coefficient"
        tcksample $tractogram_commit $dti_adc $weights_adc -stat_tck $STAT
  fi

  if [[ ! -f $weights_icvf_noddi ]]; then echo "Sampling intra-cellular volume fraction from NODDI"
        tcksample $tractogram_commit $noddi_icvfMap $weights_icvf_noddi -stat_tck $STAT
  fi

  if [[ ! -f $weights_od_noddi ]]; then echo "Sampling orientation dispersion from NODDI"
        tcksample $tractogram_commit $noddi_odMap $weights_od_noddi -stat_tck $STAT
  fi

  if [[ ! -f ${weights_icvf_commit} ]]; then echo "Sampling icvf-commit weights for subject $id"
  	tcksample $tractogram_commit $map_icvf_commit $weights_icvf_commit -stat_tck $STAT
  fi

# -----------------------------------------------------------------------------------------------
# ------------------- Connectomes computed from COMMIT-filtered Tractogram ----------------------
# -----------------------------------------------------------------------------------------------

declare -a connectomeWeightArray=("fa" "qt1" "rd" "ad" "adc" "icvf-noddi" "od-noddi" "icvf-commit")

# Build TRACTOMETRY & COMMIT connectomes

  if [[ ! -d ${dir_connectome} ]]; then mkdir -p ${dir_connectome} ; fi
  if [[ ! -d ${dir_assignments} ]]; then mkdir -p ${dir_assignments} ; fi

  parcellations=$(find ${dir_volum} -name "*.nii.gz" ! -name "*cerebellum*" ! -name "*subcortical*")
  echo "*------------------  Building full (cortical-subcortical-cerebellum) connectomes: filtered by $filtername_commit  ------------------------------*"

  for seg in $parcellations ; do
        parc_name=$(echo "${seg/.nii.gz/}" | awk -F 'nativepro_' '{print $2}')
        echo "Parcellation: ** $parc_name **"
        nom=${dir_connectome}/${id}_${tracts}_${filtername_commit}_full_${parc_name}
        assignments=${dir_assignments}/${id}_${tracts}_${filtername_commit}_full_${parc_name}_nodeassignments.txt
        lut="${util_lut}/lut_${parc_name}_mics.csv"
        dwi_all=$tmp/${id}_${parc_name}-full_dwi.nii.gz                                                                         # Segmentation in dwi space

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


      # COMMIT
        conn="${nom}_COMMIT.txt"
        if [ ! -f $conn ]; then echo "Connectome: $conn"
                tck2connectome -nthreads $threads $tractogram_commit $dwi_all $conn \
                               -tck_weights_in $weights_commit -scale_length \
                               -scale_invnodevol -quiet
                Rscript ${MICAPIPE}/functions/connectome_slicer.R --conn=$conn --lut1=${lut_sc} --lut2=${lut} --mica=${MICAPIPE}
        else echo "$conn ALREADY EXISTS!" ; fi
        if [[ -f $conn ]]; then ((Nparc++)); fi


      # LoS
        conn="${nom}_los.txt"
        if [ ! -f $conn ]; then echo "Connectome: $conn"
                tck2connectome -nthreads $threads $tractogram_commit $dwi_all $conn \
                               -scale_length -stat_edge mean -quiet
                Rscript ${MICAPIPE}/functions/connectome_slicer.R --conn=$conn --lut1=${lut_sc} --lut2=${lut} --mica=${MICAPIPE}
        else echo "$conn ALREADY EXISTS!" ; fi
	if [[ -f $conn ]]; then ((Nparc++)); fi

      # weighted connectomes
        for connectomeWeight in "${connectomeWeightArray[@]}" ; do
                WEIGHT="$STAT-$connectomeWeight"
                conn="${nom}_$WEIGHT.txt"
                edgeWeights="${dir_weights}/${id}_$WEIGHT.csv"
                if [ ! -f $conn ]; then echo "Connectome: $conn"
                        tck2connectome -nthreads $threads $tractogram_commit $dwi_all $conn \
                                       -scale_file $edgeWeights -stat_edge mean -symmetric -zero_diagonal -quiet
                        Rscript ${MICAPIPE}/functions/connectome_slicer.R --conn=$conn --lut1=${lut_sc} --lut2=${lut} --mica=${MICAPIPE}
                else echo "$conn ALREADY EXISTS!" ; fi
		if [[ -f $conn ]]; then ((Nparc++)); fi
        done
  done


# -----------------------------------------------------------------------------------------------
# ------------------- Connectomes computed from gm-filtered Tractogram --------------------------
# -----------------------------------------------------------------------------------------------

# Build NoS & SIFT2 connectomes

  parcellations=$(find ${dir_volum} -name "*.nii.gz" ! -name "*cerebellum*" ! -name "*subcortical*")
  echo "*------------------  Building full (cortical-subcortical-cerebellum) connectomes: filtered by $filtername  ------------------------------*"

  for seg in $parcellations ; do
        parc_name=$(echo "${seg/.nii.gz/}" | awk -F 'nativepro_' '{print $2}')
        echo "Parcellation: ** $parc_name **"
        nom=${dir_connectome}/${id}_${tracts}_${filtername}_full_${parc_name}
        assignments=${dir_assignments}/${id}_${tracts}_${filtername}_full_${parc_name}_nodeassignments.txt
        lut="${util_lut}/lut_${parc_name}_mics.csv"
        dwi_all=$tmp/${id}_${parc_name}-full_dwi.nii.gz                                                                         # Segmentation in dwi space

        # NoS
        conn="${nom}_nos.txt"
        if [ ! -f $conn ]; then echo "Connectome: $conn"
                tck2connectome -nthreads $threads $tractogram $dwi_all $conn \
                               -scale_invnodevol -out_assignments $assignments -quiet
                Rscript ${MICAPIPE}/functions/connectome_slicer.R --conn=$conn --lut1=${lut_sc} --lut2=${lut} --mica=${MICAPIPE}
        else echo "$conn ALREADY EXISTS!" ; fi
        if [[ -f $conn ]]; then ((Nparc++)); fi

      # SIFT2
        conn="${nom}_SIFT2.txt"
        if [ ! -f $conn ]; then echo "Connectome: $conn"
                tck2connectome -nthreads $threads $tractogram $dwi_all $conn \
                               -tck_weights_in $weights_sift2 \
                               -scale_invnodevol -quiet
                Rscript ${MICAPIPE}/functions/connectome_slicer.R --conn=$conn --lut1=${lut_sc} --lut2=${lut} --mica=${MICAPIPE}
        else echo "$conn ALREADY EXISTS!" ; fi
        if [[ -f $conn ]]; then ((Nparc++)); fi

        done
  done
  chmod 777 -R ${dir_connectome}/*


# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------


# Notification of completion
  eri=$(echo "$lopuu - $aloita" | bc)
  eri=$(echo print $eri/60 | perl)

  echo "*------------ Connectome processing for $id complete -------------*"
  echo "  Number of connectomes generated	: $Nparc"
  echo "  Proc time				: $eri minutes"
