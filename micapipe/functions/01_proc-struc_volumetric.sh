#!/bin/bash
#
# T1w Structural processing with bash:
#
# Preprocessing workflow for structural T1w.
#
# This workflow makes use of ANTS, FSL and AFNI
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

BIDS=$1
id=$2
out=$3
SES=$4
PROC=$5
nocleanup=$6
threads=$7
here=`pwd`

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

#------------------------------------------------------------------------------#
Title "Running MICA structural processing: Volumetric"
micapipe_software
# print the names on the terminal
bids_print.variables
Info "Not erasing temporal dir: $nocleanup"

# GLOBAL variables for this script
Info "ANTs will use $threads threads"

#	Timer
aloita=$(date +%s)

# if temporary directory is empty
#if [ -z ${tmp} ]; then tmp=/tmp; info "temp dir was unset: is now $tmp"; else info "temp dir was set to $tmp";  fi
# Create temporal directory
#tmp=${tmp}/${RANDOM}_micapipe_proc_struc-vol_${subject}
#if [ ! -d $tmp ]; then Do_cmd mkdir -p $tmp; fi

# Hardcode temp dir (Found this to work better)
  tmp=${proc_struct}/tmp_01volumetric
  if [ ! -d ${tmp} ]; then Do_cmd mkdir -p ${tmp} ;  fi

# TRAP in case the script fails
trap cleanup INT TERM

# BIDS T1w processing
N=${#bids_T1ws[@]} # total number of T1w
n=$((${N} - 1))

# FSL tries to submit to SGE
# unset FSLPARALLEL #unset SGE_ROOT         <<<<<<<<<<<<< This is not working

# Creates the t1w_nativepro for structural processing
if [ ! -f ${proc_struct}/${id}_t1w_*mm_nativepro.nii.gz ] || [ ! -f ${proc_struct}/${id}_t1w_*mm_nativepro_brain.nii.gz ]; then
    # Reorient  to LPI with AFNI
    # LPI is the standard 'neuroscience' orientation, where the x-axis is
    # Left-to-Right, the y-axis is Posterior-to-Anterior, and the z-axis is Inferior-to-Superior.
    for ((k=0; k<=$n; k++)); do
      nom=`echo ${bids_T1ws[k]} | awk -F '/' '{print $NF}'`
      t1_reo=${tmp}/${nom/sub/reo}
      t1_obl=${tmp}/${nom/sub/obl}
      # Deoblique
      Do_cmd 3dresample -orient LPI -prefix $t1_obl -inset ${bids_T1ws[k]}
      # Reorient to standart
      Do_cmd fslreorient2std $t1_obl $t1_reo
    done

    # If multiple T1w were provided, Register and average to the first T1w
    if [ $N -gt 1 ]; then
      ref=${bids_T1ws[0]} # reference to registration
      ref_run=`echo ${bids_T1ws[0]} | awk -F 'run-' '{print $2}'| sed 's:_T1w.nii.gz::g'`
      ref_res=`echo ${bids_T1ws[0]} | awk -F 'run-' '{print $2}'| sed 's:_T1w.nii.gz::g'`
      # Variables for N4BiasFieldCorrection & Native output
      T1str_nat=`t1w_str ${id} ${ref} nativepro`
      T1n4=${tmp}/${T1str_nat}_n4.nii.gz
      # Loop over each T1
      for ((i=1; i<=$n; i++)); do
          run=`echo ${bids_T1ws[i]} | awk -F 'run-' '{print $2}'| sed 's:_T1w.nii.gz::g'`
          t1ref=run-${ref_run}
          T1mat_str=${dir_warp}/${id}_t1w_run-${run}_to_${t1ref}_
          T1mat=${T1mat_str}0GenericAffine.mat
          T1run_2_T1=${tmp}/${id}_t1w_run-${run}_to_${t1ref}.nii.gz

          Do_cmd antsRegistrationSyN.sh -d 3 -m ${bids_T1ws[i]} -f $ref  -o $T1mat_str -t a -n $threads -p d
          Do_cmd antsApplyTransforms -d 3 -i ${bids_T1ws[i]} -r $ref -t $T1mat -o $T1run_2_T1 -v -u int
      done
      # Calculate the mean over all T1w registered to the 1st t1w run
      t1s_reg=$(find ${tmp}/*_to_${t1ref}.nii.gz)
      t1s_add=$(echo "-add $(echo ${t1s_reg} | sed 's: : -add :g')")
      Do_cmd fslmaths $ref $t1s_add -div $N $T1n4 -odt float

    # If only one T1w is provided
    elif [ $N -eq 1 ]; then
      # Variables for N4BiasFieldCorrection & Native output
      T1str_nat=`t1w_str ${id} ${t1_reo} nativepro`
      T1n4=${tmp}/${T1str_nat}_n4.nii.gz
      Do_cmd cp -v $t1_reo $T1n4
    fi

    # Output names
    T1nativepro=${proc_struct}/${T1str_nat}.nii.gz
    T1nativepro_brain=${T1nativepro/.nii.gz/_brain.nii.gz}
    T1nativepro_first=${proc_struct}/first/${T1str_nat}.nii.gz
    T1nativepro_5tt=${T1nativepro/.nii.gz/_5TT.nii.gz}

    # Intensity Non-uniform correction - N4
    Do_cmd N4BiasFieldCorrection  -d 3 -i $T1n4 -r -o $T1n4 -v

    # Rescale intensity [100,0]
    Do_cmd ImageMath 3 $T1nativepro RescaleImage $T1n4 0 100

    # If no T1native pro exit_status "something is wrong" exit
    if [ ! -f ${T1nativepro} ]; then Error "$T1str_nat was not generated"; Do_cmd exit; fi

    # Brainmask
    Do_cmd bet $T1nativepro $T1nativepro_brain  -B -f 0.25 -v

    # If no T1native pro exit_status "something is wrong" exit
    if [ ! -f ${T1nativepro_brain} ]; then Error "$T1str_nat masked was not generated"; Do_cmd exit; fi

else
    Info "Subject ${id} has a t1w_nativepro and t1w_nativepro_brain"
    # Output names get the names
    T1str_nat=`t1w_str ${id} ${proc_struct}/${id}_t1w_*mm_nativepro.nii.gz nativepro`
    T1nativepro=${proc_struct}/${T1str_nat}.nii.gz
    T1nativepro_brain=${T1nativepro/.nii.gz/_brain.nii.gz}
    T1nativepro_first=${proc_struct}/first/${T1str_nat}.nii.gz
    T1nativepro_5tt=${T1nativepro/.nii.gz/_5TT.nii.gz}
fi

# FSL first on the t1w_nativepro
 unset SGE_ROOT
 export FSLPARALLEL=0
 firstout=${T1nativepro_first/.nii.gz/_all_fast_firstseg.nii.gz}
 if [ ! -f ${firstout} ]; then
     Info "FSL first is running, output file:\n\t\t\t ${firstout}"
     Do_cmd run_first_all -v -i $T1nativepro_brain -o $T1nativepro_first -b &
    # Wait until $firstout exists, IF it was sent to SGE
     wait $!
     until [ -f $firstout ]; do sleep 120; done

     Info "Changing FIRST output names to maintain MICA-BIDS naming convention"
     for i in ${proc_struct}/first/*pro-*; do mv -v $i ${i/pro-/pro_}; done
     mv -v ${proc_struct}/${T1str_nat}_brain_to_std_sub.nii.gz ${proc_struct}/${id}_t1w_1mm_MNI152_brain_affine.nii.gz
     mv -v ${proc_struct}/${T1str_nat}_brain_to_std_sub.mat ${dir_warp}/${T1str_nat}_brain_to_1mm_MNI152_brain.mat
   else
     Info "Subject ${id} has FSL-first"
 fi

# FSL FAST on the t1w_nativepro
if [ ! -f ${T1nativepro_brain/.nii.gz/_pve_0.nii.gz} ]; then
    Do_cmd fast -N -v $T1nativepro_brain
else
    Info "Subject ${id} has FSL-fast"
fi

# Loop over all requested templates - could use a check for whether the template and the mask exists.
# mmTemplates is a fixed value=(0.8 2) of the MNI152 atlas resolution
for mm in 2 0.8; do
  # Only runs if the output doesn't exist
  if [ ! -f ${proc_struct}/${id}_t1w_${mm}mm_MNI152_brain_pve_0.nii.gz ]; then
      # ---------------------------------------------------------
      # MNI152 templates
      MNI152_brain=${util_MNIvolumes}/MNI152_T1_${mm}mm_brain.nii.gz
      MNI152_mask=${util_MNIvolumes}/MNI152_T1_${mm}mm_brain_mask.nii.gz

      # T1 to MNI152 transformation matrix and warp fields
      mat_MNI152_SyN=${dir_warp}/${T1str_nat}_brain_to_${mm}mm_MNI152_SyN_brain_
      T1_MNI152_warp=${mat_MNI152_SyN}1Warp.nii.gz
      T1_MNI152_warped=${mat_MNI152_SyN}Warped.nii.gz
      T1_MNI152_affine=${mat_MNI152_SyN}0GenericAffine.mat
      # Outputs
      T1_MNI152=${proc_struct}/${id}_t1w_${mm}mm_MNI152.nii.gz
      T1_MNI152_brain=${proc_struct}/${id}_t1w_${mm}mm_MNI152_brain.nii.gz

      # ANTs - 3 steps registration: rigid, affine and SyN
      Do_cmd antsRegistrationSyN.sh -d 3 -f $MNI152_brain -m $T1nativepro_brain -o ${mat_MNI152_SyN} -t s -n $threads

      # Nativepro to MNI152, first the brain then the full volume
      Do_cmd mv $T1_MNI152_warped $T1_MNI152_brain

  else
      Info "Subject ${id} has t1w_${mm}mm_nativepro on MNI152 space and FSL-fast "
  fi
done

# Generate a five-tissue-type image for anatomically constrained tractography
if [[ ! -f $T1nativepro_5tt ]]; then
    # --------------------------------------------------------------
    # Step by step 5tt
    # Process data in tmp directory
    cd ${tmp}
    # Convert to mif format
    Do_cmd mrconvert $T1nativepro_brain ${tmp}/input.mif
    # Change strides
    first_input='T1.nii'
    Do_cmd mrconvert input.mif $first_input -strides -1,+2,+3
    Info "Convert FIRST meshes to partial volume images:"
    firstr="${proc_struct}/first/${T1str_nat}_"
    sgm_structures=(L_Accu R_Accu L_Caud R_Caud L_Pall R_Pall L_Puta R_Puta L_Thal R_Thal)
    echo "         	${sgm_structures[@]}"
    for struct in "${sgm_structures[@]}"; do
      Info "Mesh to volume: $struct"
      pve_image_path=mesh2voxel_${struct}.mif
      vtk_in_path=${firstr}${struct}_first.vtk
      vtk_temp_path=${struct}.vtk
      Do_cmd meshconvert $vtk_in_path $vtk_temp_path -transform first2real $first_input
      Do_cmd mesh2voxel $vtk_temp_path $T1nativepro_brain $pve_image_path
    done
    mrmath mesh2voxel_* sum - | mrcalc - 1.0 -min all_sgms.mif

    Info "Generating partial volume images for SGM structures"
    T1_pve_1=${T1nativepro/.nii.gz/_brain_pve_1.nii.gz}
    T1_pve_2=${T1nativepro/.nii.gz/_brain_pve_2.nii.gz}
    T1_pve_0=${T1nativepro/.nii.gz/_brain_pve_0.nii.gz}

    mrthreshold $T1_pve_2 - -abs 0.001 | maskfilter - connect - -connectivity | mrcalc 1 - 1 -gt -sub remove_unconnected_wm_mask.mif -datatype bit
    Do_cmd mrcalc $T1_pve_0 remove_unconnected_wm_mask.mif -mult csf.mif
    Do_cmd mrcalc 1.0 csf.mif -sub all_sgms.mif -min sgm.mif
    Do_cmd mrcalc 1.0 csf.mif sgm.mif -add -sub $T1_pve_1 $T1_pve_2 -add -div multiplier.mif
    Do_cmd mrcalc multiplier.mif -finite multiplier.mif 0.0 -if multiplier_noNAN.mif
    Do_cmd mrcalc $T1_pve_1 multiplier_noNAN.mif -mult remove_unconnected_wm_mask.mif -mult cgm.mif
    Do_cmd mrcalc $T1_pve_2 multiplier_noNAN.mif -mult remove_unconnected_wm_mask.mif -mult wm.mif
    Do_cmd mrcalc 0 wm.mif -min path.mif
    mrcat cgm.mif sgm.mif wm.mif csf.mif path.mif - -axis 3 | mrconvert - combined_precrop.mif -strides +2,+3,+4,+1
    Do_cmd mv combined_precrop.mif result.mif
    Do_cmd mrconvert result.mif $T1nativepro_5tt
    Do_cmd cd $here
else
    Info "Subject ${id} has 5TT nifti"
fi

# -----------------------------------------------------------------------------------------------
# Clean temporal directory
if [[ $nocleanup == "FALSE" ]]; then Do_cmd rm -rf $tmp; else Info "Mica-pipe tmp directory was not erased: \n\t\t\t${tmp}"; fi

# QC notification of completition
lopuu=$(date +%s)
eri=$(echo "$lopuu - $aloita" | bc)
eri=`echo print $eri/60 | perl`

# Notification of completition
Title "Volumetric tructural processing ended in \033[38;5;220m `printf "%0.3f\n" ${eri}` minutes \033[38;5;141m:\n\tlogs:
`ls ${dir_logs}/proc-volumetric_*.txt`"
echo "${id}, proc_struc, COMPLETED, `whoami`, `uname -n`, $(date), `printf "%0.3f\n" ${eri}`, $PROC" >> ${out}/brain-proc.csv
bids_variables_unset
