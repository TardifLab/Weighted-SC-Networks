#!/bin/bash
#
# Main routing point to call all processing-related functions in MICAPIPE and SCRIPTS_CUSTOM folders
# Change the FUNC_ID variable to choose your function.
# All function calls are passed through the cluster (qbatch)
#  * * NOTE: THIS IS INTENDED TO BE RUN ON THE BIC AT MCGILL! PLEASE ADAPT FOR YOUR PURPOSES
#
# Logs are stored in $log_dir.
# Additional option to gauge resource needs using /usr/bin/time (see bottom)
#
# Input:
#	$1 : FUNC_ID, tag indicating which function to send to qbatch { volumetric, post_structural, dwi, SC, FC, noddi, icvf_conn, conn_slice, tck_gen }
#
#
#
# 2021 Mark C Nelson, McConnell Brain Imaging Centre, MNI, McGill
#------------------------------------------------------------------------------------------------------------------------------------

root_dir=YOUR/PATH/TO/THESE/SCRIPTS
scripts_dir=${root_dir}/scripts_custom
log_dir=${root_dir}/logs

# Choose function you wish to call for all subs
FUNC_ID=$1  						#{ volumetric, post_structural, dwi, noddi, SC, FC, noddi, icvf_conn, conn_slice, tck_gen }

# Set options based on desired function
if [ $FUNC_ID == volumetric ]        ; then VM=9  ;  Fn=1  ; script=${scripts_dir}/02_micapipe_call.sh ;
elif [ $FUNC_ID == post_structural ] ; then VM=3  ;  Fn=2  ; script=${scripts_dir}/02_micapipe_call.sh ;
elif [ $FUNC_ID == dwi ]             ; then VM=8  ;  Fn=3  ; script=${scripts_dir}/02_micapipe_call.sh ;
elif [ $FUNC_ID == noddi ]           ; then VM=10 ;  Fn=4  ; script=${scripts_dir}/03_run_noddi.sh ;
elif [ $FUNC_ID == SC ]              ; then VM=5  ;  Fn=5  ; script=${scripts_dir}/02_micapipe_call.sh ; 	# VM=5 for 5M tck;  VM=15 for 40M tck
elif [ $FUNC_ID == commit_prep ]     ; then VM=5  ;  Fn=6  ; script=${scripts_dir}/04_run_commit0.sh ;
elif [ $FUNC_ID == commit ]          ; then VM=40 ;  Fn=7  ; script=${scripts_dir}/04_run_commit1.sh ;         	# very memory intense, not ideal if tck size > 10-15M streamlines, VM=30 works most of the time
elif [ $FUNC_ID == connectomes ]     ; then VM=10 ;  Fn=8  ; script=${scripts_dir}/05_connectomes.sh ; 		# VM=3 works most of the time
elif [ $FUNC_ID == FC ]              ; then VM=10 ;  Fn=9  ; script=${scripts_dir}/02_micapipe_call.sh ;
elif [ $FUNC_ID == conn_slice ]      ; then VM=5  ;  Fn=00 ; script=${scripts_dir}/connSlicer_qbatch.sh ;
fi

# directory to store log files
log_func_dir=${log_dir}/${Fn}_${FUNC_ID}

for SUB in 13 ; do 												# Option to call single subjects
#for SUB in {02..51} ; do 											# Option to call group

	ID=HC0"$SUB"
	sub_dir=${log_func_dir}/$ID
	if [ ! -d ${sub_dir} ]; then mkdir -p ${sub_dir} ;  fi
	cd ${sub_dir} 												# cd or logs will output in cwd

	# Call to desired function
#	qbatch -verbose -l h_vmem=${VM}G -N "S${SUB}_${Fn}" ${script} $SUB $FUNC_ID 				# Standard call
	qbatch -verbose -l h_vmem=${VM}G -N "S${SUB}_${Fn}" /usr/bin/time --verbose ${script} $SUB $FUNC_ID 	# option to gauge resource allocation
done

