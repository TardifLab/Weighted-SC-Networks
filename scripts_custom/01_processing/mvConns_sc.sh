#!/bin/bash
#
# Often used to move connectomes around
#
# Inputs:
#	$1 : source directory
#	$2 : target directory
#
# note: files should conform to expected naming structure i.e.
# 	subjectid_nstreamline_filter_segmentation_parcellation_weight.txt e.g.
# 	"HC001_5M_unfiltered_full_glasser-360_median-fa.txt"
#       "HC050_5M_filt-commit_full_schaefer-100_COMMIT.txt"
# 2021 Mark C Nelson, McConnell Brain Imaging Centre, MNI, McGill
#------------------------------------------------------------------------------------------------------------------------------------

sourcedir=$1
targetdir=$2

cd "$sourcedir"
shopt -s globstar
for i in **/connectomes/*/*5M*full*.txt; do 				# careful to avoid nodeassignments files if they exist
	fnold=$(basename $i)
        id=$(echo $fnold | cut -d "_" -f1)
        tracts=$(echo $fnold | cut -d "_" -f2)
        filter=$(echo $fnold | cut -d "_" -f3)
        seg=$(echo $fnold | cut -d "_" -f4)
        parc=$(echo $fnold | cut -d "_" -f5)
        weight=$(echo $fnold | cut -d "_" -f6 | cut -d "." -f1)

	# new full file
	newdir=${targetdir}/${filter}/${parc}/${weight}
	fnnew=${newdir}/${id}_${parc}_${weight}_${tracts}_${filter}_${seg}.txt

	if [ ! -d ${newdir} ]; then mkdir -p ${newdir} ;  fi

	if [ ! -f ${fnnew} ]; then
        	cp --verbose $i $fnnew
#		mv -v $i $fnnew
	else
		echo "Skipping: $fnold"
	fi
done


