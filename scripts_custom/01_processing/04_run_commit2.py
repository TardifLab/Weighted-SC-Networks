#!/PATH/TO/YOUR/commit/env/bin/python3

import sys
import numpy as np
import commit
from commit import trk2dictionary

"""
Novice py script for running the original commmit implementation

	see Daducci et al. 2013 & 2015 IEEE Trans Med Imaging

# 2021 Mark C Nelson, McConnell Brain Imaging Centre, MNI, McGill
#------------------------------------------------------------------------------------------------------------------------------------
"""

#-----------------------------------#
#------------- SETUP ---------------#
#-----------------------------------#

# Inputs
ID  		= sys.argv[1]
tractogram 	= sys.argv[2]
out_dir 	= sys.argv[3]
#settings 	= sys.argv[4]

print(".\n *** Initializing COMMIT for: ", ID)
print(".\n *** TRACTOGRAM: ", tractogram)
print(".\n *** OUTPUT DIRECTORY: ", out_dir)
#print(".\n *** IC REGULARIZATION: ", settings)

# Dirs
proj_dir     	= "/PATH/TO/YOUR/STUFF"
subProc_dir  	= proj_dir + "/derivatives/sub-" + ID  + "/ses-01/proc_dwi"
commit_dir 	= subProc_dir + "/commit"

# Files
dwi_b0 	  	= subProc_dir + "/" + ID  + "_dwi_b0.nii.gz"
wm_fod         	= commit_dir + "/tmp/" + ID  + "_wm_fod_norm.nii.gz"
wm_mask        	= commit_dir + "/tmp/" + ID  + "_dwi_wm_mask.nii.gz"
dwi_corr       	= commit_dir + "/tmp/" + ID  + "_dwi_corr.nii.gz"
scheme         	= commit_dir + "/commit.scheme"                                 	# bvals manually adjusted

#------------------------------------
# Import usual COMMIT structure
#------------------------------------
commit.core.setup()                                                                     # precomputes the rotation matrices used internally by COMMIT
trk2dictionary.run(
        filename_tractogram     = tractogram,
        filename_peaks          = wm_fod,
        filename_mask           = wm_mask,
        TCK_ref_image           = dwi_b0,
        path_out                = out_dir ,
        fiber_shift             = 0.5,
        peaks_use_affine        = True
)

# load data
mit = commit.Evaluation( out_dir, '.' )                                                 # study_path, subject (relative to study_path)
mit.load_data(
        dwi_filename    = dwi_corr,
        scheme_filename = scheme
)

# set forward model (default params used here)
mit.set_model( 'StickZeppelinBall' )                                                    # model described in (Panagiotaki et al., NeuroImage, 2012)
d_par   = 1.7E-3                                                                        # Parallel diffusivity [mm^2/s]
d_perps = [ 0.51E-3 ]                                                                   # Perpendicular diffusivity(s) [mm^2/s]
d_isos  = [ 1.7E-3, 3.0E-3 ]                                                            # Isotropic diffusivity(s) [mm^2/s]
mit.model.set( d_par, d_perps, d_isos )
mit.generate_kernels( regenerate=True )
mit.load_kernels()

# Load dictionary (sparse data structure)
mit.load_dictionary( out_dir )

# Build linear operator A
mit.set_threads()                                                                       # use max possible; mit.set_threads( n ) to set manually
mit.build_operator()

# Fit model to data
mit.fit( tol_fun=1e-4, max_iter=1000 )                  				# stopping criterion = rel_error < 1e-4 | 1000 iterations

# Store results
mit.save_results()

#------------------------------------------------------------------------------------------------------------------
##---------------------- CODE BELOW THIS LINE COULD BE USED TO EXTEND COMMIT --------------------------------------
#------------------------------------------------------------------------------------------------------------------
#
# (saving IC coeffs to file)
# x_ic, x_ec, x_iso = mit.get_coeffs() 		# to get the estimated coeffs from all compartments
#
# save_path = out_dir + '/Results_StickZeppelinBall/coeffs_ic.csv'
#
# out_file = open(save_path, 'w') 		# opens file in write mode (caution, will clear existing files!)
#
# x_ic_list = x_ic.tolist()
# x_ic_str = str(x_ic_list)
#
# for ii in x_ic_str:
# 	out_file.write(ii)
#-----------------------------------------------------------------------------------------------------------------



#---------------------------------------------------------------------------
# Define the regularisation term  (each compartment must be done separately)

# The user can choose among the following penalties:

        # $\sum_{g\in G}w_g|x_g|_k$ : commit.solvers.group_sparsity with $k\in {2, \infty}$ (only for IC compartment)

        # $|x|_1$ : commit.solvers.norm1

        # $|x|_2$ : commit.solvers.norm2

        # $\iota_{\ge 0}(x)$ : commit.solvers.non_negative (Default for all compartments)


# If the chosen regularisation for the IC compartment is $\sum_{g\in G}|x_g|_k$,
# we can define $k$ via the group_norm field, which must be

        # $|x|_2$ : commit.solvers.norm2


# In this example we consider the following penalties:

        # Intracellular: group sparsity with 2-norm of each group

        # Extracellular: 2-norm

        # Isotropic: 1-norm
#---------------------------------------------------------------------------

#regnorms = [commit.solvers.group_sparsity,
#           commit.solvers.norm2,
#           commit.solvers.norm1
#]
#group_norm = 2                                                 # each group is penalised with its 2-norm

# Specify lambdas (regularisation params)
#lambdas = [10.,10.,10.]                                        # (do not consider this choice standard)

#------------------------------------
# Call the constructor of the data structure
#------------------------------------
#regterm = commit.solvers.init_regularisation(
#    mit,
#    regnorms    = regnorms,
#    structureIC = structureIC,
#    weightsIC   = weightsIC,
#    group_norm  = group_norm,
#    lambdas     = lambdas)

# perform optimization
#mit.fit(regularisation=regterm, max_iter=1000)

# saves to out_dir + /Results_StickZeppelinBall + path_suffix + /*
#mit.save_results(path_suffix='_AdvancedSolvers')

