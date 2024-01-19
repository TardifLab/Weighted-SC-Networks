
# <div align="center">Characterizing Weighted Structural Brain Networks</div> #
**<div align="center">by Mark Cameron Nelson</div>**  
  
This library of code can be used to   
* preprocess multi-modal MRI data,   
* compute a range of weighted connectomes and   
* run simple analyses summarized with pretty pictures.  
  
  But first, know before you start...

### How to cite this work ###
Nelson, Mark C., et al. (2023). The Human Brain Connectome Weighted by the Myelin Content and Total Intra-Axonal Cross-Sectional Area of White Matter Tracts. *Network Neuroscience*  
https://doi.org/10.1162/netn_a_00330  
  
Cruces, Raul R., et al. (2022). Micapipe: A pipeline for multimodal neuroimaging and connectome analysis. *NeuroImage*  
https://doi.org/10.1016/j.neuroimage.2022.119612  
  
Royer, Jessica, et al. (2022). An Open MRI Dataset For Multiscale Neuroscience. *Scientific Data*  
https://www.nature.com/articles/s41597-022-01682-y  
  

#### + Want to download the raw data I used? ####
  
&nbsp;&nbsp; [MICs dataset can be found here](https://portal.conp.ca/dataset?id=projects/mica-mics)  
  

#### + Prefer to just jump right in with my connectomes? ####
&nbsp;&nbsp; These can be found in the `derivatives` folder  
  

#### + Need some help using `micapipe`? ####
&nbsp;&nbsp; [micapipe.readthedocs.io](http://micapipe.readthedocs.io/en/latest/)  
  
  
  Ok let's get started!

## Step 1: data processing ##
This section can be used to preprocess raw structural and diffusion MRI data to estimate whole brain structural connectivity networks.  
It is composed primarily of a set of shell wrappers that interface with processing tools from `micapipe`
  
Note: The entire `micapipe` library is included, HOWEVER!!!   
  
**only the functions included in the processing order below were actually used i.e. some debugging may be necessary for the others**.  
  
In particular, `01_proc-struc_freesurfer.sh`, `02_proc-rsfmri.sh` & `03_FC.py` were not run, but will be necessary for structure-function investigations.  


####  +++  IMPORTANT!!  +++  ####  
The version of `micapipe` included here was derived from a beta version. The current code on the `micapipe` GitHub has undergone many changes in syntax since.  
  
**DO NOT ATTEMPT TO PULL THAT REPO OR USE CODE FROM IT DIRECTLY.**  
  
Any desired code from that repo will need to be modified manually to conform with the syntax of this version of `micapipe`.  
  

### Dependencies ###
See `scripts_custom/init.sh`  

| *Software* |    *Version*  | *Links* |
|------------|---------------|--------------|  
| dcm2niix   | v1.0.20190902 | https://github.com/rordenlab/dcm2niix |
| Freesurfer | 6.0.0         | https://surfer.nmr.mgh.harvard.edu/ |
| FSl        | 6.0           | https://fsl.fmrib.ox.ac.uk/fsl/fslwiki |
| AFNI       | 20.3.03       | https://afni.nimh.nih.gov/download |
| MRtrix3    | 3.0.1         | https://www.mrtrix.org |
| ANTs       | 2.3.4         | https://github.com/ANTsX/ANTs |
| workbench  | 1.4.2         | https://www.humanconnectome.org/software/connectome-workbench |
| FIX        | 1.06          | https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FIX |
| R          | 3.6.3         | https://www.r-project.org |
| python     | 3.7.6         | https://www.python.org/downloads/ |
  

### Order of operations ###
Processing scripts should be run in this order.  
  
  1. `scripts_custom/01_qbatch_subcall.sh`  
     - Used to run jobs on remote cluster, adapt to your needs  
  
  2. `scripts_custom/02_micapipe_call.sh`  
     - Used to call specific micapipe modules  
  
  3. `micapipe/functions/01_proc-struc_volumetric.sh`  
  
  4. `scripts_custom/seg4DConvert.sh`  
     - Workaround for occasional segmentation failure in 01_proc-struc_volumetric.sh  
  
  5. `micapipe/functions/02_post-structural.sh`  
  
  6. `micapipe/functions/02_proc-dwi.sh`  
  
  7. `scripts_custom/03_run_noddi.sh`  
  
  8. `micapipe/functions/03_SC.sh`  
  
  9. `scripts_custom/04_run_commit0.sh`  
  
  10. `scripts_custom/04_run_commit1.sh --> 04_run_commit2.py`  
  
  11. `scripts_custom/05_connectomes.sh`  
  

#### Additional comments ####
1. `scripts_custom/init.sh` must be modified to reflect your path stucture  
   - micapipe/functions/init.sh was not used  
  
2. `micapipe/functions/utilities.sh` must be modified to reflect your file naming conventions, especially for raw data  
  


## Step 2: analysis ##
This section can be used to perform simple post-processing operations on connectomes including filtering, segmenting, etc and to make some pretty plots.  
  
It is entirely based in **Matlab** and was written using the R2016a release (although it has been tested in R2022a).  
  
See `MAINSCRIPT.m` for a simple template.
  

### Dependencies ###

|    *Software*     |      *Version*      | *Links* |
|-------------------|---------------------|----------------|  
| BCT               | 2019.03.03          | https://www.nitrc.org/projects/bct |
| BrainNetViewer    | 1.7                 | https://www.nitrc.org/projects/bnv/ |
| BrainSpace        | 0.1.1               | https://github.com/MICA-MNI/BrainSpace |
| ENIGMA            | 2.0.3               | https://github.com/MICA-MNI/ENIGMA
| gifti-master      | 2.0                 | https://www.nitrc.org/projects/gifti/ |
| linspecer         | 1.4.0.0             | https://www.mathworks.com/matlabcentral/fileexchange/42673-beautiful-and-distinguishable-line-colors-colormap |
| suplabel          | 1.5.0.0             | https://www.mathworks.com/matlabcentral/fileexchange/7772-suplabel?s_tid=srchtitle |
| Violinplot-Matlab | -                   | https://github.com/bastibe/Violinplot-Matlab |
| xmltree-main      | 2.1                 | https://github.com/gllmflndn/xmltree |
  
  

