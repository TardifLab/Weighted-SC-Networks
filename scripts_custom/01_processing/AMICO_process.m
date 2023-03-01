%AMICO_process(protocol,subject,maskin,scheme)
function AMICO_process(protocol,subject,diff,maskin,scheme)

%%  annoying the way they want to set this... check if i can modify the way the paths are specified...

%% protocol: dir of specific study e.g. /data/mril/mril10/ilana/for-others/tomas (so rotation matrices only need to be computed once)
%%		
%% subject : where data is for each sujbect e.g. 602A/Diffusion_NODDI
%% diff: diffusion file (unzipped, no path)
%% maskin: mask(no path)
%% scheme: NODDI scheme file (path from AMICO_data_path) (can be created from bvecs bvals with  bvecs_bvals2camino.pl)

% Add the stuff for AMICO set-up here

global AMICO_code_path AMICO_data_path CAMINO_path CONFIG
global niiSIGNAL niiMASK
global KERNELS bMATRIX

% Path definition: adapt these to your needs
% ==========================================
AMICO_code_path = '/YOUR/PATH/TO/matlab/AMICO-master/matlab';
AMICO_data_path = protocol;
CAMINO_path     = '/YOUR/PATH/TO/spams-matlab/absolute path to CAMINO bin folder';
NODDI_path      = '/YOUR/PATH/TO/NODDI_toolbox2/NODDI_toolbox_v0.9';
SPAMS_path      = '/YOUR/PATH/TO/matlab/spams-matlab/';

if ~isdeployed
    addpath( genpath(NODDI_path) )
    addpath( fullfile(SPAMS_path,'build') )
    addpath( fullfile(AMICO_code_path,'kernels') )
    addpath( fullfile(AMICO_code_path,'models') )
    addpath( fullfile(AMICO_code_path,'optimization') )
    addpath( fullfile(AMICO_code_path,'other') )
    addpath( fullfile(AMICO_code_path,'vendor','NIFTI') )
end

AMICO_SetSubject( '', subject );
AMICO_PrecomputeRotationMatrices();
CONFIG.dwiFilename    = fullfile( CONFIG.DATA_path, diff);
CONFIG.maskFilename   = fullfile( CONFIG.DATA_path, maskin );
CONFIG.schemeFilename = fullfile(scheme);
AMICO_LoadData
AMICO_SetModel( 'NODDI' );
AMICO_GenerateKernels( false );
AMICO_ResampleKernels();
AMICO_Fit()


