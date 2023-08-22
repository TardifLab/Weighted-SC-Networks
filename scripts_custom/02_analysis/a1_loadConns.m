function a1_loadConns(varargin)
%
%  An example script for loading the shared individual connectome files and
%  resaving them as  single file.
%  Also loads shared metadata files and resaves them with a common tag
%  for easy use throughout.
% 
% Input
%                           + Optional +
%   runloc   : switch used for path setup (string) (default = local)
%   saveloc  : new derivative location (default = original location)
%   savetag  : label attached to all saved outputs
%   
%
% Output
%   P        : paths (structure)
%
% 2023 Mark C Nelson, McConnell Brain Imaging Centre, MNI, McGill
%--------------------------------------------------------------------------

%% Optional inputs
nin=max(nargin,1); 
defaults                    = {'local' [] 'allWeights'};
[runloc,saveloc,savetag]    = INhandler(varargin,nin,defaults);


%% Setup

% Set Paths
P                           = a0_pathDefs(runloc);                          % Define Paths
addpath(genpath(P.HOME))                                                    % Add to Matlab
locDeriv                    = P.dataDerivative;

if isempty(saveloc); saveloc=locDeriv; end


%% Load weights, parcellations and subjects
load([locDeriv '/1_parcellations_SC_FC_schaefer-400.mat'])                  % Ps = parcellation labels
load([locDeriv '/2_weights_SC_FC.mat'])                                     % Ws = functional & structural labels/weights
load([locDeriv '/3_subjects_SC_FC.mat']);                                   % Ss = subject labels


%% Load files for each individual weighted connectome
tmpw8vec    = {'fc' 'commit' 'sift2' 'los' 'fa' 'icvf' 'qt1' 'rd' 'nos'};   % Just used to load each individual file
Cs          = cell(numel(tmpw8vec), numel(Ps));
for dd = 1 : numel(tmpw8vec)
    disp(['Loading connectomes weighted by: ' tmpw8vec{dd}])
    load([locDeriv '/0_connectomes_SC_FC_schaefer-400_' tmpw8vec{dd} '.mat'])
    Cs{dd}  = dtmp;
end


%% Save cell arrays for easy reuse
save([saveloc '/0_connectomes_'   savetag  '.mat'],'Cs')
save([saveloc '/1_parcellations_' savetag  '.mat'],'Ps')
save([saveloc '/2_weights_'       savetag  '.mat'],'Ws')
save([saveloc '/3_subjects_'      savetag  '.mat'],'Ss')

% Alternate save option for large files (~2GB or larger usually) 
% save([P.dataDerivative '/0_connectomes_'   svtag filtname  '.mat'],'Cs','-v7.3')


end




