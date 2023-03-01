function a1_compileConns(varargin)
%
%  Compiles connectome files (.txt) into cell array and saves as .mat for
%  easier use.
%
% Notes: 
%   1. Directory structure and Filenames must follow expected pattern (or code must be adjusted accordingly)
%   2. All files will be stacked: Cant select subset until after cell array is created; 
%      shouldnt have non-connectome files in directory (else add to )
% 
% Input
%                           + Optional +
%   runloc   : switch used for path setup (string) (default = local)
%   filtname : string indicating which filtered data type to load (default='unfiltered')
%   svtag    : label used to save outputs
%   
%
% Output
%   P        : paths (structure)
%
% 2021 Mark C Nelson, McConnell Brain Imaging Centre, MNI, McGill
%--------------------------------------------------------------------------

%% Optional inputs
nin=max(nargin,1); 
defaults                = {'local' 'unfiltered' 'helpful_label_'};
[runloc,filtname,svtag] = INhandler(varargin,nin,defaults);


%% Compile connectomes into cell array

% Set Paths
P                       = a0_pathDefs(runloc);                              % Define Paths
addpath(genpath(P.HOME))                                                     % Add to Matlab

% Stack connectomes into cell array
[Cs,Ps,Ws,Ss]           = conn_stacker([P.dataSource '/' filtname]);        % Cs = weights x parcellations 

% Select subset of data
select_weights          = {'SIFT2' 'nos' 'COMMIT' 'los' 'qt1'};             % Change to suit your needs
select_parcs            = {'schaefer-100' 'schaefer-400' 'aparc'};          % Change to suit your needs
[Cs,Ws]                 = conn_selectData(Cs,Ws,select_weights,1);       
[Cs,Ps]                 = conn_selectData(Cs,Ps,select_parcs,2);


% Scale COMMIT by mean edge length if desired
% COMMIT edge weights should be: sum(Xi * Li) / Lmean
icommit                 = cellstrfind(Ws,'commit');                         % Index for COMMIT
ilos                    = cellstrfind(Ws,'los');                            % Index for edge length
Cs(icommit,:)           = commit_scaler(Cs(icommit,:), Cs(ilos,:)); 


% Save cell arrays for easy reuse
save([P.dataDerivative '/0_connectomes_'   svtag filtname  '.mat'],'Cs')
save([P.dataDerivative '/1_parcellations_' svtag filtname  '.mat'],'Ps')
save([P.dataDerivative '/2_weights_'       svtag filtname  '.mat'],'Ws')
save([P.dataDerivative '/3_subjects_'      svtag filtname  '.mat'],'Ss')

% Alternate save option for large files (~2GB or larger usually) 
% save([P.dataDerivative '/0_connectomes_'   svtag filtname  '.mat'],'Cs','-v7.3')


end




