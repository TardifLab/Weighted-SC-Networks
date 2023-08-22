%% A template analysis script

% Add path to required tools
% addpath(genpath('PATH/TO/YOUR/STUFF'))
addpath(genpath(pwd))                                                      


% Setup
runloc              = 'local';                                              % Describes desired path setup
filtname            = 'unfiltered';                                         % Describes which connectome filter type
seg                 = 'cor';                                                % Describes desired connectome segmentation
HCex                = [];                                                   % Subject to exclude
savetag             = 'helpful_label_';                                     % Label attached to derivatives
saveloc             = [];


%% Load your data (two example methods given below)

% Option 1: load the shared compiled connectome files
a1_loadConns(runloc,saveloc,savetag)


% Option 2: load your own connectome files
a1_compileConns(runloc,filtname)



%% Analysis 1: Log convert skewed data and normalize all data for these plots
% Connectomes & distributions are plot
a2_procConns1_lognorm(runloc,filtname,savetag,seg,HCex);



%% Analysis 2: no transforming or normalization of data
% All other plots generated here
a2_procConns2_nomod(runloc,filtname,savetag,seg,HCex)


