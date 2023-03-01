%% A template analysis script

addpath(genpath('PATH/TO/YOUR/STUFF'))                                      % Add path to required tools


% Setup
runloc              = 'local';                                              % Describes desired path setup
filtname            = 'unfiltered';                                         % Describes which connectome filter type
seg                 = 'cor';                                                % Describes desired connectome segmentation
HCex                = [];                                                   % Subject to exclude
svtag               = 'helpful_label_';                                     % Label attached to derivatives



%% Compile connectomes
a1_compileConns(runloc,filtname)                                            % These inputs are optional



%% Analysis 1: Log convert skewed data and normalize all data for these plots
% Connectomes & distributions are plot
a2_procConns1_lognorm(runloc,filtname,svtag,seg,HCex);



%% Analysis 2: no transforming or normalization of data
% All other plots generated here
a2_procConns2_nomod(runloc,filtname,svtag,seg,HCex)


