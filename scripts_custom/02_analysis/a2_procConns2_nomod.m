function a2_procConns2_nomod(varargin)
%
%  Processing stream 2: 
%   1. Does NOT log transform or normalize data
%   2. Computes & plots quantitative metrics
% 
% Input
%                           + Optional +
%   runloc   : switch for path setup (string) (default = local)
%   filtname : string indicating which filtered data type to load 
%              {'unfiltered' 'filt-gm' 'filt-commit'}; (default='unfiltered')
%   svtag    : label attached to derivatives
%   seg      : desired portion (segmentation) of connectomes to retain 
%             {'full' 'cor' 'sub'}; default='cor'
%   HCex     : index of subject(s) to exclude (default=[])
%
% 2021 Mark C Nelson, McConnell Brain Imaging Centre, MNI, McGill
%--------------------------------------------------------------------------

%% Optional inputs
nin=max(nargin,1); 
defaults                            = {'local' 'unfiltered' 'helpful_label_' 'cor' []};
[runloc,filtname,svtag,seg,HCex]    = INhandler(varargin,nin,defaults);

%% Load Pre-compiled connectomes
P               = a0_pathDefs(runloc);                                 
addpath(genpath(P.HOME))
locLabel        = P.dataLabel;
locDeriv        = P.dataDerivative; clear P

% load data
load([locDeriv '/0_connectomes_'   svtag filtname '.mat']);                 % Cs = functional & structural connectomes
load([locDeriv '/1_parcellations_' svtag filtname '.mat'])                  % Ps  = parcellation labels
load([locDeriv '/2_weights_'       svtag filtname '.mat'])                  % Ws = functional & structural labels/weights
load([locDeriv '/3_subjects_'      svtag filtname '.mat']); Ss(HCex)=[];    % Ss = subject labels


% Store backups (helpful while debugging/testing)
C_backup=Cs; W_backup=Ws; P_backup=Ps; S_backup=Ss;
% Cs=C_backup; Ps=P_backup; Ws=W_backup; Ss=S_backup;

% Select subset of weights or parcs
select_parcs        = {'schaefer-100'};                                     % parcellations
select_weights      = {'commit' 'qt1' 'nos' 'los'};                         % weights
[Cs,Ws]             = conn_selectData(Cs,Ws,select_weights,1);
[Cs,Ps]             = conn_selectData(Cs,Ps,select_parcs,2);

% Load parcellation info
pinfo               = conn_getParcInfo(Ps,locLabel,seg);  

%--------------------------------------------------------------------------
%% Connectome operations

% Pre-filtering
Cs                  = conn_excludeSubs(Cs,Ws,HCex);                         % Exclude subjects
Cs                  = conn_cut(Cs,Ws,seg,pinfo);                            % Trim to desired segmentation
Cs                  = conn_diagSet(Cs,0);                                   % Set value of main diagonal elements
[Cs,Ws]             = conn_sc_qt12r1(Cs,Ws);                                % Transform T1 to R1 (R1 = 1/T1)
[Cs,Ws]             = conn_fc_r2z(Cs,Ws);                                   % Fischer R-Z Transform of FC data


% Filtering 
filtstr             = 'filt-';
Cs                  = conn_sc_filter(Cs,Ws,'commit');                       % edges with commit weight <1e-12 set to 0 in all SC networks (FILTER BEFORE LOG-TRANSFORMING!!)
filtstr             = [filtstr 'commit-'];

Cs                  = conn_sc_filter(Cs,Ws,'uniform',50);                   % uniform threshold 50% filter
filtstr             = [filtstr 'Uthr50'];


% Post-filtering
Cs                  = conn_fullMat(Cs);                                     % Ensure full & symmetric matrices                             


% Manul relabeling of weights for better plot labels (Change to fit your needs)
Ws  = cellfun(@(c) regexprep(c,'qR1','R1','ignorecase'),Ws,'UniformOutput',0);
Ws  = cellfun(@(c) regexprep(c,'median-','','ignorecase'),Ws,'UniformOutput',0);
Ws  = strrep(Ws,'los', 'LoS'); Ws  = strrep(Ws,'nos', 'NoS'); Ws  = strrep(Ws,'FCz', 'FC');


%% Visualize Data

P   = 1; 
p1_plotConns2_nomod(Cs(:,P),Ws,pinfo(P),filtstr)



%% ----------------------- Network Topology ---------------------------%%
% Currently, null networks are computed within each topology script
% separately. Would be smarter to compute the null networks once, then
% derive all desired topology measures together.

% Create colormap for this section

%% Small Worldness
tmppath = [locDeriv '/NAME_FILE_' svtag '.mat'];
str     = '';                                                               % Char vector for plot titles

A       = conn_gt_smallWorld(Cs,Ws,pinfo(P),tmppath,str);                   % Compute measures (FC excluded in the function)

P_conn_gt_smallWorld([],tmppath,str,cmap)                                   % Note: need to create your colormap here


%% Rich Club
tmppath = [locDeriv '/NAME_FILE_' svtag '.mat'];
str     = '';                                                               % Char vector for plot titles

A       = conn_gt_richClub(Cs,Ws,pinfo(P),tmppath,str);                     % Compute measures: best to exclude FC, and maybe LoS too

P_conn_gt_richClub([],tmppath,str,cmap)


%% Hub-Consistency
SFN = [locDeriv '/NAME_FILE_' svtag '.mat'];
str     = '';                                                               % Char vector for plot titles

% Compute measures
A       = conn_gt_hubConsistency(Cs,Ws,pinfo(P),SFN,str);                   % Best to exclude FC and maybe LoS     

% Plot
P_conn_gt_hubConsistency(pinfo(P),[],SFN,str,cmap)



%% Additional plots using brainnet viewer (external package)
% Notes:
%   1. This package can be challenging to work with.
%   2. Requires files be created for each layer of plot i.e. edges, nodes, etc
%   3. Recommend starting with the GUI before command line use

% Setup
P       = 1;
Cst     = Cs(:,P); 
iexp    = cell2mat(cellstrfind(Ws,{'nos' 'sift2' 'commit'}));               % Skewed datasets
inrm    = setdiff(1:numel(Ws),iexp);                                        % Near-normal datasets
icommit = cellstrfind(Ws,'COMMIT');
dstr    = 'Weighted Networks';

% Plots (Skewed datasets may not plot nicely automatically here)
plot_conn_networksGroup(Cst,Ws,pinfo(P),'nz',dstr)                          % All networks

plot_conn_networksGroup(Cst(iexp),Ws(iexp),pinfo(P),'nz',dstr)              % Just networks with skewed weights

plot_conn_networksGroup(Cst(inrm),Ws(inrm),pinfo(P),'nz',dstr)              % Just networks with normally distributed weights

% Just Group COMMIT-log10
d=Cst{icommit}; d(d==0)=nan; d=log10(d); d=d+abs(min(d(:))); d(isnan(d))=0;
plot_conn_networksGroup({d},Ws(icommit),pinfo(P),'nz',dstr)


% Just networks with skewed weights: Ensuring plots are nice
% 1. Group average
% 2. Log Transform
% 3. SHIFT ALL WEIGHTS TO > 0

Cst=Cs(iexp,P); Wst=Ws(iexp);
for ii=1:numel(iexp); 
    dt              = Cst{ii}; 
    dt              = groupavg(dt,3,'nz');                                  % group averaging
    
    dt(dt==0)       = nan;
    dt              = log10(dt);                                            % log conversion
    if length(find(dt<0))>1 
        minval      = abs(min(dt(dt~=0)));                                  % Potentially not ideal that this value is dataset specific...
        dt          = dt+minval+(minval*.001);                              % shifting to > 0
    end
    dt(isnan(dt))   = 0;
    Cst{ii}         = dt;
end
plot_conn_networksGroup(Cst,Wst,pinfo(P),'nz',dstr)



%--------------------------------------------------------------------------
end




