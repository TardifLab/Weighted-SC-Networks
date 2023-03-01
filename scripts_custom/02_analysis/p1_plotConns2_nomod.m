function p1_plotConns2_nomod(Cs,Ws,varargin)
% Generates plots to summarize results from this section
%
%  Produces the following plots:
%   1. Variance
%   2. Edge Length Binned
%   3. All Weights vs Edge Length
%   4. All Weights vs FC
%   5. Control for Edge Length and Replot FC Correlations
%   6. 
%
% Input
%                           + Required +
%   Cs      : IxJ cell array, each cell containing an NxNxS stack of
%             connectivity matrices (S=subs)
%   Ws      : Cell array of strings describing edge weights
%
%                           + Optional +
%   pinfo   : 1xJ structure with parcellation info (see conn_getParcInfo.m)
%   str     : Info for plot title
%
% Outputs
%
% 2022 Mark C Nelson, McConnell Brain Imaging Centre, MNI, McGill
%--------------------------------------------------------------------------
%% Optional inputs
nin=max(nargin,1)-2; 
defaults        = {[] 'unknown data'};
[pinfo,str]     = INhandler(varargin,nin,defaults);


%% Setup
[I,J]                       = size(Cs);

% Generic parcellation names if none given
if ~isempty(pinfo); 
    Ps                      = {pinfo.name};
else
    Ps                      = genericParcNames(J);
end

% Organize and generate colormap for plots
Wiexp   = cell2mat(cellstrfind(Ws,{'nos' 'sift2' 'commit'}));               % exponentially distributed datasets
Winrm   = setdiff(1:numel(Ws),Wiexp);                                       % near normally distributed datasets
cmap    = linspecer(numel(Ws)); 
Ifc     = cell2mat(cellstrfind(Ws,{'FC'}));
Ilos    = cell2mat(cellstrfind(Ws,{'los'}));


%% PLOTS 1: VARIANCE

plot_conn_var(Cs,Ws,pinfo,'cqd',str,cmap);



%% PLOTS 2: EDGE-LENGTH BINNED

% Length binned comparison of 2 weights
Witmp   = cell2mat(cellstrfind(Ws,{'LoS' 'r1' 'nos'}));
Cst=Cs(Witmp); Wst=Ws(Witmp);cmapt=cmap(Witmp,:); bins=8;

plot_conn_dbins_corr(Cst,Wst,'los',bins,str,cmapt);


% Length binned comparison of > 2 weights
Witmp   = cell2mat(cellstrfind(Ws,{'los' 'nos' 'COMMIT' 'R1'}));
Cst=Cs(Witmp); Wst=Ws(Witmp); cmapt=cmap(Witmp,:); bins=8;

plot_conn_dbins_all(Cst,Wst,'los',bins,str,cmapt);



%% PLOTS 3: ALL WEIGHTS vs Edge Length

plot_conn_vsw8(Cs,Ws,'los',Ps,str,'group');                                 % heat scatter all weights (whole-brain, group-level)

plot_conn_corrAllW8s(Cs,Ws,'los','Spearman',pinfo,str,cmap)                 % violin all weights subject & group-level



%% PLOTS 4: ALL WEIGHTS vs FC

plot_conn_vsw8(Cs,Ws,'fc',Ps,str,'group');                                  % heat scatter all weights (whole-brain, group-level)

plot_conn_corrAllW8s(Cs,Ws,'fc','Spearman',pinfo,str,cmap)                  % violin all weights subject & group-level



%% PLOTS 5: CONTROL FOR EDGE LENGTH & REPLOT FC CORRELATIONS

[Resid_sub,Ws_rsd,Resid_grp]    = conn_controlW8(Cs,Ws,'los');

Witmp = cell2mat(cellstrfind(Wsr,{'NoS' 'SIFT2' 'COMMIT' 'R1' 'ICVF' 'FA' 'RD' 'fc'}));
plot_conn_corrAllW8s(Resid_sub,Ws_rsd,'fc','Spearman',pinfo,'Residual Edge Weights',cmap(Witmp,:)) % vs FC 


%% PLOTS 6: SPECIFIC WEIGHT PAIRS (Manual Input Required)

disp(' ')
disp('* * * * * * *     ATTENTION     * * * * * * *')
disp('Manual input required for plots!')
disp('See line 100 in p1_plotConns2_nomod.m')
keyboard

% * NOTE * Manual input required at lines 108 & 126

% Settings
selWs       = {'nos' 'r1'};                                                 % MANUALLY indicate your weights of interest
qWs         = {'LoS' 'FC' };                                                % Quantitative weights
Wsubset     = [qWs selWs];
Isubset     = cell2mat(cellstrfind(Ws,Wsubset));
Csubset     = Cs(Isubset);
cmapsubset  = cmap(Isubset,:);
Dsubset     = {Csubset Wsubset};


% Prep Data
fntsz=25;
los    = Isubset(1);            x    = Isubset(3);         y    = Isubset(4);
LOS    = Csubset{1};            X    = Csubset{3};         Y    = Csubset{4};         
LOS    = groupavg(LOS,3,'nz');  X    = groupavg(X,3,'nz'); Y    = groupavg(Y,3,'nz'); 
LOS    = LOS(Y~=0);             X    = X(Y~=0);            Y    = Y(Y~=0);            
losstr = Ws{los};               xstr = Ws{x};              ystr = Ws{y};             

% % Add log transform?
X    = log(X);  xstr = [Ws{x} ' (log)'];                                   % MANUALLY indicate if you want either weight log transformed
% Y    = log(Y);  ystr = [Wsr{y} ' (log)'];


% Linear regression of LoS
P = [ones(length(LOS),1), LOS];
[~,~,Yresidual] = regress(Y,P);
[~,~,Xresidual] = regress(X,P);


% Plot Heat Scatter
r           = corr(Xresidual,Yresidual,'Type','Pearson');                   % Plot the residuals together
rstr        = {sprintf('r = %.3f', r)};
myfig(str,[-2559 811 511 397]), heatscatter(Xresidual,Yresidual);           
xlabel([xstr '-residual']); ylabel([ystr '-residual']); title(rstr); set(gca,'FontSize',fntsz);
colormap(smartcmaps('plasma'));


plot_conn_connsGroup(Dsubset{:},pinfo,'nz',str)                             % All connectomes together   

plot_conn_var(Dsubset{:},'cqd',pinfo,str,cmapsubset);                       % Variability

%--------------------------------------------------------------------------
end