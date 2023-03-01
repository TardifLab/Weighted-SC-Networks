function p1_plotConns1_lognorm(Cs,Ws,varargin)
% 
%  Produces the following plots:
%   1. Connectomes (group)
%   2. Histograms (subject)
%   3. Edge Groups
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
Wiexp   = cell2mat(cellstrfind(Ws,{'nos' 'sift2' 'commit'}));               % skewed datasets
Winrm   = setdiff(1:numel(Ws),Wiexp);                                       % near normally distributed datasets
cmap    = linspecer(numel(Ws));                                             % External colormap package
Ifc     = cell2mat(cellstrfind(Ws,{'FC'}));
Ilos    = cell2mat(cellstrfind(Ws,{'los'}));


%% PLOTS 1: CONNECTOMES

% All connectomes together
[~,AX]      = plot_conn_connsGroup(Cs,Ws,pinfo,'nz',str);   
AX          = unicmap(AX,'plasma');


% Individual connectomes
for Wi = 1:I
    connplot(Cs{Wi},Ws{Wi},pinfo,'',0); unicmap(gca,'kindlmann');
end

%% PLOTS 2: HISTOGRAMS

% Options
ppos = [-2559 839 1276 367];                                                % Desired size & position of plots

% Histograms: Count
plot_conn_histNorm(Cs(Wiexp),Ws(Wiexp),Ps,str,cmap(Wiexp,:),ppos)           % Overlay of exponentially distributed weights
plot_conn_histNorm(Cs(Winrm),Ws(Winrm),Ps,str,cmap(Winrm,:),ppos)           % Overlay of normally distributed weights


% Histograms: PDE
plot_conn_histSmooth(Cs(Wiexp),Ws(Wiexp),Ps,str,cmap(Wiexp,:),ppos)         % Overlay of exponentially distributed weights
plot_conn_histSmooth(Cs(Winrm),Ws(Winrm),Ps,str,cmap(Winrm,:),ppos)         % Overlay of normally distributed weights



%% PLOTS 3: EDGE GROUPS
% This section can be used to look at a variety of different metrics within
% edge groups.

plot_conn_edgeGroups(Cs,Ws,pinfo,'',[],str);


%--------------------------------------------------------------------------
end








