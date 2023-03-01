function [varargout] = plot_conn_corrAllW8s(D,W8s,varargin)
%  Computes & plots a correlation coefficient between all input 
%  connectivity datasets and a target dataset at both the subject- & 
%  group-level
%
% Input:
%                           + Required +
%   D       : IxJ cell array, each cell containing an NxNxS stack of 2D
%             square, undirected connectivity matrices (S=subs)
%   W8s     : Cell array of strings describing edge weights 
%             (required for finding Wvs)
%
%                           + Optional +
%   Wvs     : Character vector label indicating which weight to run 
%             correlations against i.e. all other weights vs this one
%             (default = 1st label)
%   Type    : Character vector indicating type of correlation
%             {'Pearson' 'Spearman'} (default=Pearson)
%   pinfo   : 1xJ structure with parcellation info (see conn_getParcInfo.m)
%   str     : character vector for plot title
%   cmap    : Wx3 double used as color map
%   avg     : String indicating method for group-averaging (see groupavg.m)
%             {'all' 'nz'}; default = 'nz' (non-zero )
%   plotON  : []/1 switch to turn on plot (1=ON; default)
%
% Output:
%   V       : Plot handle
%   R       : IxJ nested-cell array of sc-fc correlations. Each cell in the
%             IxJ array is a 1x2 cell with the following data:
%               - 1xSn double of subject-level sc-fc correlations
%               - group-level sc-fc correlation
%
% 2021 Mark C Nelson, McConnell Brain Imaging Centre, MNI, McGill
%--------------------------------------------------------------------------
[I,J]                       = size(D);
fntsz                       = 25;

%% Optional inputs
nin                                     = max(nargin,1) - 2;
defaults                                = {1,'Pearson','','unknown data',[],'nz',1};
[Iwvs,Type,pinfo,str,cmap,avg,plotON]   = INhandler(varargin,nin,defaults);

if isempty(cmap); colormap(hsv(I)); end

%% Setup
if isnumeric(Iwvs)
    Wvs         = W8s{Iwvs};
else
    Wvs         = Iwvs;
end

it          = cellstrfind(W8s, Wvs);
R           = cell(I-1,J);
ppos        = [-2559 702 921 463];

% Handle potential ambiguity
if numel(it)>1
    it      = find(cellfun(@(c) strcmpi(Wvs,c),W8s));                       % more rigid-search criteria
end
if numel(it)~=1; error('ERROR: Check your input to Wvs & your labels'); end

if exist('pinfo','var')
    parcs   = {pinfo.name};                                                 % get parcellation names
else
    parcs   = genericParcNames(J);
end

switch lower(Type)
    case 'pearson'
        strtyp          = 'Pearson''s r';
    case 'spearman'
        strtyp          = 'Spearman''s \rho';
end

%% Function
for jj = 1 : J
    sp = 0;
    wii = setdiff(1:I, it);
    for ii = wii
        sp                  = sp + 1;
        dtmp1               = D{ii,jj};
        dtmp2               = D{it,jj};                                     % re-extract target dataset each time (retained edges varies)
        
        % Get subject-level corr
        Sn                  = min([size(dtmp1,3),size(dtmp1,3)]);
        r                   = nan(1,Sn);
        for ss = 1 : Sn
            dtmp1s          = dtmp1(:,:,ss);
            dtmp2s          = dtmp2(:,:,ss);
            r(ss)           = conn_corrmat(dtmp1s,dtmp2s,Type);
        end
        
        % Get group-level corr
        dtmp1               = groupavg(dtmp1,3,avg);
        dtmp2               = groupavg(dtmp2,3,avg);
        rg                  = conn_corrmat(dtmp1,dtmp2,Type);
        R{sp,jj}            = {r,rg};
    end
end
        
%% Plot
if isempty(plotON) || plotON~=1
    return
end

% violin plot
for jj = 1 : J
    dtmp    = R(:,jj);
    rsubj   = cellfun(@(c) c{1}, dtmp, 'UniformOutput',0);
    rgrp    = cellfun(@(c) c{2}, dtmp);
    
    strplt1  = ['Subject- vs group-level correlations: ' str ': ' parcs{jj}];
    strplt2  = ['Subject- & group-level correlations: ' parcs{jj}];
    strtitl  = {['\fontsize{20}' strplt2] ; '\fontsize{12} (Colored bar marks group-level)'};
    Ylbl     = [strtyp ' vs ' upper(Wvs)]; 
    
    [~,V] = plot_violin(rsubj',W8s(wii),rgrp',Ylbl,strplt1,parcs{jj},ppos,cmap(wii,:));
    title(strtitl);
end

%% Optional output
if nargout > 0
    varargout = cell(1,nargout);
    for oo = 1 : nargout
        switch oo
            case 1; varargout{oo} = V;
            case 2; varargout{oo} = R;
        end
    end
end

%--------------------------------------------------------------------------
end