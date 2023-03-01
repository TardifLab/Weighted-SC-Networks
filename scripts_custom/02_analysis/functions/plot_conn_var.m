function varargout = plot_conn_var(D,W8s,pinfo,varargin)
% Computes & plots the indicated measure of varability for all connectivity
% datasets in 2 forms:
%   1. Intra-subject: variability measure computed across edges
%      within each subject separately. Subjectwise distribution plotted.
%   2. Inter-subject: The inverse of 1 i.e. cross-subject variability
%      computed at each edge and the edge distribution is plotted.
%
% Input
%                           + Required +
%   D       : IxJ cell array, each cell containing an NxNxS stack of
%             connectivity matrices (S=subs)
%   W8s     : Cell array of strings describing edge weights
%   pinfo   : 1xJ structure with parcellation info (see conn_getParcInfo.m)
%
%                           + Optional +
%   metric  : character vector indicating desired variability metric
%             Options: {'cova' 'madmean' 'madmedian' 'cqd' 'iqr'}
%             (default = 'cova')
%   str     : Info for plot title
%   cmap    : Ix3 double used as color map
%   ppos    : plot size & position
%
% Outputs
%                           + Optional +
%   Dout    : 1X2 nested cell array containing Evar & Svar
%
% 2022 Mark C Nelson, McConnell Brain Imaging Centre, MNI, McGill
%--------------------------------------------------------------------------

[I,J]                       = size(D);

%% Optional inputs
nin                         = max(nargin,1) - 3;
defaults                    = {'cova','',colormap(hsv(I)),[-2559 507 1037 701]};
[metric,str,cmap,ppos]      = INhandler(varargin,nin,defaults);


%% Setup
fntsz                       = 25;
Nnode                       = unique(cellfun(@(c) size(c,1),D));

if numel(Nnode) > J
    error('All weighted datasets should have the same number of nodes within each parcellation')
end

if ~isempty(pinfo); 
    parcs                   = {pinfo.name};
else
    parcs                   = genericParcNames(J);
end

% Set metric specific labels
switch lower(metric)
    case 'cova'        
        metricstrS = 'CoVa';
%         metricstrL = 'Coefficient of Variation'; 
    case 'madmean'
        metricstrS = 'MAD-Mean';
%         metricstrL = 'Mean Absolute Deviation'; 
    case 'madmedian'
        metricstrS = 'MAD-Median';
%         metricstrL = 'Median Absolute Deviation';
    case 'cqd'
        metricstrS = 'CQD';
%         metricstrL = 'Quartile Coefficient of Dispersion';
    case 'iqr'
        metricstrS = 'IQR';
%         metricstrL = 'Interquartile Range';
    otherwise; error('Unable to use your input to METRIC')
end

%% Function
% Loop over datasets
Evar                = cell(I,J);                                            % edge variance
Svar                = cell(I,J);                                            % subject variance
for jj = 1 : J
    for ii = 1 : I
        
        dt              = D{ii,jj};
        
        % Prep data
        if length(size(dt)) == 3
            
            % This is edge data (connectomes)
            Nsub            = size(dt,3);                                   % Number of subjects
            ilt             = find(tril(ones(Nnode(jj)),-1)~=0);            % Indices for lower triangle LESS main diagonal
            Nedge           = length(ilt);                                  % Number of unique edges (assuming symmetry and 0 main diagonal!)

            % Stack all vectorized subject-level connectomes (assumes symmetry & a 0 main diagonal!)
            stack           = zeros(Nedge,Nsub);
            for ss = 1 : Nsub
                dst         = dt(:,:,ss);
                stack(:,ss) = dst(ilt);                                 
            end
            stack(stack==0) = nan;                                          % 0s to nans for exclusion from statistics
            Dstr            = 'Edge';
            
        elseif length(size(dt)) == 2
            
            % This is node data (assume it is already NODE X SUBJ)
            stack           = dt;
            stack(stack==0) = nan;                                          % 0s to nans for exclusion from statistics
            Dstr            = 'Node';
            
        else
            keyboard
            
            
        end
        
        %% Compute desired variance measure
        switch lower(metric)
            
            % Coefficient of Variation (std/mu)
            case 'cova'                                                     
                Estat   = nanstd(stack)      ./ nanmean(stack);             % Edge-Variability
                Sstat   = nanstd(stack,[],2)' ./ nanmean(stack,2)';         % Subject-Variability
                
            % Mean Absolute Deviation
            case 'madmean'                                                   
                Estat   = mad(stack,0,1);
                Sstat   = mad(stack,0,2)';
                
            % Median Absolute Deviation
            case 'madmedian'                                                
                Estat   = mad(stack,1,1);
                Sstat   = mad(stack,1,2)';
                
            % Quartile Coefficient of Dispersion (Q3 - Q1)/(Q3 + Q1)
            case 'cqd'           
                EQ1     = quantile(stack,.25,1);
                EQ3     = quantile(stack,.75,1);
                Estat   = (EQ3 - EQ1) ./ (EQ3 + EQ1);
                
                SQ1     = quantile(stack,.25,2);
                SQ3     = quantile(stack,.75,2);
                Sstat   = (SQ3 - SQ1) ./ (SQ3 + SQ1);
                Sstat   = Sstat';
                
            % Interquartile range
            case 'iqr'                                                      
                Estat   = iqr(stack,1);
                Sstat   = iqr(stack,2)';
        end
        
        Evar{ii,jj}     = Estat; clear Estat
        Svar{ii,jj}     = Sstat; clear Sstat
    end
end


%% VIOLIN plots of variance metrics
for jj = 1 : J
    
    % Intra-subject variability
    Estatall    = Evar(:,jj);
    [FIG,AX]    = plot_violin(Estatall',W8s,[],metricstrS,'Intra-Subject',parcs{jj},ppos,cmap);
    set(gca,'FontSize',fntsz); set(gcf,'Name',str);   
    
    % Inter-Subject variability
    Sstatall    = Svar(:,jj);
    [FIG2,AX2]  = plot_violin(Sstatall',W8s,[],metricstrS,'Inter-Subject',parcs{jj},ppos,cmap);
    set(gca,'FontSize',fntsz); set(gcf,'Name',str);
end


%% Surface plots of inter-subject variability

for jj = 1 : J
        
    % Subject variability
    Sstatall        = Svar(:,jj);
    Dplt            = cell(1,I);
    titlestr        = ['Regional Mean Inter-Subject ' metricstrS ': ' str ', ' parcs{jj}];
    
    if strncmpi(Dstr,'Edge',4)
        for ii = 1 : I
            Sstatt      = Sstatall{ii};
            Dt          = zeros(Nnode(jj));
            Dt(ilt)     = Sstatt;                                               % Project variance back onto connectome
            Dt          = Dt + Dt';
            Dt(Dt==0)   = nan;                                                  % Exclude 0s from stat
            Dtstat      = nanmean(Dt,2);                                        % Mean over edges
            Dplt{ii}    = Dtstat;
        end
    else
        for ii = 1 : I
            Dtstat          = Sstatall{ii};
            Dplt{ii}        = Dtstat;
        end
        titlestr        = ['Regional Inter-Subject ' metricstrS ': ' str ', ' parcs{jj}];
    end
    
    % Plot on surface
    if strncmpi(pinfo(jj).seg,'cor',3)
        Hcx         = plot_conn_surf(Dplt,pinfo(jj),'cortex',W8s, titlestr);
    else
        [Hcx,Hsx]   = plot_conn_surf(Dplt,pinfo(jj),'both',W8s, titlestr);
    end
    
    
    % Group all handles together
    AXall           = [];
    CBall           = [];
    for aa = 1 : length(Hcx);
        AXall       = [AXall; Hcx{aa}.axes(:)];                             % cortical handles
        CBall       = [CBall; Hcx{aa}.cb(:)];
    end
    
    if exist('Hsx','var')
        for aa = 1 : length(Hsx)
            AXall   = [AXall; Hsx{aa}.axes(:)];                             % subcortical handles
            CBall   = [CBall; Hsx{aa}.cb(:)];
        end
    end
    
    AXall           = unicb(AXall);                                         % Standardize all color axes
    caxismod        = caxis;
    AXall           = unicmap(AXall,'plasma');                              % Change colormaps
    
    % Fix colorbars
    for cc = 1 : length(CBall)
        CBall(cc).FontSize      = fntsz;
        CBall(cc).FontWeight    ='bold';
        CBall(cc).Ticks         = caxismod;
    end  
end



%% Optional output
if nargout > 0
    varargout = cell(1,nargout);
    for oo = 1 : nargout
        switch oo
            case 1
                Dout            = {Evar, Svar};
                varargout{oo}   = Dout;
        end
    end
end


%--------------------------------------------------------------------------
end