function plot_conn_edgeGroups(D,W8s,pinfo,varargin)
%  Computes & plots multiple metrics across edges grouped using multiple
%  methods
%
% * IMPORTANT * All Plots may not be desired. Code is designed to pause at
%               plotting section (~line 570) and wait for manual input. 
%
%  * Communities are defined by "pinfo"
%
%   Edge Grouping Methods:
%       - Hemisphere within/between
%       - Communities within/between
%       - Transmodal/Unimodal Within/Between
%
%   Metrics:
%       - Edge weight distributions (group & subject)
%       - connection density (group & subject)
%       - Edgewise variability in edge weights (within group & subject) 
%       - Subjectwise variability in edge weights (within edge)
%
%
% Input:
%                           + Required +
%   D       : IxJ cell array, each cell containing an NxNxS stack of 2D
%             square, undirected connectivity matrices (S=subs)
%   W8s     : 1xW Cell array of character vectors describing edge weights
%             of datasets in D
%   pinfo   : 1xJ structure with parcellation info (see conn_getParcInfo.m)
%
%
%                           + Optional +
%   target  : target dataset(s) to be correlated with all others
%             1xn cell array of strings, must match labels in W8s
%   Type    : Character vector indicating type of correlation
%             {'Pearson' 'Spearman'} (default=Pearson)
%   str     : character vector for plot title
%
% 2022 Mark C Nelson, McConnell Brain Imaging Centre, MNI, McGill
%--------------------------------------------------------------------------

%% Optional inputs
nin                 = max(nargin,1) - 3;
defaults            = {'','Pearson','unknown data'};
[target,Type,str]   = INhandler(varargin,nin,defaults);

%% Setup
TRUE            = 0;
[I,J]           = size(D);
Nmeth           = 3;                                                        % Number of edge grouping methods
Nsub            = size(D{1,1},3);                                           % Assumes all datasets have same number of subjects and subs are 3rd dim
avg             = 'nz';                                                     % {'nz' 'all'} nz excludes 0s from cross-subject averages
parcs           = {pinfo.name};                                             % parcellation names
fntsz           = 25;                                                       % font size for plots
ppos            = [1442 185 735 612];                                       % plot size & location
% [scale,ind]     = DscaleCheck(D);                                           % Can be used to determine if all datasets can be plot together

switch I
    case 1
         ppos2  = [1441 416 236 381];
    case 2
        ppos2   = [-2559 827 463 381];
    case 3
        ppos2   = [1441 340 588 457]; 
    otherwise
        ppos2   = [-2559 456 1336 752];
end

% Handle variables necessary for correlations
if ~iscell(target); target={target}; end
itrg            = cell2mat(cellstrfind(W8s, target));                       % All weights correlated with these weights
if ~isempty(itrg)
    TRUE            = 1;
    Copt            = {'Type',Type,'Rows','complete'};                      % used in correlation function calls
    W8trg           = W8s(itrg);                                            % Target weight(s)
    Icorr           = setdiff(1:I,itrg);                                    % These weights will be correlated with the target
    Nicorr          = length(Icorr);
    Ntrg            = length(W8trg);
    W8s_corr        = W8s(Icorr);
    
    % Extract target datasets
    DTRG_s          = cell(Ntrg,J);
    DTRG_g          = cell(Ntrg,J);
    for jj = 1 : J
        for tt = 1 : Ntrg
            DTRG_s{tt,jj}   = D{itrg(tt),jj};
            DTRG_g{tt,jj}   = groupavg(DTRG_s{tt,jj},3,avg);
        end
    end
    
    switch lower(Type)
        case 'pearson'
            strtyp          = 'Pearson''s r';
        case 'spearman'
            strtyp          = 'Spearman''s \rho';
    end
else
    Icorr           = [];
    Nicorr          = 1;
    Ntrg            = 1;
end

%% Setup
% Storage 
grpW8                               = cell(1,J);                            % Group level edge weights
subW8                               = cell(1,J);                            % Subject level edge weights
denstyg                             = cell(1,J);                            % Connection density group level
denstys                             = cell(1,J);                            % Connection density subject level
vrncW8_edgeg                        = cell(1,J);                            % Edge Variance within group
vrncW8_edges                        = cell(1,J);                            % Edge Variance within subject
vrncW8_subjs                        = cell(1,J);                            % Subject Variance within edge
targcorr_g                          = cell(1,J);                            % Correlations with target
targcorr_s                          = cell(1,J);                            % Correlations with target

% Loop over parcellations
for jj = 1 : J
    dparc                       = D(:,jj);                                  % all datasets for this parcellation
    N                           = length(pinfo(jj).labels);                 % node count
    iut                         = find(triu(ones(N))~=0);                   % indices for upper triangle WITH main diagonal
    
    % Masks for edge groups
    maskhemi                    = pinfo(jj).masks.hemisphere;               % 1 = within hemisphere, 2 = Between hemispheres
    maskcomm                    = pinfo(jj).masks.community;                % 1 = within community,  2 = Between communities
    maskgrad                    = pinfo(jj).masks.gradient;                 % 1 = unimodal, 2 = between, 3 = transmodal, 0 = excluded (subcortex & limbic)
    
    % Number of edges in each group
    Nedg_hemi                   = [length(find(maskhemi==1)) length(find(maskhemi==2))];
    Nedg_comm                   = [length(find(maskcomm==1)) length(find(maskcomm==2))];
    Nedg_grad                   = [length(find(maskgrad==1)) length(find(maskgrad==2)) length(find(maskgrad==3))];
  
    % Store info in these temporary cells
    vrncW8_edgegt               = cell(I,Nmeth);
    vrncW8_edgest               = cell(I,Nmeth);
    vrncW8_subjst               = cell(I,Nmeth);
    denstygt                    = cell(I,Nmeth);
    denstyst                    = cell(I,Nmeth);
    grpW8t                      = cell(I,Nmeth);
    subW8t                      = cell(I,Nmeth);
    targcorr_gt                 = cell(Nicorr,Nmeth);
    targcorr_st                 = cell(Nicorr,Nmeth);

        
    % Iterate over all weighted networks
    it=0;
    for ii = 1 : I
        ds                      = dparc{ii};                                % Stack of subject-level connectomes
        
        %% ----------------------     GROUP-LEVEL     ---------------------- %%
        
        % Prep Group data
        dg                          = groupavg(ds,3,avg);                       % Group average
        dg(isnan(dg))               = 0;                                   
        dg(iut)                     = nan;                                      % upper triangle to nans to exclude (assumes symmetry)
                
        
        %% Hemisphere within/between
        maskuse                     = maskhemi;
        IN                          = 1;
        dstore                      = cell(1,2);
        rg1                         = nan(1,Ntrg);
        rg2                         = nan(1,Ntrg);
        
        % Split data using mask
        d1                          = dg(maskuse==1);                       % within hemisphere edges group
        d2                          = dg(maskuse==2);                       % between hemispheric edges group
        d1(d1==0)                   = nan;                                  % Remaining 0s back to nans (retained for density calc)
        d2(d2==0)                   = nan;
        
        % Correlations
        if TRUE==1 && ismember(ii,Icorr)
            it=it+1;
            for tt = 1 : Ntrg
                TRG_g               = DTRG_g{tt,jj};
                TRG_d1              = TRG_g(maskuse==1);
                TRG_d2              = TRG_g(maskuse==2);
                TRG_d1(isnan(d1))   = nan;                                  % isolate non zero edges from dataset of interest in target (unnecessary with Copt)
                TRG_d2(isnan(d2))   = nan;
                
                rg1(tt)             = corr(d1,TRG_d1,Copt{:});
                rg2(tt)             = corr(d2,TRG_d2,Copt{:});
            end
            targcorr_gt{it,IN}      = {rg1 rg2};
        end
        
        % Edge weights
        dstore{1}                   = d1(~isnan(d1));
        dstore{2}                   = d2(~isnan(d2));
        grpW8t{ii,IN}               = dstore;                                   % Store all group level weights less the nans
        
        % Connection density
        dstore{1}                   = length(find(~isnan(d1)))/length(d1); 
        dstore{2}                   = length(find(~isnan(d2)))/length(d2);
        denstygt{ii,IN}             = dstore;

        % Edge weight variance: CQD
        EQ1                         = quantile(d1,.25);
        EQ3                         = quantile(d1,.75);
        dstore{1}                   = (EQ3 - EQ1) ./ (EQ3 + EQ1);
        EQ1                         = quantile(d2,.25);
        EQ3                         = quantile(d2,.75);
        dstore{2}                   = (EQ3 - EQ1) ./ (EQ3 + EQ1);
        vrncW8_edgegt{ii,IN}        = dstore;
        
        
        
        %% Community within/between
        maskuse                     = maskcomm;
        IN                          = 2;
        dstore                      = cell(1,2);
        rg1                         = nan(1,Ntrg);
        rg2                         = nan(1,Ntrg);
        
        % Split data using mask
        d1                          = dg(maskuse==1);                       % within community edges group
        d2                          = dg(maskuse==2);                       % between community edges group
        d1(d1==0)                   = nan;                                  % Remaining 0s back to nans (retained for density calc)
        d2(d2==0)                   = nan;
        
        % Correlations
        if TRUE==1 && ismember(ii,Icorr)
            for tt = 1 : Ntrg
                TRG_g               = DTRG_g{tt,jj};
                TRG_d1              = TRG_g(maskuse==1);
                TRG_d2              = TRG_g(maskuse==2);
                TRG_d1(isnan(d1))   = nan;                                  % isolate non zero edges from dataset of interest in target (unnecessary with Copt)
                TRG_d2(isnan(d2))   = nan;
                
                rg1(tt)             = corr(d1,TRG_d1,Copt{:});
                rg2(tt)             = corr(d2,TRG_d2,Copt{:});
            end
            targcorr_gt{it,IN}      = {rg1 rg2};
        end
        
        % Edge weights
        dstore{1}                   = d1(~isnan(d1));
        dstore{2}                   = d2(~isnan(d2));
        grpW8t{ii,IN}               = dstore;                                   % Store all group level weights less the nans
        
        % Connection density
        dstore{1}                   = length(find(~isnan(d1)))/length(d1); 
        dstore{2}                   = length(find(~isnan(d2)))/length(d2);
        denstygt{ii,IN}             = dstore;

        % Edge weight variance: CQD
        EQ1                         = quantile(d1,.25);
        EQ3                         = quantile(d1,.75);
        dstore{1}                   = (EQ3 - EQ1) ./ (EQ3 + EQ1);
        EQ1                         = quantile(d2,.25);
        EQ3                         = quantile(d2,.75);
        dstore{2}                   = (EQ3 - EQ1) ./ (EQ3 + EQ1);
        vrncW8_edgegt{ii,IN}        = dstore;
        
        
        
        %% Unimodal/Transmodal within/between: % 1 = unimodal, 2 = between, 3 = transmodal
        maskuse                     = maskgrad;
        IN                          = 3;
        dstore                      = cell(1,3);
        rg1                         = nan(1,Ntrg);
        rg2                         = nan(1,Ntrg);
        rg3                         = nan(1,Ntrg);
        
        % Split data using mask
        d1                          = dg(maskuse==1);                       % within UNIMODAL edges group
        d2                          = dg(maskuse==2);                       % between edges group
        d3                          = dg(maskuse==3);                       % within TRANSMODAL edges group
        d1(d1==0)                   = nan;                                  % Remaining 0s back to nans (retained for density calc)
        d2(d2==0)                   = nan;
        d3(d3==0)                   = nan;
        
        % Correlations
        if TRUE==1 && ismember(ii,Icorr)
            for tt = 1 : Ntrg
                TRG_g               = DTRG_g{tt,jj};
                TRG_d1              = TRG_g(maskuse==1);
                TRG_d2              = TRG_g(maskuse==2);
                TRG_d3              = TRG_g(maskuse==3);
                TRG_d1(isnan(d1))   = nan;                                  % isolate non zero edges from dataset of interest in target (unnecessary with Copt)
                TRG_d2(isnan(d2))   = nan;
                TRG_d3(isnan(d3))   = nan;
                
                rg1(tt)             = corr(d1,TRG_d1,Copt{:});
                rg2(tt)             = corr(d2,TRG_d2,Copt{:});
                rg3(tt)             = corr(d3,TRG_d3,Copt{:});
            end
            targcorr_gt{it,IN}      = {rg1 rg2 rg3};
        end
        
        % Edge weights
        dstore{1}               = d1(~isnan(d1));
        dstore{2}               = d2(~isnan(d2));
        dstore{3}               = d3(~isnan(d3));
        grpW8t{ii,IN}           = dstore;                                   % Store all group level weights less the nans
        
        % Connection density
        dstore{1}               = length(find(~isnan(d1)))/length(d1); 
        dstore{2}               = length(find(~isnan(d2)))/length(d2);
        dstore{3}               = length(find(~isnan(d3)))/length(d3);
        denstygt{ii,IN}         = dstore;

        % Edge weight variance: CQD
        EQ1                     = quantile(d1,.25);
        EQ3                     = quantile(d1,.75);
        dstore{1}               = (EQ3 - EQ1) ./ (EQ3 + EQ1);
        EQ1                     = quantile(d2,.25);
        EQ3                     = quantile(d2,.75);
        dstore{2}               = (EQ3 - EQ1) ./ (EQ3 + EQ1);
        EQ1                     = quantile(d3,.25);
        EQ3                     = quantile(d3,.75);
        dstore{3}               = (EQ3 - EQ1) ./ (EQ3 + EQ1);
        vrncW8_edgegt{ii,IN}    = dstore;
        


        %% ----------------     SUBJECT LEVEL     ----------------- %%
        
        
        % Prep Subject data
        ds(isnan(ds))           = 0;                                        % Ensure no nans
        for kk = 1 : Nsub
            dst                 = ds(:,:,kk);
            dst(iut)            = nan;
            ds(:,:,kk)          = dst;
        end
        
        %% Organize subject data by groups
        
        %% Hemisphere
        maskuse                         = maskhemi;
        IN                              = 1;
        stack1                          = nan(Nedg_hemi(1),Nsub);
        stack2                          = nan(Nedg_hemi(2),Nsub);
        rs1                             = nan(Nsub,Ntrg);
        rs2                             = nan(Nsub,Ntrg);
        for kk = 1 : Nsub
            dst                         = ds(:,:,kk);
            d1                          = dst(maskuse==1);                  % within hemisphere edges subject
            d2                          = dst(maskuse==2);                  % between hemispheric edges subject
            d1(d1==0)                   = nan;
            d2(d2==0)                   = nan;
            
            % Correlations
            if TRUE==1 && ismember(ii,Icorr)
                for tt = 1 : Ntrg
                    TRG_s               = DTRG_s{tt,jj}(:,:,kk);
                    TRG_d1              = TRG_s(maskuse==1);
                    TRG_d2              = TRG_s(maskuse==2);
                    TRG_d1(isnan(d1))   = nan;                              % isolate non zero edges from dataset of interest in target (unnecessary with Copt)
                    TRG_d2(isnan(d2))   = nan;

                    rs1(kk,tt)          = corr(d1,TRG_d1,Copt{:});
                    rs2(kk,tt)          = corr(d2,TRG_d2,Copt{:});
                end
            end
            
            stack1(:,kk)                = d1;                               % Stack subject edge weights
            stack2(:,kk)                = d2;
        end
        subW8t{ii,IN}                   = {stack1 stack2};
        if TRUE==1 && ismember(ii,Icorr)
            targcorr_st{it,IN}          = {rs1 rs2};
        end
        
        
        %% Communities
        maskuse                         = maskcomm;
        IN                              = 2;
        stack1                          = nan(Nedg_comm(1),Nsub);
        stack2                          = nan(Nedg_comm(2),Nsub);
        rs1                             = nan(Nsub,Ntrg);
        rs2                             = nan(Nsub,Ntrg);
        for kk = 1 : Nsub
            dst                         = ds(:,:,kk);
            d1                          = dst(maskuse==1);                  % within community edges subject
            d2                          = dst(maskuse==2);                  % between community edges subject
            d1(d1==0)                   = nan;
            d2(d2==0)                   = nan;
            
            % Correlations
            if TRUE==1 && ismember(ii,Icorr)
                for tt = 1 : Ntrg
                    TRG_s               = DTRG_s{tt,jj}(:,:,kk);
                    TRG_d1              = TRG_s(maskuse==1);
                    TRG_d2              = TRG_s(maskuse==2);
                    TRG_d1(isnan(d1))   = nan;                              % isolate non zero edges from dataset of interest in target (unnecessary with Copt)
                    TRG_d2(isnan(d2))   = nan;

                    rs1(kk,tt)          = corr(d1,TRG_d1,Copt{:});
                    rs2(kk,tt)          = corr(d2,TRG_d2,Copt{:});
                end
            end
            
            stack1(:,kk)                = d1;
            stack2(:,kk)                = d2;
        end
        subW8t{ii,IN}                   = {stack1 stack2};
        if TRUE==1 && ismember(ii,Icorr)
            targcorr_st{it,IN}          = {rs1 rs2};
        end
        
        
        %% Gradient
        maskuse                         = maskgrad;
        IN                              = 3;
        stack1                          = nan(Nedg_grad(1),Nsub);
        stack2                          = nan(Nedg_grad(2),Nsub);
        stack3                          = nan(Nedg_grad(3),Nsub);
        rs1                             = nan(Nsub,Ntrg);
        rs2                             = nan(Nsub,Ntrg);
        rs3                             = nan(Nsub,Ntrg);
        for kk = 1 : Nsub
            dst                         = ds(:,:,kk);
            d1                          = dst(maskuse==1);                  % within UNIMODAL edges subject
            d2                          = dst(maskuse==2);                  % between edges subject
            d3                          = dst(maskuse==3);                  % within TRANSMODAL edges subject
            d1(d1==0)                   = nan;
            d2(d2==0)                   = nan;
            d3(d3==0)                   = nan;
            
            % Correlations
            if TRUE==1 && ismember(ii,Icorr)
                for tt = 1 : Ntrg
                    TRG_s               = DTRG_s{tt,jj}(:,:,kk);
                    TRG_d1              = TRG_s(maskuse==1);
                    TRG_d2              = TRG_s(maskuse==2);
                    TRG_d3              = TRG_s(maskuse==3);
                    TRG_d1(isnan(d1))   = nan;                              % isolate non zero edges from dataset of interest in target (unnecessary with Copt)
                    TRG_d2(isnan(d2))   = nan;
                    TRG_d3(isnan(d3))   = nan;

                    rs1(kk,tt)          = corr(d1,TRG_d1,Copt{:});
                    rs2(kk,tt)          = corr(d2,TRG_d2,Copt{:});
                    rs3(kk,tt)          = corr(d3,TRG_d3,Copt{:});
                end
            end
            
            stack1(:,kk)                = d1;
            stack2(:,kk)                = d2;
            stack3(:,kk)                = d3;
        end
        subW8t{ii,IN}                   = {stack1 stack2 stack3};
        if TRUE==1 && ismember(ii,Icorr)
            targcorr_st{it,IN}          = {rs1 rs2 rs3};
        end
        
        %% Compute metrics
        
        %% Hemisphere
        IN                      = 1;
        stack1                  = subW8t{ii,IN}{1};
        stack2                  = subW8t{ii,IN}{2};
        dstore                  = cell(1,2);
        
        % Edge variance: CQD
        EQ1                     = quantile(stack1,.25);
        EQ3                     = quantile(stack1,.75);
        dstore{1}               = (EQ3 - EQ1) ./ (EQ3 + EQ1);
        EQ1                     = quantile(stack2,.25);
        EQ3                     = quantile(stack2,.75);
        dstore{2}               = (EQ3 - EQ1) ./ (EQ3 + EQ1);
        vrncW8_edgest{ii,IN}    = dstore;
        
        % Subject variance: CQD
        EQ1                     = quantile(stack1,.25,2);
        EQ3                     = quantile(stack1,.75,2);
        dstore{1}               = (EQ3 - EQ1) ./ (EQ3 + EQ1);
        EQ1                     = quantile(stack2,.25,2);
        EQ3                     = quantile(stack2,.75,2);
        dstore{2}               = (EQ3 - EQ1) ./ (EQ3 + EQ1);
        vrncW8_subjst{ii,IN}    = dstore;
        
        % Connection density
        stack1(~isnan(stack1))  = 1;
        stack2(~isnan(stack2))  = 1;
        sum1                    = nansum(stack1);
        sum2                    = nansum(stack2);
        dstore{1}               = sum1/Nedg_hemi(1);      
        dstore{2}               = sum2/Nedg_hemi(2);
        denstyst{ii,IN}         = dstore;
        
        
        %% Communities
        IN                      = 2;
        stack1                  = subW8t{ii,IN}{1};
        stack2                  = subW8t{ii,IN}{2};
        dstore                  = cell(1,2);
        
        % Edge variance: CQD
        EQ1                     = quantile(stack1,.25);
        EQ3                     = quantile(stack1,.75);
        dstore{1}               = (EQ3 - EQ1) ./ (EQ3 + EQ1);
        EQ1                     = quantile(stack2,.25);
        EQ3                     = quantile(stack2,.75);
        dstore{2}               = (EQ3 - EQ1) ./ (EQ3 + EQ1);
        vrncW8_edgest{ii,IN}    = dstore;
        
        % Subject variance: CQD
        EQ1                     = quantile(stack1,.25,2);
        EQ3                     = quantile(stack1,.75,2);
        dstore{1}               = (EQ3 - EQ1) ./ (EQ3 + EQ1);
        EQ1                     = quantile(stack2,.25,2);
        EQ3                     = quantile(stack2,.75,2);
        dstore{2}               = (EQ3 - EQ1) ./ (EQ3 + EQ1);
        vrncW8_subjst{ii,IN}    = dstore;
        
        % Connection density
        stack1(~isnan(stack1))  = 1;
        stack2(~isnan(stack2))  = 1;
        sum1                    = nansum(stack1);
        sum2                    = nansum(stack2);
        dstore{1}               = sum1/Nedg_comm(1);      
        dstore{2}               = sum2/Nedg_comm(2);
        denstyst{ii,IN}         = dstore;

        
        %% Gradient
        IN                      = 3;
        stack1                  = subW8t{ii,IN}{1};
        stack2                  = subW8t{ii,IN}{2};
        stack3                  = subW8t{ii,IN}{3};
        dstore                  = cell(1,3);
        
        % Edge variance: CQD
        EQ1                     = quantile(stack1,.25);
        EQ3                     = quantile(stack1,.75);
        dstore{1}               = (EQ3 - EQ1) ./ (EQ3 + EQ1);
        EQ1                     = quantile(stack2,.25);
        EQ3                     = quantile(stack2,.75);
        dstore{2}               = (EQ3 - EQ1) ./ (EQ3 + EQ1);
        EQ1                     = quantile(stack3,.25);
        EQ3                     = quantile(stack3,.75);
        dstore{3}               = (EQ3 - EQ1) ./ (EQ3 + EQ1);
        vrncW8_edgest{ii,IN}    = dstore;
        
        % Subject variance: CQD
        EQ1                     = quantile(stack1,.25,2);
        EQ3                     = quantile(stack1,.75,2);
        dstore{1}               = (EQ3 - EQ1) ./ (EQ3 + EQ1);
        EQ1                     = quantile(stack2,.25,2);
        EQ3                     = quantile(stack2,.75,2);
        dstore{2}               = (EQ3 - EQ1) ./ (EQ3 + EQ1);
        EQ1                     = quantile(stack3,.25,2);
        EQ3                     = quantile(stack3,.75,2);
        dstore{3}               = (EQ3 - EQ1) ./ (EQ3 + EQ1);
        vrncW8_subjst{ii,IN}    = dstore;
        
        % Connection density
        stack1(~isnan(stack1))  = 1;
        stack2(~isnan(stack2))  = 1;
        stack3(~isnan(stack3))  = 1;
        sum1                    = nansum(stack1);
        sum2                    = nansum(stack2);
        sum3                    = nansum(stack3);
        dstore{1}               = sum1/Nedg_grad(1);      
        dstore{2}               = sum2/Nedg_grad(2);
        dstore{3}               = sum3/Nedg_grad(3);
        denstyst{ii,IN}         = dstore;
        
    end                                                                     % over weights
    grpW8{jj}                               = grpW8t;
    subW8{jj}                               = subW8t;
    vrncW8_edgeg{jj}                        = vrncW8_edgegt;
    vrncW8_edges{jj}                        = vrncW8_edgest;
    vrncW8_subjs{jj}                        = vrncW8_subjst;
    denstyg{jj}                             = denstygt;
    denstys{jj}                             = denstyst;
    targcorr_g{jj}                          = targcorr_gt;
    targcorr_s{jj}                          = targcorr_st;
end                                                                         % over parcs

%% --------------------------- Plots ---------------------------- %%
opt = {'UniformOutput',0};

for jj = 1 : J
    
    disp(' --> Awaiting manual selection of plots in plot_conn_edgeGroups.m <-- ')
    disp('     See approx line 575 ')
    keyboard
    
    % NOTE: This 1st large section plots violins for all metrics, but a
    % better version of some plots is available below using dual yaxis
    % boxplots (~line 925)

    
    %% Hemisphere
    Ngrp            = 2;
    IN              = 1;
    Nplt            = I*Ngrp;
    dplt            = cell(1,Nplt);
    grplbls         = {'Within' 'Between'};
    methstr         = 'Hemisphere';
    
    
    % Group Edge Weights
    dt                  = grpW8{jj}(:,IN);
    dplt(1:Ngrp:Nplt)   = cellfun(@(c) c{1},dt,opt{:});
    dplt(2:Ngrp:Nplt)   = cellfun(@(c) c{2},dt,opt{:});
    
    [~,V]               = plot_violin(dplt);                                    % Violin plot
    set(gcf,'Position',ppos2);
    for vv = 1:Ngrp:Nplt
        V(vv).ViolinColor = [1,0,0];                                        % red
        V(vv+1).ViolinColor = [0,0,0];                                      % black
    end
    set(gca,'XTick',[1.5:Ngrp:(Nplt)-.5],'XTickLabels',W8s,'FontSize',fntsz);
    ylabel('Edge Weight'); title(['Group Edge Weights: ' methstr]);
    legend([V(1).ViolinPlot V(2).ViolinPlot],grplbls);
    
    
    % Subject Edge Weights
    dt                  = subW8{jj}(:,IN);
    dplt(1:Ngrp:Nplt)   = cellfun(@(c) c{1}(:),dt,opt{:});
    dplt(2:Ngrp:Nplt)   = cellfun(@(c) c{2}(:),dt,opt{:});
    
    [~,V]               = plot_violin(dplt);                                % Violin plot
    set(gcf,'Position',ppos2);
    for vv = 1:Ngrp:Nplt
        V(vv).ViolinColor   = [1,0,0];                                      % red
        V(vv+1).ViolinColor = [0,0,0];                                      % black
    end
    set(gca,'XTick',[1.5:Ngrp:(Nplt)-.5],'XTickLabels',W8s,'FontSize',fntsz);
    ylabel('Edge Weight'); title(['Subject Edge Weights: ' methstr]);
    legend([V(1).ViolinPlot V(2).ViolinPlot],grplbls);
    
    
    % Edge Variance (within group & subject)
    dt                  = vrncW8_edges{jj}(:,IN);                           % Subject level variance
    dplt(1:Ngrp:Nplt)   = cellfun(@(c) c{1},dt,opt{:});
    dplt(2:Ngrp:Nplt)   = cellfun(@(c) c{2},dt,opt{:});
    
    dtg                 = vrncW8_edgeg{jj}(:,IN);                           % Include group level variance
    dbar                = cell(1,Nplt);
    dbar(1:Ngrp:Nplt)   = cellfun(@(c) c{1},dtg,opt{:});
    dbar(2:Ngrp:Nplt)   = cellfun(@(c) c{2},dtg,opt{:});
    
    [~,V]               = plot_violin(dplt,'',dbar);                        % Violin plot
    set(gcf,'Position',ppos2);
    for vv = 1:Ngrp:Nplt
        V(vv).ViolinColor   = [1,0,0];                                      % red
        V(vv+1).ViolinColor = [0,0,0];                                      % black
    end
    set(gca,'XTick',[1.5:Ngrp:(Nplt)-.5],'XTickLabels',W8s,'FontSize',fntsz);
    ylabel('CQD'); title(['Edge Variance: ' methstr]);
    legend([V(1).ViolinPlot V(2).ViolinPlot],grplbls);
    
    
    % Subject Variance (within edge)
    dt                  = vrncW8_subjs{jj}(:,IN);                           % Subject level variance
    dplt(1:Ngrp:Nplt)   = cellfun(@(c) c{1},dt,opt{:});
    dplt(2:Ngrp:Nplt)   = cellfun(@(c) c{2},dt,opt{:});
        
    [~,V]               = plot_violin(dplt);                                % Violin plot
    set(gcf,'Position',ppos2);
    for vv = 1:Ngrp:Nplt
        V(vv).ViolinColor   = [1,0,0];                                      % red
        V(vv+1).ViolinColor = [0,0,0];                                      % black
    end
    set(gca,'XTick',[1.5:Ngrp:(Nplt)-.5],'XTickLabels',W8s,'FontSize',fntsz);
    ylabel('CQD'); title(['Subject Variance: ' methstr]);
    legend([V(1).ViolinPlot V(2).ViolinPlot],grplbls);
    
    
    % Correlations with target
    if TRUE==1
        Nplt_corr                           = Nicorr*Ngrp;
        dplt_corr_s                         = cell(1,Nplt_corr);
        dplt_corr_g                         = cell(1,Nplt_corr);
        
        dt_s                                = targcorr_s{jj}(:,IN);         % Subject level correlations
        dt_g                                = targcorr_g{jj}(:,IN);         % Group level correlations
        
        for tt = 1 : Ntrg
            dplt_corr_s(1:Ngrp:Nplt_corr)   = cellfun(@(c) c{1}(:,tt),dt_s,opt{:});
            dplt_corr_s(2:Ngrp:Nplt_corr)   = cellfun(@(c) c{2}(:,tt),dt_s,opt{:});
            
            dplt_corr_g(1:Ngrp:Nplt_corr)   = cellfun(@(c) c{1}(tt),dt_g,opt{:});
            dplt_corr_g(2:Ngrp:Nplt_corr)   = cellfun(@(c) c{2}(tt),dt_g,opt{:});
            
            % Plot
            [~,V] = plot_violin(dplt_corr_s,'',dplt_corr_g);                % Violin plot
            set(gcf,'Position',ppos2);
            for vv = 1:Ngrp:Nplt_corr
                V(vv).ViolinColor   = [1,0,0];                              % red
                V(vv+1).ViolinColor = [0,0,0];                              % black
            end
            titlstr     = ['vs ' W8trg{tt} ' (subject & group): ' methstr];
            set(gca,'XTick',[1.5:Ngrp:(Nplt_corr)-.5],'XTickLabels',W8s_corr,'FontSize',fntsz);
            ylabel(strtyp); title(titlstr);
            legend([V(1).ViolinPlot V(2).ViolinPlot],grplbls);
        end
    end
    
    %% Community
    Ngrp            = 2;
    IN              = 2;
    Nplt            = I*Ngrp;
    dplt            = cell(1,Nplt);
    grplbls         = {'Within' 'Between'};
    methstr         = 'Communities';
    
    
    % Group Edge Weights
    dt                  = grpW8{jj}(:,IN);
    dplt(1:Ngrp:Nplt)   = cellfun(@(c) c{1},dt,opt{:});
    dplt(2:Ngrp:Nplt)   = cellfun(@(c) c{2},dt,opt{:});
    
    [~,V]               = plot_violin(dplt);                                % Violin plot
    set(gcf,'Position',ppos2);
    for vv = 1:Ngrp:Nplt
        V(vv).ViolinColor = [1,0,0];                                        % red
        V(vv+1).ViolinColor = [0,0,0];                                      % black
    end
    set(gca,'XTick',[1.5:Ngrp:(Nplt)-.5],'XTickLabels',W8s,'FontSize',fntsz);
    ylabel('Edge Weight'); title(['Group Edge Weights: ' methstr]);
    legend([V(1).ViolinPlot V(2).ViolinPlot],grplbls);
    
    
    % Subject Edge Weights
    dt                  = subW8{jj}(:,IN);
    dplt(1:Ngrp:Nplt)   = cellfun(@(c) c{1}(:),dt,opt{:});
    dplt(2:Ngrp:Nplt)   = cellfun(@(c) c{2}(:),dt,opt{:});
    
    [~,V]               = plot_violin(dplt);                                % Violin plot
    set(gcf,'Position',ppos2);
    for vv = 1:Ngrp:Nplt
        V(vv).ViolinColor   = [1,0,0];                                      % red
        V(vv+1).ViolinColor = [0,0,0];                                      % black
    end
    set(gca,'XTick',[1.5:Ngrp:(Nplt)-.5],'XTickLabels',W8s,'FontSize',fntsz);
    ylabel('Edge Weight'); title(['Subject Edge Weights: ' methstr]);
    legend([V(1).ViolinPlot V(2).ViolinPlot],grplbls);
    
    
    % Edge Variance (within group & subject)
    dt                  = vrncW8_edges{jj}(:,IN);                           % Subject level variance
    dplt(1:Ngrp:Nplt)   = cellfun(@(c) c{1},dt,opt{:});
    dplt(2:Ngrp:Nplt)   = cellfun(@(c) c{2},dt,opt{:});
    
    dtg                 = vrncW8_edgeg{jj}(:,IN);                           % Include group level variance
    dbar                = cell(1,Nplt);
    dbar(1:Ngrp:Nplt)   = cellfun(@(c) c{1},dtg,opt{:});
    dbar(2:Ngrp:Nplt)   = cellfun(@(c) c{2},dtg,opt{:});
    
    [~,V]               = plot_violin(dplt,'',dbar);                        % Violin plot
    set(gcf,'Position',ppos2);
    for vv = 1:Ngrp:Nplt
        V(vv).ViolinColor   = [1,0,0];                                      % red
        V(vv+1).ViolinColor = [0,0,0];                                      % black
    end
    set(gca,'XTick',[1.5:Ngrp:(Nplt)-.5],'XTickLabels',W8s,'FontSize',fntsz);
    ylabel('CQD'); title(['Edge Variance: ' methstr]);
    legend([V(1).ViolinPlot V(2).ViolinPlot],grplbls);
    
    
    % Subject Variance (within edge)
    dt                  = vrncW8_subjs{jj}(:,IN);                           % Subject level variance
    dplt(1:Ngrp:Nplt)   = cellfun(@(c) c{1},dt,opt{:});
    dplt(2:Ngrp:Nplt)   = cellfun(@(c) c{2},dt,opt{:});
        
    [~,V]               = plot_violin(dplt);                                % Violin plot
    set(gcf,'Position',ppos2);
    for vv = 1:Ngrp:Nplt
        V(vv).ViolinColor   = [1,0,0];                                      % red
        V(vv+1).ViolinColor = [0,0,0];                                      % black
    end
    set(gca,'XTick',[1.5:Ngrp:(Nplt)-.5],'XTickLabels',W8s,'FontSize',fntsz);
    ylabel('CQD'); title(['Subject Variance: ' methstr]);
    legend([V(1).ViolinPlot V(2).ViolinPlot],grplbls);
    
    
    % Correlations with target
    if TRUE==1
        Nplt_corr                           = Nicorr*Ngrp;
        dplt_corr_s                         = cell(1,Nplt_corr);
        dplt_corr_g                         = cell(1,Nplt_corr);
        
        dt_s                                = targcorr_s{jj}(:,IN);             % Subject level correlations
        dt_g                                = targcorr_g{jj}(:,IN);             % Group level correlations
        
        for tt = 1 : Ntrg
            dplt_corr_s(1:Ngrp:Nplt_corr)   = cellfun(@(c) c{1}(:,tt),dt_s,opt{:});
            dplt_corr_s(2:Ngrp:Nplt_corr)   = cellfun(@(c) c{2}(:,tt),dt_s,opt{:});
            
            dplt_corr_g(1:Ngrp:Nplt_corr)   = cellfun(@(c) c{1}(tt),dt_g,opt{:});
            dplt_corr_g(2:Ngrp:Nplt_corr)   = cellfun(@(c) c{2}(tt),dt_g,opt{:});
            
            % Plot
            [~,V] = plot_violin(dplt_corr_s,'',dplt_corr_g);                % Violin plot
            set(gcf,'Position',ppos2);
            for vv = 1:Ngrp:Nplt_corr
                V(vv).ViolinColor   = [1,0,0];                                      % red
                V(vv+1).ViolinColor = [0,0,0];                                      % black
            end
            titlstr     = ['vs ' W8trg{tt} ' (subject & group): ' methstr];
            set(gca,'XTick',[1.5:Ngrp:(Nplt_corr)-.5],'XTickLabels',W8s_corr,'FontSize',fntsz);
            ylabel(strtyp); title(titlstr);
            legend([V(1).ViolinPlot V(2).ViolinPlot],grplbls);
        end
    end
    
    
    %% Gradient
    Ngrp            = 3;
    IN              = 3;
    Nplt            = I*Ngrp;
    dplt            = cell(1,Nplt);
    grplbls         = {'Unimodal' 'Between' 'Transmodal'};
    methstr         = 'Gradient';
    
    
    % Group Edge Weights
    dt                  = grpW8{jj}(:,IN);
    dplt(1:Ngrp:Nplt)   = cellfun(@(c) c{1},dt,opt{:});
    dplt(2:Ngrp:Nplt)   = cellfun(@(c) c{2},dt,opt{:});
    dplt(3:Ngrp:Nplt)   = cellfun(@(c) c{3},dt,opt{:});
    
    [~,V]               = plot_violin(dplt);                                % Violin plot
    set(gcf,'Position',ppos2);
    for vv = 1:Ngrp:Nplt
        V(vv).ViolinColor   = [1,0,0];                                      % red
        V(vv+1).ViolinColor = [0,1,0];                                      % green
        V(vv+2).ViolinColor = [0,0,1];                                      % blue
    end
    set(gca,'XTick',[2:Ngrp:Nplt],'XTickLabels',W8s,'FontSize',fntsz);
    ylabel('Edge Weight'); title(['Group Edge Weights: ' methstr]);
    legend([V(1).ViolinPlot V(2).ViolinPlot V(3).ViolinPlot],grplbls);
    
    
    % Subject Edge Weights
    dt                  = subW8{jj}(:,IN);
    dplt(1:Ngrp:Nplt)   = cellfun(@(c) c{1}(:),dt,opt{:});
    dplt(2:Ngrp:Nplt)   = cellfun(@(c) c{2}(:),dt,opt{:});
    dplt(3:Ngrp:Nplt)   = cellfun(@(c) c{3}(:),dt,opt{:});
    
    [~,V]               = plot_violin(dplt);                                % Violin plot
    set(gcf,'Position',ppos2);
    for vv = 1:Ngrp:Nplt
        V(vv).ViolinColor   = [1,0,0];                                      % red
        V(vv+1).ViolinColor = [0,1,0];                                      % green
        V(vv+2).ViolinColor = [0,0,1];                                      % blue
    end
    set(gca,'XTick',[2:Ngrp:Nplt],'XTickLabels',W8s,'FontSize',fntsz);
    ylabel('Edge Weight'); title(['Subject Edge Weights: ' methstr]);
    legend([V(1).ViolinPlot V(2).ViolinPlot V(3).ViolinPlot],grplbls);
    
    
    % Edge Variance (within group & subject)
    dt                  = vrncW8_edges{jj}(:,IN);                           % Subject level variance
    dplt(1:Ngrp:Nplt)   = cellfun(@(c) c{1},dt,opt{:});
    dplt(2:Ngrp:Nplt)   = cellfun(@(c) c{2},dt,opt{:});
    dplt(3:Ngrp:Nplt)   = cellfun(@(c) c{3},dt,opt{:});
    
    dtg                 = vrncW8_edgeg{jj}(:,IN);                           % Include group level variance
    dbar                = cell(1,Nplt);
    dbar(1:Ngrp:Nplt)   = cellfun(@(c) c{1},dtg,opt{:});
    dbar(2:Ngrp:Nplt)   = cellfun(@(c) c{2},dtg,opt{:});
    dbar(3:Ngrp:Nplt)   = cellfun(@(c) c{3},dtg,opt{:});
    
    [~,V]               = plot_violin(dplt,'',dbar);                        % Violin plot
    set(gcf,'Position',ppos2);
    for vv = 1:Ngrp:Nplt
        V(vv).ViolinColor   = [1,0,0];                                      % red
        V(vv+1).ViolinColor = [0,1,0];                                      % green
        V(vv+2).ViolinColor = [0,0,1];                                      % blue
    end
    set(gca,'XTick',[2:Ngrp:Nplt],'XTickLabels',W8s,'FontSize',fntsz);
    ylabel('CQD'); title(['Edge Variance: ' methstr]);
    legend([V(1).ViolinPlot V(2).ViolinPlot V(3).ViolinPlot],grplbls);
    
    
    % Subject Variance (within edge)
    dt                  = vrncW8_subjs{jj}(:,IN);                           % Subject level variance
    dplt(1:Ngrp:Nplt)   = cellfun(@(c) c{1},dt,opt{:});
    dplt(2:Ngrp:Nplt)   = cellfun(@(c) c{2},dt,opt{:});
    dplt(3:Ngrp:Nplt)   = cellfun(@(c) c{3},dt,opt{:});
        
    [~,V]               = plot_violin(dplt);                                % Violin plot
    set(gcf,'Position',ppos2);
    for vv = 1:Ngrp:Nplt
        V(vv).ViolinColor   = [1,0,0];                                      % red
        V(vv+1).ViolinColor = [0,1,0];                                      % green
        V(vv+2).ViolinColor = [0,0,1];                                      % blue
    end
    set(gca,'XTick',[2:Ngrp:Nplt],'XTickLabels',W8s,'FontSize',fntsz);
    ylabel('CQD'); title(['Subject Variance: ' methstr]);
    legend([V(1).ViolinPlot V(2).ViolinPlot V(3).ViolinPlot],grplbls);
    
    % Correlations with target
    if TRUE==1
        Nplt_corr                           = Nicorr*Ngrp;
        dplt_corr_s                         = cell(1,Nplt_corr);
        dplt_corr_g                         = cell(1,Nplt_corr);
        
        dt_s                                = targcorr_s{jj}(:,IN);             % Subject level correlations
        dt_g                                = targcorr_g{jj}(:,IN);             % Group level correlations
        
        for tt = 1 : Ntrg
            dplt_corr_s(1:Ngrp:Nplt_corr)   = cellfun(@(c) c{1}(:,tt),dt_s,opt{:});
            dplt_corr_s(2:Ngrp:Nplt_corr)   = cellfun(@(c) c{2}(:,tt),dt_s,opt{:});
            dplt_corr_s(3:Ngrp:Nplt_corr)   = cellfun(@(c) c{3}(:,tt),dt_s,opt{:});
            
            dplt_corr_g(1:Ngrp:Nplt_corr)   = cellfun(@(c) c{1}(tt),dt_g,opt{:});
            dplt_corr_g(2:Ngrp:Nplt_corr)   = cellfun(@(c) c{2}(tt),dt_g,opt{:});
            dplt_corr_g(3:Ngrp:Nplt_corr)   = cellfun(@(c) c{3}(tt),dt_g,opt{:});
            
            % Plot
            [~,V] = plot_violin(dplt_corr_s,'',dplt_corr_g);                % Violin plot
            set(gcf,'Position',ppos2);
            for vv = 1:Ngrp:Nplt_corr
                V(vv).ViolinColor   = [1,0,0];                                      % red
                V(vv+1).ViolinColor = [0,1,0];                                      % green
                V(vv+2).ViolinColor = [0,0,1];                                      % blue
            end
            titlstr     = ['vs ' W8trg{tt} ' (subject & group): ' methstr];
            set(gca,'XTick',[2:Ngrp:Nplt_corr],'XTickLabels',W8s_corr,'FontSize',fntsz);
            ylabel(strtyp); title(titlstr);
            legend([V(1).ViolinPlot V(2).ViolinPlot V(3).ViolinPlot],grplbls);
        end
    end
    
        
    disp(' --> Awaiting input in yyaxis box plotting section of plot_conn_edgeGroups.m <-- ')
    disp('     See approx line 920 ')
    keyboard
    
    %% Boxyy plots

    % Manually set split for 2 axes (leaving my code as an example)
    axind       = cell(1,2);
    axind{1}    = cell2mat(cellstrfind(W8s,{'LoS' 'R1' 'ICVF' 'FA' 'RD' 'FC'}));    % linear
    axind{2}    = cell2mat(cellstrfind(W8s,{'NoS' 'SIFT2' 'COMMIT'}));              % log
    axind{1}    = 4:9;                                                              % linear
    axind{2}    = [1 2 3];                                                          % log
    W8sr        = W8s([axind{1} axind{2}]);
    PPOS        = [-2559 839 1276 367];
    
    
    %% Group level Edge Weights
    STRTLT              = 'Group Edge Weights';
    
    % PRINCIPLE GRADIENT
    IN                  = 3;                                                % 1=hemisphere; 2=communities; 3=gradient
    Nplt                = 3;                                                % Number of plots per dataset
    DTMP                = grpW8{jj}(:,IN);
    cmapbx              = [1,0,0; 0,1,0; 0,0,1];
    binlbls             = {'Unimodal' 'Between' 'Transmodal'};
    
    dplt            = {};                                                   % Actual data to plot will be put here
    for ii = 1 : I
        dtmp           = DTMP{ii};
        for nn = 1 : length(dtmp)
            dplt       = [dplt;dtmp{nn}];                                   % append data to vector for plotting
        end
    end
    boxplotyy(dplt,Nplt,axind,cmapbx,binlbls,STRTLT)
        
    set(gcf,'Position',PPOS);
    set(gca,'FontSize',fntsz,'XTick',[2:Nplt:(I*Nplt)],'XTickLabel',W8sr,'XTickLabelRotation',20); 
    yyaxis right; ax=gca; ax.YScale='log';
    
    
    % WITHIN/BETWEEN NETWORK
    IN                  = 2;                                                % 1=within; 2=between
    Nplt                = 2;                                                % Number of plots per dataset
    DTMP                = grpW8{jj}(:,IN);
    cmapbx              = [1,0,0; 0,0,0];
    binlbls             = {'Within' 'Between'};
    
    dplt            = {};                                                   % Actual data to plot will be put here
    for ii = 1 : I
        dtmp           = DTMP{ii};
        for nn = 1 : length(dtmp)
            dplt       = [dplt;dtmp{nn}];                                   % append data to vector for plotting
        end
    end
    boxplotyy(dplt,Nplt,axind,cmapbx,binlbls,STRTLT)
        
    set(gcf,'Position',PPOS);
    set(gca,'FontSize',fntsz,'XTick',[1.5:Nplt:(I*Nplt)-.5],'XTickLabel',W8sr,'XTickLabelRotation',20); 
    yyaxis right; ax=gca; ax.YScale='log';
    
    
    
    
    %% Subject level Edge Weights
    % SUBJECT DATA IS 2D (EDGES BY SUBS)
    % NEED TO MOD BOXPLOTYY TO VECTORIZE IF 2D
    STRTLT              = 'Subject Edge Weights';
    
    % PRINCIPLE GRADIENT
    IN                  = 3;                                                % 1=hemisphere; 2=communities; 3=gradient
    Nplt                = 3;                                                % Number of plots per dataset
    DTMP                = subW8{jj}(:,IN);
    cmapbx              = [1,0,0; 0,1,0; 0,0,1];
    binlbls             = {'Unimodal' 'Between' 'Transmodal'};
    
    dplt            = {};                                                   % Actual data to plot will be put here
    for ii = 1 : I
        dtmp           = DTMP{ii};
        for nn = 1 : length(dtmp)
            dplt       = [dplt;dtmp{nn}];                                   % append data to vector for plotting
        end
    end
    boxplotyy(dplt,Nplt,axind,cmapbx,binlbls,STRTLT)
        
    set(gcf,'Position',PPOS)
    set(gca,'FontSize',fntsz,'XTick',[2:Nplt:(I*Nplt)],'XTickLabel',W8sr,'XTickLabelRotation',20); 
    yyaxis right; ax=gca; ax.YScale='log';
    
    
    % WITHIN/BETWEEN NETWORK
    IN                  = 2;                                                % 1=within; 2=between
    Nplt                = 2;                                                % Number of plots per dataset
    DTMP                = subW8{jj}(:,IN);
    cmapbx              = [1,0,0; 0,0,0];
    binlbls             = {'Within' 'Between'};
    
    dplt            = {};                                                   % Actual data to plot will be put here
    for ii = 1 : I
        dtmp           = DTMP{ii};
        for nn = 1 : length(dtmp)
            dplt       = [dplt;dtmp{nn}];                                   % append data to vector for plotting
        end
    end
    boxplotyy(dplt,Nplt,axind,cmapbx,binlbls,STRTLT)
        
    set(gcf,'Position',PPOS);
    set(gca,'FontSize',fntsz,'XTick',[1.5:Nplt:(I*Nplt)-.5],'XTickLabel',W8sr,'XTickLabelRotation',20); 
    yyaxis right; ax=gca; ax.YScale='log';

    
    
end
     
    
%--------------------------------------------------------------------------
end