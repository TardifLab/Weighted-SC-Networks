function [A] = conn_gt_coreDecomp(D,W8s,pinfo,varargin)
%  Computes k-core and s-core in binary & weighted connectivity data
%
% Input:
%                           + Required +
%   D       : IxJ cell array, each cell containing an NxNxS stack of
%             connectivity matrices (S=subs)
%   W8s     : Cell array of strings describing edge weights
%   pinfo   : 1xJ structure with parcellation info (see conn_getParcInfo.m)
%
%                           + Optional +
%   str     : character vector for plot title
%
% Output:
%   A       : 1xJ structure with fields:
%
%               A.  : 
%
% 2022 Mark C Nelson, McConnell Brain Imaging Centre, MNI, McGill
%--------------------------------------------------------------------------
%% Optional inputs
nin                             = max(nargin,1) - 3;
defaults                        = {[]};
[str]                           = INhandler(varargin,nin,defaults);                                                
                                                
%% Setup
[I,J]                           = size(D);                                  % data dimensions: weights x parcellations
parcs                           = {pinfo.name};                             % parcellation names
prctilevec_s                    = 5:.01:99;                                 % percentiles of strength to test s-core
prctilevec_k                    = 5:.1:99;                                  % percentiles of degree to test k-core

% Need to exclude FC?
ifc                             = cellstrfind(W8s,'fc');
isc                             = setdiff(1:I,ifc);
I                               = numel(isc);

% Get number of subjects
Nsub                            = unique(cellfun(@(c) size(c,3), D));       % Assume 3rd dim is subjects
if numel(Nsub) > 1
    error('The number of subjects should be equivalent in each network')
end

%% Get network attributes
A                               = [];                                       % primary storage structure
A.info.weights                  = W8s(isc);
A.info.parcs                    = parcs;
A.info.details                  = str;


% Init cell arrays to store outputs of loops (g=group & s=subject)
SCORE_g_core                    = cell(I,J);                                % Group
SCORE_g_size                    = cell(I,J);
SCORE_g_strn                    = cell(I,J);

KCORE_g_core                    = cell(1,J);
KCORE_g_size                    = cell(1,J);
KCORE_g_strn                    = cell(1,J);

SCORE_s_core                    = cell(I,J,Nsub);                           % Subject
SCORE_s_size                    = cell(I,J,Nsub);
SCORE_s_strn                    = cell(I,J,Nsub);

KCORE_s_core                    = cell(1,J,Nsub);
KCORE_s_size                    = cell(1,J,Nsub);
KCORE_s_strn                    = cell(1,J,Nsub);


% Iterate over parcellations (2nd dim of D)
for jj = 1 : J
    disp(['* * * Processing parcellation scheme ' num2str(jj) ': ' parcs{jj}])
    
    % Iterate over weights (1st dim of D)
    in=1;
    for ii = isc % 1 : I
        disp(' ')
        disp(['Computing network attributes for dataset ' num2str(in) ': ' W8s{ii}])
        
        d                       = D{ii,jj};
        
        %% ----------------------- GROUP-level ------------------------- %%
        
        dg                      = groupavg(d,3,'nz');                       % group connectivity matrix
        strength                = sum(dg);                                  % node strength
        
        % Init storage variables
        scor_cor                = {};                                       % Output cores (connectivity matrices)
        scor_siz                = [];                                       % Output core size in nodes
        scor_str                = [];                                       % Input strength used to find core
        
        % Loop over percentiles of node strength
        for pp = prctilevec_s                                                 % Sample at high density
            s_prct              = prctile(strength,pp);
            [dg_s, sn_s]        = score_wu(dg,s_prct);                      % S-core Decomposition
            
            if ~isempty(sn_s) && sn_s~=0                                    % Store if core recovered
                scor_cor        = [scor_cor dg_s];
                scor_siz        = [scor_siz sn_s];
                scor_str        = [scor_str s_prct];
            end
        end
        
        % Filter duplicate cores to minimize file size
        unique_cor              = unique(scor_siz);
        for cc = 1 : numel(unique_cor)
            trg_cor             = unique_cor(cc);
            Itrg_cor            = find(scor_siz==trg_cor);
            if numel(Itrg_cor) > 1
                Irm_cor         = Itrg_cor(2:end);                          % Just keep 1st instance of each core
            else
                continue
            end
            scor_cor(Irm_cor)   = [];
            scor_siz(Irm_cor)   = [];
            scor_str(Irm_cor)   = [];
        end
        
        % Store
        SCORE_g_core{in,jj}     = scor_cor;
        SCORE_g_size{in,jj}     = scor_siz;
        SCORE_g_strn{in,jj}     = scor_str;                                 disp('... Group level completed')
        
        
        %% --------------------- SUBJECT-level ---------------------- %%
        
        for kk = 1 : Nsub
            tic;
            dt                      = d(:,:,kk);                            % subject connectivity matrix
            strength                = sum(dt);                              % node strength
            
            % Init storage variables
            scor_cor                = {};
            scor_siz                = [];
            scor_str                = [];
            
            % Loop over percentiles of node strength
            for pp = prctilevec_s                                             % Sample at high density
                s_prct              = prctile(strength,pp);
                [dt_s, sn_s]        = score_wu(dt,s_prct);                  % S-core Decomposition

                if ~isempty(sn_s) && sn_s~=0                                % Store if core recovered
                    scor_cor        = [scor_cor dt_s];
                    scor_siz        = [scor_siz sn_s];
                    scor_str        = [scor_str s_prct];
                end
            end
            
            % Filter duplicate cores to minimize file size
            unique_cor              = unique(scor_siz);
            for cc = 1 : numel(unique_cor)
                trg_cor             = unique_cor(cc);
                Itrg_cor            = find(scor_siz==trg_cor);
                if numel(Itrg_cor) > 1
                    Irm_cor         = Itrg_cor(2:end);                      % Just keep 1st instance of each core
                else
                    continue
                end
                scor_cor(Irm_cor)   = [];
                scor_siz(Irm_cor)   = [];
                scor_str(Irm_cor)   = [];
            end
            
            % Store
            SCORE_s_core{in,jj,kk}  = scor_cor;
            SCORE_s_size{in,jj,kk}  = scor_siz;
            SCORE_s_strn{in,jj,kk}  = scor_str;
            
            % Printout status
            time=toc; 
            disp([W8s{ii} ': Subject # ' num2str(kk) ' / ' num2str(Nsub) ' completed in ' num2str(time) 's'])
        end
        disp('... Subject level completed')
                     
        
        %% store & save weighted
        
        % Store everything in structure (store & save on each iteration)
        A.grp.wei.core              = SCORE_g_core;
        A.grp.wei.size              = SCORE_g_size;
        A.grp.wei.strn              = SCORE_g_strn;
        
        A.sub.wei.core              = SCORE_s_core;
        A.sub.wei.size              = SCORE_s_size;
        A.sub.wei.strn              = SCORE_s_strn;
        

        % save on each pass (quick solution to mitigate potential crashing)
        save('~/Desktop/ConnectomeProj/myelinWeightedConnectome/0_SCFingerprinting/data/5_networkAtts_coreDecomp_schaefer400.mat','A')

        in=in+1;
    end
    
    
    %% ----------------------- Binary Group-level ----------------------- %%
    ii                          = isc(1);                                   % Use 1st structural dataset
    disp(' ')
    disp(['Computing network attributes for group level ' W8s{ii} ' BINARY'])

    d                           = D{ii,jj};
    dg                          = groupavg(d,3,'nz');                       % group connectivity matrix
    dgb                         = double(dg~=0);                            % binarize
    degree                      = sum(dgb);                                 % node degree

    % Init storage variables
    kcor_cor                    = {};                                       % Output cores (connectivity matrices)
    kcor_siz                    = [];                                       % Output core size in nodes
    kcor_deg                    = [];                                       % Input degree used to find core

    % Loop over percentiles of node degree
    for pp = prctilevec_k                                                     % Sample at high density
        k_prct                  = prctile(degree,pp);
        [dgb_k, sn_k]           = kcore_bu(dgb,k_prct);                     % K-core Decomposition

        if ~isempty(sn_k) && sn_k~=0                                        % Store if core recovered
            kcor_cor            = [kcor_cor dgb_k];
            kcor_siz            = [kcor_siz sn_k];
            kcor_deg            = [kcor_deg k_prct];
        end
    end

    % Filter duplicate cores to minimize file size
    unique_cor                  = unique(kcor_siz);
    for cc = 1 : numel(unique_cor)
        trg_cor                 = unique_cor(cc);
        Itrg_cor                = find(kcor_siz==trg_cor);
        if numel(Itrg_cor) > 1
            Irm_cor             = Itrg_cor(2:end);                          % Just keep 1st instance of each core
        else
            continue
        end
        kcor_cor(Irm_cor)       = [];
        kcor_siz(Irm_cor)       = [];
        kcor_deg(Irm_cor)       = [];
    end

    % Store
    KCORE_g_core{1,jj}          = kcor_cor;
    KCORE_g_size{1,jj}          = kcor_siz;
    KCORE_g_strn{1,jj}          = kcor_deg;                                 disp('... Binary Group level completed')
            
            
    %% ----------------------- Binary Subject-level ----------------------- %%
    
    disp(' ')
    disp('Running SUBJECT level Binary ')
    
    for kk = 1 : Nsub
        tic;
        dt                      = d(:,:,kk);                                % subject connectivity matrix
        dtb                     = double(dt~=0);                            % binarize
        degree                  = sum(dtb);                                 % node degree
        
        % Init storage variables
        kcor_cor                = {};                                       % Output cores (connectivity matrices)
        kcor_siz                = [];                                       % Output core size in nodes
        kcor_deg                = [];                                       % Input degree used to find core

        % Loop over percentiles of node degree
        for pp = prctilevec_k                                                 % Sample at high density
            k_prct              = prctile(degree,pp);
            [dtb_k, sn_k]       = kcore_bu(dtb,k_prct);                     % K-core Decomposition

            if ~isempty(sn_k) && sn_k~=0                                    % Store if core recovered
                kcor_cor        = [kcor_cor dtb_k];
                kcor_siz        = [kcor_siz sn_k];
                kcor_deg        = [kcor_deg k_prct];
            end
        end

        % Filter duplicate cores to minimize file size
        unique_cor                  = unique(kcor_siz);
        for cc = 1 : numel(unique_cor)
            trg_cor                 = unique_cor(cc);
            Itrg_cor                = find(kcor_siz==trg_cor);
            if numel(Itrg_cor) > 1
                Irm_cor             = Itrg_cor(2:end);                      % Just keep 1st instance of each core
            else
                continue
            end
            kcor_cor(Irm_cor)       = [];
            kcor_siz(Irm_cor)       = [];
            kcor_deg(Irm_cor)       = [];
        end

        % Store
        KCORE_s_core{1,jj,kk}       = kcor_cor;
        KCORE_s_size{1,jj,kk}       = kcor_siz;
        KCORE_s_strn{1,jj,kk}       = kcor_deg;

        % Printout status
        time=toc; 
        disp(['Binary: Subject # ' num2str(kk) ' / ' num2str(Nsub) ' completed in ' num2str(time) 's'])
    end
    disp('... Binary Subject level completed')

    
    
    %% TEST PLOT
% % %     ppos=[-2558 617 2559 591];
% % %     
% % %     % Compute density of subgraph across range of k
% % %     sgdensity=zeros(1,kmax);
% % %     for kk = 1 : kmax
% % %         sgdensity(kk) = Ek(kk) / (Nk(kk)* (Nk(kk)-1));
% % %     end
% % %     
% % %     % Plot Empirical only
% % %     myfig('Empirical Binary Group Structural Network',ppos); 
% % %     subplot(2,2,1); plot(1:kmax,RCCt_gb,'LineWidth',3); xlabel('k level'); ylabel('phi'); title('Rich Club Coeff'); set(gca,'FontSize',20);
% % %     subplot(2,2,2); plot(1:kmax,NDkt_gb,'LineWidth',3); xlabel('k level'); title('Remaining Nodes'); set(gca,'FontSize',20);
% % %     subplot(2,2,3); plot(1:kmax,EGkt_gb,'LineWidth',3); xlabel('k level'); title('Remaining Edges'); set(gca,'FontSize',20);
% % %     subplot(2,2,4); plot(1:kmax,sgdensity,'LineWidth',3); xlabel('k level'); title('Subgraph Density'); set(gca,'FontSize',20);
% % %     
% % %     
% % %     % Compute mean across nulls
% % %     RC_null_mean    = mean(RCCt_nb);
% % %     ND_null_mean    = mean(NDkt_nb);
% % %     EG_null_mean    = mean(EGkt_nb);
% % %     
% % %     % Compute density of NULL subgraph across range of k
% % %     sgdensity_null=zeros(1,kmax);
% % %     for kk = 1 : kmax
% % %         sgdensity_null(kk) = EG_null_mean(kk) / (ND_null_mean(kk)* (ND_null_mean(kk)-1));
% % %     end
% % %     
% % %     % Normalize empirical RC
% % %     RC_norm         = RCCt_gb ./ RC_null_mean;
% % %     
% % %     
% % %     % Plot Null & Empirical
% % %     n=2;m=3;
% % %     myfig('Null & Normalized Empirical Binary Group Structural Network',ppos); 
% % %     subplot(n,m,1); plot(1:kmax,RCCt_gb,'-b','LineWidth',3); hold on; plot(1:kmax,RC_null_mean,'-r','LineWidth',3); 
% % %                     xlabel('k level'); ylabel('phi'); title('Rich Club Coeff'); set(gca,'FontSize',20); legend({'Empr' 'Null'});
% % %                     
% % %     subplot(n,m,2); plot(1:kmax,NDkt_gb,'-b','LineWidth',3); hold on; plot(1:kmax,ND_null_mean,'-r','LineWidth',3);
% % %                     xlabel('k level'); title('Remaining Nodes'); set(gca,'FontSize',20); legend({'Empr' 'Null'});
% % %                     
% % %     subplot(n,m,3); plot(1:kmax,EGkt_gb,'-b','LineWidth',3); hold on; plot(1:kmax,EG_null_mean,'-r','LineWidth',3);
% % %                     xlabel('k level'); title('Remaining Edges'); set(gca,'FontSize',20); legend({'Empr' 'Null'});
% % %     
% % %     subplot(n,m,4); plot(1:kmax,sgdensity,'-b','LineWidth',3); hold on; plot(1:kmax,sgdensity_null,'-r','LineWidth',3);
% % %                     xlabel('k level'); title('Subgraph Density'); set(gca,'FontSize',20); legend({'Empr' 'Null'});
% % %                     
% % %     subplot(n,m,5); plot(1:kmax,RC_norm,'-b','LineWidth',3);
% % %                     xlabel('k level');  ylabel('phi norm'); title('Normalized RCC'); set(gca,'FontSize',20);
    

    %% store outputs   

    % Store everything in structure (store & save on each iteration)
    A.grp.bin.core              = KCORE_g_core;
    A.grp.bin.size              = KCORE_g_size;
    A.grp.bin.strn              = KCORE_g_strn;

    A.sub.bin.core              = KCORE_s_core;
    A.sub.bin.size              = KCORE_s_size;
    A.sub.bin.strn              = KCORE_s_strn;


    % save on each pass (quick solution to mitigate potential crashing)
    save('~/Desktop/ConnectomeProj/myelinWeightedConnectome/0_SCFingerprinting/data/5_networkAtts_coreDecomp_schaefer400.mat','A')
    
    
end
disp('* Network analysis complete! *')
keyboard

%--------------------------------------------------------------------------
end