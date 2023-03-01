function [A] = conn_gt_hubConsistency(D,W8s,pinfo,savefn,varargin)
%  Computes several node centrality measures on weighted & binary
%  connectivity data.
%
% Input:
%                           + Required +
%   D           : IxJ cell array, each cell containing an NxNxS stack of
%                 connectivity matrices (S=subs)
%   W8s         : Cell array of strings describing edge weights
%   pinfo       : 1xJ structure with parcellation info (see conn_getParcInfo.m)
%   savefn      : full path for saving outputs
%
%                           + Optional +
%   str         : character vector for plot title
%
% Output:
%   A       : 1xJ structure
%
% 2022 Mark C Nelson, McConnell Brain Imaging Centre, MNI, McGill
%--------------------------------------------------------------------------
%% Optional inputs
nin                             = max(nargin,1) - 4;
defaults                        = {[]};
[str]                           = INhandler(varargin,nin,defaults);                                                
                                                
%% Setup
[I,J]                           = size(D);                                  % data dimensions: weights x parcellations
parcs                           = {pinfo.name};                             % parcellation names
CIs                             = {pinfo.cirois};                           % community assignments
xfm                             = 'log';                                    % Wij --> Lij transform
Nnode                           = unique(cellfun(@(c) size(c,1), D),'rows');% Node Counts (assumes columns of D are parcellations)
Nsub                            = unique(cellfun(@(c) size(c,3), D));       % Assume 3rd dim is subjects

% Quick Checks
if size(Nnode,1) > 1; error('Network NODE Counts should be equivalent across weights within each parcellation'); end
if numel(Nsub) > 1;   error('The number of subjects should be equivalent in each network'); end

% Need to exclude FC?
ifc                                 = cellstrfind(W8s,'fc');
isc                                 = setdiff(1:I,ifc);
I                                   = numel(isc);


%% Get network attributes
A                                   = [];                                   % primary storage structure
A.info.weights                      = W8s(isc);
A.info.parcs                        = parcs;
A.info.details                      = str;
A.info.xfm                          = xfm;

% Test save file path before running analysis
disp(['Testing Save Path: ' savefn])
save(savefn,'A')


% Init cell arrays to store outputs of loops (g=group & s=subject)
STRENGTH_g                          = cell(I,J);                            % Group
CLOSNESS_g                          = cell(I,J);
BETWNESS_g                          = cell(I,J);
EIGNVCTR_g                          = cell(I,J);
MDDGZSCR_g                          = cell(I,J);
PRTCPACN_g                          = cell(I,J);
PAGERANK_g                          = cell(I,J);
CLUSTRNG_g                          = cell(I,J);

STRENGTH_s                          = cell(I,J);                            % Subject
CLOSNESS_s                          = cell(I,J);
BETWNESS_s                          = cell(I,J);
EIGNVCTR_s                          = cell(I,J);
MDDGZSCR_s                          = cell(I,J);
PRTCPACN_s                          = cell(I,J);
PAGERANK_s                          = cell(I,J);
CLUSTRNG_s                          = cell(I,J);


% Iterate over parcellations (2nd dim of D)
for jj = 1 : J
    disp(['* * * Processing parcellation scheme ' num2str(jj) ': ' parcs{jj}])
    CIt                             = CIs{jj};                              % Community Assignments for this Parcellation
    Nnodet                          = Nnode(jj);
    
    % Iterate over weights (1st dim of D)
    in=1;
    for ii = isc % 1 : I
        disp(' ')
        disp(['Computing network attributes for dataset ' num2str(in) ': ' W8s{ii}])
        
        d                           = D{ii,jj};                             % Stack of subject connectivity matrices
        
        %% ----------------------- GROUP-level ------------------------- %%
        
        dg                          = groupavg(d,3,'nz');                   % group connectivity matrix
        
        % Compute Length & Distance Networks
        switch xfm
            case 'log'
                if any(dg(:)<0) || any(dg(:)>1)
                    error('The -log transform requires edge weights to be in the interval [0,1)')
                else
                    Leng            = -log(dg);
                end
            case 'inv'
                Leng                = 1./dg;
            otherwise
                error('Variable XFM should be set to either LOG or INV')
        end
        [Dist,~,~]                  = distance_wei_floyd(dg,xfm);           % distance matrix (shortest weighted path)
        
                
        % Compute Centrality Measures
        STRENGTH_g_t(1,:)           = sum(dg);                              % node strength
        CLOSNESS_g_t(1,:)           = 1./mean(Dist);                        % Closeness Centrality
        BETWNESS_g_t(1,:)           = betweenness_wei(Leng);                % node betweenness centrality
        EIGNVCTR_g_t(1,:)           = eigenvector_centrality_und(dg);       % Eigenvector Centrality
        PAGERANK_g_t(1,:)           = pagerank_centrality(dg,.85);          % PageRank Centrality

        % Module-Centrality Measures
        MDDGZSCR_g_t(1,:)           = module_degree_zscore(dg,CIt,0);       % within-module degree z-score (0=undirected)
        PRTCPACN_g_t(1,:)           = participation_coef(dg,CIt,0);         % intermodular participation coeff (0=undirected)
        
        % Clustering
        CLUSTRNG_g_t(1,:)           = clustering_coef_wu(dg);               % clustering coeff
        
        % Store
        STRENGTH_g{in,jj}           = STRENGTH_g_t;
        CLOSNESS_g{in,jj}           = CLOSNESS_g_t;
        BETWNESS_g{in,jj}           = BETWNESS_g_t;
        EIGNVCTR_g{in,jj}           = EIGNVCTR_g_t;
        MDDGZSCR_g{in,jj}           = MDDGZSCR_g_t;
        PRTCPACN_g{in,jj}           = PRTCPACN_g_t;
        PAGERANK_g{in,jj}           = PAGERANK_g_t;
        CLUSTRNG_g{in,jj}           = CLUSTRNG_g_t;                         disp('... Group level completed')
        
        
        %% --------------------- SUBJECT-level ---------------------- %%
        
        % Init storage variables
        STRENGTH_s_t                = zeros(Nsub,Nnodet);
        CLOSNESS_s_t                = zeros(Nsub,Nnodet);
        BETWNESS_s_t                = zeros(Nsub,Nnodet);
        EIGNVCTR_s_t                = zeros(Nsub,Nnodet);
        PAGERANK_s_t                = zeros(Nsub,Nnodet);
        MDDGZSCR_s_t                = zeros(Nsub,Nnodet);
        PRTCPACN_s_t                = zeros(Nsub,Nnodet);
        CLUSTRNG_s_t                = zeros(Nsub,Nnodet);
        
        % Loop over subjects
        for kk = 1 : Nsub
            tic;
            dt                      = d(:,:,kk);                            % subject connectivity matrix
            
            % Compute Length & Distance Networks
            switch xfm
                case 'log'
                    if any(dt(:)<0) || any(dt(:)>1)
                        error('The -log transform requires edge weights to be in the interval [0,1)')
                    else
                        Leng        = -log(dt);
                    end
                case 'inv'
                    Leng            = 1./dt;
                otherwise
                    error('Variable XFM should be set to either LOG or INV')
            end
            [Dist,~,~]              = distance_wei_floyd(dt,xfm);           % distance matrix (shortest weighted path)

            % Compute Centrality Measures
            STRENGTH_s_t(kk,:)      = sum(dt);                              % node strength
            CLOSNESS_s_t(kk,:)      = 1./mean(Dist);                        % Closeness Centrality
            BETWNESS_s_t(kk,:)      = betweenness_wei(Leng);                % node betweenness centrality
            EIGNVCTR_s_t(kk,:)      = eigenvector_centrality_und(dt);       % Eigenvector Centrality
            PAGERANK_s_t(kk,:)      = pagerank_centrality(dt,.85);          % PageRank Centrality

            % Module-Centrality Measures
            MDDGZSCR_s_t(kk,:)      = module_degree_zscore(dt,CIt,0);       % within-module degree z-score (0=undirected)
            PRTCPACN_s_t(kk,:)      = participation_coef(dt,CIt,0);         % intermodular participation coeff (0=undirected)
            
            % Clustering
            CLUSTRNG_s_t(1,:)       = clustering_coef_wu(dt);               % clustering coeff

            % Printout status
            time=toc; 
            disp([W8s{ii} ': Subject # ' num2str(kk) ' / ' num2str(Nsub) ' completed in ' num2str(time) 's'])
        end
        
        % Store
        STRENGTH_s{in,jj}           = STRENGTH_s_t;
        CLOSNESS_s{in,jj}           = CLOSNESS_s_t;
        BETWNESS_s{in,jj}           = BETWNESS_s_t;
        EIGNVCTR_s{in,jj}           = EIGNVCTR_s_t;
        MDDGZSCR_s{in,jj}           = MDDGZSCR_s_t;
        PRTCPACN_s{in,jj}           = PRTCPACN_s_t;
        PAGERANK_s{in,jj}           = PAGERANK_s_t;
        CLUSTRNG_s{in,jj}           = CLUSTRNG_s_t;                         disp('... Subject level completed')
                     
        
        %% store & save weighted
        
        % Store everything in structure (store & save on each iteration)
        A.grp.wei.strength          = STRENGTH_g;
        A.grp.wei.closeness         = CLOSNESS_g;
        A.grp.wei.betweenness       = BETWNESS_g;
        A.grp.wei.eigenvector       = EIGNVCTR_g;
        A.grp.wei.moddegzscore      = MDDGZSCR_g;
        A.grp.wei.participation     = PRTCPACN_g;
        A.grp.wei.pagerank          = PAGERANK_g;
        A.grp.wei.clustering        = CLUSTRNG_g;
        
        A.sub.wei.strength          = STRENGTH_s;
        A.sub.wei.closeness         = CLOSNESS_s;
        A.sub.wei.betweenness       = BETWNESS_s;
        A.sub.wei.eigenvector       = EIGNVCTR_s;
        A.sub.wei.moddegzscore      = MDDGZSCR_s;
        A.sub.wei.participation     = PRTCPACN_s;
        A.sub.wei.pagerank          = PAGERANK_s;
        A.sub.wei.clustering        = CLUSTRNG_s;
        

        % save on each pass (quick solution to mitigate potential crashing)
        save(savefn,'A')

        in=in+1;
    end
    
    
    %% ----------------------- Binary Group-level ----------------------- %%
    ii                          = isc(1);                                   % Use 1st structural dataset
    disp(' ')
    disp(['Computing network attributes for group level ' W8s{ii} ' BINARY'])

    d                           = D{ii,jj};
    dg                          = groupavg(d,3,'nz');                       % group connectivity matrix
    dgb                         = double(dg~=0);                            % binarize
    Dist                        = distance_bin(dgb);                        % Shortest Paths

    % Compute Centrality Measures
    NDDEGREE_g_t(1,:)           = sum(dgb);                                 % node degree
    CLOSNESS_g_t(1,:)           = 1./mean(Dist);                            % Closeness Centrality
    BETWNESS_g_t(1,:)           = betweenness_bin(dgb);                     % node betweenness centrality
    EIGNVCTR_g_t(1,:)           = eigenvector_centrality_und(dgb);          % Eigenvector Centrality
    PAGERANK_g_t(1,:)           = pagerank_centrality(dgb,.85);             % PageRank Centrality

    % Module-Centrality Measures
    MDDGZSCR_g_t(1,:)           = module_degree_zscore(dgb,CIt,0);          % within-module degree z-score (0=undirected)
    PRTCPACN_g_t(1,:)           = participation_coef(dgb,CIt,0);            % intermodular participation coeff (0=undirected)
    
    % Clustering
    CLUSTRNG_g_t(1,:)           = clustering_coef_bu(dgb);                  % clustering coeff

    
    disp('... Binary Group level completed')
            
            
    %% ----------------------- Binary Subject-level ----------------------- %%
    
    disp(' ')
    disp('Running SUBJECT level Binary ')
    
    % Init storage variables
    NDDEGREE_s_t                = zeros(Nsub,Nnodet);
    CLOSNESS_s_t                = zeros(Nsub,Nnodet);
    BETWNESS_s_t                = zeros(Nsub,Nnodet);
    EIGNVCTR_s_t                = zeros(Nsub,Nnodet);
    PAGERANK_s_t                = zeros(Nsub,Nnodet);
    MDDGZSCR_s_t                = zeros(Nsub,Nnodet);
    PRTCPACN_s_t                = zeros(Nsub,Nnodet);
    CLUSTRNG_s_t                = zeros(Nsub,Nnodet);
    
    % Loop over subjects
    for kk = 1 : Nsub
        tic;
        dt                      = d(:,:,kk);                                % subject connectivity matrix
        dtb                     = double(dt~=0);                            % binarize
        Dist                    = distance_bin(dtb);                        % Shortest Paths

        % Compute Centrality Measures
        NDDEGREE_s_t(kk,:)      = sum(dtb);                                 % node degree
        CLOSNESS_s_t(kk,:)      = 1./mean(Dist);                            % Closeness Centrality
        BETWNESS_s_t(kk,:)      = betweenness_bin(dtb);                     % node betweenness centrality
        EIGNVCTR_s_t(kk,:)      = eigenvector_centrality_und(dtb);          % Eigenvector Centrality
        PAGERANK_s_t(kk,:)      = pagerank_centrality(dtb,.85);             % PageRank Centrality

        % Module-Centrality Measures
        MDDGZSCR_s_t(kk,:)      = module_degree_zscore(dtb,CIt,0);          % within-module degree z-score (0=undirected)
        PRTCPACN_s_t(kk,:)      = participation_coef(dtb,CIt,0);            % intermodular participation coeff (0=undirected)
        
        % Clustering
        CLUSTRNG_s_t(kk,:)      = clustering_coef_bu(dtb);                  % clustering coeff


        % Printout status
        time=toc; 
        disp(['Binary: Subject # ' num2str(kk) ' / ' num2str(Nsub) ' completed in ' num2str(time) 's'])
    end
    disp('... Binary Subject level completed')

    
        

    %% store outputs   

    % Store everything in structure (store & save on each iteration)
    A.grp.bin.degree            = NDDEGREE_g_t;
    A.grp.bin.closeness         = CLOSNESS_g_t;
    A.grp.bin.betweenness       = BETWNESS_g_t;
    A.grp.bin.eigenvector       = EIGNVCTR_g_t;
    A.grp.bin.moddegzscore      = MDDGZSCR_g_t;
    A.grp.bin.participation     = PRTCPACN_g_t;
    A.grp.bin.pagerank          = PAGERANK_g_t;
    A.grp.bin.clustering        = CLUSTRNG_g_t;

    A.sub.bin.degree            = NDDEGREE_s_t;
    A.sub.bin.closeness         = CLOSNESS_s_t;
    A.sub.bin.betweenness       = BETWNESS_s_t;
    A.sub.bin.eigenvector       = EIGNVCTR_s_t;
    A.sub.bin.moddegzscore      = MDDGZSCR_s_t;
    A.sub.bin.participation     = PRTCPACN_s_t;
    A.sub.bin.pagerank          = PAGERANK_s_t;
    A.sub.bin.clustering        = CLUSTRNG_s_t;

    % save on each pass (quick solution to mitigate potential crashing)
    save(savefn,'A')
    
    
end
disp('* Network analysis complete! *')

%--------------------------------------------------------------------------
end