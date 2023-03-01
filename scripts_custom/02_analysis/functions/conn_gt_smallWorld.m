function [A] = conn_gt_smallWorld(D,W8s,pinfo,varargin)
%  Computes clustering coefficient and lambda on connectivity data
%
% Input:
%                           + Required +
%   D           : IxJ cell array, each cell containing an NxNxS stack of
%                 connectivity matrices (S=subs)
%   W8s         : Cell array of strings describing edge weights
%   pinfo       : 1xJ structure with parcellation info (see conn_getParcInfo.m)
%
%                           + Optional +
%   savepath    : full file path to save outputs
%   str         : character vector for plot title
%
% Output:
%   A       : 1xJ structure
%
% 2022 Mark C Nelson, McConnell Brain Imaging Centre, MNI, McGill
%--------------------------------------------------------------------------
%% Optional inputs
nin                 = max(nargin,1) - 3;
defaults            = {'',[]};
[savepath, str]     = INhandler(varargin,nin,defaults);    

if isempty(savepath); disp(' * No save location provided * Results will not be saved'); end
                                                
%% Setup
[I,J]                       = size(D);                                      % data dimensions: weights x parcellations
parcs                       = {pinfo.name};                                 % parcellation names
xfm                         = 'inv';                                        % Wij --> Lij transform {'log' 'inv'}
Nnull                       = 50;                                           % number of null networks to generate
edge_swaps                  = 1e6;                                          % Number of edge swaps in null network generation

% Need to exclude FC?
ifc                         = cellstrfind(W8s,'fc');
isc                         = setdiff(1:I,ifc);
I                           = numel(isc);

%% Get network attributes
A               = [];                                                       % primary storage structure
A.info.weights  = W8s(isc);
A.info.parcs    = parcs;
A.info.details  = str;

% Init cell arrays to store outputs of loops (g=group & s=subject)
CLUg  = cell(I,J);          CLUs  = cell(I,J);                              % clustering Coefficient
LMBDg = cell(I,J);          LMBDs = cell(I,J);                              % characteristic path length
DISTg = cell(I,J);          DISTs = cell(I,J);                              % distance (shortest path) matrices
HOPSg = cell(I,J);          HOPSs = cell(I,J);                              % number of edges in shortest path
PMATg = cell(I,J);          PMATs = cell(I,J);                              % 

% For Null network measures
CLUnull  = cell(I,J);
LMBDnull = cell(I,J);

% Iterate over parcellations (2nd dim of D)
for jj = 1 : J
    disp(['* * * Processing parcellation scheme ' num2str(jj) ': ' parcs{jj}])
    
    % Iterate over weights (1st dim of D)
    in=1;
    for ii = isc % 1 : I
        disp(' ')
        disp(['Computing network attributes for dataset ' num2str(ii) ': ' W8s{ii}])
        
        d                   = D{ii,jj};
        d(d<0)              = 0;                                            % Set any negative values to 0 (necessary for some metrics)
        Sn                  = size(d,3);                                    % assume 3rd dim is subjects
        Nnode               = size(d,2);
        
        % Temporary storage arrays
        CLUgt  = zeros(1,Nnode);           
        CLUst  = zeros(Sn,Nnode);
                                        
        LMBDst = zeros(Sn,1);
        
        DISTst = zeros(Nnode,Nnode,Sn);
        HOPSst = zeros(Nnode,Nnode,Sn);
        PMATst = zeros(Nnode,Nnode,Sn);
        
        
        %% ----------------------- GROUP-level ------------------------- %%
        dg                      = groupavg(d,3,'nz');                           % group connectivity matrix
        
        % Null networks
        CLUnullt                = zeros(Nnull,Nnode);
        LMBDnullt               = zeros(Nnull,1);
        
        for nn = 1:Nnull
            tic;
            dnull               = fcn_randomize_str(dg,'maxswap',edge_swaps);   % Degree & strength preserving null network
            CLUnullt(nn,:)      = clustering_coef_wu(dnull);                    % Clustering Coeff
            
            [Dist,~,~]          = distance_wei_floyd(dnull,xfm);
            LMBDnullt(nn)       = charpath(Dist);                               % Characteristic path length
            
            % Check for infinite path length
            if isinf(LMBDnullt(nn))
                LMBDnullt(nn)   = charpath(Dist,0,0);                           % Network could contain disconnected nodes
            end
            if isinf(LMBDnullt(nn)) 
                disp(['WARNING: Infinite path length found: Null ' num2str(nn) ' ' W8s{ii}])
            end
            
            time=toc; 
            disp([W8s{ii} ' Null Model# ' num2str(nn) ' completed in ' num2str(time) 's'])
        end
        
        
        % Empirical Networks
        CLUgt(1,:)          = clustering_coef_wu(dg);                       % Clustering Coeff             
        
%         Leng                        = weight_conversion(dg,'lengths');      % length matrix conversion
%         Dist                        = distance_wei(Leng);                   % distance matrix (shortest weighted path)
        [Dist,hops,Pmat]    = distance_wei_floyd(dg,xfm);
        LMBDgt              = charpath(Dist);                               % Characteristic path length
        
        % Check for infinite path length
        if isinf(LMBDgt)
            LMBDgt          = charpath(Dist,0,0);                           % Network could contain disconnected nodes
        end
        if isinf(LMBDgt) 
            disp(['WARNING: Infinite path length found: Group ' W8s{ii}])
        end
        
        DISTgt              = Dist;
        HOPSgt              = hops;
        PMATgt              = Pmat;
        
        disp('... Group level completed')
        
        
        %% --------------------- SUBJECT-level ---------------------- %%
        for kk = 1 : Sn
            disp(['  - Running subject: ' num2str(kk)])
            dt                      = d(:,:,kk);                            % subject connectivity matrix
            
            
            % Clustering Coefficient (weighted undirected)
            CLUst(kk,:)             = clustering_coef_wu(dt);               
                        
            % Characteristic path length
%             Leng                    = weight_conversion(dt,'lengths');      % length matrix
%             Dist                    = distance_wei(Leng);                   % distance matrix (shortest weighted path)
            [Dist,hops,Pmat]        = distance_wei_floyd(dt,xfm);
            LMBDst(kk)              = charpath(Dist);
            
            % Check for infinite path length
            if isinf(LMBDst(kk))
                LMBDst(kk)          = charpath(Dist,0,0);                   % Network could contain disconnected nodes
            end
            if isinf(LMBDst(kk)) 
                disp(['WARNING: Infinite path length found: Subject ' num2str(kk) ' ' W8s{ii}])
            end
            
            DISTst(:,:,kk)          = Dist;
            HOPSst(:,:,kk)          = hops;
            PMATst(:,:,kk)          = Pmat;
                       
        end
        disp('... Subject level completed')
        
        
        %% store outputs
        CLUg{in,jj}     = CLUgt;       CLUs{in,jj}      = CLUst;
        LMBDg{in,jj}    = LMBDgt;      LMBDs{in,jj}     = LMBDst;
        DISTg{in,jj}    = DISTgt;      DISTs{in,jj}     = DISTst;
        HOPSg{in,jj}    = HOPSgt;      HOPSs{in,jj}     = HOPSst;
        PMATg{in,jj}    = PMATgt;      PMATs{in,jj}     = PMATst;
        
        CLUnull{in,jj}  = CLUnullt;    LMBDnull{in,jj}  = LMBDnullt;

        
        % Store everything in structure (store & save on each iteration)
        A.grp.clusteringcoeff                       = CLUg;
        A.grp.charpath                              = LMBDg;
        A.grp.dist.dist                             = DISTg;
        A.grp.dist.hops                             = HOPSg;
        A.grp.dist.pmat                             = PMATg;

        A.sub.clusteringcoeff                       = CLUs;
        A.sub.charpath                              = LMBDs;
        A.sub.dist.dist                             = DISTs;
        A.sub.dist.hops                             = HOPSs;
        A.sub.dist.pmat                             = PMATs;
        
        A.null.clusteringcoeff                      = CLUnull;
        A.null.charpath                             = LMBDnull;

        % save on each pass (quick solution to mitigate potential crashing)
        save(savepath,'A')

        in=in+1;
    end
end
disp('* Network analysis complete! *')

%--------------------------------------------------------------------------
end