function [A] = conn_gt_richClub(D,W8s,pinfo,varargin)
%  Computes rich club coefficients in binary & weighted connectivity data
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
[savepath,str]      = INhandler(varargin,nin,defaults);     

if isempty(savepath); disp(' * No save location provided * Results will not be saved'); end
                                                
%% Setup
[I,J]                       = size(D);                                      % data dimensions: weights x parcellations
parcs                       = {pinfo.name};                                 % parcellation names
Nnull                       = 100;                                          % number of null networks to generate
edge_swaps                  = 1e6;                                          % Number of edge swaps in null network generation

% Need to exclude FC?
ifc                         = cellstrfind(W8s,'fc');
isc                         = setdiff(1:I,ifc);
I                           = numel(isc);

%% Get network attributes
A                           = [];                                           % primary storage structure
A.info.weights              = W8s(isc);
A.info.parcs                = parcs;
A.info.details              = str;
A.info.edgeswapswei         = edge_swaps;


% Init cell arrays to store outputs of loops (g=group & s=subject)
RCC_g                       = cell(I,J);                                    % Group    
RCC_s                       = cell(I,J);                                    % Subject
RCC_n                       = cell(I,J);                                    % Null

RCC_gb                      = cell(1,J);                                    % Binary Group
RCC_nb                      = cell(1,J);                                    % Binary Null

NDk_gb                      = cell(1,J);                                    % Binary Group Node Count
NDk_nb                      = cell(1,J);                                    % Binary Null Node Count

EGk_gb                      = cell(1,J);                                    % Binary Group Edge Count
EGk_nb                      = cell(1,J);                                    % Binary Null Edge Count

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
        
        %% ----------------------- GROUP-level ------------------------- %%
        dg                  = groupavg(d,3,'nz');                           % group connectivity matrix
        degree              = sum((dg~=0));                                 % node degree
        kmax                = max(degree);                                  % Max k-level used for subject & group level 
        
        % Null networks
        RCCt_n              = zeros(Nnull,kmax);
        for nn = 1 : Nnull
            tic;
            dnull           = fcn_randomize_str(dg,'maxswap',edge_swaps);   % Degree & strength preserving null network
            RCCt_n(nn,:)    = rich_club_wu(dnull,kmax);                     % Rich Club Curve
            time=toc; 
            disp([W8s{ii} ' Null Model # ' num2str(nn) ' completed in ' num2str(time) 's'])
        end
        
        
        % Empirical Networks
        RCCt_g              = rich_club_wu(dg,kmax);                        % Rich Club Curve             
        
        
        disp('... Group level completed')
        
        
        %% --------------------- SUBJECT-level ---------------------- %%
        RCCt_s              = zeros(Sn,kmax);
        for kk = 1 : Sn
            disp(['  - Running subject: ' num2str(kk)])
            dt              = d(:,:,kk);                                    % subject connectivity matrix
            
            RCCt_s(kk,:)    = rich_club_wu(dt,kmax);                        % Rich Club Curve        
             
        end
        disp('... Subject level completed')
        
%         keyboard
             
        
        %% store outputs
        RCC_g{in,jj}            = RCCt_g;       
        RCC_s{in,jj}            = RCCt_s;
        RCC_n{in,jj}            = RCCt_n;

        
        % Store everything in structure (store & save on each iteration)
        A.grp.richclubcurve     = RCC_g;
        A.sub.richclubcurve     = RCC_s;
        A.null.richclubcurve    = RCC_n;

        % save on each pass (quick solution to mitigate potential crashing)
        save(savepath,'A')

        in=in+1;
    end
    
    
    
    %% ----------------------- Binary Group-level ----------------------- %%
    edge_swaps          = 100;                                              % no. of flips/edge
    ii                  = isc(1);                                           % Use 1st structural dataset
    disp(' ')
    disp(['Computing network attributes for group level ' W8s{ii} ' BINARY'])

    d                   = D{ii,jj};
    d(d<0)              = 0;                                                % Set any negative values to 0 (necessary for some metrics)
    dg                  = groupavg(d,3,'nz');                               % group connectivity matrix
    dgb                 = double(dg~=0);                                    % binarize
    degree              = sum(dgb);                                         % node degree
    kmax                = max(degree);                                      % Max k-level used for subject & group level 

    % Null networks
    RCCt_nb             = zeros(Nnull,kmax);
    NDkt_nb             = zeros(Nnull,kmax);                                % Number of NODES with degree > k
    EGkt_nb             = zeros(Nnull,kmax);                                % Number of EDGES with degree > k
    for nn = 1 : Nnull
        tic;
        dnull           = randmio_und(dgb,edge_swaps);                      % Degree preserving null
        [R,Nk,Ek]       = rich_club_bu(dnull,kmax);                         % Rich Club Curve
        RCCt_nb(nn,:)   = R;
        NDkt_nb(nn,:)   = Nk;
        EGkt_nb(nn,:)   = Ek;
        time=toc; 
        disp(['BINARY Null Model # ' num2str(nn) ' completed in ' num2str(time) 's'])
    end


    % Empirical
    [R,Nk,Ek]           = rich_club_bu(dgb,kmax);                           % Rich Club Curve
    RCCt_gb             = R;
    NDkt_gb             = Nk;
    EGkt_gb             = Ek;


    disp('... BINARY Group level completed')
    
   

    %% store outputs   
    RCC_gb{jj}          = RCCt_gb;
    RCC_nb{jj}          = RCCt_nb;
    NDk_gb{jj}          = NDkt_gb;
    NDk_nb{jj}          = NDkt_nb;
    EGk_gb{jj}          = EGkt_gb;
    EGk_nb{jj}          = EGkt_nb;


    % Store everything in structure
    A.info.edgeswapsbin         = edge_swaps;
    A.bin.grp.richclubcurve     = RCC_gb;
    A.bin.null.richclubcurve    = RCC_nb;
    A.bin.grp.nodecount         = NDk_gb;
    A.bin.null.nodecount        = NDk_nb;
    A.bin.grp.edgecount         = EGk_gb;
    A.bin.null.edgecount        = EGk_nb;
    

    % save on each pass (quick solution to mitigate potential crashing)
    save(savepath,'A')
    
    
end
disp('* Network analysis complete! *')

%--------------------------------------------------------------------------
end