function [A] = conn_gt_all(D,W8s,pinfo,varargin)
%  Returns a range of network attributes computed by applying graph theory
%  to connectivity data
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
%               A.degree  : 
%
% 2021 Mark C Nelson, McConnell Brain Imaging Centre, MNI, McGill
%--------------------------------------------------------------------------
%% Optional inputs
nin                         = max(nargin,1) - 3;

if nin < 1; str             = []; end

if nin > 0
    for kk = 1:nin
        switch kk
            case 1
                strin       = varargin{kk};
            case 2
                PARCS       = varargin{kk};
        end  
    end
end
                                                
                                                
%% Setup
[I,J]                       = size(D);                                      % data dimensions: weights x parcellations
parcs                       = {pinfo.name};                                 % parcellation names
CIs                         = {pinfo.cirois};                               % community assignments
fci                         = cellstrfind(W8s,'fc');
sci                         = setdiff(1:I,fci);
gamma                       = 1;                                            % modularity resolution parameter (>1=smaller; [0 1)=larger)
xfm                         = 'log';                                        % Wij --> Lij transform
Nnull                       = 50;                                           % number of null networks to generate
edge_swaps                  = 1e6;                                          % Number of edge swaps in null network generation

%% Get network attributes
A                   = [];                                                       % primary storage structure
A.info.weights      = W8s;
A.info.parcs        = parcs;
A.info.fci          = fci;
A.info.sci          = sci;
A.info.xfm          = xfm;
A.info.Nnull        = Nnull;
A.info.edgeswaps    = edge_swaps;
A.info.gamma        = gamma;
A.info.details      = strin;

% Init cell arrays to store outputs of loops (g=group & s=subject)
DEGg  = cell(I,J);          DEGs  = cell(I,J);                              % node degree  
STRg  = cell(I,J);          STRs  = cell(I,J);                              % node strength
CLUg  = cell(I,J);          CLUs  = cell(I,J);                              % clustering Coefficient
Qg    = cell(I,J);          Qs    = cell(I,J);                              % modularity
Tg    = cell(I,J);          Ts    = cell(I,J);                              % transitivity
ASTg  = cell(I,J);          ASTs  = cell(I,J);                              % assortativity
SGCg  = cell(I,J);          SGCs  = cell(I,J);                              % node subgraph centrality
RCCg  = cell(I,J);          RCCs  = cell(I,J);                              % rich club coefficients curve
ELOCg = cell(I,J);          ELOCs = cell(I,J);                              % efficiency local
NBCg  = cell(I,J);          NBCs  = cell(I,J);                              % node betweenness centrality
EBCg  = cell(I,J);          EBCs  = cell(I,J);                              % edge betweenness centrality
SIg   = cell(I,J);          SIs   = cell(I,J);                              % search information
PTVg  = cell(I,J);          PTVs  = cell(I,J);                              % path transitivity
MDZg  = cell(I,J);          MDZs  = cell(I,J);                              % module degree zscore
PNCg  = cell(I,J);          PNCs  = cell(I,J);                              % participation coefficient
LMBDg = cell(I,J);          LMBDs = cell(I,J);                              % characteristic path length
EFFg  = cell(I,J);          EFFs  = cell(I,J);                              % global efficiency
ECCg  = cell(I,J);          ECCs  = cell(I,J);                              % nodal eccentricity (max path length)
RADg  = cell(I,J);          RADs  = cell(I,J);                              % radius (min eccentricity)
DIAMg = cell(I,J);          DIAMs = cell(I,J);                              % diamater (max eccentricity)
REPg  = cell(I,J);          REPs  = cell(I,J);                              % routing efficiency pair-wise
RELg  = cell(I,J);          RELs  = cell(I,J);                              % routing efficiency local


% Iterate over parcellations (2nd dim of D)
for jj = 1 : J
    disp(['* * * Processing parcellation scheme ' num2str(jj) ': ' parcs{jj}])
    CIt                     = CIs{jj};
    
    % Iterate over weights (1st dim of D)
    for ii = 1 : I
        disp(' ')
        disp(['Computing network attributes for dataset ' num2str(ii) ': ' W8s{ii}])
        
        d                   = D{ii,jj};
        d(d<0)              = 0;                                            % Set any negative values to 0 (necessary for some metrics)
        [Nnode,~,Nsub]      = size(d);
        
        % Temporary storage arrays
        DEGgt  = zeros(1,Nnode);    DEGst  = zeros(Nsub,Nnode);
        STRgt  = zeros(1,Nnode);    STRst  = zeros(Nsub,Nnode);
        CLUgt  = zeros(1,Nnode);    CLUst  = zeros(Nsub,Nnode);
                                    Qst    = zeros(Nsub,1);
                                    Tst    = zeros(Nsub,1);
                                    ASTst  = zeros(Nsub,1);
        SGCgt  = zeros(1,Nnode);    SGCst  = zeros(Nsub,Nnode);
        
        ELOCgt = zeros(1,Nnode);    ELOCst = zeros(Nsub,Nnode);
        NBCgt  = zeros(1,Nnode);    NBCst  = zeros(Nsub,Nnode);
                                    EBCst  = zeros(Nnode,Nnode,Nsub);
                                    SIst   = zeros(Nnode,Nnode,Nsub);
                                    PTVst  = zeros(Nnode,Nnode,Nsub);
        RCCgt  = nan(1,Nnode);      RCCst  = nan(Nsub,Nnode);
        
        MDZgt  = zeros(1,Nnode);    MDZst  = zeros(Nsub,Nnode);
        PNCgt  = zeros(1,Nnode);    PNCst  = zeros(Nsub,Nnode);
        
                                    LMBDst = zeros(Nsub,1);
                                    EFFst  = zeros(Nsub,1);
        ECCgt  = zeros(1,Nnode);    ECCst  = zeros(Nsub,Nnode);
                                    RADst  = zeros(Nsub,1);
                                    DIAMst = zeros(Nsub,1);
                                        
                                    REPst  = zeros(Nnode,Nnode,Nsub);
        RELgt  = zeros(1,Nnode);    RELst  = zeros(Nsub,Nnode);
        
        
        %% ----------------------- GROUP-level ------------------------- %%
        dg                          = groupavg(d,3,'nz');                   % group connectivity matrix
        dgb                         = double(dg~=0);                        % binary
        
        % Compute length & distance matrices
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
%         Dist                        = distance_wei(Leng);                   
        [Dist,~,~]                  = distance_wei_floyd(dg,xfm);           % distance matrix (shortest weighted path)
                
        
%         % Network measures computed only in SC data (Doesn't support alternte W-->L transforms)
%         if ~ismember(ii, fci)
%             ELOCgt(1,:)             = efficiency_wei(dg,2);                 % local efficiency (internally converts Lij = 1/Wij)
%         end
        
        % Network measures derived directly from connectome data
        DEGgt(1,:)                  = degrees_und(dg);                      % node degree
        STRgt(1,:)                  = strengths_und(dg);                    % node strength
        CLUgt(1,:)                  = clustering_coef_wu(dg);               % clustering coeff
        Tgt                         = transitivity_wu(dg);                  % transitivity
        ASTgt                       = assortativity_wei(dg,0);              % assortativity weighted (0=undirected)
        SGCgt(1,:)                  = subgraph_centrality(dgb);             % node subgraph centrality (binary matrix)
%         EGLOgt                      = efficiency_wei(dg);                    % global efficiency (internally converts Lij = 1/Wij) (same as charpath-efficiency)
        
        % rich club coefficients curve (weighted undirected)
        rccst               = rich_club_wu(dg);
        lnt                 = length(rccst);
        RCCgt(1,1:lnt)      = rccst;                                        % size varies by dataset...
        
        % Network measures derived from length matrix
        NBCgt(1,:)                  = betweenness_wei(Leng);                % node betweenness centrality
        EBCgt                       = edge_betweenness_wei(Leng);           % edge betweenness centrality (matrix)
        SIgt                        = search_information(dg,Leng);          % search information, no memory (matrix)
        PTVgt                       = path_transitivity(Leng);              % path transitivity (matrix)

        % Degree centrality
        MDZgt(1,:)                  = module_degree_zscore(dg,CIt,0);       % within-module degree z-score (0=undirected)
        PNCgt(1,:)                  = participation_coef(dg,CIt,0);         % intermodular participation coeff (0=undirected)
        
        % Characteristic path length
        [lmd,eff,ecc,rad,diam]      = charpath(Dist);
        LMBDgt                      = lmd;                                  % characteristic path length
        EFFgt                       = eff;                                  % efficiency global
        ECCgt(1,:)                  = ecc;                                  % nodal eccentricity
        RADgt                       = rad;                                  % radius
        DIAMgt                      = diam;                                 % diameter
        
        % Routing efficiency
        [~,PErout,Eloc]         = rout_efficiency(Leng);
        REPgt                   = PErout;                                   % pair-wise routing efficiency
        RELgt(1,:)              = Eloc;                                     % local routing efficiency
        
        % Modularity
        %             [M,Q]=community_louvain(dt,[],CIt);
        s                           = sum(dg(:));                           % total connectome weight
        B                           = (dg-gamma*(sum(dg,2)*sum(dg,1))/s)/s; % modularity matrix (non-zero weights)
        Qgt                         = sum(B(bsxfun(@eq,CIt,CIt.')));        % compute modularity
        
        
        %% Compute Null Networks on Group Data & rerun network analysis
        
        % Temporary storage arrays
        CLU_n                       = zeros(Nnull,Nnode);
        Q_n                         = zeros(Nnull,1);
        T_n                         = zeros(Nnull,1);
        AST_n                       = zeros(Nnull,1);
        SGC_n                       = zeros(Nnull,Nnode);
%         ELOC_n                    = zeros(Nnull,Nnode);
        NBC_n                       = zeros(Nnull,Nnode);
        EBC_n                       = zeros(Nnode,Nnode,Nnull);
        SI_n                        = zeros(Nnode,Nnode,Nnull);
        PTV_n                       = zeros(Nnode,Nnode,Nnull);
        MDZ_n                       = zeros(Nnull,Nnode);
        PNC_n                       = zeros(Nnull,Nnode);
        LMBD_n                      = zeros(Nnull,1);
        EFF_n                       = zeros(Nnull,1);
        ECC_n                       = zeros(Nnull,Nnode);
        RAD_n                       = zeros(Nnull,1);
        DIAM_n                      = zeros(Nnull,1);
        REP_n                       = zeros(Nnode,Nnode,Nnull);
        REL_n                       = zeros(Nnull,Nnode);
        DNULL                       = zeros(Nnode,Nnode,Nnull);
        
        for nn = 1 : Nnull
            tic;
            dnullt                  = fcn_randomize_str(dg,'maxswap',edge_swaps);   % Degree & strength preserving nulls
            
            dnullb                  = double(dnullt~=0);                         % binary

            % Compute length & distance matrices
            switch xfm
                case 'log'
                    Leng_n          = -log(dnullt);
                case 'inv'
                    Leng_n          = 1./dnullt;
            end
%             Dist_n                  = distance_wei(Leng_n);                   % distance matrix (shortest weighted path)
            [Dist_n,~,~]            = distance_wei_floyd(dnullt,xfm);       % distance matrix (shortest weighted path)

            % Network measures derived directly from connectome data
            CLU_n(nn,:)             = clustering_coef_wu(dnullt);               % clustering coeff
            T_n(nn)                 = transitivity_wu(dnullt);                  % transitivity
            AST_n(nn)               = assortativity_wei(dnullt,0);              % assortativity weighted (0=undirected)
            SGC_n(nn,:)             = subgraph_centrality(dnullb);             % node subgraph centrality (binary matrix)

            % Network measures derived from length matrix
            NBC_n(nn,:)             = betweenness_wei(Leng_n);                % node betweenness centrality
            EBC_n(:,:,nn)           = edge_betweenness_wei(Leng_n);           % edge betweenness centrality (matrix)
            SI_n(:,:,nn)            = search_information(dnullt,Leng_n);          % search information, no memory (matrix)
            PTV_n(:,:,nn)           = path_transitivity(Leng_n);              % path transitivity (matrix)

            % Degree centrality
            MDZ_n(nn,:)             = module_degree_zscore(dnullt,CIt,0);           % within-module degree z-score (0=undirected)
            PNC_n(nn,:)             = participation_coef(dnullt,CIt,0);             % intermodular participation coeff (0=undirected)

            % Characteristic path length
            [lmd,eff,ecc,rad,diam]  = charpath(Dist_n);
            LMBD_n(nn)              = lmd;                                          % characteristic path length
            EFF_n(nn)               = eff;                                          % efficiency global
            ECC_n(nn,:)             = ecc;                                          % nodal eccentricity
            RAD_n(nn)               = rad;                                          % radius
            DIAM_n(nn)              = diam;                                         % diameter

            % Routing efficiency
            [~,PErout,Eloc]         = rout_efficiency(Leng_n);
            REP_n(:,:,nn)           = PErout;                                       % pair-wise routing efficiency
            REL_n(nn,:)             = Eloc;                                         % local routing efficiency

            % Modularity
            s                       = sum(dnullt(:));                               % total connectome weight
            B                       = (dnullt-gamma*(sum(dnullt,2)*sum(dnullt,1))/s)/s; % modularity matrix (non-zero weights)
            Q_n(nn)                 = sum(B(bsxfun(@eq,CIt,CIt.')));                % compute modularity

            
            DNULL(:,:,nn)           = dnullt;
            time=toc; 
            disp(['Network Attributes for ' W8s{ii} ' Null Model ' num2str(nn) ' / ' num2str(Nnull)  ' completed in ' num2str(time) 's'])
        end
        
        %% Normalize Group-Level Network Measures

        % Average over Null Networks
        CLU_n                       = mean(CLU_n);
        Q_n                         = mean(Q_n);
        T_n                         = mean(T_n);
        AST_n                       = mean(AST_n);
        SGC_n                       = mean(SGC_n);
%         ELOC_n                    = mean(ELOC_n);
        NBC_n                       = mean(NBC_n);
        EBC_n                       = mean(EBC_n,3);
        SI_n                        = mean(SI_n,3);
        PTV_n                       = mean(PTV_n,3);
        MDZ_n                       = mean(MDZ_n);
        PNC_n                       = mean(PNC_n);
        LMBD_n                      = mean(LMBD_n);
        EFF_n                       = mean(EFF_n);
        ECC_n                       = mean(ECC_n);
        RAD_n                       = mean(RAD_n);
        DIAM_n                      = mean(DIAM_n);
        REP_n                       = mean(REP_n,3);
        REL_n                       = mean(REL_n);
        
        % Normalize
        CLUgt                       = CLUgt./CLU_n;
        Qgt                         = Qgt/Q_n;
        Tgt                         = Tgt/T_n;
        ASTgt                       = ASTgt/AST_n;
        SGCgt                       = SGCgt./SGC_n;
%         ELOC_n
        NBCgt                       = NBCgt./NBC_n;
        EBCgt                       = EBCgt./EBC_n;
        SIgt                        = SIgt./SI_n;
        PTVgt                       = PTVgt./PTV_n;
        MDZgt                       = MDZgt./MDZ_n;
        PNCgt                       = PNCgt./PNC_n;
        LMBDgt                      = LMBDgt/LMBD_n;
        EFFgt                       = EFFgt/EFF_n;
        ECCgt                       = ECCgt./ECC_n;
        RADgt                       = RADgt/RAD_n;
        DIAMgt                      = DIAMgt/DIAM_n;
        REPgt                       = REPgt./REP_n;
        RELgt                       = RELgt./REL_n;
        
        disp('... Group level completed')

        
        %% --------------------- SUBJECT-level ---------------------- %%
        for kk = 1 : Nsub
            disp(['  - Running subject: ' num2str(kk)])
            dt                      = d(:,:,kk);                            % subject connectivity matrix
            dbt                     = double(dt~=0);                        % binary

            % Compute length & distance matrices
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
%             Dist                    = distance_wei(Leng);                   
            [Dist,~,~]              = distance_wei_floyd(dt,xfm);           % distance matrix (shortest weighted path)
            
            % Network measures computed only in SC data
            if ~ismember(ii, fci)
                DEGst(kk,:)         = degrees_und(dt);                      % node degree
                ASTst(kk)           = assortativity_wei(dt,0);              % assortativity (weighted, 0=undirected)
                SGCst(kk,:)         = subgraph_centrality(dbt);             % node subgraph centrality (binary matrix)
                PTVst(:,:,kk)       = path_transitivity(Leng);              % path transitivity (matrix)
                ELOCst(kk,:)        = efficiency_wei(dt,2);                 % local efficiency (group-level SC only?)
                
                % rich club coefficients curve (weighted undirected)
                rccst               = rich_club_wu(dt);
                lnt                 = length(rccst);
                RCCst(kk,1:lnt)     = rccst;                                % size varies by dataset...
            end
            
            % Network measures derived directly from connectome data
            STRst(kk,:)             = strengths_und(dt);                    % node strength (unsigned undirected)
            CLUst(kk,:)             = clustering_coef_wu(dt);               % clustering coefficient (weighted undirected)
            Tst(kk)                 = transitivity_wu(dt);                  % transitivity (weighted undirected)
%             EGLOst(kk)              = efficiency_wei(dt);                   % global efficiency (same as charpath-efficiency)
                        
            % Network measures derived from length matrix
            NBCst(kk,:)             = betweenness_wei(Leng);                % node betweenness centrality
            EBCst(:,:,kk)           = edge_betweenness_wei(Leng);           % edge betweenness centrality (matrix)
            SIst(:,:,kk)            = search_information(dt,Leng);          % search information, no memory (matrix)
            
            % Degree centrality
            MDZst(kk,:)             = module_degree_zscore(dt,CIt,0);       % within-module degree z-score (0=undirected)
            PNCst(kk,:)             = participation_coef(dt,CIt,0);         % intermodular participation coeff (0=undirected)

            % Characteristic path length
            [lmd,eff,ecc,rad,diam]  = charpath(Dist);
            LMBDst(kk)              = lmd;                                  % characteristic path length
            EFFst(kk)               = eff;                                  % efficiency global
            ECCst(kk,:)             = ecc;                                  % nodal eccentricity
            RADst(kk)               = rad;                                  % radius
            DIAMst(kk)              = diam;                                 % diameter
            
            % Routing efficiency
            [~,PErout,Eloc]         = rout_efficiency(Leng);
            REPst(:,:,kk)           = PErout;                               % pair-wise
            RELst(kk,:)             = Eloc;                                 % local
            
            % Modularity (weights must be positive for this implementation)
            s                       = sum(dt(:));                           % total connectome weight
            B                       = (dt-gamma*(sum(dt,2)*sum(dt,1))/s)/s; % modularity matrix (non-zero weights)
            Qst(kk)                 = sum(B(bsxfun(@eq,CIt,CIt.')));        % compute modularity
            
            
            %% Normalize Subject-Level Network Measures

            % Normalize
            CLUst(kk,:)                 = CLUst(kk,:)./CLU_n;
            Qst(kk)                     = Qst(kk)/Q_n;
            Tst(kk)                     = Tst(kk)/T_n;
            ASTst(kk)                   = ASTst(kk)/AST_n;
            SGCst(kk,:)                 = SGCst(kk,:)./SGC_n;
    %         ELOC_n
            NBCst(kk,:)                 = NBCst(kk,:)./NBC_n;
            EBCst(:,:,kk)               = EBCst(:,:,kk)./EBC_n;
            SIst(:,:,kk)                = SIst(:,:,kk)./SI_n;
            PTVst(:,:,kk)               = PTVst(:,:,kk)./PTV_n;
            MDZst(kk,:)                 = MDZst(kk,:)./MDZ_n;
            PNCst(kk,:)                 = PNCst(kk,:)./PNC_n;
            LMBDst(kk)                  = LMBDst(kk)/LMBD_n;
            EFFst(kk)                   = EFFst(kk)/EFF_n;
            ECCst(kk,:)                 = ECCst(kk,:)./ECC_n;
            RADst(kk)                   = RADst(kk)/RAD_n;
            DIAMst(kk)                  = DIAMst(kk)/DIAM_n;
            REPst(:,:,kk)               = REPst(:,:,kk)./REP_n;
            RELst(kk,:)                 = RELst(kk,:)./REL_n;
        end
        disp('... Subject level completed')
        
        
        
        %% store outputs
        DEGg{ii,jj}  = DEGgt;       DEGs{ii,jj}  = DEGst;
        STRg{ii,jj}  = STRgt;       STRs{ii,jj}  = STRst;
        CLUg{ii,jj}  = CLUgt;       CLUs{ii,jj}  = CLUst;
        Qg{ii,jj}    = Qgt;         Qs{ii,jj}    = Qst;
        Tg{ii,jj}    = Tgt;         Ts{ii,jj}    = Tst;
        ASTg{ii,jj}  = ASTgt;       ASTs{ii,jj}  = ASTst;
        SGCg{ii,jj}  = SGCgt;       SGCs{ii,jj}  = SGCst;
        RCCg{ii,jj}  = RCCgt;       RCCs{ii,jj}  = RCCst;
        
        ELOCg{ii,jj} = ELOCgt;      ELOCs{ii,jj} = ELOCst;
        
        NBCg{ii,jj}  = NBCgt;       NBCs{ii,jj}  = NBCst;        
        EBCg{ii,jj}  = EBCgt;       EBCs{ii,jj}  = EBCst;        
        SIg{ii,jj}   = SIgt;        SIs{ii,jj}   = SIst;
        PTVg{ii,jj}  = PTVgt;       PTVs{ii,jj}  = PTVst;
        
        MDZg{ii,jj}  = MDZgt;       MDZs{ii,jj}  = MDZst;
        PNCg{ii,jj}  = PNCgt;       PNCs{ii,jj}  = PNCst;
        
        LMBDg{ii,jj} = LMBDgt;      LMBDs{ii,jj} = LMBDst;
        EFFg{ii,jj}  = EFFgt;       EFFs{ii,jj}  = EFFst;
        ECCg{ii,jj}  = ECCgt;       ECCs{ii,jj}  = ECCst;
        RADg{ii,jj}  = RADgt;       RADs{ii,jj}  = RADst;
        DIAMg{ii,jj} = DIAMgt;      DIAMs{ii,jj} = DIAMst;
        
        REPg{ii,jj}  = REPgt;       REPs{ii,jj}  = REPst;
        RELg{ii,jj}  = RELgt;       RELs{ii,jj}  = RELst;
        
        % Store everything in structure (store & save on each iteration)
        A.grp.degree                                = DEGg;
        A.grp.strength                              = STRg;
        A.grp.clusteringcoeff                       = CLUg;
        A.grp.transitivity                          = Tg;
        A.grp.modularity                            = Qg;
        A.grp.assortativity                         = ASTg;
        A.grp.subgraphcentrality                    = SGCg;
        A.grp.richclubcurve                         = RCCg;
        A.grp.efficiency_local                      = ELOCg;
        A.grp.moduledegreezscore                    = MDZg;
        A.grp.participationcoeff                    = PNCg;
        A.grp.nodebetweenness                       = NBCg;
        A.grp.edgebetweenness                       = EBCg;
        A.grp.searchinformation                     = SIg;
        A.grp.pathtransitivity                      = PTVg;
        A.grp.charpath.lambda                       = LMBDg;
        A.grp.charpath.efficiency                   = EFFg;
        A.grp.charpath.eccentricity                 = ECCg;
        A.grp.charpath.radius                       = RADg;
        A.grp.charpath.diameter                     = DIAMg;
        A.grp.routeffic.pairwise                    = REPg;
        A.grp.routeffic.local                       = RELg;

        A.sub.degree                                = DEGs;
        A.sub.strength                              = STRs;
        A.sub.clusteringcoeff                       = CLUs;
        A.sub.transitivity                          = Ts;
        A.sub.modularity                            = Qs;
        A.sub.assortativity                         = ASTs;
        A.sub.subgraphcentrality                    = SGCs;
        A.sub.richclubcurve                         = RCCs;
        A.sub.efficiency_local                      = ELOCs;
        A.sub.moduledegreezscore                    = MDZs;
        A.sub.participationcoeff                    = PNCs;
        A.sub.nodebetweenness                       = NBCs;
        A.sub.edgebetweenness                       = EBCs;
        A.sub.searchinformation                     = SIs;
        A.sub.pathtransitivity                      = PTVs;
        A.sub.charpath.lambda                       = LMBDs;
        A.sub.charpath.efficiency                   = EFFs;
        A.sub.charpath.eccentricity                 = ECCs;
        A.sub.charpath.radius                       = RADs;
        A.sub.charpath.diameter                     = DIAMs;
        A.sub.routeffic.pairwise                    = REPs;
        A.sub.routeffic.local                       = RELs;

        % save on each pass (quick solution to mitigate potential crashing)
        save('~/Desktop/ConnectomeProj/myelinWeightedConnectome/0_SCFingerprinting/data/5_networkAtts_Uthr50_UniDensity_schaefer400_redux.mat','-v7.3','A')

        
    end
end
disp('* Network analysis complete! *')
keyboard
% % % % Store everything in structure
% % % A=[];
% % % A.grp.degree                                = DEGg;
% % % A.grp.strength                              = STRg;
% % % A.grp.clusteringcoeff                       = CLUg;
% % % A.grp.transitivity                          = Tg;
% % % A.grp.modularity                            = Qg;
% % % A.grp.assortativity                         = ASTg;
% % % A.grp.subgraphcentrality                    = SGCg;
% % % A.grp.richclubcurve                         = RCCg;
% % % A.grp.efficiency_local                      = ELOCg;
% % % A.grp.moduledegreezscore                    = MDZg;
% % % A.grp.participationcoeff                    = PNCg;
% % % A.grp.nodebetweenness                       = NBCg;
% % % A.grp.edgebetweenness                       = EBCg;
% % % A.grp.searchinformation                     = SIg;
% % % A.grp.pathtransitivity                      = PTVg;
% % % A.grp.charpath.lambda                       = LMBDg;
% % % A.grp.charpath.efficiency                   = EFFg;
% % % A.grp.charpath.eccentricity                 = ECCg;
% % % A.grp.charpath.radius                       = RADg;
% % % A.grp.charpath.diameter                     = DIAMg;
% % % A.grp.routeffic.pairwise                    = REPg;
% % % A.grp.routeffic.local                       = RELg;
% % % 
% % % A.sub.degree                                = DEGs;
% % % A.sub.strength                              = STRs;
% % % A.sub.clusteringcoeff                       = CLUs;
% % % A.sub.transitivity                          = Ts;
% % % A.sub.modularity                            = Qs;
% % % A.sub.assortativity                         = ASTs;
% % % A.sub.subgraphcentrality                    = SGCs;
% % % A.sub.richclubcurve                         = RCCs;
% % % A.sub.efficiency_local                      = ELOCs;
% % % A.sub.moduledegreezscore                    = MDZs;
% % % A.sub.participationcoeff                    = PNCs;
% % % A.sub.nodebetweenness                       = NBCs;
% % % A.sub.edgebetweenness                       = EBCs;
% % % A.sub.searchinformation                     = SIs;
% % % A.sub.pathtransitivity                      = PTVs;
% % % A.sub.charpath.lambda                       = LMBDs;
% % % A.sub.charpath.efficiency                   = EFFs;
% % % A.sub.charpath.eccentricity                 = ECCs;
% % % A.sub.charpath.radius                       = RADs;
% % % A.sub.charpath.diameter                     = DIAMs;
% % % A.sub.routeffic.pairwise                    = REPs;
% % % A.sub.routeffic.local                       = RELs;
% % % 
% % % 
% % % save('~/Desktop/ConnectomeProj/myelinWeightedConnectome/0_SCFingerprinting/data/5_networkAtts_Uthr50_UniDensity.mat','-v7.3','A')
% % % 
% load('~/Desktop/ConnectomeProj/myelinWeightedConnectome/0_SCFingerprinting/data/5_networkAtts_Uthr50_UniDensity.mat')

%% Plot them
% % % fntsz           = 20;
% % % cmap            = hsv(I);
% % % lnstyl          = {'-' '--' ':'};
% % % maxsp           = 20;                                                       % maximum number of subplots
% % % if (I) > maxsp
% % %     n           = 5;                                                        % dims for subplot
% % %     m           = 4;
% % %     ppos        = [1442 -539 2559 1336];                                    % plot size & location [1442 -539 2350 1336]
% % % else
% % %     n           = ceil(sqrt(I));
% % %     m           = ceil((I)/n);
% % %     ppos        = [1442 -457 2226 1254];
% % % end
% % % 
% % % %% Node Degree
% % % Dgrp        = A.grp.degree;
% % % pltstr      = ['Degree distributions: group-level (' strin ')'];
% % % 
% % % % Histograms (group)
% % % for jj = 1:J
% % %     Nt      = N(jj);
% % %     nbins   = Nt*.1;
% % %     myfig([pltstr ': ' parcs{jj}], ppos)
% % %     
% % %     % Find better x limits for visualization
% % %     xmaxsc  = max(cellfun(@(c) max(c),DEGg(sci)));
% % %     xlimssc = [0 xmaxsc+xmaxsc*.01];                                        % structural
% % %     xminfc  = min(cellfun(@(c) min(c),DEGg(fci)));
% % %     xmaxfc  = max(cellfun(@(c) max(c),DEGg(fci)));
% % %     xlimsfc = [xminfc-xminfc*.01  xmaxfc+xmaxfc*.01];                       % functional
% % %     for ii = 1 : I
% % %         d   = DEGg{ii,jj};
% % %         subplot(m,n,ii)
% % %         histogram(d,nbins,'Normalization','probability')
% % %         title(W8s{ii});
% % %         xlabel('node degree (k)'); 
% % %         ylabel('Probability (degree=k)')
% % %         set(gca,'FontSize', fntsz)
% % %         if ismember(ii,fci)
% % %             xlim(xlimsfc)
% % %         else
% % %             xlim(xlimssc)
% % %         end
% % %     end 
% % % end
% % % 
% % % 
% % % %% Node Strength
% % % Dgrp        = A.grp.strength;
% % % pltstr      = ['Node Strength: group-level (' strin ')'];
% % % 
% % % % Histograms (group)
% % % for jj = 1:J
% % %     Nt      = N(jj);
% % %     nbins   = Nt*.1;
% % %     myfig([pltstr ': ' parcs{jj}], ppos)
% % %     for ii = 1 : I
% % %         d   = Dgrp{ii,jj};
% % %         subplot(m,n,ii)
% % %         histogram(d,nbins,'Normalization','count')
% % %         title(W8s{ii});
% % %         ylim([0 50])
% % %         xlabel('Node Strength (S)'); 
% % %         ylabel('Count (strength=S)')
% % %         set(gca,'FontSize', fntsz) 
% % %     end 
% % % end
% % % 
% % % 
% % % %% ------------- Cumulative Distributions -------------- %%
% % % % Strength, Clustering Coefficient, NodeBetweenness, Participation Coefficient, Local Efficiency
% % % % datanames   = {'strength' 'clusteringcoeff' 'efficiency_local' 'nodebetweenness' 'participationcoeff' };
% % % % datanames   = {'strength' 'clusteringcoeff' 'efficiency_local'};
% % % datanames   = {'eccentricity'};
% % % isel        = [1 2 3];         
% % % % isel        = [8 10 12];         
% % % % isel        = [3 4 12];
% % % W8sel       = W8s(isel);
% % % W8sel=strrep(W8sel,'(log)','');W8sel=strrep(W8sel,'-u','');W8sel=strrep(W8sel,'-s','');W8sel=strrep(W8sel,'-c','');W8sel=strrep(W8sel,'nos','NoS');
% % % W8sel=strrep(W8sel,'icvf','ICVF');
% % % ppos        = [1442 341 594 456];
% % % 
% % %           
% % % for dns = 1:length(datanames)
% % %     metric      = datanames{dns};
% % % %     Dgrp        = A.grp.(metric)(isel,:);
% % %     Dgrp        = A.grp.charpath.(metric)(isel,:);
% % %     Dall        = cell(numel(isel),1);                                      % store all group level data for violin plot
% % %     pltstr      = [upper(strrep(metric,'_','-')) ' (group-level)'];
% % % 
% % %     % plot cdf
% % %     myfig([pltstr ': ' parcs{1}], ppos)
% % %     lnstyln             = 1;                                                % used to iterate through line styles
% % %     for ii = 1 : numel(isel)
% % %         dg              = Dgrp{ii,1};
% % %         if strncmpi(metric,'nodebe',6)
% % %             dg          = log10(dg(dg~=0));                                 % node betweenness is exp distributed
% % %         end
% % %         Dall{ii}        = dg;                                               % store for violin
% % %         h               = cdfplot(dg); hold on
% % %         h.LineWidth     = 3;
% % %         h.Color         = cmap(ii,:);
% % %         h.LineStyle     = lnstyl{lnstyln};
% % %         if lnstyln < 3
% % %             lnstyln     = lnstyln + 1;
% % %         else
% % %             lnstyln = 1;
% % %         end
% % %     end
% % %     title(''); %     title(pltstr)
% % %     legend(W8sel);
% % %     if strncmpi(metric,'nodebe',6)
% % %         xlabel(['log10-' upper(strrep(metric,'_','-'))]);
% % %     else
% % %         xlabel(upper(strrep(metric,'_','-'))); 
% % %     end
% % %     ylabel('cdf');
% % %     set(gca,'FontSize', fntsz);
% % %     
% % %     % Plot violin
% % %     T       = cell2table(Dall');
% % %     plotstr = [parcs{1} ', ' strin];
% % %     myfig(plotstr,[1442 465 608 332])
% % %     violinplot(T, 'ShowMean'); hold on
% % %     ylabel(upper(strrep(metric,'_','-')));
% % % %     title('Group-level node values')
% % %     set(gca, 'XTickLabel',W8sel, 'XTickLabelRotation',20, 'FontSize', fntsz)
% % % end
% % % 
% % % 
% % % 
% % % %% -------------------- KS statistic --------------------- %%
% % % % Strength, Clustering Coefficient, NodeBetweenness, Participation Coefficient 
% % % datanames   = {'strength' 'clusteringcoeff' 'nodebetweenness' ... 
% % %                'participationcoeff'};
% % % % datanames   = {'eccentricity'};
% % % 
% % %            
% % % isel    = [8 10 12];
% % % % isel    = [3 4 12];
% % % W8sel   = W8s(isel);
% % % ppos = [1442 341 594 456];
% % % 
% % %            
% % % for dns = 1:length(datanames)
% % %     metric      = datanames{dns};
% % %     Dgrp = A.grp.(metric)(isel,:);          Dsub = A.sub.(metric)(isel,:);
% % % %     Dgrp = A.grp.charpath.(metric)(isel,:); Dsub = A.sub.charpath.(metric)(isel,:);
% % %     ks2stat     = cell(size(Dsub));
% % %     pltstr      = upper(strrep(metric,'_','-'));
% % % 
% % %     % Get KS statistic
% % %     for jj = 1:J
% % %         for ii = 1:numel(isel)
% % %             dg                  = Dgrp{ii,jj};
% % %             ds                  = Dsub{ii,jj};
% % %             kst                 = nan(size(ds,1),1);
% % %             for ss = 1 : size(ds,1)
% % %                 [~,~,kst(ss)]   = kstest2(dg,ds(ss,:));                     % get ksstat for this subject vs group
% % %             end
% % %             ks2stat{ii,jj}      = kst; 
% % %         end
% % %     end
% % %     
% % %     % Plot KS statistic
% % %     for jj = 1:J
% % %         strtitl = {['\fontsize{18}' pltstr] ; ['\fontsize{10}' strin]};
% % % 
% % %         myfig([pltstr ': ' parcs{jj}], ppos)       
% % %         dplt    = num2cell(cell2mat(ks2stat(:,jj)'));
% % %         T       = cell2table(dplt);
% % %         violinplot(T, 'ShowMean'); hold on
% % %         ylabel('KS statistic');
% % %         title(strtitl)
% % %         ylim([0 1])
% % % 
% % %         % Rotate xlabels if they are too long
% % %         if any(cellfun(@(c) length(c),W8sel) > 10)
% % %             setopt = {'XTickLabel',W8sel,'XTickLabelRotation',25,'FontSize',fntsz};
% % %         else
% % %             setopt = {'XTickLabel',W8sel,'FontSize',fntsz};
% % %         end
% % %         set(gca, setopt{:});
% % %     end
% % % end
% % % 
% % % 
% % % 
% % % %% ------------ Violin subject & group ---------------- %%
% % % 
% % % % Modularity, Transitivity, Assortativity, Global Efficiency
% % % datanames       = {'modularity' 'transitivity' 'assortativity' 'efficiency_global'};
% % % % datanames       = {'lambda' 'efficiency' 'radius' 'diameter'};
% % % % datanames       = {'lambda' 'efficiency'};
% % % % datanames       = {'modularity' 'transitivity'};
% % % 
% % % 
% % % % isel            = [3 4 12];
% % % % isel            = [8 10 12];
% % % isel        = [1 2 3];         
% % % W8sel           = W8s(isel); 
% % % W8sel=strrep(W8sel,'(log)','');W8sel=strrep(W8sel,'-u','');W8sel=strrep(W8sel,'-s','');W8sel=strrep(W8sel,'-c','');W8sel=strrep(W8sel,'nos','NoS');
% % % W8sel=strrep(W8sel,'icvf','ICVF');
% % % ppos            = [1442 465 608 332];
% % % 
% % % for dns = 1:length(datanames)
% % %     metric      = datanames{dns};
% % % %     Dgrp        = A.grp.(metric)(isel,:);
% % % %     Dsub        = A.sub.(metric)(isel,:);
% % %     Dgrp        = A.grp.charpath.(metric)(isel,:);
% % %     Dsub        = A.sub.charpath.(metric)(isel,:);
% % %     pltstr      = 'Subject vs group';
% % % 
% % %     % plot
% % %     for jj = 1 : J
% % %         Qgt     = Dgrp(:,jj);
% % %         Qst     = num2cell(cell2mat(Dsub(:,jj)'));                          % convert to sub x weight cell array
% % %         T       = cell2table(Qst);
% % %         strtitl = {['\fontsize{18}' pltstr]};
% % %         myfig(strin,ppos)
% % %         V       = violinplot(T, 'ShowMean'); hold on
% % % 
% % %         % Add group level bar
% % %         for rr = 1:numel(Qgt)
% % %             V(rr).MeanPlot.YData        = [Qgt{rr} Qgt{rr}];
% % %             V(rr).MeanPlot.LineWidth    = 3;
% % %             V(rr).ShowMean              = 1;
% % %         end
% % %         ylabel(upper(strrep(metric,'_','-')));
% % % %         title(strtitl)
% % %         line(xlim,[0 0],'LineStyle', '--', 'Color', 'r')
% % %         set(gca,'XTickLabel',W8sel,'XTickLabelRotation',20,'FontSize',fntsz);
% % %     end
% % % end
% % % 
% % % 
% % % 
% % % % Mean strength
% % % isel            = [3 4 12];
% % % W8sel           = W8s(isel);
% % % ppos            = [1442 341 594 456];
% % % 
% % % metric      = 'strength';
% % % Dgrp        = A.grp.(metric)(isel,:);
% % % Dsub        = A.sub.(metric)(isel,:);
% % % pltstr      = 'Subject vs group';
% % % 
% % % % plot
% % % for jj = 1 : J
% % %     Dgt             = Dgrp(:,jj);
% % %     Dst             = Dsub(:,jj); 
% % %     for ii = 1 : numel(isel)
% % %         Dgt{ii,jj}  = mean(Dgt{ii,jj});                                     % compute mean over nodes
% % %         Dst{ii,jj}  = mean(Dst{ii,jj},2);
% % %     end
% % %         
% % %     Dst     = num2cell(cell2mat(Dst'));                                     % convert to sub x weight cell array
% % %     T       = cell2table(Dst);
% % %     strtitl = {['\fontsize{18}' pltstr]};
% % %     myfig(strin,ppos)
% % %     V       = violinplot(T, 'ShowMean'); hold on
% % %     % Add group level bar
% % %     for rr = 1:numel(Dgt)
% % %         V(rr).MeanPlot.YData        = [Dgt{rr} Dgt{rr}];
% % %         V(rr).MeanPlot.LineWidth    = 3;
% % %         V(rr).ShowMean              = 1;
% % %     end
% % %     ylabel('Mean node strength');
% % %     title(strtitl);
% % %     line(xlim,[0 0],'LineStyle', '--', 'Color', 'r')
% % % 
% % %     % Rotate xlabels if they are too long
% % %     if any(cellfun(@(c) length(c),W8sel) > 10)
% % %         setopt = {'XTickLabel',W8sel,'XTickLabelRotation',20,'FontSize',fntsz};
% % %     else
% % %         setopt = {'XTickLabel',W8sel,'FontSize',fntsz};
% % %     end
% % %     set(gca, setopt{:});
% % % end
% % % 
% % % 
% % % %% ------------------- Surfaces & Cross-subject variance --------------------- %%
% % % % Network properties
% % % % A.grp. { 'strength'  'clusteringcoeff'  'nodebetweenness'  'efficiency_local'  'participationcoeff' }
% % % % A.grp.charpath. { 'eccentricity' } 
% % % 
% % % % Setup
% % % nNode                       = [pinfo.n];
% % % conte69is                   = [pinfo.conte69i];
% % % parcs                       = {pinfo.name};
% % % parcs                       = strrep(parcs,'-',' ');
% % % [srf_lh,srf_rh]             = load_conte69();                               % load surface from brainspace
% % % 
% % % isel        = [1 2 3];         
% % % % isel    = [8 10 12];
% % % % isel        = [3 4 12];
% % % W8sel       = W8s(isel); 
% % % W8sel=strrep(W8sel,'(log)','');W8sel=strrep(W8sel,'-u','');W8sel=strrep(W8sel,'-s','');W8sel=strrep(W8sel,'-c','');W8sel=strrep(W8sel,'nos','NoS');
% % % W8sel=strrep(W8sel,'icvf','ICVF');
% % % metric      = 'clusteringcoeff';
% % % pltmetric  = upper(strrep(metric,'_','-'));
% % % drank       = nan(numel(isel),nNode);
% % % 
% % % % Rank data
% % % it      = 0;                                                                % just an iterator
% % % for ii = isel
% % %     it                  = it + 1;
% % %     d                   = A.grp.(metric){ii }';
% % % %     d                 = -1 * A.grp.charpath.(metric){ii}';
% % %     [drank(it,:), ~]  = tiedrank(d);                                              % rank data
% % % end
% % % 
% % % % Convert rankings to consistency measure
% % % drank(drank<350)    = 0;                                                    % focus on top 50 nodes
% % % drank(drank~=0)     = 1;                                                    % binarize top nodes
% % % drank               = sum(drank,1);                                         % sum over weights to get consistency ranking
% % % 
% % % % Convert to full surface
% % % dfull               = parcel2full(drank',conte69is(:,1));
% % % 
% % % % Plot on surface
% % % h=plot_hemispheres(dfull,{srf_lh,srf_rh}, 'labeltext',pltmetric);
% % % h.figure.NumberTitle='off';
% % % h.figure.Name=[pltmetric ' top 50 ranked nodes: ' strin ', ' parcs{1}];
% % % % h.figure.Name=[pltmetric ' reverse ranked data: ' strin ', ' parcs{1}];
% % % h.cb.Ticks  = [1 2 3];
% % % modmap      = zeros(4,3);           % custom colormap
% % % modmap(1,:) = [0.2081 0.1663 0.5292]; % dark blue
% % % modmap(2,:) = [0 1 1];              % cyan
% % % modmap(3,:) = [1 1 0];              % yellow
% % % modmap(4,:) = [1 0 0];              % red
% % % for ax =1:4
% % %     colormap(h.axes(ax),modmap)
% % % end
% % % 
% % % % Set common color mapping
% % % % cmax = max(cellfun(@(c) max(c), dfull(subsel)));
% % % % cmin = min(cellfun(@(c) min(c), dfull(subsel)));
% % % cmax = 400; cmin = 350;
% % % for ax =1:4
% % %     h.axes(ax,1).CLim=[cmin, cmax];
% % %     h.axes(ax,2).CLim=[cmin, cmax];
% % % end
% % % h.cb(1).Ticks = [cmin, cmax];
% % % h.cb(2).Ticks = [cmin, cmax];
% % % 
% % % 
% % % % A.grp. { 'strength'  'clusteringcoeff'  'nodebetweenness'  'efficiency_local'  'participationcoeff' }
% % % % A.grp.charpath. { 'eccentricity' } 
% % % % Cross-subject variance on surface
% % % isel        = [3 4 12];
% % % W8sel       = W8s(isel);
% % % metric      = 'strength';
% % % pltmetric   = upper(strrep(metric,'_','-'));
% % % dfull       = cell(numel(isel),1);
% % % dvar        = nan(400,numel(isel));
% % % Dsub        = A.sub.(metric)(isel);
% % % % Dsub        = A.sub.(metric);
% % % 
% % % 
% % % % Convert nodal data to full surface
% % % for ii = 1 : numel(isel)
% % %     d           = Dsub{ii};
% % %     dvar(:,ii)  = std(d).^2;                                                % compute cross-subject variance at each node
% % %     dfull{ii}   = parcel2full(dvar(:,ii),conte69is(:,1));
% % % end
% % % 
% % % % Plot nodal data on surface
% % % subsel = [1 2];
% % % D1=dfull{subsel(1),1};
% % % D2=dfull{subsel(2),1};
% % % h=plot_hemispheres([D1,D2],{srf_lh,srf_rh}, 'labeltext',{W8sel{subsel(1)}, W8sel{subsel(2)}});
% % % h.figure.NumberTitle='off';
% % % h.figure.Name=['Cross-subject variance in ' pltmetric ': ' strin ', ' parcs{1}];
% % % 
% % % subsel = [3 3];
% % % D1=dfull{subsel(1),1};
% % % D2=dfull{subsel(2),1};
% % % h=plot_hemispheres([D1,D2],{srf_lh,srf_rh}, 'labeltext',{W8sel{subsel(1)}, W8sel{subsel(2)}});
% % % h.figure.NumberTitle='off';
% % % h.figure.Name=['Cross-subject variance in ' pltmetric ': ' strin ', ' parcs{1}];
% % % 
% % % % violin plot of variance
% % % T       = array2table(dvar);
% % % plotstr = [parcs{1} ', ' strin];
% % % myfig(plotstr,[1442 465 608 332])
% % % violinplot(T, 'ShowMean'); hold on
% % % ylabel('Cross-subject variance');
% % % title(pltmetric)
% % % set(gca, 'XTickLabel',W8sel, 'XTickLabelRotation',20, 'FontSize', fntsz)
% % % 
% % % 
% % % % % % Box plots of ALL NODES cross-subject variance of nodal properties
% % % % % metric = 'clusteringcoeff';
% % % % % isel = [2 3 4 8 10 12];
% % % % % Dsub = A.sub.(metric);
% % % % % manlabels = {'50' '100' '150' '200' '250' '300' '350' '400'};
% % % % % 
% % % % % for ii = isel
% % % % %     d = Dsub{ii}; 
% % % % %     myfig, 
% % % % %     boxplot(d) 
% % % % %     title(['Subject distribution of nodal ' upper(metric) ': ' W8s{ii}]); 
% % % % %     xlabel('nodes'); ylabel(upper(metric));
% % % % %     set(gca,'XTick',50:50:400, 'XTickLabel',manlabels,  'XTickLabelRotation',90, 'FontSize', fntsz)
% % % % % end
% % % 
% % % 
% % % 
% % % %% ----------------- Matrices ------------------ %%
% % % isel            = [8 10 12];
% % % % EdgeBetweenness (group)
% % % for it = isel
% % %     strt = ['Edge Betweeness (group): ' W8s{it}]; 
% % %     connplot(A.grp.edgebetweenness{it}, strt, pinfo)
% % % end
% % % 
% % % % Search Information (group)
% % % for it = isel
% % %     strt = ['Search Information (group): ' W8s{it}]; 
% % %     connplot(A.grp.searchinformation{it}, strt, pinfo)
% % % end
% % % 
% % % % Path Transitivity (group)
% % % for it = isel
% % %     strt = ['Path Transitivity (group): ' W8s{it}]; 
% % %     connplot(A.grp.pathtransitivity{it}, strt, pinfo)
% % % end
% % % 
% % % % Correlate metrics with FC
% % % % {'searchinformation' 'pathtransitivity' 'edgebetweenness'}
% % % metric = 'searchinformation';
% % % % isel            = [3 4 12];
% % % % isel            = [8 10 12];
% % % isel            = [1 2 3];         
% % % W8sel           = W8s(isel);
% % % W8sel=strrep(W8sel,'(log)','');W8sel=strrep(W8sel,'-u','');W8sel=strrep(W8sel,'-s','');W8sel=strrep(W8sel,'-c','');W8sel=strrep(W8sel,'nos','NoS');
% % % W8sel=strrep(W8sel,'icvf','ICVF');
% % % Dgrp            = A.grp.(metric)(isel);
% % % Dsub            = A.sub.(metric)(isel);
% % % % ifc             = cellstrfind(W8s, 'FC');
% % % Dsfc            = D{ifc};
% % % Dgfc            = groupavg(Dsfc,3,'nz');
% % % Sn              = size(Dsfc,3);
% % % Rg              = nan(1,numel(isel));
% % % Rs              = nan(Sn,numel(isel));
% % % [it,cil,cirm]   = conn_useCIs(Dgfc,pinfo);                                  % Get node order for communities from FC
% % % Dgfc            = Dgfc(it,it);                                              % reorder FC data
% % % 
% % % % plot connectivity matrices
% % % for ii = 1 : numel(isel)
% % %     strt                        = ['FC vs Path Transitivity (' W8sel{ii} ')'];
% % %     Dg                          = Dgrp{ii}(it,it);
% % %     Dg(tril(ones(400),-1)~=0)   = Dgfc(tril(ones(400),-1)~=0);
% % %     myfig(strt);
% % %     imagesc(Dg); axis square; colorbar; title(strt); set(gca,'FontSize',20)
% % %     set(gca,'XTick',ceil(cirm'),'XTickLabel',cil','XTickLabelRotation',90)
% % %     set(gca,'YTick',ceil(cirm'),'YTickLabel',cil')
% % % end
% % % 
% % % % group-level
% % % Dgfc            = groupavg(Dsfc,3,'nz');
% % % for ii = 1 : numel(isel)
% % %     Dg  = Dgrp{ii};
% % %     connplot(Dg,[upper(metric) ': ' W8s{isel(ii)}],pinfo);
% % %     r1 = conn_corrmat(Dg,Dgfc);                                             % correlate lower tri with FC     
% % %     r2 = conn_corrmat(Dg',Dgfc);                                            % and upper tri
% % %     Rg(ii) = mean([r1 r2]);                                                 % average r of both tris
% % % end
% % % % connplot(Dgfc,'Group-level FC',pinfo);
% % % 
% % % plot_conn_vsw8([{Dgfc};Dgrp(1)], {W8s{1}(1:2), W8sel{1}}, 'fc')             % heat scatter vs FC
% % % set(gcf,'Position',[1442         377         568         420])
% % % plot_conn_vsw8([{Dgfc};Dgrp(2)], {W8s{1}(1:2), W8sel{2}}, 'fc')             % heat scatter vs FC
% % % set(gcf,'Position',[1442         377         568         420])
% % % plot_conn_vsw8([{Dgfc};Dgrp(3)], {W8s{1}(1:2), W8sel{3}}, 'fc')             % heat scatter vs FC
% % % set(gcf,'Position',[1442         377         568         420])
% % % 
% % % set(gcf,'Name', ['FC vs' upper(metric) ': group-level'])
% % % 
% % % % subject-level
% % % for ii = 1 : numel(isel)
% % %     Ds              = Dsub{ii};
% % %     for ss = 1 : Sn
% % %         ds          = Ds(:,:,ss);
% % %         dsfc        = Dsfc(:,:,ss);
% % %         r1          = conn_corrmat(ds,dsfc);                                % correlate lower tri with FC     
% % %         r2          = conn_corrmat(ds',dsfc);                               % and upper tri
% % %         Rs(ss,ii)   = mean([r1 r2]);                                        %      
% % %     end
% % % end
% % % 
% % % T       = array2table(Rs);
% % % plotstr = [parcs{1} ', ' strin];
% % % myfig(plotstr,[1442 465 608 332])
% % % V       = violinplot(T, 'ShowMean'); hold on
% % % % Add group level bar
% % % for rr = 1:numel(isel)
% % %     V(rr).MeanPlot.YData        = [Rg(rr) Rg(rr)];
% % %     V(rr).MeanPlot.LineWidth    = 3;
% % %     V(rr).ShowMean              = 1;
% % % end
% % % ylabel('Linear correlation with FC');
% % % title({['\fontsize{20}' 'Path Transitivity from SC'] ; ['\fontsize{12}' '(subject & group)']})
% % % % title({['\fontsize{20}' 'Search Information from SC'] ; ['\fontsize{12}' '(subject & group)']})
% % % line(xlim,[0 0],'LineStyle', '--', 'Color', 'r')
% % % set(gca, 'XTickLabel',W8sel, 'XTickLabelRotation',20, 'FontSize', fntsz)
% % % 
% % % 
% % % 
% % % % Difference of search information
% % % metric = 'searchinformation';
% % % % isel            = [3 4 12];
% % % % isel            = [8 10 12];
% % % isel            = [1 2 3];         
% % % W8sel           = W8s(isel);
% % % W8sel=strrep(W8sel,'(log)','');W8sel=strrep(W8sel,'-u','');W8sel=strrep(W8sel,'-s','');W8sel=strrep(W8sel,'-c','');W8sel=strrep(W8sel,'nos','NoS');
% % % W8sel=strrep(W8sel,'icvf','ICVF');
% % % Dgrp            = A.grp.(metric)(isel);
% % % Dsub            = A.sub.(metric)(isel);
% % % 
% % % ii1=3; ii2=2;
% % % test = Dgrp{ii1} - Dgrp{ii2};
% % % strtmp = ['Difference SI ' W8sel{ii1} ' - ' W8sel{ii2} ' (group)'];
% % % connplot(test,strtmp,pinfo)
% % % 
% % % ii1=3; ii2=1;
% % % test = Dgrp{ii1} - Dgrp{ii2};
% % % strtmp = ['Difference SI ' W8sel{ii1} ' - ' W8sel{ii2} ' (group)'];
% % % connplot(test,strtmp,pinfo)
% % % 
% % % ii1=2; ii2=1;
% % % test = Dgrp{ii1} - Dgrp{ii2};
% % % strtmp = ['Difference SI ' W8sel{ii1} ' - ' W8sel{ii2} ' (group)'];
% % % connplot(test,strtmp,pinfo)
% % % 
% % % 
% % % 
% % % 
% % % 
% % % %% -- testing
% % % highcclim = .6;
% % % lowcclim = .33;
% % % 
% % % ccg_commit      = A.grp.clusteringcoeff{12};
% % % myfig, histogram(ccg_commit)
% % % title('clustering coeff: commit-sc group'); set(gca, 'FontSize', 20)
% % % line([highcclim highcclim],ylim,'LineStyle', '--', 'Color', 'r')
% % % line([lowcclim lowcclim],ylim,'LineStyle', '--', 'Color', 'r')
% % % 
% % % 
% % % % get indices for high & low clustering coeff nodes
% % % highccnodes     = find(ccg_commit>highcclim);
% % % lowccnodes      = find(ccg_commit<lowcclim);
% % % 
% % % % compare other metrics for these nodes
% % % mean(A.grp.strength{12}(highccnodes))                   % high cc nodes are much weaker
% % % mean(A.grp.strength{12}(lowccnodes))                    
% % % 
% % % mean(A.grp.nodebetweenness{12}(highccnodes))            % high cc nodes have much lower betweenness centrality
% % % mean(A.grp.nodebetweenness{12}(lowccnodes))
% % % 
% % % mean(A.grp.participationcoeff{12}(highccnodes))            % high cc nodes have much lower betweenness centrality
% % % mean(A.grp.participationcoeff{12}(lowccnodes))
% % % 
% % % % Compare the commit weights for these nodes
% % % dhcc            = D{12}(highccnodes,:,:);
% % % dlcc            = D{12}(lowccnodes,:,:);
% % % mean(dhcc(dhcc~=0))                                     % edges of high cc nodes have higher commit weights
% % % mean(dlcc(dlcc~=0))
% % % cmap            = colormap(hsv(2));
% % % myfig('COMMIT-weight of edges connecting nodes with high vs low clustering coeff')
% % % histogram(dhcc(dhcc~=0),100,'FaceAlpha',.2,'FaceColor',cmap(1,:)); hold on;
% % % histogram(dlcc(dlcc~=0),100,'FaceAlpha',.2,'FaceColor',cmap(2,:));
% % % title('COMMIT-weight of edges connecting nodes with high vs low clustering coeff')
% % % legend({'high cc' 'low cc'})
% % % set(gca, 'FontSize', 20)
% % % 
% % % dhcclos            = D{2}(highccnodes,:,:);
% % % dlcclos            = D{2}(lowccnodes,:,:);
% % % mean(dhcclos(dhcclos~=0))
% % % mean(dlcclos(dlcclos~=0))
% % % myfig('Length of edges connecting nodes with high vs low clustering coeff')
% % % histogram(dhcclos(dhcclos~=0),100,'FaceAlpha',.2,'FaceColor',cmap(1,:)); hold on;
% % % histogram(dlcclos(dlcclos~=0),100,'FaceAlpha',.2,'FaceColor',cmap(2,:));
% % % title('Length of edges connecting nodes with high vs low clustering coeff')
% % % legend({'high cc' 'low cc'})
% % % set(gca, 'FontSize', 20)
% % % 
% % % length(find(dhcc~=0)) / numel(highccnodes) / 50 / 399           % have much lower density in COMMIT
% % % length(find(dlcc~=0)) / numel(lowccnodes) / 50 / 399 
% % % 
% % % dtmphcc = D{3}(highccnodes,:,:);
% % % dtmplcc = D{3}(lowccnodes,:,:);
% % % length(find(dtmphcc~=0)) / numel(highccnodes) / 50 / 399           % have much lower density
% % % length(find(dtmplcc~=0)) / numel(lowccnodes) / 50 / 399 
% % % 
% % % dtmphcc = D{4}(highccnodes,:,:);
% % % dtmplcc = D{4}(lowccnodes,:,:);
% % % length(find(dtmphcc~=0)) / numel(highccnodes) / 50 / 399           % have much lower density
% % % length(find(dtmplcc~=0)) / numel(lowccnodes) / 50 / 399 



%--------------------------------------------------------------------------
end