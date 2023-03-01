function [A] = P_conn_gt_hubConsistency(pinfo,varargin)
%  Visualizes node centrality measures computed on weighted & binary
%  connectivity data. 
%  Summarizes across measures --> Hub Consistency 
%
% Input:
%                           + Required +
%   pinfo   : 1xJ structure with parcellation info (see conn_getParcInfo.m)
%             necessary for surface plots
%
%                           + Optional +
%   A       : Structure storing network attributes
%   path    : full file path to saved A structure
%   str     : character vector for plot title  
%   cmap    : Ix3 colormap
%
% Output:
%   A       : 1xJ structure
%
% 2022 Mark C Nelson, McConnell Brain Imaging Centre, MNI, McGill
%--------------------------------------------------------------------------
%% Optional inputs
nin                 = max(nargin,1) - 1;
defaults            = {[],[],[],[]};
[A,path,str,cmap]   = INhandler(varargin,nin,defaults);
keyboard
%% SETUP

% Load A if necessary
if isempty(A)
    if ~isempty(path)
        load(path)
    else
        error('* * Must supply either A or PATH input')
    end
end

% Get info from A
AW8s                    = A.info.weights;
Aparcs                  = A.info.parcs;
I                       = numel(AW8s);
J                       = numel(Aparcs);
% Isel                        = cell2mat(cellstrfind(AW8s,{'nos' 'sift2' 'commit' 'R1' 'ICVF' 'fa' 'rd'})); % reordering
Nnode                   = unique(cellfun(@(c) length(c),A.grp.wei.strength),'rows');
centrality_measures     = {'strength' 'closeness' 'betweenness' 'eigenvector'}; % Centrality measures of interest
Ncms                    = numel(centrality_measures);

% Plot settings
if isempty(cmap); cmap  = colormap(hsv(I)); end
fntsz                   = 25;
ppos                    = [-2560 -82 2560 1290];                            % monitor


%%---------------------------- GROUP-LEVEL ------------------------------%%


%% DEFINE HUBS: Method 1
% 1. Nodes ranked on centrality measures
% 2. Arbitrary threshold identifies top ranked nodes in each metric
% 3. Arbitrary consistency measure identifies hubs as top ranked nodes across minimum fraction of centrality measures

% THRESHOLD               = .05;                                              % Percent of network nodes to treat as hubs within each centrality measure
% CONSISTENCY             = .75;                                              % Percent consistency across centrality measures required for network HUB    
% Nhub                    = ceil(THRESHOLD.*Nnode);
% 
% 
% HUBS_g                  = cell(I,J);
% 
% for jj = 1 : J
%     Nhubt               = Nhub(jj);
%     
%     for ii = 1 : I
%         hubs            = zeros(Ncms,Nhubt);
%         
%         % Rank & Threshold Centrality measures
%         for cc = 1 : Ncms
%             CM_t        = centrality_measures{cc};
%             dt          = A.grp.wei.(CM_t){ii,jj};
%             dt_rank     = tiedrank(-1.*dt);
%             hubs(cc,:)  = find(dt_rank<=Nhubt);                             % Would be better to sort nodes by ranking too
%         end
%         
%         % Reduce to indciated % consistency across centrality measures
%         hubs_superset   = unique(hubs);
%         HUBS_g_t        = [];
%         for hh = 1 : length(hubs_superset)
%             hubt        = hubs_superset(hh);
%            if length(find(hubs==hubt)) / Ncms >= CONSISTENCY
%                HUBS_g_t = [HUBS_g_t hubt]; 
%            end
%         end
%         
%         % Store Hubs for this Weight
%         HUBS_g{ii,jj}   = HUBS_g_t;
%     
%     end
%     
%     % Compute Consistency Across Weights???
%     
%     
%     % Plot Each Surface Separately???
%     
%     
% end



%% DEFINE HUBS: Method 2 
% 1. Hubs = node with at least 1 centrality measure which is at least 2 SDs > mean over all regions
% (Bassett et al. (2008) J. of Neuroscience)

% THR                     = 2;                                                % Number of standard deviations above mean to threshold
% HUBS_g                  = cell(I,J);
% 
% for jj = 1 : J
%     
%     for ii = 1 : I
%         hubs            = [];
%         
%         % Identify nodes with centrality at least 2STD > mean
%         for cc = 1 : Ncms
%             CM_t        = centrality_measures{cc};
%             dt          = A.grp.wei.(CM_t){ii,jj};
%             dt_mean     = mean(dt);
%             dt_std      = std(dt);
%             cutoff      = dt_mean + (THR*dt_std);                           % Nodes above this value will be labeled hubs
%             hubs        = [hubs find(dt >= cutoff)];                        
%             
% % % %             % TEST PLOTS
% % % %             ppos=[-2559 927 2560 281];
% % % %             myfig(['Node Centrality: ' AW8s{ii}],ppos);
% % % %             plot(1:Nnode(jj), A.grp.wei.(centrality_measures{cc}){ii,jj},'LineWidth',3); ylabel(centrality_measures{cc}); set(gca,'FontSize',fntsz); title(AW8s{ii});
% % % %             L1 = line(xlim,[dt_mean dt_mean],'Color','k','LineStyle','--','LineWidth',3);
% % % %             L2 = line(xlim,[(dt_mean+(THR*dt_std))  (dt_mean+(THR*dt_std))],'Color','r','LineStyle','--','LineWidth',3);
% % % %             legend([L1 L2],{'Mean' [num2str(THR) '*STD']});
%         end
%         
%         % Store hubs for this weight
%         HUBS_g{ii,jj}   = unique(hubs);
%     end
%     
%     
%     
%     % Plot Surfaces
%     Dplt                = zeros(Nnode(jj),I);                               % plot_conn_surf.m requires node x dataset double to plot 
%     for ii = 1 : I; Dplt(HUBS_g{ii},ii) = 1; end                            % 1 = hub node
%     titlestr            = 'Hubs (group)';
%     [Hcx,Hsx]           = plot_conn_surf(Dplt,pinfo(jj),'both',AW8s, titlestr);
%     
% 
%     % Compute Consistency Across Weights
%     D_sum               = sum(Dplt,2);
%     D_consistency       = D_sum./I;                                         % Hub frequency across datasets normalized by I
%     titlestr            = 'Hub consistency (group)';
%     [Hcx,Hsx]           = plot_conn_surf(D_consistency,pinfo(jj),'both',{'hub f / Wn'}, titlestr);
%     
%     % Mod colormaps & colorbars
%     AXall           = [];
%     CBall           = [];
%     for aa = 1 : length(Hcx);
%         AXall       = [AXall; Hcx{aa}.axes(:)];                             % cortical handles
%         CBall       = [CBall; Hcx{aa}.cb(:)];
%     end
%     for aa = 1 : length(Hsx)
%         AXall       = [AXall; Hsx{aa}.axes(:)];                             % subcortical handles
%         CBall       = [CBall; Hsx{aa}.cb(:)];
%     end
%     AXall           = unicb(AXall);                                         % Standardize all color axes
%     AXall           = unicmap(AXall,'inferno');                             % Change colormaps
%     
%     % Fix colorbars
%     clims           = AXall(1).CLim;
%     for cc = 1 : length(CBall)
%         CBall(cc).Ticks     = clims;
%         CBall(cc).FontSize  = fntsz;
%     end  
% 
%     
% end



%% DEFINE HUBS: Method 3
% 1. "Hubness" computed for all nodes (score of 0-4)
%   + 1 for top 20% strength
%   + 1 for top 20% centrality (betweenness)
%   + 1 for top 20% closeness (aka lowest 20% path length)
%   + 1 for lowest 20% clustering coefficient
% 2. Hubs = nodes with Hubness >= 2
% (Van den Heuvel et al. (2010) J. of Neuroscience)

THR                         = 80;                                           % Arbitrary threshold used to identify hubs
SCORES                      = cell(I,J);

for jj = 1 : J
    
    Nnodet                  = Nnode(jj);
        
    for ii = 1 : I
        
        scores_total        = zeros(1,Nnodet);                                  % Empty vector to keep track of overall nodal scores
        
        % Identify central nodes: +1 for each measure if in top 20%
        for cc = 1 : Ncms
            CM_t            = centrality_measures{cc};
            dt              = A.grp.wei.(CM_t){ii,jj};
            Ihubs           = find(dt > prctile(dt,THR));
            hubs            = zeros(1,Nnodet);                              % empty vector to assign scores for this centrality measure
            hubs(Ihubs)     = 1;                                            % 1 point assigned for this measure
            scores_total    = scores_total + hubs;                          % Add these scores to overall node scores
        end
        
        % Identify nodes with low clustering: +1 if in bottom 20%
        dt                  = A.grp.wei.clustering{ii,jj};
        Ihubs               = find(dt < prctile(dt,100-THR));
        hubs                = zeros(1,Nnodet);                              % empty vector to assign scores for this centrality measure
        hubs(Ihubs)         = 1;                                            % 1 point assigned for this measure
        scores_total        = scores_total + hubs;                          % Add these scores to overall node scores
                
        
        % Store Nodal scores for this Weight
        SCORES{ii,jj}       = scores_total;
    
    end
    
    
    
    % ALL SURFACES
    titlestr            = 'Hubs (group)';
    Dplt                = cell2mat(SCORES)';
    [Hcx,Hsx]           = plot_conn_surf(Dplt,pinfo(jj),'both',AW8s, titlestr);
   
    % Mod colormaps & colorbars
    AXall           = [];
    CBall           = [];
    for aa = 1 : length(Hcx);
        AXall       = [AXall; Hcx{aa}.axes(:)];                             % cortical handles
        CBall       = [CBall; Hcx{aa}.cb(:)];
    end
    for aa = 1 : length(Hsx)
        AXall       = [AXall; Hsx{aa}.axes(:)];                             % subcortical handles
        CBall       = [CBall; Hsx{aa}.cb(:)];
    end
    AXall           = unicb(AXall);                                         % Standardize all color axes
    
        modmap      = zeros(6,3);           % custom colormap
    modmap(1,:) = [.95 .95 .95];           % very faint gray
    modmap(2,:) = [0 1 1];              % cyan
    modmap(3,:) = [0 1 0];              % green
    modmap(4,:) = [0 0 1];              % blue
    modmap(5,:) = [1 1 0];              % yellow
    modmap(6,:) = [1 0 0];              % red
    for xx = 1 : length(AXall)
        colormap(AXall(xx),modmap);
    end
    
    % Fix colorbars
%     clims           = AXall(1).CLim;
    for cc = 1 : length(CBall)
        CBall(cc).Ticks     = [];   % clims
        CBall(cc).FontSize  = fntsz;
    end
    
    
    
    % AVERAGE SURFACE
    Iavg                = [1 5];                                            % Compute average across these networks
    D_sum               = sum(Dplt(:,Iavg),2);
    D_avg               = ceil(D_sum./length(Iavg));                        % Hub frequency across datasets normalized by number of datasets
    w8str               = AW8s{Iavg(1)}; 
    for xx = 2 : length(Iavg) 
        w8str = [w8str ', ' AW8s{Iavg(xx)}]; 
    end
    titlestr            = ['AVERAGE Hubs (group): ' w8str];
    [Hcx,Hsx]           = plot_conn_surf(D_avg,pinfo(jj),'both',{'hub f / Wn'}, titlestr);

    % Mod colormaps & colorbars
    AXall=[]; CBall=[];
    for aa=1:length(Hcx); AXall=[AXall; Hcx{aa}.axes(:)]; CBall=[CBall; Hcx{aa}.cb(:)]; end         % cortical handles
    for aa = 1 : length(Hsx); AXall=[AXall; Hsx{aa}.axes(:)]; CBall=[CBall; Hsx{aa}.cb(:)]; end     % subcortical handles
    for xx = 1 : length(AXall); colormap(AXall(xx),modmap); end                                     % Change colormaps
    for cc=1:length(CBall); CBall(cc).Ticks=[]; CBall(cc).FontSize=fntsz; end                       % Fix colorbars
    
    
    
    % Hub DISTANCE
    Idist                   = [7 2 1 5 4 6 3];                              % Compute consistency across these networks
%     Idist                   = [1 2];                              % Compute consistency across these networks
    Ndist                   = length(Idist);
    Ddist                   = Dplt(:,Idist);
    hubdist                 = pdist(Ddist','euclidean');
    Npwdist                 = length(hubdist);
    hubdistmat              = zeros(Ndist);
    Ilthubdist              = logical(tril(ones(Ndist),-1));
    hubdistmat(Ilthubdist)  = hubdist;
    myfig; imagesc(hubdistmat); axis square; colorbar; cm=red2bluecmap(gca,0);
    set(gca,'XTickLabel',AW8s(Idist),'XTickLabelRotation',40); set(gca,'YTickLabel',AW8s(Idist),'YTickLabelRotation',40); set(gca,'FontSize',20);
    title('Group Hub Score Distance (Euclidean)');

    
    D_consistency       = ceil(D_sum./length(Icnstncy));                    % Hub frequency across datasets normalized by number of datasets
    w8str               = AW8s{Icnstncy(1)}; 
    for xx = 2 : length(Icnstncy) 
        w8str = [w8str ', ' AW8s{Icnstncy(xx)}]; 
    end
    titlestr            = ['CONSISTENCY Hubs (group): ' w8str];
    [Hcx,Hsx]           = plot_conn_surf(D_consistency,pinfo(jj),'both',{'hub f / Wn'}, titlestr);
    
    % Mod colormaps & colorbars
    AXall=[]; CBall=[];
    for aa=1:length(Hcx); AXall=[AXall; Hcx{aa}.axes(:)]; CBall=[CBall; Hcx{aa}.cb(:)]; end         % cortical handles
    for aa = 1 : length(Hsx); AXall=[AXall; Hsx{aa}.axes(:)]; CBall=[CBall; Hsx{aa}.cb(:)]; end     % subcortical handles
    for xx = 1 : length(AXall); colormap(AXall(xx),modmap); end                                     % Change colormaps
    for cc=1:length(CBall); CBall(cc).Ticks=[]; CBall(cc).FontSize=fntsz; end                       % Fix colorbars
    
end

%--------------------------------------------------------------------------
end