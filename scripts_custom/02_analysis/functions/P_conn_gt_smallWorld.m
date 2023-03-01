function P_conn_gt_smallWorld(varargin)
% Small worldness is plot for all group-level weighted networks.
% 
% Requires the network features CLUSTERING COEFFICINT and CHARACTERISTIC 
% PATH LENGTH computed in both empirical and null networks.
% (see conn_gt_smallWorld.m & conn_gt_all.m)
%
% Null network measures are averaged and used to normalize empirical 
% measures.
%
% Input:
%                           + Optional +
%   A       : Structure storing network attributes
%   path    : full file path to saved A structure
%   str     : character vector for plot title  
%   cmap    : Ix3 colormap
%
% 2022 Mark C Nelson, McConnell Brain Imaging Centre, MNI, McGill
%--------------------------------------------------------------------------

%% Optional inputs
nin                 = max(nargin,1);
defaults            = {[],[],[],[]};
[A,path,str,cmap]   = INhandler(varargin,nin,defaults);

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
AW8s                        = A.info.weights;
Aparcs                      = A.info.parcs;
Astr                        = A.info.details;
I                           = numel(AW8s);
J                           = numel(Aparcs);

% Indicate which weights to plot
Isel                        = cell2mat(cellstrfind(AW8s,AW8s));             % All


% Plot settings
idlineopt                   = {'Color','k','LineStyle',':','LineWidth',2};  % For identity line
if isempty(cmap); cmap      = colormap(hsv(I)); end
fntsz                       = 25;
% ppos                        = [-2559 670 659 538];                        % large
ppos                        = [-2559 841 414 367];                          % small
mrkrs                       = {'s' 's' 's' '.' '.' '.' '.' 's' }; 
mrkrsize                    = [25 25 25 100 100 100 100 25]; 



%---------------------------- GROUP-LEVEL ------------------------------%

%% Normalize network measures

% Char Path
metric                  = 'charpath';
NULL_raw                = A.null.(metric);
EMPR_raw                = A.grp.(metric);
EMPR_nrm                = cell(size(EMPR_raw));

for jj = 1 : J
    for ii = 1 : I
        
        null_raw_tmp    = NULL_raw{ii,jj};
        null_raw_tmp    = mean(null_raw_tmp);                               % Measures normalized by mean across null networks
        
        emp_raw_tmp     = EMPR_raw{ii,jj};
        EMPR_nrm{ii,jj} = emp_raw_tmp/null_raw_tmp;                         % Normalize
    end
end
CharPath_nrm            = cell2mat(EMPR_nrm);


% Clustering Coeff
metric                  = 'clusteringcoeff';
NULL_raw                = A.null.(metric);
EMPR_raw                = A.grp.(metric);
EMPR_nrm                = cell(size(EMPR_raw));

for jj = 1 : J
    for ii = 1 : I
        
        null_raw_tmp    = NULL_raw{ii,jj};
        null_raw_tmp    = mean(null_raw_tmp,1);                             % Measures normalized by mean across null networks WITHIN NODE
        
        emp_raw_tmp     = EMPR_raw{ii,jj};
        EMPR_nrm{ii,jj} = emp_raw_tmp./null_raw_tmp;                        % Normalize each node separately
        EMPR_nrm{ii,jj} = nanmean(EMPR_nrm{ii,jj});                         % Mean over nodes
    end
end
ClstrCof_nrm            = cell2mat(EMPR_nrm);

SmllWrld_nrm_g          = ClstrCof_nrm./CharPath_nrm;


%% Plot

disp('* ------------------ PLOTTING: GROUP Small-Worldness -----------------*')

for jj = 1 : J
    Dx                  = CharPath_nrm(:,jj);
    Dy                  = ClstrCof_nrm(:,jj);
    pltlims             = [floor(min([Dx; Dy])) ceil(max([Dx; Dy]))];
    pltstr              = ['Normalized Small-Worldness (group, ' Aparcs{jj} ')'];
    myfig(pltstr,ppos);
    
    for ii = Isel
        plot(Dx(ii),Dy(ii), 'Marker',mrkrs{ii}, 'MarkerSize',mrkrsize(ii), 'MarkerFaceColor',cmap(ii,:), 'MarkerEdgeColor',cmap(ii,:)); hold on
    end
    
%     xlabel('Characteristic Path Length'); ylabel('Mean Clustering Coeff'); title('Normalized Small-Worldness: Group'); 
    ylabel('C/C_n_u_l_l'); xlabel('L/L_n_u_l_l'); title('Small-worldness');
    set(gca,'FontSize',fntsz); 
    xlim(pltlims); ylim(pltlims); 
    
    RL=refline(1,0); set(RL,idlineopt{:});                                  % Identity line (works in some plots...)
    legend(AW8s(Isel))
    
end



%---------------------------- SUBJECT-LEVEL ------------------------------%


%% Normalize network measures

% Char Path
metric                  = 'charpath';
NULL_raw                = A.null.(metric);
EMPR_raw                = A.sub.(metric);
EMPR_nrm                = cell(size(EMPR_raw));

for jj = 1 : J
    for ii = 1 : I
        
        null_raw_tmp    = NULL_raw{ii,jj};
        null_raw_tmp    = mean(null_raw_tmp);                               % Measures normalized by mean across null networks
        
        emp_raw_tmp     = EMPR_raw{ii,jj};
        EMPR_nrm{ii,jj} = emp_raw_tmp/null_raw_tmp;                         % Normalize
    end
end
CharPath_nrm            = EMPR_nrm;


% Clustering Coeff
metric                  = 'clusteringcoeff';
NULL_raw                = A.null.(metric);
EMPR_raw                = A.sub.(metric);
EMPR_nrm                = cell(size(EMPR_raw));

for jj = 1 : J
    for ii = 1 : I
        
        null_raw_tmp    = NULL_raw{ii,jj};
        null_raw_tmp    = mean(null_raw_tmp,1);                             % Measures normalized by mean across null networks WITHIN NODE
        
        emp_raw_tmp     = EMPR_raw{ii,jj};
        subarray        = zeros(size(emp_raw_tmp));
        
        for ss = 1 : size(emp_raw_tmp,1)
            subarray(ss,:)  = emp_raw_tmp(ss,:)./null_raw_tmp;              % Normalize each node separately
        end
            
        EMPR_nrm{ii,jj} = nanmean(subarray,2);                              % Mean over nodes
    end
end
ClstrCof_nrm            = EMPR_nrm;


%% Plot

disp('* ------------------ PLOTTING: SUBJECT Small-Worldness -----------------*')

% Scatter
for jj = 1 : J
    Dx                  = CharPath_nrm(:,jj);
    Dy                  = ClstrCof_nrm(:,jj);
    pltlims             = [floor(min(cellfun(@(c) min(c), [Dx; Dy]))) ceil(max(cellfun(@(c) max(c), [Dx; Dy])))];
    pltstr              = ['Normalized Small-Worldness (Subject, ' Aparcs{jj} ')'];
    myfig(pltstr,ppos);
    
    for ii = Isel
        plot(Dx{ii},Dy{ii}, 'k.', 'MarkerSize',25, 'MarkerEdgeColor',cmap(ii,:)); hold on
    end
    
%     xlabel('Characteristic Path Length'); ylabel('Mean Clustering Coeff'); title('Normalized Small-Worldness');
    ylabel('C/C_n_u_l_l'); xlabel('L/L_n_u_l_l'); title('Small-worldness');
    set(gca,'FontSize',fntsz);
    xlim(pltlims); ylim(pltlims);
    RL=refline(1,0); set(RL,idlineopt{:});                                  % Identity line (works in some plots...)
    legend(AW8s(Isel))
end


% Violin
% pposv                   = [1442 445 449 352];
for jj = 1 : J
    Dx                  = CharPath_nrm(:,jj);
    Dy                  = ClstrCof_nrm(:,jj);
    Dplt                = zeros(size(emp_raw_tmp,1),I);
    pltstr              = ['Normalized Small-Worldness (Subject, ' Aparcs{jj} ')'];
    
    for ii = 1 : I
        Dx_tmp          = Dx{ii};
        Dy_tmp          = Dy{ii};
        Dplt(:,ii)      = Dy_tmp./Dx_tmp;                                   % Small-Worldness Computation
    end
    [~,V]=plot_violin(Dplt(:,Isel),AW8s(Isel),SmllWrld_nrm_g(Isel),'Small-Worldness',pltstr,Aparcs{jj},ppos,cmap(Isel,:)); hold on
    for vv=1:length(V); V(vv).ViolinAlpha=.1; end
    for vv=1:length(V); V(vv).ScatterPlot.MarkerFaceAlpha=.6; end
    line(xlim,[1 1],'LineStyle', '--', 'Color', 'k','LineWidth',2);
    
end




%--------------------------------------------------------------------------
end
