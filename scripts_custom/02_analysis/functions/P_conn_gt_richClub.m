function P_conn_gt_richClub(varargin)
% Rich Club Curves are plot for all group-level weighted & binary networks.
% 
% Requires network features computed in both empirical and null networks.
% (see conn_gt_richClub.m & conn_gt_all.m)
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
% Astr                        = A.info.details;
I                           = numel(AW8s);
J                           = numel(Aparcs);

% Which weights to plot
Isel                        = cell2mat(cellstrfind(AW8s,AW8s));             % All

% Plot settings
if isempty(cmap); cmap      = colormap(hsv(I)); end
fntsz                       = 25;
ppos                        = [-2559 625 2560 583];                          % monitor


%% --------------------------- GROUP-LEVEL ----------------------------- %%

%% Normalize Rich Club Curve

% Weighted Networks
metric                  = 'richclubcurve';
NULL_raw                = A.null.(metric);
EMPR_raw                = A.grp.(metric);
EMPR_nrm                = cell(size(EMPR_raw));

for jj = 1 : J
    for ii = 1 : I
        
        null_raw_tmp    = NULL_raw{ii,jj};
        null_raw_tmp    = nanmean(null_raw_tmp);                            % Measures normalized by mean across null networks
        
        emp_raw_tmp     = EMPR_raw{ii,jj};
        EMPR_nrm{ii,jj} = emp_raw_tmp./null_raw_tmp;                        % Normalize
    end
end
% RCC_nrm                 = cell2mat(EMPR_nrm);
RCC_nrm                 = EMPR_nrm;


% Binary
metric                  = 'richclubcurve';
NULL_raw_bin            = A.bin.null.(metric);
EMPR_raw_bin            = A.bin.grp.(metric);
EMPR_nrm_bin            = cell(size(EMPR_raw_bin));

for jj = 1 : J
        
    null_raw_tmp        = NULL_raw_bin{jj};
    null_raw_tmp        = nanmean(null_raw_tmp);                            % Measures normalized by mean across null networks

    emp_raw_tmp         = EMPR_raw_bin{jj};
    EMPR_nrm_bin{jj}    = emp_raw_tmp./null_raw_tmp;                        % Normalize
end
RCC_nrm_bin             = EMPR_nrm_bin;




%% Plot All Weights Together
lnwdth                  = 4;

disp('* ------------------ PLOTTING: GROUP Rich Club Curve -----------------*')

for jj = 1 : J
    kmax                = unique(cellfun(@(c) length(c), RCC_nrm(:,jj)));
    pltstr              = ['Normalized Rich Club Curve (group, ' Aparcs{jj} ')'];
    myfig(pltstr,ppos);
    
    % Plot weights
    for ii = Isel
        Dy              = RCC_nrm{ii,jj};
        plot(1:kmax, Dy, 'LineWidth', lnwdth, 'Color',cmap(ii,:)); hold on
    end
    
    % Add Binary
    Dy                  = RCC_nrm_bin{jj};
    plot(1:kmax, Dy, 'LineWidth', 2, 'Color',[.5 .5 .5],'LineStyle','--');
    
    xlabel('k level'); ylabel('Phi norm'); 
    legend([AW8s(Isel) 'Binary']); % legend(AW8s(Isel));
    title('Normalized Rich Club Curve: Group'); set(gca,'FontSize',fntsz);
    line(xlim,[1 1],'Color','k','LineWidth',2,'LineStyle',':');
end


%% Plot Specific Weights Together
lnwdth                  = 5;


% NoS, SIFT2 & COMMIT
pltw8s      = {'nos' 'sift2' 'commit'};
ipltw8s     = cell2mat(cellstrfind(AW8s, pltw8s));

for jj = 1 : J
    kmax                = unique(cellfun(@(c) length(c), RCC_nrm(ipltw8s,jj)));
    pltstr              = ['Normalized Rich Club Curve (group, ' Aparcs{jj} ')'];
    myfig(pltstr,[-2559 774 1314 434]);  % myfig(pltstr,ppos);
    
    
    for ii = ipltw8s
        Dy              = RCC_nrm{ii,jj};
        plot(1:kmax, Dy, 'LineWidth', lnwdth, 'Color',cmap(ii,:)); hold on
    end
    
    xlabel('k level'); ylabel('Phi norm'); legend(AW8s(ipltw8s));
    title('Normalized Rich Club Curve: Group'); set(gca,'FontSize',fntsz);
    line(xlim,[1 1],'Color','k','LineWidth',3,'LineStyle','--');
end


% R1, ICVF, RD, FA + binary
pltw8s      = {'r1' 'icvf' 'rd' 'fa'};
ipltw8s     = cell2mat(cellstrfind(AW8s, pltw8s));

for jj = 1 : J
    kmax                = unique(cellfun(@(c) length(c), RCC_nrm(ipltw8s,jj)));
    pltstr              = ['Normalized Rich Club Curve (group, ' Aparcs{jj} ')'];
    myfig(pltstr,[-2559 774 1314 434]);     % myfig(pltstr,ppos);
    
    % Weighted networks
    for ii = ipltw8s
        Dy              = RCC_nrm{ii,jj};
        plot(1:kmax, Dy, 'LineWidth', lnwdth, 'Color',cmap(ii,:)); hold on
    end
    
    % Add Binary
    Dy                  = RCC_nrm_bin{jj};
    plot(1:kmax, Dy, 'LineWidth', lnwdth, 'Color',[.5 .5 .5],'LineStyle',':');
    
    xlabel('k level'); ylabel('Phi norm'); legend([AW8s(ipltw8s) 'Binary']);
    title('Normalized Rich Club Curve: Group'); set(gca,'FontSize',fntsz);
    line(xlim,[1 1],'Color','k','LineWidth',3,'LineStyle','--');
end


%% Plot individual weights with p value

metric                  = 'richclubcurve';
NULL_raw                = A.null.(metric);
EMPR_raw                = A.grp.(metric);
m=3; n=3;

for jj = 1 : J
    kmax                = unique(cellfun(@(c) length(c), EMPR_raw(:,jj)));
    Nnull               = unique(cellfun(@(c) length(c), NULL_raw(:,jj)));
    pltstr              = ['Rich Club Curve (group, ' Aparcs{jj} ')'];
    myfig(pltstr,[-2559 134 2560 1074]);
    sp=1;
    for ii = Isel
        subplot(m,n,sp)
        pval            = zeros(1,kmax);
        Dgroup          = EMPR_raw{ii,jj};
        Dnull           = NULL_raw{ii,jj};
        
        for ik = 1 : kmax
            if ~isnan(Dgroup(ik))
                pval(ik) = sum(Dnull(:,ik)>Dgroup(ik)) / Nnull;             % Compute pvalue
            else 
                pval(ik) = nan;
            end
        end
        
        plot([Dgroup; mean(Dnull); pval]','LineWidth',3);
        xlabel('k level'); legend({'phi','phi rand','p'}); 
        title(AW8s(ii)); set(gca,'FontSize',25);
        sp=sp+1;
    end
end



%% Plot Binary RCC

ppos            = [-2558 472 2559 736];

GRPrcc          = cell2mat(A.bin.grp.richclubcurve);
GRPnds          = cell2mat(A.bin.grp.nodecount);
GRPegs          = cell2mat(A.bin.grp.edgecount);

NULrcc          = cell2mat(A.bin.null.richclubcurve);
NULnds          = cell2mat(A.bin.null.nodecount);
NULegs          = cell2mat(A.bin.null.edgecount);

% Compute mean across nulls
RC_null_mean    = mean(NULrcc);
ND_null_mean    = mean(NULnds);
EG_null_mean    = mean(NULegs);


% Normalize empirical RC
RC_norm         = GRPrcc ./ RC_null_mean;


% Plot Null & Empirical
n=2;m=2;
myfig('Null & Empirical Binary Group Structural Network',ppos);

subplot(n,m,1); 
plot(1:kmax,GRPrcc,'-b','LineWidth',3); hold on; plot(1:kmax,RC_null_mean,'-r','LineWidth',3); 
xlabel('k level'); ylabel('phi'); title('Rich Club Coeff'); set(gca,'FontSize',25); legend({'Empr' 'Null'});

subplot(n,m,3); 
plot(1:kmax,RC_norm,'-b','LineWidth',3);
xlabel('k level');  ylabel('phi norm'); title('Normalized RCC'); set(gca,'FontSize',25);

subplot(n,m,2); 
plot(1:kmax,GRPnds,'-b','LineWidth',3); hold on; plot(1:kmax,ND_null_mean,'-r','LineWidth',3);
xlabel('k level'); title('Remaining Nodes'); set(gca,'FontSize',25); legend({'Empr' 'Null'});

subplot(n,m,4); 
plot(1:kmax,GRPegs,'-b','LineWidth',3); hold on; plot(1:kmax,EG_null_mean,'-r','LineWidth',3);
xlabel('k level'); title('Remaining Edges'); set(gca,'FontSize',25); legend({'Empr' 'Null'});

    


%% --------------------------- SUBJECT-LEVEL -----------------------------%


% % % %% Normalize network measures
% % % 
% % % % Char Path
% % % metric                  = 'charpath';
% % % NULL_raw                = A.null.(metric);
% % % EMPR_raw                = A.sub.(metric);
% % % EMPR_nrm                = cell(size(EMPR_raw));
% % % 
% % % for jj = 1 : J
% % %     for ii = 1 : I
% % %         
% % %         null_raw_tmp    = NULL_raw{ii,jj};
% % %         null_raw_tmp    = mean(null_raw_tmp);                               % Measures normalized by mean across null networks
% % %         
% % %         emp_raw_tmp     = EMPR_raw{ii,jj};
% % %         EMPR_nrm{ii,jj} = emp_raw_tmp/null_raw_tmp;                         % Normalize
% % %     end
% % % end
% % % RCC_nrm            = EMPR_nrm;
% % % 
% % % 
% % % % Clustering Coeff
% % % metric                  = 'clusteringcoeff';
% % % NULL_raw                = A.null.(metric);
% % % EMPR_raw                = A.sub.(metric);
% % % EMPR_nrm                = cell(size(EMPR_raw));
% % % 
% % % for jj = 1 : J
% % %     for ii = 1 : I
% % %         
% % %         null_raw_tmp    = NULL_raw{ii,jj};
% % %         null_raw_tmp    = mean(null_raw_tmp,1);                             % Measures normalized by mean across null networks WITHIN NODE
% % %         
% % %         emp_raw_tmp     = EMPR_raw{ii,jj};
% % %         subarray        = zeros(size(emp_raw_tmp));
% % %         
% % %         for ss = 1 : size(emp_raw_tmp,1)
% % %             subarray(ss,:)  = emp_raw_tmp(ss,:)./null_raw_tmp;              % Normalize each node separately
% % %         end
% % %             
% % %         EMPR_nrm{ii,jj} = mean(subarray,2);                                 % Mean over nodes
% % %     end
% % % end
% % % ClstrCof_nrm            = EMPR_nrm;
% % % 
% % % 
% % % %% Plot
% % % 
% % % disp('* ------------------ PLOTTING: SUBJECT Small-Worldness -----------------*')
% % % 
% % % % Scatter
% % % for jj = 1 : J
% % %     Dx                  = RCC_nrm(:,jj);
% % %     Dy                  = ClstrCof_nrm(:,jj);
% % %     pltlims             = [floor(min(cellfun(@(c) min(c), [Dx; Dy]))) ceil(max(cellfun(@(c) max(c), [Dx; Dy])))];
% % %     pltstr              = ['Normalized Small-Worldness (Subject, ' Aparcs{jj} ')'];
% % %     myfig(pltstr,ppos);
% % %     
% % %     for ii = Isel
% % %         plot(Dx{ii},Dy{ii}, 'k.', 'MarkerSize',25, 'MarkerEdgeColor',cmap(ii,:)); hold on
% % %     end
% % %     
% % %     xlabel('Characteristic Path Length'); ylabel('Mean Clustering Coeff'); 
% % %     title('Normalized Small-Worldness'); set(gca,'FontSize',fntsz);
% % %     xlim(pltlims); ylim(pltlims);
% % %     RL=refline(1,0); set(RL,idlineopt{:});                                  % Identity line (works in some plots...)
% % %     legend(AW8s(Isel))
% % % end
% % % 
% % % 
% % % % Violin
% % % % pposv                   = [1442 445 449 352];
% % % for jj = 1 : J
% % %     Dx                  = RCC_nrm(:,jj);
% % %     Dy                  = ClstrCof_nrm(:,jj);
% % %     Dplt                = zeros(size(emp_raw_tmp,1),I);
% % %     pltstr              = ['Normalized Small-Worldness (Subject, ' Aparcs{jj} ')'];
% % %     
% % %     for ii = 1 : I
% % %         Dx_tmp          = Dx{ii};
% % %         Dy_tmp          = Dy{ii};
% % %         Dplt(:,ii)      = Dy_tmp./Dx_tmp;                                   % Small-Worldness Computation
% % %     end
% % %     [~,V]=plot_violin(Dplt(:,Isel),AW8s(Isel),SmllWrld_nrm_g(Isel),'Small-Worldness',pltstr,Aparcs{jj},ppos,cmap(Isel,:)); hold on
% % %     for vv=1:length(V); V(vv).ViolinAlpha=.1; end
% % %     for vv=1:length(V); V(vv).ScatterPlot.MarkerFaceAlpha=.6; end
% % %     line(xlim,[1 1],'LineStyle', '--', 'Color', 'k','LineWidth',2);
% % %     
% % % end




%--------------------------------------------------------------------------
end
