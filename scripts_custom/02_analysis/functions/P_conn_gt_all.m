function P_conn_gt_all(D,W8s,pinfo,varargin)
%  Currently creates a bunch of different plots of network topological
%  measures derived from weighted connectivity data. 
%
% Input:
%                           + Required +
%   D       : IxJ cell array, each cell containing an NxNxS stack of
%             connectivity matrices (S=subs)                                @MCN Only really need to include FC & LoS
%   W8s     : Cell array of strings describing edge weights                 @MCN maybe not ideal to use this to describe input datasets & network attributes of interest???
%   pinfo   : 1xJ structure with parcellation info (see conn_getParcInfo.m)
%
%                           + Optional +
%   A       : Structure storing network attributes (see conn_gt_all.m)
%   path    : full file path to saved A structure
%   str     : character vector for plot title  
%   cmap    : Ix3 colormap
%
% 2022 Mark C Nelson, McConnell Brain Imaging Centre, MNI, McGill
%--------------------------------------------------------------------------

[I,J]               = size(D);                                      % data dimensions: weights x parcellations

%% Optional inputs
nin                 = max(nargin,1) - 3;
defaults            = {[],[],[],colormap(hsv(I))};
[A,path,str,cmap]   = INhandler(varargin,nin,defaults);

%% Setup
if isempty(A)
    if ~isempty(path)
        load(path)
    else
        error('* * Must supply either A or PATH input to plot_conn_networkAtts.m')
    end
end
% % Manual load
% load('~/Desktop/ConnectomeProj/myelinWeightedConnectome/0_SCFingerprinting/data/5_networkAtts_Uthr50_UniDensity.mat')

% Get info from A & D
parcs                       = {pinfo.name};                                 % parcellation names
CIs                         = {pinfo.cirois};                               % community assignments
fci                         = cellstrfind(W8s,'fc');
losi                        = cellstrfind(W8s,'los');
Wiint                       = setdiff(1:I,[fci losi]);                      % indices for weights of interest (excluding fc & los)
Nwint                       = numel(Wiint);
Nsubs                       = min(cellfun(@(c) size(c,3),D));              % assume 3rd dim is subjects & consistent across datasets
AW8s                        = A.info.weights;
Aparcs                      = A.info.parcs;
Astr                        = A.info.details;
AWiint                      = cellfun(@(c) cellstrfind(AW8s,c),W8s(Wiint)); % Subset of weights in A we want & the order we want them in
pltW8s                      = AW8s(AWiint);                                 % labels for those weights of interest
lslineopt                   = {'Color','r','LineWidth',4,'LineStyle','--'}; % Options for least squares line used in scatter plots
idlineopt                   = {'Color','k','LineStyle',':','LineWidth',2};% For identity line

% Ensure A & D match
disp('* ------------------ PLOTTING: Network Attributes -----------------*')
disp('*  Note: the user must ensure that the data supplied in D matches the data which was used to compute A ')
disp(' ')
disp('  These data correspond to: ')
disp(['      Connectomes      : ' str])
disp(['      Network measures : ' A.info.details])
disp(['      D Parcellations  : ' parcs(:)])
disp(['      A Parcellations  : ' parcs(:)])
disp(['      D Weights        : ' W8s(Wiint)])
disp(['      A Weights        : ' pltW8s])


%% Plot them
fntsz           = 25;
lnstyl          = {'-' '--' ':'};
% ppos            = [1442 341 594 456];                             % monitor
ppos            = [1442 445 449 352];                             % monitor only 3 weights
% ppos            = [1 424 445 398];                                  % laptop
% rot8            = 20;
cmapt           = cmap(Wiint,:);



% Don't forget to compute density!
pltstr                      = {'Subject & group' Astr};
metricstr                   = 'DENSITY';
Dsub                        = D(Wiint,:);                                   % subset of weights of interest
DNSg                        = cell(1,Nwint);
DNSs                        = cell(Nsubs,Nwint);

for jj = 1 : J
    Dsp                     = Dsub(:,jj);                                   % datasets by parcellation
    for ii = 1 : Nwint
        Dspw                = Dsp{ii,jj};                                   % by parcellation & weight
        Dgpw                = groupavg(Dspw,3,'nz'); 
        DNSg{ii}            = density_und(Dgpw);
        for ss = 1 : Nsubs
            DNSs{ss,ii}     = density_und(Dspw(:,:,ss));
        end
    end
    plot_violin(DNSs,pltW8s,DNSg,metricstr,pltstr,parcs{jj},ppos,cmapt)
end


keyboard

%% ------------------------ Global measures ---------------------- %%

% ---------------------- VIOLIN PLOTS ----------------------- %
% Modularity, Transitivity, Assortativity
datanames       = {'modularity' 'transitivity' 'assortativity'};
pltstr          = {'Subject & group' Astr};
for dns = 1:length(datanames)
    metric      = datanames{dns};
    metricstr   = upper(strrep(metric,'_','-'));
    Dgrp        = A.grp.(metric)(AWiint,:);
    Dsub        = A.sub.(metric)(AWiint,:);

    % plot
    for jj = 1 : J
        Qgt     = Dgrp(:,jj);
        Qst     = num2cell(cell2mat(Dsub(:,jj)'));                          % convert to sub x weight cell array
        plot_violin(Qst,pltW8s,Qgt,metricstr,pltstr,parcs{jj},ppos,cmapt)   % Violin plot
    end
end

% Characteristic Path Length
% datanames       = {'lambda' 'efficiency' 'radius' 'diameter'};
datanames       = {'lambda' 'efficiency' 'radius' 'diameter'};
pltstr          = {'Subject & group' Astr};
for dns = 1:length(datanames)
    metric      = datanames{dns};
    metricstr   = upper(strrep(metric,'_','-'));
    Dgrp        = A.grp.charpath.(metric)(AWiint,:);
    Dsub        = A.sub.charpath.(metric)(AWiint,:);

    % plot
    for jj = 1 : J
        Qgt     = Dgrp(:,jj);
        Qst     = num2cell(cell2mat(Dsub(:,jj)'));
        plot_violin(Qst,pltW8s,Qgt,metricstr,pltstr,parcs{jj},ppos,cmapt)
    end
end

% Mean strength
metric      = 'strength';
metricstr   = upper(strrep(metric,'_','-'));
Dgrp        = A.grp.(metric)(AWiint,:);
Dsub        = A.sub.(metric)(AWiint,:);

% plot
for jj = 1 : J
    Dgt             = Dgrp(:,jj);
    Dst             = Dsub(:,jj); 
    for ii = 1 : Nwint
        Dgt{ii,jj}  = mean(Dgt{ii,jj});                                     % compute mean over nodes
        Dst{ii,jj}  = mean(Dst{ii,jj},2);
    end
        
    Dst     = num2cell(cell2mat(Dst'));
    plot_violin(Dst,pltW8s,Dgt,['Mean ' metricstr],pltstr,parcs{jj},ppos,cmapt)
end


% Network density (currently not computed in conn_networkAtts.m)
% Note: Network Density is necessary later to normalize betweenness!
metricstr                   = 'DENSITY';
Dsub                        = D(Wiint,:);                                   % subset of weights of interest
DNSg                        = cell(1,Nwint);
DNSs                        = cell(Nsubs,Nwint);

for jj = 1 : J
    Dsp                     = Dsub(:,jj);                                   % datasets by parcellation
    for ii = 1 : Nwint
        Dspw                = Dsp{ii,jj};                                   % by parcellation & weight
        Dgpw                = groupavg(Dspw,3,'nz'); 
        DNSg{ii}            = density_und(Dgpw);
        for ss = 1 : Nsubs
            DNSs{ss,ii}     = density_und(Dspw(:,:,ss));
        end
    end
    plot_violin(DNSs,pltW8s,DNSg,metricstr,pltstr,parcs{jj},ppos,cmapt)
end

% Used to manually mod plots
% TMPX=[.4 .8]; TMPY=[.4 .8]; title('');ylabel('');xlabel('');set(gca,'XTick',TMPX,'XTickLabel',TMPX,'FontSize',40,'YTick',TMPY,'YTickLabel',TMPY)



%% Comparing subject raking across weights in Global metrics

% ---------------------- SCATTER PLOTS ----------------------- %
% First Global Topology Measures
datanames       = {'modularity' 'transitivity' 'assortativity'};
% ppost           = [1442 210 747 587];
Idx             = 2;
Idy             = 3;
jj              = 1;
pltstr          = ['Subjects ' Astr];

% plot all data
for dns = 1:length(datanames)
    
    metric          = datanames{dns};
    metricstr       = upper(strrep(metric,'_','-'));
    % Dgrp            = cell2mat(A.grp.(metric)(AWiint,jj)');
    Dsub            = cell2mat(A.sub.(metric)(AWiint,jj)');
    
    myfig(pltstr,ppos);
    plot(Dsub(:,Idx),Dsub(:,Idy),'ko','MarkerFaceColor','b','MarkerSize',15,'LineWidth',2); hold on
    ylabel(W8s(Wiint(Idy))); xlabel(W8s(Wiint(Idx))); title(metricstr); set(gca,'FontSize',fntsz);
    xlims=xlim; ylims=ylim;
    L=lsline; set(L,lslineopt{:});                                          % Least sqaures line
%     RL=refline(1,0); set(RL,idlineopt{:});                                  % Identity line (works in some plots...)
    xlim(xlims); ylim(ylims);
%     legend([L RL], {'Least Squares' 'Identity'})
    
    % Correlation
    Dt1 = Dsub(:,Idx);
    Dt2 = Dsub(:,Idy);
    Dt1(Dt1==0) = nan;
    Dt2(Dt2==0) = nan;
%     r       = corr(Dt1,Dt2,'Type','Pearson','Rows','complete');
    r       = corr(Dt1,Dt2,'Type','Spearman','Rows','complete');
    corrstr = sprintf('r = %.2f',r);
    disp(['Correlation for ' metricstr ' from ' W8s{Wiint(Idy)} ' & ' W8s{Wiint(Idx)} ': ' corrstr])
    title([metricstr ' (' corrstr ')']);
%     set(gcf,'Position',[473   428   470   369])
end

% Repeat with double y axes
datanames       = {'modularity' 'transitivity' 'assortativity'};
Idy             = [2 3];
Idx             = 1;
jj              = 1;
pltstr          = ['Subjects ' Astr];

% plot all data
for dns = 1:length(datanames)
    
    metric          = datanames{dns};
    metricstr       = upper(strrep(metric,'_','-'));
    % Dgrp            = cell2mat(A.grp.(metric)(AWiint,jj)');
    Dsub            = cell2mat(A.sub.(metric)(AWiint,jj)');
    
    % Correlation
    Dt1 = Dsub(:,Idx);
    Dt2 = Dsub(:,Idy(1));
    Dt3 = Dsub(:,Idy(2));
    Dt1(Dt1==0) = nan;
    Dt2(Dt2==0) = nan;
    Dt3(Dt3==0) = nan;
%     r1       = corr(Dt1,Dt2,'Type','Pearson','Rows','complete');
%     r2       = corr(Dt1,Dt3,'Type','Pearson','Rows','complete');
    r1       = corr(Dt1,Dt2,'Type','Spearman','Rows','complete');
    r2       = corr(Dt1,Dt3,'Type','Spearman','Rows','complete');
    corrstr1 = sprintf('r = %.2f',r1);
    corrstr2 = sprintf('r = %.2f',r2);

    % Plot
    myfig(pltstr,ppos);
    yyaxis left
    plot(Dsub(:,Idx),Dsub(:,Idy(1)),'ko','MarkerEdgeColor','b','MarkerSize',5,'LineWidth',2); hold on
    ylabel(W8s(Wiint(Idy(1))));  xlims=xlim; ylims=ylim;
    LL=lsline; set(LL,'Color','b','LineWidth',4,'LineStyle','--');
    xlim(xlims); ylim(ylims);

    yyaxis right
    plot(Dsub(:,Idx),Dsub(:,Idy(2)),'kd','MarkerEdgeColor','r','MarkerSize',5,'LineWidth',2); hold on
    ylabel(W8s(Wiint(Idy(2))));  xlims=xlim; ylims=ylim;
    LR=lsline; set(LR,'Color','r','LineWidth',4,'LineStyle','--');  
    xlim(xlims); ylim(ylims);
    
    % Set common L & R yaxes limits
    yyaxis left;  tmp = ylim;
    yyaxis right; tmp = [tmp; ylim];
    yylims = [min(tmp(:,1)) max(tmp(:,2))];
    yyaxis left; ylim(yylims); yyaxis right; ylim(yylims);
    
    xlabel(W8s(Wiint(Idx))); title(metricstr); set(gca,'FontSize',fntsz);
    legend([LL LR], {corrstr1 corrstr2})
%     set(gcf,'Position',[473   428   470   369])
end


% Char Path Metrics
datanames       = {'lambda' 'efficiency'};
% ppost           = [1442 210 747 587];
Idx             = 2;
Idy             = 3;
jj              = 1;
pltstr          = ['Subjects ' Astr];

% plot all data
for dns = 1:length(datanames)
    
    metric          = datanames{dns};
    metricstr       = upper(strrep(metric,'_','-'));
    % Dgrp            = cell2mat(A.grp.(metric)(AWiint,jj)');
    Dsub            = cell2mat(A.sub.charpath.(metric)(AWiint,jj)');

    myfig(pltstr,ppos);
    plot(Dsub(:,Idx),Dsub(:,Idy),'ko','MarkerFaceColor','b','MarkerSize',15,'LineWidth',2); hold on
    ylabel(W8s(Wiint(Idy))); xlabel(W8s(Wiint(Idx))); title(metricstr); set(gca,'FontSize',fntsz);
    xlims=xlim; ylims=ylim;
    L=lsline; set(L,lslineopt{:});
%     RL=refline(1,0); set(RL,idlineopt{:});
    xlim(xlims); ylim(ylims);
    
    % Correlation
    Dt1 = Dsub(:,Idx);
    Dt2 = Dsub(:,Idy);
    Dt1(Dt1==0) = nan;
    Dt2(Dt2==0) = nan;
%     r       = corr(Dt1,Dt2,'Type','Pearson','Rows','complete');
    r       = corr(Dt1,Dt2,'Type','S','Rows','complete');
    corrstr = sprintf('r = %.2f',r);
    disp(['Correlation for ' metricstr ' from ' W8s{Wiint(Idy)} ' & ' W8s{Wiint(Idx)} ': ' corrstr])
    title([metricstr ' (' corrstr ')']);
end


% Repeat with double y axes
datanames       = {'lambda' 'efficiency'};
Idy             = [2 3];
Idx             = 1;
jj              = 1;
pltstr          = ['Subjects ' Astr];

% plot all data
for dns = 1:length(datanames)
    
    metric          = datanames{dns};
    metricstr       = upper(strrep(metric,'_','-'));
    % Dgrp            = cell2mat(A.grp.(metric)(AWiint,jj)');
    Dsub            = cell2mat(A.sub.charpath.(metric)(AWiint,jj)');
    
    % Correlation
    Dt1 = Dsub(:,Idx);
    Dt2 = Dsub(:,Idy(1));
    Dt3 = Dsub(:,Idy(2));
    Dt1(Dt1==0) = nan;
    Dt2(Dt2==0) = nan;
    Dt3(Dt3==0) = nan;
%     r1       = corr(Dt1,Dt2,'Type','Pearson','Rows','complete');
%     r2       = corr(Dt1,Dt3,'Type','Pearson','Rows','complete');
    r1       = corr(Dt1,Dt2,'Type','Spearman','Rows','complete');
    r2       = corr(Dt1,Dt3,'Type','Spearman','Rows','complete');

    corrstr1 = sprintf('r = %.2f',r1);
    corrstr2 = sprintf('r = %.2f',r2);

    % Plot
    myfig(pltstr,ppos);
    yyaxis left
    plot(Dsub(:,Idx),Dsub(:,Idy(1)),'ko','MarkerEdgeColor','b','MarkerSize',5,'LineWidth',2); hold on
    ylabel(W8s(Wiint(Idy(1))));  xlims=xlim; ylims=ylim;
    LL=lsline; set(LL,'Color','b','LineWidth',4,'LineStyle','--');
    xlim(xlims); ylim(ylims);

    yyaxis right
    plot(Dsub(:,Idx),Dsub(:,Idy(2)),'kd','MarkerEdgeColor','r','MarkerSize',5,'LineWidth',2); hold on
    ylabel(W8s(Wiint(Idy(2))));  xlims=xlim; ylims=ylim;
    LR=lsline; set(LR,'Color','r','LineWidth',4,'LineStyle','--');  
    xlim(xlims); ylim(ylims);
    
    % Set common L & R yaxes limits
    yyaxis left;  tmp = ylim;
    yyaxis right; tmp = [tmp; ylim];
    yylims = [min(tmp(:,1)) max(tmp(:,2))];
    yyaxis left; ylim(yylims); yyaxis right; ylim(yylims);

    
    xlabel(W8s(Wiint(Idx))); title(metricstr); set(gca,'FontSize',fntsz);
    legend([LL LR], {corrstr1 corrstr2})
%     set(gcf,'Position',[473   428   470   369])
end



% Mean Strength
datanames       = {'strength'};
% ppost           = [1442 210 747 587];
% Idx             = 2;
% Idy             = 3;
% jj              = 1;
% pltstr          = ['Subjects ' Astr];

metric          = datanames{1};
metricstr       = upper(strrep(metric,'_','-'));
Dsub            = A.sub.(metric)(AWiint,jj);
Dsub            = cell2mat(cellfun(@(c) mean(c,2),Dsub,'UniformOutput',0)');

myfig(pltstr,ppos);
plot(Dsub(:,Idx),Dsub(:,Idy),'ko','MarkerFaceColor','b','MarkerSize',15,'LineWidth',2); hold on
ylabel(W8s(Wiint(Idy))); xlabel(W8s(Wiint(Idx))); title(metricstr); set(gca,'FontSize',fntsz);
xlims=xlim; ylims=ylim;
L=lsline; set(L,lslineopt{:});
RL=refline(1,0); set(RL,idlineopt{:});
xlim(xlims); ylim(ylims);

% Correlation
Dt1 = Dsub(:,Idx);
Dt2 = Dsub(:,Idy);
Dt1(Dt1==0) = nan;
Dt2(Dt2==0) = nan;
r       = corr(Dt1,Dt2,'Type','Pearson','Rows','complete');
corrstr = sprintf('r = %.2f',r);
title([metricstr ' (' corrstr ')']);
disp(['Correlation for ' metricstr ' from ' W8s{Wiint(Idy)} ' & ' W8s{Wiint(Idx)} ': ' corrstr])





% % % % Further manual exploration of Modularity
% % % % With subset of subjects marked
% % % [~,tmpIX]   = sort(Dsub(:,Idx)','descend');
% % % [~,tmpIY]   = sort(Dsub(:,Idy)','descend');
% % % tmpIX 
% % % tmpIY
% % % find(tmpIX==tmpIY(1))
% % % 
% % % Isub_loY_hiX    = 18;                                                       % Manually defined subset of subjects to focus on
% % % Isub_hiY_loX    = 9;
% % % Isub_loYonly    = 14;
% % % Isub_hiYonly    = 13;
% % % Isub_loboth     = [16 26];
% % % Isub_hiboth     = [44 40];
% % % 
% % % Xall            = Dsub(:,Idx);                  
% % % Yall            = Dsub(:,Idy);
% % % 
% % % X_loY_hiX       = Xall(Isub_loY_hiX);  
% % % Y_loY_hiX       = Yall(Isub_loY_hiX);
% % % 
% % % X_hiY_loX       = Xall(Isub_hiY_loX);
% % % Y_hiY_loX       = Yall(Isub_hiY_loX);
% % % 
% % % X_loboth        = Xall(Isub_loboth);  
% % % Y_loboth        = Yall(Isub_loboth);  
% % % 
% % % X_hiboth        = Xall(Isub_hiboth);
% % % Y_hiboth        = Yall(Isub_hiboth);
% % % 
% % % X_hiYonly        = Xall(Isub_hiYonly);
% % % Y_hiYonly        = Yall(Isub_hiYonly);
% % % 
% % % X_loYonly        = Xall(Isub_loYonly);
% % % Y_loYonly        = Yall(Isub_loYonly);
% % % 
% % % 
% % % myfig(pltstr,ppost);
% % % plot(Xall,Yall,'ko','MarkerSize',10,'LineWidth',2); hold on
% % % xlims=xlim; ylims=ylim;
% % % L=lsline; set(L, 'Color', 'k','LineWidth',2,'LineStyle','--');
% % % xlim(xlims); ylim(ylims);
% % % 
% % % XL=plot(X_loY_hiX,Y_loY_hiX,'ks','MarkerFaceColor',cmapt(Wiint(Idx),:),'MarkerSize',15,'LineWidth',2); hold on
% % % YL=plot(X_hiY_loX,Y_hiY_loX,'kh','MarkerFaceColor',cmapt(Wiint(Idy),:),'MarkerSize',15,'LineWidth',2); hold on
% % % BL=plot(X_loboth,Y_loboth,'kv','MarkerFaceColor',[0,0,1],'MarkerSize',15,'LineWidth',2); hold on
% % % BH=plot(X_hiboth,Y_hiboth,'k^','MarkerFaceColor',[0,1,0],'MarkerSize',15,'LineWidth',2); hold on
% % % YHo=plot(X_hiYonly,Y_hiYonly,'k^','MarkerFaceColor',cmapt(Wiint(Idy),:),'MarkerSize',15,'LineWidth',2); hold on
% % % YLo=plot(X_loYonly,Y_loYonly,'kv','MarkerFaceColor',cmapt(Wiint(Idy),:),'MarkerSize',15,'LineWidth',2); hold on
% % % 
% % % 
% % % ylabel(W8s(Wiint(Idy))); xlabel(W8s(Wiint(Idx))); title(metricstr); set(gca,'FontSize',fntsz);
% % % xlim(xlims); ylim(ylims);
% % % strlbl1 = strcat('Hi-',W8s(Wiint([Idx Idy])),'-Lo-',  W8s(Wiint([Idy Idx])));
% % % strlbl2 = ['Hi-' W8s{Wiint(Idy)} ' only'];
% % % strlbl3 = ['Lo-' W8s{Wiint(Idy)} ' only'];
% % % legend([XL YL BL BH YHo YLo], [strlbl1 'Lo-both' 'Hi-both' strlbl2 strlbl3])





%% ------------------ Local measures at Group Level ---------------------------- %%

% ---------------------- VIOLIN PLOTS ----------------------- %
% Strength, Clustering Coefficient, NodeBetweenness, Participation Coefficient, Local Efficiency
datanames   = {'strength' 'clusteringcoeff' 'efficiency_local' 'nodebetweenness' 'participationcoeff' };
Nnodes      = length(pinfo.labels);                                         % node count (havent extended this section to handle multiple parcs)
pltstr      = {'Group-level node values' Astr};

for jj = 1 : J
    for dns = 1:length(datanames)
        metric      = datanames{dns};
        metricstr   = upper(strrep(metric,'_','-'));
        Dgrp        = A.grp.(metric)(AWiint,jj);
        Dall        = cell(Nwint,1);                                        % store all group level data for violin plot

        % plot cdf
        myfig([pltstr{2} ': ' parcs{jj}], ppos)
        lnstyln             = 1;                                            % used to iterate through line styles
        for ii = 1 : Nwint
            dg              = Dgrp{ii,1};
            if strncmpi(metric,'nodebe',6)                                  % if metric is node betweenness
                DNSgt       = DNSg{ii};                                     % NEED TO RUN DENSITY SECTION FIRST (see lines 130-150ish)
                E           = find(dg);
                dg(E)       = dg(E)./((Nnodes-1)*(Nnodes-2)*DNSgt*.01);     % Rescale to percentage of total possible paths given network density
            end
            Dall{ii}        = dg;                                           % store for violin
            
            % plot cdf
            h               = cdfplot(dg); hold on
            h.LineWidth     = 3;
            h.Color         = cmap(ii,:);
            h.LineStyle     = lnstyl{lnstyln};
            if lnstyln < 3
                lnstyln     = lnstyln + 1;
            else
                lnstyln = 1;
            end
        end
        title(pltstr{1}); legend(pltW8s); xlabel(metricstr); ylabel('cdf');
        set(gca,'FontSize', fntsz);

        % Plot violin
        plot_violin(Dall',pltW8s,[],metricstr,pltstr,parcs{jj},ppos,cmapt)
    end


    % Eccentricity
    datanames   = {'eccentricity'};

    for dns = 1:length(datanames)
        metric      = datanames{dns};
        metricstr   = upper(strrep(metric,'_','-'));
        Dgrp        = A.grp.charpath.(metric)(AWiint,jj);
        Dall        = cell(Nwint,1);                                        % store all group level data for violin plot

        % plot cdf
        myfig([pltstr{2} ': ' parcs{jj}], ppos)
        lnstyln             = 1;                                            % used to iterate through line styles
        for ii = 1 : Nwint
            dg              = Dgrp{ii,1};
            Dall{ii}        = dg;                                           % store for violin
            h               = cdfplot(dg); hold on
            h.LineWidth     = 3;
            h.Color         = cmap(ii,:);
            h.LineStyle     = lnstyl{lnstyln};
            if lnstyln < 3
                lnstyln     = lnstyln + 1;
            else
                lnstyln = 1;
            end
        end
        title(pltstr{1}); legend(pltW8s); xlabel(metricstr); ylabel('cdf');
        set(gca,'FontSize', fntsz);

        % Plot violin
        plot_violin(Dall',pltW8s,[],metricstr,pltstr,parcs{jj},ppos,cmapt)
    end
end


%% ------------------ Compare Node raking across weights --------------- %%

% ---------------------- SCATTER PLOTS ----------------------- %
% Group Level
datanames       = {'strength' 'clusteringcoeff' 'efficiency_local' 'nodebetweenness' 'participationcoeff' };
Ndatanames      = length(datanames);
% ppost           = [1442 210 747 587];
Idx             = 1;
Idy             = 3;
jj              = 1;
pltstr          = ['Nodes ' Astr];

for dns = 1:Ndatanames
    metric          = datanames{dns};
    metricstr       = upper(strrep(metric,'_','-'));
    Xgrp            = A.grp.(metric){AWiint(Idx),jj};
    Ygrp            = A.grp.(metric){AWiint(Idy),jj};
    
    if strncmpi(metric,'nodebetweenness',6)
        Irm         = find(Xgrp==0 | Ygrp==0);
        Xgrp(Irm)   = [];
        Ygrp(Irm)   = [];
        Xgrp        = log(Xgrp);                
        Ygrp        = log(Ygrp);
        metricstr   = [metricstr ' (log)'];
    end
    
    % plot all data
    myfig(pltstr,ppos);
    plot(Xgrp,Ygrp,'ko','MarkerFaceColor','b','MarkerSize',15,'LineWidth',2); hold on
    ylabel(W8s{Wiint(Idy)}); xlabel(W8s{Wiint(Idx)}); title(['Group Node ' metricstr]); set(gca,'FontSize',fntsz);
    xlims=xlim; ylims=ylim;
    L=lsline; set(L,lslineopt{:});
%     RL=refline(1,0); set(RL,idlineopt{:});
    xlim(xlims); ylim(ylims);
    
    % Correlation
    Xgrp(Xgrp==0)   = nan;
    Ygrp(Ygrp==0)   = nan;
    r               = corr(Xgrp',Ygrp','Type','Pearson','Rows','complete');
    corrstr = sprintf('r=%.2f',r);
    title([metricstr ' (' corrstr ')']);
    disp(['Correlation for ' metricstr ' from ' W8s{Wiint(Idy)} ' & ' W8s{Wiint(Idx)} ': ' corrstr])
    
    % Alternate plot option: HEAT SCATTER
    myfig(pltstr,ppos); 
    heatscatter(Xgrp',Ygrp');
    ylabel(W8s{Wiint(Idy)}); xlabel(W8s{Wiint(Idx)}); title(['Group Node ' metricstr]); set(gca,'FontSize',fntsz);
    xlims=xlim; ylims=ylim;
    text(xlims(1)+1,ylims(2)-1,corrstr,'FontSize',30);
end


% Repeat with double y axes
datanames       = {'strength' 'clusteringcoeff'}; % 'efficiency_local' 'nodebetweenness' 'participationcoeff' };
Idy             = [2 3];
Idx             = 1;
jj              = 1;
pltstr          = ['Nodes ' Astr];

for dns = 1:length(datanames)
    
    metric          = datanames{dns};
    metricstr       = upper(strrep(metric,'_','-'));
    Dgrp            = cell2mat(A.grp.(metric)(AWiint,jj));
%     Dsub            = cell2mat(A.sub.(metric)(AWiint,jj)');
    
    % Correlation
    Dt1 = Dgrp(Idx,   :)';
    Dt2 = Dgrp(Idy(1),:)';
    Dt3 = Dgrp(Idy(2),:)';
    Dt1(Dt1==0) = nan;
    Dt2(Dt2==0) = nan;
    Dt3(Dt3==0) = nan;
    r1       = corr(Dt1,Dt2,'Type','Pearson','Rows','complete');
    r2       = corr(Dt1,Dt3,'Type','Pearson','Rows','complete');
    corrstr1 = sprintf('r = %.2f',r1);
    corrstr2 = sprintf('r = %.2f',r2);

    % Plot
    myfig(pltstr,ppos);
    yyaxis left
    plot(Dgrp(Idx,:),Dgrp(Idy(1),:),'k+','MarkerEdgeColor','b','MarkerSize',6,'LineWidth',1); hold on
    ylabel(W8s(Wiint(Idy(1))));  xlims=xlim; ylims=ylim;
    LL=lsline; set(LL,'Color','b','LineWidth',4,'LineStyle','--');
    xlim(xlims); ylim(ylims);

    yyaxis right
    plot(Dgrp(Idx,:),Dgrp(Idy(2),:),'k.','MarkerEdgeColor','r','MarkerSize',6,'LineWidth',1); hold on
    ylabel(W8s(Wiint(Idy(2))));  xlims=xlim; ylims=ylim;
    LR=lsline; set(LR,'Color','r','LineWidth',4,'LineStyle','--');  
    xlim(xlims); ylim(ylims);
    
    xlabel(W8s(Wiint(Idx))); title(metricstr); set(gca,'FontSize',fntsz);
    legend([LL LR], {corrstr1 corrstr2})
%     set(gcf,'Position',[473   428   470   369])
end



datanames       = {'eccentricity'};
Ndatanames      = length(datanames);

for dns = 1:Ndatanames
    metric          = datanames{dns};
    metricstr       = upper(strrep(metric,'_','-'));
    Xgrp            = A.grp.charpath.(metric){AWiint(Idx),jj};
    Ygrp            = A.grp.charpath.(metric){AWiint(Idy),jj};
    
    % plot all data
    myfig(pltstr,ppos);
    plot(Xgrp,Ygrp,'ko','MarkerFaceColor','b','MarkerSize',15,'LineWidth',2); hold on
    ylabel(W8s{Wiint(Idy)}); xlabel(W8s{Wiint(Idx)}); title(['Group Node ' metricstr]); set(gca,'FontSize',fntsz);
    xlims=xlim; ylims=ylim;
    L=lsline; set(L,lslineopt{:});
    RL=refline(1,0); set(RL,idlineopt{:});
    xlim(xlims); ylim(ylims);
    
    % Correlation
    Xgrp(Xgrp==0)   = nan;
    Ygrp(Ygrp==0)   = nan;
    r               = corr(Xgrp',Ygrp','Type','Pearson','Rows','complete');
    corrstr         = sprintf('r = %.2f',r);
    title([metricstr ' (' corrstr ')']);
    disp(['Correlation for ' metricstr ' from ' W8s{Wiint(Idy)} ' & ' W8s{Wiint(Idx)} ': ' corrstr])
end



%% ------------------ Local measures at Subject Level ---------------------------- %%

% ---------------------- VIOLIN PLOTS ----------------------- %
% Strength, Clustering Coefficient, NodeBetweenness, Participation Coefficient, Local Efficiency
datanames   = {'strength' 'clusteringcoeff' 'efficiency_local' 'nodebetweenness' 'participationcoeff' };
Nnodes      = length(pinfo.labels);                                         % node count (havent extended this section to handle multiple parcs)
pltstr      = {'Subject-level node values' Astr};

for jj = 1 : J
    for dns = 1:length(datanames)
        metric      = datanames{dns};
        metricstr   = upper(strrep(metric,'_','-'));
        Dsub        = A.sub.(metric)(AWiint,jj);
        Dall        = cell(Nwint,1);                                        % store all data for violin plot

        % plot cdf
        myfig([pltstr{2} ': ' parcs{jj}], ppos)
        lnstyln             = 1;                                            % used to iterate through line styles
        for ii = 1 : Nwint
            ds              = Dsub{ii,1}(:)';
            if strncmpi(metric,'nodebe',6)                                  % if metric is node betweenness
                ds          = Dsub{ii,1};                                   % Re-extract data because need nodes organized by subject
                DNSst       = cell2mat(DNSs(:,ii));                         % NEED TO RUN DENSITY SECTION FIRST (see lines 130-150ish)
                for kk = 1 : size(ds,1)
                    dst     = ds(kk,:);
                    E       = find(dst);
                    dst(E)  = dst(E)./((Nnodes-1)*(Nnodes-2)*DNSst(kk)*.01);% Rescale to percentage of total possible paths given network density
                    ds(kk,:)= dst;
                end
                ds          = ds(:)';
            end
            Dall{ii}        = ds;                                           % store for violin
            
            % plot cdf
            h               = cdfplot(ds); hold on
            h.LineWidth     = 3;
            h.Color         = cmap(ii,:);
            h.LineStyle     = lnstyl{lnstyln};
            if lnstyln < 3
                lnstyln     = lnstyln + 1;
            else
                lnstyln = 1;
            end
        end
        title(pltstr{1}); legend(pltW8s); xlabel(metricstr); ylabel('cdf');
        set(gca,'FontSize', fntsz);

        % Plot violin
        plot_violin(Dall',pltW8s,[],metricstr,pltstr,parcs{jj},ppos,cmapt)
    end


    % Eccentricity
    datanames   = {'eccentricity'};

    for dns = 1:length(datanames)
        metric      = datanames{dns};
        metricstr   = upper(strrep(metric,'_','-'));
        Dsub        = A.sub.charpath.(metric)(AWiint,jj);
        Dall        = cell(Nwint,1);                                        % store all group level data for violin plot

        % plot cdf
        myfig([pltstr{2} ': ' parcs{jj}], ppos)
        lnstyln             = 1;                                            % used to iterate through line styles
        for ii = 1 : Nwint
            ds              = Dsub{ii,1}(:)';
            Dall{ii}        = ds;                                           % store for violin
            h               = cdfplot(ds); hold on
            h.LineWidth     = 3;
            h.Color         = cmap(ii,:);
            h.LineStyle     = lnstyl{lnstyln};
            if lnstyln < 3
                lnstyln     = lnstyln + 1;
            else
                lnstyln = 1;
            end
        end
        title(pltstr{1}); legend(pltW8s); xlabel(metricstr); ylabel('cdf');
        set(gca,'FontSize', fntsz);

        % Plot violin
        plot_violin(Dall',pltW8s,[],metricstr,pltstr,parcs{jj},ppos,cmapt)
    end
end


%% ------------------ Compare Node ranking across weights --------------- %%

% ---------------------- SCATTER PLOTS ----------------------- %
% Subject Level
datanames       = {'strength' 'clusteringcoeff' 'efficiency_local' 'nodebetweenness' 'participationcoeff' };
Ndatanames      = length(datanames);
% ppost           = [1442 210 747 587];
Idx             = 1;
Idy             = 2;
jj              = 1;
pltstr          = ['Subject Nodes ' Astr];

for dns = 1:Ndatanames
    metric          = datanames{dns};
    metricstr       = upper(strrep(metric,'_','-'));
    Xsub            = A.sub.(metric){AWiint(Idx),jj}(:);
    Ysub            = A.sub.(metric){AWiint(Idy),jj}(:);
    
    if strncmpi(metric,'nodebetweenness',6)
        Irm         = find(Xsub==0 | Ysub==0);
        Xsub(Irm)   = [];
        Ysub(Irm)   = [];
        Xsub        = log(Xsub);                
        Ysub        = log(Ysub);
        metricstr   = [metricstr ' (log)'];
    end
    
    % Correlation
    Xsub(Xsub==0)   = nan;
    Ysub(Ysub==0)   = nan;
    r               = corr(Xsub,Ysub,'Type','Pearson','Rows','complete');
    corrstr         = sprintf('r = %.2f',r);
    
    % plot all data
    myfig(pltstr,ppos);
    plot(Xsub,Ysub,'ko','MarkerEdgeColor','k','MarkerFaceColor','b','MarkerSize',4,'LineWidth',.5); hold on
    ylabel(W8s{Wiint(Idy)}); xlabel(W8s{Wiint(Idx)}); title(['Subject Node ' metricstr ' (' corrstr ')']); set(gca,'FontSize',fntsz);
    xlims=xlim; ylims=ylim;
    L=lsline; set(L,lslineopt{:});
%     RL=refline(1,0); set(RL,idlineopt{:});
    xlim(xlims); ylim(ylims);
        
    % Alternate plot option: HEAT SCATTER
    myfig(pltstr,ppos); 
    heatscatter(Xsub,Ysub);
    ylabel(W8s{Wiint(Idy)}); xlabel(W8s{Wiint(Idx)}); title(['Subject Node ' metricstr]); set(gca,'FontSize',fntsz);
    xlims=xlim; ylims=ylim;
    text(xlims(1)+1,ylims(2)-1,corrstr,'FontSize',30);
%     set(gcf,'Position',[473   428   470   369])
end


datanames       = {'eccentricity'};
Ndatanames      = length(datanames);

for dns = 1:Ndatanames
    metric          = datanames{dns};
    metricstr       = upper(strrep(metric,'_','-'));
    Xsub            = A.sub.charpath.(metric){AWiint(Idx),jj}(:);
    Ysub            = A.sub.charpath.(metric){AWiint(Idy),jj}(:);
    
    % plot all data
    myfig(pltstr,ppos);
    plot(Xsub,Ysub,'ko','MarkerEdgeColor','k','MarkerFaceColor','b','MarkerSize',4,'LineWidth',.5); hold on
    ylabel(W8s{Wiint(Idy)}); xlabel(W8s{Wiint(Idx)}); title(['Subject Node ' metricstr]); set(gca,'FontSize',fntsz);
    xlims=xlim; ylims=ylim;
    L=lsline; set(L,lslineopt{:});
    RL=refline(1,0); set(RL,idlineopt{:});
    xlim(xlims); ylim(ylims);
    
    % Correlation
    Xsub(Xsub==0)   = nan;
    Ysub(Ysub==0)   = nan;
    r               = corr(Xsub,Ysub,'Type','Pearson','Rows','complete');
    corrstr         = sprintf('r = %.2f',r);
    title([metricstr ' (' corrstr ')']);
    disp(['Correlation for ' metricstr ' from ' W8s{Wiint(Idy)} ' & ' W8s{Wiint(Idx)} ': ' corrstr])
end


%% Now plot these top nodes on the surface WITH SUBCORTEX
% {'strength' 'clusteringcoeff' 'efficiency_local' 'nodebetweenness' 'participationcoeff' }

% Setup
jj                      = 1;
isel                    = [1 3];
Nnodecor                = [pinfo.ncor];
Nnodesub                = [pinfo.nsub];
Nnodeall                = Nnodecor + Nnodesub;
metric                  = 'participationcoeff';
pltmetric               = upper(strrep(metric,'_','-'));

% Rank data
drank                   = nan(numel(isel),Nnodeall(jj));
it                      = 0;                                                                % just an iterator
for ii = isel
    it                  = it + 1;
    d                   = A.grp.(metric){AWiint(ii),jj}';
    d                   = -1 * d;                                             % to get highest value as rank 1
    [drank(it,:), ~]    = tiedrank(d);                                        % rank data
end

% Convert rankings to consistency measure
cutoff                  = 50;                                               % Number of nodes to plot
drank(drank>cutoff)     = 0;
drank(1,drank(1,:)~=0)  = 1;                                                % binarize top nodes
drank(2,drank(2,:)~=0)  = 2;                                                % binarize top nodes
drank                   = sum(drank,1);                                     % sum over weights to get consistency ranking

% Plot surface data
titlestr    = [pltmetric ' top ' num2str(cutoff) ' ranked nodes: ' str ', ' parcs{jj}];
[Hcx,Hsx]   = plot_conn_surf(drank',pinfo(jj),'both',{pltmetric}, titlestr);

% Manual mods to colormaps
h=Hcx{1};
h.cb.Ticks  = [0 1 2 3];
modmap      = zeros(4,3);           % custom colormap
modmap(1,:) = [0.2081 0.1663 0.5292]; % dark blue
modmap(2,:) = [0 1 1];              % cyan
modmap(3,:) = [1 1 0];              % yellow
modmap(4,:) = [1 0 0];              % red
for ax =1:4
    colormap(h.axes(ax),modmap)
end

H=Hsx{1}; 
CB=H.cb;
CB.FontSize=fntsz;
CB.Ticks= [0 1 2 3];
for AXIS = 1 : 4
    colormap(H.axes(AXIS),modmap)
end


%% Absolute values of measures & Variance on surface
% {'strength' 'clusteringcoeff' 'efficiency_local' 'nodebetweenness' 'participationcoeff' }

jj                      = 1;
isel                    = AWiint;
Nnodecor                = [pinfo.ncor];
Nnodesub                = [pinfo.nsub];
Nnodeall                = Nnodecor + Nnodesub;
metric                  = 'clusteringcoeff';
pltmetric               = upper(strrep(metric,'_','-'));

% Extract data
dgrp                    = cell2mat(A.grp.(metric)(AWiint,jj));

% Plot surface data
titlestr    = [pltmetric ': ' str ', ' parcs{jj}];
[Hcx,Hsx]   = plot_conn_surf(dgrp',pinfo(jj),'both',AW8s(AWiint), titlestr);

% Standardize colormaps & colorbars across cortical & subcortical plots
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
AXall           = unicb(AXall);                                         % Standardize all colormaps
clims           = AXall(1).CLim;

for cc = 1 : length(CBall)
    CBall(cc).Ticks     = clims;
    CBall(cc).FontSize  = fntsz;
end 

% Variability in subject level data
dsub = A.sub.(metric)(AWiint,jj);
dsub = cellfun(@(c) c', dsub, 'UniformOutput', 0);
plot_conn_var(dsub,pltW8s,'cqd',pinfo,[upper(metric) '; ' str],cmapt);
% set(gcf,'Position',ppos)


%% ----------------- Matrices ------------------ %%         *** @MCN Could do separate modules exploring these further?
% ppos            = [1442 1 1010 796];                                        % Small number of plots
ppos            = [1441 362 520 435];                                       % large number of plots
Nnodes          = length(pinfo.labels);                                     % node count (havent extended this section to handle multiple parcs)
axs             = cell(Nwint,1);                                            % Temporary storage for axes

% Plot EdgeBetweenness (group)                      @MCN Separate module exploring this further?
metric          = 'edgebetweenness';
for ii = 1 : Nwint
    itA         = AWiint(ii);
    strt        = ['Edge Betweeness (group): ' AW8s{itA}];
    DNSgt       = DNSg{ii};                                                 % density of this network at group level
    dt          = A.grp.(metric){itA};
    E           = find(dt);
    dt(E)       = dt(E)./((Nnodes-1)*(Nnodes-2)*DNSgt*.01);                 % Rescale to percentage of total possible paths relative to density
    connplot(dt, strt, pinfo, ppos)
    axs{ii}     = gca;
end
unicb(axs);                                                                 % enforce uniform colorbar


% Search Information (group)                        @MCN need to average triangles
metric          = 'searchinformation';
for ii = 1 : Nwint
    itA         = AWiint(ii);
    strt        = ['Search Information (group): ' AW8s{itA}]; 
    connplot(A.grp.(metric){itA}, strt, pinfo, ppos)
    axs{ii}     = gca;
end
unicb(axs);


% Path Transitivity (group)
metric          = 'pathtransitivity';
for ii = 1 : Nwint
    itA         = AWiint(ii);
    strt        = ['Path Transitivity (group): ' AW8s{itA}]; 
    connplot(A.grp.(metric){itA}, strt, pinfo, ppos);
    axs{ii}     = gca;
end
unicb(axs);


%% Relationships with FC
jj                  = 1;
Nctx                = pinfo(jj).ncor;
Nstx                = pinfo(jj).nsub;
Nnode               = Nctx + Nstx;
Dsfc                = D{fci,jj};
Dgfc                = groupavg(Dsfc,3,'nz');
[it,cil,cirm]       = conn_useCIs(Dgfc,pinfo);                              % Get node order for communities from FC
Dgfcr               = Dgfc(it,it);                                          % reorder FC data
Dgfcrz              = zscore(Dgfcr);                                        % zscore so we can show in same matrix with metric
ltmask              = tril(ones(Nnode),-1)~=0;                              % mask for lower triangle
mdmask              = logical(eye(Nnode));                                  % mask for main diagonal
Ystr                = 'Linear correlation with FC';

datanames           = {'searchinformation' 'pathtransitivity' 'edgebetweenness'};
for dns = 1:length(datanames)
    metric          = datanames{dns};
    Dgrp            = A.grp.(metric)(AWiint,:);
    Dsub            = A.sub.(metric)(AWiint,:);
    Rg              = nan(1,Nwint);
    Rs              = nan(Nsubs,Nwint);
    pltstr          = {{['\fontsize{20}' upper(metric) ' from SC'] ; ['\fontsize{12}' '(subject & group)']}, Astr};
    
    % Distributions
    plot_conn_histNorm(Dsub,AW8s(AWiint),parcs(jj),str,cmapt,[ 1465 -310 1272 316])
    title([upper(metric) ' from SC (subject)']);

        
    % Visual comparison of Metric & FC in adjacency matrices
    for ii = 1 : Nwint
        itA             = AWiint(ii);                                       % Index for this weight in A
        strt            = ['FC vs ' upper(metric) '-' AW8s{itA} ' (zscore)'];
        Dg              = Dgrp{ii}(it,it);
        if strncmpi(metric,'edgebetweenness',6)
            Dg(Dg==0)   = nan;
            Dg          = log(Dg);
            metric      = [metric ' (log)'];
        end
        
        Dg(isnan(Dg))   = 0;
        Dg              = zscore(Dg);
        Dg(mdmask)      = nan;
        Dg(ltmask)      = Dgfcrz(ltmask);
        myfig(strt);
        imagesc(Dg); axis square; colorbar; title(strt); set(gca,'FontSize',fntsz); caxis([-4 4]);
        set(gca,'XTick',ceil(cirm'),'XTickLabel',cil','XTickLabelRotation',40)
        set(gca,'YTick',ceil(cirm'),'YTickLabel',cil','YTickLabelRotation',40)
    end

    % Compute group-level correlation with FC
    for ii = 1 : Nwint
        Dg              = Dgrp{ii};
        r1              = conn_corrmat(Dg,  Dgfc);                          % correlate lower tri with FC   @MCN BETTER TO AVERAGE TRIS IN SI???  
        r2              = conn_corrmat(Dg', Dgfc);                          % and upper tri
        Rg(ii)          = mean([r1 r2]);                                    % average r of both tris
    end
    
    % subject-level correlations
    for ii = 1 : Nwint
        Ds              = Dsub{ii};
        for ss = 1 : Nsubs
            ds          = Ds(:,:,ss);
            dsfc        = Dsfc(:,:,ss);
            r1          = conn_corrmat(ds,dsfc);                                % correlate lower tri with FC     
            r2          = conn_corrmat(ds',dsfc);                               % and upper tri
            Rs(ss,ii)   = mean([r1 r2]);                                        %      
        end
    end
    
    % Subject & Group violin
    plot_violin(Rs,pltW8s,Rg,Ystr,pltstr,parcs{1},[1442 465 608 332],cmapt) 
    
    
    
    % Heat scatters vs FC at group level
    for ii = 1 : Nwint
        itA             = AWiint(ii);
        W8stmp          = {W8s{fci}(1:2), AW8s{itA}};
        plot_conn_vsw8([{Dgfc};Dgrp(ii)], W8stmp, 'fc');                    % heat scatter vs FC
        set(gcf,'Position',[1442 377 568 420])
        set(gcf,'Name', ['FC vs ' upper(metric) ': group-level'])
    end 
   

    % Heat scatters vs FC at subject level (EXPENSIVE!)
    ii=2;
    Dtsub = Dsub{ii};
    Dtsubvec=[]; Dfcsubvec=[];
    for ss = 1 :50;
    Dtsc=Dtsub(:,:,ss); Dtfc=Dsfc(:,:,ss);
    Dtsubvec=[Dtsubvec Dtsc(ltmask)']; Dfcsubvec=[Dfcsubvec Dtfc(ltmask)'];
    end
    Dtsubvec(Dtsubvec==0)=nan; R=corr(Dtsubvec',Dfcsubvec','Type','Pearson','Rows','complete');
    Dtsubvec(Dtsubvec==0)=nan; Rho=corr(Dtsubvec',Dfcsubvec','Type','Spearman','Rows','complete');
    corrstr = sprintf('r = %.2f',R);
    
    myfig('',[1442 377 568 420]), heatscatter(Dtsubvec',Dfcsubvec')
    xlabel(AW8s(AWiint(ii))); ylabel('FC'); set(gca,'FontSize',25); 
    title({['\fontsize{25}' 'SC-' upper(metric) ' vs FC'], ['\fontsize{12}' '(pooled subject; ' corrstr ')']}); 
    % set(gcf,'Position',ppos);                     % Smaller
    % set(gcf,'Position',[1441 304 628 493]);       % larger
%     text(70,1.2,corrstr,'FontSize',25);
%     myfig, heatscatter(Dtsubvec(~isnan(Dtsubvec))',Dfcsubvec(~isnan(Dtsubvec))')
    
    
    % Variance
    Dstack = nan(Nnode,Nnode,length(Dgrp));
    for ii = 1 : Nwint
        Dstack(:,:,ii) = Dgrp{ii};                                          % Stack group level data to compare to subject    
    end
%     plot_conn_var({Dstack},{'Group'},'cqd',pinfo,[upper(metric) '; ' str]);
%     set(gcf,'Position',ppos)

    Dtmp    = Dsub;
    Dtmp    = [Dtmp; Dstack];                                               % Combine subject & Group level data
    lbltmp  = pltW8s;
    lbltmp  = [lbltmp 'Group'];                                             % Combine subject & Group level data
    cmaptmp = cmapt;
    cmaptmp = [cmaptmp; [0 0 0]];                                           % Add black color for group level
    plot_conn_var(Dtmp,lbltmp,'cova',pinfo,[upper(metric) '; ' str],cmaptmp);
    plot_conn_var(Dtmp,lbltmp,'cqd',pinfo,[upper(metric) '; ' str],cmaptmp);
    set(gcf,'Position',ppos)
end



% Difference of search information
metric = 'searchinformation';
% isel            = [3 4 12];
% isel            = [8 10 12];
isel            = [1 2 3];         
W8sel           = W8s(isel);
W8sel=strrep(W8sel,'(log)','');W8sel=strrep(W8sel,'-u','');W8sel=strrep(W8sel,'-s','');W8sel=strrep(W8sel,'-c','');W8sel=strrep(W8sel,'nos','NoS');
W8sel=strrep(W8sel,'icvf','ICVF');
Dgrp            = A.grp.(metric)(isel);
Dsub            = A.sub.(metric)(isel);

ii1=3; ii2=2;
test = Dgrp{ii1} - Dgrp{ii2};
strtmp = ['Difference SI ' W8sel{ii1} ' - ' W8sel{ii2} ' (group)'];
connplot(test,strtmp,pinfo)

ii1=3; ii2=1;
test = Dgrp{ii1} - Dgrp{ii2};
strtmp = ['Difference SI ' W8sel{ii1} ' - ' W8sel{ii2} ' (group)'];
connplot(test,strtmp,pinfo)

ii1=2; ii2=1;
test = Dgrp{ii1} - Dgrp{ii2};
strtmp = ['Difference SI ' W8sel{ii1} ' - ' W8sel{ii2} ' (group)'];
connplot(test,strtmp,pinfo)


%% ------------------- Surfaces & Cross-subject variance --------------------- %%
% Network properties
% A.grp. { 'strength'  'clusteringcoeff'  'nodebetweenness'  'efficiency_local'  'participationcoeff' }
% A.grp.charpath. { 'eccentricity' } 

% Setup
isel                        = [1 3];
W8sel                       = W8s(isel);
nNode                       = [pinfo.ncor];
conte69is                   = [pinfo.conte69i];
parcs                       = {pinfo.name};
parcs                       = strrep(parcs,'-',' ');
% [srf_lh,srf_rh]             = load_conte69();                               % load surface from brainspace

metric      = 'strength';
pltmetric  = upper(strrep(metric,'_','-'));
drank       = nan(numel(isel),nNode);

% Rank data
it      = 0;                                                                % just an iterator
for ii = isel
    it                  = it + 1;
    d                   = A.grp.(metric){ii}';
    d(1:14)             = [];
%     d                 = -1 * A.grp.charpath.(metric){ii}';
    [drank(it,:), ~]  = tiedrank(d);                                              % rank data
end

% Convert rankings to consistency measure
drank(drank<390)    = 0;                                                    % focus on top 10 nodes
drank(drank~=0)     = 1;                                                    % binarize top nodes
drank               = sum(drank,1);                                         % sum over weights to get consistency ranking

% Convert to full surface
dfull               = parcel2full(drank',conte69is(:,1));

% Plot on surface
h=plot_hemispheres(dfull,{srf_lh,srf_rh}, 'labeltext',pltmetric);
h.figure.NumberTitle='off';
h.figure.Name=[pltmetric ' top 50 ranked nodes: ' strin ', ' parcs{1}];
% h.figure.Name=[pltmetric ' reverse ranked data: ' strin ', ' parcs{1}];
h.cb.Ticks  = [1 2 3];
modmap      = zeros(4,3);           % custom colormap
modmap(1,:) = [0.2081 0.1663 0.5292]; % dark blue
modmap(2,:) = [0 1 1];              % cyan
modmap(3,:) = [1 1 0];              % yellow
modmap(4,:) = [1 0 0];              % red
for ax =1:4
    colormap(h.axes(ax),modmap)
end

% Set common color mapping
% cmax = max(cellfun(@(c) max(c), dfull(subsel)));
% cmin = min(cellfun(@(c) min(c), dfull(subsel)));
cmax = 400; cmin = 350;
for ax =1:4
    h.axes(ax,1).CLim=[cmin, cmax];
    h.axes(ax,2).CLim=[cmin, cmax];
end
h.cb(1).Ticks = [cmin, cmax];
h.cb(2).Ticks = [cmin, cmax];


% A.grp. { 'strength'  'clusteringcoeff'  'nodebetweenness'  'efficiency_local'  'participationcoeff' }
% A.grp.charpath. { 'eccentricity' } 
% Cross-subject variance on surface
isel        = [3 4 12];
W8sel       = W8s(isel);
metric      = 'strength';
pltmetric   = upper(strrep(metric,'_','-'));
dfull       = cell(numel(isel),1);
dvar        = nan(400,numel(isel));
Dsub        = A.sub.(metric)(isel);
% Dsub        = A.sub.(metric);


% Convert nodal data to full surface
for ii = 1 : numel(isel)
    d           = Dsub{ii};
    dvar(:,ii)  = std(d).^2;                                                % compute cross-subject variance at each node
    dfull{ii}   = parcel2full(dvar(:,ii),conte69is(:,1));
end

% Plot nodal data on surface
subsel = [1 2];
D1=dfull{subsel(1),1};
D2=dfull{subsel(2),1};
h=plot_hemispheres([D1,D2],{srf_lh,srf_rh}, 'labeltext',{W8sel{subsel(1)}, W8sel{subsel(2)}});
h.figure.NumberTitle='off';
h.figure.Name=['Cross-subject variance in ' pltmetric ': ' strin ', ' parcs{1}];

subsel = [3 3];
D1=dfull{subsel(1),1};
D2=dfull{subsel(2),1};
h=plot_hemispheres([D1,D2],{srf_lh,srf_rh}, 'labeltext',{W8sel{subsel(1)}, W8sel{subsel(2)}});
h.figure.NumberTitle='off';
h.figure.Name=['Cross-subject variance in ' pltmetric ': ' strin ', ' parcs{1}];

% violin plot of variance
T       = array2table(dvar);
plotstr = [parcs{1} ', ' strin];
myfig(plotstr,[1442 465 608 332])
violinplot(T, 'ShowMean'); hold on
ylabel('Cross-subject variance');
title(pltmetric)
set(gca, 'XTickLabel',W8sel, 'XTickLabelRotation',20, 'FontSize', fntsz)


%% TESTING
% High level stats for some metric
in=7;
means=nan(1,in); maxs=nan(1,in); stds=nan(1,in);
metric = 'edgebetweenness';
Nnodes = length(pinfo.labels);                                     % node count (havent extended this section to handle multiple parcs)
for ii = 1:in 
    itA     = AWiint(ii); 
    dt      = A.grp.(metric){itA};
    E       = find(dt);
    dt(E)   = dt(E)./((Nnodes-1)*(Nnodes-2));                               % Rescale to proportion of total possible paths
    means(ii)=mean(dt(:));
    maxs(ii)=max(dt(:));
    stds(ii)=std(dt(:));
end

myfig, 

%--------------------------------------------------------------------------
end
