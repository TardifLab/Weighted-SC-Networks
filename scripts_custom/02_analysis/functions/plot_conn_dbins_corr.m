function plot_conn_dbins_corr(D,W8s,varargin)
% Computes & plots a range of metrics across binned data.
% This variant of plot_conn_dbins_* only compares 2 datasets and runs
% correlations.
%
% * WARNING * Produces a large number of plots!
%
% Input
%                           + Required +
%   D       : 1X3 cell array, each cell containing an NxNxS stack of
%             connectivity matrices (S=subs) 
%   W8s     : 1X3 cell array of strings describing edge weights in D
%
%                           + Optional +
%   dbin    : character vector indicating which dataset in D to use to
%             discretize the remaining (default = 1st dataset)
%   bins    : number of bins for discretizing data (integer, default=10)
%   str     : Info for plot title
%   cmap    : Ix3 double used as color map
%   ppos    : plot size & position
%
%
% 2022 Mark C Nelson, McConnell Brain Imaging Centre, MNI, McGill
%--------------------------------------------------------------------------

[I,J]                       = size(D);

%% Optional inputs
nin                         = max(nargin,1) - 2;
defaults                    = {W8s{1},10,'',[],[-2559 797 494 411]};
[dbin,bins,str,cmap,ppos]   = INhandler(varargin,nin,defaults);

if isempty(cmap); cmap=colormap(hsv(I)); end

%% Setup
% opt                         = {'rows','complete','type','Spearman'};        % options used in correlations
opt                         = {'rows','complete','type','Pearson'};        % options used in correlations
fntsz                       = 25;
Idbin                       = cellstrfind(W8s,dbin);                        % index of dataset used to discretize
Iothr                       = setdiff(1:I,Idbin);
cmapdbin                    = cmap(Idbin,:);      
cmapothr                    = cmap(Iothr,:);
bincountg                   = zeros(1,bins);                                % Count data points in each bin
bincounts                   = zeros(1,bins);
bcstrg                      = cell(1,bins);                                 % store strings of bin counts for legend
bcstrs                      = cell(1,bins);
cmapbins                    = colormap(hsv(bins));

% String used in plots depends on correlation type
if strncmpi(opt{4},'spear',5)
    corrstr = '\rho';
else
    corrstr = 'r';
end
                    

%% Discretize the indicated dataset
Dbin                        = D{Idbin};
mask                        = logical(tril(ones(size(Dbin,1)),-1));         % Lower tri LESS main diagonal (assuming symmetry!)

% Group level
Dbing                       = groupavg(Dbin,3,'nz');                        % Group average (assumes 3rd dimension is subjects!)
Dbingvec                    = Dbing(mask);                                  % Vectorize
Dbingmin                    = min(Dbingvec(Dbingvec~=0));
Dbingmax                    = max(Dbingvec(Dbingvec~=0));
edgesg                      = linspace(Dbingmin,Dbingmax,bins+1);           % edges = bins + 1
Ibing                       = discretize(Dbingvec,edgesg);                  % bin indices used to discretize data (0s left out)


% Subject level
Nsub                        = size(Dbin,3);
Dbinsvec                    = zeros(length(Dbingvec),Nsub);
for ss = 1 : Nsub
    dts                     = Dbin(:,:,ss);
    Dbinsvec(:,ss)          = dts(mask);
end
Dbinsvec                    = Dbinsvec(:);
Dbinsmin                    = min(Dbinsvec(Dbinsvec~=0));
Dbinsmax                    = max(Dbinsvec(Dbinsvec~=0));
edgess                      = linspace(Dbinsmin,Dbinsmax,bins+1);
Ibins                       = discretize(Dbinsvec,edgess);                  % bin indices used to discretize data

% Bin labels used in plots
binlblg                             = {};
binlbls                             = {};
if strncmpi(dbin,'los',3)
    switch bins
        case 2
            binlblg{1}              = 'shorter'; 
            binlblg{2}              = 'longer';
        case 3
            binlblg{1}              = 'short'; 
            binlblg{2}              = 'med'; 
            binlblg{3}              = 'long';
        case 4
            binlblg{1}              = 'short';   
            binlblg{2}              = 'med-short'; 
            binlblg{3}              = 'med-long';       
            binlblg{4}              = 'long';
        case 5
            binlblg{1}              = 'short';   
            binlblg{2}              = 'med-short'; 
            binlblg{3}              = 'medium'; 
            binlblg{4}              = 'med-long';       
            binlblg{5}              = 'long';
        case 6
            binlblg{1}              = 'short';   
            binlblg{2}              = 'med-short'; 
            binlblg{3}              = 'medium'; 
            binlblg{4}              = 'medium'; 
            binlblg{5}              = 'med-long';
            binlblg{6}              = 'long';
        otherwise
            binlblg{1}              = 'short';   
            binlblg{bins}           = 'long';
                for bb = 2 : bins-1                        
                    binlblg{bb}     = ['med ' num2str(bb)];
                end
    end
    binlbls                         = binlblg;
else
    % Just use the edge values and number the bins
    for bb = 1 : bins
        binlblg{bb} = ['bin ' num2str(bb) ' ' sprintf('(%.3f)', edgesg(bb))];
        binlbls{bb} = ['bin ' num2str(bb) ' ' sprintf('(%.3f)', edgess(bb))];
    end
end


%% Vectorize other datasets
Dgall                   = nan(length(Ibing),2);
Dsall                   = nan(length(Ibins),2);
for ii = 1 : 2
    Dt                  = D{Iothr(ii)};
    
    % Group Level
    Dgt                 = groupavg(Dt,3,'nz');
    Dgtvec              = Dgt(mask);
    Dgtvec(Dgtvec==0)   = nan;
    Dgall(:,ii)         = Dgtvec;
    
    % Subject Level
    Dtsvec              = zeros(length(Dbingvec),Nsub);
    for ss = 1 : Nsub
        dts             = Dt(:,:,ss);
        Dtsvec(:,ss)    = dts(mask);
    end
    Dtsvec(Dtsvec==0)   = nan;
    Dsall(:,ii)         = Dtsvec(:);
end


%% Get some info about these datasets across all edges
allmedng            = nan(1,2);
allmedns            = nan(1,2);
allcqdg             = nan(1,2);
allcqds             = nan(1,2);


for ii = 1 : 2    
    % Median
    allmedng(ii)     = nanmedian(Dgall(:,ii));                              % Group
    allmedns(ii)     = nanmedian(Dsall(:,ii));                              % Subject
    
    % CQD
    EQ1             = quantile(Dgall(:,ii),.25,1); 
    EQ3             = quantile(Dgall(:,ii),.75,1); 
    allcqdg(ii)      = (EQ3-EQ1) ./ (EQ3+EQ1);                              % Group
    EQ1             = quantile(Dsall(:,ii),.25,1); 
    EQ3             = quantile(Dsall(:,ii),.75,1); 
    allcqds(ii)      = (EQ3-EQ1) ./ (EQ3+EQ1);                              % Subject
end

% Correlation with each other
allcorrg            = corr(Dgall(:,1),Dgall(:,2),opt{:});                   % Group
allcorrs            = corr(Dsall(:,1),Dsall(:,2),opt{:});                   % Subject

% Correlation with the dataset used to bin the data (typically edge length)
allcorrg_Dbin       = zeros(1,I-1);
allcorrs_Dbin       = zeros(1,I-1);
allcorrg_Dbin(1)    = corr(Dgall(:,1),Dbingvec,opt{:});                     % Group
allcorrg_Dbin(2)    = corr(Dgall(:,2),Dbingvec,opt{:});
allcorrs_Dbin(1)    = corr(Dsall(:,1),Dbinsvec,opt{:});                     % Subject
allcorrs_Dbin(2)    = corr(Dsall(:,2),Dbinsvec,opt{:});


%% Compute some metrics within these bins
bincorrg                 = zeros(1,bins);
bincorrg_Dbin            = zeros(2,bins);
binmedng                 = zeros(2,bins);
bincqdg                  = zeros(2,bins);
bincorrs                 = zeros(1,bins);
bincorrs_Dbin            = zeros(2,bins);
binmedns                 = zeros(2,bins);
bincqds                  = zeros(2,bins);

D1g=Dgall(:,1); D2g=Dgall(:,2);
D1s=Dsall(:,1); D2s=Dsall(:,2);
for bb = 1 : bins
    % Group Level
    bincountg(bb)           = length(find(Ibing==bb));                      % number of data points in each bin
    bcstrg{bb}              = num2str(bincountg(bb));                       % used in plot legends
    if ~isempty(find(Ibing==bb,1))                                          % if bin is not empty
        % Data for this bin
        dt1                 = D1g(Ibing==bb);
        dt2                 = D2g(Ibing==bb);
        dt_Dbin             = Dbingvec(Ibing==bb);
        
        % Quantiative metrics
        bincorrg(bb)        = corr(dt1,dt2,opt{:});                         % Corr
        bincorrg_Dbin(1,bb) = corr(dt1,dt_Dbin,opt{:}); 
        bincorrg_Dbin(2,bb) = corr(dt2,dt_Dbin,opt{:});
        binmedng(1,bb)      = nanmedian(dt1);                               % Median
        binmedng(2,bb)      = nanmedian(dt2);
        EQ1=quantile(dt1,.25,1); EQ3=quantile(dt1,.75,1); bincqdg(1,bb)=(EQ3-EQ1) ./ (EQ3+EQ1); % Edge Weight CQD
        EQ1=quantile(dt2,.25,1); EQ3=quantile(dt2,.75,1); bincqdg(2,bb)=(EQ3-EQ1) ./ (EQ3+EQ1);
    end
    
    % Subject Level
    bincounts(bb)           = length(find(Ibins==bb));
    bcstrs{bb}              = num2str(bincounts(bb));
    if ~isempty(find(Ibins==bb,1))                                          % if bin is not empty
        % Data for this bin
        dt1                 = D1s(Ibins==bb);
        dt2                 = D2s(Ibins==bb);
        dt_Dbin             = Dbinsvec(Ibins==bb);
        
        % Metrics
        bincorrs(bb)        = corr(dt1,dt2,opt{:});                         % Corr
        bincorrs_Dbin(1,bb) = corr(dt1,dt_Dbin,opt{:}); 
        bincorrs_Dbin(2,bb) = corr(dt2,dt_Dbin,opt{:});
        binmedns(1,bb)      = nanmedian(dt1);                               % Median
        binmedns(2,bb)      = nanmedian(dt2);
        EQ1=quantile(dt1,.25,1); EQ3=quantile(dt1,.75,1); bincqds(1,bb)=(EQ3-EQ1) ./ (EQ3+EQ1); % CQD
        EQ1=quantile(dt2,.25,1); EQ3=quantile(dt2,.75,1); bincqds(2,bb)=(EQ3-EQ1) ./ (EQ3+EQ1);
    end

end

%% ----------------- Plots --------------- %%

%% Histograms of binned dataset
% Group
myfig(['BinnedData-Group-' str],ppos); histogram(Dbingvec(Dbingvec~=0),100); hold on
for ee=1:bins;k=edgesg(ee);L(ee)=line([k k],ylim,'LineStyle', '--', 'Color', cmapbins(ee,:),'LineWidth',3); end
ylabel('Count'); xlabel(W8s{Idbin}); title('Binned data (group)');
set(gca,'FontSize',fntsz); legend(L,strcat('N = ',bcstrg));

% Subject
myfig(['BinnedData-Subject-' str],ppos); histogram(Dbinsvec(Dbinsvec~=0),100); hold on
for ee=1:bins;k=edgess(ee);L(ee)=line([k k],ylim,'LineStyle', '--', 'Color', cmapbins(ee,:),'LineWidth',3); end
ylabel('Count'); xlabel(W8s{Idbin}); title('Binned data (subject)');
set(gca,'FontSize',fntsz); legend(L,strcat('N = ',bcstrs)); 

%% Heat scatter plots of binned data vs LoS
 
for ii = 1 : I-1
    Y   = Dgall(:,ii); 
    X   = Dbingvec;
    Y   = Y(X~=0);   
    X   = X(X~=0);
    X   = X(~isnan(Y));
    Y   = Y(~isnan(Y));   
    

    % Log transform data if skewed
    if abs(skewness(Y)) > 3    
        Y       = log(Y);
        ylbl    = [W8s{Iothr(ii)} ' (log)'];
    else
        ylbl    = W8s{Iothr(ii)};
    end

    % Plot Group level
    r       = corr(X,Y,opt{:}); 
    rstr    = sprintf(' = %.3f', r);                                        % Compute correlation to display on plot
    myfig(['Group-' str],ppos); heatscatter(X,Y); hold on; xlims=ylim;
    for ee = 1 : bins
        k   = edgesg(ee);
        L(ee)=line([k k],ylim,'LineStyle', '--', 'Color', cmapbins(ee,:),'LineWidth',3); % Mark position of bins
    end 
    ylim(xlims); xlabel(W8s{Idbin}); ylabel(ylbl); title([corrstr rstr]); set(gca,'FontSize',fntsz);

    % Repeat at subject level (SLOW)
% % %     Y           = Dsall(:,ii); 
% % %     X           = Dbinsvec;
% % %     Y           = Y(X~=0);    
% % %     X           = X(X~=0);
% % %     if abs(skewness(Y))>5; Y=log(Y); end
% % %     r           = corr(X,Y,opt{:}); 
% % %     rstr        = sprintf(' = %.3f', r);
% % %     myfig(['Subject-' str],ppos), heatscatter(X,Y); hold on; xlims=ylim;
% % %     for ee = 1 : bins
% % %         k       = edgess(ee);
% % %         L(ee)   = line([k k],ylim,'LineStyle', '--', 'Color', cmapbins(ee,:),'LineWidth',3); 
% % %     end 
% % %     ylim(xlims); xlabel(W8s{Idbin}); ylabel(ylbl); title([corrstr rstr]); set(gca,'FontSize',fntsz); xlim([0 1.2]);
    
end


%% Scatter plots of binned data vs each other with points colored by binID
% Group
X           = Dgall(:,2);       
Y           = Dgall(:,1);

if abs(skewness(Y)) > 3    
    Y       = log(Y);
    ylbl    = [W8s{Iothr(1)} ' (log)'];
else
    ylbl    = W8s{Iothr(1)};
end
if abs(skewness(X)) > 3    
    X       = log(X);
    xlbl    = [W8s{Iothr(2)} ' (log)'];
else
    xlbl    = W8s{Iothr(2)};
end

r=corr(X,Y,opt{:}); rstr=sprintf(' = %.3f', r);
myfig(['Group-' str],[-2559 434 993 774]);
for bb = 1 : bins
    bcstr{bb} = num2str(bincountg(bb));
    plot(X(Ibing==bb),Y(Ibing==bb),'k.','MarkerEdgeColor',cmapbins(bb,:),'MarkerSize',12); hold on
end
ylabel(ylbl); xlabel(xlbl); title(['Group: ' corrstr rstr]); set(gca,'FontSize',fntsz); legend(binlblg);



% Subject (SLOW!)
% % % X   = Dsall(:,2);
% % % Y   = Dsall(:,1);
% % % if abs(skewness(Y))>3; Y=log(Y); end
% % % if abs(skewness(X))>3; X=log(X); end
% % % r=corr(X,Y,opt{:}); rstr=sprintf(' = %.3f', r);
% % % myfig(['Subject-' str],[-2559 434 993 774]); cmapbins=colormap(hsv(bins));
% % % for bb = 1 : bins
% % %     bcstr{bb} = num2str(bincounts(bb));
% % % %     plot(X(Ibins==bb),Y(Ibins==bb),'ko','MarkerEdgeColor', cmapbins(bb,:),'MarkerFaceColor',cmapbins(bb,:),'MarkerSize',10); hold on
% % % %     plot(X(Ibins==bb),Y(Ibins==bb),'ko','MarkerEdgeColor', cmapbins(bb,:),'MarkerSize',10,'LineWidth',2); hold on
% % %     plot(X(Ibing==bb),Y(Ibing==bb),'k.','MarkerEdgeColor',cmapbins(bb,:),'MarkerSize',12); hold on
% % % end
% % % ylabel(ylbl); xlabel(xlbl); title(['Subject: ' corrstr rstr]); set(gca,'FontSize',fntsz); legend(binlbls);



%% -------------------- LINE PLOTS ------------------------- %%
ppost   = [-2559 694 783 514];
xlims   = [.5 bins+.5];
pltopt  = {'LineWidth',.1,'MarkerSize',20};
lineopt = {'LineStyle', '--','LineWidth',3};
gcaopt  = {'FontSize',fntsz,'XTick',1:bins,'XTickLabel',binlbls,'XTickLabelRotation',30};


%% Correlations

% Plot correlation between weights
dplt    = bincorrg; dplt2   = bincorrs;
dline   = allcorrg; dline2  = allcorrs;
pltstr  = [W8s{Iothr(1)} ' vs ' W8s{Iothr(2)}];
lgndstr = {'bin-grp' 'bin-sub' 'all-grp' 'all-sub'};
% Plot
myfig([pltstr '-' str],ppost); 
P1=plot(1:bins,dplt, ':d','Color','b','MarkerFaceColor','b',pltopt{:}); hold on
P2=plot(1:bins,dplt2,':s','Color','g','MarkerFaceColor','g',pltopt{:}); 
xlim(xlims);
Rg=line(xlim,[dline  dline], 'Color', 'b',lineopt{:});
Rs=line(xlim,[dline2 dline2],'Color', 'g',lineopt{:});
title(pltstr); ylabel([opt{4} '''s ' corrstr]); set(gca,gcaopt{:});
legend([P1 P2 Rg Rs],lgndstr);
ylims = ylim;
if all(0>ylims(1) & 0<ylims(2)) || all(0<ylims(1) & 0>ylims(2))             % 0 lies between color limits
    line(xlim,[0 0],'LineStyle', '--', 'Color', 'r','LineWidth',2);         % 0 line
end


% GROUP correlation between both weights & the data used to define bins
dplt    = bincorrg_Dbin;
dline   = allcorrg_Dbin;
pltstr  = ['Correlation with ' W8s{Idbin} ' (group)'];
lgndstr = {[W8s{Iothr(1)} ' (bin)'] [W8s{Iothr(2)} ' (bin)'] [W8s{Iothr(1)} ' (all)'] [W8s{Iothr(2)} ' (all)']};
% Plot
myfig([pltstr '-' str],ppost); 
P1=plot(1:bins,dplt(1,:),':d','Color',cmapothr(1,:),'MarkerFaceColor',cmapothr(1,:),pltopt{:}); hold on
P2=plot(1:bins,dplt(2,:),':s','Color',cmapothr(2,:),'MarkerFaceColor',cmapothr(2,:),pltopt{:}); 
xlim(xlims);
Rw1=line(xlim,[dline(1) dline(1)],'Color', cmapothr(1,:),lineopt{:});
Rw2=line(xlim,[dline(2) dline(2)],'Color', cmapothr(2,:),lineopt{:});
title(pltstr); ylabel([opt{4} '''s ' corrstr]); set(gca,gcaopt{:});
legend([P1 P2 Rw1 Rw2],lgndstr);
ylims = ylim;
if all(0>ylims(1) & 0<ylims(2)) || all(0<ylims(1) & 0>ylims(2))
    line(xlim,[0 0],'LineStyle', '-', 'Color', 'k','LineWidth',2);
end

% SUBJECT correlation between both weights & the data used to define bins
dplt    = bincorrs_Dbin;
dline   = allcorrs_Dbin;
pltstr  = ['Correlation with ' W8s{Idbin} ' (subject)'];
lgndstr = {[W8s{Iothr(1)} ' (bin)'] [W8s{Iothr(2)} ' (bin)'] [W8s{Iothr(1)} ' (all)'] [W8s{Iothr(2)} ' (all)']};
% Plot
myfig([pltstr '-' str],ppost); 
P1=plot(1:bins,dplt(1,:),':d','Color',cmapothr(1,:),'MarkerFaceColor',cmapothr(1,:),pltopt{:}); hold on
P2=plot(1:bins,dplt(2,:),':s','Color',cmapothr(2,:),'MarkerFaceColor',cmapothr(2,:),pltopt{:}); 
xlim(xlims);
Rw1=line(xlim,[dline(1) dline(1)],'Color', cmapothr(1,:),lineopt{:});
Rw2=line(xlim,[dline(2) dline(2)],'Color', cmapothr(2,:),lineopt{:});
title(pltstr); ylabel([opt{4} '''s ' corrstr]); set(gca,gcaopt{:});
legend([P1 P2 Rw1 Rw2],lgndstr);
ylims = ylim;
if all(0>ylims(1) & 0<ylims(2)) || all(0<ylims(1) & 0>ylims(2))
    line(xlim,[0 0],'LineStyle', '-', 'Color', 'k','LineWidth',2);
end



%% Variance

% Plot GROUP variability in both weights
dplt    = bincqdg;
dline   = allcqdg;
pltstr  = ['Group Edge Weight Variance across '  W8s{Idbin} ' bins'];
lgndstr = {[W8s{Iothr(1)} ' (bin)'] [W8s{Iothr(2)} ' (bin)'] [W8s{Iothr(1)} ' (all)'] [W8s{Iothr(2)} ' (all)']};
ylbl    = 'CQD';
C1      = cmapothr(1,:);
C2      = cmapothr(2,:);
myfig([pltstr '-' str],ppost);
P1=plot(1:bins,dplt(1,:),':d','Color',C1,'MarkerFaceColor',C1,pltopt{:}); hold on
P2=plot(1:bins,dplt(2,:),':s','Color',C2,'MarkerFaceColor',C2,pltopt{:}); 
xlim(xlims);
Rw1=line(xlim,[dline(1) dline(1)],'Color',C1,lineopt{:});
Rw2=line(xlim,[dline(2) dline(2)],'Color',C2,lineopt{:});
title(pltstr); ylabel(ylbl); set(gca,gcaopt{:});
legend([P1 P2 Rw1 Rw2],lgndstr);
ylims = ylim;
if all(0>ylims(1) & 0<ylims(2)) || all(0<ylims(1) & 0>ylims(2))
    line(xlim,[0 0],'LineStyle', '-', 'Color', 'k','LineWidth',2);
end


% Plot SUBJECT variability in both weights
dplt    = bincqds;
dline   = allcqds;
pltstr  = ['Subject Edge Weight Variance across '  W8s{Idbin} ' bins'];
lgndstr = {[W8s{Iothr(1)} ' (bin)'] [W8s{Iothr(2)} ' (bin)'] [W8s{Iothr(1)} ' (all)'] [W8s{Iothr(2)} ' (all)']};
ylbl    = 'CQD';
myfig([pltstr '-' str],ppost); 
P1=plot(1:bins,dplt(1,:),':d','Color',cmapothr(1,:),'MarkerFaceColor',cmapothr(1,:),pltopt{:}); hold on
P2=plot(1:bins,dplt(2,:),':s','Color',cmapothr(2,:),'MarkerFaceColor',cmapothr(2,:),pltopt{:}); 
xlim(xlims);
Rw1=line(xlim,[dline(1) dline(1)],'Color', cmapothr(1,:),lineopt{:});
Rw2=line(xlim,[dline(2) dline(2)],'Color', cmapothr(2,:),lineopt{:});
title(pltstr); ylabel(ylbl); set(gca,gcaopt{:});
legend([P1 P2 Rw1 Rw2],lgndstr);
ylims = ylim;
if all(0>ylims(1) & 0<ylims(2)) || all(0<ylims(1) & 0>ylims(2))
    line(xlim,[0 0],'LineStyle', '-', 'Color', 'k','LineWidth',2);
end


%--------------------------------------------------------------------------
end