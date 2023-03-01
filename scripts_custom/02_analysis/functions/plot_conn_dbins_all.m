function plot_conn_dbins_all(D,W8s,varargin)
% Computes & plots a range of metrics across binned data.
% This variant of plot_conn_dbins_* compares all datasets given and thus
% does not run correlations, but computes all metrics within each network. 
%
% * WARNING * Produces a large number of plots!
%
% Input
%                           + Required +
%   D       : 1XN cell array, each cell containing an NxNxS stack of
%             connectivity matrices (S=subs) 
%   W8s     : 1XN cell array of strings describing edge weights in D
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
defaults                    = {W8s{1},10,'',colormap(hsv(I)),[-2559 797 494 411]};
[dbin,bins,str,cmap,ppos]   = INhandler(varargin,nin,defaults);

%% Setup
metric                      = 'cqd';
opt                         = {'rows','complete','type','Spearman'};        % options used in correlations
fntsz                       = 25;
Idbin                       = cellstrfind(W8s,dbin);                        % index of dataset used to discretize
Iothr                       = setdiff(1:I,Idbin);
N                           = numel(Iothr);
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

% Set metric specific labels
switch lower(metric)
    case 'cova'        
        metricstrS = 'CoVa';
    case 'cqd'
        metricstrS = 'CQD';
    otherwise; error('Unable to use your input to METRIC')
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
Dgall                   = nan(length(Ibing),N);
Dsall                   = nan(length(Ibins),N);
for ii = 1 : N
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
allmedng            = nan(1,N);
allmedns            = nan(1,N);
allvarg             = nan(1,N);
allvars             = nan(1,N);


for ii = 1 : N    
    % Median
    allmedng(ii)     = nanmedian(Dgall(:,ii));                              % Group
    allmedns(ii)     = nanmedian(Dsall(:,ii));                              % Subject
        
    
    switch lower(metric)
        % Coefficient of Variation (std/mu)
        case 'cova'                                                     
            allvarg(ii)     = nanstd(Dgall(:,ii)) ./ nanmean(Dgall(:,ii));
            allvars(ii)     = nanstd(Dsall(:,ii)) ./ nanmean(Dsall(:,ii));

        % Quartile Coefficient of Dispersion (Q3 - Q1)/(Q3 + Q1)
        case 'cqd'           
            EQ1             = quantile(Dgall(:,ii),.25,1); 
            EQ3             = quantile(Dgall(:,ii),.75,1); 
            allvarg(ii)     = (EQ3-EQ1) ./ (EQ3+EQ1);                       % Group
            EQ1             = quantile(Dsall(:,ii),.25,1); 
            EQ3             = quantile(Dsall(:,ii),.75,1); 
            allvars(ii)     = (EQ3-EQ1) ./ (EQ3+EQ1);                       % Subject
    end
end



%% Compute some metrics within these bins
binmedng                    = zeros(N,bins);
binvarg                     = zeros(N,bins);
binmedns                    = zeros(N,bins);
binvars                     = zeros(N,bins);

for bb = 1 : bins
    
    % Group Level
    bincountg(bb)           = length(find(Ibing==bb));                      % number of data points in each bin
    bcstrg{bb}              = num2str(bincountg(bb));                       % used in plot legends
    if ~isempty(find(Ibing==bb,1))                                          % if bin is not empty
                
        for ii = 1 : N
            
            % Get data for this bin
            DGALL                   = Dgall(:,ii);
            dt                      = DGALL(Ibing==bb);

            % Median
            binmedng(ii,bb)         = nanmedian(dt);
            
            % Variance
            switch lower(metric)
                % Coefficient of Variation
                case 'cova'                                                     
                    binvarg(ii,bb)  = nanstd(dt) ./ nanmean(dt);

                % Quartile Coefficient of Dispersion (Q3 - Q1)/(Q3 + Q1)
                case 'cqd'           
                    EQ1=quantile(dt,.25,1); EQ3=quantile(dt,.75,1); binvarg(ii,bb)=(EQ3-EQ1) ./ (EQ3+EQ1);
            end
        end
    end
    
    % Subject Level
    bincounts(bb)               = length(find(Ibins==bb));
    bcstrs{bb}                  = num2str(bincounts(bb));
    if ~isempty(find(Ibins==bb,1))                                          % if bin is not empty
        
        for ii = 1 : N
            
            % Data for this bin
            DSALL                   = Dsall(:,ii);
            dt                      = DSALL(Ibins==bb);

            % Median
            binmedns(ii,bb)         = nanmedian(dt);
            
            % Variance
            switch lower(metric)
                % Coefficient of Variation
                case 'cova'                                                     
                    binvars(ii,bb)  = nanstd(dt) ./ nanmean(dt);

                % Quartile Coefficient of Dispersion (Q3 - Q1)/(Q3 + Q1)
                case 'cqd'           
                    EQ1=quantile(dt,.25,1); EQ3=quantile(dt,.75,1); binvars(ii,bb)=(EQ3-EQ1) ./ (EQ3+EQ1);
            end
        end
    end
end

%% ----------------- Plots --------------- %%

%% Histograms of binned dataset
% Group
myfig(['BinnedData-Group-' str],ppos); histogram(Dbingvec(Dbingvec~=0),100); hold on
for ee=1:bins;k=edgesg(ee);L(ee)=line([k k],ylim,'LineStyle', '--', 'Color', cmapbins(ee,:),'LineWidth',3); end
ylabel('Count'); xlabel(W8s{Idbin}); title('Binned data (group)');
set(gca,'FontSize',fntsz); legend(L,strcat('N = ',bcstrg)); % xlim([0 1]);

% Subject
myfig(['BinnedData-Subject-' str],ppos); histogram(Dbinsvec(Dbinsvec~=0),100); hold on
for ee=1:bins;k=edgess(ee);L(ee)=line([k k],ylim,'LineStyle', '--', 'Color', cmapbins(ee,:),'LineWidth',3); end
ylabel('Count'); xlabel(W8s{Idbin}); title('Binned data (subject)');
set(gca,'FontSize',fntsz); legend(L,strcat('N = ',bcstrs));
%  xlim([0 1]);

%% Heat scatter plots of binned data vs LoS
 
for ii = 1 : N
    Y   = Dgall(:,ii); 
    X   = Dbingvec;
    Y   = Y(X~=0);   
    X   = X(X~=0);

    % Log transform data if skewed
    if abs(skewness(Y)) > 5    
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
    
end



%% -------------------- LINE PLOTS ------------------------- %%
ppost   = [-2559 694 783 514];
xlims   = [.5 bins+.5];
pltopt  = {'LineWidth',2,'MarkerSize',15};
% lineopt = {'LineStyle', '--','LineWidth',3};
gcaopt  = {'FontSize',fntsz,'XTick',1:bins,'XTickLabel',binlbls,'XTickLabelRotation',30};


%% Variance

% Plot GROUP variability in all weights
dplt    = binvarg;
pltstr  = ['Group Edge Weight Variance across '  W8s{Idbin} ' bins'];
ylbl    = ['Percent Max ' metricstrS];
myfig([pltstr '-' str],ppost);
for ii = 1 : N
    dt  = dplt(ii,:);
    dt  = dt./max(dt);                                                      % Convert to proportion to handle differences in scale
    plot(1:bins,dt,':d','Color',cmapothr(ii,:),'MarkerFaceColor',cmapothr(ii,:),pltopt{:}); hold on    
end
xlim(xlims); title(pltstr); ylabel(ylbl); set(gca,gcaopt{:});
legend(W8s(Iothr));
ylims = ylim;
if all(0>ylims(1) & 0<ylims(2)) || all(0<ylims(1) & 0>ylims(2))
    line(xlim,[0 0],'LineStyle', '-', 'Color', 'k','LineWidth',2);
end


% Plot SUBJECT variability in both weights
dplt    = binvars;
pltstr  = ['Subject Edge Weight Variance across '  W8s{Idbin} ' bins'];
ylbl    = ['Percent Max ' metricstrS];
myfig([pltstr '-' str],ppost); 
for ii = 1 : N
    dt  = dplt(ii,:);
    dt  = dt./max(dt);
    plot(1:bins,dt,':d','Color',cmapothr(ii,:),'MarkerFaceColor',cmapothr(ii,:),pltopt{:}); hold on
end
xlim(xlims); title(pltstr); ylabel(ylbl); set(gca,gcaopt{:});
legend(W8s(Iothr));
ylims = ylim;
if all(0>ylims(1) & 0<ylims(2)) || all(0<ylims(1) & 0>ylims(2))
    line(xlim,[0 0],'LineStyle', '-', 'Color', 'k','LineWidth',2);
end



% Plot SUBJECT & GROUP together
pltstr      = ['Edge Weight Variance across '  W8s{Idbin} ' bins'];
ylbl    = ['Percent Max ' metricstrS];
myfig([pltstr '-' str],ppost);

dplt        = binvarg;
pp          = 1;
for ii = 1 : N
    dt      = dplt(ii,:);
    dt      = dt./max(dt);                                                      % Convert to proportion to handle differences in scale
    PLT(pp) = plot(1:bins,dt,':d','Color',cmapothr(ii,:),'MarkerFaceColor',cmapothr(ii,:),pltopt{:}); hold on   
    pp      = pp + 1;
end

dplt        = binvars;
for ii = 1 : N
    dt      = dplt(ii,:);
    dt      = dt./max(dt);
    PLT(pp) = plot(1:bins,dt,':x','Color',cmapothr(ii,:),'MarkerFaceColor',cmapothr(ii,:),pltopt{:}); hold on
    pp      = pp + 1;
end
xlim(xlims); title(pltstr); ylabel(ylbl); set(gca,gcaopt{:});
ylims = ylim;
if all(0>ylims(1) & 0<ylims(2)) || all(0<ylims(1) & 0>ylims(2))
    line(xlim,[0 0],'LineStyle', '-', 'Color', 'k','LineWidth',2);
end
legend(PLT,[ strcat('group-',W8s(Iothr)) strcat('subject-',W8s(Iothr)) ]);


%--------------------------------------------------------------------------
end