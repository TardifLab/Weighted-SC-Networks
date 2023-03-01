function varargout = plot_violin(Dm,varargin)
%  Creates violin plots using external software 
%
% Input:
%                           + Required +
%   Dm      : Main data input; SxL cell array; S=samples, L=labels. 
%             Contents of Dm will be plotted such that S data points are 
%             grouped under L columns (Cells should contain numeric data) 
%
%
%                           + Optional +
%   Dlbl    : 1xL cell array of character vectors describing columns (xaxis)                
%   Do      : Optional additional data input; 1xL cell array; Contents of
%             Do will be plotted in each column using a colored bar.
%   Ylbl    : Character vector describing data being plotted (yaxis)
%   infostr : 1x2 cell array of character vectors for plot & fig titles
%   parc    : Character vector describing parcellation scheme
%   ppos    : 1x4 numeric array indicating size & position of plot
%   cmap    : Lx3 colormap
%
% Outputs
%                           + Optional +
%   fig     : figure handle
%   V       : plot handle
%
% 2022 Mark C Nelson, McConnell Brain Imaging Centre, MNI, McGill
%--------------------------------------------------------------------------

%% Handle inputs
nin                                     = max(nargin,1) - 1;
defaults                                = {'',[],'unknown data',{'' ''},'',[1442 341 594 456],[]};
[Dlbl,Do,Ylbl,infostr,parc,ppos,cmap]   = INhandler(varargin,nin,defaults);


%% Setup
rot8                                    = 20;                               % used to rotate x-axis labels
fntsz                                   = 25;
if length(infostr)==2 & iscell(infostr)
    [pltstr,figstr]                     = deal(infostr{1}, infostr{2});
else
    figstr                              = infostr;
    pltstr                              = infostr;
end


%% Plot violin
fig                                     = myfig([figstr ': ' parc],ppos);

% Violinplot.m is picky about the shape of the input data
% [a,b]                                   = size(Dm);
% if a > b 
%     Dm                                  = Dm';end
% if ~isempty(Do)
%     [a,b]                               = size(Do);
%     if a > b 
%         Do                              = Do';
%     end
% end
L                                       = size(Dm,2);


% Plot
try
    if iscell(Dm)
        T                               = cell2table(Dm);
        V                               = violinplot(T,'ShowMean');hold on
    else
        T                               = array2table(Dm);
        V                               = violinplot(T,'ShowMean');hold on
    end
catch
    Dm                                  = cellfun(@(c) c',Dm,'UniformOutput',0); % Try just flipping the vectors 1st
    try
        T                               = cell2table(Dm);
        V                               = violinplot(T,'ShowMean');hold on
    catch
            disp('Violin fail. See plot_violin.m line 50')
            keyboard
    end
end

% Modify violins
for vv = 1:L
    V(vv).ViolinAlpha                   = .5;                               % make violin darker
    V(vv).ScatterPlot.SizeData          = 50;                              % make data points larger
    V(vv).ScatterPlot.MarkerFaceAlpha   = .3;                               % And darker
    V(vv).ViolinPlot.LineWidth          = 3;
    
    V(vv).MedianPlot.Marker             = 'd';
    V(vv).MedianPlot.SizeData           = 100;
    V(vv).ShowNotches                   = 1;
    V(vv).ShowMean                      = 1;                                % Make mean more distinct
    V(vv).MeanPlot.LineWidth            = 8;
    V(vv).MeanPlot.Color                = [1 1 1];
    
    % Make mean bar wider
    mw                                  = V(vv).MeanPlot.XData;  
    V(vv).MeanPlot.XData                = [mw(1)-.1 mw(2)+.1];

    % Add colored bar for Do if given
    if ~isempty(Do)
        if iscell(Do)
            V(vv).MeanPlot.YData        = [Do{vv} Do{vv}];                  % Manually change median line to your input
        else
            V(vv).MeanPlot.YData        = [Do(vv) Do(vv)];
        end
    end
    
    if ~isempty(cmap)
        V(vv).ViolinColor               = cmap(vv,:);
        V(vv).ViolinPlot.EdgeColor      = cmap(vv,:);
    end
end

ylabel(Ylbl); title(pltstr); xlim([.5 L+.5]);
set(gca,'XTickLabel',Dlbl,'XTickLabelRotation',rot8,'FontSize',fntsz);

% Add 0 line
ylims = ylim;
if ylims(1)<0 && ylims(2)>0
    line(xlim,[0 0],'LineStyle', '--', 'Color', 'r','LineWidth',3);
end

%% Optional output
if nargout > 0
    varargout = cell(1,nargout);
    for oo = 1 : nargout
        switch oo
            case 1
                varargout{1}    = fig;
            case 2
                varargout{2}    = V;
        end
    end
end


%--------------------------------------------------------------------------
end
