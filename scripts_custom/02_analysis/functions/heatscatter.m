function [varargout] = heatscatter(X, Y, numbins, markersize, marker, plot_colorbar, plot_lsf, xlab, ylab, title)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% heatscatter(X, Y, outpath, outname, numbins, markersize, marker, plot_colorbar, plot_lsf, xlab, ylab, title)
% mandatory:
%            X                  [x,1] array containing variable X
%            Y                  [y,1] array containing variable Y
% optional:
%            numbins            [double], default 50
%                                number if bins used for the
%                                heat3-calculation, thus the coloring
%            markersize         [double], default 10
%                                size of the marker used in the scatter-plot
%            marker             [char], default 'o'
%                                type of the marker used in the scatter-plot
%            plot_colorbar      [double], boolean 0/1, default 1
%                                set whether the colorbar should be plotted
%                                or not
%            plot_lsf           [double], boolean 0/1, default 1
%                                set whether the least-square-fit line
%                                should be plotted or not (together with
%                                the correlation/p-value of the data
%            xlab               [char], default ''
%                                lable for the x-axis
%            ylab               [char], default ''
%                                lable for the y-axis
%            title              [char], default ''
%                                title of the figure
%
% Output:
%                           + Optional +
%   SH     : handle to scatter plot
%   LH     : handle to lsline
%
%
% Modification history:
%
% Mark C Nelson, 2020, MNI & McGill University
%   1. remove save file related input variables & comment out functionality
%   2. comment out force double from char inputs (caused errors)
%   3. Remove display of r & p values with text (doesn't work with subplot)
%   4. Add optional output of scatter plot handle
%--------------------------------------------------------------------------


    %%%% mandatory
    if ~exist('X','var') || isempty(X)
        error('Param X is mandatory! --> EXIT!');
    end
    if ~exist('Y','var') || isempty(Y)
        error('Param Y is mandatory! --> EXIT!');
    end
%     if ~exist('outpath','var')
%         error('Param outpath is mandatory! --> EXIT!');
%     end
%     if ~exist('outname','var') || isempty(outname)
%         error('Param outname is mandatory! --> EXIT!');
%     end
    
    %%%% optional
    if ~exist('numbins','var') || isempty(numbins)
        numbins = 50;
    else
%         % force number, not char input
%         numbins = double(numbins);
    end
    if ~exist('markersize','var') || isempty(markersize)
        markersize = 10;
    else
        % force number, not char input
        markersize = str2double(markersize);
    end
    if ~exist('marker','var') || isempty(marker)
        marker = 'o';
    end
    if ~exist('plot_colorbar','var') || isempty(plot_colorbar)
        plot_colorbar = 1;
    end
    if ~exist('plot_lsf','var') || isempty(plot_lsf)
        plot_lsf = 1;
    end
    if ~exist('xlab','var') || isempty(xlab)
        xlab = '';
    end
    if ~exist('ylab','var') || isempty(ylab)
        ylab = '';
    end
    if ~exist('title','var') || isempty(title)
        title = '';
    end
    
    [values, centers] = hist3([X Y], [numbins numbins]);

    centers_X = centers{1,1};
    centers_Y = centers{1,2};

    binsize_X = abs(centers_X(2) - centers_X(1)) / 2;
    binsize_Y = abs(centers_Y(2) - centers_Y(1)) / 2;
    bins_X = zeros(numbins, 2);
    bins_Y = zeros(numbins, 2);

    for i = 1:numbins
        bins_X(i, 1) = centers_X(i) - binsize_X;
        bins_X(i, 2) = centers_X(i) + binsize_X;
        bins_Y(i, 1) = centers_Y(i) - binsize_Y;
        bins_Y(i, 2) = centers_Y(i) + binsize_Y;
    end

    scatter_COL = zeros(length(X), 1);

    onepercent = round(length(X) / 100);
    
    
    for i = 1:length(X)

        last_lower_X = NaN;
        last_higher_X = NaN;
        id_X = NaN;

        c_X = X(i);
        last_lower_X = find(c_X >= bins_X(:,1));
        if (~isempty(last_lower_X))
            last_lower_X = last_lower_X(end);
        else
            last_higher_X = find(c_X <= bins_X(:,2));
            if (~isempty(last_higher_X))
                last_higher_X = last_higher_X(1);
            end
        end
        if (~isnan(last_lower_X))
            id_X = last_lower_X;
        else
            if (~isnan(last_higher_X))
                id_X = last_higher_X;
            end
        end

        last_lower_Y = NaN;
        last_higher_Y = NaN;
        id_Y = NaN;

        c_Y = Y(i);
        last_lower_Y = find(c_Y >= bins_Y(:,1));
        if (~isempty(last_lower_Y))
            last_lower_Y = last_lower_Y(end);
        else
            last_higher_Y = find(c_Y <= bins_Y(:,2));
            if (~isempty(last_higher_Y))
                last_higher_Y = last_higher_Y(1);
            end
        end
        if (~isnan(last_lower_Y))
            id_Y = last_lower_Y;
        else
            if (~isnan(last_higher_Y))
                id_Y = last_higher_Y;
            end
        end

        scatter_COL(i) = values(id_X, id_Y);
    
    end
    
    SH=scatter(X, Y, markersize, scatter_COL, marker);
    
    if (plot_colorbar)
        colorbar;
    end
    
   if (plot_lsf) 													% display corr (ON) & p value (OFF) loc doesn't work with subplot
       LH = lsline;
       set(LH, 'Color', 'k');
   end
    
    if (~isempty(xlab))
        xlabel(xlab);
    end
    if (~isempty(ylab))
        ylabel(ylab);
    end
    if (~isempty(title))
        title(title);
    end
    
    
    %% Optional output
nout = nargout;
if nout > 0
    varargout = cell(1,nout);
    for oo = 1 : nout
        switch oo
            case 1; varargout{oo} = SH;
            case 2; varargout{oo} = LH;    
        end
    end
end
    
       
end