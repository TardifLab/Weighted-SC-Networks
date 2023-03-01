function varargout = boxplotyy(D,N,axind,varargin)
%  Plots colored boxes with data points overlaid and 2 y axes
%
%   Note: This code is a little complicated. Please double check that it is
%   working correctly for your application
%
% Input:
%                           + Required +
%   D       : Data to plot in cell array
%   N       : Number of bins for each dataset (data binned prior to function call)
%   axind   : 1x2 cell array of indices used to map datasets to 1 of the 2
%             y axes 
%             e.g. {[1 2] [ 3 4]} -> [1 2] to LEFT & [3 4] to RIGHT
%
%                           + Optional +
%   cmap    : Nx3 colormap used to color data points & boxes
%   binlbls : 1xN cell array of labels for data bins
%   str     : character vector describing data
%
% Output:
%
% 2022 Mark C Nelson, McConnell Brain Imaging Centre, MNI, McGill
%--------------------------------------------------------------------------
T       = length(D);
ppos    = [-2559 596 2135 612];                                             % plot size & location

%% Optional inputs
nin                     = max(nargin,1) - 3;
defaults                = {colormap(hsv(N)),'','unknown data'};
[cmap,binlbls,str]      = INhandler(varargin,nin,defaults);


% Ensure Contents of D are Correct Shape
dims    = cell2mat(cellfun(@(c)  size(c), D,'UniformOutput',0));
if any(dims(:,2)>1)
    D   = cellfun(@(c)  c(:), D,'UniformOutput',0);
end


%% Boxplot
ax=myfig(str,ppos);

DPLT                = D;                                                    % Store as master data to reuse later

% Get indices for data in split which will be mapped to each axis
switch N
    case 1
        indL    = sort(axind{1}*N);
        indR    = sort(axind{2}*N);
    case 2
        indL    = sort([sort(axind{1})*N-1 sort(axind{1})*N]);
        indR    = sort([sort(axind{2})*N-1 sort(axind{2})*N]);
    case 3
        indL    = sort([sort(axind{1})*N-2 sort(axind{1})*N-1 sort(axind{1})*N]);
        indR    = sort([sort(axind{2})*N-2 sort(axind{2})*N-1 sort(axind{2})*N]);
    case 4
        indL    = sort([axind{1}*N-3 axind{1}*N-2 axind{1}*N-1 axind{1}*N]);
        indR    = sort([axind{2}*N-3 axind{2}*N-2 axind{2}*N-1 axind{2}*N]);
    otherwise
        keyboard
end
NindL           = length(indL);
NindR           = length(indR);

%% Left Y axis

% Sort out data for LEFT axis
D(1:NindL)          = DPLT(indL);                                           % Map these datasets to 1st position to use left axis
it                  = NindL;
for nn = 1 : NindR
    it              = it+1;
    indtmp          = indR(nn);
    dtmp            = nan(size(DPLT{indtmp}));                              % This is a place holder because boxplot is dumb
    D{it}           = dtmp;
end

% Make g vector to match new data order (finicky, requires contents of D to be column vector)
g                   = [];
it                  = 0;
for dd = 1 : T
    Ndatapoints         = numel(D{dd});
    gt                  = ones(Ndatapoints,1)*it; it=it+1;
    g                   = [g;gt];
%     it                  = it+1;
end

% % % % Modify g vector to control data location on x axis
% % % switch N
% % %     
% % %     case 1
% % %         keyboard
% % %     case 2
% % %         keyboard
% % %     case 3
% % %         
% % %         for vv = N-3:N:it;  g(g==vv) = vv+0.5; end
% % %         for vv = N-1:N:it;  g(g==vv) = vv-0.5; end
% % %         
% % %     case 4
% % %         keyboard
% % %     otherwise
% % %         keyboard
% % % end

% Plot left axis datasets
yyaxis left
boxplot(cell2mat(D),g,'symbol', ''); hold on;

% Figure out ideal size of data points
Nmaxdps                     = max(cellfun(@(c) length(c), D));
mrkrclr                     = cmap; 
if Nmaxdps <= 100
    mrkrsz                  = 6;
    bxfcalpha               = .4;
    mrkrclr(mrkrclr==0)     = .2;
elseif Nmaxdps <= 500
    mrkrsz                  = 4;
    bxfcalpha               = .4;
    mrkrclr(mrkrclr==0)     = .3;
elseif Nmaxdps <= 5000
    mrkrsz                  = 3;
    bxfcalpha               = .4;
    mrkrclr(mrkrclr==0)     = .4;
elseif Nmaxdps <= 10000
    mrkrsz                  = 2;
    bxfcalpha               = .5;
    mrkrclr(mrkrclr==0)     = .5;
else 
    mrkrsz                  = 1;
    bxfcalpha               = .6;
    mrkrclr(mrkrclr==0)     = .6;
end

% Add LEFT data points
spread  = 0.5;                                                              % Spread of data points from center
for nn = 1 : N                                                              % Iterate over data bins
    for dd=nn:N:NindL                                                       % Over datasets in D
        plot(rand(size(D{dd}))*spread-(spread/2)+dd,D{dd},...
        'o','linewidth',1,'MarkerFaceColor',mrkrclr(nn,:),...
        'MarkerEdgeColor',mrkrclr(nn,:),'MarkerSize',mrkrsz);
    end
end


% Color LEFT boxes
L       = cell(1,N);                                                        % Used to create legend
h       = findobj(gca,'Tag','Box');
for nn = 1 : N                                                              % Iterate over data bins
    for bb = T-nn+1 : -N : T-NindL
        L{nn}=patch(get(h(bb),'XData'),get(h(bb),'YData'),cmap(nn,:),...
                    'FaceAlpha',bxfcalpha);
        h(bb).Color     = [0 0 0];                                          % Make boxes black
    end
end
boxplot(cell2mat(D),g,'symbol', ''); hold on;                               % Replot boxes over data points (easier to see)


% Make median lines bigger and easier to see
lines               = findobj(gcf, 'type', 'line', 'Tag', 'Median');
lines(T+1:end)      = [];                                                   % delete 1st boxes (they are underneath)
for bb = T : -1 : T-NindL
    lines(bb).Color      = [0 0 0];
    lines(bb).LineWidth  = 2;
end

%% Right Y axis

% Sort out RIGHT axis
D(NindL+1:end)       = DPLT(indR);                                          % Map these datasets to 2nd position to use right axis
for nn = 1 : NindL
    indtmp      = indL(nn);
    dtmp        = nan(size(DPLT{indtmp}));                                  % This is a place holder because boxplot is dumb
    D{nn}       = dtmp;
end

% Make g vector again to match new data order
g                   = [];
it                  = 0;
for dd = 1 : T
    Ndatapoints         = length(D{dd});
    gt                  = ones(Ndatapoints,1)*it; it=it+1;
    g                   = [g;gt];
    it                  = it+1;
end


% Plot datasets for right axis
yyaxis right
boxplot(cell2mat(D),g,'symbol', ''); hold on;                               % Plot right axis elements here


% Add RIGHT data points
spread  = 0.5;                                                              % Spread of data points from center
for nn = 1 : N                                                              % Iterate over data bins
    for dd = nn+NindL : N : T                                               % Over datasets in D
        plot(rand(size(D{dd}))*spread-(spread/2)+dd,D{dd},...
        'o','linewidth',1,'MarkerFaceColor',mrkrclr(nn,:),...
        'MarkerEdgeColor',mrkrclr(nn,:),'MarkerSize',mrkrsz);
    end
end


% Color RIGHT Boxes
h       = findobj(gca,'Tag','Box');
for nn = 1 : N                                                              % Iterate over data bins
    for bb = NindR-nn+1 : -N : 1
        L{nn}=patch(get(h(bb),'XData'),get(h(bb),'YData'),cmap(nn,:),...
                    'FaceAlpha',bxfcalpha);
        h(bb).Color     = [0 0 0];                                          % Make boxes black
    end
end
boxplot(cell2mat(D),g,'symbol', ''); hold on;                               % Replot boxes


% Make median lines bigger and easier to see
lines               = findobj(gcf, 'type', 'line', 'Tag', 'Median');
lines(T+1:end)      = [];                                                   % delete 1st boxes (they are underneath)
for bb = T-NindL : -1 : 1
    lines(bb).Color      = [0 0 0];
    lines(bb).LineWidth  = 2;
end

% Final Mods
lineloc = length(axind{1})*N+.5;
linelims=[-1e10 1e10];
line([lineloc lineloc],linelims,'LineStyle', '--', 'Color', [.25 .25 .25],'LineWidth',1); % Line to separate data by axis
yyaxis left
line([lineloc lineloc],linelims,'LineStyle', '--', 'Color', [.25 .25 .25],'LineWidth',1); % Better to add line to both y axes
set(gca,'FontSize',25);
title(str);
legend([L{:}],binlbls,'FontSize',25);

%--------------------------------------------------------------------------
end
