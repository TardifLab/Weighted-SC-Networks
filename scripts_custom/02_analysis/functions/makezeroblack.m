function [varargout] = makezeroblack(ax,varargin)
%  Modifys the colormap for an existing plot to show 0 values in black
%
% Input:
%                           + Required +
%   ax      : plot handle
%
%                           + Optional +
%   cm      : colormap
%
% 2022 Mark C Nelson, McConnell Brain Imaging Centre, MNI, McGill
%--------------------------------------------------------------------------
% Optional input
if ~isempty(varargin)
    cm                      = varargin{1};
end

% Create colormap if none supplied
if ~exist('cm','var')
    cm                      = colormap;
end

% Function
lssd                        = 1e5;                                          % linspace sampling density (precision of 0 value marking)
cax                         = caxis;                                        % range of current colormap
ds                          = linspace(min(cax), max(cax), lssd);           % linearly spaced bins over color range
ncol                        = size(cm,1);                                   % range of colors

% Find 0 in colormap
if ~any([all(cax<=0) all(cax>=0)])                                          % if 0 lies between color limits
    dpn                     = find(ds==max(ds(ds<0)))/lssd;                 % proportion of data below 0
    zpos                    = floor(dpn * ncol);                            % find row corresponding to 0
    cm(zpos:zpos+1,:)       = 0;                                            % set that position +1 to black
elseif any(cax==0)                                                          % if 0 is a limit in caxis
    % Find 0 bound & set to 0
    if all(cax>=0)                                                          % if caxis is positive
        cm(1,:)             = 0;
    elseif all(cax<=0)                                                      % if caxis is negative
        cm(ncol,:)          = 0;
    end
end

try
    colormap(ax,cm)                                                         % activate it
catch
    disp('Unable to apply your colormap directly. Return as output & apply manually')
end


% Optional output
if nargout > 0
    varargout                           = cell(1,nargout);
    for kk = 1 : nargout
        switch kk
            case 1; varargout{kk}       = cm;
        end
    end
end

%--------------------------------------------------------------------------
end
