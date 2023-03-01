function [varargout] = red2bluecmap(ax,varargin)
%  Joins hot & cold colormaps approximatley at the data value in 
%  "d" indicated by "val"
%
% Input:
%                           + Required +
%   ax      : plot handle
%
%                           + Optional +
%   val     : target value (default=0)
%
% Output:
%                           + Optional +
%   cm      : colormap handle
%
% 2022 Mark C Nelson, McConnell Brain Imaging Centre, MNI, McGill
%--------------------------------------------------------------------------

%% Optional input
if ~isempty(varargin)
    val                     = varargin{1};
else
    val                     = 0;                                            % assume 0 is target if no value is given
end

%% Setup
lssd                        = 1e5;                                          % linspace sampling density (precision of target value marking)
cax                         = caxis;                                        % range of current colormap
ds                          = linspace(min(cax), max(cax), lssd);           % linearly spaced bins over range
ncol                        = 64;                                           % assume standard size colormap
cm                          = zeros(ncol,3);                                % init blank colormap

%% Build colormap
if ~any([all(cax<=val) all(cax>=val)])                                        % if val is contained within your color axis
    % Find target value position in data
    dpn                     = find(ds==max(ds(ds<val)))/lssd;               % proportion of data below target value
    zpos                    = ceil(dpn * ncol);                             % find approximate row corresponding to target value in colormap
    cm(zpos,:)              = 1;                                            % Set middle row to white

     % Create color gradient  
     maxseg                 = max([zpos, ncol-zpos]);                       % number of rows in cmap for segment of data containing maxval
     cgrad                  = linspace(.98,0,maxseg);

    % Build BLUE gradient for portion of data below target
    it                      = zpos-1:-1:1;                                  % span of rows for this color gradient (reverse order)
    itn                     = length(it);
    cm(it,1)                = cgrad(1:itn);
    cm(it,2)                = cgrad(1:itn);
    cm(it,3)                = 1;

    % Build RED gradient to portion of data above target
    it                      = zpos+1:ncol;
    itn                     = length(it);
    cm(it,1)                = 1;
    cm(it,2)                = cgrad(1:itn);
    cm(it,3)                = cgrad(1:itn);
    
elseif all(cax<=0)                                                           % If all color limits are negative
    cr      = linspace(0,min(cax),100);                                     % linearly spaced range from 0 to min value in caxis
    cpn     = find(cr==max(cr(cr<max(cax))));                               % proportion of range [min 0] corresponding to this caxis
    cgtmp   = linspace(.98,0,100);                                          % linearly spaced range over color gradient
    cgmin   = cgtmp(cpn);                                                   % minimum value for color gradient

    % Create blue gradient from 0
    cgrad                   = linspace(0,cgmin,ncol);
    cm(:,1)                 = cgrad;
    cm(:,2)                 = cgrad;
    cm(:,3)                 = 1;
    
elseif all(cax>=0)                                                           % If all color limits are positive
    % Create red gradient from 0
    cr      = linspace(0,max(cax),100);                                     % linearly spaced range from 0 to max value in caxis
    cpn     = find(cr==min(cr(cr>min(cax))));                               % proportion of range [0 max] corresponding to this caxis
    cgtmp   = linspace(.98,0,100);                                          % linearly spaced range over color gradient
    cgmin   = cgtmp(cpn);                                                   % minimum value for color gradient
    
    cgrad                   = linspace(cgmin,0,ncol);
    cm(:,1)                 = 1;
    cm(:,2)                 = cgrad;
    cm(:,3)                 = cgrad;
else
    error('Unable to build your colormap. Check your ax input')
end

% Activate colormap
try
    colormap(ax,cm)
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
