function varargout = unicmap(axin,cmname)
% Changes colormap of all axes in axin
%
% Input:
%   axin    : 1xN array of axis handles
%   cmname  : string indicating the desired colormap (default=inferno)
%             (must be a matlab built in or map supported by smartcmaps.m
%
% Output:
%                           + Optional +
%   axcell  : input cell reflecting changes
%
% 2022 Mark C Nelson, McConnell Brain Imaging Centre, MNI, McGill
%--------------------------------------------------------------------------

cmname      = lower(cmname);

% Get colormap if not a matlab built in
cstmcmaps   = {'inferno' 'kindlmann' 'plasma' 'viridis' 'smoothcoolwarm' 'bentcoolwarm'};
if ismember(cmname,cstmcmaps)
    cmap    = smartcmaps(cmname);
else
    cmap    = cmname;
end

if iscell(axin)
    n       = length(axin);
    
    for ii = 1 : n
        colormap(axin{ii},cmap);
    end
    
else
    % Quick solution to handle multiple cases
    axin                    = axin(:);
    n                       = length(axin);
    
    for ii = 1 : n
        colormap(axin(ii),cmap);
    end
end



%% Optional output
if nargout > 0
    varargout = cell(1,nargout);
    for oo = 1 : nargout
        switch oo
            case 1
                varargout{1} = axin;
        end
    end
end

%--------------------------------------------------------------------------
end