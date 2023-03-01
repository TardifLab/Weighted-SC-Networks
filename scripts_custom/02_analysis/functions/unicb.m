function varargout = unicb(axin)
% Applies uniform color limits to indicated axes
%
% Input:
%   axin    : 1xN array of axis handles
%
% Output:
%                           + Optional +
%   axcell  : input cell reflecting changes
%
% 2022 Mark C Nelson, McConnell Brain Imaging Centre, MNI, McGill
%--------------------------------------------------------------------------

if iscell(axin)
    n                       = length(axin);
    clims                   = [min(cellfun(@(c) c.CLim(:,1), axin)) ... 
                               max(cellfun(@(c) c.CLim(:,2), axin))];         % Assumes color axis not inverted
    for ii = 1 : n
        axin{ii}.CLim       = clims;
    end
    
else
    % Quick solution to handle multiple cases
    axin                    = axin(:);
    n                       = length(axin);
    
    clims                   = [min(arrayfun(@(c) c.CLim(:,1), axin)) ... 
                               max(arrayfun(@(c) c.CLim(:,2), axin))];
                           
    for ii = 1 : n
        axin(ii).CLim       = clims;
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