function varargout = myfig(varargin)
% Quick open figure with custom settings
%
% Optional inputs:
%   name    : character vector for plot title
%   ppos    : 1x4 double array for plot size & position
%
% Optional outputs:
%   ax      : handle to figure
%
% 2021 Mark C Nelson, McConnell Brain Imaging Centre, MNI, McGill
%--------------------------------------------------------------------------
%% Optional inputs
nin             = nargin;
defaults        = {'',[-2558 617 658 591]};
[str,ppos]      = INhandler(varargin,nin,defaults);


%% Figure
ax=figure('Name',str,'NumberTitle','off','color','w','Position',ppos);

%% Optional output
if nargout > 0
    varargout = cell(1,nargout);
    for oo = 1 : nargout
        switch nargout
            case 1
                varargout{oo} = ax;
        end
    end
end

%--------------------------------------------------------------------------
end