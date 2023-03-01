function varargout = INhandler(Vin,nin,D)
%  Handles optional inputs. Defaults are assigned where no input is given.
%
% Input:
%   Vin     : VARARGIN
%   nin     : number of optional inputs supplied
%   D       : Default Values
%
% Output:
%   Vout    : VARARGOUT
%
% 2022 Mark C Nelson, McConnell Brain Imaging Centre, MNI, McGill
%--------------------------------------------------------------------------

N       = length(D);
Vout    = cell(1,N);

% Assign defaults where no optional input is given
for kk = 1 : N
    if nin<kk || isempty(Vin{kk}); Vout{kk} = D{kk}; end
end

% Extract supplied optional inputs
if nin > 0
    for kk = 1 : nin
        if ~isempty(Vin{kk});
            Vout{kk} = Vin{kk};
        end
    end
end

varargout=Vout;

%--------------------------------------------------------------------------
end
