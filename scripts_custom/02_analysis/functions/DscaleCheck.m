function [scale,ind] = DscaleCheck(D)
%  Rough check if all datasets are similarly distributed and can be plotted
%  on same axis
%
% Input:
%   D       : Cell array of input connectomes
%
% Output:
%   scale   : scale required to visualize these data 
%             (cell array of strings)
%   ind     : maps datasets to their scale
%
% 2022 Mark C Nelson, McConnell Brain Imaging Centre, MNI, McGill
%--------------------------------------------------------------------------

N                   = numel(D);
skewcheck           = zeros(1,N);

% Determine skewness of each dataset
for nn = 1 : N
    Dt              = D{nn};
    Dt(isnan(Dt))   = 0;
    Dt(isinf(Dt))   = 0;
    Dt              = Dt(Dt~=0);
    
    skewcheck(nn)   = abs(skewness(Dt));
end

% Create output variables
if all(skewcheck>5)
    scale   = {'log'};
    ind     = {1:N};
elseif all(skewcheck<5)
    scale   = {'linear'};
    ind     = {1:N};
else 
    scale   = {'linear' 'log'};
    ind     = cell(1,2);
    ind{1}  = find(skewcheck<5);
    ind{2}  = find(skewcheck>5);
end
    

%--------------------------------------------------------------------------
end
