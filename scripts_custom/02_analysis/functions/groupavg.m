function [G] = groupavg(D,varargin)
%  Computes group-average data 
%
% Input:
%                           + Required +
%   D       : N dimensional numeric array
%
%                           + Optional +
%   dim     : Integer indicating dimension for averaging i.e. subjects 
%             (default is 3rd dimension if 3D, or smallest)
%   avg     : Character vector indicating which method to use in averaging.
%             (default is 'nz'). Options include:
%                   - 'all' : mean of all array elements (typical mean)
%                   - 'nz'  : mean of non-zero array elements only
%             
%
% Output:
%   G       : N-1 dimensional averaged numeric array
%
% 2021 Mark C Nelson, McConnell Brain Imaging Centre, MNI, McGill
%--------------------------------------------------------------------------

%% Optional inputs
nin=max(nargin,1) - 1; 
defaults            = {[] 'nz'};
[dim,avg]   = INhandler(varargin,nin,defaults);


%% Setup
if ~isnumeric(D)                                                            % Input array should be numeric
    error('Input array should be numeric')
end

if isempty(dim)                                                             % If no dimension provided for averaging
    if length(size(D)) == 3                                                 % If data is 3D
        dim         = 3;                                                    % assume 3rd dimension is subjects
    else
        dim         = min(size(D));                                         % assume smallest dimension is subjects
    end
end

if length(size(D)) < dim || size(D,dim) == 1                                % Indicated dimension should exist & be non-singleton
    error('The indicated dimension does not exist in your input array')
end

avg                 = lower(avg);

%% GROUP-AVERAGE
D(isinf(D))         = nan;
switch avg
    case 'all'
        D(isnan(D)) = 0;
        G           = mean(D,dim);
    case 'nz'
        G           = nzmean(D,dim);
    otherwise
        error('Check your avg input!')
end

%--------------------------------------------------------------------------
end