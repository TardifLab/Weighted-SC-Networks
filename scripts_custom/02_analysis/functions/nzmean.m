function [Am] = nzmean(A,dim)
% Computes mean of nonzero-elements of N dimensional array
%
% Input:
%   A     : N dimensional numeric array
%   dim   : dimension for computation
%
% Output:
%   A     : N-1 dimensional numeric array (mean along given dim)
%
% 2021 Mark C Nelson, McConnell Brain Imaging Centre, MNI, McGill
%--------------------------------------------------------------------------

%% Ensure data is in correct format
if length(size(A)) < dim || size(A,dim) == 1                                % Indicated dimension should exist & be non-singleton
    error('The indicated dimension does not exist in your input array')
end
if ~isnumeric(A)                                                            % Input array should be numeric
    error('Input array should be numeric')
end
A(isinf(A))         = nan;

%% MEAN OF NON-ZERO ARRAY ELEMENTS
% Compute mean
A(A==0)             = nan;
Am                  = nanmean(A,dim);
Am(isnan(Am))       = 0;

%--------------------------------------------------------------------------
end