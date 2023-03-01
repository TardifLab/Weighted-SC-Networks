function [D,W8s] = conn_fc_r2z(D,W8s)
%  Fisher R-Z transform functional connectivity data
%
% Input:
%   D       : IxJ cell array, each cell containing an NxNxS stack of
%             2D connectivity matrices
%   W8s     : Cell array of strings describing edge weights
%
% Output:
%   D       : Same cell array except now with R-Z transformed  FC data
%   W8s     : modified to reflect change 
%
% 2021 Mark C Nelson, McConnell Brain Imaging Centre, MNI, McGill
%--------------------------------------------------------------------------

%Setup
[I,J]   = size(D);
ii      = cellstrfind(W8s, 'fc');

% Main Function
if ~isempty(ii) && ismember(ii, 1:I)
    disp(['Fisher r-z at index: ' num2str(ii')])
    for kk = 1 : numel(ii)
        if isempty(regexp(W8s{ii(kk)}, 'fcz','ignorecase'))                 % make sure it hasn't already been transformed
            for jj = 1 : J
                dims            = size(D{ii(kk),jj});
                Z               = nan(dims);
                for dd = 1 : dims(3)
                    Rcorr       = D{ii(kk),jj}(:,:,dd);
                    Z(:,:,dd)   = 0.5 * log( (1+Rcorr) ./ (1-Rcorr) );      % Fisher r-z
                end
                D{ii(kk),jj}    =   Z;
            end
            W8s{ii(kk)} = regexprep(W8s{ii(kk)},'fc','FCz','ignorecase');
        else
            disp(['ERROR: data already r-z transformed at index: ' num2str(ii(kk))])
            continue
        end
    end
else
    disp('** ERROR: Unable to Fisher r-z transform FC data... **')
end
%--------------------------------------------------------------------------
end