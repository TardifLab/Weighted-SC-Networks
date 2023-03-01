function [D,W8s] = conn_sc_qt12r1(D,W8s)
%  converts qT1 (longitudinal relaxation rate) to R1
%
% Input:
%   D       : IxJ cell array, each cell containing an NxNxS stack of
%             matrices
%   W8s     : Cell array of strings describing edge weights
%
% Output:
%   D       : modified data
%   W8s     : modified data labels 
%
% 2021 Mark C Nelson, McConnell Brain Imaging Centre, MNI, McGill
%--------------------------------------------------------------------------

%% setup
[I,J]   = size(D);
ii      = cellstrfind(W8s, 't1');                                           % get indices for T1 data

% some checks
if iscolumn(ii); ii = ii'; end                                              % ensure row vector to avoid failure
if isempty(ii)
    error('* ERROR: Unable to find indices for t1 using your labels')
end
if any(~ismember(ii, 1:I))
    error('* ERROR: mismatch between dimensions of yoru labels & input data')
end

%% Compute inverse
disp(['Converting T1 to R1 at indices: ' num2str(ii)])
for kk = 1 : numel(ii)
    for jj = 1 : J
            Dtmp                = D{ii(kk),jj};
            Dtmp                = 1./Dtmp;
            Dtmp(isinf(Dtmp))   = 0;
            D{ii(kk),jj}        = Dtmp; 
    end
    W8s{ii(kk)}                 = regexprep(W8s{ii(kk)},'t1','R1','ignorecase');
end
%--------------------------------------------------------------------------
end
