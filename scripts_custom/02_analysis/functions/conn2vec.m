function [V,varargout] = conn2vec(Cm)
%  Accepts a cell array of 2D connectivity matrices & vectorizes them
%  such that only non-zero elements common to all matrices are retained.
%
%  Input matrices should be symmetric, and uniform in size.
%
% Input:
%   Cm      : N dimensional cell array, each cell containing a 2D square,
%             symmetric, nxn connectivity matrix.
%
% Outputs:
%   V       : Cell array of size(Cm), each cell containing the vectorized
%             form of the corresponding input matrix (non-zero entries only)
%
%                           + Optional +
%   idx     : indices of non-zero values common to all input matrices
%   masko   : nxn binary matrix with idx=1 in lower triangle
%             
%
% 2021 Mark C Nelson, McConnell Brain Imaging Centre, MNI, McGill
%--------------------------------------------------------------------------
opt = 'UniformOutput';

%% Check inputs
if ~iscell(Cm)                                                              % enforce cell class
    Cm = {Cm};
end
if any(cellfun(@(c) ~isnumeric(c),Cm))                                      % input matrices should be numeric (wont catch nans)
    error('Input matrices should be numeric')
end
if length(unique(cellfun(@(c) numel(c),Cm))) ~= 1                           % ensure uniform size
    error('Input matrices should be uniform in size')
end
if any(cellfun(@(c) any(any(c~=c')),Cm))                                    % ensure symmetry
    warning('Check symmetry of input matrices! Only lower tri is vectorized...')
end

%% Vectorize data
Do              = Cm{1};                                                    % use 1st cell as starting point
Do(isnan(Do))   = 0;
masko           = tril(Do,-1)~=0;                                           % lower tri non-zero elemenst of 1st dataset

% Iteratively modify initial mask to exclude zeros from other datasets
for ii = 2 : numel(Cm)
    Dt                          = Cm{ii};
    maski                       = tril(Dt,-1)~=0;
    masko(masko==1 & maski==0)  = 0;
end

% Vectorize all matrices using these indices
V                               = cell(size(Cm));
idx                             = find(masko);
for ii = 1:numel(Cm)
    V{ii}                       = Cm{ii}(idx);
end

%% Optional outputs
nout                            = max(nargout,1) - 1;
if nout > 0
    optout                      = {idx, masko};
    varargout                   = cell(1,nout);
    for kk = 1:nout
        varargout{kk}           = optout{kk};
    end
end

%--------------------------------------------------------------------------
end