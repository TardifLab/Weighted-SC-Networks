function [C,varargout] = conn_corrmat(D1,D2,varargin)
% Vectorizes & computes correlations for 2 input matrices.
% All zero valued entries are excluded.
%
%   Assumes:
%       - matrices are same size
%       - matrices represent undirected 2D connectivity (symmetric)
%
%
% INPUT:
%                           + Required +
%   D1      : NxN undirected connectivity matrix
%   D2      : NxN undirected connectivity matrix
%
%                           + Optional +
%   Type    : Character vector indicating type of correlation
%             {'Pearson' 'Spearman'} (default=Pearson)
%
%
% OUTPUT:
%                           + Required +
%   C       : correlation coefficient
%
%                           + Optional +
%   V       : Lx2 vectorized elements of D1 & D2
%
% 2021 Mark C Nelson, McConnell Brain Imaging Centre, MNI, McGill
%--------------------------------------------------------------------------

%% Optional inputs
nin             = max(nargin,1) - 2;
defaults        = {'Pearson'};
Type            = INhandler(varargin,nin,defaults);


%% Function
Copt            = {'Type',Type,'Rows','complete'};                          % Correlation call flags

N               = length(D1);                                               % Node count
mask            = logical(tril(ones(N),-1));                                % Lower triangle mask LESS main diagonal

D1(isinf(D1))   = 0;                                                        % Ensure no infs
D2(isinf(D2))   = 0;

D1(D1==0)       = nan;                                                      % Zeros to nans for exclusion
D2(D2==0)       = nan;

V1              = D1(mask);                                                 % vectorize both mats
V2              = D2(mask);

C               = corr(V1,V2,Copt{:});                                      % run correlation

%% Optional output
nout = nargout - 1;
if nout > 0
    varargout = cell(1,nout);
    for oo = 1 : nout
        switch oo
            case 1; varargout{oo} = [V1,V2];
        end
    end
end

%--------------------------------------------------------------------------
end