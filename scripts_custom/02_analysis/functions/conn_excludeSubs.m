function [D] = conn_excludeSubs(D,W8s,HCi,varargin)
%  Cuts out the indicated subject data from STRUCTURAL connectomes only.
%  Removed datasets correspond to 2D matrix slices in the overall n x n x
%  subjects 3D shape. 
%
%  NOTE: Recommend excluding before consistency-based thresholding, as it 
%        relies on cross-subject variance.
%
% Input:
%                           + Required +
%   D    : IxJ cell array, each cell containing an NxNxS stack of
%          2D connectivity matrices
%   W8s  : Cell array of edge weights (to differentiate SC & FC data)
%   HCi  : Indices of subjects to be excluded
%
%                           + Optional +
%   FCrm : []/1 option to exclude subs in FC data too (default=exclude)
%
% Output:
%   D   : Dataset minus exclusions
%
% 2021 Mark C Nelson, McConnell Brain Imaging Centre, MNI, McGill
%--------------------------------------------------------------------------

%% Optional inputs
nin=max(nargin,1)-3; 
defaults        = {1};
[FCrm]          = INhandler(varargin,nin,defaults);

%% Main Function
[I,J]           = size(D);
Sn              = size(D{1,1},3);

if ~isempty(HCi) && round(HCi)==HCi && HCi>0 && HCi<= Sn                    % is not empty & is valid subject index
    disp(['EXCLUDING the following subjects: ' num2str(HCi)])
    if ~isempty(FCrm); disp('EXCLUDING from FC data too!'); end
    for jj = 1 : J
        for ii = 1 : I
            if isempty(regexp(W8s{ii},'fc','ignorecase')) || ~isempty(FCrm)
                Dtmp              = D{ii,jj};
                Dtmp(:,:,HCi)     = [];                                     % remove subs (mats are shape n x n x subs)
                D{ii,jj}          = Dtmp;
            end
        end
    end
else
    disp('HOORAY, NO SUBS to exclude!')
end
%--------------------------------------------------------------------------
end