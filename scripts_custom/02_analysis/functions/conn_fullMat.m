function D = conn_fullMat(D)
% Ensures matrices are full (assumes symmetry)
%
% Input:
%   D       : IxJ cell array, each cell containing an NxNxS stack of
%             2D connectivity matrices
%
% Output:
%   D       : Cell array reflecting desired modification   
%
% 2021 Mark C Nelson, McConnell Brain Imaging Centre, MNI, McGill
%--------------------------------------------------------------------------
[I,J]=size(D);

% Transform to full mats
for jj = 1 : J
    for ii = 1 : I
        dtmp    = D{ii,jj};
        K       = size(dtmp,3);
        for kk = 1 : K
            if ~isequal(dtmp(:,:,kk), dtmp(:,:,kk).')                       % check for symmetry first
                dtmp(:,:,kk) = dtmp(:,:,kk) + dtmp(:,:,kk).';
            end
        end
        D{ii,jj} = dtmp;
    end
end
%--------------------------------------------------------------------------
end