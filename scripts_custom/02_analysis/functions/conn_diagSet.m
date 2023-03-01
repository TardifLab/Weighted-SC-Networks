function D = conn_diagSet(D,value)
% Changes value on main diagoinal of square 2D matrices in a cell array.
% Matrices may be stacked in 3rd dimension i.e, n x n x s
%
% Input:
%   D       : IxJ cell array, each cell containing an NxNxS stack of
%             2D connectivity matrices
%   value   : Desired value on main diagonal
%
% Output:
%   D       : Cell array reflecting desired modification   
%
% 2021 Mark C Nelson, McConnell Brain Imaging Centre, MNI, McGill
%--------------------------------------------------------------------------

% Get dims for loop
[I,J]=size(D);

% loop over datasets
if ~isempty(value)
    for ii = 1 : I
        for jj = 1 : J

            % Get number of matrices in stack
            K = size(D{ii,jj},3);                                           % should be cell array

            for kk = 1 : K
                Dtmp                                = D{ii,jj}(:,:,kk);
                Dtmp(logical(eye(size(Dtmp,1))))    = value;
                D{ii,jj}(:,:,kk)                    = Dtmp;
            end
        end
    end
else
    disp('ERROR: Value parameter missing. Diagonals not changed...')
end
%--------------------------------------------------------------------------
end