function labels = genericParcNames(N)
% Creates generic parcellation labels if none are provided 
%
% Input
%   N       : integer value indicating the number of parcellations
%
% Outputs
%   labels  : 1xN cell array of strings, generic parcellation names
%
% 2022 Mark C Nelson, McConnell Brain Imaging Centre, MNI, McGill
%--------------------------------------------------------------------------

labels                  = cell(1,N);
for nn = 1 : N 
    labels{nn}          = ['Parc ' num2str(nn)];
end


%--------------------------------------------------------------------------
end