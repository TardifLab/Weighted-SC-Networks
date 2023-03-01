function [D] = commit_scaler(D,LoS)
%  Scales COMMIT weights by mean edge length 
%
% Input:
%   D       : NxNxS stack of subject-level connectomes with edges weighted
%             by: SUM(Xk * Lk) over all streamlines connecting nodes i & j;
%             Where Xk is the COMMIT weight & Lk is the length of 
%             streamline k.
%   LoS     : NxNxS stack of subject-level connectomes with edges weighted
%             by the mean edge length over streamlines.
%
% Output:
%   D       : NxNxS stack of subject-level connectomes with edges weighted
%             by: SUM(Xk * Lk)/Lmean_ij 
%
% 2021 Mark C Nelson, McConnell Brain Imaging Centre, MNI, McGill
%--------------------------------------------------------------------------

[~,J]                       = size(D);                                      % data dimensions: weights x parcellations

if ~iscell(D);   D          = {D};   end
if ~iscell(LoS); LoS        = {LoS}; end


for jj = 1 : J
    dcommit                 = D{jj};
    dlos                    = LoS{jj};
    dlos(dlos==0)           = nan;
    dcommit(dcommit==0)     = nan;
    dt                      = dcommit ./ dlos;                              % Normalized by LoS
    dt(isnan(dt))           = 0;
    D{jj}                   = dt;
end

end

