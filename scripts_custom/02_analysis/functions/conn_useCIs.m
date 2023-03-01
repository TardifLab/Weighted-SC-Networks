function [it,cil,cirb,cirm] = conn_useCIs(D,pinfo)
%  Returns indices which can be used to reorder connectivity matrices
%  according to input community assignments. 
%  Nodes additionally ordered by strength within each community. 
%
% Input:
%   D       : N x N undirected connectivity matrix
%   pinfo   : Structure with parcellation info (see conn_getParcInfo.m)
%
% Output:
%   it      : 1 x N vector of node indices reordered by input CIs
%   cil     : 1 x COMMUNITIES cell array of character labels
%   cirb    : 1 x COMMUNITIES double indicating the boundary between
%             communities
%   cirm    : 1 x COMMUNITIES double indicating tick locations for adding
%             labels to axes (the midpoint of each community)
%
% 2021 Mark C Nelson, McConnell Brain Imaging Centre, MNI, McGill
%--------------------------------------------------------------------------

cil             = pinfo.clabels;                                            % community labels by network
civ             = pinfo.cirois';                                            % nodes labeled by community assignments
cic             = pinfo.cis;                                                % indices for community assignments by network

cir             = cellfun(@(c)  numel(c),cic);                              % number of nodes in each community


% Get border of each communitie
cirb            = cir;                                                      
for cc = 2 : numel(cir)
    cirb(cc)    = cirb(cc) + cirb(cc-1);                                    
end

% Get midpoint of each community
cirm            = cirb;
for cc = 1 : numel(cir)
    cirm(cc)    = ceil(cirm(cc) - (cir(cc)/2));                             
end

strength        = nansum(abs(D))';                                          % Note absolute value; doesn't matter if norm [0 1]
[~,it]          = sortrows([civ strength],[1 -2]);


%--------------------------------------------------------------------------
end