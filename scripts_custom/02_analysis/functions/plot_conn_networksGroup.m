function plot_conn_networksGroup(D,W8s,pinfo,varargin)
%  Plots 3D brain networks of both nodes & edges at group level using
%  external software package.
%
% Notes:
%   1. This package can be challenging to work with.
%   2. Requires files be created for each layer of plot i.e. edges, nodes, etc
%   3. Recommend starting with the GUI before command line use
%
% Input:
%                           + Required +
%   D       : IxJ cell array, each cell containing an NxNxS stack of
%             connectivity matrices (S=subs)
%   W8s     : Cell array of strings describing edge weights
%   pinfo   : 1xJ structure with parcellation info (see conn_getParcInfo.m)
%
%                           + Optional +
%   avg     : String indicating method for group-averaging (see groupavg.m)
%             (default = non-zero)
%   str     : character vector for plot title
%
% 2022 Mark C Nelson, McConnell Brain Imaging Centre, MNI, McGill
%--------------------------------------------------------------------------
%% Optional inputs
nin                 = max(nargin,1) - 3;
defaults            = {'nz','unknown data'};
[avg,str]           = INhandler(varargin,nin,defaults);

%% Setup
[I,J]               = size(D);
fntsz               = 25;
parcs               = {pinfo.name};
ppos                = [-2559 556 1000 652];
if ~iscell(D); D={D}; end

% Files used by BrainNet Viewer
rootDir     = '/PATH/TO/YOUR/bin'; 
configdir   = [rootDir '/BrainNetViewer_20191031/Data/Configurations'];
fileSurf    = [rootDir '/BrainNetViewer_20191031/Data/SurfTemplate/BrainMesh_ICBM152_smoothed.nv']; % Alternate surfaces can be be found in package BRAPH-2-Matlab-master

% My config files (included as examples)
% INFERNO Config Files
% cnfg          = [configdir '/cnfg_inferno.mat'];                  
% cnfg_nodes    = [configdir '/cnfg-nodes_inferno_full-view.mat'];
% cnfg_edges    = [configdir '/cnfg-edges_inferno_full-view.mat'];
% cnfg_edges    = [configdir '/cnfg-edges_inferno_axial-view.mat'];

% KINDLMANN Config File
% cnfg_edges    = [configdir '/cnfg-edges_kindlmann_axial-view.mat'];
% cnfg_edges    = [configdir '/cnfg-edges_kindlmann_med-lat-dors-view.mat'];
cnfg_edges    = [configdir '/cnfg-edges_kindlmann_full-view.mat'];


%% Data Prep
for jj = 1 : J  
    strplt          = ['Group ' str ': ' parcs{jj}];
    
    % Initialize Node file for this parcellation
    N               = pinfo(jj).ncor;                                       % node count
    cstmNode        = nan(N,6);
    cstmNode(:,1:3) = pinfo(jj).coor;
    cstmNode(:,4)   = 1;                                                    % Node color (currently set with cnfg file)
    
    for ii = 1 : I
        d           = D{ii,jj};

                
        % Group Average 
        if length(size(d)) == 3                                             % if 3D, need to average across subjects
            d               = groupavg(d,3,avg);
        end
        
        %% Edges
        
        % Save Custom Node File
        cstmNode(:,5)       = 1;                                            % All nodes same size
        fileNode            = [rootDir '/BrainNetViewer_20191031/Data/Nodes/' parcs{jj} '.node'];
        save(fileNode,'cstmNode','-ascii')
        
        % Save Edge File for BrainNet Viewer
        fileEdge    = [rootDir '/BrainNetViewer_20191031/Data/Edges/' W8s{ii} '_group.edge'];
        save(fileEdge,'d','-ascii')
        
        % Plot edges with diameter as a function of weight & nodes of equivalent and small size
        FIG         = BrainNet_MapCfg(fileSurf, fileNode, fileEdge, cnfg_edges);
        figsavename = [rootDir '/BrainNetViewer_20191031/Data/Networks/' parcs{jj} '_' W8s{ii} '_group-edges.png'];
        saveas(FIG,figsavename)
        close(FIG)
                

    end
end
%--------------------------------------------------------------------------
end