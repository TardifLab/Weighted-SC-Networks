function varargout = plot_conn_surf(D,pinfo,varargin)
%  Plots cortical & subcortical node data on the surface.
%
%                           *Usage Notes*
%  1. Currently only conte69 cortical surface is supported.
%  2. Currently does not support Cerebellum!
%  3. Uses BrainSpace & ENIGMA toolboxes
%
% INPUT:
%                           + Required +
%   D       : NxI double of node data to plot (N=Nodes, I=Datasets)
%             (Also supports 1xI cell arrays with Nx1 size cells)
%   pinfo   : 1xJ structure with parcellation info (see conn_getParcInfo.m)
%
%                           + Optional +
%   Type    : String indicating desired plot type 
%             {'cortex' 'subcortex' 'both'}; default = 'cortex'
%   lbls    : Cell array of plot labels
%   str     : String for plot title
%
%
% OUTPUT:
%                           + Optional +
%   Hcx     : 1xFigures cell array of Handles to cortical surface plots. 
%            (Includes figure, axes, colorbar handles) (size depends on I)
%   Hsx     : 1xFigures cell array of Handles to subcortical surface plots.
%            (Includes figure, axes, colorbar handles) (size depends on I)
%
% 2022 Mark C Nelson, McConnell Brain Imaging Centre, MNI, McGill
%--------------------------------------------------------------------------

if iscell(D); D     = cell2mat(D); end
errstr              = '* * Surface Plot Fail * *     Check input: ';
[N,I]               = size(D);

% Reshape data if necessary (assumes Ndatasets < Nnodes)
if I > N
    D               = D';
    [N,I]           = size(D);
end

%% Optional inputs
nin                 = max(nargin,1) - 2;
defaults            = {'cortex','',''};
[Type,lbls,str]     = INhandler(varargin,nin,defaults);


%% Setup
Type                        = lower(Type);
parc                        = pinfo.name;
parc                        = strrep(parc,'-',' ');
Ncx                         = pinfo.ncor;
Nsx                         = pinfo.nsub;
Nall                        = Ncx + Nsx;

% Prep to split data if necessary
if strncmpi(Type,'both',4) || strncmpi(Type,'subcortex',4)
    Isx_lbls    = cellstrfind(pinfo.clabels,'subcortex');                   % Location of subcortical indices in pinfo structure
    Isx_nodes   = pinfo.cis{Isx_lbls};                                      % Indices for subcortical nodes in D
    Icx_nodes   = setdiff([1:N]',Isx_nodes);                                % Indices for cortical nodes in D
end

% Double check inputs and prep data according to the desired plot
switch Type
    
    case 'both'
        
        if isempty(Ncx) || isempty(Nsx); error([errstr ' Type & pinfo ']); end
        if Nall ~= N;                    error([errstr ' pinfo & D ']);    end
        
        % Split cortical & Subcortical data
        Dsx     = D(Isx_nodes,:);
        Dcx     = D(Icx_nodes,:);
%         SUBCTX  = true;
        
    case 'cortex'
        
        if isempty(Ncx); error([errstr ' Type & pinfo ']); end
        if Ncx ~= N;     error([errstr ' pinfo & D ']);    end
        
        Dcx     = D;
        
    case 'subcortex'
        
        if isempty(Nsx); error([errstr ' Type & pinfo ']); end
        if Nsx ~= N;     error([errstr ' pinfo & D ']);    end
        
        Dsx     = D;
%         SUBCTX  = true;
        
    otherwise
        
        error([errstr ' Type '])
end


% Reorder subcortical data to match ENIGMA setup (ventricles removed)
if exist('Dsx','var');
    enigmasxlbls    = { 'Left-Accumbens-area'  'Left-Amygdala'  'Left-Caudate'  'Left-Hippocampus'  'Left-Pallidum'  'Left-Putamen'  'Left-Thalamus-Proper' ... 
                       'Right-Accumbens-area' 'Right-Amygdala' 'Right-Caudate' 'Right-Hippocampus' 'Right-Pallidum' 'Right-Putamen' 'Right-Thalamus-Proper'};
    mysxlbls        = pinfo.labels(Isx_nodes)';
    Isxr            = cell2mat(cellstrfind(enigmasxlbls,mysxlbls));         % Indices to reorder subcortical labels
    Dsxr            = zeros(length(enigmasxlbls),I);
    Dsxr(Isxr,:)    = Dsx;                                                  % Reorder subcotical nodes
    Dsx             = Dsxr; clear Dsxr
end


%% Plots

% Plot cortical data on surface (Brainspace)
if exist('Dcx','var')
    
    conte69is               = pinfo.conte69i;                               % indices for mapping node data to cortical surface
    [srf_lh,srf_rh]         = load_conte69();                               % load surface from brainspace
    SURF                    = {srf_lh,srf_rh};
    Dcxf                    = cell(1,I);
    
    for ii = 1 : I
        Dcxf{ii}            = parcel2full(Dcx(:,ii),conte69is);             % Map parcellated cortical data to the surface
    end
    
    % If more than 4 datasets plot them all separately (only 4 can be plot together)
    if I > 4
        Hcx     = cell(1,I);
        for ii = 1 : I
                Dt          = Dcxf{ii};
                opt         = {'labeltext',lbls{ii}};
                Hcxt        = plot_hemispheres(Dt,SURF,opt{:});
                Hcx{ii}     = Hcxt;
                Hcxt.figure.NumberTitle='off'; Hcxt.figure.Name=str;
        end
    % Else plot them together
    else
        Dt          = cell2mat(Dcxf);
        opt         = {'labeltext',lbls};
        Hcx{1}      = plot_hemispheres(Dt,SURF,opt{:});
        Hcx{1}.figure.NumberTitle='off'; Hcx{1}.figure.Name=str;
    end
end


% Plot subcortical data on surface (ENIGMA)
if exist('Dsx','var')
    
    Hsx                 = cell(1,I);
    for ii = 1 : I
        Figsx           = myfig(str,[1442 210 747 587]);
        opt             = {'label_text',lbls{ii},'ventricles','False'};
        [AXsx,CBsx]     = plot_subcortical(Dsx(:,ii),opt{:});

        % Package handles
        Hsxt            = [];
        Hsxt.figure     = Figsx;
        Hsxt.axes       = AXsx;
        Hsxt.cb         = CBsx;
        Hsx{ii}         = Hsxt;
    end
end


%% Optional output
nout = nargout;
if nout > 0
    varargout = cell(1,nout);
    for oo = 1 : nout
        switch oo
            case 1; varargout{oo} = Hcx;
            case 2; varargout{oo} = Hsx;
        end
    end
end

%--------------------------------------------------------------------------
end
