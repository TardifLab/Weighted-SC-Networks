function varargout = connplot(D,varargin)
% Quick plot of 2D square matrix datasets
%
% Input
%                           + Required +
%   D       : Data to plot
%
%                           + Optional +
%   str     : Info for plot title
%   pinfo   : Structure with parcellation info
%   ppos    : plot size & position
%   cmmod   : modify colormap to make 0 values black (1/0, default = 0)
%
% Outputs
%                           + Optional +
%   fig     : handle to figure
%   plt     : handle to plot
%   cb      : handle to colorbar
%
% 2021 Mark C Nelson, McConnell Brain Imaging Centre, MNI, McGill
%--------------------------------------------------------------------------
%% Optional inputs
nin                     = max(nargin,1) - 1;
defaults                = {'',[],[-2558 788 560 420],0};
[str,pinfo,ppos,cmmod]  = INhandler(varargin,nin,defaults);

%% Function

% Check dimensionality and average over 3rd dim if necessary
if length(size(D))==3
    D                   = groupavg(D,3,'nz');                               % Assumes 3rd dim is subjects!
end

D(D==0) = nan;

% Get community info for plot
if ~isempty(pinfo)
    [it,cil,cirb,cirm]  = conn_useCIs(D,pinfo);
else
    it                  = 1:size(D,1);
end

fig=myfig(str,ppos);
imagesc(D(it,it)); hold on;
axis square; cb=colorbar; title(str); set(gca,'FontSize', 25);

% Smart clims
% clims=caxis; lims=ylim;
clims=[prctile(D(D~=0),1) prctile(D(D~=0),99)]; lims=ylim; 

% dens                = round(density_und(d)*100);                          % network density

try
    caxis(clims);
catch
    disp('* * * NOTE * * * unable to change coloraxis')
end

if exist('cil','var')
    set(gca,'XTick',ceil(cirm'),'XTickLabel',cil','XTickLabelRotation',40);
    set(gca,'YTick',ceil(cirm'),'YTickLabel',cil','YTickLabelRotation',40);
    Mticks = [0 cirb]+lims(1);
    M=mesh(Mticks, Mticks, zeros(numel(Mticks)));                           % Overlay grid separating communities
    M.FaceColor='none'; M.EdgeColor='k';
end


% Modify color mapping if 0 is not an extreme in data
if cmmod == 1
    if all(0>clims(1) & 0<clims(2)) || all(0<clims(1) & 0>clims(2))         % 0 lies between color limits
        
        % External function
%         makezeroblack(D,ax)

        % Manual version
        d                   = D(:);
        lssd                = 1e5;                                          % linspace sampling density (precision of 0 value marking)
        ds                  = linspace(min(clims),max(clims),lssd);         % linearly spaced range of data
        dpn                 = find(ds==max(ds(ds<0)))/lssd;                 % proportion of data below 0
        modmap              = colormap;                                     % starting map
        ncol                = size(modmap,1);                               % range of colors
        zpos                = 1 + floor(dpn * ncol);                        % find row corresponding to 0
        modmap(zpos,:)      = [0 0 0];                                      % set that position to black
        colormap(modmap)                                                    % activate it
    end
end

%% Optional output
if nargout > 0
    varargout = cell(1,nargout);
    for oo = 1 : nargout
        switch oo
            case 1
                varargout{1}    = fig;
            case 2
                plt             = gca;
                varargout{2}    = plt;
            case 3
                varargout{3}    = cb;
        end
    end
end

% --------------- Additional unused colormap mod code --------------- %
%% Setting 0s to black
% Set 0 to black at existing location in caxis (not centered on 0)
% D                   = D(:);
% lssd                = 1000000;                                              % linspace sampling density (precision of 0 value marking)
% ds                  = linspace(min(D), max(D), lssd);                       % linearly spaced range of data
% dpn                 = find(ds==max(ds(ds<0)))/lssd;                         % proportion of data below 0
% 
% modmap              = colormap;                                             % starting map
% ncol                = size(modmap,1);                                       % range of colors
% zpos                = 1 + floor(dpn * ncol);                                % find row corresponding to 0
% modmap(zpos,:)      = [0 0 0];                                              % set that position to black
% colormap(modmap);                                                           % activate it
% 
% 
% % Centering caxis on 0, then set to black
% D                   = D(:);
% caxis(max(abs(D)) * [-1 1])                                                 % center color axis on 0
% modmap              = colormap;                                             % starting map
% ncol                = size(modmap,1);                                       % range of colors
% zpos                = 1+floor(.5 * ncol);                                   % find this row in colormap
% modmap(zpos,:)      = [0 0 0];                                              % set that position to black
% colormap(modmap);                                                           % activate it


%--------------------------------------------------------------------------
end