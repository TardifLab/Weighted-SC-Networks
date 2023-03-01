function [varargout] = plot_conn_connsGroup(D,W8s,pinfo,varargin)
%  Plots group level connectivity data in 2D square matrix form.
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
%   sparse  : Sets 0 values to black if 0 is contained within dataset  
%             (default = 0)
%
% Output:
%                           + Optional +
%   FIG     : Figure handles
%   AX      : Axes handles
%
% 2021 Mark C Nelson, McConnell Brain Imaging Centre, MNI, McGill
%--------------------------------------------------------------------------
%% Optional inputs
nin                 = max(nargin,1) - 3;
defaults            = {'nz','Connectomes',0};
[avg,str,sparse]    = INhandler(varargin,nin,defaults);

%% Setup
[I,J]                       = size(D);
fntsz                       = 25;
parcs                       = {pinfo.name};

% Get shape & number of plots
maxsp                       = 6;                                            % max number of subplots
if I > maxsp                                                                % if number of weights > 
    n                       = 3;                                            % dims for subplot
    m                       = 2;
    ppos                    = [-2559 -16 2220 1224];                        % plot size & location
elseif I*J < maxsp && J > 1                                                 % if all parcs could fit on same plot
    n                       = ceil(sqrt(I*J));
    m                       = ceil(I*J/n);
    ppos                    = [-2559 -16 2220 1224];
    
    % Reconfigure variables to plot all parcs together
    Dnew                    = cell(I*J,1);
    lblsnew                 = cell(1,I*J);                                  % will combine 
    Piter                   = nan(1,J);                                     % new parc iterator
    kk                      = 0;
    for ii = 1 : I 
        W8tmp               = W8s{ii};
        for jj = 1 : J
            Ptmp            = parcs{jj};
            parcs{jj}       = '';                                           % null original parcs labels (need variable)
            kk              = kk+1;
            Dnew(kk)        = D(ii,jj);
            lblsnew{kk}     = [W8tmp '(' Ptmp ')'];                         % join labels for weights & parcs
            Piter(kk)       = jj;                                           % Store index for this parc in new shape
        end
    end
    D                       = Dnew;                                         % Replace original data
    W8s                     = lblsnew;                                      
    [I,J]                   = size(D);
else
    n                       = ceil(sqrt(I));
    m                       = ceil(I/n);
    ppos                    = [-2559 -16 2220 1224];
end
AX                          = cell(1,I*J);                                  % Total number of plots
FIG                         = cell(1,ceil((I*J)/(n*m)));                    % Total number of figures
iFIG=1; iAX=1;                                                              % Iterators for storing figure & plot handles

%% Plot
for jj = 1 : J  
    strplt                  = ['Group-level ' str ': ' parcs{jj}];
    FIG{iFIG}=myfig(strplt,ppos);
%     H=suptitle (strplt);
%     H.FontSize=20; H.FontWeight='bold';
    sp                      = 0;
    for ii = 1 : I
        if exist('Piter','var')
            pi              = Piter(ii);
        else
            pi              = jj;
        end
        d                   = D{ii,jj};
        N                   = length(pinfo(pi).labels);                     % node count

        % If you'd like to automatically transform skewed data
%         if abs(skewness(d(:))) > 5                                          % if data skewed (5 arbitrarily chosen)
%             d               = conn_sc_logweight(d);                         % log convert
%             d               = conn_normalize(d,{'range', [1e-6 1]});        % normalize to interval [0 1]
%             d               = d{1};
%         end
        
        if length(size(d)) == 3                                             % if 3D, need to average across subjects
            d               = groupavg(d,3,avg);
        end
        clims               = [prctile(d(d~=0),.01) prctile(d(d~=0),99.99)];% Determine limits for this weight
        dens                = round(density_und(d)*100);                    % network density
        
        % Use community assignments to reorder data if you have them
        if ~isempty(pinfo(pi).cirois)
            [it,cil,cirb,cirm]  = conn_useCIs(d,pinfo(pi));
        else
            it                  = [1:N]';
        end
                        
        % Plot
        sp                  = sp + 1;
        if sp > maxsp
            sp              = 1;
            iFIG            = iFIG+1;
            FIG{iFIG}=myfig(strplt,ppos);
%             H=suptitle (strplt);
%             H.FontSize=20; H.FontWeight='bold';
        end
        AX{iAX}=subplot(m,n,sp); iAX=iAX+1;
        imagesc(d(it,it)); hold on; axis square; colorbar; caxis(clims); % caxis([0 1]);
        title([W8s{ii} ' (' num2str(dens) '%)']); 
        set(gca,'FontSize', fntsz); lims=ylim;
        if exist('cil','var')
            set(gca,'XTick',ceil(cirm')+lims(1),'XTickLabel',cil','XTickLabelRotation',40);
            set(gca,'YTick',ceil(cirm')+lims(1),'YTickLabel',cil','YTickLabelRotation',40);
            Mticks = [0 cirb]+lims(1);
            M=mesh(Mticks, Mticks, zeros(numel(Mticks)));                         % Overlay grid separating communities
            M.FaceColor='none'; M.EdgeColor='k';
        end
        
        % Option to map values to custom blue2red colormap
%         cm=red2bluecmap(ax,0);                                              % Positive values to red & negative values to blue         

        
        % Modify colormap if 0 is within the color axis
        if sparse == 1
            if ~any([all(clims<=0) all(clims>=0)])
                cm=colormap;
                makezeroblack(AX{iAX-1},cm);                                % 0 to black
            end
        end
    end
end


%% Optional output
if nargout > 0
    varargout = cell(1,nargout);
    for oo = 1 : nargout
        switch oo
            case 1; varargout{oo} = FIG;
            case 2; varargout{oo} = AX;
        end
    end
end

%--------------------------------------------------------------------------
end