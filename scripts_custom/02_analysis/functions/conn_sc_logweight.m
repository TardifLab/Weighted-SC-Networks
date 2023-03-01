function [D,varargout] = conn_sc_logweight(D,varargin)
%  log  transforms the datasets in D indicated by TARGET.
%  If no target given, datasets are log transformed based on skewness.
%
% Input:
%                           + Required +
%   D       : IxJ cell array, each cell containing an NxNxS stack of
%             connectivity matrices
%
%                           + Optional +
%   labels  : Cell array of strings describing edge weights in D
%   target  : cell array of strings describing weights to transform 
%
% Output:
%   D       : modified data
%   labels  : modified to reflect change (optional)
%
% 2021 Mark C Nelson, McConnell Brain Imaging Centre, MNI, McGill
%--------------------------------------------------------------------------

%% Optional inputs
nin=max(nargin,1) -1; 
defaults            = {{} {}};
[labels,target]     = INhandler(varargin,nin,defaults);


%% setup
if ~iscell(D);      D       = {D};                                      end
if ~iscell(labels); labels  = {labels};                                 end
if ~iscell(target); target  = {target};                                 end
[~,J]                       = size(D);


% Get target indices
if ~isempty(target)                                                         % if target weight given
    itrg                    = cellstrfind(labels, target);
    itrg                    = unique(cell2mat(itrg));
else                                                                        % else use skewness of dataset to find targets (arbitrary value of 5 chosen here!)
    itrg                    = find(cellfun(@(c) abs(skewness(c(:)))>5, D(:,1))); % assumes all columns contain same weights
    itrg                    = itrg';
end
 
% Ensure indices properly set
if isempty(itrg)
    error('* ERROR: Unable to log convert your data!')
end

%% log conversion
for kk = itrg
    
    % Make sure this dataset hasn't already been transformed
    if ~isempty(labels) && ~isempty(regexp(labels{kk}, 'log','ignorecase'))
        disp(['ERROR: data already log transformed at index: ' num2str(kk)])
        continue
    end
        
    for jj = 1 : J
        dt                  = D{kk,jj};
        dt(dt==0)           = nan;                                          % avoids -inf values
        dt                  = log10(dt);
        dt(dt==0)           = 1e-30;                                        % differentiate between 0s pre & post log (shouldnt be any)
        dt(isnan(dt))       = 0;                                            % restore original 0s 
        D{kk,jj}            = dt;
    end
    if ~isempty(labels)
        labels{kk}          = [labels{kk} '(log10)'];
    end
end

%% Optional outputs
nout                            = max(nargout,1) - 1;
if nout ~= 0
    varargout                   = cell(1,nout);
    for kk = 1:nout
        switch kk
            case 1
                varargout{kk}   = labels;
        end  
    end
end

%--------------------------------------------------------------------------
end