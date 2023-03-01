function P = a0_pathDefs(runloc)
% Defines path structure
%
% Input
%   runloc  : switch for path setup (string)
%
% Output
%   P       : paths (structure)
%
% 2020 Mark C Nelson, McConnell Brain Imaging Centre, MNI, McGill
%--------------------------------------------------------------------------
%%

switch runloc
    
    case 'local'
        
        disp('Getting LOCAL path defs...')
        P.root              = '~/Desktop';
        P.HOME              = [P.root '/01_aim1'];
        P.MICAPIPE          = [P.HOME '/micapipe'];
        P.project           = [P.HOME '/someProject'];
        P.dataSource        = [P.HOME '/sourceData'];
        P.dataSC            = [P.dataSource '/sc'];
        P.dataFC            = [P.dataSource '/fc'];                         % If stored separately
        P.dataLabel         = [P.MICAPIPE '/parcellations'];
        P.dataDerivative    = [P.project '/derivatives'];
        
    case 'remote'
        
        disp('Getting REMOTE path defs...')
        P.root              = '/PATH/TO/YOUR/REMOTE/STUFF';
        P.HOME              = [P.root '/01_aim1'];
        P.project           = [P.HOME '/someProject'];
        P.MICAPIPE          = [P.HOME '/micapipe'];
        P.dataSource        = [P.HOME '/sourceData'];
        P.dataSC            = [P.dataSource '/sc'];
        P.dataFC            = [P.dataSource '/fc'];
        P.dataLabel         = [P.MICAPIPE '/parcellations'];
        P.dataDerivative    = [P.project '/derivatives'];
        
        
    otherwise
        
        disp('ERROR: to get path defs, supply keywords: local, remote, etc')
end

% Create output dir if necessary
if ~isdir(P.dataDerivative)
    disp(['Creating output directory: ' P.dataDerivative])
    mkdir(P.dataDerivative)
end

% print output
fprintf('\n*--------------------------------  PATH DEFS  --------------------------------*\n')
disp([ 'root        : '  P.root           ])
disp([ 'project     : '  P.project        ])
disp([ 'SC data     : '  P.dataSC         ])
disp([ 'FC data     : '  P.dataFC         ])
disp([ 'labels      : '  P.dataLabel      ])
disp([ 'output      : '  P.dataDerivative ])
fprintf('*-----------------------------------------------------------------------------*\n')

return;