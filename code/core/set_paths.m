function paths = set_paths()
    % Initialize the paths struct with empty fields
    paths = struct( ...
        'project_path', '', ...      
        'cplex', '', ...
        'data', '', ...
        'models', '', ...
        'toolbox_path', '', ... % This is the correct key for the toolbox path
        'matTFA', '', ...
        'redgem', '', ...
        'redhuman', '', ...
        'remi', '', ...
        'BCT', '', ...
        'imm', '', ... 
        'annotater', '', ...
        'spams', '', ...
        'fastica', '');

    % The starting path for the uigetdir dialogs can be set to the user's home directory
    start_path = getenv('HOME');

    % Ask the user to select the project path
    paths.project_path = uigetdir(start_path, 'Select the Project Path');
    if paths.project_path == 0
        error('User cancelled the project path selection.');
    end

    % Ask the user to select the toolbox path
    paths.toolbox_path = uigetdir(start_path, 'Select the parent folder of matTFA');
    if paths.toolbox_path == 0
        error('User cancelled the toolbox path selection.');
    end

    % Ask the user to select the solver path
    paths.cplex = uigetdir(start_path, 'Select the folder that contains IBM ILOG CPLEX');
    if paths.cplex == 0
        error('User cancelled the solver path selection.');
    end

    % Update the dependent paths
    paths.data = fullfile(paths.project_path, 'data');
    paths.models = fullfile(paths.project_path, 'models');
    paths.matTFA = fullfile(paths.toolbox_path, 'matTFA'); % Corrected the path for matTFA
    paths.redgem = fullfile(paths.toolbox_path, 'redgem'); 
    paths.redhuman = fullfile(paths.toolbox_path, 'redhuman');
    paths.remi = fullfile(paths.toolbox_path, 'lcsb_remi'); 
    paths.BCT = fullfile(paths.toolbox_path, 'Brain_Connectivity_Toolbox');
    paths.imm = fullfile(paths.toolbox_path, 'imm'); 
    paths.annotater = fullfile(paths.toolbox_path, 'omar&denis'); 
    paths.spams = fullfile(paths.toolbox_path, 'spams-matlab-v2.6'); 
    paths.fastica = fullfile(paths.toolbox_path, 'FastICA_25');

    % Confirm with the user that these are the correct paths before proceeding
    disp('The following paths have been set:');
    field_names = fieldnames(paths);
    for i = 1:numel(field_names)
        disp([field_names{i} ': ' paths.(field_names{i})]);
    end

    % Ask the user which paths are incorrect, if any
    incorrectPath = true;
    while incorrectPath
        choice = questdlg('Are these paths correct?', ...
            'Confirm Paths', ...
            'Yes', 'No', 'Yes');
        
        if strcmp(choice, 'Yes')
            incorrectPath = false;
        else
            % Allow the user to correct paths
            for i = 1:numel(field_names)
                choice = questdlg(['Is the path for ' field_names{i} ' correct?'], ...
                    'Path Check', 'Yes', 'No', 'Yes');
                if strcmp(choice, 'No')
                    paths.(field_names{i}) = uigetdir(start_path, ['Select the correct path for ' field_names{i}]);
                    if paths.(field_names{i}) == 0
                        error(['User cancelled the selection for ' field_names{i} '.']);
                    end
                    % Update dependent paths if they are based on the path just corrected
                    if strcmp(field_names{i}, 'project_path')
                        paths.data = fullfile(paths.project_path, 'data');
                        paths.models = fullfile(paths.project_path, 'models');
                    elseif strcmp(field_names{i}, 'toolbox_path')
                        paths.matTFA = fullfile(paths.toolbox_path, 'matTFA');
                        paths.redgem = fullfile(paths.toolbox_path, 'redgem');
                        paths.redhuman = fullfile(paths.toolbox_path, 'redhuman');
                        paths.remi = fullfile(paths.toolbox_path, 'lcsb_remi');
                        paths.BCT = fullfile(paths.toolbox_path, 'Brain_Connectivity_Toolbox');
                        paths.imm = fullfile(paths.toolbox_path, 'imm');
                        paths.annotater = fullfile(paths.toolbox_path, 'omar&denis');
                        paths.spams = fullfile(paths.toolbox_path, 'spams-matlab-v2.6');
                        paths.fastica = fullfile(paths.toolbox_path, 'FastICA_25');
                    end
                end
            end
        end
    end

    % Remove the toolboxes path
     paths = rmfield(paths,"toolbox_path"); 

    % Return the struct only if all paths are confirmed by the user
    return;
end
