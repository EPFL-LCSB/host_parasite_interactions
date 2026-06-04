function solver_config = set_solver()
    % Set default values for name and type
    default_name = 'cplex_direct';
    default_type = 'LP';

    % Ask the user for the name and type of the solver using inputdlg
    prompt = {'Enter the name of the solver:', 'Enter the type of the solver:'};
    dlgtitle = 'Set Solver';
    dims = [1 50];
    definput = {default_name, default_type};
    answer = inputdlg(prompt, dlgtitle, dims, definput);

    % Get the user's answers from the answer cell array
    name = answer{1};
    type = answer{2};

    % Create the solver_config structure with the name and type fields
    solver_config.name = name;
    solver_config.type = type;
end
