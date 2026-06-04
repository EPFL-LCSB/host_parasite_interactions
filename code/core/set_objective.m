function model = set_objective(model,objective_rxns,coeffs,direction,TFA)

%% Handle the missing arguments
if ~exist('objective_rxns', 'var') || isempty(objective_rxns)
    objective_rxns = {'biomass'};
end
if ~exist('coeffs', 'var') || isempty(coeffs)
    coeffs = 1;
end
if ~exist('direction', 'var') || isempty(direction)
    direction = 'max';
end
if ~exist('TFA', 'var') || isempty(TFA)
    TFA = 1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('##########################\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FBA level 

% Reset the objective vector
model.c = zeros(size(model.c));

% Find the index of the objective rxns
[~,ind_objs] = ismember(objective_rxns,model.rxns);

% Set the coeffs individually
for i = 1:length(objective_rxns)
    
    ind_obj = ind_objs(i);
    coeff = coeffs(i);
    
    if strcmp(direction,"max")
        dir= 1;
    elseif strcmp(direction,"min")
        dir = -1;
    else
        error("wrong direction");
    end

    model.c(ind_obj) = dir*coeff;
        
end

%fprintf(formatSpec,string(rxn),model.lb(ind),model.ub(ind));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TFA level
if TFA
    
    % Reset the objective vector
    model.f = zeros(size(model.f));
    
    % Find the index of the objective vars
    [~,ind_objs_NF] = ismember(strcat('NF_',objective_rxns),model.varNames);
    
    % Set the coeffs individually
    for i = 1:length(ind_objs_NF)
    
        ind_obj_NF = ind_objs_NF(i);
        coeff = coeffs(i);

        model.f(ind_obj_NF) = coeff;

        if strcmp(direction,"max")
            model.objtype = -1;
        elseif strcmp(direction,"min")
            model.objtype = 1;
        else
            error("wrong direction");
        end
    
    end

    
end
