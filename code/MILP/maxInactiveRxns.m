function [rxns2add,uptakes_pb,secretions_pb] = maxInactiveRxns(model_red,model_full,time,model_prep_save_file,uptakes,uptakes_save_file,secretions,secretions_save_file,thermo_db)

%% Identify the difference in reactions and update the bounds of the full models with the bounds of the reduced model
if ~exist(model_prep_save_file,'file')

    rxns_same = intersect(model_red.rxns,model_full.rxns);
    
    [~,id_rxns_same_full] = ismember(rxns_same,model_full.rxns);
    [~,id_rxns_same_r] = ismember(rxns_same,model_red.rxns);
    
    model_full_w_updated_bounds = model_full;
    
    model_full_w_updated_bounds.lb(id_rxns_same_full) = model_red.lb(id_rxns_same_r);
    model_full_w_updated_bounds.ub(id_rxns_same_full) = model_red.ub(id_rxns_same_r);
    
    solveFBAmodel_selections(model_full_w_updated_bounds)
    
    model_full_w_updated_bounds = FBA2TFA(model_full_w_updated_bounds,thermo_db,1e3,5e4,'DGo');
    
    solveAndCheckGrowthFBAandTFAmodel(model_full_w_updated_bounds);
    

    %% Add variables and constraints to maximise the number of inactive different rxns
    
    % Get the rxns that are different and add the constraints
    
    rxns_diff = setdiff(model_full_w_updated_bounds.rxns,model_red.rxns);
    rxns_diff = setdiff(rxns_diff,model_full_w_updated_bounds.rxns(find(endsWith(model_full_w_updated_bounds.rxns,'_test'))));
    model_full_w_updated_bounds = addBinaryVarAndCons(model_full_w_updated_bounds,rxns_diff);
    
    %% Impose minimal flux through the drains of interest, depending on uptakes or secretions
    
    % Save the original model
    model_full_w_updated_bounds_original = model_full_w_updated_bounds;
    inds_bfunuse_original = find(startsWith(model_full_w_updated_bounds_original.varNames,'BF_UNUSED_'));    

    save(model_prep_save_file);
else
    load(model_prep_save_file);    
end

% Uptakes
[rxns2add_uptakes,uptakes_pb] = Rxns2Add(model_full_w_updated_bounds_original,inds_bfunuse_original,uptakes,'uptake',time);
save(uptakes_save_file)

% Secretions
[rxns2add_secretions,secretions_pb] = Rxns2Add(model_full_w_updated_bounds_original,inds_bfunuse_original,secretions,'secretion',time);
save(secretions_save_file)


%% Final processing
rxns2add = union(rxns2add_uptakes,rxns2add_secretions);

end


%% Helper functions

function model = addBinaryVarAndCons(model,rxns)

    [num_constr,num_vars] = size(model.A);

    for i = 1:length(rxns)

        % Add the variables
        model.varNames{num_vars+i,1}    = strcat('BF_UNUSED_',rxns{i});
        model.vartypes{num_vars+i,1}    = 'B';    
        model.var_lb(num_vars+i,1)      = 0;
        model.var_ub(num_vars+i,1)      = 1;
    
        % Get the indices for the constraints
        ind_bfunuse     = find(ismember(model.varNames,strcat('BF_UNUSED_',rxns{i})));
        ind_buse        = find(ismember(model.varNames,strcat('BU_',rxns{i})));
        ind_fuse        = find(ismember(model.varNames,strcat('FU_',rxns{i})));
    
        % Add the constraints: BFUNUSED + BUSE + FUSE = 1 
        model.constraintNames{num_constr+i,1}   = strcat('cons_UNUSED_',rxns{i});
        model.A(num_constr+i, ind_bfunuse)      = 1;
        model.A(num_constr+i, ind_buse)         = 1;
        model.A(num_constr+i, ind_fuse)         = 1;
        model.constraintType{num_constr+i,1}    = '=';
        model.rhs(num_constr+i,1)               = 1;    
    
        % Objective function
        model.f(num_vars+i,1)=1; 

    end

end

function [rxns2add,exc_pb] = Rxns2Add(model_original,inds_bfunuse,exc_rxns,exc_type,time)

    if strcmp(exc_type,'uptake')
        F_lb_bound = 1e-6;
        F_ub_bound = 50;
        R_lb_bound = 0;
        R_ub_bound = 0;
    elseif strcmp(exc_type,'secretion')
        F_lb_bound = 0;
        F_ub_bound = 0;
        R_lb_bound = 1e-6;
        R_ub_bound = 50;
    else
        error("wrong type of exchange")
    end

    rxns2add = [];
    exc_pb = [];

    for i = 1:length(exc_rxns)
    
        % Initiate variable:

        % Restart with a fresh model
        model = model_original;     
    
        % Identify the metabolite we want to produce/expel in/from cytosol
        met         = strcat(exc_rxns(i),'_test');
        dmd         = strcat('DM_',met);
        F_dm        = strcat('F_',dmd);
        R_dm        = strcat('R_',dmd);
        NF_dm       = strcat('NF_',dmd);
        ind_F_dm    = find(ismember(model.varNames,F_dm));
        ind_R_dm    = find(ismember(model.varNames,R_dm));
        ind_NF_dm   = find(ismember(model.varNames,NF_dm));

        % Block the other test rxns
        ind_test_variables = find(endsWith(model.varNames,'_test') & ~contains(model.varNames,met));
        model.var_lb(ind_test_variables) = 0;
        model.var_ub(ind_test_variables) = 0;

    
        % Impose minimal flux through the secretion of the metabolites
        model.var_lb(ind_F_dm) = F_lb_bound;
        model.var_ub(ind_F_dm) = F_ub_bound;
        model.var_lb(ind_R_dm) = R_lb_bound;
        model.var_ub(ind_R_dm) = R_ub_bound;
        model.var_lb(ind_NF_dm) = -50;
        model.var_ub(ind_NF_dm) = 50;

        % Maximise the number of inactive rxns
        sol = solveTFAmodel_selections(model,'TimeInSec',time,'timePolishing',time*2/3);

        % Look for alternatives
       if ~(isempty(sol.x) || sum(isnan(sol.x)) > 0)

            % Look for alternatives
            [DPs,~,~,~] = find_alt(model,10,sol,inds_bfunuse,time,0.5,sol.val);

            % Select the union of alternatives
            mins_2add = inds_bfunuse(DPs(inds_bfunuse,1) == 0);
            if size(DPs,2) > 1
                for j = 2:size(DPs,2)
                    mins_2add = union(mins_2add,inds_bfunuse(DPs(inds_bfunuse,j) == 0));
                end
            end
            mins_2add = strrep(model.varNames(mins_2add),'BF_UNUSED_','');
    
            % Update the pool of rxns
            rxns2add = union(rxns2add,mins_2add);

       else

           exc_pb = [exc_pb;met];

       end
        
    end

end

% Seeking alternatives
function [DPs, objectives, model, count_DP_YIELD] = find_alt(model, maxNumEPs, sol, use_vectors_indices, TimeLimitSecs, SizeLUMPNetwork, stop_search_value)

    num_sols = 0;
    count = 1;
    smallest_subnetwork = Inf;
    actUSEvecs = {};

    while (num_sols < maxNumEPs) && ~isempty(sol)

        fprintf('#######################################\n')
        fprintf('Looking for DP number %d\n', count+1)

        [~,num_vars] = size(model.A);


        if ~isempty(sol.x)
            num_sols = num_sols + 1;
            objectives(num_sols, 1) = sol.val;
            DPs(:, num_sols) = sol.x;

    
            % we find all the use vectors and formulate them into a new integer cut
            % constraint
            use_vectors_sol = ones(num_vars,1);
            use_vectors_sol(use_vectors_indices) = sol.x(use_vectors_indices);
            actUSEvec = find(use_vectors_sol<1e-9);
            actUSEvecs = [actUSEvecs,{actUSEvec}];

            if count == 1
                smallest_subnetwork = length(actUSEvec);
            end

            % Update the model with a new integer cut constraint
            model = addIntegerCutConstraint(model, actUSEvec, num_sols);

            % Solve the TFA model
            sol = solveTFAmodel_selections(model,'TimeInSec',TimeLimitSecs ,'stop_search_value', stop_search_value);

            if isempty(sol.x)
                break;
            end

            length_sol = nnz(sol.x(use_vectors_indices) < 1e-9);

            if abs(length_sol - smallest_subnetwork) > SizeLUMPNetwork
                fprintf('In findDP_YIELD_TimeLimit.m: break because size of subnetwork > size of smallest subnetwork by %d variable(s)\n: ', length_sol - smallest_subnetwork); 
                break;
            end

            fprintf('Number of DPs:\t%d\n', num_sols);


            % Sanity check 
            if rem(num_sols,3) == 0
                 
                % Sanity check 1: make sure all the elements are same size or smaller size
                vectorSizes = cellfun(@numel, actUSEvecs);
                max_size = size(drains,1) - stop_search;
                isValid = all(vectorSizes <= max_size);
                
                if ~isValid
                    error('Error: At least one vector in the cell array has a size larger than %d.\n',max_size);
                else
                    fprintf('All vectors in the cell array have a size equal to or less than %d.\n',max_size);
                end
        
                % Sanity check 2 to ensure actUSEvec is different from previous iterations
                if size(unique(DPs(use_vectors_indices,:)','rows'),1) < size(DPs,2)
                    error('New actUSEvec is the same as a previous one. Aborting search.');
                end
        
            end


        end
        count = count + 1;
    end

    fprintf('Analysis for DPs ends with %d alternatives\n', count);
    count_DP_YIELD = count;
    fprintf('#################################################################################\n')

end


function updatedModel = addIntegerCutConstraint(model, active_use_vectors, num_sols)
    [num_cons, num_vars] = size(model.A);
    updatedModel = model;
    updatedModel.constraintNames{num_cons + 1} = ['CUT_' num2str(num_sols)];
    new_constraint = zeros(num_vars, 1);
    new_constraint(active_use_vectors) = 1;
    updatedModel.constraintType{num_cons + 1} = '>';
    updatedModel.A(num_cons + 1, :) = new_constraint;
    updatedModel.rhs(num_cons + 1) = 0.5;
end
