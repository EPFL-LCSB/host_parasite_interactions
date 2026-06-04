function [DPs, objectives, model, count_DP_YIELD] = findDP_YIELD_TimeLimit_CallBack(model, maxNumEPs, sol, use_vectors_indices, TimeLimitSecs, SizeLUMPNetwork, stop_search_value)

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
                 
                % Sanity check to ensure actUSEvec is different from previous iterations
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

%%
% Helper function
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
