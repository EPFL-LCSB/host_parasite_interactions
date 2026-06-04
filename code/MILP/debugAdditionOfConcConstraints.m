function [constraintsToApply,constraintsToRelax,model,sol, IndOfDebugVars] = debugAdditionOfConcConstraints(model,ConstraintTable,timeLimit)
% When introducing constraints, find which ones have to be relaxed in order
% to get a feasible solution.
%
% INPUTS:
% - model: Input TFBA mode
% - ConstraintTable: Set of variable names for Log concentrations, or
%         reaction fluxes. This should be a cell with three columns:
%         1st column: variable names
%         2nd column: lower bound
%         3rd column: upper bound
% 
% OUTPUTS:
% - constraintsToApply: Which constraints can be immediately applied
% - constraintsToRelax: which of the constraints should be relaxed

% Initialize a variable to store all the indices of the new variables for
% debugging
IndOfDebugVars = [];

% Reset the vector of the objective function (zeros)
model.f = zeros(length(model.f),1);

% Find the indices of the variables to be constrained in the varNames
[~, Constr_idx] = ismember(ConstraintTable(:,1), model.varNames);
[~,num_vars] = size(model.A);

% Some contraints of the conc not present in smaller models so remove zeros
Constr_idx = Constr_idx(find(Constr_idx));
[~, id_conc_add]=ismember(model.varNames(Constr_idx),ConstraintTable(:,1));
ConstraintTable=ConstraintTable(id_conc_add,:);

% Create new binary variable LCUSE,  
for i=1:length(ConstraintTable(:,1))

    % Create a binary variable DebugRelax_LC_met
    model.varNames{num_vars+i,1} = strcat('DebugRelax_',model.varNames{Constr_idx(i)});
    model.var_ub(num_vars+i,1) = 1;
    model.var_lb(num_vars+i,1) = 0;
    model.vartypes{num_vars+i,1} = 'B';

    % set the vector of the objective function for the newly generated
    % variables to 1
    model.f(num_vars+i)=1; 
    IndOfDebugVars(i,1)=num_vars+i;

end

% Create a new constraint for the newly generated variable
[num_constr,~] = size(model.A);

% from left to right:
% LC_met_i < new_ub_LC_met_i (negative) + (diff_btwn_old_and_new_ub_for_LC_met)*binary_var_DebuG_to_minimize
for i=1:length(ConstraintTable(:,1))

    % Add the name for the constraint
    model.constraintNames{num_constr+i} = strcat('DebugRelax_ub_', model.varNames{Constr_idx(i)});

    % Link the variables LC_met with this particular constraints
    % in the A matrix and put a weight of 1
    model.A(num_constr+i, Constr_idx(i))=1;

    % Link the variable DebugRelax_LC_met with this particular constraints
    % in the A matrix and put a weight of the difference between the new
    % and old bounds
    model.A(num_constr+i, IndOfDebugVars(i))= (-model.var_ub(Constr_idx(i))) + ConstraintTable{i,3};

    model.constraintType{num_constr+i}='<';
    model.rhs(num_constr+i)= ConstraintTable{i,3}; % new upper bound for the LConcentration

end
 
% Create one more constraint involving the newly generated variable
[num_constr,~] = size(model.A);

for i=1:length(ConstraintTable(:,1))

    model.constraintNames{num_constr+i} = strcat('DebugRelax_lb_', model.varNames{Constr_idx(i)});
    model.A(num_constr+i, Constr_idx(i))=1;
    model.A(num_constr+i, IndOfDebugVars(i)) = (-model.var_lb(Constr_idx(i))) + ConstraintTable{i,2};
    model.constraintType{num_constr+i}='>';
    model.rhs(num_constr+i)= ConstraintTable{i,2};


end

% Set the solver to minimization of the number of new binary variables: if
% the DebugRelax_LC_met is 1, then this meant that the LC did not change.
model.objtype = +1;

% Solve the system to find the minimal set of new variables that can be
% added
sol = solveTFAmodel_selections(model,'TimeInSec', timeLimit);

% Get the newly added binary variables that were not able to be set to 0
ValuesOfDebugVars = sol.x(IndOfDebugVars);

ids_VarToRelax = find(ValuesOfDebugVars==1);
ids_VarToApply = find(ValuesOfDebugVars==0);

% Those concentrations must be relaxed to get a feasible solution
constraintsToRelax = regexprep(model.varNames(IndOfDebugVars(ids_VarToRelax)),'^DebugRelax_','');
constraintsToApply = regexprep(model.varNames(IndOfDebugVars(ids_VarToApply)),'^DebugRelax_','');

% Check
%sol.x(strcmp(model.varNames,'F_biomass'))
%sol.x(startsWith(model.varNames,'DebugRelax_'))
%checkTable = [ConstraintTable(:,1) cell2mat(ConstraintTable(:,2)),sol.x(Constr_idx),cell2mat(ConstraintTable(:,3))];

end
