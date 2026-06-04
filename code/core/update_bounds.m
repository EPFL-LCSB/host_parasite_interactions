function model = update_bounds(model,rxns,lbs,ubs,TFA,verbose)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INPUTS: 
% - rxns: string array containing the name(s) of the rxn(s) for which we
%  want to update the bounds
% - lbs: lower bound for the reaction(s), in a FBA format 
% - ubs: upper bound for the reaction(s), in a FBA format 
% - TFA: Boolean flag indicating if one should consider updating the bounds
% at the TFA levels

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Allows for flexibility

if length(lbs) == 1
    lbs = repmat(lbs,length(rxns),1);
end

if length(ubs) == 1
    ubs = repmat(ubs,length(rxns),1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loop through the reactions

fprintf('##########################\n');

for i = 1:length(rxns)

    %% Get more practical names
    rxn = rxns(i);
    ind = strcmp(model.rxns,rxn);
    if ind == 0
        error("reaction not found")
    end
    ind_F = strcmp(model.varNames,strcat('F_',rxn));
    ind_R = strcmp(model.varNames,strcat('R_',rxn));
    ind_NF = strcmp(model.varNames,strcat('NF_',rxn));
    lb = lbs(i);
    ub = ubs(i);
    
    %% Fill with the actual values if it is empty
    if isnan(lb)
       lb = model.lb(ind);
    end        
    if isnan(ub) 
        ub = model.ub(ind);
    end
        
    %% Impose the bounds at the FBA level
    model.lb(ind) = lb;     
    model.ub(ind) = ub;
              
    %% Check 
    formatSpec = '(FBA) For the reaction %s, lb = %f and ub = %f \n';
    fprintf(formatSpec,string(rxn),model.lb(ind),model.ub(ind));
    
    %% Impose the bounds at the TFA level
    if TFA

        % Impose the bounds by considering the different possibilities
        if lb < 0 && ub < 0
            model.var_ub(ind_R) = abs(lb); model.var_lb(ind_R) = abs(ub);
            model.var_ub(ind_F) = 0; model.var_lb(ind_F) = 0;
        elseif lb <= 0 && ub >= 0
            model.var_ub(ind_R) = abs(lb); model.var_lb(ind_R) = 0;
            model.var_ub(ind_F) = ub; model.var_lb(ind_F) = 0;
        elseif lb > 0 && ub > 0
            model.var_ub(ind_R) = 0; model.var_lb(ind_R) = 0;
            model.var_ub(ind_F) = ub; model.var_lb(ind_F) = lb;
        else
            error("PROBLEM WHEN IMPOSING BOUNDS AT TFA LEVEL!");
        end

        % Impose the constraints on the NF
        model.var_lb(ind_NF) = lb;
        model.var_ub(ind_NF) = ub;

        %% Check 
        formatSpec = '(TFA) For the variable %s, lb = %f and ub = %f \n(TFA) For the variable %s, lb = %f and ub = %f \n';
        fprintf(formatSpec,string(model.varNames(ind_R)),model.var_lb(ind_R), model.var_ub(ind_R),string(model.varNames(ind_F)),model.var_lb(ind_F), model.var_ub(ind_F));

    end
    
    %% Separator
    fprintf('##########################\n');
    
end


%% Overwrite the function fprintf
function fprintf(varargin)
    if verbose
      builtin('fprintf', varargin{:});
    end
  end



end