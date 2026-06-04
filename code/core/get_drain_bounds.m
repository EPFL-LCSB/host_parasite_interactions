function [drains_F, drains_R, drains_NF, drains_F_lb, drains_F_ub, drains_R_lb, drains_R_ub, drains_NF_lb, drains_NF_ub, drains_F_ind, drains_R_ind, drains_NF_ind] = get_drain_bounds(model)

% extract drains and remove any variable containing 'biomass' (ignoring case)
drains = extract_drains(model);
drains(contains(lower(drains), 'biomass')) = [];

% create drain names with prefixes
drains_F = strcat('F_', drains);
drains_R = strcat('R_', drains);
drains_NF = strcat('NF_', drains);

% get indices of drain variables in model.varNames
[~, drain_ind] = ismember([drains_F, drains_R, drains_NF], model.varNames);
drains_F_ind = drain_ind(1:numel(drains_F));
drains_R_ind = drain_ind(numel(drains_F)+1:numel(drains_F)+numel(drains_R));
drains_NF_ind = drain_ind(numel(drains_F)+numel(drains_R)+1:end);

% get lower and upper bounds of drain variables
drains_F_lb = model.var_lb(drains_F_ind);
drains_F_ub = model.var_ub(drains_F_ind);

drains_R_lb = model.var_lb(drains_R_ind);
drains_R_ub = model.var_ub(drains_R_ind);

drains_NF_lb = model.var_lb(drains_NF_ind);
drains_NF_ub = model.var_ub(drains_NF_ind);