function step1_context_specific_host_GEM_reconstruction()

%% Prepare the environment

% Load required data
load('./data/thermoData_Recon3.mat');  % thermo database for Recon3D model
load('./data/PmetsRecon3_redHUMAN.mat');  % Load table describing phosphate group metabolites

% Load thermocurated Recon3D model
model = load('./models/GEM_Recon3_thermocurated_redHUMAN.mat').model;

% Set time limit
timeLimit = 2*60;

%% Load the tables

% Load fluxomics and metabolomics data
Tuptake = readtable('./data/Uptakes.csv', 'ReadVariableNames', true);
Tconc = readtable('./data/MetabolomicsData.csv', 'ReadVariableNames', true);

% Find unique elements in Tconc.nutrients and their first occurrence and
% keep only unique
[~, uniqueIdx] = unique(Tconc.nutrients, 'stable');
Tconc = Tconc(uniqueIdx, :);

% Analyze fluxomics and metabolomics tables
medium_mets_uptake = strrep(table2cell(Tuptake(:, 1)), 'EX_', '');
medium_rxns_uptake = Tuptake.nutrients;
medium_rxns_uptake_all_related = unique(model.rxns(any(model.S(ismember(model.mets, medium_mets_uptake), :), 1)));

%% Update key bounds 

% Set constants
hepato_growth_rate_lb = 0; 
hepato_growth_rate_ub = 0.014;
ATPM_rate = 1.07;

% Update bounds
ATPm_rxn = 'ATPM';
biomass_rxns = model.rxns(contains(model.rxns, 'biomass'));
rxns = [ATPm_rxn; biomass_rxns];
lbs = [ATPM_rate, hepato_growth_rate_lb, 0, 0];
ubs = [NaN, hepato_growth_rate_ub, 0, 0];
TFA = 1;
verbose = 1;
model = update_bounds(model, rxns, lbs, ubs, TFA, verbose);

% Remove sink and Pmets reactions
DM_sink_rxns = DM_sink_rxnsANALYSIS(model);
DM_sink_rxns_2rm = setdiff(DM_sink_rxns, union('EX_pi_e', medium_rxns_uptake_all_related));

[Prxns, ~] = biggMetsANALYSIS(model, Pmets);
Prxns_2rm = setdiff(Prxns, union('EX_pi_e', medium_rxns_uptake_all_related));

model = removeRxns(model, setdiff([DM_sink_rxns_2rm'; Prxns_2rm], medium_rxns_uptake_all_related));

% Set bounds for medium reactions
model = update_bounds(model, medium_rxns_uptake, Tuptake.lb, Tuptake.ub, 0, 1);

% Block other uptakes
drains = extract_drains(model);
other_uptks = setdiff(drains, Tuptake.nutrients);
model = update_bounds(model, other_uptks, 0, NaN, 0, 1);

% Set bounds for key inorganic compound exchanges
openflux = {'EX_k_e', 'EX_na1_e', 'EX_pi_e'};
model = update_bounds(model, openflux, -100, NaN, 0, 1);

%% Conversion model to thermo format

% Test Feasibility
solTFA = solveTFAmodel_selections(model);

% Conversion
model = FBA2TFA(model,thermoData_Recon3,1e3,1e6);

% Check for growth

[solFBA, solTFA] = solveAndCheckGrowthFBAandTFAmodel(model);
disp(solFBA.f);
disp(solTFA.val);

%% Impose constraints for metabolite concentration ranges

% Create log-scaled constraint table
Tconc_cell = table2cell(Tconc);
Tconc_cell(:, 2:3) = num2cell(log(cell2mat(Tconc_cell(:, 2:3))));
ConstraintTable = Tconc_cell;

% Impose a minimum growth of 90% of the max observed growth rate
model_test = update_bounds(model, "biomass", hepato_growth_rate_ub*0.90, hepato_growth_rate_ub, 1, 1);

% Determine constraints to apply and relax
[~, constraintsToRelax, modelRelaxed, solRelaxed, ind_bfuse] = debugAdditionOfConcConstraints(model_test, ConstraintTable, timeLimit);

% Find alternative sets of constraints to relax
[maxNumEPs, sol, use_vectors_indices, TimeLimitSecs, SizeLUMPNetwork] = deal(10, solRelaxed, ind_bfuse, 10*60, 0.5);
[DPs, ~, ~, ~] = findDP_YIELD_TimeLimit_CallBack(modelRelaxed, maxNumEPs, sol, use_vectors_indices, TimeLimitSecs, SizeLUMPNetwork, sol.val);

% Analyze results from the alternatives
if size(DPs, 2) > 1
    for i = 1:size(DPs, 2)
        bfuse = DPs(use_vectors_indices, i);
        constr2relax_DP = string(strrep(modelRelaxed.varNames(use_vectors_indices(bfuse > 0.5)), 'DebugRelax_', ''));
        disp(constr2relax_DP);
        constr2relax_DPs(1, i) = {constr2relax_DP};
    end
end

% User has to make a choice in case there are several alternatives
constraintsToRelax_choice = constr2relax_DPs(1,1);

% Impose the updated constraints
constraintsToRelax_choice = constraintsToRelax;
ConstraintTable_up = Tconc(~ismember(Tconc.nutrients, constraintsToRelax_choice), :);

% Get indices and check alignment
[~, ind_mets_model] = ismember(ConstraintTable_up.nutrients, model.varNames);

% Update log concentration
model_up = model;
model_up.var_lb(ind_mets_model(ind_mets_model > 0)) = log(ConstraintTable_up.lb(ind_mets_model > 0));
model_up.var_ub(ind_mets_model(ind_mets_model > 0)) = log(ConstraintTable_up.ub(ind_mets_model > 0));

% Test model solvability
[~, ~] = solveAndCheckGrowthFBAandTFAmodel(model);

% Update model
model = model_up;

%% Adjust model genes
recon3d = load('./models/GEM_Recon3_thermocurated_redHUMAN.mat', 'model');
model = genes_fromGEM(model, recon3d.model);
model = genes_HGNC(model);
clear recon3d;

%% Find the reactions supposed to be active based on transcriptomics data
rxns_w_GPR_expressed = HPA_processing(model);

%% Run a minmax to remove rxns not carrying fluxes

% Conduct the minmax analysis
minmax = runMinMax_selections(model, model.rxns, true);

% Reactions to remove: rxns carrying zero fluxes based on minmax
ind_rxns2rmv = find(abs(minmax(:, 1)) < 1e-9 & abs(minmax(:, 2)) < 1e-9);
rxns2rmv_minmax = model.rxns(ind_rxns2rmv);

% Remove the reactions from the model based on minmax
model = removeRxns(model, model.rxns(ind_rxns2rmv));

% Solve the model
sol_test = solveFBAmodel_selections(model);

%% Conversion model to thermo

% Conversion
model = FBA2TFA(model,thermoData_Recon3,1e3,1e6);

% Check for growth
[solFBA, solTFA] = solveAndCheckGrowthFBAandTFAmodel(model);


%% Max number of simultaneously active reactions based on transcriptomics data

% Update the rxns to min due to the minmax
rxns_w_GPR_expressed_u = setdiff(rxns_w_GPR_expressed, rxns2rmv_minmax);

% Maximise the number of simultaneously active reactions
[model_minRxns, expressedRxns, nonCoreRxnsToAdd, coreRxnsUpdated,coreRxnsOriginal, nonCoreRxnsUpdated, nonCoreRxnsOriginal, problematicRxns,nonCoreRxnsAddedRecords] = ...
    minNumRxns2Add2Core(model,rxns_w_GPR_expressed_u,{'biomass'}, 1e-6, ... 
    15 * 60,[],'./workspaces/step1_minRxns_identification_in_progress.mat');

% See if the model composed of rxnsC grows
rxns2rmv = setdiff(model.rxns, coreRxnsUpdated);
trimmedModel = removeRxns(model, rxns2rmv);

% Check growth ability
trimmedModel.c = zeros(size(trimmedModel.c));
trimmedModel.c(find(ismember(trimmedModel.rxns,'biomass'))) = 1;
solveFBAmodel_selections(trimmedModel)

% Conversion to thermo
model = FBA2TFA(trimmedModel,thermoData_Recon3,1e3,5e4);

% Check growth ability 
[~, ~] = solveAndCheckGrowthFBAandTFAmodel(model);

%% Impose the metabolomics constraints

% Check if all of the constraints related to concentrations can be imposed
ConstraintTable = Tconc(ismember(Tconc.nutrients, model.varNames), :);
ConstraintTable.lb = log(ConstraintTable.lb);
ConstraintTable.ub = log(ConstraintTable.ub);

model_test = update_bounds(model, "biomass", hepato_growth_rate_ub * 0.90, hepato_growth_rate_ub, 1, 1);
solveTFAmodel_selections(model_test,'TimeInSec',2*60,'emphPar',0)

[~, constraintsToRelax, modelRelaxed, solRelaxed, ind_bfuse] = debugAdditionOfConcConstraints(model_test, table2cell(ConstraintTable), timeLimit);

[maxNumEPs, TimeLimitSecs, SizeLUMPNetwork] = deal(100, 10 * 60, 0.5);
[DPs, ~, ~, ~] = findDP_YIELD_TimeLimit_CallBack(modelRelaxed, maxNumEPs, solRelaxed, ind_bfuse, TimeLimitSecs, SizeLUMPNetwork, length(constraintsToRelax));

for i = 1:size(DPs, 2)
    bfuse = DPs(ind_bfuse, i);
    constr2relax_DP = strrep(modelRelaxed.varNames(ind_bfuse(bfuse > 0.5)), 'DebugRelax_', '');
    constr2relax_DPs(1, i) = {string(constr2relax_DP)};
end

% Impose the updated constraints
constraintsToRelax_choice = constr2relax_DPs{1, 2};
ConstraintTable_up = Tconc(~ismember(Tconc.nutrients, constraintsToRelax_choice), :);

[~, ind_mets_model] = ismember(ConstraintTable_up.nutrients, model.varNames);
model_up = model;
model_up.var_lb(ind_mets_model(ind_mets_model > 0)) = log(ConstraintTable_up.lb(ind_mets_model > 0));
model_up.var_ub(ind_mets_model(ind_mets_model > 0)) = log(ConstraintTable_up.ub(ind_mets_model > 0));

model = model_up;

%% Save the final model & workspace

save('./models/hepatocyte_thermocurated.mat', 'model');
save('./workspaces/step1_final.mat')

end