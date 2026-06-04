function step2_parasitosomes_generation(medium_type)

%% Prepare the environment

% Load model
model = load('./models/iPfa.mat').ipfa4liver;

% Set maximum growth rate
max_growth = log(90000) / (6.5 * 24);
ind_bm = find(contains(model.rxns, 'biomass_c', 'IgnoreCase', true));
ind_bm_F = find(ismember(model.varNames, 'F_EXC_BOTH_Biomass_c'));
model.ub(ind_bm) = max_growth;
model.var_ub(ind_bm_F) = max_growth;

% Load medium and thermo database
ipfaMedium = load('ipfaMedium.mat').ipfaMedium;
load('thermo_data.mat');

% Set time limit
timeLimit = 2 * 60;

%% Conversion to TFA

% Conversion
[model,~] = FBA2TFA(model,DB_AlbertyUpdate,1e3,1e6,'DGo');

% Check growth ability
[solFBA, solTFA] = solveAndCheckGrowthFBAandTFAmodel(model,{'EXC_BOTH_Biomass_c'});
disp(solFBA.f);
disp(solTFA.val);

%% Minmax to identify blocked rxns and clarify the effective size of the model

ind = getAllVar(model,{'NF'});
minmax = runTMinMax_save_selections(model,model.varNames(ind),3*60);
 
blocked_rxns = find(abs(minmax(:,1))<1e-9 & abs(minmax(:,2))<1e-9);
 
model_effective = removeRxns(model,strrep(model.varNames(ind(blocked_rxns)),'NF_',''));

disp(numel(model_effective.rxns));
disp(numel(model_effective.mets));

% Resulting model: 737 metabolites and 889 reactions

%% Enumerate the in silico Minimal Media (iMMs) if the medium is compromised

if strcmp(medium_type,'prototrophic') && ~exist('./workspaces/step2_iMM_recap_table.mat','file')

    % Set the parameters and load data
    sol = solveTFAmodel_selections(model);
    params = struct('flagUpt', true, ...
                    'minObj', 0.99 * sol.val, ...
                    'maxObj', sol.val, ...
                    'NumAlt', 1e5, ...
                    'drainsForiMM', {{}}, ...
                    'rxnNoThermo', [], ...
                    'ReactionDB', load('thermo_data.mat'), ...
                    'metabData', [], ...
                    'save_path','./workspaces/step2_iMM_enumeration.mat');
    
    
    % Enumeration iMM
    [model_iMM, drains_iMM, MinMets_iMM, cDPs_iMM] = ...
        analysisIMM_save(model, params.flagUpt, params.minObj, params.maxObj, ...
                    params.NumAlt, params.drainsForiMM, params.rxnNoThermo, ...
                    params.ReactionDB, params.metabData,params.save_path, 'EXC_BOTH_Biomass_c');
    
    
    %% Analyse the iMMs
     
    % Get the indices of the BF_UNUSE
    ind_BFuse = startsWith(model_iMM.varNames, 'BFUSE_');
    
    % Initialize the exchange_table
    DPs = cDPs_iMM;
    exchange_table = table('VariableTypes', ["cell", "cell"], ...
                           'VariableNames', ["Uptakes", "Secretions"], ...
                           'Size', [size(DPs, 2), 2]);
    
    % Loop through the Distribution Profiles
    for i = 1:size(DPs, 2)
    
        % Get active uptake reactions
        active_uptakes_BFUSE = model_iMM.varNames(ind_BFuse & (DPs(:, i) == 0));
        uptakes = strrep(active_uptakes_BFUSE, 'BFUSE_', '');
        uptakes_flux = DPs(ismember(model_iMM.varNames,uptakes), i);
    
        % Check 
        ind = startsWith(model_iMM.varNames, 'R_EXC_BOTH_');
        uptakes_test = model_iMM.varNames(ind & (DPs(:, i) > 1e-9));
        uptakes_test_flux = DPs(ismember(model_iMM.varNames, uptakes_test), i);
    
        if ~all(isequal(uptakes_flux,uptakes_test_flux))
            error('problem in uptakes fluxes')
        end
    
        % Create and fill the uptakes table
        uptakes_table = table(uptakes, uptakes_flux);
        exchange_table.Uptakes{i} = uptakes_table;
    
        % Check the associated secretions for each set of active uptake reactions
        ind = startsWith(model_iMM.varNames, 'F_EXC_BOTH_');
        secretions = model_iMM.varNames(ind & (DPs(:, i) > 0));
        secretions_flux = DPs(ismember(model_iMM.varNames, secretions), i);
    
        % Create and fill the secretions table
        secretions_table = table(secretions, secretions_flux);
        exchange_table.Secretions{i} = secretions_table;
    end
    
    % Perform the union of uptakes and secretions
    secretions_union = unique(vertcat(exchange_table.Secretions{:}).secretions);
    uptakes_union = unique(vertcat(exchange_table.Uptakes{:}).uptakes);
    
    % Extract the first columns of secretions and uptakes tables
    secretions = arrayfun(@(x) x{1}.secretions, exchange_table.Secretions, 'UniformOutput', false);
    uptakes = arrayfun(@(x) x{1}.uptakes, exchange_table.Uptakes, 'UniformOutput', false);
    
    % Compute the intersection using the first cell as a starting point
    secretions_intersection = secretions{1};
    uptakes_intersection = uptakes{1};
    
    for i = 2:numel(secretions)
        secretions_intersection = intersect(secretions_intersection, secretions{i});
        uptakes_intersection = intersect(uptakes_intersection, uptakes{i});
    end
    
    
    % Save the results
    save('./workspaces/step2_iMM_recap_table.mat')

else 
    load('./workspaces/step2_iMM_recap_table.mat')
end

%% Preparation for LumpGEM 

% Rename the biomass exchange reaction to DM
biomassRxnIdx = find(ismember(model.rxns, "EXC_BOTH_Biomass_c"));
model.rxns(biomassRxnIdx) = {'DM_Biomass_c'};
model.varNames(ismember(model.varNames, "F_EXC_BOTH_Biomass_c")) = {'F_DM_Biomass_c'};
model.varNames(ismember(model.varNames, "R_EXC_BOTH_Biomass_c")) = {'R_DM_Biomass_c'};
model.varNames(ismember(model.varNames, "NF_EXC_BOTH_Biomass_c")) = {'NF_DM_Biomass_c'};

if strcmp(medium_type,'prototrophic')

    % Prepare the model for lumping
    [~, ~, EXC_rxns, ~] = putDrainsForward(model);
    EXC_rxns(contains(lower(EXC_rxns), 'biomass')) = [];
    biomass_rxns = model.rxns(contains(lower(model.rxns), 'biomass'));
    sol = solveTFAmodel_selections(model);
    
    % Loop through the different biomass growth rate
    rxns = biomass_rxns(1);
    uptakes_union = strrep(uptakes_union, 'R_', '');
    [~, uptakes_union_ind] = ismember(uptakes_union, model.rxns);
    uptakes_block = setdiff(EXC_rxns, uptakes_union);
    [~, uptakes_block_ind] = ismember(uptakes_block, model.rxns);
    rxns = [rxns; uptakes_union; uptakes_block];
    
    lbs = [sol.val * 0.99; repmat(-50, length(uptakes_union), 1); zeros(length(uptakes_block), 1)];
    ubs = [sol.val; model.ub(uptakes_union_ind); model.ub(uptakes_block_ind)];
    check_table = table(rxns,lbs,ubs);
    disp(size(check_table,1))
    model = update_bounds(model, rxns, lbs, ubs, 1, 1);
    
    % Change the objective function
    model = set_objective(model, biomass_rxns(2), 1,"max", 1);
    
    % Test ability of model to grow
    [~, ~] = solveAndCheckGrowthFBAandTFAmodel(model,biomass_rxns(2));
    
    num_drains_active = length(uptakes_union);
        
    % Lumping
    [activeRxns, LumpedRxnFormulas, bbbNames, DPsAll, IdNCNTNER, relaxedDGoVarsValues_ForEveryLumpedRxn] = ...
        lumpGEM_as_prototrophic(model, 1:length(model.rxns), DB_AlbertyUpdate, [], ...
                   'mmol', 'Smin', 'LCSBDefault', ...
                    [],'./MIPstart/lumpGEM_as_prototrophic.mst','./workspaces/step2_lumping_subnetworks_generated_as_prototrophic.mat',max_growth,num_drains_active, ...
                    "./workspaces/step2_lumping_ongoing_as_prototrophic.mat");
    
    % Save
    save('./workspaces/step2_lumping_done_as_prototrophic.mat');

elseif strcmp(medium_type,'auxotrophic')

    % Prepare the model for lumping
    biomass_rxns = model.rxns(contains(lower(model.rxns), 'biomass'));
    sol = solveTFAmodel_selections(model);
    
    % Impose minimal biomass growth rate
    rxns = biomass_rxns(1);
    lbs = sol.val * 0.99;
    ubs = sol.val;
    model = update_bounds(model, rxns, lbs, ubs, 1, 1);
    
    % Change the objective function
    model = set_objective(model, biomass_rxns(2), 1,"max", 1);
    
    % Test ability of model to grow
    [~, ~] = solveAndCheckGrowthFBAandTFAmodel(model,biomass_rxns(2));
    
    % Lumping
    [activeRxns, LumpedRxnFormulas, bbbNames, DPsAll, IdNCNTNER, relaxedDGoVarsValues_ForEveryLumpedRxn] = ...
        lumpGEM_as_auxotrophic(model, 1:length(model.rxns), DB_AlbertyUpdate, [], ...
                   'mmol', 'Smin', 'LCSBDefault', ...
                    [],'./MIPstart/lumpGEM_as_auxotrophic.mst','./workspaces/step2_lumping_subnetworks_generated_as_auxotrophic.mat',max_growth, ...
                    "./workspaces/step2_lumping_ongoing_as_auxotrophic.mat");

    % Save 
    save('./workspaces/step2_lumping_done_as_auxotrophic.mat');

end

%% Lumped reactions analysis

% Model
model_original = model;

% Identify the unique reactions
LumpRxns_unique = unique(LumpedRxnFormulas);

% Store the uptakes and secretions
Uptakes     = table('VariableTypes',["cell","cell"],'VariableNames',["id","coeff"],'Size',[size(LumpRxns_unique,1),2]);
Secretions  = table('VariableTypes',["cell","cell"],'VariableNames',["id","coeff"],'Size',[size(LumpRxns_unique,1),2]);
    
% Prepare the variables
Uptakes_union       = {};
Uptakes_inter       = {};
Secretions_union    = {};
Secretions_inter    = {};

% Index biomass parasite
[~,ind_BM_para]    = ismember('Biomass_c',model.mets);

% Restart the model 
model = model_original;

% Loop through the lumped rxns
for j = 1:size(LumpRxns_unique,1)

    model           = addReaction(model,strcat('parasitosome_N_',num2str(j)),LumpRxns_unique{j});
    [~,ind_para]    = ismember(strcat('parasitosome_N_',num2str(j)),model.rxns);

    ind_mets            = find(model.S(:,ind_para) < -1e-9);
    [~,ind_rxns]         = ismember(strcat('EXC_BOTH_',model.mets(ind_mets)),model.rxns);
    Uptakes.id{j}       = model.rxns(ind_rxns);
    Uptakes.coeff{j}    = full(model.S(ind_mets,ind_para))/full(model.S(ind_BM_para,ind_para));

    ind_mets            = find(model.S(:,ind_para) > 1e-9);
    [~,ind_rxns]         = ismember(strcat('EXC_BOTH_',model.mets(ind_mets)),model.rxns);
    ind_mets(ind_mets==find(contains(model.mets,'biomass_c','IgnoreCase',true))) = [];
    ind_rxns(ind_rxns==0) = [];
    Secretions.id{j}    = model.rxns(ind_rxns);
    Secretions.coeff{j} = full(model.S(ind_mets,ind_para))/full(model.S(ind_BM_para,ind_para));

    if j == 1
        Uptakes_union       = Uptakes.id{j};
        Uptakes_inter       = Uptakes.id{j};
        Secretions_union    = Secretions.id{j};
        Secretions_inter    = Secretions.id{j};
    else
        Uptakes_union       = union(Uptakes_union,Uptakes.id{j});
        Uptakes_inter       = intersect(Uptakes_inter,Uptakes.id{j});
        Secretions_union    = union(Secretions_union,Secretions.id{j});
        Secretions_inter    = intersect(Secretions_inter,Secretions.id{j});     
    end
 
end

Lump_exchanges_table = table('VariableTypes',repmat("cell",1,6),'VariableNames',["Uptakes","Uptakes_union","Uptakes_inter","Secretions","Secretions_union","Secretions_inter"],'Size',[1,6]);

Lump_exchanges_table.Uptakes            = {Uptakes};
Lump_exchanges_table.Uptakes_union      = {Uptakes_union};
Lump_exchanges_table.Uptakes_inter      = {Uptakes_inter};
Lump_exchanges_table.Secretions         = {Secretions};
Lump_exchanges_table.Secretions_union   = {Secretions_union};
Lump_exchanges_table.Secretions_inter   = {Secretions_inter};

%% Save the final workspace

save(strcat('./workspaces/step2_final_as_',medium_type,'.mat'));

end
