function step3_3_integrating_parasitosomes(medium_type)

    %% Loading the data 
    
    % load parasite model
    load('./models/iPfa.mat');
    
    % load human model
    load('./data/thermoData_Recon3.mat');
    load('./data/PmetsRecon3_redHUMAN.mat');
    
    % Set the thermo values
    bigM = 1e3;
    bigM_thermo = 5e4;
    
    %% Prepare the environment
    
    % Load the host with required demand rxns for parasite growth
    host = load(strcat('./models/hepatocyte_with_drains_for_para_as_',medium_type,'.mat'),'host_red_m').host_red_m;

    para = load('./models/iPfa.mat').ipfa4liver;
    
    % Load the parasitosome reactions
    para_rxns = load(strcat('./workspaces/step2_final_as_',medium_type,'.mat'),'LumpedRxnFormulas').LumpedRxnFormulas;
    
    % Load the matches table
    matchesTable = load(strcat('./workspaces/step3_2_final_as_',medium_type,'.mat'),'matchesTable').matchesTable;
    
    %% Test growth of each parasitosome in the host
    
    % Description: Test minimal growth of each parasitosome across each yield cut
    
    % Add demand rxns
    biomass_para = 'biomass_para';
    host = addReaction(host,['DM_',biomass_para],[biomass_para,' <=> ']);
    host.lb(end) = 0;
    host.ub(end) = 50;
    ind_para = find(ismember(host.rxns,'DM_biomass_para'));
    ind_matches = size(matchesTable,1)+1;
    matchesTable.Names(ind_matches) = "Biomass_parasite";
    matchesTable.Source(ind_matches) = {{biomass_para}};
    matchesTable.Is_Matched(ind_matches) = 1;
    matchesTable.KEGG_ID(ind_matches) = {'Biomass_c'};
    
    % Create a table to store the results
    growth_rate = table('VariableTypes',["double","double"],'VariableNames',["parasitosome","DM_para"],'Size',[size(para_rxns,1),2]);
    
    % Loop through the parasitosome
    for i = 1:size(para_rxns,1)
    
        % Get the parasitosome reaction
        rxn = para_rxns{i};
        rxnName = strcat('parasitosome_N_',num2str(i));
    
       % Add the reaction
        host = addParaRxn(host,para,matchesTable,rxn,rxnName);
    
        % Test for non null growth (i.e., above 1e-9?)
        ind_parasitosome = length(host.rxns);
        host.c = zeros(length(host.rxns),1);
        host.c(ind_parasitosome) = 1;
    
        sol = solveFBAmodel_selections(host);
        disp(sol.f);
        
        if ~isnan(sol.f)
            growth_rate.parasitosome(i) = sol.f;
        else
            growth_rate.parasitosome(i) = nan;
        end
    
        host.c = zeros(length(host.rxns),1);
        host.c(ind_para) = 1;
    
        sol = solveFBAmodel_selections(host);
        disp(sol.f);
    
    
        if ~isnan(sol.f)
            growth_rate.DM_para(i) = sol.f;
        else
            growth_rate.DM_para(i) = nan;
        end
    
        host.lb(ind_parasitosome) = 0;
        host.ub(ind_parasitosome) = 0;
    
    end
    
    % Check if parasite cannot grow
    ind_para_no_growth = find(growth_rate.parasitosome == 0);
    
    if ind_para_no_growth
        error("parasite not growing");
    end
       
    %% Convert the host to thermo
    
    % Convert the host model to thermo 
    host = FBA2TFA(host,thermoData_Recon3,bigM,bigM_thermo);
    
    % Test for feasibility in reasonable time
    solveAndCheckGrowthFBAandTFAmodel(host);
    
    %% Save the integrated model and the workspace
    
    % Save the model
    save(strcat('./models/integrated_hepatocyte_model_for_para_as_',medium_type,'.mat'),'host');

    % Save the workspace    
    save(strcat('./workspaces/step3_3_host_with_thermo_for_para_as_',medium_type,'.mat'));
          
end