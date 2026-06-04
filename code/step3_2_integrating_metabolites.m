function step3_2_integrating_metabolites(medium_type)

    %% Loading the data and models
    
    % load ifpa
    load('./models/iPfa.mat');
        
    % load model
    load('./data/thermoData_Recon3.mat');
    load('./data/PmetsRecon3_redHUMAN.mat');
    host_red = load('./models/hepatocyte_thermocurated.mat').model;
    host_full = load('./models/GEM_Recon3_thermocurated_redHUMAN.mat').model;

    % Load the appropriate data
    matchesTable = load('./data/matches_table_para_host.mat').matches; 
    load(strcat('./workspaces/step2_final_as_',medium_type,'.mat'),'Lump_exchanges_table');
    
    % Set the bigM for the thermo 
    bigM = 1e3;
    bigM_thermo = 5e4;
    
    %% Find the missing metabolites based on KEGG ids

    % Extract the KEGG ids
    uptakes_KEGG    = extractBetween(Lump_exchanges_table.Uptakes_union{1},'EXC_BOTH_','_e');
    secretions_KEGG = extractBetween(Lump_exchanges_table.Secretions_union{1},'EXC_BOTH_','_e');
    
    % Check the coverages from host_full
    uptakes_KEGG_missing    = uptakes_KEGG(ismember(uptakes_KEGG,matchesTable.KEGG_ID(cellfun(@(x) isempty(x), matchesTable.Source))));
    secretions_KEGG_missing = secretions_KEGG(ismember(secretions_KEGG,matchesTable.KEGG_ID(cellfun(@(x) isempty(x), matchesTable.Source))));
    
    % Get the names 
    uptakes_KEGG_missing_names      = queryKEGGIDS(uptakes_KEGG_missing);
    secretions_KEGG_missing_names   = queryKEGGIDS(secretions_KEGG_missing);
    
    % Display the names 
    disp(uptakes_KEGG_missing_names);
    disp(secretions_KEGG_missing_names);
    
    %% Loop through the missing metabolites and handle the special cases
    
    % host_full_m stands for host full modified
    host_full_m = host_full;

    % Save the name of the txt file containing the biGG database
    bigg_metabolites_file = './data/bigg_models_metabolites.txt'; % Replace this with the correct file path

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Uptakes: add drains or supplementary reactions for the nutrients supposed
    % to be taken up by the parasites

    for i = 1:length(uptakes_KEGG_missing)

        % Get the KEGG id for the metabolite
        metKEGGID = uptakes_KEGG_missing{i};

        % Special cases

        % Special case #1: Methylcobalamin
        if strcmp(metKEGGID,'C06453')

            % Update the matches tables
            matchesTable.Source{strcmp(matchesTable.KEGG_ID,metKEGGID)} = {'cbl2_c'};
            matchesTable.Is_Matched(strcmp(matchesTable.KEGG_ID,metKEGGID)) = 1;
    
            % Test for feasibility
            solveFBAmodelCplex(host_full_m)

        % Special case #2: Sulfur donor   
        elseif strcmp(metKEGGID,'C17023')
    
            % We chose to select L-methionine to act as a sulfur donor
            host_full_m = addReaction(host_full_m,'conv_sulfur_donor','met_L_c <=> sulfur_donor_c ');
            host_full_m.lb(end) = 0;
            host_full_m.ub(end) = 50;
    
            % Update the matches tables
            matchesTable.Source{strcmp(matchesTable.KEGG_ID,metKEGGID)} = {'sulfur_donor_c'};
            matchesTable.Is_Matched(strcmp(matchesTable.KEGG_ID,metKEGGID)) = 1;
            
            % Test for feasibility
            solveFBAmodelCplex(host_full_m)

        % Special case #3: Apoprotein
        elseif strcmp(metKEGGID,'C16240')

            metsN6 = [{'M00213_c'};{'Ndmelys_c'};{'Ntmelys_c'}];
            for j = 1:length(metsN6)
                host_full_m = addReaction(host_full_m,strcat(strrep(metsN6{j},'_c',''),'_to_Nnmelys'),strcat(metsN6{j},' <=> Nnmelys_c '));
                host_full_m.lb(end) = 0;
                host_full_m.ub(end) = 50;
            end

            % Update the matches tables
            matchesTable.Source{strcmp(matchesTable.KEGG_ID,metKEGGID)} = {'Nnmelys_c'};
            matchesTable.Is_Matched(strcmp(matchesTable.KEGG_ID,metKEGGID)) = 1;

            % Test for feasibility            
            solveFBAmodelCplex(host_full_m)

        elseif strcmp(metKEGGID,'H_ETC')
            host_full_m = addReaction(host_full_m,'conv_h_H_passed','H_ETC <=> h_c ');
            host_full_m.lb(end) = -50;
            host_full_m.ub(end) = 50;

            matchesTable.Source{strcmp(matchesTable.KEGG_ID,metKEGGID)}  = {'H_ETC'};
            matchesTable.Is_Matched(strcmp(matchesTable.KEGG_ID,metKEGGID)) = 1;

        % Generic case
        else
            % Query the bigg database for a good name as a metabolite
            bigg_metabolites = find_bigg_metabolites_by_kegg_id(metKEGGID, bigg_metabolites_file);
            met = strcat(bigg_metabolites{1},'_c');
            
            % Add demand reactions
            host_full_m = addReaction(host_full_m,strcat('DM_',met),strcat(met,' <=> '));
            host_full_m.lb(end) = -50;
            host_full_m.ub(end) = 0;
            
            % Update the matches tables
            matchesTable.Source{strcmp(matchesTable.KEGG_ID,metKEGGID)} = {met};
            matchesTable.Is_Matched(strcmp(matchesTable.KEGG_ID,metKEGGID)) = 1;
            
            solveFBAmodel_selections(host_full_m)            

        end

    end    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Secretions

     for i = 1:length(secretions_KEGG_missing)

        % Get the KEGG id for the metabolite
        metKEGGID = secretions_KEGG_missing{i};

        % Special cases

        % Special case #1: H_plus
        if strcmp(metKEGGID,'C00080')

            % Add demand reactions    
            host_full_m = addReaction(host_full_m,'conv_h_H_plus','h_c <=> H_plus_c ');
            host_full_m.lb(end) = -50;
            host_full_m.ub(end) = 50;
        
            % Fill in the tables
            matchesTable.Source{strcmp(matchesTable.KEGG_ID,metKEGGID)}  = {'H_plus_c'};
            matchesTable.Is_Matched(strcmp(matchesTable.KEGG_ID,metKEGGID)) = 1;
    
        % Generic case
        else
            % Query the bigg database for a good name as a metabolite
            bigg_metabolites = find_bigg_metabolites_by_kegg_id(metKEGGID, bigg_metabolites_file);
            met = strcat(bigg_metabolites{1},'_c');
            
            % Add demand reactions
            host_full_m = addReaction(host_full_m,strcat('DM_',met),strcat(met,' <=> '));
            host_full_m.lb(end) = 0;
            host_full_m.ub(end) = 50;
            
            % Update the matches tables
            matchesTable.Source{strcmp(matchesTable.KEGG_ID,metKEGGID)} = {met};
            matchesTable.Is_Matched(strcmp(matchesTable.KEGG_ID,metKEGGID)) = 1;
            
            % Test for feasibility            
            solveFBAmodelCplex(host_full_m)            

        end

    end    
 
    %% Sanity check 
    
    [producable,uptakes_pb,secretable,secretions_pb] = testForParaMetsInHost(host_full_m,matchesTable,uptakes_KEGG,uptakes_KEGG_missing,secretions_KEGG, secretions_KEGG_missing,10);
     
    % Add uptakes
    for i = 1:length(uptakes_pb)    
        met = matchesTable.Source{contains(matchesTable.KEGG_ID,secretions_pb{i})}{1};
        rxn = strcat('DM_',met);
        formula = strcat(met,' <=> ');
        [host_full_m,~] = addReaction(host_full_m,rxn,formula);
        host_full_m.lb(end) = -50;
        host_full_m.ub(end) = 0;    
    end

    % Add secretions
    for i = 1:length(secretions_pb)    
        met = matchesTable.Source{contains(matchesTable.KEGG_ID,secretions_pb{i})}{1};
        rxn = strcat('DM_',met);
        formula = strcat(met,' <=> ');
        [host_full_m,~] = addReaction(host_full_m,rxn,formula);
        host_full_m.lb(end) = 0;
        host_full_m.ub(end) = 50;    
    end    
    
    %% PRELIMINARY STEP 1: Add the required demand reactions before the minimization step
    
    host_full_m_w_dmd_test = host_full_m;
    
    for i = 1:size(uptakes_KEGG,1)
    
        % Identify the met to produce
        met = matchesTable.Source(find(contains(matchesTable.KEGG_ID,uptakes_KEGG(i))));
        met = met{1}{1};
        met_test = strcat(met,'_test');
    
        % Add conversion reaction
        rxn_conversion = strcat(met,'_conv_',met_test);
        formula_conversion = [met, ' <=> ', met_test];
        [host_full_m_w_dmd_test,~] = addReaction(host_full_m_w_dmd_test,rxn_conversion,formula_conversion);
        host_full_m_w_dmd_test.lb(end) = -50;
        host_full_m_w_dmd_test.ub(end) = 50;
    
        % Add drain reaction
        rxn_drain = strcat('DM_',met_test);
        formula_drain = strcat(met_test,' <=> ');
        [host_full_m_w_dmd_test,~] = addReaction(host_full_m_w_dmd_test,rxn_drain,formula_drain);   
        host_full_m_w_dmd_test.lb(end) = 1e-6;
        host_full_m_w_dmd_test.ub(end) = 50;
    
    end
     
    %%%%%%%%%%%%%%%%%%%%%%%%
    % Secretions of parasites
        
    for i = 1:size(secretions_KEGG,1)   
    
        % Identify the met to produce
        met = matchesTable.Source(find(contains(matchesTable.KEGG_ID,secretions_KEGG(i))));
        met = met{1}{1};
    
        % Add conversion reaction
        rxn_conversion = strcat(met,'_conv_',met_test);
        formula_conversion = [met, ' <=> ', met_test];
        [host_full_m_w_dmd_test,~] = addReaction(host_full_m_w_dmd_test,rxn_conversion,formula_conversion);
        host_full_m_w_dmd_test.lb(end) = -50;
        host_full_m_w_dmd_test.ub(end) = 50;
    
        % Add drain reaction
        rxn_drain = strcat('DM_',met_test);
        formula_drain = strcat(met_test,' <=> ');
        [host_full_m_w_dmd_test,~] = addReaction(host_full_m_w_dmd_test,rxn_drain,formula_drain);   
        host_full_m_w_dmd_test.lb(end) = -50;
        host_full_m_w_dmd_test.ub(end) = -1e-6;
    
    end
    
    solveFBAmodelCplex(host_full_m_w_dmd_test)
    
    %% PRELIMINARY STEP 2: Convert full host model to thermo and check for growth before minimization
    
    host_full_m_w_dmd_test = FBA2TFA(host_full_m_w_dmd_test,thermoData_Recon3,bigM,bigM_thermo);
    
    solveAndCheckGrowthFBAandTFAmodel(host_full_m_w_dmd_test);
    
    save(strcat('./workspaces/step3_2_host_full_thermo_for_para_as_',medium_type,'.mat'))
    
    %% Add the minimal number of rxns to the context-specific host model
    
    % Preparation for maximization of inactive rxns
    uptakes = cellfun(@(x) x{1},matchesTable.Source(find(contains(matchesTable.KEGG_ID,uptakes_KEGG))),'UniformOutput',false);
    secretions = cellfun(@(x) x{1},matchesTable.Source(find(contains(matchesTable.KEGG_ID,secretions_KEGG))),'UniformOutput',false);   
    dir_path = strcat('./workspaces/maxInactiveRxns_on_',medium_type,'_medium/');
    [model_prep_save_file,uptakes_save_file,secretions_save_file] = deal(strcat(dir_path,'prep_model_',medium_type,'.mat'),strcat(dir_path,'uptakes_rxns_',medium_type,'.mat'),strcat(dir_path,'secretions_rxns_',medium_type,'.mat'));
   
    % Max the number of inactive different reactions
    [rxns2add,uptakes_pb,secretions_pb] = maxInactiveRxns(host_red,host_full_m_w_dmd_test,15*60,model_prep_save_file,uptakes,uptakes_save_file,secretions,secretions_save_file,thermoData_Recon3);
        
    % Save 
    save(strcat('./workspaces/step3_2_rxns2add_thermo_analysis_for_para_as_',medium_type,'.mat'));
    
    %% Add the minimal number of rxns to the context-specific host model and test for the prodcution/secretions of metabolites of interest
    
    % Preserve original model 
    host_red_m = host_red;
    
    % Add the reactions
    host_red_m = addReactionFromGEM(host_red_m,host_full_m,rxns2add);
        
    % Test for growth
    solveFBAmodel_selections(host_red_m)
    
    % Test for the metabolites
    [producable,uptakes_pb,secretable,secretions_pb] = testForParaMetsInHost(host_red_m,matchesTable,uptakes_KEGG,uptakes_KEGG_missing,secretions_KEGG, secretions_KEGG_missing, 1e-6);
    
    % Add drains for problematic uptakes
    for i = 1:length(uptakes_pb)
    
        met = matchesTable.Source{contains(matchesTable.KEGG_ID,uptakes_pb{i})}{1};
        rxn = strcat('DM_',met);
        formula = strcat(met,' <=> ');
        [host_red_m,~] = addReaction(host_red_m,rxn,formula);
        host_red_m.lb(end) = -50;
        host_red_m.ub(end) = 0;
    
    end
    
    % Add drains for problematic secretions
    for i = 1:length(secretions_pb)
    
        met = matchesTable.Source{contains(matchesTable.KEGG_ID,secretions_pb{i})}{1};
        rxn = strcat('DM_',met);
        formula = strcat(met,' <=> ');
        [host_red_m,~] = addReaction(host_red_m,rxn,formula);
        host_red_m.lb(end) = 0;
        host_red_m.ub(end) = 50;
    
    end
  
    %% Convert reduced host model to thermo and check for growth
    
    % Convert the model to thermo
    host_red_m = FBA2TFA(host_red_m,thermoData_Recon3,bigM,bigM_thermo);
    
    % Check for the ability of the model to grow
    solveAndCheckGrowthFBAandTFAmodel(host_red_m);
    
    %% Save the workspace as well as the model 
    
    % Save the workspace
    save(strcat('./workspaces/step3_2_final_as_',medium_type,'.mat'));

    % Save the model
    save(strcat('./models/hepatocyte_with_drains_for_para_as_',medium_type,'.mat'),'host_red_m');

end