function step4_gene_essentiality_analyses(medium_type)

%% Prepare the environment

% Load the integrated host 
host = load(strcat('./models/integrated_hepatocyte_model_for_para_as_',medium_type,'.mat'),'host').host;

% Load the parasitosomes reactions
para_rxns = load(strcat('./workspaces/step2_final_as_',medium_type,'.mat'),'LumpedRxnFormulas').LumpedRxnFormulas;

% Identify the indices of the parasitosomes reactions
ind_parasitosomes_F = find(startsWith(host.varNames,'F_parasitosome'));

%% Update the gene-related field in the host

% find indices 
row2add = length(host.rxns) - size(host.rxnGeneMat,1);
ind_new = size(host.rxnGeneMat,1)+1:size(host.rxnGeneMat,1)+row2add;

% complete the matrix
host.rxnGeneMat(ind_new,:) = zeros(row2add,size(host.rxnGeneMat,2));
host.grRules(ind_new) = {''};
host.rules(ind_new) = {''};

% Adjust the csense
host.csense = repmat('E',length(host.mets),1);


%% Find the relaxation necesarry for each parasite to grow to its maximal biomass grwoth rate (for each parasite)

% Save the original model
host_original = host;

% Number of parasitosomes
num_parasitosomes = size(para_rxns,1);

% Close the parasitosome
host.var_lb(ind_parasitosomes_F) = 0;
host.var_ub(ind_parasitosomes_F) = 0;

% Get the drains and their bounds
[drains_F, drains_R, drains_NF, drains_F_lb, drains_F_ub, drains_R_lb, drains_R_ub, drains_NF_lb, drains_NF_ub, drains_F_ind, drains_R_ind, drains_NF_ind] = get_drain_bounds(host);

% Filter the drains to keep only the ones which can carry flux
ind_F = find(drains_F_ub>0);
ind_R = find(drains_R_ub>0);

% Create constraint cell
variableNames = [drains_F(ind_F); drains_R(ind_R)];

% Prepare variables to store the results

% Save the original medium 
drainNames = [drains_F;drains_R;drains_NF];
drainLBs = [drains_F_lb;drains_R_lb;drains_NF_lb];
drainUBs = [drains_F_ub;drains_R_ub;drains_NF_ub];
orig_medium = table(drainNames, drainLBs, drainUBs, 'VariableNames', {'Drain', 'LowerBound', 'UpperBound'});

% Creating the table for the results
gene10s = cell(num_parasitosomes,0);
updated_media = cell(num_parasitosomes,0);
orig_media = repmat({orig_medium}, num_parasitosomes, 1);
results_table = table(orig_media, updated_media, gene10s, 'VariableNames', {'OriginalMedium', 'UpdatedMedium', 'Genes_10PercentThreshold'});

% Load the gene database
gene_db = load('./data/humanGenome_hgnc_entrez.mat').humanGenome_hgnc_entrez;

% Run gene essentiality for the liver model

% Conduct gene essentiality 10% threshold for the given parasitosome
[growth_boolean_10_healthy,~,~] = geneEssentialityThermoInfectedHost(host, [], 0, 0.1, 0, [],0);

% Save the results in table
genes_ENTREZ = host.genes(find(growth_boolean_10_healthy<0.5));
genes_HGNC = cellfun(@(x) string(x),convertGeneNames(genes_ENTREZ,gene_db),'UniformOutput',false);    
genes_10_table_healthy = table(genes_ENTREZ,genes_HGNC,'VariableNames',["genes_ENTREZ","genes_HGNC"]);

save(strcat('./workspaces/step4_gene_essentiality_healthy_10_percent_for_para_as_',medium_type,'_done.mat'));

%% Relax the bounds in the model and conduct gene essentiality with active parasites

% Loop through the parasitosome
for i = 1:size(para_rxns,1)

    % Select a fresh model 
    host = host_original;

    % Open one parasitosome
    host.var_ub(ind_parasitosomes_F(i)) = 50;

    % Impose max growth of parasite and host
    host.var_lb(ind_biomass_NF) = 0.99*max_growth_host;
    host.var_lb(ind_para_NF)    = 0.99*max_growth_para;  
    
    % Find how much to relax
    [model_Slack,sol_Slack, slackValues, new_bounds] = MaxSlackVariables2ReduceBounds(host,variableNames,10*60);
    
    % Find the indices of the row for which the bounds are different
    ind_diff = find(new_bounds.Old_ub ~= new_bounds.New_ub);
    drain_diff = new_bounds.Reaction(ind_diff);

    % Reset the bounds
    host.var_lb(ind_biomass_NF) = 0;
    host.var_lb(ind_para_NF)    = 0;  

    % Save the updated medium in the results table
    updated_medium = orig_medium;
    [~,ind_medium] = ismember(drain_diff,updated_medium.Drain);
    updated_medium.UpperBound(ind_medium) = new_bounds.New_ub(ind_diff);
    results_table.UpdatedMedium{i} = updated_medium;

    % Adjust the medium in the host
    [~,ind_host] = ismember(drain_diff,host.varNames);
    host.var_ub(ind_host) = new_bounds.New_ub(ind_diff);

    % Conduct gene essentiality 10% threshold for the given parasitosome
    [growth_boolean_10,~,~] = geneEssentialityThermoInfectedHost(host, [], 1, 0.1, 0.1, [],0);

    % Save the results in table
    genes_ENTREZ = host.genes(find(growth_boolean_10<0.5));
    genes_HGNC = cellfun(@(x) string(x),convertGeneNames(genes_ENTREZ,gene_db),'UniformOutput',false);    
    genes_10_table = table(genes_ENTREZ,genes_HGNC,'VariableNames',["genes_ENTREZ","genes_HGNC"]);
    results_table.Genes_10PercentThreshold{i} = genes_10_table;

    save(strcat('./workspaces/step4_gene_essentiality_infected_10_percent_for_para_as_',medium_type,'_in_progress.mat'));

end

save(strcat('./workspaces/step4_gene_essentiality_infected_10_percent_for_para_as_',medium_type,'_done.mat'));


%% Union and Intersection

% Find the union
union_essential_genes_host_para_10 = vertcat(results_table.Genes_10PercentThreshold{:});
union_essential_genes_host_para_ENTREZ_10 = unique(union_essential_genes_host_para_10.genes_ENTREZ);
union_essential_genes_host_para_HGNC_10 = unique(cellfun(@(x) char(x), union_essential_genes_host_para_10.genes_HGNC, 'UniformOutput',false));

% Find the intersection
inter_essential_genes_host_para_ENTREZ_10 = results_table.Genes_10PercentThreshold{1}.genes_ENTREZ;
for i = 2:size(results_table,1)
    inter_essential_genes_host_para_ENTREZ_10 = intersect(inter_essential_genes_host_para_ENTREZ_10, results_table.Genes_10PercentThreshold{i}.genes_ENTREZ);
end
inter_essential_genes_host_para_HGNC_10 = cellfun(@(x) char(x),convertGeneNames(inter_essential_genes_host_para_ENTREZ_10,gene_db),'UniformOutput',false);  



%% Identify essential genes exclusively for the parasite

for i = 1:size(results_table,1)

    % Find the essential genes
    genes_ENTREZ    = setdiff(results_table.Genes_10PercentThreshold{i}.genes_ENTREZ ,genes_10_table_healthy.genes_ENTREZ);
    genes_HGNC = cellfun(@(x) string(x),convertGeneNames(genes_ENTREZ,gene_db),'UniformOutput',false);    
    genes_10_ess_table = table(genes_ENTREZ,genes_HGNC,'VariableNames',["genes_ENTREZ","genes_HGNC"]);

    % Fill the result table
    results_table.Genes_10Essential_Para{i} = genes_10_ess_table;

end

% Find the union
union_essential_genes_only_para_10 = vertcat(results_table.Genes_10Essential_Para{:});
union_essential_genes_only_para_ENTREZ_10 = unique(union_essential_genes_only_para_10.genes_ENTREZ);
union_essential_genes_only_para_HGNC_10 = unique(cellfun(@(x) char(x), union_essential_genes_only_para_10.genes_HGNC, 'UniformOutput',false));

% Find the intersection
inter_essential_genes_only_para_ENTREZ_10 = results_table.Genes_10Essential_Para{1}.genes_ENTREZ;
for i = 2:size(results_table,1)
    inter_essential_genes_only_para_ENTREZ_10 = intersect(inter_essential_genes_only_para_ENTREZ_10, results_table.Genes_10PercentThreshold{i}.genes_ENTREZ);
end
inter_essential_genes_only_para_HGNC_10 = cellfun(@(x) char(x),convertGeneNames(inter_essential_genes_only_para_ENTREZ_10,gene_db),'UniformOutput',false);  


%% Save the final results

% Save
save(strcat('./workspaces/step4_final_as_',medium_type,'.mat'));