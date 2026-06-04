function rxns_w_GPR_expressed = HPA_processing(model)

	%% Load liver transcriptomics data from Human Protein Atlas
	HPA = load('./data/transcriptomics_hepatocyte.mat');
	HPA = HPA.transcriptomics;

	%% Match genes from HPA to the genes in the model

	% Map HPA genes to model genes
	[genes_HPA_in_model, ind_model, ind_HPA] = intersect(model.genesHGNC, HPA.Gene_name);

	% Check for duplicates
	genes_model_w_data = model.genesHGNC(contains(model.genesHGNC, genes_HPA_in_model));

	% Check distribution of expressed genes in the model
	[~, ind_genes_HPA_in_model] = ismember(genes_HPA_in_model, HPA.Gene_name);

	% Iterate through the genes in the Human Protein Atlas dataset
	for i = 1:length(genes_HPA_in_model)
	    % Prepare the main field name
	    gene_name = string(genes_HPA_in_model(i));
	    gene = regexprep(gene_name, '-', '_'); % Replace '-' with '_'

	    % Add Ensembl ID
	    [~, ind_Ensembl] = ismember(HPA.Gene_name, gene);
	    ind_Ensembl = find(ind_Ensembl);

	    % Store names, indices, and NX values for the corresponding Ensembl ID
	    genesTable.(gene).Ensembl.name = string(HPA.Gene_name(ind_Ensembl));
	    genesTable.(gene).Ensembl.ind = ind_Ensembl;
	    genesTable.(gene).Ensembl.NX = HPA.NX(ind_Ensembl);

	    % Deduce the expression level (conservative assumption)
	    genesTable.(gene).NX_max = max(genesTable.(gene).Ensembl.NX);
	    genesTable.(gene).NX_min = min(genesTable.(gene).Ensembl.NX);
	    genesTable.(gene).NX_mean = mean(genesTable.(gene).Ensembl.NX);

	    % Add NCBI Entrez ID
	    [~, ind_NCBI] = ismember(model.genesHGNC, gene);
	    ind_NCBI = find(ind_NCBI);

	    % Store names, indices, and NX values for the corresponding NCBI ID
	    genesTable.(gene).NCBI.name = string(model.genes(ind_NCBI));
	    genesTable.(gene).NCBI.ind = ind_NCBI;
	end

	% Match the genes from HPA to the reactions in the model

	% Define parameters for active gene selection
	detection_threshold = 1;

	% Match genes and find expressed and non-expressed reactions with GPR
	varNames = ["Gene_Entrez_id", "Gene_HGNC_id", "mean_NX", "min_NX", "max_NX"];
	varTypes = ["string", "string", "double", "double", "double"];
	sz = [length(model.genes), length(varNames)];
	table_genes2rxns = table('VariableNames', varNames, 'VariableTypes', varTypes, 'Size', sz);

	% Iterate through the genes in the model
	for i = 1:length(model.genes)
	    % Fill the table with gene IDs
	    table_genes2rxns.Gene_Entrez_id(i) = model.genes(i);
	    table_genes2rxns.Gene_HGNC_id(i) = model.genesHGNC(i);

	    % Get the name of the gene of interest
	    gene = string(model.genesHGNC(i));

	    % Test if the gene is present in the data
	    if isfield(genesTable, gene)
		% Fill the table with appropriate values
		table_genes2rxns.mean_NX(i) = genesTable.(gene).NX_mean;
		table_genes2rxns.min_NX(i) = genesTable.(gene).NX_min;
		table_genes2rxns.max_NX(i) = genesTable.(gene).NX_max;
	    else
		% Fill the table with NaN values
		table_genes2rxns.mean_NX(i) = NaN;
		table_genes2rxns.min_NX(i) = NaN;
		table_genes2rxns.max_NX(i) = NaN;
	    end
	end

	% Map genes to reactions
	genes = model.genes;
	geneRatio = table_genes2rxns.max_NX;
	f_and = @min;
	f_or = @max;
	[rxns_w_GPR, rxnNXs] = evaluateGPR(model, genes, geneRatio, f_and, f_or);

	% Evaluate whether a reaction is active based on expression data
	rxns_w_GPR_null = rxns_w_GPR(rxnNXs < detection_threshold)';
	rxns_w_GPR_no_data = rxns_w_GPR(isnan(rxnNXs))';
	rxns_w_GPR_expressed = rxns_w_GPR(rxnNXs >= detection_threshold)';
	rxns_no_GPR = setdiff(model.rxns, rxns_w_GPR);

end
