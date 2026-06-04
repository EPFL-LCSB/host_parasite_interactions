function model = genes_HGNC(model,humanGenome_hgnc_entrez)

if nargin<2
    load('humanGenome_hgnc_entrez.mat')
end

genes_model = [];
    for this_gene = model.genes'    
        this_gene = strsplit(this_gene{1},'.');
        genes_model = [genes_model;this_gene(1)];
    end
    genes_recon3_unique = unique(genes_model);

    ncbi_ids = humanGenome_hgnc_entrez.entrez_id; %arrayfun(@num2str, humanGenome_hgnc_entrez.entrez_id, 'un', 0);
    table = [ncbi_ids humanGenome_hgnc_entrez.symbol];
    [~,bb] = ismember(genes_recon3_unique, table(:,1));

    table_genes_recon3 = [genes_recon3_unique(find(bb)) table(bb(find(bb)),:)];

    model.genesHGNC={};
    for this_gene = genes_model'
        id = find(ismember(table_genes_recon3(:,1),this_gene));
        if ~isempty(id)
            model.genesHGNC = [model.genesHGNC;table_genes_recon3(id,3)];
        elseif isempty(id)
            model.genesHGNC = [model.genesHGNC;{''}];
        end
    end
end