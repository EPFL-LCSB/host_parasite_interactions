function geneNames_HGNC = convertGeneNames(geneNames, geneMapping)

%% Description
%
% convertGeneNames converts a list of gene names from Entrez ID format to HGNC
% symbol format using a provided gene mapping dataset.
%
% Usage:
%   geneNames_HGNC = convertGeneNames(geneNames, geneMapping)
%
% Inputs:
%   geneNames      : A cell array of gene names in Entrez ID format.
%   geneMapping    : A table containing the mapping between Entrez IDs and
%                    HGNC symbols.
%
% Output:
%   geneNames_HGNC : A cell array of gene names in HGNC symbol format.
%
% Notes:
%   The function iterates through each gene in the input list, splitting
%   the Entrez ID from the version number and finding the corresponding
%   HGNC symbol in the geneMapping table.
%
% Example:
%   load('./data/humanGenome_hgnc_entrez.mat');
%   geneNames_Entrez = {'1.1', '2.1', '3.1'};
%   geneNames_HGNC = convertGeneNames(geneNames_Entrez, humanGenome_hgnc_entrez);
%
% This will convert the gene names from Entrez ID format to HGNC symbol format.


geneNames_HGNC = cell(size(geneNames));

for i = 1:length(geneNames)
    gene = geneNames{i};
    gene = strsplit(gene, '.');
    gene = gene{1};
    gene_HGNC = geneMapping.symbol(strcmp(geneMapping.entrez_id, gene));
    geneNames_HGNC{i} = gene_HGNC;
end

end