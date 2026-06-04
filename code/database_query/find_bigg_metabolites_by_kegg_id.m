function [bigg_metabolites] = find_bigg_metabolites_by_kegg_id(kegg_id, bigg_metabolites_file)
    % Read the BiGG metabolites file
    bigg_data = readtable(bigg_metabolites_file, 'Delimiter', '\t', 'ReadVariableNames', true);

    % Initialize the list of BiGG metabolites
    bigg_metabolites = {};

    % Iterate through the rows and extract the BiGG metabolite IDs
    for i = 1:height(bigg_data)
        database_links = bigg_data{i, 'database_links'}{1};

        kegg_pattern = strcat('KEGG Compound: http://identifiers.org/kegg.compound/', kegg_id);

        if contains(database_links, kegg_pattern)
            bigg_metabolite_id = bigg_data{i, 'universal_bigg_id'}{1};
            bigg_metabolites = [bigg_metabolites; bigg_metabolite_id];
        end
    end
end
