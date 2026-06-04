function step3_1_matching_metabolites()

%% Loading the data 

% load models
para = load('./models/iPfa.mat').ipfa4liver;
recon3d = load('./models/GEM_Recon3_thermocurated_redHUMAN.mat').model;

% load human metabolites list of equivalences
recon3dTable = readtable('./data/Recon3D.xls','Sheet','chemicals');

%% Get the unique exchanged metabolites KEGG id

exchange_rxns = startsWith(para.rxns, 'EXC_BOTH_') & endsWith(para.rxns, '_e');
exchange_mets = sum(para.S(:, exchange_rxns), 2) ~= 0;
rowNames = string(para.metNames(exchange_mets));
Size = [length(rowNames), 2];
ipfaTable = table('Size', Size, 'VariableNames', ["KEGG_ID", "SEED_ID"], 'VariableTypes', ["string", "string"], 'RowNames', rowNames);
ipfaTable.KEGG_ID = strrep(para.mets(exchange_mets), '_e', '');
ipfaTable.SEED_ID = strrep(para.metSEEDID(exchange_mets), '_e', '');

% Extract KEGG IDs, SEED IDs, Xref, and Source columns
kegg_ids = ipfaTable.KEGG_ID;
seed_ids = ipfaTable.SEED_ID;
xref_col_original = recon3dTable.Xref;
source_col = cellfun(@(src) strjoin(cellfun(@(s) extractAfter(s, 'M_'), strsplit(src, ';'), 'UniformOutput', false), ';'), recon3dTable.Source, 'UniformOutput', false);
source_col = cellfun(@(cell_str) strjoin(regexp(cell_str, '[\w]*_c', 'match'), ';'), source_col, 'UniformOutput', false);
ind_c = find(cellfun(@(x) ~isempty(x), source_col));
source_col = source_col(ind_c);
xref_col = xref_col_original(ind_c);

% Functions to find specific KEGG or SEED ID in the Xref column
find_match_KEGG_ID = @(id) find(cellfun(@(xref) contains(xref, ['keggC:' id]), xref_col), 1, 'first');
find_match_SEED_ID = @(id) find(cellfun(@(xref) contains(xref, ['seedM:' id]), xref_col), 1, 'first');

% Find matching KEGG IDs and SEED IDs
kegg_matches_idx = cellfun(find_match_KEGG_ID, kegg_ids, 'UniformOutput', false);
seed_matches_idx = cellfun(find_match_SEED_ID, seed_ids, 'UniformOutput', false);

% Combine KEGG and SEED matches
combined_matches_idx = cellfun(@union, kegg_matches_idx, seed_matches_idx, 'UniformOutput', false);

% Get the matched sources
matched_sources = cellfun(@(idx) source_col(idx), combined_matches_idx, 'UniformOutput', false);

% Create a table with KEGG IDs, SEED IDs, corresponding cleaned Sources, and Is_Matched column
matches = table(rowNames, kegg_ids, seed_ids, matched_sources, ~cellfun(@isempty, matched_sources), 'VariableNames', {'Names', 'KEGG_ID', 'SEED_ID', 'Source', 'Is_Matched'});

% Number of matched metabolites
disp(length(find(matches.Is_Matched)));

%% For the missing matches, try to operate based on name, and ask the user to accept the change or not

idx_missing = find(matches.Is_Matched == 0);
recon3dTable_cytosolic = recon3dTable(ind_c,:);


for i = 1:length(idx_missing)
    met = ipfaTable.Properties.RowNames{idx_missing(i)};
    idx_potential = find(contains(lower(recon3dTable_cytosolic.Name), lower(met)));

    if ~isempty(idx_potential)
        cell_array = strcat(recon3dTable_cytosolic.Name(idx_potential),' | ', recon3dTable_cytosolic.Source(idx_potential));
        fprintf('%s\n', repmat('#', 1, 20));
        fprintf('We found match(es) for %s:\n\n', met);

        for j = 1:length(cell_array)
            fprintf('%d. %s\n', j, cell_array{j});
        end

        user_choice = input('\nWhich option is the best? (answer 0 if none is suitable): ');

        if user_choice > 0
            % Fill the matches table with the value
            % Assuming both tables have the same variable names
            source = recon3dTable_cytosolic.Source(idx_potential(user_choice));
            source = cellfun(@(src) strjoin(cellfun(@(s) extractAfter(s, 'M_'), strsplit(src, ';'), 'UniformOutput', false), ';'), source, 'UniformOutput', false);
            source = cellfun(@(cell_str) strjoin(regexp(cell_str, '[\w]*_c', 'match'), ';'), source, 'UniformOutput', false);

            matches.Source(idx_missing(i)) = {source};
            matches.Is_Matched(idx_missing(i)) = true;

            disp('Value added to the matches table.');
        else
            disp('Value not added to the matches table.');
        end
    end
end

%% Initialize variables
idx_missing = find(matches.Is_Matched == 0);

% Iterate through the missing matches
for i = 1:length(idx_missing)
    met = ipfaTable.Properties.RowNames{idx_missing(i)};
    KEGG_ID = ipfaTable.KEGG_ID{idx_missing(i)};
    idx_potential_KEGG = find_met_id(recon3d.metKEGGID, KEGG_ID);
    SEED_ID = ipfaTable.SEED_ID{idx_missing(i)};
    idx_potential_SEED = find_met_id(recon3d.metSEEDID, SEED_ID);

    idx_potential = union(idx_potential_KEGG, idx_potential_SEED);

    if ~isempty(idx_potential)
        cell_array = strcat(recon3d.metNames(idx_potential), ' | ', recon3d.mets(idx_potential));
        idx_potential(~endsWith(cell_array, '_c')) = [];
        cell_array(~endsWith(cell_array, '_c')) = [];

        if ~isempty(cell_array)
            fprintf('%s\n', repmat('#', 1, 20));
            fprintf('We found match(es) for %s:\n\n', met);

            for j = 1:length(cell_array)
                fprintf('%d. %s\n', j, cell_array{j});
            end

            user_choice = input('\nWhich option is the best? (answer 0 if none is suitable): ');

            if user_choice > 0
                source = recon3d.mets(idx_potential(user_choice));
                matches.Source(idx_missing(i)) = {source};
                matches.Is_Matched(idx_missing(i)) = true;

                disp('Value added to the matches table.');
            else
                disp('Value not added to the matches table.');
            end
        end
    end
end

%% Save
save('./data/matches_table_para_host.mat','matches');

%% Save
save('./workspaces/step_3_1_final.mat');

end
