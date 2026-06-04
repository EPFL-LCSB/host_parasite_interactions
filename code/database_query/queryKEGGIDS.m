function compound_names = queryKEGGIDS(compound_ids)

% Initialize an empty cell array to store compound names
compound_names = cell(size(compound_ids));

% KEGG API base URL
kegg_api_url = 'http://rest.kegg.jp/get/';

% Loop through the compound_ids and query the KEGG database
for i = 1:length(compound_ids)
    compound_id = compound_ids{i};
    url = [kegg_api_url, compound_id];
    
    % Fetch compound information
    try
        compound_info = webread(url);
        
        % Extract compound names using regular expressions
        match = regexp(compound_info, 'NAME\s+((?:[^\r\n]+\r?\n?)+)', 'tokens');
        if ~isempty(match)
            compound_name_raw = match{1}{1};
            % Replace line breaks and whitespaces with semicolons and single spaces, respectively
            compound_name = regexprep(compound_name_raw, '\s*\n\s*', '; ');
            compound_name = regexprep(compound_name, ';$', '');
            if contains(compound_name,'FORMULA')
                compound_name = extractBefore(compound_name,'FORMULA');
            elseif contains(compound_name,'COMMENT')
                compound_name = extractBefore(compound_name,'COMMENT');
            end
            compound_names{i} = compound_name;
        else
            compound_names{i} = 'Not found';
        end
    catch
        if strcmp(compound_id,'H_ETC')
            compound_names{i} = 'H_ETC';
        else
            compound_names{i} = 'Error in fetching';
        end
    end
end

compound_names = strrep(compound_names,';;',';');
