function idx = find_met_id(met_id_col, met_id)
    if ~strcmp(met_id, 'NA') && ~isempty(met_id)
        idx = find(contains(lower(met_id_col), lower(met_id)));
    else
        idx = [];
    end
end