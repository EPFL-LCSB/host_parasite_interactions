function [producable,uptakes_pb,secretable,secretions_pb] = testForParaMetsInHost(model,matchesTable,uptakes_KEGG,uptakes_KEGG_missing,secretions_KEGG, secretions_KEGG_missing,threshold)

%%%%%%%%%%%%%%%%%%%%%%%%
% Uptakes of parasites

% Array to store the results
producable = nan(size(uptakes_KEGG));

for i = 1:size(producable,1)

    % Identify the met to produce
    met = matchesTable.Source(find(contains(matchesTable.KEGG_ID,uptakes_KEGG(i))));
    met = met{1}{1};

    % Add a drain reaction
    rxn = strcat('DM_',met);
    formula = strcat(met,' <=> ');
    [model_i,rxnIDexists] = addReaction(model,rxn,formula);
    if rxnIDexists
        ind = rxnIDexists;
    else
        ind = length(model_i.rxns);
    end
    model_i.lb(ind) = 1e-6;
    model_i.ub(ind) = 50;

    % Impose secretion to force production
    model_i.c = zeros(length(model_i.rxns),1);
    model_i.c(ind) = 1;

    sol = solveFBAmodel_selections(model_i);

    if ~isnan(sol.f)
        producable(i) = sol.f;
    end

end

% Find potential problem 
uptakes_pb = uptakes_KEGG(producable<=threshold);
uptakes_pb = setdiff(uptakes_pb,uptakes_KEGG_missing);

%%%%%%%%%%%%%%%%%%%%%%%%
% Secretions of parasites

% Array to store the results
secretable = nan(size(secretions_KEGG));

for i = 1:size(secretable,1)   

    % Identify the met to produce
    met = matchesTable.Source(find(contains(matchesTable.KEGG_ID,secretions_KEGG(i))));
    met = met{1}{1};

    % Add a drain reaction
    rxn = strcat('DM_',met);
    formula = strcat(met,' <=> ');
    [model_i,rxnIDexists] = addReaction(model,rxn,formula);
    if rxnIDexists
        ind = rxnIDexists;
    else
        ind = length(model_i.rxns);
    end    
    model_i.lb(ind) = -50;
    model_i.ub(ind) = -1e-6;

    % Impose secretion to force production
    model_i.c = zeros(length(model_i.rxns),1);
    model_i.c(ind) = -1;

    sol = solveFBAmodel_selections(model_i);

    if ~isnan(sol.f)
        secretable(i) = sol.f;
    end

end

% Should exclude from the problematic reactions the reactions that were
% added as specific in DM

% Find potential problem 
secretions_pb = secretions_KEGG(secretable<=threshold);
secretions_pb = setdiff(secretions_pb,secretions_KEGG_missing);
