% OBTAIN GENES FROM ORIGINAL MODEL
function model = genes_fromGEM(model,gem)
    
    if nargin<2
        gem = model.OriginalGEM;
    end

    model.lb(find(model.c)) = 0;
%     model.var_lb(find(model.f)+1) = 0;

    [~,ind_rxns] = ismember(model.rxns, gem.rxns);
    ind_rxns = ind_rxns(find(ind_rxns)); % because of the lumps
    rxnGeneMatrix = gem.rxnGeneMat(ind_rxns,:);

    ind_genes = find(sum(rxnGeneMatrix,1));

    model.rxnGeneMat = rxnGeneMatrix(:,ind_genes);
    model.genes = gem.genes(ind_genes);
    model.rules = gem.rules(ind_rxns);
    model.grRules = gem.grRules(ind_rxns);

    rules = model.rules;

    for i = 1:length(ind_genes)
        rules = strrep(rules,['(',num2str(ind_genes(i)),')'],['(_',num2str(ind_genes(i)),'_)']);
    end

    for i = 1:length(ind_genes)
       rules = strrep(rules,['(_',num2str(ind_genes(i)),'_)'],['(',num2str(i),')']);
    end

    model.rules = rules;

    % include the lumps (zeros everywhere)
    if any(contains(model.rxns,'LMPD_'))
        ind_lumps = find(contains(model.rxns,'LMPD_'));
        model.rules(ind_lumps) = {''};
        model.grRules(ind_lumps) = {''};
        model.rxnGeneMat(ind_lumps,:) = zeros(length(ind_lumps),length(model.genes));
    end
    
    if any(contains(model.rxns,'Parasitosome'))
        ind_para = find(contains(model.rxns,'Parasitosome'));
        model.rules = [model.rules(1:ind_para(1)-1);repmat({''},length(ind_para),1);model.rules(ind_para(1):end)];
        model.grRules = [model.grRules(1:ind_para(1)-1);repmat({''},length(ind_para),1);model.grRules(ind_para(1):end)];
        matrix = zeros(length(ind_para),length(model.genes));%MODIFY THIS TOO!
        model.rxnGeneMat = [model.rxnGeneMat(1:ind_para(1)-1,:);matrix;model.rxnGeneMat(ind_para(1):end,:)];
    end
    
    if any(contains(model.rxns,'DM_'))
        ind_para = find(contains(model.rxns,'DM_'));
        model.rules = [model.rules(1:ind_para(1)-1);repmat({''},length(ind_para),1);model.rules(ind_para(1):end)];
        model.grRules = [model.grRules(1:ind_para(1)-1);repmat({''},length(ind_para),1);model.grRules(ind_para(1):end)];
        matrix = zeros(length(ind_para),length(model.genes));%MODIFY THIS TOO!
        model.rxnGeneMat = [model.rxnGeneMat(1:ind_para(1)-1,:);matrix;model.rxnGeneMat(ind_para(1):end,:)];
    end
end


  