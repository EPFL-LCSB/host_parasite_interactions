function [Prxns,T_Prxns] = biggMetsANALYSIS(model,Pmets)
% Function to extract the reactions where the Pmets_e participate

% Now we find the reactions in which these [e] mets participate. 

[~,indPM] = ismember(Pmets, model.mets);

% What is the point of indPM(find(indPM)) ? Only select in indPM the values
% that are nonzeros ? 
indPM = indPM(find(indPM));
Prxns = [];
for i = 1:length(indPM)
    indPR = find(model.S(indPM(i),:));
    aux_r = model.rxns(indPR);
    for j = 1:length(indPR)
        if ~isempty(setdiff(aux_r(j),Prxns))
            Prxns = [Prxns; aux_r(j)];
        end
    end
end

[~,indPR] = ismember(Prxns, model.rxns);
T_Prxns = [Prxns printRxnFormula(model, Prxns, 0) num2cell(model.lb(indPR)) num2cell(model.ub(indPR))];

end