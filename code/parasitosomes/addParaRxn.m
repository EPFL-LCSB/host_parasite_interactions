function host = addParaRxn(host,parasite,matches,rxn,rxnName)
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Add reaction to parasite GEM to extract stoichiometric coefficients easily

% get the equations for the parasite
parasite = addReaction(parasite,rxnName,rxn);
ind_bm_parasite = find(contains(parasite.mets,'biomass_c','IgnoreCase',true));

% get the normalized coefficients
coeffs_ = full(parasite.S(find(parasite.S(:,end)),end))/full(parasite.S(ind_bm_parasite,end));
mets    = parasite.mets(find(parasite.S(:,end)));
mets    = cellfun(@(x) extractBefore(x,length(x)-1), mets,'UniformOutput',false);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create a new format, proper for integration in host

rxn_host = '<=>';

for k = 1:length(coeffs_)

    coeff = coeffs_(k);
    met = matches.Source(find(contains(matches.KEGG_ID,mets(k),'IgnoreCase',true)));
    met = met{1};
    
    if coeff < 0
        text = strcat({' + '},num2str(abs(coeff)),{' '},met);
        rxn_host = insertBefore(rxn_host,'<=>',text);
    elseif coeff > 0
        text = strcat({' '},num2str(coeff),{' '},met,{' +'});
        rxn_host = insertAfter(rxn_host,'<=>',text);                
    end

end

rxn_host(1:3) = [];
rxn_host(end) = [];
rxn_host = insertBefore(rxn_host,'<=>',{' '});

if ~contains(rxn_host,'1 biomass_para')
    error("Stoichiometric coefficients in rxn host not normalized");
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Add the reaction in the host
% Add the reaction in the right format
host = addReaction(host,rxnName,rxn_host);

end