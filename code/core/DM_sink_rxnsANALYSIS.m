function DM_sink_rxns = DM_sink_rxnsANALYSIS(model)

exchgRxns = extract_drains(model);
DM_sink_rxns=[];
for i=1:length(exchgRxns)
    if ~isempty(strfind(exchgRxns{i},'DM')) || ~isempty(strfind(exchgRxns{i},'sink'))
        DM_sink_rxns=[DM_sink_rxns,exchgRxns(i)];
    end
    
end
end