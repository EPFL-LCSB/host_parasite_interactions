function [model, expressedRxns, nonCoreRxnsToAdd, coreRxnsUpdated, nonCoreRxnsUpdated, coreRxnsOriginal, nonCoreRxnsOriginal, problematicRxns,nonCoreRxnsAddedRecords] = minNumRxns2Add2Core_updated(model, expressedRxns,biomassRxn, minFlux, timeLimit, directionalityFlag, save_path)
% minNumRxns2Add2Core_updated - Identifies the minimal number of non-core reactions to be added to the core reactions set to enable growth.
% 	Inputs:
% 	model (struct): A TFA model structure.
% 	expressedRxns (cell array): A list of core reactions that are expressed.
% 	minFlux (double): A minimal flux threshold for the reactions to be tested.
% 	timeLimit (double): Time limit (in seconds) for the solver.
% 	directionalityFlag (string, optional): Specifies directionality of the reactions ('forward', 'reverse', or []).
% 	save_path (string, optional): File path to save intermediate results.
% 
% 	Outputs:
% 	model (struct): The updated TFA model structure.
% 	expressedRxns (cell array): The updated list of core reactions that are expressed.
% 	nonCoreRxnsToAdd (cell array): The minimal set of non-core reactions to add to the core reactions.
% 	coreRxnsUpdated (cell array): The updated core reactions set.
% 	nonCoreRxnsUpdated (cell array): The updated non-core reactions set.
% 	coreRxnsOriginal (cell array): The original core reactions set.
% 	nonCoreRxnsOriginal (cell array): The original non-core reactions set.
% 	problematicRxns (cell array): The list of problematic reactions for which no solution was found.
% 	nonCoreRxnsAddedRecords (cell array): A record of non-core reactions added at each iteration.
% 	numRxnsToAddRecords (double array): A record of the number of reactions to add at each iteration.
% 
% 	Notes:
% 	- The basal flux of 1e-5 caused the solver to crash, which is why it was not used.
% 	- This function uses an updated algorithm to identify the minimal number of non-core reactions to be added.
%
%   Edits:
%   - Add the biomass rxn as part of the core rxn
%   - Add the newly found rxns to the core at each iteration
%   - Test for feasibility before going to the loop, faster


    if ~exist(save_path,'file')

    if nargin < 5
        directionalityFlag = [];
    end

    % Add biomass rxn to the core rxn
    expressedRxns = [biomassRxn;expressedRxns];

    coreRxnsOriginal = expressedRxns;
    coreRxnsUpdated = coreRxnsOriginal;
    nonCoreRxnsOriginal = setdiff(model.rxns, union(coreRxnsUpdated,{'biomass'}));
    nonCoreRxnsUpdated = nonCoreRxnsOriginal;
    problematicRxns = {};
    nonCoreRxnsAddedRecords = [];

    [numConstraints, numVars] = size(model.A);

    for i = 1:length(nonCoreRxnsOriginal)
        model.varNames{numVars + i, 1} = strcat('BF_UNUSE_', nonCoreRxnsOriginal{i});
        model.var_lb(numVars + i, 1) = 0;
        model.var_ub(numVars + i, 1) = 1;
        model.vartypes{numVars + i, 1} = 'B';
        model.A(:, numVars + i) = zeros(numConstraints, 1);
    end
    
    % Reset objective function
    model.f = zeros(size(model.A,2),1);

    [numConstraints, ~] = size(model.A);

    for i = 1:length(nonCoreRxnsOriginal)
        ind_bfunuse = find(ismember(model.varNames, strcat('BF_UNUSE_', nonCoreRxnsOriginal{i})));
        ind_buse = find(ismember(model.varNames, strcat('BU_', nonCoreRxnsOriginal{i})));
        ind_fuse = find(ismember(model.varNames, strcat('FU_', nonCoreRxnsOriginal{i})));

        model.constraintNames{numConstraints + i, 1} = strcat('BF_UNUSE_CONS_', nonCoreRxnsOriginal{i});
        model.A(numConstraints + i, ind_bfunuse) = 1;
        model.A(numConstraints + i, ind_buse) = 1;
        model.A(numConstraints + i, ind_fuse) = 1;
        model.constraintType{numConstraints + i, 1} = '=';
        model.rhs(numConstraints + i, 1) = 1;
    end

    [~, ind_ncRxns_Binary] = ismember(strcat('BF_UNUSE_', nonCoreRxnsOriginal), model.varNames);
    ind_ncRxns_Binary_updated = ind_ncRxns_Binary;
    ind_nonCoreRxnsAdded_record = [];
    model.objtype = -1;
     
    rxnsToTest = expressedRxns;
    
    model_original = model;
    start = 1;

    else
        load(save_path)
        start = i+1;
    end

    for i = start:length(rxnsToTest)

        model = model_original;
        model.var_ub(ind_nonCoreRxnsAdded_record) = 0;       
        
        rxn = rxnsToTest{i};
        ind_rxn_F = find(ismember(model.varNames, strcat('F_', rxn)));
        ind_rxn_R = find(ismember(model.varNames, strcat('R_', rxn)));
        
        [numConstraints, ~] = size(model.A);
        
        model.constraintNames{numConstraints + 1, 1} = strcat('flux_min_', rxn);
        if isempty(directionalityFlag)
            model.A(numConstraints + 1, ind_rxn_F) = 1;
            model.A(numConstraints + 1, ind_rxn_R) = 1;
        elseif strcmp(directionalityFlag, 'forward')
            model.A(numConstraints + 1, ind_rxn_F) = 1;
        elseif strcmp(directionalityFlag, 'reverse')
            model.A(numConstraints + 1, ind_rxn_R) = 1;
        end
        
        model.constraintType{numConstraints + 1, 1} = '>';
        model.rhs(numConstraints + 1, 1) = minFlux;  

        % Test for feasibility
        model.f = zeros(size(model.f));
        model.var_lb(ind_ncRxns_Binary_updated) = 1;
        sol_feas = solveTFAmodel_selections(model,'emphPar',0);


        % If not feasible, find the reactions one needs to add
        if isempty(sol_feas.val)

            % Readjust the bounds
            model.var_lb(ind_ncRxns_Binary_updated) = 0;           

            % Set the objective function
            model.f = zeros(size(model.f));            
            model.f(ind_ncRxns_Binary_updated) = 1;

            % Solve
            sol = solveTFAmodel_selections(model,'TimeInSec',timeLimit,'timePolishing',timeLimit*3/4);                 

            % If solution found, update the results
            if ~isempty(sol.val)                
                ind_nonCoreRxnsAdded = ind_ncRxns_Binary_updated(sol.x(ind_ncRxns_Binary_updated) <= 1e-9);
                ind_nonCoreRxnsAdded_record = union(ind_nonCoreRxnsAdded_record,ind_nonCoreRxnsAdded);
                ind_ncRxns_Binary_updated = setdiff(ind_ncRxns_Binary_updated,ind_nonCoreRxnsAdded);
                nonCoreRxnsAdded = strrep(model.varNames(ind_nonCoreRxnsAdded), 'BF_UNUSE_', '');            
                coreRxnsUpdated = union(coreRxnsUpdated, nonCoreRxnsAdded);
                nonCoreRxnsUpdated = setdiff(nonCoreRxnsUpdated, nonCoreRxnsAdded);                   
            else
                problematicRxns = [problematicRxns; expressedRxns(i)];
                nonCoreRxnsAdded = [];
                numRxnsToAdd = nan;
            end
                   
            nonCoreRxnsAddedRecords = [nonCoreRxnsAddedRecords; {nonCoreRxnsAdded}];
    
            if rem(i,round(length(rxnsToTest)/100)) == 0, save(save_path); end

        end
            
    end
   

    nonCoreRxnsToAdd = setdiff(coreRxnsUpdated, expressedRxns);

end

