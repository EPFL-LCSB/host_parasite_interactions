function [growth_boolean,growth_ratio_host,growth_ratio_para] = geneEssentialityThermoInfectedHost(model, geneList, para_active, grKO_percent_host, grKO_percent_para, timeSOLVE,ratioFlag)
% geneEssentialityThermoInfectedHost evaluates the feasibility of single-gene
% knockout models in a host-parasite metabolic system.
%
% USAGE:
%   [growth_boolean, growth_ratio_host, growth_ratio_para] = geneEssentialityThermoInfectedHost(model, geneList, para_active, grKO_percent_host, grKO_percent_para, timeSOLVE, ratioFlag)
%
% INPUTS:
%   model: A TFA model structure.
%   geneList (optional): A cell array of gene IDs to be tested for knockout. Default is all genes in the model.
%   para_active (optional): A binary flag indicating if the parasite is active (1) or not (0). Default is 0.
%   grKO_percent_host (optional): A scalar representing the percentage of host growth rate to be used as a lower bound. Default is 0.1.
%   grKO_percent_para (optional): A scalar representing the percentage of parasite growth rate to be used as a lower bound. Default is 0.1.
%   timeSOLVE (optional): A scalar representing the time limit in seconds for the solver. Default is 2*60.
%   ratioFlag (optional): A binary flag indicating whether or not to compute the growth ratios (host growth rate / maximal host growth rate, parasite growth rate / maximal parasite growth rate). Default is 0.
%
% OUTPUTS:
%   growth_boolean: A binary vector of length equal to the number of genes in geneList, with 1 indicating that the model is feasible upon gene knockout and 0 otherwise.
%   growth_ratio_host: A vector of host growth ratios (host growth rate / maximal host growth rate) for each feasible knockout model. If ratioFlag is set to 0, this will be an empty vector.
%   growth_ratio_para: A vector of parasite growth ratios (parasite growth rate / maximal parasite growth rate) for each feasible knockout model. If para_active is set to 0 or if ratioFlag is set to 0, this will be an empty vector.
%
% NOTES:
%   The function first calculates the maximal biomass of the host and parasite (if active).
%   Then, it iterates through the genes in geneList, performing single-gene knockouts and
%   assessing the feasibility of the resulting models. If the gene knockout has an effect,
%   the feasibility of the model is tested, and the growth vector is updated accordingly.
%



%% Default argument values
if nargin < 2 || isempty(geneList), geneList = model.genes; end
if nargin < 3 || isempty(para_active), para_active = 0; end
if nargin < 4 || isempty(grKO_percent_host), grKO_percent_host = 0.1; end
if nargin < 5 || isempty(grKO_percent_para), grKO_percent_para = 0.1; end
if nargin < 6 || isempty(timeSOLVE), timeSOLVE = 2*60; end
if nargin < 7 || isempty(ratioFlag), ratioFlag = 0; end

%% Prepare the environment

% Fix values for key variables
nDelGenes = length(geneList);
growth_boolean = zeros(nDelGenes, 1);
growth_ratio_host = zeros(nDelGenes, 1);
growth_ratio_para = zeros(nDelGenes, 1);
ind_bm_host = find(ismember(model.varNames, 'F_biomass'));
ind_bm_para = find(ismember(model.varNames, 'F_DM_biomass_para'));

% Compute host maximal biomass
model = setModelObjective(model, ind_bm_host);
solWT_max = solveTFAmodelCplex_selections(model, 'TimeInSec', timeSOLVE);
grRateWTmax = solWT_max.val;

% Compute parasite maximal biomass and fix it to a given value
if para_active
    model = setModelObjective(model, ind_bm_para);
    solPARA_max = solveTFAmodelCplex_selections(model, 'TimeInSec', timeSOLVE);
    grRatePARAmax = solPARA_max.val;
    model.var_lb(ind_bm_para) = grKO_percent_para * grRatePARAmax;
end

% Fix host maximal value to the value of choice
model.var_lb(ind_bm_host) = grKO_percent_host * grRateWTmax;

% Set objective function to zero, just testing for feasibility
model = setModelObjective(model, []);

%% Conduct the gene essentiality

% Deleting genes one at a time
for i = 1:nDelGenes
    gene = geneList{i};
    [modelDel, hasEffect, ~, ~] = thermoDeleteModelGenes(model, gene);

    if hasEffect
        solFEAS = solveTFAmodelCplex_selections(modelDel, 'TimeInSec', timeSOLVE, 'emphPar' ,0);
        if ~isempty(solFEAS.x)
            growth_boolean(i) = 1;
            if ratioFlag == 1
                % test host growth
                modelDel = setModelObjective(modelDel, ind_bm_host);
                grHOST = solveTFAmodelCplex_selections(model, 'TimeInSec', timeSOLVE).val;
                growth_ratio_host(i) = grHOST/grRateWTmax;
                % test parasite growth
                modelDel = setModelObjective(modelDel, ind_bm_para);
                grPARA = solveTFAmodelCplex_selections(model, 'TimeInSec', timeSOLVE).val;
                growth_ratio_para(i) = grPARA/grRatePARAmax;
            end
        end
    else
        growth_boolean(i) = 1;
    end
end

%% Helper function

function model = setModelObjective(model, objIdx)
    model.f = zeros(size(model.f));
    model.f(objIdx) = 1;
    model.objtype = -1;
end

end