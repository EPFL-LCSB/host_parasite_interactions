%% Prepare & Clean the environment

% Clean the environment
clear; clc; restoredefaultpath

% Define paths and parameters
clear paths;
paths = {
    './';    
    '/PATH/TO/CPLEX';
    '/PATH/TO/MATTFA';
    '/PATH/TO/REGEM';
    };

% Add the folders to the paths
for i = 1:length(paths)    
    addpath(genpath(paths{i}));
end

% Check for available solvers
cplex_available = ~isempty(which('cplex.p'));

if cplex_available
    solver = 'cplex_direct';
else
    error('No available solver found. Please install CPLEX and add it to the MATLAB path.');
end

% Set the chosen solver
global TFA_MILP_SOLVER
TFA_MILP_SOLVER = solver;

global TFA_LP_SOLVER
TFA_LP_SOLVER = solver;

CplexParameters = 'LCSBDefault';
GurobiParameters = 'LCSBDefault';

% Define the type of medium to study
medium_type = 'prototrophic'; % select{'auxotrophic','prototrophic'}

%% Generating host, parasite and host-parasite metabolic models and performing gene essentiality analysis.

% Step 1: generation of a hepatocyte model
step1_context_specific_host_GEM_reconstruction()

% Step 2: generation of parasitosomes
step2_parasitosomes_generation(medium_type)

% Step 3: reconstruction of an integrated host-parasite model
step3_1_matching_metabolites()
step3_2_integrating_metabolites(medium_type)
step3_3_integrating_parasitosomes(medium_type)

% Step 4: perform gene essentiality analysis
step4_gene_essentiality_analyses(medium_type)
