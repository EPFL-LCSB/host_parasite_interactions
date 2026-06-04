function model = addActiveParasitosomeConstraint(model, ind_FU_para, ind_BU_para)

%% Description
%
% addActiveParasitosomeConstraint adds a constraint to the model that ensures
% the growth of only one parasitosome.
%
% Usage:
%   model = addActiveParasitosomeConstraint(model, ind_FU_para, ind_BU_para)
%
% Inputs:
%   model        : A TFA model structure.
%   ind_FU_para  : A vector of indices for 'FU_parasitosome' variables in the model.
%   ind_BU_para  : A vector of indices for 'BU_parasitosome' variables in the model.
%
% Output:
%   model        : The updated TFA model structure with the added constraint.
%
% Notes:
%   The function adds a constraint to the model in the form of
%   BFUNUSED + BUSE + FUSE = 1, which ensures that only one parasitosome is
%   allowed to grow. The constraint is added to the model by updating the
%   constraint names, A matrix, constraint types, and rhs values.
%
% Example:
%   ind_FU_para = find(contains(model.varNames, 'FU_parasitosome'));
%   ind_BU_para = find(contains(model.varNames, 'BU_parasitosome'));
%   model = addActiveParasitosomeConstraint(model, ind_FU_para, ind_BU_para);
%
% This will add a constraint to the model to ensure that only one parasitosome grows.

num_constr = size(model.A, 1);
constraintName = 'constr_active_parasitosomes';
constraintType = '=';
rhsValue = 1;

model.constraintNames{num_constr + 1, 1} = constraintName;
model.A(num_constr + 1, ind_FU_para) = 1;
model.A(num_constr + 1, ind_BU_para) = 1;
model.constraintType{num_constr + 1, 1} = constraintType;
model.rhs(num_constr + 1, 1) = rhsValue;

end