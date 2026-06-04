function model = setParasitosomeBounds(model, ind_parasitosomes_F, lb, ub)
    model.var_lb(ind_parasitosomes_F) = lb;
    model.var_ub(ind_parasitosomes_F) = ub;
end