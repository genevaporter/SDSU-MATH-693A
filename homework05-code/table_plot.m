function nice_table = table_plot(info_matrix)

format short

iterations          = info_matrix(:,1);
x_1                 = info_matrix(:,2);
x_2                 = info_matrix(:,3);
alpha               = info_matrix(:,4);
p_k1                = info_matrix(:,5);
p_k2                = info_matrix(:,6);
fx_k                = info_matrix(:,7);
inner_iterations    = info_matrix(:,8);
convergence         = info_matrix(:,9);

nice_table = table(iterations, x_1, x_2, alpha, p_k1, p_k2, fx_k, ...
             inner_iterations, convergence);

end

