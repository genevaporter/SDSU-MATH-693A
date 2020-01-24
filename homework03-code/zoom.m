function a = zoom(q_eval, q_grad, p_k, alpha_0, alpha_1)

    %%% ZOOM Introduction
    %       Interpolates previous two estimates for alpha, reducing the
    %       search interval for alpha.
    %
    %   q_eval: function handle with respect to alpha for
    %               q_eval(alpha) = f_eval(x_k + alpha * p_k)
    %
    %   q_grad: function handle with respect to alpha for
    %               q_grad(alpha) = f_grad(x_k + alpha * p_k)
    %
    %   a_0:    previous estimate for alpha
    %
    %   a_1:    current estimate for alpha
    
    a_0 = alpha_0;
    a_1 = alpha_1;
    
    d1 = p_k' * (q_grad(a_0) + q_grad(a_1)) - 3 * ...
         (q_eval(a_0) - q_eval(a_1)) / (a_0 - a_1);
     
    d2 = (a_1 - a_0) / abs(a_1 - a_0) * ...
         (d1^2 - q_grad(a_0)' * q_grad(a_1))^(1/2);
    
    a = a_1 - (a_1 - a_0) * (p_k'*q_grad(a_1) + d1 + d2) / ...
        (p_k' * (q_grad(a_1) - q_grad(a_0)) + 2 * d2);
    
end
