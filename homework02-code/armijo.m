function boolean = armijo(q_eval, q_grad, p_k, a_k, c1)

    %%% ARMIJO Introduction
    %       Returns true if Armijo conditio is met, returns false if the 
    %       Armijo condition is not met.
    %
    %   q_eval: function handle with respect to alpha for
    %               q_eval(alpha) = f_eval(x_k + alpha * p_k)
    %
    %   q_grad: function handle with respect to alpha for
    %               q_grad(alpha) = f_grad(x_k + alpha * p_k)
    %
    %   a_k:    current estimate for alpha
    %
    %   c1:     given parameter

    if q_eval(a_k) <= (q_eval(0) + c1*a_k*p_k'*q_grad(0))
        boolean = 1;
    else
        boolean = 0;
    end

end

