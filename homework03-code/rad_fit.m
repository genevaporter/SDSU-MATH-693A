function [p_k, rho] = rad_fit(f, m_k, x_k, delta)

    %%% RAD_FIT Introduction
    %                   Establishes the scalar for the trust region radius
    %
    %   f:              The function beig evaluated
    %
    %   m_k:            The quadratic function used for optimization
    %
    %   delta:          The current trust radius

    
    
    %%% SETTING UP FUNCTIONS
    
    x          = sym('x', [2,1]);
    f_grad(x)  = gradient(f);
    f_hess(x)  = hessian(f);

    f_eval     = matlabFunction(f, 'Vars', {x});
    f_grad     = matlabFunction(f_grad, 'Vars', {x});
    f_hess     = matlabFunction(f_hess, 'Vars', {x});

    
    
    %%% ESTABLISHING TRUST REGION SCALAR
    
    tau = min(1, norm(f_grad(x_k))^3 / delta * ...
          f_grad(x_k)' * f_hess(x_k) * f_grad(x_k));
      
    p_k = -tau * delta * f_grad(x_k) / norm(f_grad(x_k));
    
    rho = (f_eval(x_k) - f_eval(x_k + p_k))/...
          (m_k([0;0]) - m_k(p_k));

end

