function info_matrix = linesearch5(method, x_0, f, param)

    %%% BFGS Introduction    
    %           Returns all relevant values for optimizatin of the given
    %           function, using either the Newton or BFGS Method.
    %
    %   method: The method used, either "Newton" or "BFGS"
    %
    %   f:      The function to be evaluated, using symbolic variables
    %
    %   x_0:    The initial point to begin the minimization search, given
    %           as a vertical vector.
    %
    %   param:  The basic parameters, given in a vector with the format:
    %           [alpha_0, tolerance]
    
    %%% Establishing Parameters
    
    x_k = x_0;
    
    I = eye(2);
    H_k = I;
    outer_iterations = 0;
    inner_iterations = 0;
    max_outer_iterations = 20;
    max_inner_iterations = 1000;
    
    %%% Setting Up Functions
    
    x                   = sym('x', [2,1]);
    
    f_evaluated(x)      = f;
    f_gradient(x)       = gradient(f_evaluated);
    f_hessian(x)        = hessian(f_evaluated);
    
    f_eval              = matlabFunction(f_evaluated, 'Vars', {x});
    f_grad              = matlabFunction(f_gradient, 'Vars', {x});
    f_hess              = matlabFunction(f_hessian, 'Vars', {x});
    
    
    %%% Line Search Algorithm
    
    i = 1;
    convergence = 0;
    info_matrix(1,:) = [outer_iterations, x_k', param(1), ...
                        0, 0, f_eval(x_k), ...
                        inner_iterations, convergence];
    
    while f_eval(x_k) > param(2)
        
        outer_iterations = outer_iterations + 1;
        if outer_iterations > max_outer_iterations
            break;
        end
        
        % Choosing search direction
        
        if method == "BFGS"
            p_k = -H_k * f_grad(x_k) / norm(f_grad(x_k));            
            convergence = norm(H_k-f_hess(x_k)*p_k) / norm(p_k);
        elseif method == "Newton"
            p_k = -inv(f_hess(x_k))' * f_grad(x_k) / norm(f_grad(x_k));
            convergence = NaN;
        end
        
        inner_iterations = 0;
        
        % Alpha refinement algorithm
        
        while 1
            
            inner_iterations = inner_iterations + 1;
            if inner_iterations > max_inner_iterations
                param(1) = 0;
                break;
            end
            
            delta = (1 + abs(param(1))) * 0.01;
            check1 = (f_eval(x_k + (param(1) + delta)*p_k) - ...
                f_eval(x_k + (param(1) - delta)*p_k)) / (2*delta);
            check2 = (f_eval(x_k + (param(1) + delta)*p_k) - ...
                2 * f_eval(x_k + param(1)*p_k) + ...
                f_eval(x_k + (param(1) - delta)*p_k) ) / (delta^2);
            
            if check2 == 0
                break;
            end
            
            condition = check1 / check2;
            param(1) = param(1) - condition;
            
            if abs(condition) < param(2)
                break
            end
            
        end
        
        % Updating Hessian approximation
        
        x_k = x_k + param(1) * p_k;
        
        s_k = param(1) * p_k;
        y_k = f_grad(x_k + param(1) * p_k) - f_grad(x_k);
        
        rho = 1 / (y_k'*s_k);
        H_k = (I-rho*s_k*y_k')*H_k*(I - rho*y_k*s_k') + rho*(s_k*s_k');  
        
        
                        
        % Storing values
        
        i = i+1;
        info_matrix(i,:) = [outer_iterations, x_k', param(1), ...
                            p_k', f_eval(x_k), ...
                            inner_iterations, convergence]; 
    
    end

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    

