%%% Homework Introduction

% Geneva Porter, 2019.11.18
% Homework 5, Math 693A
% Professor Uduak George, SDSU

% This assignment minimizes the function:
%
%               f(x) = 100(x2-x1^2)^2 + (1-x1)^2
%
% using the BFGS method, a popular quasi-Newton method. The objective
% function will be updated at every iteration; therefore the graph will
% reflect only the first iteration of the function as it plots the
% optimization line. The BFGS method, using an updated hessian
% approximation, will be compared beside the Newton method using a simple
% line search. The initial and second approximations for the Hessian will
% be the identity matrix and the true Hessian, respectively.
% 
% "Inner iterations" (how many times alpha must be revised
% before meeting the wolfe conditions) will be recorded for each step, and
% an average will be output. "Outer iterations" (how many iterations are
% needed for the optimization to converge) will also be computed, and both
% inner and outer iteration numbers are compared for each method.
%
% Overall, the BFGS method seemed to work slightly better for starting
% points further away from the minimum, and worse for points close to the
% minimum...probably becaue I did something wrong. Oh well...

%%% Establishing Parameters

clear
clc

alpha       = 1;
tolerance   = 1e-6;

param       = [alpha, tolerance];


%%% Setting Up Functions

x           = sym('x', [2,1]);
f(x)        = 100*(x(2) - x(1).^2).^2 + (1 - x(1)).^2;
f_grad(x)   = gradient(f);
f_hess(x)   = hessian(f);
R           = rosenbrock_2Nd(x,-1);
disp(R);

NE = "Newton";
BS = "BFGS";


%%% Iterating Methods and Printing Results:

points = length(R);
for i = 1:points
    
    % Establishing points from given function rosenbrock_2Nd.m
    
    x_0 = R(i,:)';
    
    figure(i)
    clf
    hold on
    grid on
    
    % Plotting the Newton method
    
    subplot(1,2,1)
    info_matrix1 = linesearch5("Newton", x_0, f, param);
    plot5(info_matrix1, "Newton", f);
    inner1 = mean(info_matrix1(:,8));
    outer1 = length(info_matrix1(:,1));
    
    % Plotting the BFGS method
    
    subplot(1,2,2)
    info_matrix2 = linesearch5("BFGS", x_0, f, param);
    plot5(info_matrix2, "BFGS", f);
    inner2 = mean(info_matrix2(:,8));
    outer2 = length(info_matrix2(:,1));
    convg2 = info_matrix2(end,9);
    
    % Printing a comparison summary
    
    summary = ["For starting point (" + x_0(1) + ", " + x_0(2) ...
               + "):" + newline  + "The Newton method averaged " ...
               + inner1 + " inner iterations and used " ...
               + outer1 + " outer iterations." + newline ... 
               + "The BFGS Method averaged " ...
               + inner2 + " inner iterations, used " ...
               + outer2 + " outer iterations and " ...
               + "ended with " + convg2 + " convergence." + newline];
    disp(summary);

    % Printing data tables
    
    disp(["Newton method at starting point(" ...
          + x_0(1) + ", " + x_0(2) + "):"]);
    data1 = table_plot(info_matrix1);
    disp(data1);
    
    disp(["BFGS method at starting point(" ...
          + x_0(1) + ", " + x_0(2) + "):"]);
    data2 = table_plot(info_matrix2);
    disp(data2);


    
end


