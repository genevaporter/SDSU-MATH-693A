%%% Homework Introduction

% Geneva Porter, 2019.09.26
% Homework 1, Math 693A
% Professor Uduak George, SDSU

% This assignment programs the steepest descent and Newton algorithms using
% the backtracking line search to minimize the Rosenbrock function:
%
%               f(x) = 100(x2-x1^2)^2 + (1-x1)^2
%
% The initial step length is alpha_0 = 1, and each step length used by each
% method is reported at each iteration. First, the initial point (1.2,1.2)
% is used, then the more difficult point (-1.2, 1) is used for each method.
% The suggested values of alpha, rho, and c are used, as shown below. Since
% we know that the minimum of the function is (0,0) from straightforward
% analysis, the iteration stops when the absolute value of f or the norm of 
% the gradient is less than our tolerance, 10^(-8). Only the initial
% values, the first 10 iterations, and the last iteration are shown in the
% output. The function backtrack_line_search contains the algorithm that
% produces these results.


%%% Establishing Parameters

clear
clc

alpha       = 1.0;
rho         = 0.5;
c           = 1e-4;

tolerance   = 10e-8;

param = [alpha, rho, c, tolerance];

p1          = [ 1.2; 1.2];
p2          = [-1.2; 1.0];



%%% Setting Up Functions

x    = sym('x', [2,1]);
f(x) = 100*(x(2) - x(1)^2)^2 + (1 - x(1))^2;

NE = "Newton";
SD = "steepest descent";



%%% Iterating Methods and Printing Results:

SD_point1 = backtrack_line_search(SD, f, p1, param);
figure(1)
plot_search(SD_point1, SD, f, [1 1.25 1 1.3]);

SD_point2 = backtrack_line_search(SD, f, p2, param);
figure(2)
plot_search(SD_point2, SD, f, [-1.5 1 0 1.5]);

newton_point1 = backtrack_line_search(NE, f, p1, param);
figure(3)
plot_search(newton_point1, NE, f, [1 1.3 1 1.5]);

newton_point2 = backtrack_line_search(NE, f, p2, param);
figure(4)
plot_search(newton_point2, NE, f, [-2 2 0 1.5]);


