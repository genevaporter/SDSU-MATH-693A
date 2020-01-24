%%% Homework 4 #2 Introduction

% Geneva Porter, 2019.11.07
% Homework 4, Problem 2, Math 693A
% Professor Uduak George, SDSU

% This assignment applies the standard conjugate gradient method to solve
% the linear systems describing the Hilbert matrix in 5, 8, 12, and 20
% dimensions. The residual norms and eigenvalue spreads are plotted.
% Condition factors for CG and steepest descent methods are compared.

%%% Establishing variables and parameters

clear
clc
format short

Dimension = [5, 8, 12, 20]';
Condition_Number = zeros(4,1);

%%% Forming each n-dimensional matrix A

for k=1:4
    
    figure(k)
    clf
    hold on
    
    A = zeros(Dimension(k));
    b = ones(Dimension(k),1);
    x = zeros(Dimension(k),1);
    
    for i=1:Dimension(k)
        for j=1:Dimension(k)
            A(i,j) = 1/(i+j-1);
        end
    end
    
    % Compute condition number
    Condition_Number(k) = cond(A);
    
    % Plotting the norm of the residual
    subplot(1,2,1);
    r_norm = cg_standard(A,b,x,"no");
    plot_cg(r_norm, Dimension(k));
    
    % Plotting the spread of the eigenvalues
    subplot(1,2,2);
    plot_eg(A);
    
    
end

%%% Computing other values

% Steepest descent comparison

CG_Condition_Factor = (Condition_Number.^(1/2) - 1) ...
                       ./ (Condition_Number.^(1/2) + 1);
SD_Condition_Factor = (Condition_Number - 1) ...
                       ./ (Condition_Number + 1);

Data = table(Dimension, Condition_Number, ...
       CG_Condition_Factor, SD_Condition_Factor);
   
disp(Data);



