%%% Homework 4 #1 Introduction

% Geneva Porter, 2019.11.07
% Homework 4, Problem 1, Math 693A
% Professor Uduak George, SDSU

% This assignment applies the standard conjugate gradient method to solve
% the linear systems describing the Helical Coordinate Preconditioner for
% the Laplacian in one, two, and three dimensions. The variable "n" on line
% 20 can be changed to represent the size of the n x n matrix A. The
% function I created for using the CG method is cg_standard.m.



clear
clc

%%% Establishing Parameters

n = 3;

%%% One-dimensional solution

d1 = ones(n,1);
x1 = zeros(size(d1));
b1 = ones(size(d1));
A1 = spdiags([d1, -2*d1, d1], [-1 0 1], n, n);
r_norm_1 = cg_standard(A1,b1,x1,"yes");

%%% Two-dimensional solution

d2 = ones(n^2, 1);
x2 = zeros(size(d2));
b2 = ones(size(d2));
A2 = spdiags([d2 d2 -4*d2 d2 d2], [-n -1 0 1 n], n^2, n^2);
r_norm_2 = cg_standard(A2,b2,x2,"yes");

%%% Three-dimensional solution

d3 = ones(n^3, 1);
x3 = zeros(size(d3));
b3 = ones(size(d3));
A3 = spdiags([d3 d3 d3 -6*d3 d3 d3 d3], [-n^2 -n -1 0 1 n n^2], n^3, n^3);
r_norm_3 = cg_standard(A3,b3,x3,"yes");

