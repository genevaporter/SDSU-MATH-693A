function plot_eg(A)

    %%% PLOT_EG Introduction
    %                   Plots the eigenvalue spread of a matrix A

    
%%% Establishing Plot Vectors

e = eig(A);

%%% Establishing Graph Components

xlabel('Real');
ylabel('Imaginary');
title("Eigenvalues of the Hilbert Matrix");

%%% Printing Results

hold on
box on
grid on

Eigenvalues = plot(real(e), imag(e), 'bo');
Eigenvalues.MarkerSize = 5;
Eigenvalues.MarkerFaceColor = 'b';

axis tight

end

