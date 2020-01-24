function plot_cg(residue_vector, dimension)

    %%% PLOT_CG Introduction
    %                   Plots the residual norm over the number of
    %                   iterations needed for convergence in the CG method.
    %
    %   residue_vector: the norm of the residue at every iteration of the
    %                   CG method.
    %
    %   dimension:      the dimension of the nxn matrix A



%%% Establishing Plot Vectors

x1 = 0:1:length(residue_vector)-1;
x2 = residue_vector;



%%% Establishing Graph Components

header = ["CG method residual in " + dimension + " dimensions"];
     
xlab = "Iterations";
ylab = "Residual 2-Norm";

%%% Printing Results

hold on
box on
grid on

title(header);
xlabel(xlab);
ylabel(ylab);

minline = plot(x1,x2);
minline.Color = 'b';
minline.LineWidth = 1.5;
axis tight


end
