function plot_search(info_matrix, method, f, limits)

%%%     PLOT_SEARCH Introduction
%                       Plots the 3D surface plot and the minimizing path.
%
%       info_matrix:    the matrix produced by iterating a backtrack line
%                       search.
%
%       method:         the method used for the backtrack line search



%%% Establishing Plot Vectors

x1 = info_matrix(:,1);
x2 = info_matrix(:,2);
fx = info_matrix(:,3);

x0 = info_matrix(1,1:2);
xval = num2str(x0(1));
yval = num2str(x0(2));

%%% Establishing Graph Components

header = ["Backtracking line search using " + method + " method"...
          "for f(x_1, x_2) = 100(x_2 - x_1^2)^2 + (1 - x_1)^2" + ...
          " starting at (" + xval + ", " + yval + ")"];
     
xlab = "x_1";
ylab = "x_2";
zlab = "f(x_1, x_2)";

%%% Printing Results

hold on
box on
grid on
colormap cool

surface = fsurf(f, limits);
surface.EdgeColor = 'interp';
surface.FaceColor = 'none';
surface.MeshDensity = 20;
surface.LineWidth = 1.5;
title(header);
xlabel(xlab);
ylabel(ylab);
zlabel(zlab);

minline = plot3(x1,x2,fx);
minline.Color = 'k';
minline.LineWidth = 3;
minline.Marker = '.';

view(-45,60);

end

