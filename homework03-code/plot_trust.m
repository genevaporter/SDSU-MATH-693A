function plot_trust(info_matrix, method, f)

    %%% PLOT_SEARCH Introduction
    %                   Plots the 3D surface plot and the minimizing path.
    %
    %   info_matrix:    the matrix produced by iterating a backtrack line
    %                   search.
    %
    %   method:         the method used for the backtrack line search
    %
    %   f               the function to be evaluated



%%% Establishing Plot Vectors

x1 = info_matrix(:,1);
x2 = info_matrix(:,2);
fx = info_matrix(:,3);

x0 = info_matrix(1,1:2);
xval = num2str(x0(1));
yval = num2str(x0(2));



%%% Establishing Graph Components

header = ["Backtracking line search using " + method + " method"...
          "for f(x_1, x_2) = 10(x_2 - x_1^2)^2 + (1 - x_1)^2" + ...
          " starting at (" + xval + ", " + yval + ")"];
     
xlab = "x_1";
ylab = "x_2";
zlab = "f(x_1, x_2)";

%%% Printing Results

hold on
box on
grid on
colormap jet

contours = fsurf(f, [-0.2 1 -1 1],'ShowContours', 'on');
contours.EdgeColor = 'none';
contours.FaceColor = 'none';
contours.MeshDensity = 50;
contours.LineWidth = 1;
title(header);
xlabel(xlab);
ylabel(ylab);
zlabel(zlab);

minline = plot3(x1,x2,fx);
minline.Color = 'k';
minline.LineWidth = 1.5;
minline.Marker = 'o';

axis tight

end
