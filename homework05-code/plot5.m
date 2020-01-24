function plot5(info_matrix, method, f)

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

x1 = info_matrix(:,2);
x2 = info_matrix(:,3);
fx = info_matrix(:,7);

padx1 = (max(x1) - min(x1))/5;
padx2 = (max(x2) - min(x2))/5;

x0 = info_matrix(1,2:3);
xval = num2str(x0(1));
yval = num2str(x0(2));



%%% Establishing Graph Components

header = ["Optimization using " + method + " method " + newline + ...
          "for f (x_1, x_2) = 100(x_2 - x_1^2)^2 + (1 - x_1)^2" + ...
           newline + " starting at (" + xval + ", " + yval + ")"];
     
xlab = "x_1";
ylab = "x_2";
zlab = "f(x_1, x_2)";

%%% Printing Results

hold on
box on
grid on
colormap jet

contours = fsurf(f, [min(x1)-padx1 max(x1)+padx1 ...
                     min(x2)-padx2 max(x2)+padx2],...
                     'ShowContours', 'on');
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
minline.Marker = 'o';
minline.MarkerFaceColor = 'k';
minline.LineWidth = 1.5;
minline.MarkerSize = 3;

start = plot3(x1(1),x2(1),100);
start.Color = 'g';
start.Marker = 'o';
start.MarkerFaceColor = 'g';
start.MarkerSize = 6;

goal = plot3(1,1,100);
goal.Color = 'r';
goal.Marker = 'o';
goal.MarkerFaceColor = 'r';
goal.MarkerSize = 6;

axis tight

end
