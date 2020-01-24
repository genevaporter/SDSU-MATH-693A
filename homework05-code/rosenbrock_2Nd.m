function [R] = rosenbrock_2Nd(x,order)

if( nargin<2 )
  order=0;
end

%
% Initial condition in R18
%
if( order == - 1 )
  xN     = [  1 , 1 ];
  x0easy = [  1.2 , 1.2 ];
  x0e2   = (xN + x0easy) / 2;
  x0e3   = (xN + x0e2)  / 2;
  x0e4   = (xN + x0e3) / 2;
  x0hard = [ -1.2 , 1.0 ];
  x0h2   = (xN + x0hard) / 2;
  x0h3   = (xN + x0h2)  / 2;
  x0h4   = (xN + x0h3) / 2;
  x0h5   = 2*x0hard;
  R      = [x0easy ; x0e2 ; x0e3 ; x0e4 ; ...
	    x0hard ; x0h2 ; x0h3 ; x0h4; x0h5];
  return
end

nx = length(x);

% 1D versions
%
% The function and derivatives needed to compute the
% gradient and the Hessian
%
rb2d      = @(x) ( 100*(x(2)-x(1).^2).^2+(1-x(1)).^2 );
rb2d_x    = @(x) ( -400*(x(2)-x(1).^2).*x(1)-2+2*x(1) );
rb2d_xx   = @(x) ( 1200*x(1).^2-400*x(2)+2 );
rb2d_xy   = @(x) ( -400*x(1) );
rb2d_y    = @(x) ( 200*x(2)-200*x(1).^2 );
rb2d_yy   = @(x) ( 200 );
rb2d_grad = @(x) ( [ rb2d_x(x) ; rb2d_y(x) ] );
rb2d_hess = @(x) ( [ rb2d_xx(x) rb2d_xy(x) ; rb2d_xy(x) rb2d_yy(x)] );

switch order
 case 0
  R = zeros(1,1);
  for k = 1:2:nx
    R = R + rb2d(x(k:(k+1)));
  end
 case 1
  R = zeros(length(x),1);
  for k = 1:2:nx
    R(k:(k+1)) = rb2d_grad(x(k:(k+1)));
  end
 case 2
  R = zeros(length(x),length(x));
  for k = 1:2:nx
    R(k:(k+1),k:(k+1)) = rb2d_hess(x(k:(k+1)));
  end
 otherwise
  error(sprintf('\nCannot compute derivatives of order %d.\n',order));
end
