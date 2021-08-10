r=28; sigma=10; beta=8/3; 
x1=0; y1=0; z1=0;
x2=sqrt(beta*(r-1)); y2=sqrt(beta*(r-1)); z2=r-1;
x3=-sqrt(beta*(r-1)); y3=-sqrt(beta*(r-1)); z3=r-1;
nx=500; nz=500;
xmin=-40; xmax=40; zmin=-40; zmax=40;
x_grid=linspace(xmin,xmax,nx); z_grid=linspace(zmin,zmax,nz);
[X,Z]=meshgrid(x_grid,z_grid);

% Write Newton's method using every grid point as the initial condition
% Perform enough iterations that every initial condition converges to a root 
% Save the x-values of the converged roots in the matrix X
% To pass the assessment, every pixel in the figure must be correctly colored

%!!!!!!!!! Set initial value y=3*sqrt(2) for all values (x,z) on the grid !!!!!!!!!!
RelTol=1.e-06; AbsTol=1.e-09;
y0=3*sqrt(2);
for i=1:nx    
  for j=1:nz    
    error=inf;    
    Y=y0;               
    while error > max(RelTol ,AbsTol)
      J=[-sigma, sigma, 0; r-Z(i, j), -1, -X(i, j); Y, X(i, j), -beta];
      rhs = -[sigma*(Y-X(i, j)); r*X(i, j)-Y-X(i, j)*Z(i, j); X(i, j)*Y-beta*Z(i, j)];
      delta_xyz=J\rhs;
      X(i, j) = X(i, j) + delta_xyz(1);
      Y = Y + delta_xyz(2);
      Z(i, j) = Z(i, j) + delta_xyz(3);
      error=max(abs(delta_xyz));
    end
  end
end

eps=1.e-03;
X1 = abs(X-x1) < eps; X2 = abs(X-x2) < eps; X3 = abs(X-x3) < eps;
X4 = ~(X1+X2+X3);
figure; 
map = [1 0 0; 0 1 0; 0 0 1; 0 0 0]; colormap(map); %[red;green;blue;black]
X=(X1+2*X2+3*X3+4*X4); 
image([xmin xmax], [zmin zmax], X); set(gca,'YDir','normal');
xlabel('$x$', 'Interpreter', 'latex', 'FontSize',14);
ylabel('$z$', 'Interpreter', 'latex', 'FontSize',14);
title('Fractal from the Lorenz Equations', 'Interpreter', 'latex','FontSize', 16)  
