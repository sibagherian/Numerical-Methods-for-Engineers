%%%%% Define the square and grid parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
L=1;  %square is 2L x 2L 
N=100; %# of intervals in x and y directions
n=N+1; %# of gridpoints in x,y directions including boundaries
h=2*L/N;  %grid size in x,y directions
x=-L + (0:N)*h; %x values on the grid
y=-L + (0:N)*h; %y values on the grid
[X,Y]=meshgrid(x,y);
%%%%% Define the indices associated with the boundaries %%%%%%%%%%%%%%%%%%%
% boundary_index = [bottom, left, top, right]
boundary_index=[          1:n,       1:n:1+(n-1)*n, ...
                1+(n-1)*n:n*n,   n:n:n*n           ]; 
%%%%% Diffusion constant and time-step parameters
D=1;
dt=h^2/(2*D); %borderline stability of FTCS scheme
alpha=dt*D/h^2; %equation parameter
nsteps=1000; %number of time steps
%%%%% CONSTRUCT THE MATRIX AND COMPUTE LU DECOMPOSITION %%%%%%%%%%%%%%%%%%%%
diagonals = [2*(1+2*alpha)*ones(n*n,1), -alpha*ones(n*n,4)];
A=spdiags(diagonals,[0 -1 1 -n n], n*n, n*n);
I=speye(n*n);
A(boundary_index,:)=I(boundary_index,:);
[L1,U1]=lu(A);
%%%%% Define initial conditions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
u=zeros(n,n,nsteps);
sigma=L/4;
u(:,:,1)=1/(2*pi*sigma^2)*exp(-0.5*(X.^2+Y.^2)/sigma^2); 
u(1,:,1)=0; u(n,:,1)=0; u(:,1,1)=0; u(:,n,1)=0; %b.c.
%%%%% ADVANCE SOLUTION u %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
b=zeros(n*n,n*n);
bdiagonals = [2*(1-2*alpha)*ones(n*n,1), alpha*ones(n*n,4)];
b=spdiags(bdiagonals,[0 -1 1 -n n], n*n, n*n);
b(boundary_index,:)=I(boundary_index,:);
for m=2:nsteps
  s = b*reshape(u(:,:,m-1), n*n, 1);
  u(:,:,m)=reshape(U1\(L1\s), n,n);
end
% %%%%% Plot with animation:  UNCOMMENT TO RUN ON MATLAB ONLINE OR DESKTOP %%%
% figure('units','normalized','outerposition',[0 0 1 1])
% s=surf(X,Y,u(:,:,1)); zlim([0, 2.6]);
% xlabel('$x$','Interpreter','latex','FontSize',14); 
% ylabel('$y$','Interpreter','latex','FontSize',14); 
% zlabel('$u(x,y,t)$','Interpreter','latex','FontSize',14); 
% title('Solution of the 2D diffusion equation','Interpreter','latex','FontSize',16);
% pause(1)
% for j=2:nsteps
%     s.ZData=u(:,:,j); pause(0.01);
% end
