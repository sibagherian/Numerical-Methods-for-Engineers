mu_min=2.4; mu_max=4; %range of mu values
n_mu=500; %number of mu pixels
n_x=400; %number of x pixels
mu_edges=linspace(mu_min,mu_max,n_mu+1); %edges of mu pixels
mu=(mu_edges(1:n_mu)+mu_edges(2:n_mu+1))/2; %values of mu on which to perform computation
x_edges=linspace(0,1,n_x+1); %edges of x pixels

n_trans=20000; %transient iterations
n_data=10000;  %number of x values per mu value

x_data=zeros(n_data,n_mu); %x-data used to construct figure

x_0=0.5; %initial condition

% WRITE THE COMPUTATIONAL ENGINE OF THE CODE.
% USE THE ALREADY DEFINED PARAMETERS AND VARIABLES: n_mu, mu, x_0, n_trans, n_data.
% YOUR FINAL RESULT WILL BE THE VARIABLE x_data, and this variable will be assessed.

x_trans = x_0*ones(1,500);
for a=1:n_trans
    x_trans = mu .* x_trans .* ( 1 - x_trans );
end
for b=1:n_mu
    x_data(1,b)= mu(b) * x_trans(b) * ( 1 - x_trans(b) );
    for c=2:n_data
        x_data(c,b) = mu(b) * x_data(c-1,b) * ( 1 - x_data(c-1,b) );    
    end
end

%%%%% bin data and plot image %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x_histogram=zeros(n_x,n_mu); %binned values of x
for i=1:n_mu
x_histogram(:,i)=histcounts(x_data(:,i),x_edges);
x_histogram(:,i)=255*x_histogram(:,i)/max(x_histogram(:,i));
end
colormap(flipud(gray(256))); brighten(-0.8); cmap=colormap;
im=image([mu_edges(1) mu_edges(end)], [x_edges(1) x_edges(end)], x_histogram);
set(gca,'YDir','normal');
xlabel('$\mu$','Interpreter','latex','FontSize',14);
ylabel('$x\;\;$','Interpreter','latex','FontSize',14);
title('Logistic Map Bifurcation Diagram','Interpreter','latex','FontSize',16)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%