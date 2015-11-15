% Solver based on characteristic curves
% Solves the system Vt + l1*Vx = 0
%                   Wt + l2*Wx = 0
close all
clear all

T_max = 0.1;
X_max = 10;
N_space = 100;
N_time = 100;
delt = T_max/N_time;
t = linspace(0,T_max,N_time);
u = zeros(N_space,N_time);
eta = zeros(N_space,N_time);

c1 = @(x,t) u(x,t) + sqrt(eta(x,t));
c2 = @(x,t) u(x,t) - sqrt(eta(x,t));
u0 = @(x) 1;
eta0 = @(x) sin(x)+1;

xq = linspace(0,X_max,N_space);

x_V = linspace(0,X_max,N_space);
x_W = linspace(0,X_max,N_space);

u(:,1) = u0(linspace(0,X_max,N_space));
eta(:,1) = eta0(linspace(0,X_max,N_space));

V_leftboundary = @(t) 1;
W_rightboundary = @(t) 1;

V = zeros(N_space,1);
W = zeros(N_space,1);

for k = 1:N_space
    V(k) = u0(x_V(k)) + 2*sqrt(eta0(x_V(k)));
    W(k) = u0(x_W(k)) - 2*sqrt(eta0(x_W(k)));
end

V(1) = V_leftboundary(0);
W(end) = W_rightboundary(0);


for i = 2:N_time

    lambda1 = u(1:N_space,i-1) + sqrt(eta(1:N_space,i-1));
    lambda2 = u(1:N_space,i-1) - sqrt(eta(1:N_space,i-1));
    
    x_V = x_V + lambda1'*delt;
    x_W = x_W + lambda2'*delt;
    
    % TO DO: Sort V according to x-values
    
    [x_V, x_Vindices] = sort(x_V);
    V_sorted = V(x_Vindices);
    
    [x_W, x_Windices] = sort(x_W);
    W_sorted = W(x_Windices);
    
    
    V = interp1(x_V,V_sorted,xq);
    W = interp1(x_W,W_sorted,xq);
    
    V(1) = V_leftboundary(t(i));
    W(end) = W_rightboundary(t(i));
    
    x_V = xq;
    x_W = xq;
    
    u(:,i) = (V+W)/2;
    eta(:,i) = (V-W).^2/16;
  
end

figure
mesh(t,xq,u)
title('Velocity')
xlabel('t')
ylabel('x')


figure
mesh(t,xq,eta)
title('Surface')
xlabel('t')
ylabel('x')