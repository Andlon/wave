% Solver based on characteristic curves
% Solves the system Vt + l1*Vx = 0
%                   Wt + l2*Wx = 0
close all
clear all
clc

% choice of initial conditions
% 1 - sinus 
% 2 - cylinder
% 3 - Gauss
% 4 - nice
InitialValues = 4;


T_max = 70;
a = -0.1;
X_max = 100;
N_space = 100;
N_time = 100;
delt = T_max/N_time;
t = linspace(0,T_max,N_time);
u = zeros(N_space,N_time);
eta = zeros(N_space,N_time);


xq = linspace(0,X_max,N_space);
x_V = linspace(0,X_max,N_space);
x_W = linspace(0,X_max,N_space);


c1 = @(x,t) u(x,t) + sqrt(eta(x,t));
c2 = @(x,t) u(x,t) - sqrt(eta(x,t));
% u0 = @(x) 0.5;
u0 = @(x) 1;
switch InitialValues
    case 1
        eta0 = @(x) 0.1*sin(x)+1;
        u(:,1) = u0(0);
        eta(:,1)=1;
        eta(N_space/2-10:N_space/2+10,1,1) = eta0(linspace(0,pi,21));
    case 2
        eta0 = @(x) 1.1;
        u(:,1) = u0(0);
        eta(:,1)=1;
        eta(N_space/2-10:N_space/2+10,1,1) = eta0(linspace(0,pi,21));
    case 3
        eta0 = @(x) 0.1*exp(-(x-50).^2/(2*10^2));
        eta(:,1) = eta0(xq);
    case 4
        eta0 = @(x) -0.1*cos(x)+1.1;
        u(:,1) = u0(0);
        eta(:,1)=1;
        eta(N_space/2-10:N_space/2+10,1) = eta0(linspace(0,2*pi,21));
end

V = u(:,1) + 2*sqrt(eta(:,1));
W = u(:,1) - 2*sqrt(eta(:,1));
% eta(:,1)=eta(:,1)-1/X_max*xq'

for i = 2:N_time

    lambda1 = u(1:N_space,i-1) + sqrt(eta(1:N_space,i-1));
    lambda2 = u(1:N_space,i-1) - sqrt(eta(1:N_space,i-1));
  
    x_V = x_V + lambda1'*delt;
    x_W = x_W + lambda2'*delt;
    
    % TO DO: Sort V according to x-values
    
    [x_V, x_Vindices] = sort(x_V);
    V_sorted = V(x_Vindices) + a*delt;
    
    [x_W, x_Windices] = sort(x_W);
    W_sorted = W(x_Windices)+ a*delt;
    
    
    V = interp1(x_V,V_sorted,xq,'linear','extrap');
    W = interp1(x_W,W_sorted,xq,'linear','extrap');
    
%     V(1) = V_leftboundary(t(i));
%     W(end) = W_rightboundary(t(i));
    
    x_V = xq;
    x_W = xq;
    
    u(:,i) = (V+W)/2;
    eta(:,i) = (V-W).^2/16;
  
end

figure(1)
mesh(t,xq,u)
title('Velocity')
xlabel('t')
ylabel('x')



figure(2)
% mesh(t,xq,eta - a*xq'/X_max*ones(1,N_time))
mesh(t,xq,eta)

title('Surface')
xlabel('t')
ylabel('x')