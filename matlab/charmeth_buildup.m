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
InitialValues = 5;


T_max = 1;
a = -0.5;
X_max = 10;
N_space = 100;
N_time = 100;
delt = T_max/N_time;
t = linspace(0,T_max,N_time);
u = zeros(N_space,N_time);
eta = zeros(N_space,N_time);
wave = zeros(N_space,N_time);

xq = linspace(0,X_max,N_space);
x_V = linspace(0,X_max,N_space);


c1 = @(x,t) u(x,t) + sqrt(eta(x,t));
c2 = @(x,t) u(x,t) - sqrt(eta(x,t));
% u0 = @(x) 0.5;
u0 = @(x) x;
switch InitialValues
    case 1
        eta0 = @(x) 0.1*sin(x)+1;
        u(:,1) = u0(xq);
        eta(:,1)=1;
        eta(N_space/2-10:N_space/2+10,1,1) = eta0(linspace(0,pi,21));
    case 2
        eta0 = @(x) 1.1;
        u(N_space/2-10:N_space/2+10,1) = u0(linspace(0,pi,21));
        eta(:,1)=1;
        eta(N_space/2-10:N_space/2+10,1,1) = eta0(linspace(0,pi,21));
    case 3
        eta0 = @(x) 0.1*exp(-(x-50).^2/(2*10^2));
        eta(:,1) = eta0(xq);
    case 4
        eta0 = @(x) -0.1*cos(x)+1.1;
        u(:,1) = u0(xq);
        eta(:,1)=1;
        eta(N_space/2-10:N_space/2+10,1) = eta0(linspace(0,2*pi,21));
    case 5
        eta0 = @(x) a*(x-X_max);
        eta(:,1)= eta0(xq);
        eta(N_space/2-10:N_space/2+10,1) = eta(N_space/2-10:N_space/2+10,1)-0.1*cos((linspace(0,2*pi,21)))'+0.1;
        u(:,1) = 1+2*sqrt(eta(:,1));
end

V_mat = zeros(N_space,N_time);
W_mat = zeros(N_space,N_time);

V = u(:,1) + 2*sqrt(eta(:,1));
W = u(:,1) - 2*sqrt(eta(:,1));
W0 = @(x) u0(x)-2*sqrt(eta0(x));
% eta(:,1)=eta(:,1)-1/X_max*xq'

V_mat(:,1) = V;
W_mat(:,1) = W;

for i = 2:N_time

    %lambda1 = a*t(i)+(3*V_mat(:,1)+W_mat(:,1))/4;
    
    lambda1 = u(1:N_space,i-1) + sqrt(eta(1:N_space,i-1));
  
    x_V = x_V + lambda1'*delt;
    
    [x_V, x_Vindices] = sort(x_V);
    V_sorted = V(x_Vindices) + a*delt;
    
    W = W + a*delt;
    
    V = interp1(x_V,V_sorted,xq,'linear','extrap')';
    x_V = xq;
    x_W = xq;

    u(:,i) = (V+W)/2;
    eta(:,i) = (V-W).^2/16;
    
    V_mat(:,i) = V;
    W_mat(:,i) = W;
    
    wave(:,i) = eta(:,i)+a*xq';
  
end

figure(1)
mesh(t,xq,u)
title('Velocity')
xlabel('t')
ylabel('x')


figure(2)
% mesh(t,xq,eta - a*xq'/X_max*ones(1,N_time))
mesh(t,xq,wave)

title('Surface')
xlabel('t')
ylabel('x')

figure
mesh(t,xq,V_mat)
title('V')
xlabel('t')
ylabel('x')


figure
mesh(t,xq,W_mat)
title('W')
xlabel('t')
ylabel('x')