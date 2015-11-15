% Solver based on characteristic curves
% Solves the system Vt + l1*Vx = 0
%                   Wt + l2*Wx = 0

T_max = 100;
X_max = 100;
N_space = 100;
N_time = 100;
delt = T_max/N_time;
t = linspace(0,T_max,N_time);
u = zeros(N_space,N_time);
eta = zeros(N_space,N_time);

c1 = @(x,t) u(x,t) + sqrt(eta(x,t));
c2 = @(x,t) u(x,t) - sqrt(eta(x,t));
%V =  @(x,t) u(x,t) + 2*sqrt(eta(x,t));
%W = @(x,t) u(x,t) - 2*sqrt(eta(x,t));
u0 = @(x) 1;
eta0 = @(x) sin(x)+1;

xq = linspace(0,X_max,N_space);

x_V = linspace(0,X_max,N_space);
x_W = linspace(0,X_max,N_space);

u(:,1) = u0(linspace(0,X_max,N_space));
eta(:,1) = eta0(linspace(0,X_max,N_space));

V_leftboundary = @(t) 1;
W_rightboundary = @(t) 1;



V = u0(x_V) + 2*sqrt(eta0(x_V));
W = u0(x_W) - 2*sqrt(eta0(x_W));

V(1) = V_leftboundary(0);
W(end) = W_rightboundary(0);

for i = 2:N_time
    lambda1 = c1(x_V,t(i-1));
    lambda2 = c2(x_W,t(i-1));
    
    x_V = x_V + lambda1*delt;
    x_W = x_W + lambda2*delt;
    
    % TO DO: Sort V according to x-values
    
    [x_V, x_Vindices] = sort(x_V);
    Vted_sor = V(x_Vindices);
    
    [x_W(i,:), x_Windices] = sort(x_W(i,:));
    W_sorted = V(x_Windices);
    
    
    
    
    
    
    % Interpolate V and W values to find u and eta
    % 
    
end

