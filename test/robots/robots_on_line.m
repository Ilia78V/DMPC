%% Two-Cart “Maintain 0.5 m Gap” DMPC Test
clear; clc;

%% (Optional) Add YALMIP/IPOPT paths
if false
    originalPath = 'M:\DMPC proj\matlab';
    addpath(genpath('M:\DMPC proj\matlab\YALMIP-master'));
    optiPath = 'M:\DMPC proj\matlab\OPTI-master';
    cd(optiPath);
    run('opti_Install.m');
    cd(originalPath);
end

%% Physical params
m1 = 1;   % kg
m2 = 1;   % kg
b1 = 0.1; % friction coeff Cart1
b2 = 0.1; % friction coeff Cart2


%% Horizon & sim
t0        = 0;
T_h       = 2;    % 2 s to go from 2→3
N         = 11;   %  dt = 0.2 s
dt        = T_h/(N-1);
T_sim     = 2;    % run longer to see gap hold
dt_sample = dt;

%% State/control dims
n_x = 2;  % [position; velocity]
n_u = 1;  % scalar force

%% Initial & local refs
% Cart1: start at 2, want position=3
x0_1   = [2;0];
x_ref1 = [3;0];  
% Cart2: start at 0, no personal ref
x0_2   = [0;0];
x_ref2 = [0; 0];  

%% Bounds
x_min = [-5; 0];  x_max = [ 5;  5];
u_min = -5;        u_max =  5;

rho_init = 100;

%% approximation flags
approx = containers.Map('KeyType','char','ValueType','logical');
appr = true;
approx('cost')        = appr;   % werll use cost coupling
approx('dynamics')    = appr;
approx('constraints') = appr;

%% Create Agent_data
agentData1 = Agent_data(1,n_x,n_u,t0,T_h,N,x0_1,x_ref1,x_min,x_max,u_min,u_max,rho_init,approx);
agentData2 = Agent_data(2,n_x,n_u,t0,T_h,N,x0_2,x_ref2,x_min,x_max,u_min,u_max,rho_init,approx);

%% Cost weights
Q1  = diag([250, 50]);    % strongly penalize position error for Cart1
R1  = 10;
Qf1 = diag([250, 50]);

Q2  = diag([0, 50]);      % no personal tracking cost for Cart2
R2  = 0.01;
Qf2 = diag([0, 50]);

%% Local stage & terminal costs
l1 = @(x,u,t) (x - x_ref1)'*Q1*(x - x_ref1) + R1*u^2;
V1 = @(x,T)       (x - x_ref1)'*Qf1*(x - x_ref1);
l_12 = @(x1,u1,x2,u2,t) 0;  
V_12 = @(x1,x2,T) 0;

l2 = @(x,u,t) (x - x_ref2)'*Q2*(x - x_ref2) + R2*u^2;   % Cart2 only pays control effort
V2 = @(x,T)   (x - x_ref2)'*Qf2*(x - x_ref2);
l_21 = @(x2,u2,x1,u1,t) 1e4* ( (x2(1)-(x1(1)-0.5))^2 ) + 10*( x2(2) - x1(2) )^2; % + 200*( (t - T_sim)^2 ); *(1 + 5*t)  
V_21 = @(x2,x1,T)       1e4* ( (x2(1)-(x1(1)-0.5))^2 ); %+ 200*( (T - T_sim)^2 );

%% Dynamics
f1 = @(x,u) [ x(2); (u - b1*x(2))/m1 ];
f2 = @(x,u) [ x(2); (u - b2*x(2))/m2 ];

%% Zero constraints
g_i   = @(x,u,t) 0;   h_i   = @(x,u,t) 0;
g_ij  = @(x,u,xn,un,t) 0;   h_ij  = @(x,u,xn,un,t) 0;
g_i_N = @(x,t) 0;     h_i_N = @(x,t) 0;
g_ij_N= @(x,xn,t) 0;  h_ij_N= @(x,xn,t) 0;

%% Build Agents
agent1 = Agent(1, f1, @(x,T)V1(x,T), @(x,u,t) l1(x,u,t), ...
               g_i, g_i_N, h_i, h_i_N, agentData1);
agent2 = Agent(2, f2, @(x,T)V2(x,T), @(x,u,t) l2(x,u,t), ...
               g_i, g_i_N, h_i, h_i_N, agentData2);

%% Neighbor for cost only (no dynamic coupling)
nbrData1 = Neighbor_data(1, n_x, n_u, agent2, rho_init);
nbr1     = Neighbor(1, agent2, true, true, ...
                   @(x,u,xn,un) [0;0], ...   % no dynamic coupling
                   g_ij, g_ij_N, h_ij, h_ij_N, V_21, l_21, nbrData1);

nbrData2 = Neighbor_data(2, n_x, n_u, agent1, rho_init);
nbr2     = Neighbor(2, agent1, true, true, ...
                   @(x,u,xn,un) [0;0], ...   % no dynamic coupling
                   g_ij, g_ij_N, h_ij, h_ij_N, V_12, l_12, nbrData2);

agent1.register_neighbors({nbr2});       % Cart1 no coupling
agent2.register_neighbors({nbr1});   % Cart2 cost-coupled to Cart1

%% Solution containers & solver
sol1     = Solution(agent1, dt_sample);
sol2     = Solution(agent2, dt_sample);
solutions = {sol1, sol2};

rho_init  = 1.0;
max_iters = 30;
tol       = 1e-3;

solver = ADMM_Solver('ipopt', {agent1,agent2}, solutions, max_iters, tol, approx);

%% Run distributed MPC
for t_sim = 0:dt_sample:T_sim
    fprintf('t = %.2f  ',t_sim);
    solver.solve();
    solver.shift(dt_sample);
end

%% Plot the Results
figure;
subplot(2,2,1);
plot(sol1.t, sol1.x(1,:), 'b-o','LineWidth',1.5); hold on;
plot(sol2.t, sol2.x(1,:), 'r-o','LineWidth',1.5);
xlabel('Time (s)'); ylabel('Position');
title('Mass Position Trajectories');
legend('Cart 1','Cart 2');

subplot(2,2,2);
plot(sol1.t, sol1.x(2,:), 'b-o','LineWidth',1.5); hold on;
plot(sol2.t, sol2.x(2,:), 'r-o','LineWidth',1.5);
xlabel('Time (s)'); ylabel('velocity');
title('Mass velocity');
legend('Cart 1','Cart 2');

subplot(2,2,3);
plot(sol1.t(1:end-1), sol1.u, 'b-o','LineWidth',1.5); hold on;
plot(sol2.t(1:end-1), sol2.u, 'r-o','LineWidth',1.5);
xlabel('Time (s)'); ylabel('Control Input');
title('Control Inputs');
legend('Cart 1','Cart 2');

subplot(2,2,4);
plot(sol1.t(1:end-1), sol1.cost, 'b-o','LineWidth',1.5); hold on;
plot(sol2.t(1:end-1), sol2.cost, 'r-o','LineWidth',1.5);
xlabel('Time (s)'); ylabel('Cost');
title('Cost Trajectories');
legend('Cart 1','Cart 2');

%% (Optional) Save the results
%save('M:\DMPC proj\matlab\test\robots\ws_rl3_cdc.mat');                                                     
