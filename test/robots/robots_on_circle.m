%% Two-Mass “Opposite-Side” Test for DMPC Framework
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

%% Physical Parameters
m = 1; k = 1; c = 0.1;

%% Horizon & Simulation
t0        = 0;
T         = 2;      N = 11;    dt = T/(N-1);
T_sim     = 5;      dt_sample = dt;

%% Agent dims
n_x = 2; n_u = 1;

%% Initial/desired states
x0_1   = [1;0];    x_ref1 = [0.5;0];    u_ref1 = -1.5;
x0_2   = [2;0];    x_ref2 = [2.5;0];    u_ref2 =  2.0;

%% Bounds
x_min = [-10; -10];  x_max = [10; 10];
u_min = -5;          u_max =  5;

%% ADMM + approx flags
rho_init  = 0.1;
max_iters = 50;
tol       = 1e-3;
approx = containers.Map('KeyType','char','ValueType','logical');
approx('cost')        = false;   % <-- turn ON neighbor‐cost coupling
approx('dynamics')    = false;
approx('constraints') = false;

%% Create Agent_data
agentData1 = Agent_data(1,n_x,n_u,t0,T,N,x0_1,x_ref1,x_min,x_max,u_min,u_max,rho_init,approx);
agentData2 = Agent_data(2,n_x,n_u,t0,T,N,x0_2,x_ref2,x_min,x_max,u_min,u_max,rho_init,approx);

%% Path parameters for Agent 1
R_circle = 2;        % radius
omega     = pi/5;    % angular speed => period 10 s


%% Quadratic weights
Q  = diag([1,1]);
Rr = 0.1;
Qf = diag([30,30]);

%% Local dynamics (integrator + spring/damper)
f_i1 = @(x,u) [ x(2); u - k*x(1) - c*x(2) ];
f_i2 = @(x,u) [ x(2); u ];

%% Coupling “dynamics” = zero (we only cost‐couple)
v_ij = @(x,u,xn,un) [0;0];

%% Local stage & terminal costs
l1 = @(x,u,t)       (x - x_ref1)'*Q*(x - x_ref1) + (u-u_ref1)'*Rr*(u-u_ref1);
V1 = @(x)           (x - x_ref1)'*Qf*(x - x_ref1);

l2 = @(x,u,t)       (x - x_ref2)'*Q*(x - x_ref2) + (u-u_ref2)'*Rr*(u-u_ref2);
V2 = @(x)           (x - x_ref2)'*Qf*(x - x_ref2);

%% Coupling cost: penalize being “aligned” instead of “opposite”
l_ij  = @(xi,ui,xj,uj)  (xi + xj)'*Q*(xi + xj);
V_ij  = @(xi,xj)        (xi + xj)'*Qf*(xi + xj);

%% No constraints
g_i   = @(x,u,t) 0;  h_i   = @(x,u,t) 0;
g_ij  = @(x,u,xn,un,t) 0;  h_ij  = @(x,u,xn,un,t) 0;
g_i_N = @(x,t) 0;  h_i_N = @(x,t) 0;
g_ij_N= @(x,xn,t) 0; h_ij_N= @(x,xn,t) 0;

%% Build Agents
agent1 = Agent(1, f_i1, @(x)V1(x), ...
               @(x,u,t) l1(x,u,t), ...
               g_i, g_i_N, h_i, h_i_N, agentData1);
agent2 = Agent(2, f_i2, @(x)V2(x), ...
               @(x,u,t) l2(x,u,t), ...
               g_i, g_i_N, h_i, h_i_N, agentData2);

%% Neighbor_data + Neighbor (pure cost‐coupling)
nbrData1 = Neighbor_data(1,n_x,n_u,agent2,rho_init);
nbr1     = Neighbor(1, agent2.id, true, true, ...
                   v_ij, ...    % no dynamic coupling
                   [], [], [], [], ...  % no constraints
                   nbrData1);
% attach cost‐coupling functions into neighbor object for your code to pick up later:
nbr1.l_ij = l_ij;
nbr1.V_ij = V_ij;

nbrData2 = Neighbor_data(2,n_x,n_u,agent1,rho_init);
nbr2     = Neighbor(2, agent1.id, true, true, ...
                   v_ij, [], [], [], [], nbrData2);
nbr2.l_ij = l_ij;
nbr2.V_ij = V_ij;

agent1.register_neighbors({});      % Agent1 has no neighbor‐cost
agent2.register_neighbors({nbr2});  % Agent2 will incur l_ij w.r.t Agent1

%% Solutions & Solver
sol1 = Solution(agent1, dt_sample);
sol2 = Solution(agent2, dt_sample);
solutions = {sol1, sol2};

solver = ADMM_Solver({agent1,agent2}, solutions, ...
                     max_iters, tol, approx);

%% Run DMPC loop
for t = 0:dt_sample:T_sim
  fprintf('t=%.2f  ',t);
  solver.solve();
  solver.shift(dt_sample);
end

%% Plot
figure; hold on; axis equal;
theta = linspace(0,2*pi,200);
plot(R*cos(theta),R*sin(theta),'k--','LineWidth',1);
plot(sol1.x(1,:), sol1.x(2,:),'b-','LineWidth',1.5);
plot(sol2.x(1,:), sol2.x(2,:),'r-','LineWidth',1.5);
legend('Circle','Agent1','Agent2');
xlabel('x'); ylabel('y');
title('Agent1 Tracks Circle • Agent2 Opposite-Side');
