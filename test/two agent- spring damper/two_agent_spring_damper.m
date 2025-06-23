%% Two-Agent Distributed MPC Test for Spring–Damper Coupled Agents
clear; clc;

%% Add required paths
if false
    originalPath = 'M:\DMPC proj\matlab';
    addpath(genpath('M:\DMPC proj\matlab\YALMIP-master'));
    optiPath = 'M:\DMPC proj\matlab\OPTI-master';
    cd(optiPath);
    run('opti_Install.m');
    cd(originalPath);
end

%% Physical Parameters
m = 1;          % Mass (kg)
k = 1;          % Spring constant (N/m)
c = 0.1;        % Damping coefficient (N/(m/s))

%% Time Horizon and Discretization
t0 = 0;
T  = 2;         % Length of the horizon (seconds)
N  = 5;        % Number of discretization steps (dt = T/(N-1))
dt = T/(N-1);

%% Simulation Parameters
T_sim    = 5;    % Total simulation time (seconds)
dt_sample = dt;  % Sampling time (seconds)

%% Agent Specifications
n_x = 2;        % State dimension: [position; velocity]
n_u = 1;        % Control dimension: [force]

% Agent 1: initial state and reference
x0_1   = [0; 0];
x_ref1 = [1; 0];

% Agent 2: initial state and reference
x0_2   = [1; 0];
x_ref2 = [0; 0];

% State and control bounds (same for both agents)
x_min = [-10; -10];  x_max = [10; 10];
u_min = -5;          u_max = 5;
rho_init = 10;       % Initial penalty parameter

%% Create Agent Data Objects
agentData1 = Agent_data(1, n_x, n_u, t0, T, N, x0_1, x_ref1, x_min, x_max, u_min, u_max, rho_init);
agentData2 = Agent_data(2, n_x, n_u, t0, T, N, x0_2, x_ref2, x_min, x_max, u_min, u_max, rho_init);

%% Define Cost Weights (Quadratic)
Q   = diag([1, 1]);       % Stage state weight
R   = 0.1;                % Stage control weight
Qf  = diag([10, 10]);     % Terminal state weight

%% Define Dynamics and Coupling Functions
% Local dynamics for each agent (a double integrator):
%   x1_dot = x2,   x2_dot = u
f_i = @(x,u) [ x(2); u ];

% Coupling dynamics (spring–damper effect from the neighbor):
%   The acceleration is influenced by the difference in positions and velocities.
f_ij = @(x, u, x_neighbor, u_neighbor) [ 0; k*(x_neighbor(1)-x(1)) + c*(x_neighbor(2)-x(2)) ];

%% Define Stage and Terminal Cost Functions (without augmentation)
l_i = @(x,u,x_ref,t) (x - x_ref)'*Q*(x - x_ref) + u'*R*u;
V_i = @(x,x_ref) (x - x_ref)'*Qf*(x - x_ref);

%% Define Constraint Functions (no additional constraints in this demo)
g_i  = @(x,u,t) 0;
h_i  = @(x,u,t) 0;
g_ij = @(x,u,x_neighbor,u_neighbor,t) 0;
h_ij = @(x,u,x_neighbor,u_neighbor,t) 0;

g_i_N  = @(x,t) 0;
h_i_N  = @(x,t) 0;
g_ij_N = @(x,x_neighbor,t) 0;
h_ij_N = @(x,x_neighbor,t) 0;

%% Create Agent Objects
% Wrap the cost functions to capture each agent's own reference.
agent1 = Agent(1, f_i, @(x)V_i(x,x_ref1), @(x,u,t)l_i(x,u,x_ref1,t), g_i, g_i_N, h_i, h_i_N, agentData1);
agent2 = Agent(2, f_i, @(x)V_i(x,x_ref2), @(x,u,t)l_i(x,u,x_ref2,t), g_i, g_i_N, h_i, h_i_N, agentData2);

%% Create Neighbor Objects for Coupling
% For agent 1, the neighbor copy from agent 2:
neighborData1 = Neighbor_data(1, n_x, n_u, agent2, rho_init);
neighbor1 = Neighbor(1, agent2.id, true, true, f_ij, g_ij, g_ij_N, h_ij, h_ij_N, neighborData1);

% For agent 2, the neighbor copy from agent 1:
neighborData2 = Neighbor_data(2, n_x, n_u, agent1, rho_init);
neighbor2 = Neighbor(2, agent1.id, true, true, f_ij, g_ij, g_ij_N, h_ij, h_ij_N, neighborData2);

%% Register Neighbors with Each Agent
agent1.register_neighbors({neighbor2});
agent2.register_neighbors({neighbor1});

%% Create Solution Containers
solution1 = Solution(agent1, dt_sample);
solution2 = Solution(agent2, dt_sample);
solutions = {solution1, solution2};

%% Create the ADMM Solver for the Two Agents
agents = {agent1, agent2};
max_iterations = 100;
convergence_tolerance = 1e-3;
solver = ADMM_Solver(agents, solutions, max_iterations, convergence_tolerance);

%% Run the Distributed MPC Simulation Loop
for t = 0:dt_sample:T_sim
    solver.solve();
    solver.shift(dt_sample);
end

%% Plot the Results
figure;
subplot(3,1,1);
plot(solution1.t, solution1.x(1,:), 'b-o','LineWidth',1.5); hold on;
plot(solution2.t, solution2.x(1,:), 'r-o','LineWidth',1.5);
xlabel('Time (s)'); ylabel('Position');
title('Agent Position Trajectories');
legend('Agent 1','Agent 2');

subplot(3,1,2);
plot(solution1.t(:, 1:end-1), solution1.u, 'b-o','LineWidth',1.5); hold on;
plot(solution2.t(:, 1:end-1), solution2.u, 'r-o','LineWidth',1.5);
xlabel('Time (s)'); ylabel('Control Input');
title('Agent Control Inputs');
legend('Agent 1','Agent 2');

subplot(3,1,3);
plot(solution1.t(:, 1:end-1), solution1.cost, 'b-o','LineWidth',1.5); hold on;
plot(solution2.t(:, 1:end-1), solution2.cost, 'r-o','LineWidth',1.5);
xlabel('Time (s)'); ylabel('Cost');
title('Cost');
legend('Agent 1','Agent 2');