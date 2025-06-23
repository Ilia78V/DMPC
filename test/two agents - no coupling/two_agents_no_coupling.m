%% Two-Agent Distributed MPC Test
clear; clc;

%% Add paths for YALMIP/OPTI if needed
if true
    originalPath = 'M:\DMPC proj\matlab';
    addpath(genpath('M:\DMPC proj\matlab\YALMIP-master'));
    optiPath = 'M:\DMPC proj\matlab\OPTI-master';
    cd(optiPath);
    run('opti_Install.m');
    cd(originalPath);
end

%% Time Horizon and Discretization Parameters
t0 = 0;
T  = 1;           % Horizon length: 1 second
N  = 6;           % Discretization steps (dt = T/(N-1) = 1/5 = 0.2 sec)
dt = T/(N-1);

%% Simulation Parameters
T_sim     = 5;    % Total simulation time: 5 seconds
dt_sample = 0.5;  % Sampling time: 0.5 seconds

%% Agent Specifications (Simple Single Integrator)
n_x = 1;        % State dimension (scalar state)
n_u = 1;        % Control dimension (scalar control)

% Agent 1: initial state 0, desired reference 1
x0_1   = 0;
x_ref1 = 1;

% Agent 2: initial state 1, desired reference 0
x0_2   = 1;
x_ref2 = 0;

% State and control bounds (same for both agents)
x_min = -10;
x_max = 10;
u_min = -1;
u_max = 1;

rho_init = 10;  % Initial penalty parameter

%% Create Agent Data Objects
agentData1 = Agent_data(1, n_x, n_u, t0, T, N, x0_1, x_ref1, x_min, x_max, u_min, u_max, rho_init);
agentData2 = Agent_data(2, n_x, n_u, t0, T, N, x0_2, x_ref2, x_min, x_max, u_min, u_max, rho_init);

%% Define Cost Weights (Scalars)
Q   = 1;      % Stage state weight
R_c = 0.1;    % Stage control weight
Qf  = 10;     % Terminal state weight

%% Define Dynamics
%   x[k+1] = x[k] + dt * u[k]
f_i = @(x, u) u;

% Coupling dynamics
f_ij = @(x, u, x_neighbor, u_neighbor) 0;

%% Define Cost Functions
% Stage cost: quadratic penalty on state deviation and control effort
l_i = @(x, u, x_ref, t) (x - x_ref).^2 * Q + (u).^2 * R_c;
% Terminal cost: quadratic penalty on final state error
V_i = @(x, x_ref) (x - x_ref).^2 * Qf;

%% Define Dummy Constraint Functions (no additional constraints)
g_i  = @(x, u, t) 0;
h_i  = @(x, u, t) 0;
g_ij = @(x, u, x_neighbor, u_neighbor, t) 0;
h_ij = @(x, u, x_neighbor, u_neighbor, t) 0;

%% Create Agent Objects
% Wrap the cost functions to capture each agent's own reference.
agent1 = Agent(1, f_i, @(x) V_i(x, x_ref1), @(x, u, t) l_i(x, u, x_ref1, t), g_i, h_i, agentData1);
agent2 = Agent(2, f_i, @(x) V_i(x, x_ref2), @(x, u, t) l_i(x, u, x_ref2, t), g_i, h_i, agentData2);

%% Create Neighbor Objects for Coupling (with f_ij set to zero)
neighborData1 = Neighbor_data(1, n_x, n_u, agent2, rho_init);
neighbor1     = Neighbor(1, agent2.id, true, true, f_ij, g_ij, h_ij, neighborData1);

neighborData2 = Neighbor_data(2, n_x, n_u, agent1, rho_init);
neighbor2     = Neighbor(2, agent1.id, true, true, f_ij, g_ij, h_ij, neighborData2);

%% Register Neighbors
agent1.register_neighbors({neighbor2});
agent2.register_neighbors({neighbor1});

%% Create Solution Containers
solution1 = Solution(agent1, dt_sample);
solution2 = Solution(agent2, dt_sample);

%% Create the ADMM Solver for the Two Agents
agents = {agent1, agent2};
max_iterations = 50;
convergence_tolerance = 1e-3;
solver = ADMM_Solver(agents, max_iterations, convergence_tolerance);

%% Simulation Loop
for t = 0:dt_sample:T_sim
    solver.solve();
    solution1.update_solution();
    solution2.update_solution();
    solver.shift(dt_sample);
end

%% Plot the Results
figure;
subplot(3,1,1);
plot(solution1.t, solution1.x, 'b-o','LineWidth',1.5); hold on;
plot(solution2.t, solution2.x, 'r-o','LineWidth',1.5);
xlabel('Time (s)'); ylabel('Position');
title('Agent Position Trajectories');
legend('Agent 1','Agent 2');

subplot(3,1,2);
plot(solution1.t, solution1.u, 'b-o','LineWidth',1.5); hold on;
plot(solution2.t, solution2.u, 'r-o','LineWidth',1.5);
xlabel('Time (s)'); ylabel('Control Input');
title('Agent Control Inputs');
legend('Agent 1','Agent 2');

subplot(3,1,3);
plot(solution1.t, solution1.cost, 'b-o','LineWidth',1.5); hold on;
plot(solution2.t, solution2.cost, 'r-o','LineWidth',1.5);
xlabel('Time (s)'); ylabel('Cost');
title('Cost');
legend('Agent 1','Agent 2');
