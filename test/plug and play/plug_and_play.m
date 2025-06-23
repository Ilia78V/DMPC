%% SmartGrid DMPC Simulation
% This script sets up a distributed MPC simulation for a smart grid example.
% Two agents are defined with smart grid dynamics (adapted from smartGrid_agentModel.cpp)
% and are coupled via a sinusoidal coupling (from smartGrid_couplingModel.cpp).
%
% The simulation uses ADMM to solve the distributed OCP.
%
% Author: [Your Name]
% Date: [Today's Date]

%% Clear Workspace and Setup Environment
clear; clc; close all;

%% Simulation Parameters
Tsim = 20;            % Total simulation time [s]
T    = 12;            % Prediction horizon [s]
dt   = 0.1;           % Time step [s]
N    = round(T/dt) + 1; % Number of discretization points in the horizon
maxIter = 20;         % Maximum number of ADMM iterations
tol     = 0.0005;     % Convergence tolerance for ADMM

%% Smart Grid Model Parameters
% Common parameters for both agents:
I      = 1;           % Inertia
Omega  = 1;           % Frequency
kappa  = 1e-3;        % Friction term
P_max  = 0.1;         % Maximum power exchange (coupling)

%% Cost Parameters
% Cost parameters for the agent models:
% Format: {terminal_weight_state1, terminal_weight_state2, stage_weight_state1, stage_weight_state2, stage_weight_control}
P_cost = 0.1;
Q_cost = 1;
R_cost = 0.01;

%% Initial and Desired Conditions
x_init = [0; 0];      % Initial state: [x1; x2]
x_des  = [0; 0];      % Desired state
u_init = 0;           % Initial control input
u_des  = 0;           % Desired control input

%% State and Control Bounds
% (Using very generous bounds as in the C++ implementation)
x_min = [-1e6; -1e6];
x_max = [ 1e6;  1e6];
u_min = -1e6;
u_max =  1e6;
rho_init = 1;

%% Create Agent Data Objects
agentData0 = Agent_data(0, 2, 1, 0, T, N, x_init, x_des, x_min, x_max, u_min, u_max, rho_init);
agentData1 = Agent_data(1, 2, 1, 0, T, N, x_init, x_des, x_min, x_max, u_min, u_max, rho_init);

%% Define Agent Dynamics and Cost Functions
% Each agent follows smart grid dynamics. Their differences lie in the parameters:
%   Agent 0: P0 = 0, p = 1
%   Agent 1: P0 = -0.01, p = 0

% Agent 0 Dynamics:
P0_0 = 0;
p0   = 1;
f_i0 = @(x,u) [ x(2);
                (1/(I*Omega))*(p0*u + P0_0 - kappa*Omega^2) - (2*kappa/I)*x(2) ];

% Agent 1 Dynamics:
P0_1 = -0.01;
p1   = 0;
f_i1 = @(x,u) [ x(2);
                (1/(I*Omega))*(p1*u + P0_1 - kappa*Omega^2) - (2*kappa/I)*x(2) ];

% Stage Cost Function:
% Penalizes deviation of the second state and control effort.
l_i = @(x,u,t) Q_cost * (x(2) - x_des(2))^2 + R_cost * u^2;

% Terminal Cost Function:
V_i = @(x) P_cost * (x(2) - x_des(2))^2;

% No additional constraints:
g_i   = @(x,u,t) 0;
g_i_N = @(x,t)   0;
h_i   = @(x,u,t) 0;
h_i_N = @(x,t)   0;

%% Create Agent Objects
agent0 = Agent(0, f_i0, V_i, l_i, g_i, g_i_N, h_i, h_i_N, agentData0);
agent1 = Agent(1, f_i1, V_i, l_i, g_i, g_i_N, h_i, h_i_N, agentData1);

%% Define Coupling (Neighbor) Dynamics and Constraints
% The coupling model applies a sinusoidal coupling to the second state:
%   f_ij = [ 0; (P_max/(I*Omega)) * sin( x1 - x_neighbor1 ) ]
f_ij = @(x,u,x_neighbor,u_neighbor) [ 0;
               (P_max/(I*Omega)) * sin( x(1) - x_neighbor(1) ) ];

% No coupling constraints are defined:
g_ij   = @(x,u,x_neighbor,u_neighbor,t) 0;
h_ij   = @(x,u,x_neighbor,u_neighbor,t) 0;
g_ij_N = @(x,x_neighbor,t) 0;
h_ij_N = @(x,x_neighbor,t) 0;

%% Create Neighbor (Coupling) Objects
% For Agent 0, register Agent 1 as its neighbor:
neighborData0 = Neighbor_data(1, 2, 1, agent1, rho_init);
neighbor0 = Neighbor(1, agent1.id, true, true, f_ij, g_ij, g_ij_N, h_ij, h_ij_N, neighborData0);

% For Agent 1, register Agent 0 as its neighbor:
neighborData1 = Neighbor_data(0, 2, 1, agent0, rho_init);
neighbor1 = Neighbor(0, agent0.id, true, true, f_ij, g_ij, g_ij_N, h_ij, h_ij_N, neighborData1);

%% Register Neighbors with Agents
agent0.register_neighbors({neighbor0});
agent1.register_neighbors({neighbor1});

%% Create Solution Containers
solution0 = Solution(agent0, dt);
solution1 = Solution(agent1, dt);
solutions = {solution0, solution1};

%% Create ADMM Solver
agents = {agent0, agent1};
solver = ADMM_Solver(agents, solutions, maxIter, tol);

%% Run the DMPC Simulation Loop
for t = 0:dt:Tsim
    solver.solve();
    solver.shift(dt);
end

%% Plot Simulation Results
figure;
subplot(3,1,1);
plot(solution0.t, solution0.x(1,:), 'b-o','LineWidth',1.5); hold on;
plot(solution1.t, solution1.x(1,:), 'r-o','LineWidth',1.5);
xlabel('Time [s]'); ylabel('State x_1');
title('State Trajectories');
legend('Agent 0', 'Agent 1');
grid on;

subplot(3,1,2);
plot(solution0.t(1:end-1), solution0.u, 'b-o','LineWidth',1.5); hold on;
plot(solution1.t(1:end-1), solution1.u, 'r-o','LineWidth',1.5);
xlabel('Time [s]'); ylabel('Control Input u');
title('Control Inputs');
legend('Agent 0', 'Agent 1');
grid on;

subplot(3,1,3);
plot(solution0.t(1:end-1), solution0.cost, 'b-o','LineWidth',1.5); hold on;
plot(solution1.t(1:end-1), solution1.cost, 'r-o','LineWidth',1.5);
xlabel('Time [s]'); ylabel('Cost');
title('Cost Trajectories');
legend('Agent 0', 'Agent 1');
grid on;

%% Save Results (Optional)
%save('smartGrid_solution.mat');
