%% Two-Mass Spring-Damper Chain Test for DMPC Framework
% Mass 1 is attached to the wall; Mass 2 is attached to Mass 1.
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
m = 1;      % Mass (kg)
k = 1;      % Spring constant (N/m)
c = 0.1;    % Damping coefficient (N/(m/s))

%% Time Horizon and Discretization
t0 = 0;
T_ho  = 2;             % Prediction horizon (seconds)
N  = 11;            % Number of discretization points (dt = T/(N-1))
dt = T_ho/(N-1);

%% Simulation Parameters
T_sim = 5;         % Total simulation time (seconds)
dt_sample = dt;    % Sampling time (seconds)

%% Agent Specifications
n_x = 2;   % Each agent has state [position; velocity]
n_u = 1;   % Control is scalar (force)

% Agent 1 (Mass 1, attached to wall)
x0_1   = [1; 0];         % Initial state: position 1, velocity 0
x_ref1 = [0.5; 0];       % Desired state: move toward 0.5
u_ref1 = [-1.5];

% Agent 2 (Mass 2, attached to Mass 1)
x0_2   = [2; 0];         % Initial state: position 2, velocity 0
x_ref2 = [2.5; 0];       % Desired state: move toward 2.5
u_ref2 = [2];

% State and control bounds
x_min = [-10; -10];  
x_max = [10; 10];
u_min = -5;  
u_max = 5;
rho_init = 50;  % Penalty parameter

% Neighbor approximation
appr = true;
approximation = containers.Map('KeyType', 'char', 'ValueType', 'logical');
approximation('cost') = appr;
approximation('dynamics') = appr;
approximation('constraints') = appr;

%% Create Agent Data Objects
agentData1 = Agent_data(1, n_x, n_u, t0, T_ho, N, x0_1, x_ref1, x_min, x_max, u_min, u_max, rho_init, approximation);
agentData2 = Agent_data(2, n_x, n_u, t0, T_ho, N, x0_2, x_ref2, x_min, x_max, u_min, u_max, rho_init, approximation);

%tunning rho
agentData1.rho_u_i = agentData1.rho_u_i/5;
agentData2.rho_u_i = agentData2.rho_u_i/5;
%% Define Cost Weights (Quadratic)
Q   = diag([1, 1]);       % Stage state weight
R   = 0.1;                % Stage control weight
Qf  = diag([30, 30]);     % Terminal state weight

%% Define Local Dynamics
% Agent 1 dynamics: mass attached to wall with spring/damper forces
f_i1 = @(x,u) [ x(2); u - k*x(1) - c*x(2) ];
% Agent 2 dynamics: free double integrator
f_i2 = @(x,u) [ x(2); u ];

%% Define Coupling Dynamics
% For both agents, define the coupling force based on the relative state.
% For Agent 1, the coupling term is:
%   f_ij^(1) = [0; - k*(x1 - x2) - c*(v1 - v2)]
% For Agent 2, the coupling term is:
%   f_ij^(2) = [0;  k*(x1 - x2) + c*(v1 - v2)]
f_ij1 = @(x, u, x_neighbor, u_neighbor) [ 0; - k*(x(1) - x_neighbor(1)) - c*(x(2) - x_neighbor(2)) ];
f_ij2 = @(x, u, x_neighbor, u_neighbor) [ 0;  k*(x_neighbor(1) - x(1)) + c*(x_neighbor(2) - x(2)) ];

%% Define Stage and Terminal Cost Functions
% l_i = @(x, u, x_ref, t) (x - x_ref)'*Q*(x - x_ref) + u'*R*u;
l_i = @(x, u, x_ref, u_ref, t) (x - x_ref)'*Q*(x - x_ref) + (u - u_ref)'*R*(u - u_ref);
V_i = @(x, x_ref, T) (x - x_ref)'*Qf*(x - x_ref);
l_ij = @(x,u,xn,un,t) 0;
V_ij = @(x,xn,T) 0;
% V_i = @(x, x_ref, u, u_ref) (x - x_ref)'*Qf*(x - x_ref) + (u - u_ref)'*Qf*(u - u_ref);

%% Define Constraint Functions (no additional constraints for this example)
g_i  = @(x, u, t) 0;
h_i  = @(x, u, t) 0;
g_ij = @(x, u, x_neighbor, u_neighbor, t) 0;
h_ij = @(x, u, x_neighbor, u_neighbor, t) 0;
g_i_N = @(x, t) 0;
h_i_N = @(x, t) 0;
g_ij_N = @(x, x_neighbor, t) 0;
h_ij_N = @(x, x_neighbor, t) 0;

%% Create Agent Objects
agent1 = Agent(1, f_i1, @(x,T)V_i(x,x_ref1,T), @(x,u,t) l_i(x,u,x_ref1,u_ref1,t), ...
    g_i, g_i_N, h_i, h_i_N, agentData1);
agent2 = Agent(2, f_i2, @(x,T)V_i(x,x_ref2,T), @(x,u,t) l_i(x,u,x_ref2,u_ref2,t), ...
    g_i, g_i_N, h_i, h_i_N, agentData2);

%% Create Neighbor Objects for Coupling
% For Agent 1: register Agent 2 as neighbor with coupling f_ij1
neighborData1 = Neighbor_data(1, n_x, n_u, agent2, rho_init);
%tuning rho
neighborData1.rho_u_ij = neighborData1.rho_u_ij/5;
neighborData1.rho_u_ji = neighborData1.rho_u_ji/5;

neighbor1 = Neighbor(1, agent2, true, true, f_ij1, g_ij, g_ij_N, h_ij, h_ij_N, V_ij, l_ij, neighborData1);


% For Agent 2: register Agent 1 as neighbor with coupling f_ij2
neighborData2 = Neighbor_data(2, n_x, n_u, agent1, rho_init);
%tuning rho
neighborData2.rho_u_ij = neighborData2.rho_u_ij/5;
neighborData2.rho_u_ji = neighborData2.rho_u_ji/5;

neighbor2 = Neighbor(2, agent1, true, true, f_ij2, g_ij, g_ij_N, h_ij, h_ij_N, V_ij, l_ij, neighborData2);

%% Register Neighbors with Each Agent
agent1.register_neighbors({neighbor2});
agent2.register_neighbors({neighbor1});

%% Create Solution Containers
solution1 = Solution(agent1, dt_sample);
solution2 = Solution(agent2, dt_sample);
solutions = {solution1, solution2};

%% Create the ADMM Solver for the Two Agents
agents = {agent1, agent2};

max_iterations = 50;
convergence_tolerance = 1e-3;

solver = ADMM_Solver('ipopt', agents, solutions, max_iterations, convergence_tolerance, approximation);

agent1.register_solver(solver);
agent2.register_solver(solver);

%% Run the Distributed MPC Simulation Loop
for t = 0:dt_sample:T_sim
    fprintf('t_sim = %.2f\n', t);
    solver.solve();
    solver.shift(dt_sample);
end

%% Plot the Residuals
figure;
colors = 'br';
subplot(2,1,1); hold on;
plot(agentData1.primal_residual, ['-' colors(1) 'o'],'LineWidth',1.5);
ylabel('residual'); title(['agent' num2str(1)]);
plot(agentData1.dual_residual, ['-' colors(2) 'o'],'LineWidth',1.5);
legend('primal','dual');

subplot(2,1,2); hold on;
plot(agentData2.primal_residual, ['-' colors(1) 'o'],'LineWidth',1.5);
ylabel('residual'); title(['agent' num2str(2)]);
plot(agentData2.dual_residual, ['-' colors(2) 'o'],'LineWidth',1.5);
legend('primal','dual');

%% Plot the Results
figure;
subplot(3,1,1);
plot(solution1.t, solution1.x(1,:), 'b-o','LineWidth',1.5); hold on;
plot(solution2.t, solution2.x(1,:), 'r-o','LineWidth',1.5);
xlabel('Time (s)'); ylabel('Position');
title('Mass Position Trajectories');
legend('Mass 1 (Agent 1)','Mass 2 (Agent 2)');

subplot(3,1,2);
plot(solution1.t(:, 1:end-1), solution1.u, 'b-o','LineWidth',1.5); hold on;
plot(solution2.t(:, 1:end-1), solution2.u, 'r-o','LineWidth',1.5);
xlabel('Time (s)'); ylabel('Control Input');
title('Control Inputs');
legend('Mass 1 (Agent 1)','Mass 2 (Agent 2)');

subplot(3,1,3);
plot(solution1.t(:, 1:end-1), solution1.cost, 'b-o','LineWidth',1.5); hold on;
plot(solution2.t(:, 1:end-1), solution2.cost, 'r-o','LineWidth',1.5);
xlabel('Time (s)'); ylabel('Cost');
title('Cost Trajectories');
legend('Mass 1 (Agent 1)','Mass 2 (Agent 2)');

%% (Optional) Save the results
%save('M:\DMPC proj\matlab\test\two agent- spring damper\cost\workspace.mat');