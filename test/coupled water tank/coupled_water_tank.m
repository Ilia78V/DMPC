%% Water Tank DMPC Test File in MATLAB
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

%% Time Horizon and Simulation Parameters
t0        = 0;
T         = 2;           % Prediction horizon length [s]
dt        = 0.1;         % Discretization time step [s]
N         = round(T/dt) + 1;  % Number of discretization steps
T_sim     = 6;         % Total simulation time [s]
dt_sample = dt;         % Sampling time [s] (you can adjust as needed)

%% Water Tank Model Parameters
Ai = 0.1;     % Tank area
d  = 0.01;    % Outflow for agent 4

% For agent 1, the inflow is nonzero; for agents 2 and 3, no inflow;
% for agent 4, outflow is d.
% Define model parameters for each agent as: [Ai, ci, di]
params1 = [Ai, 1, 0];      % Agent 1: inflow = 1, outflow = 0
params2 = [Ai, 0, 0];      % Agent 2: no inflow/outflow
params3 = [Ai, 0, 0];      % Agent 3: no inflow/outflow
params4 = [Ai, 0, d];      % Agent 4: outflow = d

% Cost parameters (each as [terminal state weight, stage state weight, stage control weight])
% For agents 1–3, cost is set to zero (except agent 1 penalizes control effort);
% Agent 4 penalizes state deviation.
R_val = 0.1; P_val = 1; Q_val = 1;

cost1 = [0, 0, R_val];     % Agent 1: only control cost
cost2 = [0, 0, 0];         % Agent 2: no cost
cost3 = [0, 0, 0];         % Agent 3: no cost
cost4 = [P_val, Q_val, 0]; % Agent 4: penalize water level error

%% Initial and Desired Conditions (scalar water level)
x0 = 0.5;    % Initial water level
%u0 = 0.0;    % Initial control input
xdes = 2.0;  % Desired water level
%udes = 0.0;  % Desired control (typically zero)

%% Constraint Function (e.g. water level must be below 3)
h_fun = @(x,u,t) x - 3;   % Inequality: x - 3 <= 0
g_fun = @(x,u,t) 0;       % No equality constraint

%% Water Tank Agent Dynamics
% The water tank dynamic (continuous-time) is:
%   x_dot = (ci * u - di) / Ai
% Discretized using Euler:
%   x[k+1] = x[k] + dt * ((ci*u[k] - di) / Ai)
% We define separate function handles for each agent (capturing their model parameters)
f_agent = cell(4,1);
f_agent{1} = @(x,u) (params1(2) * u - params1(3)) / params1(1);
f_agent{2} = @(x,u) (params2(2) * u - params2(3)) / params2(1);
f_agent{3} = @(x,u) (params3(2) * u - params3(3)) / params3(1);
f_agent{4} = @(x,u) (params4(2) * u - params4(3)) / params4(1);

%% Water Tank Cost Functions for Each Agent
% Stage cost: l(x,u) = (state error cost) + (control cost)
% Terminal cost: V(x) = (terminal state cost)
l_cost = cell(4,1);
V_cost = cell(4,1);
% Agent 1: only control cost
l_cost{1} = @(x,u,t) u' * R_val * u;
V_cost{1} = @(x) 0;
% Agents 2 and 3: zero cost (for testing)
l_cost{2} = @(x,u,t) 0;
V_cost{2} = @(x) 0;
l_cost{3} = @(x,u,t) 0;
V_cost{3} = @(x) 0;
% Agent 4: penalize deviation of water level from xdes
l_cost{4} = @(x,u,t) (x - xdes)' * Q_val * (x - xdes);
V_cost{4} = @(x) (x - xdes)' * P_val * (x - xdes);

%% Create Agent Data Objects
% Assuming Agent_data(id, n_x, n_u, t0, T, N, x0, xdes, x_min, x_max, u_min, u_max, rho_init)
rho_init = 1;
n_x = 1; n_u = 1;
% x_min = -Inf; x_max = Inf;
% u_min = -Inf;  u_max = Inf;
x_min = -1e4; x_max = 1e4;
u_min = -1e4;  u_max = 1e4;

agentData1 = Agent_data(1, n_x, n_u, t0, T, N, x0, xdes, x_min, x_max, u_min, u_max, rho_init);
agentData2 = Agent_data(2, n_x, n_u, t0, T, N, x0, xdes, x_min, x_max, u_min, u_max, rho_init);
agentData3 = Agent_data(3, n_x, n_u, t0, T, N, x0, xdes, x_min, x_max, u_min, u_max, rho_init);
agentData4 = Agent_data(4, n_x, n_u, t0, T, N, x0, xdes, x_min, x_max, u_min, u_max, rho_init);

%% Create Agent Objects
% Assuming Agent(id, f_i, V_i, l_i, g_i, h_i, agent_data)
agent1 = Agent(1, f_agent{1}, @(x)V_cost{1}(x), @(x,u,t)l_cost{1}(x,u,t), g_fun, g_fun, @(x,u,t)h_fun(x,u,t), @(x,t)h_fun(x,0,t), agentData1);
agent2 = Agent(2, f_agent{2}, @(x)V_cost{2}(x), @(x,u,t)l_cost{2}(x,u,t), g_fun, g_fun, @(x,u,t)h_fun(x,u,t), @(x,t)h_fun(x,0,t), agentData2);
agent3 = Agent(3, f_agent{3}, @(x)V_cost{3}(x), @(x,u,t)l_cost{3}(x,u,t), g_fun, g_fun, @(x,u,t)h_fun(x,u,t), @(x,t)h_fun(x,0,t), agentData3);
agent4 = Agent(4, f_agent{4}, @(x)V_cost{4}(x), @(x,u,t)l_cost{4}(x,u,t), g_fun, g_fun, @(x,u,t)h_fun(x,u,t), @(x,t)h_fun(x,0,t), agentData4);

%% Define Water Tank Coupling Dynamics
% Coupling dynamics (continuous-time) as defined in the C++ model:
%   Let dx = x_j - x_i.
%   if |dx| < eps, then coupling = poly_param2 * dx^3 + poly_param1 * dx;
%   else coupling = (aij / Ai) * sign(dx) * sqrt(2 * g * |dx|).
eps_val = 0.01;
g_val   = 9.81;
aij     = 0.005;  % pipe area
% Compute polynomial fitting parameters:
q = (aij / Ai) * sqrt(2 * g_val * eps_val);
dqdeps = -(sqrt(2) * aij * g_val) / (2 * Ai * sqrt(g_val * eps_val));
poly_param1 = (3 * q) / (2 * eps_val) - dqdeps / 2;
poly_param2 = dqdeps / (2 * eps_val^2) - q / (2 * eps_val^3);

% Define the coupling function as a function handle:
% f_coupling = @(xi, ui, xj, uj) ...
%     (heaviside(eps_val - abs(xj - xi))) .* (poly_param2 * (xj - xi).^3 + poly_param1 * (xj - xi)) + ...
%     (heaviside(abs(xj - xi) - eps_val)) .* ((aij / Ai) * sign(xj - xi) .* sqrt(2 * g_val * abs(xj - xi)));

% f_coupling = @(xi,ui,xj,uj)...
% (abs(xj-xi) <= eps_val)*(poly_param2*(xj-xi)^3 + poly_param1*(xj-xi)) +...
% (abs(xj-xi) > eps_val)*((aij/Ai)*sign(xj-xi)*sqrt(2*g_val*abs(xj-xi)));


% f_coupling = @(xi,ui,xj,uj) ((aij/Ai)*sign(xj-xi).*sqrt(2*g_val*abs(xj-xi)));
% f_coupling = @(xi,ui,xj,uj) ((xj-xi));

delta = 1e-4;  % Smoothing parameter
k_val = 100;   % Parameter for smooth sign approximation
%w = 1./(1 + exp((abs(xj - xi) - eps_val)/delta));  % Smooth transition weight

% %smooth_sign = tanh(k_val * (xj - xi));              % Smooth approximation of sign
smooth_abs = @(z) (z.^2 )^0.5; %+ 1e-6);
% f_coupling = @(xi,ui,xj,uj) ...
%     (1./(1 + exp((smooth_abs(xj-xi)-eps_val)/delta))) .* (poly_param2 * (xj-xi).^3 + poly_param1 * (xj-xi)) + ...
%     (1 - 1./(1 + exp((smooth_abs(xj-xi)-eps_val)/delta))) .* ((aij/Ai) * tanh(k_val*(xj-xi)) .* (2*g_val*smooth_abs(xj-xi))^0.5);

f_coupling = @(xi,ui,xj,uj) ((aij/Ai) * tanh(k_val*(xj-xi)) .* (2*g_val*smooth_abs(xj-xi))^0.5);

%f_coupling = @(xi,ui,xj,uj) ceil(xj-xi);


%% Dummy Coupling Constraint Functions (set to zero)
g_coup = @(xi,ui,xj,uj,t) 0;
h_coup = @(xi,ui,xj,uj,t) 0;

%% Create Neighbor Objects for Couplings
% For each coupling, create a Neighbor_data object and a Neighbor.
% Format: Neighbor_data(neighbor_id, n_x, n_u, associated_agent, rho_init)
% Then: Neighbor(neighbor_id, associated_agent.id, sending_flag, receiving_flag, f_coupling, g_coup, h_coup, neighbor_data)
%
% According to the C++ file, couplings are:
% (1,2), (1,3); (2,1), (2,4); (3,1), (3,4); (4,2), (4,3)

% Agent 1 couplings:
neighbor_1_from2 = Neighbor(1, agent2.id, true, true, f_coupling, g_coup, g_coup, h_coup, h_coup, Neighbor_data(1, n_x, n_u, agent2, rho_init));
neighbor_1_from3 = Neighbor(1, agent3.id, true, true, f_coupling, g_coup, g_coup, h_coup, h_coup, Neighbor_data(1, n_x, n_u, agent3, rho_init));
% Agent 2 couplings:
neighbor_2_from1 = Neighbor(2, agent1.id, true, true, f_coupling, g_coup, g_coup, h_coup, h_coup, Neighbor_data(2, n_x, n_u, agent1, rho_init));
neighbor_2_from4 = Neighbor(2, agent4.id, true, true, f_coupling, g_coup, g_coup, h_coup, h_coup, Neighbor_data(2, n_x, n_u, agent4, rho_init));
% Agent 3 couplings:
neighbor_3_from1 = Neighbor(3, agent1.id, true, true, f_coupling, g_coup, g_coup, h_coup, h_coup, Neighbor_data(3, n_x, n_u, agent1, rho_init));
neighbor_3_from4 = Neighbor(3, agent4.id, true, true, f_coupling, g_coup, g_coup, h_coup, h_coup, Neighbor_data(3, n_x, n_u, agent4, rho_init));
% Agent 4 couplings:
neighbor_4_from2 = Neighbor(4, agent2.id, true, true, f_coupling, g_coup, g_coup, h_coup, h_coup, Neighbor_data(4, n_x, n_u, agent2, rho_init));
neighbor_4_from3 = Neighbor(4, agent3.id, true, true, f_coupling, g_coup, g_coup, h_coup, h_coup, Neighbor_data(4, n_x, n_u, agent3, rho_init));

%% Register Neighbors with Each Agent
% For Agent 1, register neighbors from agents 2 and 3.
agent1.register_neighbors({neighbor_2_from1, neighbor_3_from1});
% For Agent 2, register neighbors from agents 1 and 4.
agent2.register_neighbors({neighbor_1_from2, neighbor_4_from2});
% For Agent 3, register neighbors from agents 1 and 4.
agent3.register_neighbors({neighbor_1_from3, neighbor_4_from3});
% For Agent 4, register neighbors from agents 2 and 3.
agent4.register_neighbors({neighbor_2_from4, neighbor_3_from4});

%% Create Solution Containers for Each Agent
solution1 = Solution(agent1, dt_sample);
solution2 = Solution(agent2, dt_sample);
solution3 = Solution(agent3, dt_sample);
solution4 = Solution(agent4, dt_sample);

%% Create the ADMM Solver
% Assuming ADMM_Solver takes a cell array of agents, a maximum iteration count, and a convergence tolerance.
max_iterations = 25;
convergence_tolerance = 1e-3;
agents = {agent1, agent2, agent3, agent4};
solutions = {solution1, solution2, solution3, solution4};
solver = ADMM_Solver(agents, solutions, max_iterations, convergence_tolerance);

%% Simulation Loop
for t = 0:dt_sample:T_sim
    solver.solve();
    solver.shift(dt_sample);
end

%% Plot the Results
figure;
subplot(3,1,1);
plot(solution1.t, solution1.x, 'b-o', 'LineWidth',1.5); hold on;
plot(solution2.t, solution2.x, 'r-o', 'LineWidth',1.5);
plot(solution3.t, solution3.x, 'g-o', 'LineWidth',1.5);
plot(solution4.t, solution4.x, 'k-o', 'LineWidth',1.5);
xlabel('Time (s)'); ylabel('Water Level');
title('Water Tank Level Trajectories');
legend('Agent 1','Agent 2','Agent 3','Agent 4');

subplot(3,1,2);
plot(solution1.t(:, 1:end-1), solution1.u, 'b-o', 'LineWidth',1.5); hold on;
plot(solution2.t(:, 1:end-1), solution2.u, 'r-o', 'LineWidth',1.5);
plot(solution3.t(:, 1:end-1), solution3.u, 'g-o', 'LineWidth',1.5);
plot(solution4.t(:, 1:end-1), solution4.u, 'k-o', 'LineWidth',1.5);
xlabel('Time (s)'); ylabel('Control Input');
title('Control Input Trajectories');
legend('Agent 1','Agent 2','Agent 3','Agent 4');

subplot(3,1,3);
plot(solution1.t(:, 1:end-1), solution1.cost, 'b-o', 'LineWidth',1.5); hold on;
plot(solution2.t(:, 1:end-1), solution2.cost, 'r-o', 'LineWidth',1.5);
plot(solution3.t(:, 1:end-1), solution3.cost, 'g-o', 'LineWidth',1.5);
plot(solution4.t(:, 1:end-1), solution4.cost, 'k-o', 'LineWidth',1.5);
xlabel('Time (s)'); ylabel('Cost');
title('Cost Trajectories');
legend('Agent 1','Agent 2','Agent 3','Agent 4');
