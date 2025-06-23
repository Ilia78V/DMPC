%% Two-Agent Inverted Pendulum DMPC Test Case (Nonlinear Dynamics with Cart Coupling)
% Each agent is a cart–inverted pendulum system.
% Initially, Agent 1’s pendulum is at +15° and Agent 2’s at -10°.
% A spring now couples the carts (i.e. the first state of each agent).
% The goal is to drive both pendulums to the upright position (θ = 0).

clear; clc; close all;

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
M = 1;             % Mass of the cart (kg)
m = 0.1;           % Mass of the pendulum bob (kg)
L = 1;             % Length of the pendulum (m)
g = 9.81;          % Gravitational acceleration (m/s^2)
k_c = 0.5;           % Coupling spring constant (N/m) for the carts

%% Simulation and Horizon Settings
t0 = 0;            % Initial time
T  = 3;            % Prediction horizon (seconds)
N  = 21;            % Number of discretization points (dt = T/(N-1))
dt = T/(N-1);
T_sim = 9;         % Total simulation time (seconds)
dt_sample = 0.1;    % Sampling time

%% Agent Specifications
n_x = 4;           % State dimension: [x; theta; x_dot; theta_dot]
n_u = 1;           % Single control input: force on cart

% Define state bounds
x_min = [-100; -pi; -50; -50];
x_max = [ 100;  pi;  50;  50];
u_min = -20;       % Force bounds (N)
u_max = 20;
rho_init = 100;     % DMPC penalty parameter

%% Initial and Desired States
% Agent 1: initial cart at 0, pendulum at +15° (unstable)
x0_1   = [0; 5*pi/180; 0; 0];
x_ref1 = [0; 0; 0; 0];  % Desired: upright pendulum, cart at 0

% Agent 2: initial cart at 0, pendulum at -10° (unstable)
x0_2   = [0; -3*pi/180; 0; 0];
x_ref2 = [0; 0; 0; 0];  % Desired: upright pendulum, cart at 0

uinit = 0;         % Initial control input
udes  = 0;         % Desired control input

%% Create Agent Data Objects
agentData1 = Agent_data(1, n_x, n_u, t0, T, N, x0_1, x_ref1, x_min, x_max, u_min, u_max, rho_init);
agentData2 = Agent_data(2, n_x, n_u, t0, T, N, x0_2, x_ref2, x_min, x_max, u_min, u_max, rho_init);

%% Define Nonlinear Dynamics Function for the Inverted Pendulum on a Cart
% The dynamics are given by:
%   x_dot(1) = x3,
%   x_dot(2) = x4,
%   x_dot(3) = ( u - m*g*sin(x2)*cos(x2) + m*L*x4^2*sin(x2) ) / (M + m - m*cos(x2)^2),
%   x_dot(4) = ( g*sin(x2) - cos(x2)*( u - m*g*sin(x2)*cos(x2) + m*L*x4^2*sin(x2) )/(M + m - m*cos(x2)^2) ) / L.
f_i = @(x, u) [ x(3);
                x(4);
                ( u - m*g*sin(x(2))*cos(x(2)) + m*L*x(4)^2*sin(x(2)) ) / (M + m - m*cos(x(2))^2 );
                ( g*sin(x(2)) - cos(x(2))*( u - m*g*sin(x(2))*cos(x(2)) + m*L*x(4)^2*sin(x(2)) )/(M + m - m*cos(x(2))^2 ) ) / L ];

% Linearized dynamics for the inverted pendulum on a cart:
%   x1_dot = x3,
%   x2_dot = x4,
%   x3_dot = (1/M)*u - (m*g/M)*x2,
%   x4_dot = ((M+m)*g/(M*L))*x2 - (1/(M*L))*u.
% f_i = @(x, u) [ x(3);
%                 x(4);
%                 (1/M)*u - (m*g/M)*x(2);
%                 ((M+m)*g/(M*L))*x(2) - (1/(M*L))*u ];

%% Define Stage and Terminal Cost Functions
% Penalize deviation from desired state and control effort.
Q   = diag([0, 50, 10, 10]);
R   = 0.1;
Qf  = diag([0, 70, 10, 10]);

l_i = @(x, u, x_ref, t) (x - x_ref)' * Q * (x - x_ref) + u' * R * u;
V_i = @(x, x_ref) (x - x_ref)' * Qf * (x - x_ref);

%% Define Constraint Functions (none additional in this example)
g_i   = @(x, u, t) 0;
h_i   = @(x, u, t) 0;
g_i_N = @(x, t)   0;
h_i_N = @(x, t)   0;

%% Define Coupling Constraint Functions
g_ij   = @(x, u, x_neighbor, u_neighbor, t) 0;
h_ij   = @(x, u, x_neighbor, u_neighbor, t) 0;
g_ij_N = @(x, x_neighbor, t) 0;
h_ij_N = @(x, x_neighbor, t) 0;

%% Define Coupling Function (Cart Coupling)
% The coupling is now between the cart positions (state 1):
%   F_spring = k_c*(x_i(1) - x_j(1))
% This force is applied to the cart acceleration (third state).
f_ij = @(x, u, x_neighbor, u_neighbor) ...
    [ 0;
      0;
      - (k_c/M)*(x(1) - x_neighbor(1));
      0 ];

%% Create Agent Objects Using Nonlinear Dynamics
agent1 = Agent(1, f_i, @(x) V_i(x, x_ref1), @(x,u,t) l_i(x,u,x_ref1,t), ...
    g_i, g_i_N, h_i, h_i_N, agentData1);
agent2 = Agent(2, f_i, @(x) V_i(x, x_ref2), @(x,u,t) l_i(x,u,x_ref2,t), ...
    g_i, g_i_N, h_i, h_i_N, agentData2);

%% Create Neighbor (Coupling) Objects
% For Agent 1, register Agent 2 as neighbor.
neighborData1 = Neighbor_data(1, n_x, n_u, agent2, rho_init);
neighbor1 = Neighbor(1, agent2.id, true, true, f_ij, g_ij, g_ij_N, h_ij, h_ij_N, neighborData1);

% For Agent 2, register Agent 1 as neighbor.
neighborData2 = Neighbor_data(2, n_x, n_u, agent1, rho_init);
neighbor2 = Neighbor(2, agent1.id, true, true, f_ij, g_ij, g_ij_N, h_ij, h_ij_N, neighborData2);

%% Register Neighbors with Each Agent
agent1.register_neighbors({neighbor2});
agent2.register_neighbors({neighbor1});

%% Create Solution Containers
solution1 = Solution(agent1, dt_sample);
solution2 = Solution(agent2, dt_sample);
solutions = {solution1, solution2};

%% Create the ADMM Solver for DMPC
agents = {agent1, agent2};
max_iterations = 40;
convergence_tolerance = 1e-2;
solver = ADMM_Solver(agents, solutions, max_iterations, convergence_tolerance);

%% Run the Distributed MPC Simulation Loop
for t = 0:dt_sample:T_sim
    fprintf('t_sim = %.2f\n', t);
    solver.solve();
    solver.shift(dt_sample);
end

%% Plot the Results
figure;
subplot(3,2,1);
plot(solution1.t, solution1.x(1,:), 'b-o','LineWidth',1.5); hold on;
plot(solution2.t, solution2.x(1,:), 'r-o','LineWidth',1.5);
xlabel('Time (s)'); ylabel('Cart Position (m)');
title('Cart Position Trajectories');
legend('Agent 1','Agent 2');
grid on;

subplot(3,2,2);
plot(solution1.t, solution1.x(2,:)*180/pi, 'b-o','LineWidth',1.5); hold on;
plot(solution2.t, solution2.x(2,:)*180/pi, 'r-o','LineWidth',1.5);
xlabel('Time (s)'); ylabel('Pendulum Angle (°)');
title('Pendulum Angle Trajectories');
legend('Agent 1','Agent 2');
grid on;

subplot(3,2,3);
plot(solution1.t, solution1.x(3,:), 'b-o','LineWidth',1.5); hold on;
plot(solution2.t, solution2.x(3,:), 'r-o','LineWidth',1.5);
xlabel('Time (s)'); ylabel('Velocity');
title('Cart velocity');
legend('Agent 1','Agent 2');
grid on;

subplot(3,2,4);
plot(solution1.t, solution1.x(4,:)*180/pi, 'b-o','LineWidth',1.5); hold on;
plot(solution2.t, solution2.x(4,:)*180/pi, 'r-o','LineWidth',1.5);
xlabel('Time (s)'); ylabel('Angular velocity');
title('Pendulum angular velocity');
legend('Agent 1','Agent 2');
grid on;

subplot(3,2,5);
plot(solution1.t(1:end-1), solution1.u, 'b-o','LineWidth',1.5); hold on;
plot(solution2.t(1:end-1), solution2.u, 'r-o','LineWidth',1.5);
xlabel('Time (s)'); ylabel('Control Force (N)');
title('Control Inputs');
legend('Agent 1','Agent 2');
grid on;

subplot(3,2,6);
plot(solution1.t(:, 1:end-1), solution1.cost, 'b-o','LineWidth',1.5); hold on;
plot(solution2.t(:, 1:end-1), solution2.cost, 'r-o','LineWidth',1.5);
xlabel('Time (s)'); ylabel('Cost');
title('Cost Trajectories');
legend('Agent 1','Agent 2');
grid on;


%% (Optional) Save the results
%save('M:\DMPC proj\matlab\test\Two connected inverted pendulums on carts\taste cases\workspace.mat');
