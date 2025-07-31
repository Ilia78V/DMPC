%% Water Tank DMPC Test (new API)
clear; clc;

%% (Optional) Add YALMIP/IPOPT paths
if false
    %originalPath = 'W:\proj\DMPC proj\matlab';
    originalPath = 'C:\Users\ej62ixyv\proj\DMPC proj\matlab';
    addpath(genpath([originalPath, '\YALMIP-master']));
    optiPath = [originalPath, '\OPTI-master'];
    cd(optiPath);
    run('opti_Install.m');
    cd(originalPath);
end

%% Simulation & horizon
t0        = 0;
T         = 2;          % MPC horizon [s]
dt        = 0.1;        % discretization [s]
N         = round(T/dt)+1;
T_sim     = 6;          % total sim time [s]
dt_sample = dt;

%% Tank parameters
Ai = 0.1;   % cross‐section area
d  = 0.01;  % outflow coefficient for agent 4

params = {
    [Ai, 1, 0],    % agent 1: inflow=1, outflow=0
    [Ai, 1, 0],    % agent 2: no in/out
    [Ai, 1, 0],    % agent 3: no in/out
    [Ai, 1, d]     % agent 4: outflow=d
};

%% Cost weights
R_val = 0.1;  P_val = 1;  Q_val = 1;
cost_params = {
    [0, 0,    R_val],  ... % agent 1: control only
    [0, 0,    0    ],  ... % agent 2: none
    [0, 0,    0    ],  ... % agent 3: none
    [P_val, Q_val, 0 ]     % agent 4: state only
};

x0   = 0.5;    % all tanks start at 0.5 m
xdes = 2.0;    % desired for agent 4

%% Approximation flags
appr = false;
approx = containers.Map({'cost','dynamics','constraints'},{appr,appr,appr});

%% Build Agent_data objects
rho_init = 100;
n_x = 1; n_u = 1;
x_min = -1e3;  x_max = +1e3;
u_min = -1e3;  u_max = +1e3;

agentData = cell(1,4);
for i=1:4
    cp = cost_params{i};
    xr = xdes;  % only agent 4 truly tracks xdes
    agentData{i} = Agent_data( ...
        i, n_x, n_u, t0, T, N, x0, xr, ...
        x_min, x_max, u_min, u_max, rho_init, approx );
    %tunning rho
    % agentData{i}.rho_u_i = agentData{i}.rho_u_i/5;
end

%% Define dynamics and costs
f_agent = cellfun(@(p) @(x,u) (p(2)*u - p(3))/p(1), params, 'UniformOutput',false);
l_cost  = {
    @(x,u,t) Q_val*(x-xdes).^2 + R_val*u^2, ...
    @(x,u,t) Q_val*(x-xdes).^2 + R_val*u^2, ...
    @(x,u,t) Q_val*(x-xdes).^2 + R_val*u^2, ...
    @(x,u,t) Q_val*(x-xdes).^2 + R_val*u^2
};
V_cost  = {
    @(x,Th) P_val*(x-xdes).^2, ...
    @(x,Th) P_val*(x-xdes).^2, ...
    @(x,Th) P_val*(x-xdes).^2, ...
    @(x,Th) P_val*(x-xdes).^2
};

l_ij = @(x,u,xn,un,t) 0;  
V_ij = @(x,xn,T) 0;

%% Create Agent objects
agents = cell(1,4);
for i=1:4
    agents{i} = Agent( ...
        i, ...
        f_agent{i}, ...
        V_cost{i}, ...
        l_cost{i}, ...
        @(x,u,t)0, @(x,T)0, ...  % g_i, g_i_N
        @(x,u,t)x-3, @(x,T)x-3, ... % h_i, h_i_N: enforce x≤3
        agentData{i} );
end

%% Coupling: smooth pipe flow between neighbors
f_coup = @(xi,ui,xj,uj) 0;

g_ij = @(x, u, xn, un, t) 0;
h_ij = @(x, u, xn, un, t) 0;

g_ij_N = @(x, xn, t) 0;
h_ij_N = @(x, xn, t) 0;

%% Construct Neighbor objects
% Define neighbor pairs as [agent_id, neighbor_id]
neighbor_pairs = [1 2; 1 3; 2 1; 2 4; 3 1; 3 4; 4 2; 4 3];
neighbor_data = cell(size(neighbor_pairs,1), 1);
neighbors = cell(size(neighbor_pairs,1), 1);  % to store neighbor objects

for i = 1:size(neighbor_pairs,1)
    neighbor_id   = neighbor_pairs(i,1);
    agent_id = neighbor_pairs(i,2);

    neighbor_data{i} = Neighbor_data(neighbor_id, n_x, n_u, agents{agent_id}, rho_init, approx);
    %tunning rho
    % neighbor_data{i}.rho_u_ij = neighbor_data{i}.rho_u_ij/5;
    % neighbor_data{i}.rho_u_ji = neighbor_data{i}.rho_u_ji/5;

    neighbors{i} = Neighbor(neighbor_id, ...
                            agents{agent_id}, ...
                            true, true, ...
                            f_coup, ...
                            g_ij, g_ij_N, ...
                            h_ij, h_ij_N, ...
                            V_ij, l_ij, ...
                            neighbor_data{i});
    % Register neighbors
    agents{agent_id}.register_neighbors({neighbors{i}});
end

%% Register neighbors
% agents{1}.register_neighbors({neighbor_2_from1, neighbor_3_from1});
% agents{2}.register_neighbors({neighbor_1_from2, neighbor_4_from2});
% agents{3}.register_neighbors({neighbor_1_from3, neighbor_4_from3});
% agents{4}.register_neighbors({neighbor_2_from4, neighbor_3_from4});

%% Create solutions & solver
sols = cellfun(@(a) Solution(a, dt_sample), agents,'UniformOutput',false);
solver = ADMM_Solver('ipopt', agents, sols, 30, 1e-2, approx);

%% Main loop
for t_sim = 0:dt_sample:T_sim
    fprintf('t = %.2f  ',t_sim);
    solver.solve();
    solver.shift(dt_sample);
end

%% Plot residuals
figure;
colors = 'br';
for i=1:4
    subplot(4,1,i); hold on;
    plot(agentData{i}.primal_residual, ['-' colors(1) 'o'],'LineWidth',1.5);
    ylabel('residual'); title(['agent' num2str(i)]);

    plot(agentData{i}.dual_residual, ['-' colors(2) 'o'],'LineWidth',1.5);
end
legend('primal','dual');

%% Plot
figure;
% 1) Levels
subplot(3,1,1); hold on;
colors = 'brgk';
for i=1:4
    plot(sols{i}.t, sols{i}.x, ['-' colors(i) 'o'],'LineWidth',1.5);
end
ylabel('Level'); title('Tank levels');
legend('1','2','3','4');

% 2) Inputs
subplot(3,1,2); hold on;
for i=1:4
    plot(sols{i}.t(1:end-1), sols{i}.u, ['-' colors(i) 'o'],'LineWidth',1.5);
end
ylabel('u'); title('Control inputs');

% 3) Costs
subplot(3,1,3); hold on;
for i=1:4
    plot(sols{i}.t(1:end-1), sols{i}.cost, ['-' colors(i) 'o'],'LineWidth',1.5);
end
ylabel('Cost'); xlabel('Time (s)');
title('Cost per agent');

%% (Optional) Save the results
%save('W:\proj\DMPC proj\matlab\test\coupled water tank\New folder\w_decoup_poly9_it30_rho100.mat');