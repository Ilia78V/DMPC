%% Plot GRAMPC-D outputs for 4 agents
% Files: WT_chain_1.txt ... WT_chain_4.txt

clear; clc;

nAgents = 4;
data    = cell(nAgents,1);

dl = 251; %data length

for k = 1:nAgents
    fname   = sprintf('orgish/CW_swchApprox_%d.txt', k);

    % Let MATLAB guess the format; it will give Var1, Var2, ...
    T       = readtable(fname);
    %res = readtable('orgish/Residual.txt');

    % Force meaningful names, assuming the column order is:
    % t, x0, u0, v0, cost, debug_cost
    T.Properties.VariableNames = { ...
        't', ...
        'x0', ...
        'u0', ...
        'v0', ...
        'cost', ...
        'debugCost' ...
    };

    data{k} = T;
end
figure;
%% Figure 1 – State x0 vs time
subplot(2,2,1); hold on;
for k = 1:nAgents
    t  = data{k}.t(1:dl);
    x0 = data{k}.x0(1:dl);
    plot(t, x0, 'DisplayName', sprintf('Agent %d', k), 'LineWidth', 2);
end
grid on;
xlabel('Time');
ylabel('State x_0');
title('States of all agents');
legend('show','Location','best');

%% Figure 2 – Control u0 vs time
subplot(2,2,2); hold on;
for k = 1:nAgents
    t  = data{k}.t(1:dl);
    u0 = data{k}.u0(1:dl);
    plot(t, u0, 'DisplayName', sprintf('Agent %d', k), 'LineWidth', 2);
end
grid on;
xlabel('Time');
ylabel('Control u_0');
title('Controls of all agents');
legend('show','Location','best');

%% Figure 3 – Cost vs time
subplot(2,2,3); hold on;
for k = 1:nAgents
    t    = data{k}.t(1:dl);
    J    = data{k}.cost(1:dl);
    plot(t, J, 'DisplayName', sprintf('Agent %d', k), 'LineWidth', 2);
end
grid on;
xlabel('Time');
ylabel('Cost');
title('Costs of all agents');
legend('show','Location','best');

%% Figure 4 – Debug cost vs iteration index (NO time axis)
subplot(2,2,4); hold on;
for k = 1:nAgents
    dbg = data{k}.debugCost;
    it  = 1:numel(dbg);          % iteration index
    plot(it, dbg, 'DisplayName', sprintf('Agent %d', k), 'LineWidth', 2);
end
grid on;
xlabel('Iteration index');
ylabel('Debug cost');
title('Debug costs of all agents');
legend('show','Location','best');
