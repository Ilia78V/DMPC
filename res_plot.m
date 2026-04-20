folder_path = 'D:\studies\FAU\programming project\report - Copy\resources\figures\reg figures\';

files = {
    'res_noApprox_p20_noAdapt_agent1.fig'
    'res_Approx_p20_noAdapt_agent1.fig'
    'res_regApprox_p20_noAdapt_agent1.fig'
};

files = strcat(folder_path, files);

labels = {
    'without neighborhood approximation'
    'with neighborhood approximation'
    'with regional neighborhood approximation'
};

primalData = cell(numel(files),1);
dualData   = cell(numel(files),1);

for i = 1:numel(files)
    fig = openfig(files{i}, 'invisible');
    ax = findall(fig, 'Type', 'axes');
    ax = ax(1);

    lines = findall(ax, 'Type', 'line');

    for k = 1:length(lines)
        name = get(lines(k), 'DisplayName');

        if contains(lower(name), 'primal')
            primalData{i}.x = get(lines(k), 'XData');
            primalData{i}.y = get(lines(k), 'YData');
        elseif contains(lower(name), 'dual')
            dualData{i}.x = get(lines(k), 'XData');
            dualData{i}.y = get(lines(k), 'YData');
        end
    end

    close(fig);
end

% ---- Primal figure ----
figure;
hold on; grid on;
for i = 1:numel(files)
    plot(primalData{i}.x, primalData{i}.y, 'LineWidth', 1.5);
end
xlabel('ADMM iterations');
ylabel('Primal residual');
legend(labels, 'Location', 'best');
title('Primal residual');

% ---- Dual figure ----
figure;
hold on; grid on;
for i = 1:numel(files)
    plot(dualData{i}.x, dualData{i}.y, 'LineWidth', 1.5);
end
xlabel('ADMM iterations');
ylabel('Dual residual');
legend(labels, 'Location', 'best');
title('Dual residual');

% % Files
% folder_path = 'D:\studies\FAU\programming project\report - Copy\resources\figures\reg figures\';
% 
% files = {
%     'res_noApprox_p20_noAdapt_agent1.fig'
%     'res_Approx_p20_noAdapt_agent1.fig'
%     'res_regApprox_p20_noAdapt_agent1.fig'
% };
% 
% labels = {
%     'without neighborhood approximation'
%     'with neighborhood approximation'
%     'with regional neighborhood approximation'
% };
% 
% % Containers
% primalData = cell(3,1);
% dualData   = cell(3,1);
% 
% for i = 1:3
%     fig = openfig(strcat(folder_path,files{i}), 'invisible');
% 
%     % Find axes
%     ax = findall(fig, 'Type', 'axes');
%     ax = flipud(ax);   % usually fixes MATLAB reversed order
% 
%     % ---- IMPORTANT ----
%     % Assume:
%     % ax(1) = primal residual subplot
%     % ax(2) = dual residual subplot
%     %
%     % If results come out swapped, just swap ax(1) and ax(2) below.
% 
%     % Get line(s) from primal axis
%     linesP = findall(ax(1), 'Type', 'line');
%     linesP = flipud(linesP);
% 
%     % If there is only one line, this still works
%     primalData{i}.x = get(linesP(1), 'XData');
%     primalData{i}.y = get(linesP(1), 'YData');
% 
%     % Get line(s) from dual axis
%     linesD = findall(ax(2), 'Type', 'line');
%     linesD = flipud(linesD);
% 
%     dualData{i}.x = get(linesD(1), 'XData');
%     dualData{i}.y = get(linesD(1), 'YData');
% 
%     close(fig);
% end
% 
% %% Figure 1: Primal residual
% figure;
% hold on; grid on;
% for i = 1:3
%     plot(primalData{i}.x, primalData{i}.y, 'LineWidth', 1.5);
% end
% xlabel('ADMM iterations');
% ylabel('Primal residual');
% legend(labels, 'Location', 'best');
% title('Primal residual');
% set(gca, 'YScale', 'log');   % remove this if you do not want log scale
% 
% %% Figure 2: Dual residual
% figure;
% hold on; grid on;
% for i = 1:3
%     plot(dualData{i}.x, dualData{i}.y, 'LineWidth', 1.5);
% end
% xlabel('ADMM iterations');
% ylabel('Dual residual');
% legend(labels, 'Location', 'best');
% title('Dual residual');
% set(gca, 'YScale', 'log');   % remove this if you do not want log scale