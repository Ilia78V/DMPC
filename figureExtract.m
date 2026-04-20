folder_path = 'D:\studies\FAU\programming project\report - Copy\resources\figures\reg figures\';


% Open the .fig without displaying weird callbacks
fig = openfig(strcat(folder_path, 'res_regApprox_p20_noAdapt.fig'), 'invisible');

% Find all axes in the figure
ax = findall(fig, 'Type', 'axes');

% Usually findall returns them in reverse order, so flip if needed
ax = flipud(ax);

% Export each subplot into its own figure
for k = 1:length(ax)
    newFig = figure;
    newAx = copyobj(ax(k), newFig);
    set(newAx, 'Position', get(groot, 'defaultAxesPosition'));
    
    % optional: save each one
    saveas(newFig, sprintf('plot_%d.png', k));
    % savefig(newFig, sprintf('plot_%d.fig', k));
end