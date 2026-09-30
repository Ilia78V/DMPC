function cb_export_diagnostics(observer,out)
local = cell2table(observer.local,'VariableNames', ...
    {'step','iteration','agent_id','problem','duration_s','message'});
writetable(local,fullfile(out,'local_solves.csv'));
residuals = array2table(observer.residual,'VariableNames', ...
    {'step','iteration','agent_id','primal','dual','normalized_primal'});
writetable(residuals,fullfile(out,'residuals.csv'));
end
