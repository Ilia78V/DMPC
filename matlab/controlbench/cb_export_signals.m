function cb_export_signals(states,inputs,steps,out)
t = array2table(states,'VariableNames', ...
    {'step','time_s','agent_id','value','reference','lower','upper'});
t.state_name = repmat("level",height(t),1); t.unit = repmat("m",height(t),1);
writetable(t,fullfile(out,'states.csv'));
t = array2table(inputs,'VariableNames', ...
    {'step','time_s','agent_id','commanded','applied','lower','upper'});
t.input_name = repmat("flow",height(t),1); t.unit = repmat("m3/s",height(t),1);
writetable(t,fullfile(out,'inputs.csv'));
t = cell2table(steps,'VariableNames',{'step','solve_time_s','plant_time_s', ...
    'iterations','converged','termination','stopping_criterion'});
writetable(t,fullfile(out,'steps.csv'));
end
