function cb_run(configPath,out)
% Batch entry point: exact sample grid, independent plant, portable artifacts.
cb_setup();
c = jsondecode(fileread(configPath));
rng(c.experiment.seed,'twister');
model = cb_model(c); n = c.plant.agents;
metadata = cb_environment();
metadata.model = model;
metadata.simulation = 'independent_unclipped_euler';
fid = fopen(fullfile(out,'matlab.json'),'w');
assert(fid~=-1,'ControlBench:IO','Cannot write metadata.');
fprintf(fid,'%s',jsonencode(metadata));
fclose(fid);
dt = c.controller.sample_time_s; K = round(c.scenario.duration_s/dt);
assert(abs(K*dt-c.scenario.duration_s)<1e-9,'ControlBench:Config','Invalid duration.');
observer = ControlBenchDiagnostics();
cleanup = onCleanup(@()cb_export_diagnostics(observer,out)); %#ok<NASGU>
x = c.plant.initial_level_m*ones(n,1);
states = zeros((K+1)*n,7); inputs = zeros(K*n,7); steps = cell(K,7);
isDmpc = strcmp(c.controller.type,'matlab_dmpc');
if isDmpc, [solver,agents] = cb_build_dmpc(c,model,observer); end
try
for step = 0:K
    idx = step*n+(1:n);
    states(idx,:) = [step*ones(n,1),step*dt*ones(n,1),(1:n)',x, ...
        c.plant.reference_level_m*ones(n,1),c.constraints.level_m.min*ones(n,1), ...
        c.constraints.level_m.max*ones(n,1)];
    if step == K, break; end
    observer.beginStep(step); timer = tic;
    if isDmpc
        solver.solve();
        u = cellfun(@(a)value(a.data.u(:,1)),agents)';
    else
        u = cb_centralized(x,c,model,observer);
    end
    solveTime = toc(timer);
    assert(all(isfinite(u)),'ControlBench:Nonfinite','Nonfinite control.');
    timer = tic;
    next = x+dt*cb_dynamics(x,u,model);
    plantTime = toc(timer);
    assert(all(isfinite(next)),'ControlBench:Nonfinite','Nonfinite physical state.');
    inputs(idx,:) = [step*ones(n,1),step*dt*ones(n,1),(1:n)',u,u, ...
        c.constraints.input.min*ones(n,1),c.constraints.input.max*ones(n,1)];
    termination = 'iteration_limit';
    if observer.converged, termination = 'converged'; end
    criterion = 'local_optimizer_success';
    if isDmpc, criterion = 'normalized_primal_v1'; end
    steps(step+1,:) = {step,solveTime,plantTime,observer.iteration, ...
        double(observer.converged),termination,criterion};
    if isDmpc && step < K-1
        % Update all initial states only after the simultaneous plant step.
        for i = 1:n, agents{i}.data.initialize(next(i),1); end
        for i = 1:n
            for neighbor = agents{i}.neighbors
                neighbor{1}.data.initialize(1);
            end
        end
    end
    x = next;
end
catch failure
    cb_export_signals(states(1:(step+1)*n,:),inputs(1:step*n,:),steps(1:step,:),out);
    writetable(cell2table({step,'execution_failure',failure.message}, ...
        'VariableNames',{'step','type','message'}),fullfile(out,'events.csv'));
    rethrow(failure);
end
cb_export_signals(states,inputs,steps,out);
writetable(cell2table(cell(0,3),'VariableNames',{'step','type','message'}), ...
    fullfile(out,'events.csv'));
end
