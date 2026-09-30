function u = cb_centralized(x0,c,model,observer)
% Same Euler dynamics, bounds and rectangular stage objective as the DMPC core.
n = c.plant.agents; H = c.controller.horizon_steps; dt = c.controller.sample_time_s;
x = sdpvar(n,H+1,'full'); u = sdpvar(n,H,'full');
% Initial state is a measured constant: admit the same numerical tolerance
% used by validation, without projecting it. Future bounds remain unchanged.
tol = c.constraints.tolerance;
constraints = [x(:,1)==x0, ...
    c.constraints.level_m.min-tol<=x(:,1)<=c.constraints.level_m.max+tol, ...
    c.constraints.level_m.min<=x(:,2:end)<=c.constraints.level_m.max, ...
    c.constraints.input.min<=u<=c.constraints.input.max];
cost = 0;
for k = 1:H
    constraints = [constraints,x(:,k+1)==x(:,k)+dt*cb_dynamics(x(:,k),u(:,k),model)]; %#ok<AGROW>
    error1 = x(:,k)-c.plant.reference_level_m;
    cost = cost+dt*(sum(model.weights.*error1.^2) ...
        +c.plant.input_weight*sum(u(:,k).^2));
end
cost = cost+c.plant.terminal_weight*sum(model.weights.* ...
    (x(:,end)-c.plant.reference_level_m).^2);
% Start near a physically meaningful trajectory, rather than YALMIP's zeros.
assign(x,repmat(x0,1,H+1));
assign(u,zeros(n,H));
t = tic;
d = optimize(constraints,cost,sdpsettings('solver','ipopt','verbose',0, ...
    'usex0',1,'ipopt.max_iter',2000,'ipopt.tol',1e-5,'ipopt.bound_relax_factor',0));
observer.iteration = 1;
observer.recordSolve(0,d.problem,toc(t),d.info);
assert(d.problem==0,'ControlBench:LocalSolveFailed','%s',d.info);
observer.converged = true;
u = value(u(:,1));
end
