function [solver,agents] = cb_build_dmpc(c,model,observer)
n = c.plant.agents; ctl = c.controller; a = ctl.admm;
H = ctl.horizon_steps; dt = ctl.sample_time_s;
agents = cell(1,n);
zero = @(varargin) 0;
for i = 1:n
    data = Agent_data(i,1,1,0,H*dt,H+1,c.plant.initial_level_m, ...
        c.plant.reference_level_m,c.constraints.level_m.min,c.constraints.level_m.max, ...
        c.constraints.input.min,c.constraints.input.max,a.initial_penalty);
    act = model.actuation(i); w = model.weights(i); r = c.plant.input_weight;
    ref = c.plant.reference_level_m; terminal = c.plant.terminal_weight;
    agents{i} = Agent(i,@(x,u)act*u,@(x,t)terminal*w*(x-ref)^2, ...
        @(x,u,t)w*(x-ref)^2+r*u^2,zero,zero,zero,zero,data);
end
for i = 1:n
    for j = 1:n
        if model.edges(i,j) == 0, continue; end
        enabled = strcmp(a.approximation,'all');
        flags = containers.Map({'cost','dynamics','constraints'},{enabled,enabled,enabled});
        data = Neighbor_data(j,1,1,agents{i},a.initial_penalty,flags);
        gain = model.edges(i,j); cubic = model.cubic; linear = model.linear;
        f = @(x,u,y,v)gain*(cubic*(y-x)^3+linear*(y-x));
        neighbor = Neighbor(j,agents{i},true,true,f,zero,zero,zero,zero,zero,zero,data);
        agents{i}.register_neighbors({neighbor});
    end
end
solutions = cellfun(@(agent)Solution(agent,dt),agents,'UniformOutput',false);
solver = ADMM_Solver(ctl.optimizer,agents,solutions,a.max_iterations,a.convergence_tolerance);
solver.ADMM_penaltyAdapt = a.adaptive_penalty;
solver.diagnostics = observer; solver.verbosity = 0;
solver.boundary_tol = c.constraints.tolerance;
end
