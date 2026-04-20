classdef ADMM_Solver
    properties (Access = public)
        
        optimizer;

        % Agent and data
        agents;
        agent_map;

        % Solutions
        solutions;
        solution_map;

        % Solver parameters
        max_iterations;
        convergence_tolerance;
        slackWeight
        boundary_tol

        % Residuals for convergence check
        primal_residual;
        dual_residual;
        N_pr;

        approximation;
        % % ADMM states
        % x_opt;
        % u_opt;
        % x_neighbors_opt;
        % u_neighbors_opt;

        % Penalty parameters
        ADMM_penaltyAdapt;
        ADMM_PenaltyIncreaseFactor;
        ADMM_PenaltyDecreaseFactor;
        ADMM_PenaltyMin;
        ADMM_PenaltyMax;
        ADMM_PenaltyInit;
    end

    methods
        %% Constructor
        function obj = ADMM_Solver(optimizer, agents, solutions, max_iterations, convergence_tolerance)
            obj.optimizer = optimizer;
            
            obj.agents = agents;
            obj.solutions = solutions;
            obj.max_iterations = max_iterations;
            obj.convergence_tolerance = convergence_tolerance;
            obj.slackWeight = 1e8;   % tune this
            obj.boundary_tol= 1e-2;% tune this

            obj.agent_map = containers.Map('KeyType', 'double', 'ValueType', 'any');
            for i = 1:length(obj.agents)
                obj.agent_map(obj.agents{i}.id) = obj.agents{i};
            end

            obj.solution_map = containers.Map('KeyType', 'double', 'ValueType', 'any');
            for i = 1:length(obj.solutions)
                obj.solution_map(obj.solutions{i}.id) = obj.solutions{i};
            end

            % obj.approximation = approximation;
            obj.set_boarder_agents();

            % obj.approximation = containers.Map('KeyType', 'char', 'ValueType', 'logical');
            % obj.approximation('cost') = false;
            % obj.approximation('dynamics') = false;
            % obj.approximation('constraints') = false;

            % Residuals for convergence check
            obj.N_pr = zeros(1, length(agents));
            obj = obj.count_N_pr();

            % Penalty parameters
            obj.ADMM_penaltyAdapt = false;
            obj.ADMM_PenaltyIncreaseFactor = 1.5;
            obj.ADMM_PenaltyDecreaseFactor = 0.75;
            obj.ADMM_PenaltyMin            = 1e-4;
            obj.ADMM_PenaltyMax            = 1e4;
            obj.ADMM_PenaltyInit           = 1;

        end

        %% Set border agents
        function set_boarder_agents(obj)
            for agent = obj.agents
                d = agent{1}.data;
                d.border = 0;
                if d.approximation("dynamics")
                    if size(agent{1}.dyn_approx_neighbors, 2) ~= size(agent{1}.neighbors, 2)
                        d.border = 1;
                    end
                else
                    for n = agent{1}.neighbors
                        nd = n{1}.data;
                        if obj.agent_map(nd.id).data.approximation("dynamics")
                            d.border = 2;
                            break;
                        end
                    end
                end
            end
        end

        %% ADMM solver function
        function solve(obj)
            for q = 1:obj.max_iterations
                % Step 1: Compute local variables
                obj.compute_local_variables();

                % Step 2: Send local copies to sending neighbors
                obj.send_local_copies();

                % Step 3: Compute coupling variables
                obj.compute_coupling_variables();

                % Step 4: Send coupling variables to receiving neighbors
                obj.send_coupling_variables();

                % Step 5: Compute Lagrange multipliers
                obj.compute_multipliers();

                % Step 6: Send Lagrange multipliers to sending neighbors
                obj.send_multipliers();

                % Step 7: Compute residuals and adapt penalty parameters
                obj.update_residual();

                % Step 8: Send penalty parameters to sending neighbors
                obj.send_penalties();

                % Step 9: Check convergence
                if obj.check_convergence()
                    obj.update_previous_data();
                    fprintf('ADMM solver converged at %.2f\n', q)
                    break;
                elseif q == obj.max_iterations
                    warning('ADMM solver terminated\n')
                end
                obj.update_previous_data();
            end
        end

        %% Compute local variables (Step 1)
        function compute_local_variables(obj)
            for i = 1:length(obj.agents)
                agent = obj.agents{i};
                d = agent.data;

                % s_x = sdpvar(size(d.x,1), size(d.x,2));
                % s0 = zeros(size(d.x));

                % -----Define constraints-----
                constraints = [];
                constraints = [constraints, d.x(:,1) == d.x0]; % Initial state condition
                if ~isempty(d.u0)
                    constraints = [constraints, d.u(:,1) == d.u0]; % Initial control condition
                end
                
                if d.approximation('dynamics') || d.border
                    for neighbor = agent.sending_neighbors
                        nd = neighbor{1}.data;
                        constraints = [constraints, nd.x_ji(:,1) == obj.agent_map(neighbor{1}.id).data.x0]; % Initial state condition
                    end
                end

                for k = 1:d.N-1
                    % Dynamics constraint (Euler discretization)
                    dynamics = agent.f_i(d.x(:,k), d.u(:,k));
                    for j = 1:length(agent.sending_neighbors)
                        nd = agent.sending_neighbors{j}.data;
                        dynamics = dynamics + agent.sending_neighbors{j}.f_ij(d.x(:,k), d.u(:,k), nd.x_ji(:,k), nd.u_ji(:,k));
                    end

                    constraints = [constraints, d.x(:,k+1) == d.x(:,k) + d.dt * dynamics];
                    
                    %%%%%%%%%%%%%%%% Approximation of neighbor dynamics %%%%%%%%%%%%%%%%                    
                    if d.approximation('dynamics') || d.border
                        for neighbor = agent.sending_neighbors
                            nd = neighbor{1}.data;
                            neighbor_dynamics = obj.agent_map(neighbor{1}.id).f_i(nd.x_ji(:,k), nd.u_ji(:,k));
                            neighbor_dynamics = neighbor_dynamics + obj.agent_map(neighbor{1}.id).neighbor_map(agent.id).f_ij(nd.x_ji(:,k), nd.u_ji(:,k), d.x(:,k), d.u(:,k));
                            neighbor_dynamics = neighbor_dynamics + nd.v_ji(:,k);
                            constraints = [constraints, nd.x_ji(:,k+1) == nd.x_ji(:,k) + d.dt * neighbor_dynamics];
                        end
                    end
                    
                    %%%%%%%%%%%%%%%% Approximation of neighbor dynamics %%%%%%%%%%%%%%%%
                    % compute the external influence
                    if d.approximation('dynamics') || d.border
                        for neighbor = agent.receiving_neighbors
                            nd = neighbor{1}.data;

                            v = sdpvar(size(nd.v_i,1), 1);
                            for se_neighbor = agent.sending_neighbors
                                sd = se_neighbor{1}.data;
                                if neighbor{1}.id ~= se_neighbor{1}.id
                                    v = v + se_neighbor{1}.f_ij(d.x(:,k), d.u(:,k), sd.x_ji(:,k), sd.u_ji(:,k));
                                end
                            end
                            constraints = [constraints, nd.v_i(:,k) == v];
                        end
                    end
                    

                    % State and control bounds
                    % constraints = [constraints, s_x(:,k) >= s0(:,k)];
                    % constraints = [constraints, d.x_min - s_x(:,k) <= d.x(:,k) <= d.x_max + s_x(:,k)];
                    constraints = [constraints, d.x_min <= d.x(:,k) <= d.x_max ];
                    constraints = [constraints, d.u_min <= d.u(:,k) <= d.u_max];
                    
                    % Local equality and inequality constraints
                    constraints = [constraints, agent.g_i(d.x(:,k), d.u(:,k), d.t(k)) == 0];
                    constraints = [constraints, agent.h_i(d.x(:,k), d.u(:,k), d.t(k)) <= 0];

                    %%%%%%%%%%%%%%%% Approximation of neighbor constraints %%%%%%%%%%%%%%%%
                    for neighbor = agent.const_approx_neighbors
                        nd = neighbor{1}.data;
                        ag = obj.agent_map(neighbor{1}.id);
                        nda = ag.data; % Neighbor data as an agent
                        % State and control bounds
                        constraints = [constraints, nda.u_min <= nd.u_ji(:,k) <= nda.u_max];
                        % Local equality and inequality constraints
                        constraints = [constraints, ag.g_i(nd.x_ji(:,k), nd.u_ji(:,k), d.t(k)) == 0];
                        constraints = [constraints, ag.h_i(nd.x_ji(:,k), nd.u_ji(:,k), d.t(k)) <= 0];
                        % Neighbor equality and inequality constraints
                        constraints = [constraints, ag.neighbor_map(agent.id).g_ij(nd.x_ji(:,k), nd.u_ji(:,k), d.x(:,k), d.u(:,k), d.t(k)) == 0];
                        constraints = [constraints, ag.neighbor_map(agent.id).h_ij(nd.x_ji(:,k), nd.u_ji(:,k), d.x(:,k), d.u(:,k), d.t(k)) <= 0];
                    end
                    
                    % Neighbor equality and inequality constraints
                    for j = 1:length(agent.sending_neighbors)
                        nd = agent.sending_neighbors{j}.data;
                        constraints = [constraints, agent.sending_neighbors{j}.g_ij(d.x(:,k), d.u(:,k), nd.x_ji(:,k), nd.u_ji(:,k), d.t(k)) == 0];
                        constraints = [constraints, agent.sending_neighbors{j}.h_ij(d.x(:,k), d.u(:,k), nd.x_ji(:,k), nd.u_ji(:,k), d.t(k)) <= 0];
                    end
                end

                % FINAL STATE BOUNDS
                % constraints = [constraints, s_x(:,d.N) >= s0(:,d.N)];
                % constraints = [constraints, d.x_min - s_x(:,d.N) <= d.x(:,d.N) <= d.x_max + s_x(:,d.N)];
                constraints = [constraints, d.x_min <= d.x(:,d.N) <= d.x_max]; 
                constraints = [constraints, agent.g_i_N(d.x(:,d.N), d.t(d.N)) == 0];
                constraints = [constraints, agent.h_i_N(d.x(:,d.N), d.t(d.N)) <= 0];
                for j = 1:length(agent.sending_neighbors)
                    nd = agent.sending_neighbors{j}.data;
                    constraints = [constraints, agent.sending_neighbors{j}.g_ij_N(d.x(:,d.N), nd.x_ji(:,d.N), d.t(d.N)) == 0];
                    constraints = [constraints, agent.sending_neighbors{j}.h_ij_N(d.x(:,d.N), nd.x_ji(:,d.N), d.t(d.N)) <= 0];
                end
                
                    %%%%%%%%%%%%%%%% Approximation of neighbor constraints %%%%%%%%%%%%%%%%
                for neighbor = agent.const_approx_neighbors   %sending_neighbors
                    nd = neighbor{1}.data;
                    ag = obj.agent_map(neighbor{1}.id);
                    % Local equality and inequality constraints
                    constraints = [constraints, ag.g_i_N(nd.x_ji(:,d.N), d.t(d.N)) == 0];
                    constraints = [constraints, ag.h_i_N(nd.x_ji(:,d.N), d.t(d.N)) <= 0];
                    % Neighbor equality and inequality constraints
                    constraints = [constraints, ag.neighbor_map(agent.id).g_ij_N(nd.x_ji(:,d.N), d.x(:,d.N), d.t(d.N)) == 0];
                    constraints = [constraints, ag.neighbor_map(agent.id).h_ij_N(nd.x_ji(:,d.N), d.x(:,d.N), d.t(d.N)) <= 0];
                end

                % Construct cost
                cost = obj.cost_constructor(agent);
                %%%%%%%%%%%%%%%% Approximation of neighbor cost %%%%%%%%%%%%%%%%
                if ~isempty(agent.cost_approx_neighbors)
                    cost = cost * (1/(1+length(agent.cost_approx_neighbors)));
                    for neighbor = agent.cost_approx_neighbors
                        cost = cost + (1/(1+length(obj.agent_map(neighbor{1}.id).cost_approx_neighbors))) * obj.cost_constructor(neighbor{1});
                    end
                end
                % cost = cost + obj.slackWeight * sum(s_x, 'all');

                % % -----Define the ADMM augmented cost function-----
                % cost = agent.V_i(d.x(:,d.N)); % Terminal cost
                % for k = 1:d.N-1
                %     % Local stage cost
                %     cost = cost + d.dt * agent.l_i(d.x(:,k), d.u(:,k), d.t(k));
                % 
                %     % Local coupling terms and Lagrange multipliers (scaled by d.dt for integral form)
                %     z_local = [d.x(:,k); d.u(:,k)];
                %     z_coupling = [d.z_x(:,k); d.z_u(:,k)];
                %     mu_local = [d.mu_x(:,k); d.mu_u(:,k)];
                %     cost = cost + d.dt * mu_local' * (z_coupling - z_local);
                %     cost = cost + (d.dt/2) * (z_coupling - z_local)' * d.C_i * (z_coupling - z_local);
                % 
                %     % Neighbor coupling terms (scaled by d.dt for integral form)
                %     for j = 1:length(agent.sending_neighbors)
                %         nd = agent.sending_neighbors{j}.data;
                %         z_neighbor = [nd.x_ji(:,k); nd.u_ji(:,k)];
                %         z_coupling_neighbor = [nd.z_x_j(:,k); nd.z_u_j(:,k)];
                %         mu_neighbor = [nd.mu_x_ji(:,k); nd.mu_u_ji(:,k)];
                %         cost = cost + d.dt * mu_neighbor' * (z_coupling_neighbor - z_neighbor);
                %         cost = cost + (d.dt/2) * (z_coupling_neighbor - z_neighbor)' * nd.C_ji * (z_coupling_neighbor - z_neighbor);
                %     end
                % end
                % % Final state costs
                % % Local coupling terms and Lagrange multipliers (scaled by d.dt for integral form)
                % z_local = [d.x(:,d.N)];
                % z_coupling = [d.z_x(:,d.N)];
                % mu_local = [d.mu_x(:,d.N)];
                % cost = cost + d.dt * mu_local' * (z_coupling - z_local);
                % cost = cost + (d.dt/2) * (z_coupling - z_local)' * diag(d.rho_x_i) * (z_coupling - z_local);
                % 
                % % Neighbor coupling terms (scaled by d.dt for integral form)
                % for j = 1:length(agent.sending_neighbors)
                %     nd = agent.sending_neighbors{j}.data;
                %     z_neighbor = [nd.x_ji(:,d.N)];
                %     z_coupling_neighbor = [nd.z_x_j(:,d.N)];
                %     mu_neighbor = [nd.mu_x_ji(:,d.N)];
                %     cost = cost + d.dt * mu_neighbor' * (z_coupling_neighbor - z_neighbor);
                %     cost = cost + (d.dt/2) * (z_coupling_neighbor - z_neighbor)' * diag(nd.rho_x_ji) * (z_coupling_neighbor - z_neighbor);
                % end

                % % -----Update the previous data-----
                % agent.previous_data.x = d.x; 
                % agent.previous_data.u = d.u;
                % 
                % for j = 1:length(agent.sending_neighbors)
                %     agent.sending_neighbors{j}.previous_data.x_ji = agent.sending_neighbors{j}.data.x_ji;
                %     agent.sending_neighbors{j}.previous_data.u_ji = agent.sending_neighbors{j}.data.u_ji;
                % end

                % -----Solve the optimization problem-----

                options = sdpsettings('solver', obj.optimizer, 'verbose', 1, 'debug', 0);

                
                options.ipopt.max_iter = 2000;           % Set max iterations
                options.ipopt.tol = 1e-5;               % Set convergence tolerance
               
                
                sol = optimize(constraints, cost, options);
                % sol = optimize(constraints);

                options.ipopt.warm_start_init_point = 'yes';
   
                if sol.problem == 0
                    disp('Solver successfully found an optimal solution.');
                elseif sol.problem == 1
                    disp('Solver failed: Infeasible problem.');
                elseif sol.problem == 2
                    disp('Solver failed: Unbounded problem.');
                elseif sol.problem == 3
                    disp('Solver failed: Numerical issues.');
                else
                    disp(['Solver returned an unknown status: ', yalmiperror(sol.problem)]);
                end

                
                % -----Update cost-----
                % cost = sdpvar(1,1);
                % assign(cost, [0]);
                d.cost = [d.cost, value(cost)];

                % Approximation of neighbor dynamics (external influence)
                % if obj.approximation('dynamics')
                %     for neighbor = agent.receiving_neighbors
                %         nd = neighbor{1}.data;
                %         for se_neighbor = agent.sending_neighbors
                %             if neighbor{1}.id ~= se_neighbor{1}.id
                %                 nd.v_i = nd.v_i + value(neighbor{1}.f_ij(nd.x_ji, nd.u_ji, d.x, d.u));
                %             end
                %         end
                %     end
                % end


                if d.approximation('dynamics') || d.border
                    for neighbor = agent.sending_neighbors
                        nd = neighbor{1}.data;
                        if isnan(value(nd.v_ji(:, end)))
                            assign(nd.v_ji(:, end), zeros(size(nd.v_ji(:, end))));
                        end
                    end
                end
                
                % x_opt = value(d.x);
                % u_opt = value(d.u);
                % x_neighbors_opt = cellfun(@value, d.x_ji, 'UniformOutput', false);
                % u_neighbors_opt = cellfun(@value, d.u_ji, 'UniformOutput', false);
                % 
                % agent.update_agentState(x_opt, u_opt);
                % for j = 1:length(agent.sending_neighbors)
                %     agent.sending_neighbors{j}.update_local_copies(x_neighbors_opt{j}, u_neighbors_opt{j});
                % end
            end
        end

        %% Cost Constructor
        function cost = cost_constructor(obj, agent)
            if isa(agent, 'Agent')
               isagent = true;
               d = agent.data;
               x = d.x;
               u = d.u;
               border = d.border;
            elseif isa(agent, 'Neighbor')
               isagent = false;
               d = agent.data;
               x = d.x_ji;
               u = d.u_ji;
               agent = obj.agent_map(agent.id);
               border = agent.data.border;
            end
            
            if d.approximation('dynamics') || border
                % -----Define the ADMM augmented cost function-----
                cost = agent.V_i(x(:,d.N), d.N); % Terminal cost
                for neighbor = agent.neighbors
                    nd = neighbor{1}.data;
                    cost = cost + neighbor{1}.V_ij(x(:,d.N), nd.x_ji(:,d.N), d.N);
                end
                    
                for k = 1:d.N-1
                    % Local stage cost
                    cost = cost + d.dt * agent.l_i(x(:,k), u(:,k), d.t(k));
                    for neighbor = agent.neighbors
                        nd = neighbor{1}.data;
                        cost = cost + d.dt * neighbor{1}.l_ij(x(:,k), u(:,k), nd.x_ji(:,k), nd.u_ji(:,k), d.t(k));
                    end

                    if isagent
                        % Local coupling terms and Lagrange multipliers (scaled by d.dt for integral form)
                        cost = cost + d.dt * d.mu_u(:,k)' * (d.z_u(:,k) - d.u(:,k));
                        cost = cost + (d.dt/2) * (d.z_u(:,k) - d.u(:,k))' * diag(d.rho_u_i(:,k)) * (d.z_u(:,k) - d.u(:,k));
                        for neighbor = agent.receiving_neighbors
                            nd = neighbor{1}.data; 
                            cost = cost + d.dt * nd.mu_v_i(:,k)' * (nd.z_v_i(:,k) - nd.v_i(:,k));
                            cost = cost + (d.dt/2) * (nd.z_v_i(:,k) - nd.v_i(:,k))' * diag(nd.rho_v_i(:,k)) * (nd.z_v_i(:,k) - nd.v_i(:,k));
                        end
                        
                        % Neighbor coupling terms (scaled by d.dt for integral form)
                        for neighbor = agent.sending_neighbors
                            nd = neighbor{1}.data;
                            z_neighbor = [nd.u_ji(:,k); nd.v_ji(:,k)];
                            z_coupling_neighbor = [nd.z_u_j(:,k); nd.z_v_j(:,k)];
                            mu_neighbor = [nd.mu_u_ji(:,k); nd.mu_v_ji(:,k)];
                            cost = cost + d.dt * mu_neighbor' * (z_coupling_neighbor - z_neighbor);
                            cost = cost + (d.dt/2) * (z_coupling_neighbor - z_neighbor)' * diag([nd.rho_u_ji(:,k); nd.rho_v_ji(:,k)]) * (z_coupling_neighbor - z_neighbor);
                        end
                    end
    
                end
    
                if isagent
                    % Final state costs
                    % Local coupling terms and Lagrange multipliers (scaled by d.dt for integral form)
                    for neighbor = agent.receiving_neighbors
                        nd = neighbor{1}.data; 
                        cost = cost + d.dt * nd.mu_v_i(:,d.N)' * (nd.z_v_i(:,d.N) - nd.v_i(:,d.N));
                        cost = cost + (d.dt/2) * (nd.z_v_i(:,d.N) - nd.v_i(:,d.N))' * diag(nd.rho_v_i(:,d.N)) * (nd.z_v_i(:,d.N) - nd.v_i(:,d.N));
                    end

                    % Neighbor coupling terms (scaled by d.dt for integral form)
                    for neighbor = agent.sending_neighbors
                        nd = neighbor{1}.data;                         
                        cost = cost + d.dt * nd.mu_v_ji(:,d.N)' * (nd.z_v_j(:,d.N) - nd.v_ji(:,d.N));
                        cost = cost + (d.dt/2) * (nd.z_v_j(:,d.N) - nd.v_ji(:,d.N))' * diag(nd.rho_v_ji(:,d.N)) * (nd.z_v_j(:,d.N) - nd.v_ji(:,d.N));
                    end
                end
            
            else
                % DEFAULT CASE (without neighborhood approximation)
            
                % -----Define the ADMM augmented cost function-----
                cost = agent.V_i(x(:,d.N), d.N); % Terminal cost
                for neighbor = agent.neighbors
                    % d = obj.agent_map(neighbor{1}.id).data;
                    nd = neighbor{1}.data;
                    cost = cost + neighbor{1}.V_ij(x(:,d.N), nd.x_ji(:,d.N), d.N);
                end

                for k = 1:d.N-1
                    % Local stage cost
                    cost = cost + d.dt * agent.l_i(x(:,k), u(:,k), d.t(k));
                    for neighbor = agent.neighbors
                        nd = neighbor{1}.data;
                        cost = cost + d.dt * neighbor{1}.l_ij(x(:,k), u(:,k), nd.x_ji(:,k), nd.u_ji(:,k), d.t(k));
                    end
                    
                    if isagent
                        % Local coupling terms and Lagrange multipliers (scaled by d.dt for integral form)
                        z_local = [x(:,k); u(:,k)];
                        z_coupling = [d.z_x(:,k); d.z_u(:,k)];
                        mu_local = [d.mu_x(:,k); d.mu_u(:,k)];
                        cost = cost + d.dt * mu_local' * (z_coupling - z_local);
                        cost = cost + (d.dt/2) * (z_coupling - z_local)' * diag([d.rho_x_i(:,k); d.rho_u_i(:,k)]) * (z_coupling - z_local);
                        
                        % Neighbor coupling terms (scaled by d.dt for integral form)
                        for j = 1:length(agent.sending_neighbors)
                            nd = agent.sending_neighbors{j}.data;
                            z_neighbor = [nd.x_ji(:,k); nd.u_ji(:,k)];
                            z_coupling_neighbor = [nd.z_x_j(:,k); nd.z_u_j(:,k)];
                            mu_neighbor = [nd.mu_x_ji(:,k); nd.mu_u_ji(:,k)];
                            cost = cost + d.dt * mu_neighbor' * (z_coupling_neighbor - z_neighbor);
                            cost = cost + (d.dt/2) * (z_coupling_neighbor - z_neighbor)' * diag([nd.rho_x_ji(:,k); nd.rho_u_ji(:,k)]) * (z_coupling_neighbor - z_neighbor);
                        end
                    end
    
                end
    
                if isagent
                    % Final state costs
                    % Local coupling terms and Lagrange multipliers (scaled by d.dt for integral form)
                    z_local = [x(:,d.N)];
                    z_coupling = [d.z_x(:,d.N)];
                    mu_local = [d.mu_x(:,d.N)];
                    cost = cost + d.dt * mu_local' * (z_coupling - z_local);
                    cost = cost + (d.dt/2) * (z_coupling - z_local)' * diag(d.rho_x_i(:,d.N)) * (z_coupling - z_local);
                    
                    % Neighbor coupling terms (scaled by d.dt for integral form)
                    for j = 1:length(agent.sending_neighbors)
                        nd = agent.sending_neighbors{j}.data;
                        z_neighbor = [nd.x_ji(:,d.N)];
                        z_coupling_neighbor = [nd.z_x_j(:,d.N)];
                        mu_neighbor = [nd.mu_x_ji(:,d.N)];
                        cost = cost + d.dt * mu_neighbor' * (z_coupling_neighbor - z_neighbor);
                        cost = cost + (d.dt/2) * (z_coupling_neighbor - z_neighbor)' * diag(nd.rho_x_ji(:,d.N)) * (z_coupling_neighbor - z_neighbor);
                    end
                end
            end

        end

        %% Send local copies to sending neighbors (Step 2)
        function send_local_copies(obj)
            for agent = obj.agents
                d = agent{1}.data;
                for neighbor = agent{1}.sending_neighbors
                    nd = neighbor{1}.data;
                    ag = obj.agent_map(neighbor{1}.id).neighbor_map(agent{1}.id);
                    ad = obj.agent_map(neighbor{1}.id).data;
                    ag.data.u_ij = value(neighbor{1}.data.u_ji);
                    
                    if d.border == 0
                        if d.approximation('dynamics')
                            ag.data.v_ij = value(neighbor{1}.data.v_ji);
                        else
                            ag.data.x_ij = value(neighbor{1}.data.x_ji);
                        end

                    elseif d.border == 1
                        if ad.border
                            ag.data.x_ij = value(neighbor{1}.data.x_ji);
                            ag.data.v_ij = value(neighbor{1}.data.v_ji);
                        else
                            ag.data.v_ij = value(neighbor{1}.data.v_ji);
                        end

                    elseif d.border == 2
                        if ad.border
                            ag.data.x_ij = value(neighbor{1}.data.x_ji);
                            ag.data.v_ij = value(neighbor{1}.data.v_ji);
                        else
                            ag.data.x_ij = value(neighbor{1}.data.x_ji);                            
                        end
                    end

                    % if obj.approximation('dynamics') || 
                    %     ag.data.v_ij = value(neighbor{1}.data.v_ji);
                    % else
                    %     ag.data.x_ij = value(neighbor{1}.data.x_ji);
                    % end
                end
            end
        end

        %% Compute coupling variables (Step 3)
        function compute_coupling_variables(obj)
            for agent = obj.agents
                d = agent{1}.data;
                
                %%%%%%%% input coupling variables %%%%%%%%
                z_local_u = value(d.u) - (1 ./ d.rho_u_i) .* d.mu_u;
                z_neighbor_u = zeros(size(z_local_u));
                for neighbor = agent{1}.receiving_neighbors
                    nd = neighbor{1}.data;
                    z_neighbor_u = z_neighbor_u + nd.u_ij - (1 ./ nd.rho_u_ij) .* nd.mu_u_ij;
                end
                d.z_u = (1 / (1 + length(agent{1}.receiving_neighbors))) * (z_local_u + z_neighbor_u);

                %%%%%%%% external influence coupling variables %%%%%%%%
                if d.approximation('dynamics') || d.border
                    z_local_v = cell(size(agent{1}.receiving_neighbors));
                    for i = 1:length(agent{1}.receiving_neighbors)
                        nd = agent{1}.receiving_neighbors{i}.data;
                        z_local_v{i} = value(nd.v_i) - (1 ./ nd.rho_v_i) .* nd.mu_v_i;
                    end
    
                    % Compute neighbor contribution
                    z_neighbor_v = cellfun(@(m) zeros(size(z_local_v{1})), z_local_v, 'UniformOutput', false);
                    cell(size(z_local_v));
                    for i = 1:length(agent{1}.receiving_neighbors)
                        nd = agent{1}.receiving_neighbors{i}.data;
                        delta_v = nd.v_ij - (1 ./ nd.rho_v_ij) .* nd.mu_v_ij;
                        z_neighbor_v = cellfun(@(z) z + delta_v, z_neighbor_v, 'UniformOutput', false);
                    end
                    
                    % Compute final coupling variable update
                    for i = 1:length(agent{1}.receiving_neighbors)
                        nd = agent{1}.receiving_neighbors{i}.data;
                        nd.z_v_i = (1 / (1 + length(agent{1}.receiving_neighbors))) * (z_local_v{i} + z_neighbor_v{i});
                    end
                
                end
                
                %%%%%%%% state coupling variables %%%%%%%%
                if ~d.approximation('dynamics') || d.border
                    z_local_x = value(d.x) - (1 ./ d.rho_x_i) .* d.mu_x;
    
                    % Compute neighbor contribution
                    z_neighbor_x = zeros(size(z_local_x));
                    for neighbor = agent{1}.receiving_neighbors
                        nd = neighbor{1}.data;
                        z_neighbor_x = z_neighbor_x + nd.x_ij - (1 ./ nd.rho_x_ij) .* nd.mu_x_ij;
                    end
                    
                    % Compute final coupling variable update
                    d.z_x = (1 / (1 + length(agent{1}.receiving_neighbors))) * (z_local_x + z_neighbor_x);
                end

            end
        end

        %% Send coupling variables to receiving neighbors (Step 4)
        function send_coupling_variables(obj)
            for agent = obj.agents
                d = agent{1}.data;
                for neighbor = agent{1}.receiving_neighbors
                    nd = neighbor{1}.data;
                    ag = obj.agent_map(neighbor{1}.id).neighbor_map(agent{1}.id);
                    ad = obj.agent_map(neighbor{1}.id).data;
                    ag.data.z_u_j = agent{1}.data.z_u;
                    
                    if d.border == 0
                        if d.approximation('dynamics')
                            ag.data.z_v_j = neighbor{1}.data.z_v_i;
                        else
                            ag.data.z_x_j = agent{1}.data.z_x;
                        end

                    elseif d.border == 1
                        if ad.border
                            ag.data.z_x_j = agent{1}.data.z_x;
                            ag.data.z_v_j = neighbor{1}.data.z_v_i;
                        else
                            ag.data.z_v_j = neighbor{1}.data.z_v_i;                            
                        end

                    elseif d.border == 2
                        if ad.border
                            ag.data.z_x_j = agent{1}.data.z_x;
                            ag.data.z_v_j = neighbor{1}.data.z_v_i;
                        else
                            ag.data.z_x_j = agent{1}.data.z_x;                            
                        end
                    end
                end
            end
        end

        %% Compute Lagrange multipliers (Step 5)
        function compute_multipliers(obj)
            for agent = obj.agents
                d = agent{1}.data;
                
                %%%%%%%% input Lagrange multipliers %%%%%%%%
                d.mu_u = d.mu_u + (d.rho_u_i) .* (d.z_u - value(d.u));                
                for neighbor = agent{1}.sending_neighbors
                    nd = neighbor{1}.data;
                    nd.mu_u_ji = nd.mu_u_ji + (nd.rho_u_ji) .* (nd.z_u_j - value(nd.u_ji));
                end
        
                %%%%%%%% external influence Lagrange multipliers %%%%%%%%
                if d.approximation('dynamics') || d.border
                    for i = 1:length(agent{1}.receiving_neighbors)
                        nd = agent{1}.receiving_neighbors{i}.data;
                        nd.mu_v_i = nd.mu_v_i + (nd.rho_v_i) .* (nd.z_v_i - value(nd.v_i));
                    end
                    
                    % Neighbor contribution
                    for neighbor = agent{1}.sending_neighbors
                        nd = neighbor{1}.data;
                        nd.mu_v_ji = nd.mu_v_ji + (nd.rho_v_ji) .* (nd.z_v_j - value(nd.v_ji));
                    end
                end
        
                %%%%%%%% state Lagrange multipliers %%%%%%%%
                if ~d.approximation('dynamics') || d.border
                    d.mu_x = d.mu_x + (d.rho_x_i) .* (d.z_x - value(d.x));
                    
                    % Neighbor contribution
                    for neighbor = agent{1}.sending_neighbors
                        nd = neighbor{1}.data;
                        nd.mu_x_ji = nd.mu_x_ji + (nd.rho_x_ji) .* (nd.z_x_j - value(nd.x_ji));
                    end
                end
            end
        end

        %% Send Lagrange multipliers to sending neighbors (Step 6)
        function send_multipliers(obj)
            for agent = obj.agents
                d = agent{1}.data;
                for neighbor = agent{1}.sending_neighbors
                    nd = neighbor{1}.data;
                    ag = obj.agent_map(neighbor{1}.id).neighbor_map(agent{1}.id);
                    ad = obj.agent_map(neighbor{1}.id).data;
                    ag.data.mu_u_ij = neighbor{1}.data.mu_u_ji;
                    
                    if d.border == 0
                        if d.approximation('dynamics')
                            ag.data.mu_v_ij = neighbor{1}.data.mu_v_ji;
                        else
                            ag.data.mu_x_ij = neighbor{1}.data.mu_x_ji;
                        end

                    elseif d.border == 1
                        if ad.border
                            ag.data.mu_x_ij = neighbor{1}.data.mu_x_ji;
                            ag.data.mu_v_ij = neighbor{1}.data.mu_v_ji;
                        else
                            ag.data.mu_v_ij = neighbor{1}.data.mu_v_ji;                            
                        end

                    elseif d.border == 2
                        if ad.border
                            ag.data.mu_x_ij = neighbor{1}.data.mu_x_ji;
                            ag.data.mu_v_ij = neighbor{1}.data.mu_v_ji;
                        else
                            ag.data.mu_x_ij = neighbor{1}.data.mu_x_ji;                            
                        end
                    end
                end
            end
        end

        %% Check convergence (Step 7)
        function is_converged = check_convergence(obj) 
            is_converged = true;
            for i = 1: length(obj.agents)
                agent = obj.agents{i};
                d = agent.data;
                pd = agent.previous_data;

                % if d.approximation('dynamics')
                %     error_vec = vecnorm([(d.z_u - pd.z_u); (d.mu_u - pd.mu_u)], 2, 2); % Computes the 2-norm (Euclidean norm) of each row
                %     for neighbor = agent{1}.sending_neighbors
                %         nd = neighbor{1}.data;
                %         npd = neighbor{1}.previous_data;
                %         [error_vec; vecnorm([(nd.z_v_i - npd.z_v_i); (nd.mu_v_i - npd.mu_v_i)], 2, 2)];
                %     end
                %     error = norm(error_vec);
                % else
                %     error = norm([(d.z_x(:,1:end-1) - pd.z_x(:,1:end-1)); ...
                %                   (d.z_u - pd.z_u); ...
                %                   (d.mu_x(:,1:end-1) - pd.mu_x(:,1:end-1)); ...
                %                   (d.mu_u - pd.mu_u)]);
                % end 
                % is_converged = is_converged && (error <= obj.convergence_tolerance);
                is_converged = is_converged && ((d.primal_residual(end)/(obj.N_pr(i)^0.5)) <= obj.convergence_tolerance);
            end
        end
        %% Update previous data
        function update_previous_data(obj)
            for agent = obj.agents
                agent{1}.previous_data = copy(agent{1}.data);
                for neighbor = agent{1}.neighbors
                    neighbor{1}.previous_data = copy(neighbor{1}.data);
                end
            end
        end

        %% Calculate primal and dual residual/ Update penalty parameters;
        function update_residual(obj)
            for agent = obj.agents
                d = agent{1}.data;
                pd = agent{1}.previous_data;  
                
                if d.approximation('dynamics') || d.border
                    % ===== Input =====
                    %primal residual 
                    local_pr_u = abs(value(d.u) - d.z_u);
                    ADMM_local_pr_u = norm(local_pr_u, "fro");                    
                    %dual residual
                    local_dr_u = abs(pd.rho_u_i .* (d.z_u - pd.z_u));
                    ADMM_local_dr_u = norm(local_dr_u, "fro");
                    % Update rho
                    if obj.ADMM_penaltyAdapt, d.rho_u_i = adaptPenaltyParameter(obj, local_pr_u, local_dr_u, d.rho_u_i); end

                    % ===== External influence =====
                    %primal residual 
                    local_pr_v = cell(1, length(agent{1}.receiving_neighbors));
                    ADMM_local_pr_v = zeros(length(agent{1}.receiving_neighbors),1);
                    for i = 1: length(agent{1}.receiving_neighbors)
                        nd = agent{1}.receiving_neighbors{i}.data;
                        local_pr_v{i} = abs(value(nd.v_i) - nd.z_v_i);
                        ADMM_local_pr_v(i) = norm(local_pr_v{i}, "fro");
                    end

                    %dual residual
                    local_dr_v = cell(1, length(agent{1}.receiving_neighbors));
                    ADMM_local_dr_v = zeros(length(agent{1}.receiving_neighbors),1);
                    for i = 1: length(agent{1}.receiving_neighbors)
                        nd = agent{1}.receiving_neighbors{i}.data;
                        npd = agent{1}.receiving_neighbors{i}.previous_data;
                        local_dr_v {i} = abs(npd.rho_v_i .* (nd.z_v_i - npd.z_v_i));
                        ADMM_local_dr_v(i) = norm(local_dr_v{i}, "fro");
                        % Update rho
                        if obj.ADMM_penaltyAdapt, nd.rho_v_i = adaptPenaltyParameter(obj, local_pr_v{i}, local_dr_v{i}, nd.rho_v_i); end
                    end
                    
                    % initial = true;
                    % **Ensure neighbor variables are initialized** before the loop
                    neighbor_pr_v = cell(1, length(agent{1}.sending_neighbors));
                    neighbor_pr_u = cell(1, length(agent{1}.sending_neighbors));
                    neighbor_dr_v = cell(1, length(agent{1}.sending_neighbors));
                    neighbor_dr_u = cell(1, length(agent{1}.sending_neighbors));

                    ADMM_neighbor_pr_v = zeros(length(agent{1}.sending_neighbors),1);
                    ADMM_neighbor_pr_u = zeros(length(agent{1}.sending_neighbors),1);
                    ADMM_neighbor_dr_v = zeros(length(agent{1}.sending_neighbors),1);
                    ADMM_neighbor_dr_u = zeros(length(agent{1}.sending_neighbors),1);
    
                    % Neighbor contribution
                    for i = 1: length(agent{1}.sending_neighbors)
                        nd = agent{1}.sending_neighbors{i}.data;
                        npd = agent{1}.sending_neighbors{i}.previous_data;
                        
                        % if initial == true
                        %     neighbor_pr_x = zeros(size(nd.z_x_j));
                        %     neighbor_pr_u = zeros(size(nd.z_u_j));
                        %     neighbor_dr_x = zeros(size(nd.z_x_j));
                        %     neighbor_dr_u = zeros(size(nd.z_u_j));
                        %     initial = false;
                        % end

                        % ===== Input =====
                        %primal residual
                        neighbor_pr_u{i} = abs(value(nd.u_ji) - nd.z_u_j);
                        ADMM_neighbor_pr_u(i) = norm(neighbor_pr_u{i}, "fro");
                        %dual residual
                        neighbor_dr_u{i} = abs(npd.rho_u_ji .* (nd.z_u_j - npd.z_u_j));
                        ADMM_neighbor_dr_u(i) = norm(neighbor_dr_u{i}, "fro");
                        % Update rho
                        if obj.ADMM_penaltyAdapt, nd.rho_u_ji = adaptPenaltyParameter(obj, neighbor_pr_u{i}, neighbor_dr_u{i}, nd.rho_u_ji); end
                        
                        % ===== External influence =====
                        %primal residual
                        neighbor_pr_v{i} = abs(value(nd.v_ji) - nd.z_v_j);
                        ADMM_neighbor_pr_v(i) = norm(neighbor_pr_v{i}, "fro");
                        %dual residual
                        neighbor_dr_v{i} = abs(npd.rho_v_ji .* (nd.z_v_j - npd.z_v_j));
                        ADMM_neighbor_dr_v(i) = norm(neighbor_dr_v{i}, "fro");
                        % Update rho
                        if obj.ADMM_penaltyAdapt, nd.rho_v_ji = adaptPenaltyParameter(obj, neighbor_pr_v{i}, neighbor_dr_v{i}, nd.rho_v_ji); end
                    end
    
                    % neighbor_pr_x = norm(vecnorm(neighbor_pr_x, 2, 1), 1)/d.N;
                    % neighbor_pr_u = norm(vecnorm(neighbor_pr_u, 2, 1), 1)/(d.N - 1);
                    % neighbor_dr_x = norm(vecnorm(neighbor_dr_x, 2, 1), 1)/d.N;
                    % neighbor_dr_u = norm(vecnorm(neighbor_dr_u, 2, 1), 1)/(d.N - 1);
    
                    d.primal_residual = [d.primal_residual, norm([ADMM_local_pr_v; ADMM_local_pr_u; ADMM_neighbor_pr_v; ADMM_neighbor_pr_u], 2)];
                    d.dual_residual = [d.dual_residual, norm([ADMM_local_dr_v; ADMM_local_dr_u; ADMM_neighbor_dr_v; ADMM_neighbor_dr_u], 2)];                 
                else

                    %DEFAULT CASE
                    
                    % ===== Input =====
                    %primal residual
                    local_pr_u = abs(value(d.u) - d.z_u);
                    ADMM_local_pr_u = norm(local_pr_u, "fro");
                    %dual residual
                    local_dr_u = abs(pd.rho_u_i .* (d.z_u - pd.z_u));
                    ADMM_local_dr_u = norm(local_dr_u, "fro");                  
                    % Update rho
                    if obj.ADMM_penaltyAdapt, d.rho_u_i = adaptPenaltyParameter(obj, local_pr_u, local_dr_u, d.rho_u_i); end

                    % ===== State =====
                    %primal residual
                    local_pr_x = abs(value(d.x) - d.z_x);
                    ADMM_local_pr_x = norm(local_pr_x, "fro");
                    %dual residual
                    local_dr_x = abs(pd.rho_x_i .* (d.z_x - pd.z_x));
                    ADMM_local_dr_x = norm(local_dr_x, "fro");
                    % Update rho
                    if obj.ADMM_penaltyAdapt, d.rho_x_i = adaptPenaltyParameter(obj, local_pr_x, local_dr_x, d.rho_x_i); end

                    % initial = true;
                    % % **Ensure neighbor variables are initialized** before the loop
                    % neighbor_pr_x = 0;
                    % neighbor_pr_u = 0;
                    % neighbor_dr_x = 0;
                    % neighbor_dr_u = 0;
                    neighbor_pr_x = cell(1, length(agent{1}.sending_neighbors));
                    neighbor_pr_u = cell(1, length(agent{1}.sending_neighbors));
                    neighbor_dr_x = cell(1, length(agent{1}.sending_neighbors));
                    neighbor_dr_u = cell(1, length(agent{1}.sending_neighbors));

                    ADMM_neighbor_pr_x = zeros(length(agent{1}.sending_neighbors),1);
                    ADMM_neighbor_pr_u = zeros(length(agent{1}.sending_neighbors),1);
                    ADMM_neighbor_dr_x = zeros(length(agent{1}.sending_neighbors),1);
                    ADMM_neighbor_dr_u = zeros(length(agent{1}.sending_neighbors),1);
    
                    % Neighbor contribution
                    for i = 1: length(agent{1}.sending_neighbors)
                        neighbor = agent{1}.sending_neighbors{i};
                        nd = neighbor.data;
                        npd = neighbor.previous_data;
                        
                        % if initial == true
                        %     neighbor_pr_x = zeros(size(nd.z_x_j));
                        %     neighbor_pr_u = zeros(size(nd.z_u_j));
                        %     neighbor_dr_x = zeros(size(nd.z_x_j));
                        %     neighbor_dr_u = zeros(size(nd.z_u_j));
                        %     initial = false;
                        % end
    
                        %  %primal residual
                        %  neighbor_pr_x = neighbor_pr_x + (value(nd.x_ji) - nd.z_x_j);
                        %  neighbor_pr_u = neighbor_pr_u + (value(nd.u_ji) - nd.z_u_j);
                        % 
                        % %dual residual
                        %  neighbor_dr_x = neighbor_dr_x + diag(npd.rho_x_ji) * (nd.z_x_j - npd.z_x_j);
                        %  neighbor_dr_u = neighbor_dr_u + diag(npd.rho_u_ji) * (nd.z_u_j - npd.z_u_j); 


                        % ===== Input =====
                        %primal residual
                        neighbor_pr_u{i} = abs(value(nd.u_ji) - nd.z_u_j);
                        ADMM_neighbor_pr_u(i) = norm(neighbor_pr_u{i}, "fro");
                        %dual residual
                        neighbor_dr_u{i} = abs(npd.rho_u_ji .* (nd.z_u_j - npd.z_u_j));
                        ADMM_neighbor_dr_u(i) = norm(neighbor_dr_u{i}, "fro");
                        % Update rho
                        if obj.ADMM_penaltyAdapt, nd.rho_u_ji = adaptPenaltyParameter(obj, neighbor_pr_u{i}, neighbor_dr_u{i}, nd.rho_u_ji); end                        
                        
                        % ===== State =====
                        %primal residual
                        neighbor_pr_x{i} = abs(value(nd.x_ji) - nd.z_x_j);
                        ADMM_neighbor_pr_x(i) = norm(neighbor_pr_x{i}, "fro");
                        %dual residual
                        neighbor_dr_x{i} = abs(npd.rho_x_ji .* (nd.z_x_j - npd.z_x_j));
                        ADMM_neighbor_dr_x(i) = norm(neighbor_dr_x{i}, "fro");
                        % Update rho
                        if obj.ADMM_penaltyAdapt, nd.rho_x_ji = adaptPenaltyParameter(obj, neighbor_pr_x{i}, neighbor_dr_x{i}, nd.rho_x_ji); end                        

                        % neighbor_pr_u(i) = norm(vecnorm( (value(nd.u_ji) - nd.z_u_j) , 2, 1), 1)/(d.N - 1);
                        % neighbor_pr_v(i) = norm(vecnorm( (value(nd.v_ji) - nd.z_v_j) , 2, 1), 1)/d.N;
                    end
                    % 
                    % neighbor_pr_x = norm(vecnorm(neighbor_pr_x, 2, 1), 1)/d.N;
                    % neighbor_pr_u = norm(vecnorm(neighbor_pr_u, 2, 1), 1)/(d.N - 1);
                    % neighbor_dr_x = norm(vecnorm(neighbor_dr_x, 2, 1), 1)/d.N;
                    % neighbor_dr_u = norm(vecnorm(neighbor_dr_u, 2, 1), 1)/(d.N - 1);
    
                    d.primal_residual = [d.primal_residual, norm([ADMM_local_pr_x; ADMM_local_pr_u; ADMM_neighbor_pr_x; ADMM_neighbor_pr_u], 2)];
                    d.dual_residual = [d.dual_residual, norm([ADMM_local_dr_x; ADMM_local_dr_u; ADMM_neighbor_dr_x; ADMM_neighbor_dr_u], 2)];
                end
            end
        end

        %% Count N_primalResidual
        function obj = count_N_pr(obj)
            for i = 1: length(obj.agents)
                agent = obj.agents{i};
                d = agent.data;

                if d.approximation('dynamics') || d.border
                    % ===== Input =====
                    obj.N_pr(i) = obj.N_pr(i) + (d.n_u * (d.N-1));

                    % ===== External influence =====
                    for j = 1: length(agent.receiving_neighbors)
                        obj.N_pr(i) = obj.N_pr(i) + (d.n_x * (d.N));
                    end
    
                    % Neighbor contribution
                    for j = 1: length(agent.sending_neighbors)
                        nd = agent.sending_neighbors{j}.data;
                        % ===== Input =====
                        obj.N_pr(i) = obj.N_pr(i) + (nd.n_u * (d.N-1));

                        % ===== External influence =====
                        obj.N_pr(i) = obj.N_pr(i) + (nd.n_x * (d.N));
                    end
             
                else

                    %DEFAULT CASE
                    
                    % ===== Input =====
                    obj.N_pr(i) = obj.N_pr(i) + (d.n_u * (d.N-1));

                    % ===== State =====
                    obj.N_pr(i) = obj.N_pr(i) + (d.n_x * (d.N));

                    for j = 1: length(agent.sending_neighbors)
                        nd = agent.sending_neighbors{j}.data;
                        % ===== Input =====
                        obj.N_pr(i) = obj.N_pr(i) + (nd.n_u * (d.N-1));
                 
                        % ===== State =====
                        obj.N_pr(i) = obj.N_pr(i) + (nd.n_x * (d.N));
                    end
                end
            end
        end

        %% Adapt penalty parameter
        function adaptedPenalty = adaptPenaltyParameter(obj, primal_residuum, dual_residuum, penalty)
            % Elementwise version of
            % SolverLocalADMM::adaptPenaltyParameter(primal, dual, penalty)
            %
            % All inputs are matrices of identical size.
        
            % ADMM constants (can also be taken from obj properties if you stored them there)
            incFactor = obj.ADMM_PenaltyIncreaseFactor;  % 1.5
            decFactor = obj.ADMM_PenaltyDecreaseFactor;  % 0.75
            penMin    = obj.ADMM_PenaltyMin;             % 1e-4
            penMax    = obj.ADMM_PenaltyMax;             % 1e4
        
            % ---------- compute factor elementwise ----------
            % default factor = 1
            factor = ones(size(dual_residuum));
        
            % where dual_residuum > 1e-10, use primal/dual
            mask = dual_residuum > 1e-10;
            factor(mask) = primal_residuum(mask) ./ dual_residuum(mask);
        
            % clamp factor between decrease and increase factors, elementwise
            factor = min(factor, incFactor);  % elementwise min with scalar
            factor = max(factor, decFactor);  % elementwise max with scalar
        
            % ---------- scale penalty elementwise ----------
            adaptedPenalty = factor .* penalty;
        
            % clamp penalty between global min/max, elementwise
            adaptedPenalty = min(adaptedPenalty, penMax);
            adaptedPenalty = max(adaptedPenalty, penMin);
        end

        %% Send Lagrange penalty parameters to sending neighbors (Step 6)
        function send_penalties(obj)
            for agent = obj.agents
                d = agent{1}.data;
                for neighbor = agent{1}.sending_neighbors
                    nd = neighbor{1}.data;
                    ag = obj.agent_map(neighbor{1}.id).neighbor_map(agent{1}.id);
                    ad = obj.agent_map(neighbor{1}.id).data;
                    ag.data.rho_u_ij = neighbor{1}.data.rho_u_ji;
                    
                    if d.border == 0
                        if d.approximation('dynamics')
                            ag.data.rho_v_ij = neighbor{1}.data.rho_v_ji;
                        else
                            ag.data.rho_x_ij = neighbor{1}.data.rho_x_ji;
                        end

                    elseif d.border == 1
                        if ad.border
                            ag.data.rho_x_ij = neighbor{1}.data.rho_x_ji;
                            ag.data.rho_v_ij = neighbor{1}.data.rho_v_ji;
                        else
                            ag.data.rho_v_ij = neighbor{1}.data.rho_v_ji;                            
                        end

                    elseif d.border == 2
                        if ad.border
                            ag.data.rho_x_ij = neighbor{1}.data.rho_x_ji;
                            ag.data.rho_v_ij = neighbor{1}.data.rho_v_ji;
                        else
                            ag.data.rho_x_ij = neighbor{1}.data.rho_x_ji;                            
                        end
                    end
                end
            end
        end

        %% Shift and initialize
        function shift(obj, dt_sample)
            for agent = obj.agents
                % Find the step
                %k = ceil(dt_sample / agent{1}.data.dt) + (ceil(dt_sample / agent{1}.data.dt) == (dt_sample / agent{1}.data.dt));
                k = ceil(dt_sample / agent{1}.data.dt);

                %s = solutions_map(agent.id);
                d = agent{1}.data;
                dynamics = agent{1}.f_i(value(d.x(:,1)), value(d.u(:,1)));
                for neighbour = agent{1}.sending_neighbors
                    %ns = solutions_map(sending_neighbors{j}.id);
                    nd = obj.agent_map(neighbour{1}.id).data;
                    dynamics = dynamics + neighbour{1}.f_ij(value(d.x(:,1)), value(d.u(:,1)), value(nd.x(:,1)), value(nd.u(:,1)));
                end

                x = zeros(size(d.x(:,1)));
                x = value(d.x(:,1)) + dt_sample * dynamics;
                %agent{1}.data.initialize(k);

                % control the state x bounds
                for i=1:size(x,1)
                    if (x(i) < d.x_min)
                        if  false %(x(i) < d.x_min - obj.boundary_tol)
                            error('ADMM_Solver:BoundsViolation', ...
                                'x is outside [x_min, x_max] for this agent.');
                        else
                            warning('ADMM_Solver:BoundsNearViolation', ...
                                'x is outside bounds, projecting back to [x_min, x_max].');
                            x(i) = d.x_min;
                        end
                    elseif  (d.x_max < x(i))
                        if  false %(d.x_max + obj.boundary_tol < x(i))
                            error('ADMM_Solver:BoundsViolation', ...
                                'x is outside [x_min, x_max] for this agent.');
                        else
                            warning('ADMM_Solver:BoundsNearViolation', ...
                                'x is outside bounds, projecting back to [x_min, x_max].');
                            x(i) = d.x_max;
                        end
                    end
                end

                % Update solutions
                s = obj.solution_map(agent{1}.id);
                s.update(x);
                
                % Initialize agent
                d.initialize(x, k)
                
                % Shift 
                % agent{1}.data.shift(k);

                % Initialize and shift neighbors
                for neighbor = agent{1}.neighbors
                    neighbor{1}.data.initialize(k);
                    % neighbor{1}.data.shift(k);
                end
            end
        end

        %% Update solutions
        function update_solutions(obj, x)
            for solution = obj.solutions
                solution{1}.update(x);
            end 
        end

    end
end
