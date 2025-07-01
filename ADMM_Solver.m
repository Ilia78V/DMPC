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

        % Residuals for convergence check
        primal_residual;
        dual_residual;

        approximation;
        % % ADMM states
        % x_opt;
        % u_opt;
        % x_neighbors_opt;
        % u_neighbors_opt;
    end

    methods
        %% Constructor
        function obj = ADMM_Solver(optimizer, agents, solutions, max_iterations, convergence_tolerance, approximation)
            obj.optimizer = optimizer;
            
            obj.agents = agents;
            obj.solutions = solutions;
            obj.max_iterations = max_iterations;
            obj.convergence_tolerance = convergence_tolerance;

            obj.agent_map = containers.Map('KeyType', 'double', 'ValueType', 'any');
            for i = 1:length(obj.agents)
                obj.agent_map(obj.agents{i}.id) = obj.agents{i};
            end

            obj.solution_map = containers.Map('KeyType', 'double', 'ValueType', 'any');
            for i = 1:length(obj.solutions)
                obj.solution_map(obj.solutions{i}.id) = obj.solutions{i};
            end

            obj.approximation = approximation;

            % obj.approximation = containers.Map('KeyType', 'char', 'ValueType', 'logical');
            % obj.approximation('cost') = false;
            % obj.approximation('dynamics') = false;
            % obj.approximation('constraints') = false;

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

                obj.update_residual();
                % Step 7: Check convergence
                if obj.check_convergence()
                    obj.update_previous_data();
                    fprintf('ADMM solver converged at %.2f\n', q)
                    break;
                elseif q == obj.max_iterations
                    fprintf('ADMM solver terminated\n')
                end
                obj.update_previous_data();
            end
        end

        %% Compute local variables (Step 1)
        function compute_local_variables(obj)
            for i = 1:length(obj.agents)
                agent = obj.agents{i};
                d = agent.data;

                % -----Define constraints-----
                constraints = [];
                constraints = [constraints, d.x(:,1) == d.x0]; % Initial state condition
                if ~isempty(d.u0)
                    constraints = [constraints, d.u(:,1) == d.u0]; % Initial control condition
                end
                
                for k = 1:d.N-1
                    % Dynamics constraint (Euler discretization)
                    dynamics = agent.f_i(d.x(:,k), d.u(:,k));
                    for j = 1:length(agent.sending_neighbors)
                        nd = agent.sending_neighbors{j}.data;
                        dynamics = dynamics + agent.sending_neighbors{j}.f_ij(d.x(:,k), d.u(:,k), nd.x_ji(:,k), nd.u_ji(:,k));
                    end

                    constraints = [constraints, d.x(:,k+1) == d.x(:,k) + d.dt * dynamics];
                    
                    % Approximation of neighbor dynamics
                    if obj.approximation('dynamics')
                        for neighbor = agent.sending_neighbors
                            nd = neighbor{1}.data;

                            constraints = [constraints, nd.x_ji(:,1) == obj.agent_map(neighbor{1}.id).data.x0]; % Initial state condition

                            neighbor_dynamics = obj.agent_map(neighbor{1}.id).f_i(nd.x_ji(:,k), nd.u_ji(:,k));
                            neighbor_dynamics = neighbor_dynamics + obj.agent_map(neighbor{1}.id).neighbor_map(agent.id).f_ij(nd.x_ji(:,k), nd.u_ji(:,k), d.x(:,k), d.u(:,k));
                            neighbor_dynamics = neighbor_dynamics + nd.v_ji(:,k);
                            constraints = [constraints, nd.x_ji(:,k+1) == nd.x_ji(:,k) + d.dt * neighbor_dynamics];
                        end
                    end
                    % comput the external influence
                    if obj.approximation('dynamics')
                        for neighbor = agent.receiving_neighbors
                            nd = neighbor{1}.data;
                            v = sdpvar(size(nd.v_i,1), 1);
                            for se_neighbor = agent.sending_neighbors
                                if neighbor{1}.id ~= se_neighbor{1}.id
                                    v = v + neighbor{1}.f_ij(nd.x_ji(:,k), nd.u_ji(:,k), d.x(:,k), d.u(:,k));
                                end
                            end
                        end
                        constraints = [constraints, nd.v_i(:,k) == v];
                    end
                    

                    % State and control bounds
                    constraints = [constraints, d.x_min <= d.x(:,k) <= d.x_max];
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

                % Final state bounds
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
                % if obj.approximation('cost')
                %     cost = cost * (1/(1+length(agent.neighbors)));
                %     for neighbor = agent.neighbors
                %         cost = cost + (1/(1+length(obj.agent_map(neighbor{1}.id).neighbors))) * obj.cost_constructor(neighbor{1});
                %     end
                % end
                %%%%%%%%%%%%%%%% Approximation of neighbor cost %%%%%%%%%%%%%%%%
                if ~isempty(agent.cost_approx_neighbors)
                    cost = cost * (1/(1+length(agent.cost_approx_neighbors)));
                    for neighbor = agent.cost_approx_neighbors
                        cost = cost + (1/(1+length(obj.agent_map(neighbor{1}.id).cost_approx_neighbors))) * obj.cost_constructor(neighbor{1});
                    end
                end

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


                if obj.approximation('dynamics')
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
            elseif isa(agent, 'Neighbor')
               isagent = false;
               d = agent.data;
               x = d.x_ji;
               u = d.u_ji;
               agent = obj.agent_map(agent.id);
            end
            
            if obj.approximation('dynamics')
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
                        cost = cost + (d.dt/2) * (d.z_u(:,k) - d.u(:,k))' * diag(d.rho_u_i) * (d.z_u(:,k) - d.u(:,k));
                        for neighbor = agent.receiving_neighbors
                            nd = neighbor{1}.data; 
                            cost = cost + d.dt * nd.mu_v_i(:,k)' * (nd.z_v_i(:,k) - nd.v_i(:,k));
                            cost = cost + (d.dt/2) * (nd.z_v_i(:,k) - nd.v_i(:,k))' * diag(nd.rho_v_i) * (nd.z_v_i(:,k) - nd.v_i(:,k));
                        end
                        
                        % Neighbor coupling terms (scaled by d.dt for integral form)
                        for neighbor = agent.sending_neighbors
                            nd = neighbor{1}.data;
                            z_neighbor = [nd.u_ji(:,k); nd.v_ji(:,k)];
                            z_coupling_neighbor = [nd.z_u_j(:,k); nd.z_v_j(:,k)];
                            mu_neighbor = [nd.mu_u_ji(:,k); nd.mu_v_ji(:,k)];
                            cost = cost + d.dt * mu_neighbor' * (z_coupling_neighbor - z_neighbor);
                            cost = cost + (d.dt/2) * (z_coupling_neighbor - z_neighbor)' * diag([nd.rho_u_ji; nd.rho_v_ji]) * (z_coupling_neighbor - z_neighbor);
                        end
                    end
    
                end
    
                if isagent
                    % Final state costs
                    % Local coupling terms and Lagrange multipliers (scaled by d.dt for integral form)
                    for neighbor = agent.receiving_neighbors
                        nd = neighbor{1}.data; 
                        cost = cost + d.dt * nd.mu_v_i(:,d.N)' * (nd.z_v_i(:,d.N) - nd.v_i(:,d.N));
                        cost = cost + (d.dt/2) * (nd.z_v_i(:,d.N) - nd.v_i(:,d.N))' * diag(nd.rho_v_i) * (nd.z_v_i(:,d.N) - nd.v_i(:,d.N));
                    end

                    % Neighbor coupling terms (scaled by d.dt for integral form)
                    for neighbor = agent.sending_neighbors
                        cost = cost + d.dt * nd.mu_v_ji(:,d.N)' * (nd.z_v_j(:,d.N) - nd.v_ji(:,d.N));
                        cost = cost + (d.dt/2) * (nd.z_v_j(:,d.N) - nd.v_ji(:,d.N))' * diag(nd.rho_v_ji) * (nd.z_v_j(:,d.N) - nd.v_ji(:,d.N));
                    end
                end
            
            else
                % Default case (without neighborhood approximation)
            
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
                        cost = cost + (d.dt/2) * (z_coupling - z_local)' * d.C_i * (z_coupling - z_local);
                        
                        % Neighbor coupling terms (scaled by d.dt for integral form)
                        for j = 1:length(agent.sending_neighbors)
                            nd = agent.sending_neighbors{j}.data;
                            z_neighbor = [nd.x_ji(:,k); nd.u_ji(:,k)];
                            z_coupling_neighbor = [nd.z_x_j(:,k); nd.z_u_j(:,k)];
                            mu_neighbor = [nd.mu_x_ji(:,k); nd.mu_u_ji(:,k)];
                            cost = cost + d.dt * mu_neighbor' * (z_coupling_neighbor - z_neighbor);
                            cost = cost + (d.dt/2) * (z_coupling_neighbor - z_neighbor)' * nd.C_ji * (z_coupling_neighbor - z_neighbor);
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
                    cost = cost + (d.dt/2) * (z_coupling - z_local)' * diag(d.rho_x_i) * (z_coupling - z_local);
                    
                    % Neighbor coupling terms (scaled by d.dt for integral form)
                    for j = 1:length(agent.sending_neighbors)
                        nd = agent.sending_neighbors{j}.data;
                        z_neighbor = [nd.x_ji(:,d.N)];
                        z_coupling_neighbor = [nd.z_x_j(:,d.N)];
                        mu_neighbor = [nd.mu_x_ji(:,d.N)];
                        cost = cost + d.dt * mu_neighbor' * (z_coupling_neighbor - z_neighbor);
                        cost = cost + (d.dt/2) * (z_coupling_neighbor - z_neighbor)' * diag(nd.rho_x_ji) * (z_coupling_neighbor - z_neighbor);
                    end
                end
            end

        end

        %% Send local copies to sending neighbors (Step 2)
        function send_local_copies(obj)
            for agent = obj.agents
                for neighbor = agent{1}.sending_neighbors
                    ag = obj.agent_map(neighbor{1}.id).neighbor_map(agent{1}.id);
                    ag.data.u_ij = value(neighbor{1}.data.u_ji);
                    
                    if obj.approximation('dynamics')
                        ag.data.v_ij = value(neighbor{1}.data.v_ji);
                    else
                        ag.data.x_ij = value(neighbor{1}.data.x_ji);
                    end
                end
            end
        end

        %% Compute coupling variables (Step 3)
        function compute_coupling_variables(obj)
            for agent = obj.agents
                d = agent{1}.data;
                %pd = agent{1}.previous_data;

                if obj.approximation('dynamics')
                    % Compute local term: (x_i, u_i) - C_i^{-1} * mu_i
                    % z_local_x = value(d.x) - diag(1 ./ d.rho_x_i) * d.mu_x;
                    z_local_u = value(d.u) - diag(1 ./ d.rho_u_i) * d.mu_u;

                    z_local_v = cell(size(agent{1}.receiving_neighbors));
                    for i = 1:length(agent{1}.receiving_neighbors)
                        nd = agent{1}.receiving_neighbors{i}.data;
                        z_local_v{i} = value(nd.v_i) - diag(1 ./ nd.rho_v_i) * nd.mu_v_i;
                    end
    
                    % Compute neighbor contribution
                    z_neighbor_u = zeros(size(z_local_u));
                    z_neighbor_v = cellfun(@(m) zeros(size(z_local_v{1})), z_local_v, 'UniformOutput', false);
                    cell(size(z_local_v));
                    for i = 1:length(agent{1}.receiving_neighbors)
                        nd = agent{1}.receiving_neighbors{i}.data;
                        z_neighbor_u = z_neighbor_u + nd.u_ij - diag(1 ./ nd.rho_u_ij) * nd.mu_u_ij;
                        delta_v = nd.v_ij - diag(1 ./ nd.rho_v_ij) * nd.mu_v_ij;
                        z_neighbor_v = cellfun(@(z) z + delta_v, z_neighbor_v, 'UniformOutput', false);
                    end
                    
                    % Compute final coupling variable update
                    d.z_u = (1 / (1 + length(agent{1}.receiving_neighbors))) * (z_local_u + z_neighbor_u);
                    for i = 1:length(agent{1}.receiving_neighbors)
                        nd = agent{1}.receiving_neighbors{i}.data;
                        nd.z_v_i = (1 / (1 + length(agent{1}.receiving_neighbors))) * (z_local_v{i} + z_neighbor_v{i});
                    end
                
                else
                    % Compute local term: (x_i, u_i) - C_i^{-1} * mu_i
                    z_local_x = value(d.x) - diag(1 ./ d.rho_x_i) * d.mu_x;
                    z_local_u = value(d.u) - diag(1 ./ d.rho_u_i) * d.mu_u;
    
                    % Compute neighbor contribution
                    z_neighbor_x = zeros(size(z_local_x));
                    z_neighbor_u = zeros(size(z_local_u));
                    for neighbor = agent{1}.receiving_neighbors
                        nd = neighbor{1}.data;
                        z_neighbor_x = z_neighbor_x + nd.x_ij - diag(1 ./ nd.rho_x_ij) * nd.mu_x_ij;
                        z_neighbor_u = z_neighbor_u + nd.u_ij - diag(1 ./ nd.rho_u_ij) * nd.mu_u_ij;
                    end
                    
                    % Compute final coupling variable update
                    d.z_x = (1 / (1 + length(agent{1}.receiving_neighbors))) * (z_local_x + z_neighbor_x);
                    d.z_u = (1 / (1 + length(agent{1}.receiving_neighbors))) * (z_local_u + z_neighbor_u);
                end

            end
        end


        %% Send coupling variables to receiving neighbors (Step 4)
        function send_coupling_variables(obj)
            for agent = obj.agents
                for neighbor = agent{1}.receiving_neighbors
                    ag = obj.agent_map(neighbor{1}.id).neighbor_map(agent{1}.id);
                    ag.data.z_u_j = agent{1}.data.z_u;
                    
                    if obj.approximation('dynamics')
                        ag.data.z_v_j = ag.data.z_v_i;
                    else
                        ag.data.z_x_j = agent{1}.data.z_x;
                    end
                end
            end
        end

        %% Compute Lagrange multipliers (Step 5)
        function compute_multipliers(obj)
            for agent = obj.agents
                d = agent{1}.data;

                if obj.approximation('dynamics')
                    d.mu_u = d.mu_u + diag(d.rho_u_i) * (d.z_u - value(d.u));

                    for i = 1:length(agent{1}.receiving_neighbors)
                        nd = agent{1}.receiving_neighbors{i}.data;
                        nd.mu_v_i = nd.mu_v_i + diag(nd.rho_v_i) * (nd.z_v_i - value(nd.v_i));
                    end
                    
                    % Neighbor contribution
                    for neighbor = agent{1}.sending_neighbors
                        nd = neighbor{1}.data;
                        nd.mu_u_ji = nd.mu_u_ji + diag(nd.rho_u_ji) * (nd.z_u_j - value(nd.u_ji));
                        nd.mu_v_ji = nd.mu_v_ji + diag(nd.rho_v_ji) * (nd.z_v_j - value(nd.v_ji));
                    end
                
                else
                    d.mu_x = d.mu_x + diag(d.rho_x_i) * (d.z_x - value(d.x));
                    d.mu_u = d.mu_u + diag(d.rho_u_i) * (d.z_u - value(d.u));                
                    
                    % Neighbor contribution
                    for neighbor = agent{1}.sending_neighbors
                        nd = neighbor{1}.data;
                        nd.mu_x_ji = nd.mu_x_ji + diag(nd.rho_x_ji) * (nd.z_x_j - value(nd.x_ji));
                        nd.mu_u_ji = nd.mu_u_ji + diag(nd.rho_u_ji) * (nd.z_u_j - value(nd.u_ji));
                    end
                end
            end
        end
        %% Send Lagrange multipliers to sending neighbors (Step 6)
        function send_multipliers(obj)
            for agent = obj.agents
                for neighbor = agent{1}.sending_neighbors
                    ag = obj.agent_map(neighbor{1}.id).neighbor_map(agent{1}.id);
                    ag.data.mu_u_ij = neighbor{1}.data.mu_u_ji;
                    
                    if obj.approximation('dynamics')
                        ag.data.mu_v_ij = neighbor{1}.data.mu_v_ji;
                    else
                        ag.data.mu_x_ij = neighbor{1}.data.mu_x_ji;
                    end
                end
            end
        end

        %% Check convergence (Step 7)
        function is_converged = check_convergence(obj) 
            is_converged = true;
            for agent = obj.agents
                d = agent{1}.data;
                pd = agent{1}.previous_data;

                if obj.approximation('dynamics')
                    error_vec = vecnorm([(d.z_u - pd.z_u); (d.mu_u - pd.mu_u)], 2, 2); % Computes the 2-norm (Euclidean norm) of each row
                    for neighbor = agent{1}.sending_neighbors
                        nd = neighbor{1}.data;
                        npd = neighbor{1}.previous_data;
                        [error_vec; vecnorm([(nd.z_v_i - npd.z_v_i); (nd.mu_v_i - npd.mu_v_i)], 2, 2)];
                    end
                    error = norm(error_vec);
                else
                    error = norm([(d.z_x(:,1:end-1) - pd.z_x(:,1:end-1)); ...
                                  (d.z_u - pd.z_u); ...
                                  (d.mu_x(:,1:end-1) - pd.mu_x(:,1:end-1)); ...
                                  (d.mu_u - pd.mu_u)]);
                end 
                
                is_converged = is_converged && (error <= obj.convergence_tolerance);
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

        %% Calculate primal and dual residual;
        function update_residual(obj)
            for agent = obj.agents
                d = agent{1}.data;
                pd = agent{1}.previous_data;  
                
                if obj.approximation('dynamics')
                    %primal residual
                    local_pr_u = norm(vecnorm(value(d.u) - d.z_u, 2, 1), 1)/(d.N - 1);
                    local_pr_v = zeros(length(agent{1}.receiving_neighbors),1);
                    for i = 1: length(agent{1}.receiving_neighbors)
                        nd = agent{1}.receiving_neighbors{i}.data;
                        local_pr_v(i) = norm(vecnorm(value(nd.v_i) - nd.z_v_i, 2, 1), 1)/d.N;
                    end

                    
                    %dual residual
                    local_dr_u = norm(vecnorm((diag(pd.rho_u_i) * (d.z_u - pd.z_u)), 2, 1), 1)/(d.N - 1);
                    local_dr_v = zeros(length(agent{1}.receiving_neighbors),1);
                    for i = 1: length(agent{1}.receiving_neighbors)
                        nd = agent{1}.receiving_neighbors{i}.data;
                        npd = agent{1}.receiving_neighbors{i}.previous_data;
                        local_dr_v(i) = norm(vecnorm((diag(npd.rho_v_i) * (nd.z_v_i - npd.z_v_i)), 2, 1), 1)/d.N;
                    end
                    
                    % initial = true;
                    % **Ensure neighbor variables are initialized** before the loop
                    neighbor_pr_v = zeros(length(agent{1}.sending_neighbors),1);
                    neighbor_pr_u = zeros(length(agent{1}.sending_neighbors),1);
                    neighbor_dr_v = zeros(length(agent{1}.sending_neighbors),1);
                    neighbor_dr_u = zeros(length(agent{1}.sending_neighbors),1);
    
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
    
                        %primal residual
                        neighbor_pr_u(i) = norm(vecnorm( (value(nd.u_ji) - nd.z_u_j) , 2, 1), 1)/(d.N - 1);
                        neighbor_pr_v(i) = norm(vecnorm( (value(nd.v_ji) - nd.z_v_j) , 2, 1), 1)/d.N;
                        
                        %dual residual
                        neighbor_dr_u(i) = norm(vecnorm( (diag(npd.rho_u_ji) * (nd.z_u_j - npd.z_u_j)) , 2, 1), 1)/(d.N - 1);
                        neighbor_dr_v(i) = norm(vecnorm( (diag(npd.rho_v_ji) * (nd.z_v_j - npd.z_v_j)) , 2, 1), 1)/d.N;
                    end
    
                    % neighbor_pr_x = norm(vecnorm(neighbor_pr_x, 2, 1), 1)/d.N;
                    % neighbor_pr_u = norm(vecnorm(neighbor_pr_u, 2, 1), 1)/(d.N - 1);
                    % neighbor_dr_x = norm(vecnorm(neighbor_dr_x, 2, 1), 1)/d.N;
                    % neighbor_dr_u = norm(vecnorm(neighbor_dr_u, 2, 1), 1)/(d.N - 1);
    
                    d.primal_residual = [d.primal_residual, norm([local_pr_v; local_pr_u; neighbor_pr_v; neighbor_pr_u], 2)];
                    d.dual_residual = [d.dual_residual, norm([local_dr_v; local_dr_u; neighbor_dr_v; neighbor_dr_u], 2)];
                
                else
                    %DEFAULT CASE

                    %primal residual
                    local_pr_x = norm(vecnorm(value(d.x) - d.z_x, 2, 1), 1)/d.N;
                    local_pr_u = norm(vecnorm(value(d.u) - d.z_u, 2, 1), 1)/(d.N - 1);
                    
                    %dual residual
                    local_dr_x = norm(vecnorm((diag(pd.rho_x_i) * (d.z_x - pd.z_x)), 2, 1), 1)/d.N;
                    local_dr_u = norm(vecnorm((diag(pd.rho_u_i) * (d.z_u - pd.z_u)), 2, 1), 1)/(d.N - 1);
                    
                    % initial = true;
                    % % **Ensure neighbor variables are initialized** before the loop
                    % neighbor_pr_x = 0;
                    % neighbor_pr_u = 0;
                    % neighbor_dr_x = 0;
                    % neighbor_dr_u = 0;

                    neighbor_pr_x = zeros(length(agent{1}.sending_neighbors),1);
                    neighbor_pr_u = zeros(length(agent{1}.sending_neighbors),1);
                    neighbor_dr_x = zeros(length(agent{1}.sending_neighbors),1);
                    neighbor_dr_u = zeros(length(agent{1}.sending_neighbors),1);
    
                    % Neighbor contribution
                    tmp = 0;
                    for neighbor = agent{1}.sending_neighbors
                        tmp = tmp + 1;
                        nd = neighbor{1}.data;
                        npd = neighbor{1}.previous_data;
                        
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

                       %primal residual
                        neighbor_pr_x(tmp) = norm(vecnorm( (value(nd.x_ji) - nd.z_x_j), 2, 1), 1)/(d.N - 1);
                        neighbor_pr_u(tmp) = norm(vecnorm( (value(nd.u_ji) - nd.z_u_j) , 2, 1), 1)/d.N;
    
                       %dual residual
                        neighbor_dr_x(tmp) = norm(vecnorm( diag(npd.rho_x_ji) * (nd.z_x_j - npd.z_x_j) , 2, 1), 1)/(d.N - 1);
                        neighbor_dr_u(tmp) = norm(vecnorm( diag(npd.rho_u_ji) * (nd.z_u_j - npd.z_u_j) , 2, 1), 1)/d.N ;

                        % neighbor_pr_u(i) = norm(vecnorm( (value(nd.u_ji) - nd.z_u_j) , 2, 1), 1)/(d.N - 1);
                        % neighbor_pr_v(i) = norm(vecnorm( (value(nd.v_ji) - nd.z_v_j) , 2, 1), 1)/d.N;
                    end
                    % 
                    % neighbor_pr_x = norm(vecnorm(neighbor_pr_x, 2, 1), 1)/d.N;
                    % neighbor_pr_u = norm(vecnorm(neighbor_pr_u, 2, 1), 1)/(d.N - 1);
                    % neighbor_dr_x = norm(vecnorm(neighbor_dr_x, 2, 1), 1)/d.N;
                    % neighbor_dr_u = norm(vecnorm(neighbor_dr_u, 2, 1), 1)/(d.N - 1);
    
                    d.primal_residual = [d.primal_residual, norm([local_pr_x; local_pr_u; neighbor_pr_x; neighbor_pr_u], 2)];
                    d.dual_residual = [d.dual_residual, norm([local_dr_x; local_dr_u; neighbor_dr_x; neighbor_dr_u], 2)];
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
