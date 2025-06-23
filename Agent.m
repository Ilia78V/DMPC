classdef Agent < handle 
    properties (Access = public)
        % Agent Parameters
        id;

        %Model
        f_i;
        V_i;
        l_i;
        g_i;
        g_i_N;
        h_i;
        h_i_N;

        % Vector for neighbors
        neighbors;
        neighbor_map;
        receiving_neighbors;
        sending_neighbors;
        approx_neighbors;

        % Agent data
        data;
        previous_data;
        
        % Solver
        solver;

        % 
        % agentState;
        % previous_agentState;
        % 
        % couplingState;
        % previous_couplingState;
        % 
        % multiplierState;
        % previous_multiplierState;
        % 
        % penaltyState;
        % initial_penalty;
        % 
        % % Optimization parameters
        % optimizationInfo;
        % solution;
%% to do
        % % Neighbor approximation
        % is_approximating;
        % is_approximatingCost;
        % is_approximatingConstraints;
        % is_approximatingDynamics;
        % 
        % % Neighbor approximation switch parameters
        % delay_c;
        % thr;
        % needs_approx_switch_state;
        % if_init;
        % relative_ADMM_PrimalResiduals;
        % relative_ADMM_DualResiduals;
        % envelope;
        % residual_data;
        % 
        % decay_factor;
        % decay_floor;
    end

    methods
        %% Constructor
        function obj = Agent(id, f_i, V_i, l_i, g_i, g_i_N, h_i, h_i_N, agent_data)
            % Initialize key properties
            obj.id = id;

            obj.f_i = f_i;
            obj.V_i = V_i;
            obj.l_i = l_i;
            obj.g_i = g_i;
            obj.g_i_N = g_i_N;
            obj.h_i = h_i;
            obj.h_i_N = h_i_N;

            obj.data = agent_data;
            obj.previous_data = copy(agent_data);

            obj.data.agent = obj;
            obj.previous_data.agent = obj;            
        end

        %% Register neighbors
        % function register_neighbors(obj, neighbors)
        %     obj.neighbors = neighbors;
        %     obj.neighbor_map = containers.Map('KeyType', 'double', 'ValueType', 'any');
        %     obj.receiving_neighbors = {};
        %     obj.sending_neighbors = {};
        %     for i = 1:length(obj.neighbors)
        %         obj.neighbor_map(obj.neighbors{i}.id) = obj.neighbors{i};
        %         if obj.neighbors{i}.sending == true
        %             obj.sending_neighbors{end+1} = obj.neighbors{i};
        %         end
        %         if obj.neighbors{i}.receiving == true
        %             obj.receiving_neighbors{end+1} = obj.neighbors{i};
        %         end
        %     end
        % end

        function register_neighbors(obj, new_neighbors)
            % Initialize once if not already done
            if isempty(obj.neighbors)
                obj.neighbors = {};
                obj.neighbor_map = containers.Map('KeyType', 'double', 'ValueType', 'any');
                obj.receiving_neighbors = {};
                obj.sending_neighbors = {};
                obj.approx_neighbors = {};
            end
        
            for i = 1:length(new_neighbors)
                new_id = new_neighbors{i}.id;
        
                % If this neighbor is already registered, skip or update
                if isKey(obj.neighbor_map, new_id)
                    warning('Neighbor with ID %d already exists. Skipping.', new_id);
                    continue;
                end
        
                % Register new neighbor
                obj.neighbors{end+1} = new_neighbors{i};
                obj.neighbor_map(new_id) = new_neighbors{i};
        
                if new_neighbors{i}.sending
                    obj.sending_neighbors{end+1} = new_neighbors{i};
                end
                if new_neighbors{i}.receiving
                    obj.receiving_neighbors{end+1} = new_neighbors{i};
                end
                if new_neighbors{i}.approx
                    obj.approx_neighbors{end+1} = new_neighbors{i};
                end
            end
        end
        %% Update agent state
        function update_agentState(obj, x_opt, u_opt)
            obj.previous_agentState = obj.agentState;
            obj.agentState.x = x_opt;
            obj.agentState.u = u_opt;
        end

        % %% Get approximation data
        % function get_approx_data(obj, approx)
        %     obj.data.approximation = approx;
        % end
    
        %%
        function register_solver(obj, solver)
            obj.solver = solver;
        end
    end
end
