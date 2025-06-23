classdef Neighbor < handle 
    properties (Access = public)

        % Neighbor parameters
        id;
        agent;
        sending;
        receiving;
        approx;
        numberOfNeighbors;

        % Coupling model
        f_ij;
        g_ij;
        g_ij_N;
        h_ij;
        h_ij_N;
        V_ij;
        l_ij;

        % Neighbor data
        data;
        previous_data;

        % % States for neighbor
        % local_copies;
        % coupled_multiplierState;
        % coupled_penaltyState;
        % 
        % neighbors_localCopies;
        % neighbors_couplingState;
        % neighbors_coupled_multiplierState;
        % neighbors_coupled_penaltyState;
        % 
        % previous_couplingState;
        % previous_multiplierState;
        % 
        % initial_penalty;
%% to do
        % externalInfluence_agentState;
        % externalInfluence_couplingState;
        % externalInfluence_multiplierState;
        % externalInfluence_penaltyState;
        % 
        % neighbors_externalInfluence_couplingState;
        % neighbors_externalInfluence_multiplierState;
        % neighbors_externalInfluence_penaltyState;
        % 
        % previous_externalInfluence_couplingState;
        % previous_externalInfluence_multiplierState;
        % previous_neighbors_externalInfluence_couplingState;
        % previous_neighbors_couplingState;
        % 
        % % States for cost approximation
        % neighbors_desiredAgentState;
        % 
        % % Neighbor approximation
        % approximate_neighbor;
        % 
        % is_approximating;
        % is_approximatingCost;
        % is_approximatingConstraints;
        % is_approximatingDynamics;
    end

    methods
        %% Constructor
        function obj = Neighbor(id, agent, sending, receiving, f_ij, g_ij, g_ij_N, h_ij, h_ij_N, V_ij, l_ij, neighbor_data)
            obj.id = id;
            obj.agent = agent;
            obj.sending = sending;
            obj.receiving = receiving;
            obj.approx = false;
            if neighbor_data.approx('cost') || neighbor_data.approx('dynamics') || neighbor_data.approx('constraints') 
                obj.approx = true; 
            end

            obj.f_ij = f_ij;
            obj.g_ij = g_ij;
            obj.g_ij_N = g_ij_N;
            obj.h_ij = h_ij;
            obj.h_ij_N = h_ij_N;
            obj.V_ij = V_ij;
            obj.l_ij = l_ij;
            
            obj.data = neighbor_data;
            obj.previous_data = copy(neighbor_data);            
        end

        %% Update local copies for neighbor
        function update_local_copies (obj, x_neighbors_opt, u_neighbors_opt)
            obj.local_copies.x = x_neighbors_opt;
            obj.local_copies.u = u_neighbors_opt;
        end
        %% Update neighbor local copies
        function update_neighbor_local_copies()
            obj.local_copies.x = x_neighbors_opt;
            obj.local_copies.u = u_neighbors_opt;
        end
    end
end
