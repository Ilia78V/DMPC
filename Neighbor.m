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
    end

    methods
        %% Constructor
        function obj = Neighbor(id, agent, sending, receiving, f_ij, g_ij, g_ij_N, h_ij, h_ij_N, V_ij, l_ij, neighbor_data)
            obj.id = id;
            obj.agent = agent;
            obj.sending = sending;
            obj.receiving = receiving;
            obj.approx = false;
            if neighbor_data.approximation('cost') || neighbor_data.approximation('dynamics') || neighbor_data.approximation('constraints') 
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
