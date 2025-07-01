classdef Neighbor_data < handle & matlab.mixin.Copyable
    properties
        % Neighbor parameters
        id;             % Neighbor ID
        agent_data;
        approximation;
        
        % Time Variable
        t0;             % Initial time step
        t;              % Time vector
        T;              % Time horizon
        dt;             % Time step size
        N;              % Number of time steps
               
        % Neighbor parameters       
        n_x;            % Neighbor state dimension
        n_u;            % Neighbor control dimension

        % Define rho values for neighbor penalties
        rho_x_ij;       
        rho_u_ij;
        rho_v_ij;
        rho_x_ji;       % Neighbor state penalties
        rho_u_ji;       % Neighbor control penalties
        rho_v_ji;
        rho_v_i;

        % Coupling parameters based on rho values
        C_ij;
        C_ji;           % Neighbor coupling penalties

        % Neighbor variables
        x_ij;           % recevied Neighbors state local copies
        u_ij;           % Neighbors control local copies
        v_ij;

        x_ji;            % Neighbor state trajectory that supposed to be sent
        u_ji;            % Neighbor control trajectory
        v_ji;

        v_i;
        
        z_x_j;          % Neighbor state coupling variables
        z_u_j;          % Neighbor control coupling variables
        z_v_j;

        z_v_i;
        
        mu_x_ij;        
        mu_u_ij;
        mu_v_ij;
        mu_x_ji;        % Neighbor Lagrange multiplier for state
        mu_u_ji;        % Neighbor Lagrange multiplier for control
        mu_v_ji;
        mu_v_i;
       
        %% to be added
        % % Vector containing external influence
        % v;
    end
    
    methods
        %% Constructor
        function obj = Neighbor_data(id, n_x, n_u, agent, rho_init, approximation)
            if nargin > 0
                obj.id = id;
                obj.agent_data = agent.data;
                obj.approximation = approximation; %approx = containers.Map({'cost','dynamics','constraints'},{flag,flag,flag});

                obj.t0 = obj.agent_data.t0;
                obj.t = linspace(obj.agent_data.t0, obj.agent_data.T, obj.agent_data.N);   
                obj.T = obj.agent_data.T;
                obj.dt = obj.agent_data.T / (obj.agent_data.N-1);
                obj.N = obj.agent_data.N;
                
                obj.n_x = n_x;            
                obj.n_u = n_u;      

                obj.rho_x_ij = rho_init * ones(obj.agent_data.n_x, 1);       
                obj.rho_u_ij = rho_init * ones(obj.agent_data.n_u, 1);
                obj.rho_v_ij = rho_init * ones(obj.agent_data.n_x, 1);
                obj.rho_x_ji = rho_init * ones(n_x, 1);       
                obj.rho_u_ji = rho_init * ones(n_u, 1);
                obj.rho_v_ji = rho_init * ones(n_x, 1);
                obj.rho_v_i = rho_init * ones(obj.agent_data.n_x, 1);

                if obj.agent_data.approximation('dynamics') == false
                    obj.C_ij = diag([obj.rho_x_ij; obj.rho_u_ij]);                  
                    obj.C_ji = diag([obj.rho_x_ji; obj.rho_u_ji]);                  
    
                    obj.x_ij = zeros(obj.agent_data.n_x, obj.N);
                    obj.u_ij = zeros(obj.agent_data.n_u, obj.N-1); 
                    obj.x_ji = sdpvar(obj.n_x, obj.N); 
                    obj.u_ji = sdpvar(obj.n_u, obj.N-1); 
                    
                    obj.z_x_j = zeros(obj.n_x, obj.N); 
                    obj.z_u_j = zeros(obj.n_u, obj.N-1); 
                    
                    obj.mu_x_ij = zeros(obj.agent_data.n_x, obj.N);
                    obj.mu_u_ij = zeros(obj.agent_data.n_u, obj.N-1);
                    obj.mu_x_ji = zeros(obj.n_x, obj.N);
                    obj.mu_u_ji = zeros(obj.n_u, obj.N-1);
                else
                    obj.C_ij = diag([obj.rho_u_ij; obj.rho_v_ij]);                  
                    obj.C_ji = diag([obj.rho_u_ji; obj.rho_v_ji]);                  
    
                    obj.x_ij = zeros(obj.agent_data.n_x, obj.N);
                    obj.u_ij = zeros(obj.agent_data.n_u, obj.N-1);
                    obj.v_ij = zeros(obj.agent_data.n_x, obj.N);

                    obj.x_ji = sdpvar(obj.n_x, obj.N); 
                    obj.u_ji = sdpvar(obj.n_u, obj.N-1); 
                    obj.v_ji = sdpvar(obj.n_x, obj.N);

                    obj.v_i = sdpvar(obj.agent_data.n_x, obj.N);

                    obj.z_x_j = zeros(obj.n_x, obj.N); 
                    obj.z_u_j = zeros(obj.n_u, obj.N-1); 
                    obj.z_v_j = zeros(obj.n_x, obj.N);

                    obj.z_v_i = zeros(obj.agent_data.n_x, obj.N);

                    obj.mu_x_ij = zeros(obj.agent_data.n_x, obj.N);
                    obj.mu_u_ij = zeros(obj.agent_data.n_u, obj.N-1);
                    obj.mu_v_ij = zeros(obj.agent_data.n_x, obj.N);

                    obj.mu_x_ji = zeros(obj.n_x, obj.N);
                    obj.mu_u_ji = zeros(obj.n_u, obj.N-1);
                    obj.mu_v_ji = zeros(obj.n_x, obj.N);
                    
                    obj.mu_v_i = zeros(obj.agent_data.n_x, obj.N);
                end                
            end
        end
        
        %% Initialize
        function initialize(obj, k)
            x_ji0 = value(obj.x_ji(:, k+1:end));
            u_ji0 = value(obj.u_ji(:, k+1:end));
            
            obj.u_ij   = [obj.u_ij(:, k+1:end)];
            obj.z_u_j  = [obj.z_u_j(:, k+1:end)];
            obj.mu_u_ij = [obj.mu_u_ij(:, k+1:end)];
            obj.mu_u_ji = [obj.mu_u_ji(:, k+1:end)];
            
            if obj.agent_data.approximation('dynamics')
                v_ji0 = value(obj.v_ji(:, k+1:end));

                obj.v_ij   = [obj.v_ij(:, k+1:end)];
                obj.z_v_j  = [obj.z_v_j(:, k+1:end)];
                obj.mu_v_ij = [obj.mu_v_ij(:, k+1:end)];
                obj.mu_v_ji = [obj.mu_v_ji(:, k+1:end)];
            else
                obj.x_ij   = [obj.x_ij(:, k+1:end)];
                obj.z_x_j  = [obj.z_x_j(:, k+1:end)];
                obj.mu_x_ij = [obj.mu_x_ij(:, k+1:end)];
                obj.mu_x_ji = [obj.mu_x_ji(:, k+1:end)];
            end

            for i=1:k
                x_ji0 = [x_ji0, x_ji0(:, end)];
                u_ji0 = [u_ji0, u_ji0(:, end)];
                
                obj.u_ij   = [obj.u_ij,   obj.u_ij(:, end)];
                obj.z_u_j  = [obj.z_u_j,   obj.z_u_j(:, end)];
                obj.mu_u_ij = [obj.mu_u_ij, obj.mu_u_ij(:, end)];
                obj.mu_u_ji = [obj.mu_u_ji,  obj.mu_u_ji(:, end)];

                if obj.agent_data.approximation('dynamics')
                    v_ji0 = [v_ji0, v_ji0(:, end)];

                    obj.v_ij   = [obj.v_ij,  obj.v_ij(:, end)];
                    obj.z_v_j  = [obj.z_v_j,   obj.z_v_j(:, end)];
                    obj.mu_v_ij = [obj.mu_v_ij, obj.mu_v_ij(:, end)];
                    obj.mu_v_ji = [obj.mu_v_ji, obj.mu_v_ji(:, end)];
                else
                    obj.x_ij   = [obj.x_ij,  obj.x_ij(:, end)];
                    obj.z_x_j  = [obj.z_x_j,   obj.z_x_j(:, end)];
                    obj.mu_x_ij = [obj.mu_x_ij, obj.mu_x_ij(:, end)];
                    obj.mu_x_ji = [obj.mu_x_ji, obj.mu_x_ji(:, end)];
                end
            end

            assign(obj.x_ji, x_ji0);
            assign(obj.u_ji, u_ji0);
            if obj.agent_data.approximation('dynamics')
                assign(obj.v_ji, v_ji0);
            end

            %obj.x_ji = sdpvar(obj.n_x, obj.N); 
            %obj.u_ji = sdpvar(obj.n_u, obj.N-1);
        end
        
        %% Shift
        % function shift(obj, k)
        %     obj.x_ij   = [obj.x_ij(:, k+1:end),  obj.x_ij(:, end)];
        %     obj.u_ij   = [obj.u_ij(:, k+1:end),   obj.u_ij(:, end)];
        %     obj.z_x_j  = [obj.z_x_j(:, k+1:end),   obj.z_x_j(:, end)];
        %     obj.z_u_j  = [obj.z_u_j(:, k+1:end),   obj.z_u_j(:, end)];
        %     obj.mu_x_ij = [obj.mu_x_ij(:, k+1:end), obj.mu_x_ij(:, end)];
        %     obj.mu_u_ij = [obj.mu_u_ij(:, k+1:end), obj.mu_u_ij(:, end)];
        %     obj.mu_x_ji = [obj.mu_x_ji(:, k+1:end), obj.mu_x_ji(:, end)];
        %     obj.mu_u_ji = [obj.mu_u_ji(:, k+1:end),  obj.mu_u_ji(:, end)];
        % end
    end
end
