classdef Agent_data < handle & matlab.mixin.Copyable
    properties
        agent;

        % Time Variable
        t0;             % Initial time step
        t;              % Time vector
        T;              % Time horizon
        dt;             % Time step size
        N;              % Number of time steps

        % Agent parameters
        id;             % Agent ID
        n_x;            % Dimension of state vector for the agent
        n_u;            % Dimension of control vector for the agent
               
        % Define state and control bounds for the agent
        x_min;          % State bounds
        x_max;
        u_min;          % Control bounds
        u_max;

        % Initial and final conditions
        x0;
        x_ref;
        u0;
        u_ref;

        % Define rho values for local penalties
        rho_x_i;        % Local state penalty. Example:  = 10 * ones(n_x, 1)
        rho_u_i;        % Local control penalty
        % rho_v_i;

        % Coupling parameters based on rho values
        C_i;            % Local coupling penalties
        
        % Agent variables
        x;              % State trajectory
        u;              % Control trajectory
        % v;
        z_x;            % Local coupling variable for state
        z_u;            % Local coupling variable for control
        % z_v;
        mu_x;           % Lagrange multiplier for state
        mu_u;           % Lagrange multiplier for control 
        % mu_v;
       
        %ADMM variables
        primal_residual;
        dual_residual;
        cost;
        approximation;
        
        %% to be added
        % % Vector containing external influence
        % v;
    end
    
    methods
        %% Constructor
        function obj = Agent_data(id, n_x, n_u, t0, T, N, x0, x_ref, x_min, x_max, u_min, u_max, rho_init, approximation)
            if nargin > 0
                obj.id = id;

                obj.t0 = t0;
                obj.T = T;
                obj.t = linspace(t0, T, N);
                obj.dt = T / (N-1);
                obj.N = N;
                
                obj.n_x = n_x;            
                obj.n_u = n_u;            
                
                obj.x_min = x_min;        
                obj.x_max = x_max;
                obj.u_min = u_min;         
                obj.u_max = u_max;
        
                obj.x0 = x0;
                obj.x_ref = x_ref;
        
                obj.rho_x_i = rho_init * ones(obj.n_x, 1);       
                obj.rho_u_i = rho_init * ones(obj.n_u, 1);
                % obj.rho_v_i = rho_init * ones(obj.n_x, 1);

                obj.approximation = approximation;
                
                if approximation('dynamics') == false
                    obj.C_i = diag([obj.rho_x_i; obj.rho_u_i]);
                    
                    obj.x = sdpvar(obj.n_x, N);             
                    obj.u = sdpvar(obj.n_u, N-1);            
                    obj.z_x = zeros(obj.n_x, N);          
                    obj.z_u = zeros(obj.n_u, N-1);          
                    obj.mu_x = zeros(obj.n_x, N);           
                    obj.mu_u = zeros(obj.n_u, N-1);
                else
                    obj.C_i = diag([obj.rho_u_i]);

                    obj.x = sdpvar(obj.n_x, N);             
                    obj.u = sdpvar(obj.n_u, N-1);
                    % obj.v = repmat({zeros(obj.n_x, N)}, 1, length(agent.neighbors));          
                    obj.z_u = zeros(obj.n_u, N-1);
                    % obj.z_v = zeros(obj.n_x, N);           
                    obj.mu_u = zeros(obj.n_u, N-1);
                    % obj.mu_v = zeros(obj.n_x, N);
                end

                obj.primal_residual = [];
                obj.dual_residual= [];
                obj.cost = [];
            end
        end

        %% Initialize
        function initialize(obj, x, k)
            obj.x0 = x; 
            
            x0 = value(obj.x(:, k+1:end));
            u0 = value(obj.u(:, k+1:end));
            obj.z_u = [obj.z_u(:, k+1:end)];
            obj.mu_u = [obj.mu_u(:, k+1:end)];

            if obj.approximation('dynamics')
                for neighbor = obj.agent.receiving_neighbors
                    nd = neighbor{1}.data;
                    nd.z_v_i = [nd.z_v_i(:, k+1:end)]; % zeros(size(obj.z_x,1), k)];
                    nd.mu_v_i = [nd.mu_v_i(:, k+1:end)];
                end
            else
                obj.z_x = [obj.z_x(:, k+1:end)]; % zeros(size(obj.z_x,1), k)];
                obj.mu_x = [obj.mu_x(:, k+1:end)];
            end
            
            
            for i=1:k
                x0 = [x0, x0(:, end)];
                u0 = [u0, u0(:, end)];
                obj.z_u = [obj.z_u, obj.z_u(:, end)];   
                obj.mu_u = [obj.mu_u, obj.mu_u(:, end)];
                
                if obj.approximation('dynamics')
                    for neighbor = obj.agent.receiving_neighbors
                        nd = neighbor{1}.data;
                        nd.z_v_i = [nd.z_v_i, nd.z_v_i(:, end)]; % zeros(size(obj.z_x,1), k)];
                        nd.mu_v_i = [nd.mu_v_i, nd.mu_v_i(:, end)];
                    end
                else
                    obj.z_x = [obj.z_x, obj.z_x(:, end)]; % zeros(size(obj.z_x,1), k)];
                    obj.mu_x = [obj.mu_x, obj.mu_x(:, end)];
                end 
            end

            assign(obj.x, x0);
            assign(obj.u, u0);
            
            %value(obj.x(:, k));
            %obj.u0 = value(obj.u(:, k));

            % Reinitialize the decision variable for the new horizon:
            %obj.x = sdpvar(obj.n_x, obj.N);
            % Optionally, reinitialize d.u as well if needed:
            %obj.u = sdpvar(obj.n_u, obj.N-1);
        end

        %% Shift
        % function shift(obj, k)
        % 
        %     obj.z_x = [obj.z_x(:, k+1:end), obj.z_x(:, end)]; % zeros(size(obj.z_x,1), k)];
        %     obj.z_u = [obj.z_u(:, k+1:end), obj.z_u(:, end)];
        %     obj.mu_x = [obj.mu_x(:, k+1:end), obj.mu_x(:, end)];
        %     obj.mu_u = [obj.mu_u(:, k+1:end), obj.mu_u(:, end)];
        % 
        % 
        % end
        % 
    end
end
