classdef Solution < handle
    properties
        id;
        dt_sample;
        t;
        x;
        u;
        cost;
        agent;
    end
    
    methods
        %% Constructor
        function obj = Solution(agent, dt_sample)
            obj.id = agent.id;
            obj.agent = agent;
            obj.t = [0];
            obj.x = agent.data.x0;
            obj.u = [];
            obj.cost = [];
            obj.dt_sample = dt_sample;
        end

        %% Set solution
        function update(obj, x)
            % Update time
            obj.t = [obj.t, obj.t(end) + obj.dt_sample];
            % Find the step
            %k = ceil(obj.dt_sample / obj.agent.data.dt);
            % Update state
            %new_x = value(obj.agent.data.x(:,k));
            obj.x = [obj.x, x];
            % Update control            
            new_u = value(obj.agent.data.u(:,1));
            obj.u = [obj.u, new_u];
            % Update cost            
            new_cost = value(obj.agent.data.cost(end));
            obj.cost = [obj.cost, new_cost];
        end
    end
end