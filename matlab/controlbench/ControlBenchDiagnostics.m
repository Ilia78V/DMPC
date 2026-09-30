classdef ControlBenchDiagnostics < handle
    % Collect diagnostics without changing the legacy solver value semantics.
    properties
        step = 0
        iteration = 0
        converged = false
        local = cell(0,6)
        residual = zeros(0,6)
    end
    methods
        function beginStep(obj, step)
            obj.step = step;
            obj.iteration = 0;
            obj.converged = false;
        end
        function recordSolve(obj, agent, problem, duration, message)
            obj.local(end+1,:) = {obj.step,obj.iteration,agent,problem,duration,message};
        end
        function recordResiduals(obj, agents, normalization)
            for i = 1:numel(agents)
                d = agents{i}.data;
                obj.residual(end+1,:) = [obj.step,obj.iteration,agents{i}.id, ...
                    d.primal_residual(end),d.dual_residual(end), ...
                    d.primal_residual(end)/sqrt(normalization(i))];
            end
        end
    end
end
