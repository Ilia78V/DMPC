classdef TestControlBenchIntegration < matlab.unittest.TestCase
    properties
        config
        output
    end
    properties (TestParameter)
        approximation = {'none','all'}
    end
    methods (TestMethodSetup)
        function setup(testCase)
            root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(root,'matlab','controlbench')));
            cb_setup();
            testCase.config = jsondecode(fileread(fullfile(root,'matlab','tests','fixture.json')));
            testCase.output = tempname;
            mkdir(testCase.output);
            testCase.addTeardown(@()rmdir(testCase.output,'s'));
        end
    end
    methods (Test)
        function exportGrid(testCase)
            path = fullfile(testCase.output,'config.json');
            writelines(jsonencode(testCase.config),path);
            cb_run(path,testCase.output);
            states = readtable(fullfile(testCase.output,'states.csv'));
            inputs = readtable(fullfile(testCase.output,'inputs.csv'));
            testCase.verifyEqual(height(states),6);
            testCase.verifyEqual(height(inputs),4);
            testCase.verifyEqual(max(states.time_s),0.2,'AbsTol',1e-12);
        end
        function infeasibleSolveRaises(testCase)
            c = testCase.config;
            c.plant.initial_level_m = 4;
            observer = ControlBenchDiagnostics(); observer.beginStep(0);
            testCase.verifyError(@()cb_centralized([4;4],c,cb_model(c),observer), ...
                'ControlBench:LocalSolveFailed');
            testCase.verifyNotEqual(observer.local{1,4},0);
        end
        function boundaryRoundoffRemainsObservable(testCase)
            c = testCase.config;
            observer = ControlBenchDiagnostics(); observer.beginStep(0);
            state = [3+2e-8;2.8];
            input = cb_centralized(state,c,cb_model(c),observer);
            testCase.verifyTrue(all(isfinite(input)));
            testCase.verifyGreaterThan(state(1),3);
        end
        function iterationLimitRecorded(testCase)
            c = testCase.config; c.controller.admm.max_iterations = 1;
            c.controller.admm.convergence_tolerance = 1e-12;
            observer = ControlBenchDiagnostics(); observer.beginStep(0);
            [solver,~] = cb_build_dmpc(c,cb_model(c),observer);
            solver.solve();
            testCase.verifyFalse(observer.converged);
            testCase.verifyEqual(observer.iteration,1);
            testCase.verifySize(observer.local,[2 6]);
        end
        function convexControllersAgree(testCase,approximation)
            c = testCase.config;
            c.controller.admm.max_iterations = 400;
            c.controller.admm.convergence_tolerance = 1e-7;
            c.controller.admm.initial_penalty = 1;
            c.controller.admm.approximation = approximation;
            observer = ControlBenchDiagnostics(); observer.beginStep(0);
            model = cb_model(c);
            centralized = cb_centralized([0.5;0.5],c,model,observer);
            observer.beginStep(0);
            [solver,agents] = cb_build_dmpc(c,model,observer);
            solver.solve();
            distributed = cellfun(@(a)value(a.data.u(:,1)),agents)';
            testCase.verifyTrue(observer.converged);
            testCase.verifyEqual(distributed,centralized,'AbsTol',0.005);
        end
    end
end
