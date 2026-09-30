classdef TestControlBench < matlab.unittest.TestCase
    methods (TestClassSetup)
        function setupPath(testCase)
            root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(root,'matlab','controlbench')));
        end
    end
    methods (Test)
        function doubleIntegratorConstantInput(testCase)
            actual = cb_double_integrator([2;3],4,0.5);
            testCase.verifyEqual(actual,[4;5],'AbsTol',1e-12);
        end
        function doubleIntegratorEquilibrium(testCase)
            testCase.verifyEqual(cb_double_integrator([2;0],0,0.1),[2;0]);
        end
        function invalidInterval(testCase)
            testCase.verifyError(@()cb_double_integrator([0;0],0,-1), ...
                'MATLAB:validators:mustBePositive');
        end
        function conservativeCoupling(testCase)
            model = struct('actuation',[0;0;0],'edges',[0 1 0;1 0 2;0 2 0], ...
                'cubic',-0.01,'linear',0.2);
            testCase.verifyEqual(sum(cb_dynamics([1;2;0.5],[0;0;0],model)),0,'AbsTol',1e-12);
        end
        function permutationInvariantPlant(testCase)
            model = struct('actuation',[1;0;2],'edges',[0 1 0;1 0 2;0 2 0], ...
                'cubic',-0.01,'linear',0.2);
            permuted = model;
            permuted.actuation = model.actuation([3 1 2]);
            permuted.edges = model.edges([3 1 2],[3 1 2]);
            first = cb_dynamics([1;2;0.5],[0.2;0;0.3],model);
            second = cb_dynamics([0.5;1;2],[0.3;0.2;0],permuted);
            testCase.verifyEqual(second,first([3 1 2]),'AbsTol',1e-12);
        end
        function noStateProjection(testCase)
            model = struct('actuation',[1;1],'edges',zeros(2),'cubic',0,'linear',1);
            actual = [2.9;2.9]+0.1*cb_dynamics([2.9;2.9],[10;10],model);
            testCase.verifyGreaterThan(actual,3*ones(2,1));
        end
    end
end
