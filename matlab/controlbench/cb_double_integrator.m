function next = cb_double_integrator(state,input,dt)
% Exact zero-order-hold double integrator used as an independent test model.
arguments
    state (2,1) double {mustBeFinite}
    input (1,1) double {mustBeFinite}
    dt (1,1) double {mustBePositive,mustBeFinite}
end
next = [1 dt;0 1]*state + [dt^2/2;dt]*input;
end
