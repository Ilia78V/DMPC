function dx = cb_dynamics(x,u,model)
% Uses one immutable state vector for every agent. No projection or clipping.
n = numel(model.actuation);
dx = model.actuation.*u;
for i = 1:n
    for j = 1:n
        if model.edges(i,j) ~= 0
            d = x(j)-x(i);
            dx(i) = dx(i)+model.edges(i,j)*(model.cubic*d^3+model.linear*d);
        end
    end
end
end
