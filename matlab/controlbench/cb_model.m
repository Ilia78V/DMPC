function model = cb_model(c)
% Explicit scalar tank-chain model, including endpoint actuation and edge gains.
p = c.plant; n = p.agents;
model.actuation = zeros(n,1); model.actuation([1 n]) = 1/p.area_m2;
model.weights = p.state_weight*ones(n,1);
model.weights([1 n]) = p.endpoint_state_weight;
model.edges = zeros(n);
for i = 1:n-1
    gain = p.coupling_gain;
    if i == n-1, gain = gain*p.endpoint_gain; end
    model.edges(i,i+1)=gain; model.edges(i+1,i)=gain;
end
model.cubic = -0.0110497924; model.linear = 0.2131101400;
if strcmp(p.coupling_model,'linear')
    model.cubic = 0; model.linear = 1;
end
model.state_unit = 'm'; model.input_unit = 'm3/s';
model.coupling_valid_difference_m = [-3 3];
end
