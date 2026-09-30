function root = cb_setup()
% Configure project dependencies for this MATLAB process only.
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(root);
paths = {getenv('CONTROLBENCH_YALMIP'),getenv('CONTROLBENCH_OPTI')};
fallback = {fullfile(root,'YALMIP-master'),fullfile(root,'OPTI-master')};
for i = 1:2
    if isempty(paths{i}), paths{i} = fallback{i}; end
    if isfolder(paths{i}), addpath(genpath(paths{i})); end
end
assert(exist('sdpvar','file') ~= 0,'ControlBench:Dependency','YALMIP is unavailable.');
assert(exist('ipopt','file') ~= 0,'ControlBench:Dependency','IPOPT is unavailable.');
end
