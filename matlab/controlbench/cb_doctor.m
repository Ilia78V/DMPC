function cb_doctor()
cb_setup();
x = sdpvar(1);
d = optimize(x>=1,x^2,sdpsettings('solver','ipopt','verbose',0));
assert(d.problem==0,'ControlBench:Dependency','IPOPT smoke test failed: %s',d.info);
assert(abs(value(x)-1)<1e-4,'ControlBench:Dependency','Incorrect solver smoke result.');
environment = cb_environment();
out = getenv('CONTROLBENCH_OUTPUT');
if ~isempty(out)
    fid = fopen(fullfile(out,'matlab.json'),'w');
    assert(fid~=-1,'ControlBench:IO','Cannot write environment.');
    cleanup = onCleanup(@()fclose(fid)); %#ok<NASGU>
    fprintf(fid,'%s',jsonencode(environment));
end
disp(jsonencode(struct('matlab',environment.matlab,'yalmip',environment.yalmip,'status','PASS')));
end
