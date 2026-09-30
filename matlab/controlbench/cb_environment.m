function environment = cb_environment()
solver = which('ipopt');
fid = fopen(solver,'rb');
assert(fid~=-1,'ControlBench:IO','Cannot read solver binary.');
cleanup = onCleanup(@()fclose(fid)); %#ok<NASGU>
bytes = fread(fid,Inf,'*uint8');
hasher = java.security.MessageDigest.getInstance('SHA-256');
hasher.update(typecast(bytes,'int8'));
hash = lower(reshape(dec2hex(typecast(hasher.digest(),'uint8'),2)',1,[]));
environment = struct('matlab',version,'yalmip',yalmip('version'),'ipopt',solver, ...
    'ipopt_sha256',hash,'toolboxes',{ver});
end
