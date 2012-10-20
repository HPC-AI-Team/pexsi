% addpath('~/Software/hb_to_msm');
% %input_file='g20.rua';
% input_file='big.rua';
% 
% A = hb_to_msm ( input_file );

[fid,msg]=fopen('ROWPERM');
assert(isempty(msg));
N=fscanf(fid,'%d',1);
rowperm = fscanf(fid,'%d',inf)+1;
fclose(fid);

[fid,msg]=fopen('COLPERM');
assert(isempty(msg));
N=fscanf(fid,'%d',1);
colperm = fscanf(fid,'%d',inf)+1;
fclose(fid);

B=sparse(size(A));
B(colperm,colperm) = A;
