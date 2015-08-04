function A = ReadSparseMatrixFormatted(filename)
% Read a sparse matrix in compressed column storage as produced by
% SIESTA.

disp('Reading matrix from file...');
tic
fid = fopen(filename,'r');
N = fscanf(fid, '%g', 1);
N = fscanf(fid, '%g', 1);
Annz = fscanf(fid, '%g', 1);
dummy  = fscanf(fid, '%g', 1);A=
colptr = fscanf(fid, '%g', N+1);
rowind = fscanf(fid, '%g', Annz);
Aval   = fscanf(fid, '%g', Annz);
fclose(fid);
toc

tic
disp('Converting column index...');
colind = zeros(Annz, 1);
cnt = 1;
for i =  1 : Annz
	if( i == colptr( cnt+1 ) )
		cnt = cnt + 1;
	end
  colind(i) = cnt;
end
toc

tic
disp('Generating A matrix...');
A = sparse(rowind, colind, Aval, N, N);
toc
