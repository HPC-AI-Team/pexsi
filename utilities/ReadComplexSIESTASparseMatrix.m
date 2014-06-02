function A = ReadComplexSIESTASparseMatrix(filename)
% Read a sparse matrix in compressed column storage as produced by
% SIESTA.

disp('Reading matrix from file...');
tic
fid = fopen(filename,'r');
N = fscanf(fid, '%g', 1);
N = fscanf(fid, '%g', 1);
Annz = fscanf(fid, '%g', 1);
colptr = fscanf(fid, '%g', N+1);
drowind = fscanf(fid, '%g', Annz);
Aval   = fscanf(fid, '%g', 2*Annz);
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
ACpxval = Aval(1:2:end) + 1i*Aval(2:2:end);
A = sparse(rowind, colind, ACpxval, N, N);
toc
