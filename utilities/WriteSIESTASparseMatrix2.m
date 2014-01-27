function nzval = WriteSIESTASparseMatrix2(A, filename, colptr, rowind)
% Write a sparse matrix A in compressed sparse column format that is
% compatible with SIESTA format. Both upper and lower triangular
% matrices are saved. The sparsity pattern is given by colptr, rowind
%
% Lin Lin
% 02/11/2012

disp('Computing the compressed sparse column format...');
tic
disp('norm(A-A'')_{inf} = ');
norm(A-A',inf)
A = (A+A')/2; % Make sure symmetry
N = length(A);
nnzA = length(rowind);
disp('Preparing nzval according to the sparsity pattern');
tic
	nzval = zeros( nnzA, 1 );
	cnt = 1;
	for i = 1 : nnzA 
		if( i == colptr( cnt + 1 ) )
			cnt = cnt + 1;
		end
		nzval(i) = A(rowind(i), cnt);
	end
toc
	

disp('Writing the matrix to file (SIESTA text format)...');
tic
	fid = fopen(filename,'w');
	fprintf(fid, '%d %d %d\n', [N, N, nnzA]);
	fprintf(fid, '%d ', colptr);
	fprintf(fid, '\n');
	fprintf(fid, '%d ', rowind);
	fprintf(fid, '\n');
	fprintf(fid, '%25.15f \n', nzval);
	fclose(fid);
toc
