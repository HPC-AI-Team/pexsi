function nzval = WriteSparseMatrixLU2(A, filename, colptr, rowind)
% Write a sparse matrix A in compressed sparse column format, for given
% colptr and rowind. It is used for LU factorization
%
% Lin Lin
% 12/11/2012

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
	

if(0)
	disp('Writing the matrix to file (text format)...');
	tic
		fid = fopen(filename,'w');
		fprintf(fid, '%d %d ', [N, nnzA]);
		fprintf(fid, '\n');
		fprintf(fid, '%d ', colptr);
		fprintf(fid, '\n');
		fprintf(fid, '%d ', rowind);
		fprintf(fid, '\n');
		fprintf(fid, '%g ', nzval);
		fclose(fid);
	toc
else
	disp('Writing the matrix to file (binary format)...');
	tic
		fid = fopen(filename,'wb');
		serialize(fid, N, {'int'} );
		serialize(fid, nnzA, {'int'});
		serialize(fid, colptr, {'IntNumVec'});
		serialize(fid, rowind, {'IntNumVec'});
		serialize(fid, nzval,  {'DblNumVec'});
		fclose(fid);
	toc
end
