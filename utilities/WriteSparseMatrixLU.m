function [colptr, rowind, nzval] = ...
  WriteSparseMatrixLU(A, filename,formatted)
% Write a sparse matrix A in compressed sparse column format. It is used
% for LU factorization
%
% Lin Lin
% 10/22/2012

disp('Computing the compressed sparse column format...');
tic
% disp('norm(A-A'')_{inf} = ');
% norm(A-A',inf)
% A = (A+A')/2; % Make sure symmetry
N = length(A);
nnzA = nnz( A );
[rowind, colind, nzval] = find( A );
[~,colptr,~] = unique( colind, 'first' );
colptr = [ colptr; nnzA+1 ];
toc

if(formatted)
	disp('Writing the matrix to file (text format)...');
	tic
		fid = fopen(filename,'w');
		fprintf(fid, '%d %d %d %d', [N, N, nnzA, 0]);
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
		serialize(fid, N, {'int'} );
		serialize(fid, nnzA, {'int'});
		serialize(fid, 0, {'int'} );
		serialize(fid, colptr, {'IntNumVec'});
		serialize(fid, rowind, {'IntNumVec'});
		serialize(fid, nzval,  {'DblNumVec'});
		fclose(fid);
	toc
end
