function [colptr, rowind, nzval] = WriteSIESTASparseMatrix(A, filename)
% Write a sparse matrix A in compressed sparse column format that is
% compatible with the SIESTA format. Both upper and lower triangular
% matrices are saved.
%
% Lin Lin
% 02/11/2012

disp('Computing the compressed sparse column format...');
tic
disp('norm(A-A'')_{inf} = ');
norm(A-A',inf)
A = (A+A')/2; % Make sure symmetry
N = length(A);
nnzA = nnz( A );
[rowind, colind, nzval] = find( A );
[~,colptr,~] = unique( colind, 'first' );
colptr = [ colptr; nnzA+1 ];
toc

disp('Writing the matrix to file (SIESTA text format)...');
tic
	fid = fopen(filename,'w');
	fprintf(fid, '%d %d %d\n', [N, N, nnzA]);
	fprintf(fid, '%d ', colptr);
	fprintf(fid, '\n');
	fprintf(fid, '%d ', rowind);
	fprintf(fid, '\n');
	fprintf(fid, '%25.15f\n ', nzval);
	fclose(fid);
toc
