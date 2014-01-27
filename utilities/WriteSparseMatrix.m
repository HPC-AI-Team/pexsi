function [colptr, rowind, nzval] = WriteSparseMatrix(A, filename)
% Write a sparse symmetric matrix A in compressed column storage. Only
% the lower triangular part (L) is saved 

disp('Computing the compressed column storage...');
tic
Alower = tril(A);
N = length(A);
nnzAlower = nnz( Alower );
[rowind, colind, nzval] = find( Alower );
[~,colptr,~] = unique( colind, 'first' );
colptr = [ colptr; nnzAlower+1 ];
toc

if(1)
	disp('Writing the matrix to file (text format)...');
	tic
		fid = fopen(filename,'w');
		fprintf(fid, '%d %d ', [N, nnzAlower]);
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
		serialize(fid, nnzAlower, {'int'});
		serialize(fid, colptr, {'IntNumVec'});
		serialize(fid, rowind, {'IntNumVec'});
		serialize(fid, nzval,  {'DblNumVec'});
		fclose(fid);
	toc
end
