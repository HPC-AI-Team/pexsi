function [colptr, rowind, nzval] = WriteSparseMatrixUnsym(A, filename)
% Write a sparse matrix A in compressed column storage. 

disp('Computing the compressed column storage...');
tic
N = length(A);
nnzA = nnz( A );
[rowind, colind, nzval] = find( A );
[~,colptr,~] = unique( colind, 'first' );
colptr = [ colptr; nnzA+1 ];
toc

if(1)
	disp('Writing the matrix to file (text format)...');
	tic
		fid = fopen(filename,'w');
		fprintf(fid, '%d %d %d %d ', [N,N, nnzA,0]);
		fprintf(fid, '\n');
		fprintf(fid, '%d ', colptr);
		fprintf(fid, '\n');
		fprintf(fid, '%d ', rowind);
		fprintf(fid, '\n');
    if(isreal(nzval))
  		fprintf(fid, '%g ', nzval);
    else
	    disp('Writing complex numbers...');
      tmp = [real(nzval.');imag(nzval.')];
      tmp=tmp(:)';
  		fprintf(fid, '%g ', tmp);
    end
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
