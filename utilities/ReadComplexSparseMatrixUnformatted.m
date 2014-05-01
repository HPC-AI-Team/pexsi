function A = ReadComplexSparseMatrixUnformatted(filename)
% Deserialize a SparseMatrix structure in csc format in complex
% arithmetic 

disp('Reading matrix from file (NOTE: ONLY SERIAL FORMAT WITH COMPLEX ARITHMETIC IS SUPPORTED)...');
disp('Reading matrix from file (The file can also be generated from the ParaWriteSparseMatrix procedure)...'); 
tic
fid = fopen(filename,'r');
N         = deserialize( fid, {'int'} );
Annz      = deserialize( fid, {'int'} );

colptr    = deserialize( fid, {'IntNumVec'} );
rowind    = deserialize( fid, {'IntNumVec'} );
Aval      = deserialize( fid, {'CpxNumVec'} );
fclose( fid );
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
