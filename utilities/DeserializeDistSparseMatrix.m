function A = DeserializeDistSparseMatrix(filename)
% Deserialize a DistSparseMatrix structure generated from PEXSI
% (currently only serial format is supported).

disp('Reading matrix from file (NOTE: ONLY SERIAL FORMAT WITH COMPLEX ARITHMETIC IS SUPPORTED)...');
tic
fid = fopen(filename,'r');
N         = deserialize( fid, {'int'} );
Annz      = deserialize( fid, {'int'} );
AnnzLocal = deserialize( fid, {'int'} );
if( Annz ~= AnnzLocal )
	error('Must be generated from a serial calculation.');
end

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
