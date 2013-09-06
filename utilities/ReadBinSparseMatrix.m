function A = ReadBinSparseMatrix( fileName, uplo, startIdx )

fh = fopen(fileName, 'rb');
sizeA = fread(fh,1,'int');
sizeTriplet = fread(fh, 1, 'int');
RowVec = fread(fh, sizeTriplet, 'int');
ColVec = fread(fh, sizeTriplet, 'int');
ValVec = fread(fh, sizeTriplet, 'double');
fclose(fh);
disp('Finish reading');


switch(startIdx)
	case 0
		idxPlus = 1;
	case 1
		idxPlus = 0;
	otherwise
		error('starting index can only be 0 or 1.');
end

A = sparse(RowVec + idxPlus, ColVec + idxPlus, ValVec, sizeA, sizeA);

switch(uplo)
	case {'u', 'U', 'l', 'L'}
		% If triangular matrix, get the full matrix
		A = A + A' - spdiags(diag(A), 0, sizeA, sizeA);
	case {'a','A'}
		% Do nothing
	otherwise
		error('Wrong uplo options.');
end


