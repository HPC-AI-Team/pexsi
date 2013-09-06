fileName = 'DGMAT_FULL_UPPER';
fh = fopen(fileName, 'rb');
sizeA = fread(fh,1,'int');
sizeTriplet = fread(fh, 1, 'int');
RowVec = fread(fh, sizeTriplet, 'int');
ColVec = fread(fh, sizeTriplet, 'int');
ValVec = fread(fh, sizeTriplet, 'double');

A = sparse(RowVec, ColVec, ValVec, sizeA, sizeA);

fclose(fh);

