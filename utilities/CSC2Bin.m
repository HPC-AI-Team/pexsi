function A = CSC2Bin(infile, outfile)
% Read a sparse matrix in compressed sparse column format as produced by
% SIESTA, and then output it in binary format.

disp('Reading matrix from file...');
tic
fid = fopen(infile,'r');
N = fscanf(fid, '%g', 1);
N = fscanf(fid, '%g', 1);
Annz = fscanf(fid, '%g', 1);
colptr = fscanf(fid, '%g', N+1);
rowind = fscanf(fid, '%g', Annz);
Aval   = fscanf(fid, '%g', Annz);
fclose(fid);


disp('Writing the matrix to file (binary format)...');
tic
	fid = fopen(outfile,'wb');
	serialize(fid, N, {'int'} );
	serialize(fid, Annz, {'int'});
	serialize(fid, colptr, {'IntNumVec'});
	serialize(fid, rowind, {'IntNumVec'});
	serialize(fid, Aval,   {'DblNumVec'});
	fclose(fid);
toc
