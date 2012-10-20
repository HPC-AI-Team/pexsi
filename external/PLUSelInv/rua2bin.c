#include <stdio.h>

void dreadhb_serial(FILE *fp, int *nrow, int *ncol, int *nonz,
               double **nzval, int **rowind, int **colptr);

int main(int argc, char *argv[])
{
    char asc_fname[100], bin_fname[100];
    FILE *fasc, *fbin;
    int  nrow,    ncol,   nnz, n, i;
    int  *rowind, *colptr;
    double *nzvals;
 
    /* parse the command line, extract
       input and output file names */

    sscanf(*(argv+1),"%s", asc_fname);
    sscanf(*(argv+2),"%s", bin_fname);

    fasc = fopen(asc_fname,"r");
    dreadhb_serial(fasc, &nrow, &ncol, &nnz, &nzvals, &rowind, &colptr);
/*
    for (i=0; i<nrow+1; i++) colptr[i]++;
    for (i=0; i<nnz; i++) rowind[i]++;
*/
    fclose(fasc);
    n = nrow;

    fbin = fopen(bin_fname,"w");
    fwrite(&n,           4, 1,   fbin);
    fwrite(&nnz,         4, 1,   fbin);
    fwrite((int*)colptr, 4, n+1, fbin);
    for (i = 0; i < n; i++) {
       nrow = colptr[i+1]-colptr[i]; 
       fwrite(&rowind[colptr[i]], 4, nrow, fbin);
    }
    for (i = 0; i < n; i++) {
       nrow = colptr[i+1]-colptr[i]; 
       fwrite(&nzvals[colptr[i]], 8, nrow, fbin);
    }
    fclose(fbin);
    return 0;
}
