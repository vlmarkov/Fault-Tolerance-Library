#include <stdio.h>
#include <stdlib.h>

int main () 
{
    FILE * pFile;
    double *buffer = NULL;

    double ttotal, thalo, treduce;
    int niters;

    pFile = fopen("1_7_0.161894", "rb");
    if (!pFile) {
        fputs ("File error",stderr); 
        exit (1);
    }

    buffer = (double*)malloc(sizeof(double) * ((512 + 2) * (256 + 2)));
    if (!buffer) {
        fprintf(stderr, "can't allocate\n");
        exit(1);
    }

    // copy the file into the buffer:
    fread(buffer, sizeof(double), ((512 + 2) * (256 + 2)), pFile);
    fread(&ttotal, sizeof(double), 1, pFile);
    fread(&thalo, sizeof(double), 1, pFile);
    fread(&treduce, sizeof(double), 1, pFile);
    fread(&niters, sizeof(int), 1, pFile);


    for (int i = 0; i < ((512 + 2) * (256 + 2)); i++) {
        printf("%f ", buffer[i]);
    }

    printf("\n");


    printf("ttotal  %f\n", ttotal);
    printf("thalo   %f\n", thalo);
    printf("treduce %f\n", treduce);
    printf("niters  %d\n", niters);

    free(buffer);
    fclose (pFile);
    return 0;
}
