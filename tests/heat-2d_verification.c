#include <stdio.h>
#include <stdlib.h>

double *read_data(char *file_name)
{
    FILE * pFile;
    double *buffer = NULL;

    double ttotal, thalo, treduce;
    int niters;

    pFile = fopen(file_name, "rb");
    if (!pFile) {
        fprintf(stderr, "Can't open file %s\n", file_name);
        exit(1);
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

/*
    for (int i = 0; i < ((512 + 2) * (256 + 2)); i++) {
        printf("%f ", buffer[i]);
    }

    printf("\n");


    printf("ttotal  %f\n", ttotal);
    printf("thalo   %f\n", thalo);
    printf("treduce %f\n", treduce);
    printf("niters  %d\n", niters);
*/
    fclose (pFile);

    return buffer;
}

void data_verification(double *a, double *b)
{
    for (int i = 0; i < ((512 + 2) * (256 + 2)); i++) {
        if (a[i] != b[i]) {
            fprintf(stderr, "Verification failed\n");
        }
    }
    printf("Verification pass\n"); 
}

int main () 
{
    double *a;
    double *b;

    a = read_data("result/w_error/1_0_0.876316");
    b = read_data("result/wo_error/1_0_1.406019");

    data_verification(a, b);

    free(a);
    free(b);

    return 0;
}
