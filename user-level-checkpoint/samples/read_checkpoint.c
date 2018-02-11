#include <stdio.h>
#include <stdlib.h>

#include <string.h>

#define INTEGRITY_SNAPSHOT_FILE "integrity_file.txt"

int get_last_snapshot_(char *last_checkpoint)
{
    //int myrank = get_comm_rank__();
    int myrank = 1;

    FILE *file = fopen(INTEGRITY_SNAPSHOT_FILE, "r");
    if (!file) {
        fprintf(stderr, "Can't read from %s\n", INTEGRITY_SNAPSHOT_FILE);
        exit(1);
    }

    // TO DO
    char * line = NULL;
    size_t len = 0;

    while ((getline(&line, &len, file)) != -1) {
        printf("%s\n", line);
        char tmp[10] = { 0 };
        sprintf(tmp, "%c", line[0]);
        if (atoi(tmp) == myrank) {
            printf("*%s\n", line);
            strcpy(last_checkpoint, line);
        }
    }

    fclose(file);

    int i;
    for (i = 0; i < strlen(last_checkpoint); i++) {
        if (last_checkpoint[i] == '=') {
            break;
        }
    }

    i += 1;

    char tmp_checkpoint[256] = { 0 };

    for (int j = 0; i < strlen(last_checkpoint); j++, i++) {
        tmp_checkpoint[j] = last_checkpoint[i];
    }

    strcpy(last_checkpoint, tmp_checkpoint);

    last_checkpoint[strlen(last_checkpoint) - 1] = '\0';

    printf("Rankd %d, file %s, phase %d\n", myrank, last_checkpoint, last_checkpoint[0] - '0');

    return last_checkpoint[0] - '0';
}

int main () 
{

    char last_chechkpoint[256] = { 0 };

    int phase = get_last_snapshot_(last_chechkpoint);

    printf("%d %s\n", phase, last_chechkpoint);


/*
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
*/
    return 0;
}
