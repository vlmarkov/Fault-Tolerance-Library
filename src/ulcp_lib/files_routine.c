#include "ulcp.h"

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <sys/types.h>

void ulcp_open_file(MPI_File* file, int snapshot_phase)
{
    static int counter  = 0;
    char file_name[256] = { 0 };
    char file_path[256] = { 0 };
    /* 
     * 1_1 => [PHASE_OF_CALCULATION]_[COUNTER]
     *
     * PHASE_OF_CALCULATION - global state, need for 'goto'
     * COUNTER              - each PHASE_OF_CALCULATION could reach many times
     */

    sprintf(file_path,"%s/%d", ULCP_SNAPSHOT_DIR_NAME, ulcp_get_comm_rank());

    mkdir(file_path, 0777);

    sprintf(file_name,"%s/%d/%d_%d", ULCP_SNAPSHOT_DIR_NAME,
            ulcp_get_comm_rank(), snapshot_phase, counter++);

    MPI_File_open(MPI_COMM_WORLD, file_name,
                  MPI_MODE_CREATE|MPI_MODE_WRONLY,
                  MPI_INFO_NULL, file);

    ulcp_save_time_local = ulcp_wtime();
}

void ulcp_close_file(MPI_File* file)
{
    MPI_Status status;

    ulcp_snapshot_counter += 1;
    ulcp_save_time_local   = ulcp_wtime() - ulcp_save_time_local;
    ulcp_save_time        += ulcp_save_time_local;

    double elapsed_time = ulcp_wtime() - ulcp_start_time_local;

    // Write meta-data
    MPI_File_write(*file, &elapsed_time, 1, MPI_DOUBLE, &status);
    MPI_File_write(*file, &ulcp_save_time, 1, MPI_DOUBLE, &status);
    MPI_File_write(*file, &ulcp_snapshot_counter, 1, MPI_INT, &status);
    MPI_File_write(*file, &ULCP_SNAPSHOT_INTEGRITY_FILE_MARK, 
                    strlen(ULCP_SNAPSHOT_INTEGRITY_FILE_MARK), MPI_CHAR, &status);

    MPI_File_close(file);
}

void ulcp_save_file(MPI_File file, void* data, int size, MPI_Datatype type)
{
    MPI_Status status;
    MPI_File_write(file, data, size, type, &status);
}

FILE* ulcp_open_snapshot(char* file_name, char* mode)
{
    FILE *file = fopen(file_name, mode);
    if (!file) {
        fprintf(stderr, "Can't read %s\n", file_name);
        exit(1);
    }

    long offset = sizeof(ulcp_start_time) + sizeof(ulcp_save_time) \
                  + sizeof(ulcp_snapshot_counter) + strlen(ULCP_SNAPSHOT_INTEGRITY_FILE_MARK); 

    // Move file to the end of file
    fseek(file, -offset, SEEK_END);

    // copy the file into the buffer:
    fread(&ulcp_start_time, sizeof(double), 1, file);
    fread(&ulcp_save_time, sizeof(double), 1, file);
    fread(&ulcp_snapshot_counter, sizeof(int), 1, file);

    // Move file ptr to the begin of file
    rewind(file);

    return file;
}

int ulcp_get_snapshot(char* last_snapshot)
{
    int myrank = ulcp_get_comm_rank();

    char tmp_checkpoint[256] = { 0 };

    char * line = NULL;
    size_t len = 0;
    int i, j;

    FILE *file = fopen(ULCP_SNAPSHOT_INTEGRITY_FILE_NAME, "r");
    if (!file)
    {
        fprintf(stderr, "Can't read from %s\n", ULCP_SNAPSHOT_INTEGRITY_FILE_NAME);
        exit(1);
    }

    while ((getline(&line, &len, file)) != -1)
    {
        char tmp[10] = { 0 };
        sprintf(tmp, "%c", line[0]);

        if (atoi(tmp) == myrank)
        {
            strcpy(last_snapshot, line);
        }
    }

    fclose(file);

    for (i = 0; i < strlen(last_snapshot); i++)
    {
        if (last_snapshot[i] == '=')
        {
            break;
        }
    }

    i += 1;

    for (j = 0; i < strlen(last_snapshot); j++, i++)
    {
        tmp_checkpoint[j] = last_snapshot[i];
    }

    strcpy(last_snapshot, tmp_checkpoint);

    last_snapshot[strlen(last_snapshot) - 1] = '\0';

    printf("Rank %d, file %s, phase %d\n", myrank, last_snapshot, last_snapshot[0] - '0');

    return last_snapshot[0] - '0';
}

int ulcp_get_snapshot_idx_by_name(void **table, int size, void *name)
{
    int i;
    for (i = 0; i < size; i++)
    {
        if (table[i] == name)
        {
            break;
        }
    }
    return i;
}
