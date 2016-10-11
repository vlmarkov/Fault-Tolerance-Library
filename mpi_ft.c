#include "mpi_ft.h"

int get_comm_rank()
{
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    return rank;
}

int get_comm_size()
{
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    return size;
}

int get_snapshot()
{
    int myrank = get_comm_rank();
    char buf[256] = { 0 };
    char file_name[256] = { 0 };

    sprintf(file_name,"%s_%d_snapshot.txt",__FILE__, myrank);

    FILE *snapshot = fopen(file_name, "r");

    if (snapshot) {
        while (fgets(buf, sizeof(buf), snapshot)) {
            printf("Rank %d restore: %s\n", myrank, buf);
        }
    }
    printf("Rankd %d sucsefully get data\n", myrank);
    fclose(snapshot);
    return atoi(buf);
}

void make_snapshot(int value)
{
    int rc;
    char file_name[256] = { 0 };
    int myrank = get_comm_rank();

    sprintf(file_name,"%s_%d_snapshot.txt",__FILE__, myrank);

    MPI_Status status;
    MPI_File snapshot;

    rc = MPI_File_open( MPI_COMM_WORLD,
                        file_name,
                        MPI_MODE_CREATE|MPI_MODE_WRONLY,
                        MPI_INFO_NULL,
                        &snapshot);

    if (rc == MPI_SUCCESS) {
        char buf[256] = { 0 };

        sprintf(buf, "%d", value);
        printf("%s:%d\n", buf, strlen(buf));

        //MPI_File_set_view(snapshot, 0,  MPI_CHAR, buf, "native", MPI_INFO_NULL);
        MPI_File_write(snapshot, buf, strlen(buf), MPI_CHAR, &status);
        MPI_File_sync(snapshot);

    } else {
        printf("ERROR: Rank %d can't save in file\n", myrank);
        return;
    }

    MPI_File_close(&snapshot);
}