#ifndef _MPI_FT_H_
#define _MPI_FT_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#include <signal.h>
#include <unistd.h>

#include <mpi.h>

int get_comm_rank();
int get_comm_size();
int get_snapshot();
void make_snapshot(int value);

#endif /* _MPI_FT_H_ */