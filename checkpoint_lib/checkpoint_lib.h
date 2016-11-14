#ifndef _CHECKPOINT_LIB_H_
#define _CHECKPOINT_LIB_H_

#include <stdio.h>
#include <sys/time.h>

//#include <mpi.h>

//#define TABLE_SIZE 2

/*****************************************************************************/
/* Global variables                                                          */
/*****************************************************************************/
extern double GLOBAL_START_TIME;

//void*  GLOBAL_JUMP_TABLE[TABLE_SIZE];
int TIME;
int COUNTER;

/*****************************************************************************/
/* Initialize global start(time) point                                        */
/*****************************************************************************/
#define CHECKPOINT_LIB_INIT(size, time)         \
    TIME = time;                                \
    COUNTER = 0;                                \
    GLOBAL_START_TIME = wtime_();               \
    void*  GLOBAL_JUMP_TABLE[size];             \


/*****************************************************************************/
/* Assing checkpoint macros                                                  */
/*****************************************************************************/
#define CHECKPOINT_ASSIGN(label) GLOBAL_JUMP_TABLE[COUNTER++] = label;


/*****************************************************************************/
/* Control flow-macros                                                       */
/*****************************************************************************/
#define CHECKPOINT_GOTO(idx) goto *GLOBAL_JUMP_TABLE[idx];

#define CHECKPOINT_SET(name) name :


#define CHECKPOINT_HANDLER(func) signal(SIGALRM, func);


#define CHECKPOINT_TIMER() timer_init_();

#define CHECKPOINT_GET(local_grid, ttotal, thalo, treduce, niters)                                  \
    get_snapshot(double *local_grid, double *ttotal, double *thalo, double *treduce, int *niters);  \

/*****************************************************************************/
/* Checkpoint-save macros                                                    */
/*****************************************************************************/
//#define CHECKPOINT_SAVE(data, n, dtype, phase)    \
    make_snapshot(data, n, dtype, phase);           \


/*****************************************************************************/
/* Prototypes                                                                */
/*****************************************************************************/
//void make_snapshot(void *data, int count, MPI_Datatype datatype, int phase);

double wtime_();
void timer_init_();

#endif /* _CHECKPOINT_LIB_H_ */
