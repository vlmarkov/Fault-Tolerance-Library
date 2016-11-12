#ifndef _CHECKPOINT_LIB_H_
#define _CHECKPOINT_LIB_H_

#include <stdio.h>
#include <sys/time.h>

#include <mpi.h>


/*****************************************************************************/
/* Global variables                                                          */
/*****************************************************************************/
extern double GLOBAL_START_TIME;

/*****************************************************************************/
/* Initialize globa start(time) point                                        */
/*****************************************************************************/
#define CHECKPOINT_LIB_INIT() GLOBAL_START_TIME = wtime_();


/*****************************************************************************/
/* Control flow-macros                                                       */
/*****************************************************************************/
#define CHECKPOINT_SET(name) name ## _checkpoint:
#define CHECKPOINT_GOTO(name) goto name ## _checkpoint


/*****************************************************************************/
/* Checkpoint-save macros                                                    */
/*****************************************************************************/
#define CHECKPOINT_SAVE(data, n, dtype, phase) make_snapshot(data, n, dtype, phase)


/*****************************************************************************/
/* Prototypes                                                                */
/*****************************************************************************/
void make_snapshot(void *data, int count, MPI_Datatype datatype, int phase);

double wtime_();

#endif /* _CHECKPOINT_LIB_H_ */
