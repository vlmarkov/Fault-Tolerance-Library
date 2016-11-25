#ifndef _CHECKPOINT_LIB_H_
#define _CHECKPOINT_LIB_H_

#include <stdio.h>
#include <signal.h>
#include <sys/time.h>

#include <mpi.h>


/*****************************************************************************/
/* Global variables                                                          */
/*****************************************************************************/
extern double GLOBAL_START_TIME;

int TIME;
int COUNTER;

enum {
    CHECKPOINT_MODE = 0,
    RECOVERY_MODE   = 1
};


/*****************************************************************************/
/* Initializing checkpoint library macros                                    */
/*****************************************************************************/

/* 
 * Decription:
 * size - checkpoints numbers
 * time - in seconds for timer 
 * func - handler function
 */

#define CHECKPOINT_LIB_INIT(size, time, func)   \
    signal(SIGALRM, func);                      \
    TIME = time;                                \
    COUNTER = 0;                                \
    GLOBAL_START_TIME = wtime_();               \
    void*  GLOBAL_JUMP_TABLE[size];             \


/*****************************************************************************/
/* Assigning checkpoint macros                                               */
/*****************************************************************************/

/*
 * Description:
 * name - checkpoint_one, checkpoint_two, etc
 *      - phase_one, phase_two, etc
 *      - one, two, three, etc 
 */

#define DECLARATE_CHECKPOINT(name)            \
    GLOBAL_JUMP_TABLE[COUNTER++] = name;      \


/*****************************************************************************/
/* Control flow-macros                                                       */
/*****************************************************************************/
#define GO_TO_CHECKPOINT(idx)                 \
    goto *GLOBAL_JUMP_TABLE[idx];             \


#define SET_CHECKPOINT(checkpoint_name)       \
    checkpoint_name :                         \


#define SAVE_STATE(name)                                            \
    get_checkpoint_idx_by_name(GLOBAL_JUMP_TABLE, COUNTER, &&name); \


#define SAVE_STATE_AND_SET_CHECKPOINT(name)                         \
    get_checkpoint_idx_by_name(GLOBAL_JUMP_TABLE, COUNTER, &&name); \
    name :                                                          \


/*****************************************************************************/
/* Timer                                                                     */
/*****************************************************************************/
#define CHECKPOINT_TIMER_INIT() timer_init_();
#define CHECKPOINT_TIMER_STOP() timer_stop_();


/*****************************************************************************/
/* Checkpoint-save macros                                                    */
/*****************************************************************************/
#define CHECKPOINT_FILE_OPEN(file, phase) open_checkpoint_file(file, phase);
#define CHECKPOINT_FILE_CLOSE(file) close_checkpoint_file(file);

#define CHECKPOINT_SAVE(file, data, n, type) \
    make_snapshot(file, data, n, type);      \


#define CHECKPOINT_GET(last_chechkpoint) get_lastcheckpoint(last_chechkpoint);



/*****************************************************************************/
/* Prototypes                                                                */
/*****************************************************************************/
void make_snapshot(MPI_File file, void *data, int n, MPI_Datatype type);

int get_lastcheckpoint(char *last_chechkpoint);

double wtime_();
void timer_init_();
void timer_stop_();

void close_checkpoint_file(MPI_File *snapshot);
void open_checkpoint_file(MPI_File *snapshot, int phase);

int get_checkpoint_idx_by_name(void **table, int size, void *name);

#endif /* _CHECKPOINT_LIB_H_ */
