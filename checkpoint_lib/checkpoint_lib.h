#ifndef _CHECKPOINT_LIB_H_
#define _CHECKPOINT_LIB_H_

#include <stdio.h>
#include <signal.h>
#include <sys/time.h>

#include <mpi.h>

/*****************************************************************************/
/* C - check                                                                 */
/* P - point                                                                 */
/* L - library                                                               */
/*****************************************************************************/


/*****************************************************************************/
/* Global variables                                                          */
/*****************************************************************************/
extern double CPL_GLOBAL_START_TIME;
void   **CPL_GLOBAL_JUMP_TABLE;

int CPL_SIZE;
int CPL_TIME;
int CPL_COUNTER;

enum {
    CPL_CHECKPOINT_MODE = 0,
    CPL_RECOVERY_MODE   = 1
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

#define CPL_INIT(size, time, func)                                            \
    signal(SIGALRM, func);                                                    \
    CPL_TIME              = time;                                             \
    CPL_COUNTER           = 0;                                                \
    CPL_SIZE              = size;                                             \
    CPL_GLOBAL_START_TIME = wtime_();                                         \
    CPL_GLOBAL_JUMP_TABLE = init_table_(CPL_SIZE);                            \

//#define CPL_DEINIT() deinit_table_();

/*****************************************************************************/
/* Assigning checkpoint macros                                               */
/*****************************************************************************/

/*
 * Description:
 * name - checkpoint_one, checkpoint_two, etc
 *      - phase_one, phase_two, etc
 *      - one, two, three, etc 
 */

#define CPL_DECLARATE_CHECKPOINT(name)                                        \
    CPL_GLOBAL_JUMP_TABLE[CPL_COUNTER++] = name;                              \


/*****************************************************************************/
/* Control flow-macros                                                       */
/*****************************************************************************/
#define CPL_GO_TO_CHECKPOINT(idx)                                             \
    goto *CPL_GLOBAL_JUMP_TABLE[idx];                                             \


#define CPL_SET_CHECKPOINT(checkpoint_name)                                   \
    checkpoint_name :                                                         \


#define CPL_SAVE_STATE(name, callback)                                       \
    int _i_ = get_checkpoint_idx_by_name_(CPL_GLOBAL_JUMP_TABLE, CPL_SIZE, name); \
    callback(_i_);\


/*****************************************************************************/
/* Timer                                                                     */
/*****************************************************************************/
#define CPL_TIMER_INIT() timer_init_();                                       \

#define CPL_TIMER_STOP() timer_stop_();                                       \


/*****************************************************************************/
/* Checkpoint-save macros                                                    */
/*****************************************************************************/
#define CPL_FILE_OPEN(file, phase)                                            \
    open_checkpoint_file(file, phase);                                        \


#define CPL_FILE_CLOSE(file)                                                  \
    close_checkpoint_file(file);                                              \

#define CHECKPOINT_SAVE(file, data, n, type)                                  \
    make_snapshot(file, data, n, type);                                       \


#define CHECKPOINT_GET(last_chechkpoint) get_lastcheckpoint(last_chechkpoint);



/*****************************************************************************/
/* Prototypes                                                                */
/*****************************************************************************/
void make_snapshot(MPI_File file, void *data, int n, MPI_Datatype type);

int get_lastcheckpoint(char *last_chechkpoint);

void   **init_table_(int size);

double wtime_();

void   timer_init_();
void   timer_stop_();

void   close_checkpoint_file(MPI_File *snapshot);
void   open_checkpoint_file(MPI_File *snapshot, int phase);

int    get_checkpoint_idx_by_name_(void **table, int size, void *name);

#endif /* _CHECKPOINT_LIB_H_ */
