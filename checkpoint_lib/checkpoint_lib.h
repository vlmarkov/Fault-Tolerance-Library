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
/* Initializing checkpoint macros                                            */
/*****************************************************************************/
#define CHECKPOINT_LIB_INIT(size, time)         \
    TIME = time;                                \
    COUNTER = 0;                                \
    GLOBAL_START_TIME = wtime_();               \
    void*  GLOBAL_JUMP_TABLE[size];             \


/*****************************************************************************/
/* Assigning checkpoint macros                                                  */
/*****************************************************************************/
#define CHECKPOINT_ASSIGN(label) GLOBAL_JUMP_TABLE[COUNTER++] = label;


/*****************************************************************************/
/* Control flow-macros                                                       */
/*****************************************************************************/
#define CHECKPOINT_GOTO(idx) goto *GLOBAL_JUMP_TABLE[idx];
#define CHECKPOINT_SET(name) name :


/*****************************************************************************/
/* Timer                                                                     */
/*****************************************************************************/
#define CHECKPOINT_SIGNAL_HANDLER(func) signal(SIGALRM, func);
#define CHECKPOINT_TIMER_INIT() timer_init_();


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

void close_checkpoint_file(MPI_File *snapshot);
void open_checkpoint_file(MPI_File *snapshot, int phase);

#endif /* _CHECKPOINT_LIB_H_ */
