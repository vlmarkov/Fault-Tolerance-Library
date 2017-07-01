#ifndef _ULCP_LIB_H_
#define _ULCP_LIB_H_

#include <mpi.h>
#include <mpi-ext.h> // ULFM support

#define ULCP_SNAPSHOT_DIR_NAME "snapshot"
#define ULCP_SNAPSHOT_INTEGRITY_FILE_MARK "\n=end_of_file="
#define ULCP_SNAPSHOT_INTEGRITY_FILE_NAME "integrity_file.txt"

enum {
    ULCP_SAVE_MODE    = 1,
    ULCP_RECOVER_MODE = 2
};

typedef struct ulcp_snapshot
{
    int            block;  // checkpoint block number
    int            offset;
    int            size;   // data size
    MPI_Datatype   type;   // data type
    void           *data;  // data itself
    struct DeltaCP *next;
} ulcp_t;

extern void** ulcp_checkpoint_table;
extern double* ulcp_base_snapshot;
extern double ulcp_start_time;
extern double ulcp_start_time_local;
extern double ulcp_save_time;
extern double ulcp_save_time_local;
extern int ulcp_snapshot_counter;
extern int ulcp_size;
extern int ulcp_time;
extern int ulcp_counter;
extern int ulcp_run_options;

/*****************************************************************************/
/* Declaration/Initialization checkpoint macros                              */
/*****************************************************************************/

/*
 * Description:
 * name - checkpoint_one, checkpoint_two, etc
 *      - phase_one, phase_two, etc
 *      - one, two, three, etc 
 */
#define ulcp_init_checkpoit(name)                                             \
    ulcp_checkpoint_table[ulcp_counter++] = name;                             \


/*****************************************************************************/
/* Control flow-macros                                                       */
/*****************************************************************************/
#define ulpc_goto_checkpoint(idx)                                             \
    goto *ulcp_checkpoint_table[idx];                                         \


#define ulpc_set_checkpoint(checkpoint_name)                                  \
    checkpoint_name :                                                         \


#define ulcp_save_data(checkpoint, user_save_callback)                        \
    user_save_callback(ulcp_get_snapshot_idx_by_name(ulcp_checkpoint_table,   \
                                                      ulcp_size, checkpoint));\


/*****************************************************************************/
/* Main functions                                                            */
/*****************************************************************************/
void ulcp_init(int size, int argc, char *argv[]);
void ulcp_finalize();

void** ulcp_init_table(int size);

/*****************************************************************************/
/* SAVE/RECOVER mode function                                                */
/*****************************************************************************/
int ulcp_is_recovery_mode();

#endif /* _ULCP_LIB_H_ */
