/*****************************************************************************/
/* U - user                                                                  */
/* L - level                                                                 */
/* C - check                                                                 */
/* P - point                                                                 */
/*****************************************************************************/

#include "ulcp_header.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>

void** ulcp_checkpoint_table = NULL;
double* ulcp_base_snapshot   = NULL;
double ulcp_start_time       = 0.0;
double ulcp_start_time_local = 0.0;
double ulcp_save_time        = 0.0;
double ulcp_save_time_local  = 0.0;
int ulcp_snapshot_counter    = 0;
int ulcp_size                = 0;
int ulcp_time                = 0;
int ulcp_counter             = 0;
int ulcp_run_options         = 0;


void** ulcp_init_table(int size)
{
    void **jump_table = (void**)malloc(sizeof(void*) * size);
    if (!jump_table)
    {
        fprintf(stderr, "[%s] Can't allocate memory\n", __FUNCTION__);
        exit(1);
    }

    mkdir(SNAPSHOT_DIR_NAME, 0777);

    return jump_table;
}

void ulcp_init(int size, int argc, char *argv[])
{
    ulcp_counter          = 0;
    ulcp_size             = size;
    ulcp_start_time_local = ulcp_wtime();
    ulcp_checkpoint_table = ulcp_init_table(ulcp_size);

    if (ulcp_get_comm_rank() == 0)
    {
        printf("\n[CPL_LIBRARY] note, use: %s [args] [recovery]\n", argv[0]);
    }

    while (argc-->0)
    {
        if (strcmp(argv[argc], "recovery") == 0)
        {
            if (ulcp_get_comm_rank() == 0)
            {
                printf("[CPL_LIBRARY] running options 'recovery'\n");
            }
            ulcp_run_options = ULCP_RECOVER_MODE;
            break;
        }
    }
}

void ulcp_finalize()
{
    int i;

    // Zero process print statistic
    if (ulcp_get_comm_rank() == 0) {
        printf("\n");
        // Print up border
        for (i = 0; i < 80; i++) {
            printf("*");
        }
        printf("\n");
        ulcp_start_time += ulcp_wtime() - ulcp_start_time_local;
        printf("[CPL_LIBRARY] common elapsed time   = %f sec\n", ulcp_start_time);

        double time_one_save = ulcp_save_time / (double)ulcp_snapshot_counter;
        ulcp_save_time = time_one_save * (double)ulcp_snapshot_counter;

        printf("[CPL_LIBRARY] elapsed time one save = %f sec\n", time_one_save);
        printf("[CPL_LIBRARY] elapsed time all save = %f sec\n", ulcp_save_time);
        printf("[CPL_LIBRARY] amount of snapshots   = %d\n", ulcp_snapshot_counter);

        // Print down border
        for (i = 0; i < 80; i++)
        {
            printf("*");
        }
        printf("\n");
    }

    free(ulcp_base_snapshot);
}


int ulcp_is_recovery_mode()
{
    if (ulcp_run_options == ULCP_RECOVER_MODE)
    {
        return 1;
    }
    else
    {
        return 0;
    }
}
