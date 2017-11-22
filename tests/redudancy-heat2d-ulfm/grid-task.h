#ifndef _GRID_TASK_H_
#define _GRID_TASK_H_

#include "mpi.h"

typedef enum
{
    DEAD_PROC = -2,
    BORDER    = MPI_PROC_NULL, /* -1 */
} grid_task_e;

typedef struct task task_t;

struct task
{
    int     rank;          // Assigned proc rank
    int     top;           // Neighbor top
    int     bottom;        // Neighbor down
    int     left;          // Neighbor left
    int     right;         // Neighbor right
    int     redundancy_capacity;
    int     redundancy_counter;
    int     real_capacity;
    int     real_counter;
    task_t **real_task;
    task_t **redundancy_task;
    double *local_grid;
    double *local_newgrid;
}t;


typedef struct
{
    int      cols;          // Amount of cells by x (for whole grid)
    int      rows;          // Amount of cells by y (for whole grid)
    int      cols_per_proc; // Amount of processes by x (for whole grid)
    int      rows_per_proc; // Amount of processes by y (for whole grid)
    int      cols_per_task; // Amount of cells by x (per task)
    int      rows_per_task; // Amount of cells by y (per task)
    int      proc_per_node; // Amount of processes per node
    task_t **tasks;         // Grid[rows_per_proc][cols_per_proc]
} grid_task_t;

/*****************************************************************************/
/* Main destruction function                                                 */
/*****************************************************************************/
void grid_task_free(grid_task_t *grid);

/*****************************************************************************/
/* Main construction function                                                */
/*****************************************************************************/
grid_task_t *grid_task_allocate(const int cols,
                                const int rows,
                                const int nx,
                                const int ny,
                                const int commsize);

/*****************************************************************************/
/* Main intialize function                                                   */
/*****************************************************************************/
void grid_task_init(grid_task_t *grid);

/*****************************************************************************/
/* Redundancy ranks setter                                                   */
/*****************************************************************************/
void grid_task_redundancy_ranks_set(grid_task_t *grid,
                                    const int    rank,
                                    const int    row,
                                    const int    col);

/*****************************************************************************/
/* Grid getter                                                               */
/*****************************************************************************/
double *grid_task_local_grid_get(const grid_task_t *grid, const int rank);
double *grid_task_local_newgrid_get(const grid_task_t *grid, const int rank);

/*****************************************************************************/
/* Neighbor setters                                                          */
/*****************************************************************************/
void grid_task_neighbor_top_set(grid_task_t *grid, const int x, const int y);
void grid_task_neighbor_bottom_set(grid_task_t *grid, const int x, const int y);
void grid_task_neighbor_left_set(grid_task_t *grid, const int x, const int y);
void grid_task_neighbor_right_set(grid_task_t *grid, const int x, const int y);

/*****************************************************************************/
/* Neighbor getters                                                          */
/*****************************************************************************/
int grid_task_neighbor_top_get(const grid_task_t *grid, const int rank);
int grid_task_neighbor_bottom_get(const grid_task_t *grid, const int rank);
int grid_task_neighbor_left_get(const grid_task_t *grid, const int rank);
int grid_task_neighbor_right_get(const grid_task_t *grid, const int rank);

void grid_task_neighbors_get(const grid_task_t *grid,
                             const int          rank,
                             int               *top,
                             int               *bottom,
                             int               *left,
                             int               *right);


/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
void grid_task_redundancy_ranks_send_show(const grid_task_t *grid, const int rank);
void grid_task_redundancy_ranks_receive_show(const grid_task_t *grid, const int rank);


void grid_task_show(const grid_task_t *grid);
void grid_task_redundancy_show(const grid_task_t *grid);

int grid_task_redundancy_ranks_get(const grid_task_t *grid, const int rank, int *redundancy);

double *grid_task_redundancy_local_grid_get(const grid_task_t *grid, const int rank);
double *grid_task_redundancy_local_newgrid_get(const grid_task_t *grid, const int rank);

task_t *grid_task_get(const grid_task_t *grid, const int rank);

int grid_task_redundancy_task_get(const task_t *task,
                                  int          *ranks,
                                  double      **grid,
                                  double      **newgrid);

int grid_task_real_task_get(const task_t *my_task, task_t **tasks);

#endif /* _GRID_TASK_H_ */
