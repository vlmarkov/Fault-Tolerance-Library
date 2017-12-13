#include "grid-task.h"
#include "utils.h"

#include <cmath>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>


GridTask::GridTask(const int cols,
                   const int rows,
                   const int nx,
                   const int ny,
                   const int commsize)
{
    const int n = std::sqrt(commsize);
    const int m = std::sqrt(commsize);
    const int p = std::sqrt(commsize) / 2;

    if ((n < 4) || (m < 4) || (p < 2))
    {
        fprintf(stderr, "<%s> Not supproted grid size\n", __FUNCTION__);
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    this->cols          = cols;
    this->rows          = rows;
    this->cols_per_task = nx;
    this->rows_per_task = ny;
    this->cols_per_proc = n;
    this->rows_per_proc = m;
    this->proc_per_node = p;

    // Allocate memory for each row in grid-tasks
    for (int i = 0; i < n; i++)
    {
        this->tasks.push_back(std::vector<task_t>(m));
    }

    /*
     * Allocate memory for local_grid and
     * local_newgrid and redundancy_task in each task
     */
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < m; j++)
        {
            task_t *task = &(this->tasks[i][j]);

            task->local_grid = (double *)
                xcalloc((ny + 2) * (nx + 2), sizeof(*task->local_grid));

            task->local_newgrid = (double *)
                xcalloc((ny + 2) * (nx + 2), sizeof(*task->local_newgrid));

            task->redundancy_task = std::vector<task_t *>(commsize, NULL);
            task->real_task = std::vector<task_t *>(commsize, NULL);
        }
    }
}

GridTask::~GridTask() { }

void GridTask::init(grid_task_e mode)
{
    int counter = 0;

    this->mode = mode;

    int commsize = (this->cols_per_proc * this->rows_per_proc);

    for (int i = 0; i < this->cols_per_proc; i++)
    {
        for (int j = 0; j < this->rows_per_proc; j++)
        {
            this->tasks[i][j].redundancy_counter  = 0;
            this->tasks[i][j].redundancy_capacity = commsize;
            this->tasks[i][j].real_capacity       = commsize;
            this->tasks[i][j].rank                = counter++;
            this->tasks[i][j].x                   = i;
            this->tasks[i][j].y                   = j;
            this->tasks[i][j].real_counter        = 1;
            this->tasks[i][j].real_task[0]        = &(this->tasks[i][j]);
        }
    }

    counter = 0;

    for (int i = 0; i < this->cols_per_proc; i++)
    {
        for (int j = 0; j < this->rows_per_proc; j++)
        {
            // Each struct task contains redudancy list
            this->redundancyTaskSet(counter++, i, j);
        }
    }

    // Neighbors rank set
    for (int i = 0; i < this->cols_per_proc; i++)
    {
        for (int j = 0; j < this->rows_per_proc; j++)
        {
            this->neighborTopSet(i, j);
            this->neighborBottomSet(i, j);
            this->neighborLeftSet(i, j);
            this->neighborRightSet(i, j);
        }
    }
}

/*****************************************************************************/
/* Neighbor setters                                                          */
/*****************************************************************************/
void GridTask::neighborTopSet(const int x, const int y)
{
    if (x == 0)
    {
        this->tasks[x][y].top = BORDER;
    }
    else if (x > 0)
    {
        this->tasks[x][y].top = this->tasks[x - 1][y].rank;
    }
}

void GridTask::neighborBottomSet(const int x, const int y)
{
    if (x == (this->cols_per_proc - 1))
    {
        this->tasks[x][y].bottom = BORDER;
    }
    else if (x < (this->cols_per_proc - 1))
    {
        this->tasks[x][y].bottom = this->tasks[x + 1][y].rank;
    }
}

void GridTask::neighborLeftSet(const int x, const int y)
{
    if (y == 0)
    {
        this->tasks[x][y].left = BORDER;
    }
    else if (y > 0)
    {
        this->tasks[x][y].left = this->tasks[x][y - 1].rank;
    }
}

void GridTask::neighborRightSet(const int x, const int y)
{
    if (y == (this->rows_per_proc - 1))
    {
        this->tasks[x][y].right = BORDER;
    }
    else if (y < (this->rows_per_proc - 1))
    {
        this->tasks[x][y].right = this->tasks[x][y + 1].rank;
    }
}

/*****************************************************************************/
/* Redundancy ranks setter                                                   */
/*****************************************************************************/
void GridTask::redundancyTaskSet(const int rank, const int row, const int col)
{
    const int x_times = this->cols_per_proc / this->proc_per_node;
    const int y_times = this->rows_per_proc / this->proc_per_node;

    int ri = row + this->proc_per_node;
    for (int x = 0; x < x_times; x++)
    {
        ri = check_overflow(ri, this->cols_per_proc);

        int rj = col + this->proc_per_node;
        for (int y = 0; y < y_times; y++)
        {
            rj = check_overflow(rj, this->rows_per_proc);

            const int capacity = this->tasks[ri][rj].redundancy_capacity;
            const int idx      = this->tasks[ri][rj].redundancy_counter;

            if (idx == capacity)
            {
                fprintf(stderr, "<%s> Bad redundancy-counter\n", __FUNCTION__);
                MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
            }

            this->tasks[ri][rj].redundancy_task[idx] = this->taskGet(rank);
            this->tasks[ri][rj].redundancy_counter  += 1;
            
            if (this->mode == COMPUTE_REDUNDANCY)
            {
                const int capacity = this->tasks[ri][rj].real_capacity;
                const int idx      = this->tasks[ri][rj].real_counter;

                if (idx == capacity)
                {
                    fprintf(stderr, "<%s> Bad redundancy-counter\n", __FUNCTION__);
                    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
                }

                /*
                 * Do not include self
                 */
                if (this->tasks[ri][rj].real_task[0]->rank != rank)
                {
                    this->tasks[ri][rj].real_task[idx] = this->taskGet(rank);
                    this->tasks[ri][rj].real_counter  += 1;
                }
            }

            rj += this->proc_per_node;
        }

        ri += this->proc_per_node;
    }
}

/*****************************************************************************/
/* Task getters                                                              */
/*****************************************************************************/
task_t *GridTask::taskGet(const int rank)
{
    for (int i = 0; i < this->cols_per_proc; i++)
    {
        for (int j = 0; j < this->rows_per_proc; j++)
        {
            if (rank == this->tasks[i][j].rank)
            {
                return &this->tasks[i][j];
            }
        }
    }

    fprintf(stderr, "<%s> Can't find task\n", __FUNCTION__);
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    return NULL; // Never reached
}

int GridTask::realTaskGet(const task_t *my_task, task_t **tasks)
{
    if (!my_task)
    {
        fprintf(stderr, "<%s> Bad task pointer\n", __FUNCTION__);
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    int size = my_task->real_counter;

    for (int i = 0; i < size; i++)
    {
        tasks[i] = my_task->real_task[i];
    }

    return size;
}

void GridTask::show()
{
    for (int i = 0; i < this->cols_per_proc; i++)
    {
        for (int j = 0; j < this->rows_per_proc; j++)
        {
            printf("Task->coordinates   : %04d:%04d\n", i, j);
            printf("Task->rank          : %04d\n", this->tasks[i][j].rank);
            printf("Task->top           : %04d\n", this->tasks[i][j].top);
            printf("Task->bottom        : %04d\n", this->tasks[i][j].bottom);
            printf("Task->left          : %04d\n", this->tasks[i][j].left);
            printf("Task->right         : %04d\n", this->tasks[i][j].right);
            printf("Task->red_capacity  : %04d\n", this->tasks[i][j].redundancy_capacity);
            printf("Task->red_counter   : %04d\n", this->tasks[i][j].redundancy_counter);
            printf("Task->real_capacity : %04d\n", this->tasks[i][j].real_capacity);
            printf("Task->real_counter  : %04d\n", this->tasks[i][j].real_counter);

            printf("Task->red_task      : { ");
            for (int r = 0; r < this->tasks[i][j].redundancy_counter; r++)
                printf("%04d ", this->tasks[i][j].redundancy_task[r]->rank);
            printf("}\n");

            printf("Task->real_task     : { ");
            for (int r = 0; r < this->tasks[i][j].real_counter; r++)
                printf("%04d ", this->tasks[i][j].real_task[r]->rank);
            printf("}\n"); 

            printf("Task->real_task*    : { ");
            for (int r = 0; r < this->tasks[i][j].real_counter; r++)
            {
                printf("%04d(%04d:%04d) ", this->tasks[i][j].real_task[r]->rank,
                    this->tasks[i][j].real_task[r]->x, this->tasks[i][j].real_task[r]->y);
            }
            printf("}\n"); 
            printf("\n");
        }
        printf("------------------------\n");
    }
    printf("\n");
}

void GridTask::repair()
{
    for (int i = 0; i < this->cols_per_proc; i++)
    {
        for (int j = 0; j < this->rows_per_proc; j++)
        {
            if (this->tasks[i][j].rank == DEAD_PROC)
            {
                task_t *dead_task = &this->tasks[i][j];
                for (int r = 0; r < dead_task->redundancy_counter; r++)
                {
                    task_t *repair_task = dead_task->redundancy_task[r];

                    if (repair_task->rank != DEAD_PROC)
                    {
                        const int counter  = repair_task->real_counter;
                        const int capacity = repair_task->real_capacity;
                        if (counter == capacity)
                        {
                            fprintf(stderr, "<%s> Repair failed\n", __FUNCTION__);
                            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
                        }

                        dead_task->rank = repair_task->rank;

                        if (this->mode == DATA_REDUNDANCY)
                        {
                            repair_task->real_task[counter] = dead_task;
                            repair_task->real_counter++;
                        }

                        break;
                    }
                }
            }
        }
    }

    // Neighbors rank set
    for (int i = 0; i < this->cols_per_proc; i++)
    {
        for (int j = 0; j < this->rows_per_proc; j++)
        {
            this->neighborTopSet(i, j);
            this->neighborBottomSet(i, j);
            this->neighborLeftSet(i, j);
            this->neighborRightSet(i, j);
        }
    }
}

void GridTask::killRank(const int killed)
{
    for (int i = 0; i < this->cols_per_proc; i++)
    {
        for (int j = 0; j < this->rows_per_proc; j++)
        {
            if (this->tasks[i][j].rank == killed)
            {
                this->tasks[i][j].rank = DEAD_PROC;
            }
        }
    }
}
