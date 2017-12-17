#include "grid-task.h"
#include "utils.h"

#include <cmath>
#include <cstdlib>

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include <iostream>

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
        std::cerr << "<" << __FUNCTION__ << ">"
                  << "  Not supproted grid size"
                  << std::endl;

#ifdef MPI_SUPPORT
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
#else
        exit(EXIT_FAILURE);
#endif /* MPI_SUPPORT */

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
     * Allocate memory for local_grid and local_newgrid
     */
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < m; j++)
        {
            task_t *task = &(this->tasks[i][j]);

            task->local_grid = (double *)
                xCalloc((ny + 2) * (nx + 2), sizeof(*task->local_grid));

            task->local_newgrid = (double *)
                xCalloc((ny + 2) * (nx + 2), sizeof(*task->local_newgrid));
        }
    }
}

GridTask::~GridTask()
{
    // TODO: memory free
}

void GridTask::init(grid_task_e mode)
{
    int addRank = 0;

    this->mode = mode;

    // Init each task
    for (int i = 0; i < this->cols_per_proc; i++)
    {
        for (int j = 0; j < this->rows_per_proc; j++)
        {
            this->tasks[i][j].rank = addRank++;
            this->tasks[i][j].x    = i;
            this->tasks[i][j].y    = j;

            this->tasks[i][j].redundancy.addReal(&(this->tasks[i][j]));
            this->tasks[i][j].redundancy.addRedundancy(&(this->tasks[i][j]));
        }
    }

    addRank = 0;

    // Set redundancy ranks for each task
    for (int i = 0; i < this->cols_per_proc; i++)
    {
        for (int j = 0; j < this->rows_per_proc; j++)
        {
            this->redundancyTaskSetBalanced(i, j, addRank++);
        }
    }

    addRank = 0;

    // Reorder redundancy ranks for each task
    for (int i = 0; i < this->cols_per_proc; i++)
    {
        for (int j = 0; j < this->rows_per_proc; j++)
        {
            this->redundancyTaskBalancedRedored(i, j, addRank++);
        }
    }

    // Set neighbors rank for each task
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
        this->tasks[x][y].top = GRID_TASK_BORDER;
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
        this->tasks[x][y].bottom = GRID_TASK_BORDER;
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
        this->tasks[x][y].left = GRID_TASK_BORDER;
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
        this->tasks[x][y].right = GRID_TASK_BORDER;
    }
    else if (y < (this->rows_per_proc - 1))
    {
        this->tasks[x][y].right = this->tasks[x][y + 1].rank;
    }
}

/*********************************************************************/
/* Neighbor getters                                                  */
/*********************************************************************/
int GridTask::neighborGet(const int rank, const int step)
{
    if (rank == GRID_TASK_BORDER)
    {
        return GRID_TASK_BORDER;
    }

    task_t *task = this->taskGet(rank);

    if (step > task->redundancy.getRealSize())
    {
        std::cerr << "<" << __FUNCTION__ << ">"
              << "Can't get neighbor" << std::endl;

#ifdef MPI_SUPPORT
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
#else
        exit(EXIT_FAILURE);
#endif /* MPI_SUPPORT */

    }

    std::vector<task_t *> neighbors = task->redundancy.getReal();

    return neighbors[step]->rank;
}

int GridTask::neighborGetTop(const task_t *t, const int step)
{
    int x = t->x;
    int y = t->y;

    if (x == 0)
    {
        return GRID_TASK_BORDER;
    }

    task_t *top = &this->tasks[x - 1][y];

    if (!top)
    {
        return GRID_TASK_BORDER;
    }
    else
    {
        if (step > top->redundancy.getRealSize())
        {
            std::cerr << "<" << __FUNCTION__ << ">"
                  << "Can't get neighbor" << std::endl;

#ifdef MPI_SUPPORT
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
#else
            exit(EXIT_FAILURE);
#endif /* MPI_SUPPORT */

        }

        std::vector<task_t *> neighbors = top->redundancy.getReal();

        return neighbors[step]->rank;
    }
}

int GridTask::neighborGetBottom(const task_t *t, const int step)
{
    int x = t->x;
    int y = t->y;

    if (x + 1 >= (int)this->tasks.size())
    {
        return GRID_TASK_BORDER;
    }

    task_t *bottom = &this->tasks[x + 1][y];

    if (!bottom)
    {
        return GRID_TASK_BORDER;
    }
    else
    {
        if (step > bottom->redundancy.getRealSize())
        {
            std::cerr << "<" << __FUNCTION__ << ">"
                  << "Can't get neighbor" << std::endl;

#ifdef MPI_SUPPORT
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
#else
            exit(EXIT_FAILURE);
#endif /* MPI_SUPPORT */

        }

        std::vector<task_t *> neighbors = bottom->redundancy.getReal();

        return neighbors[step]->rank;
    }
}

int GridTask::neighborGetRight(const task_t *t, const int step)
{
    int x = t->x;
    int y = t->y;

    if (y + 1 >= (int)this->tasks[x].size())
    {
        return GRID_TASK_BORDER;
    }

    task_t *right = &this->tasks[x][y + 1];

    if (!right)
    {
        return GRID_TASK_BORDER;
    }
    else
    {
        if (step > right->redundancy.getRealSize())
        {
            std::cerr << "<" << __FUNCTION__ << ">"
                  << "Can't get neighbor" << std::endl;

#ifdef MPI_SUPPORT
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
#else
            exit(EXIT_FAILURE);
#endif /* MPI_SUPPORT */

        }

        std::vector<task_t *> neighbors = right->redundancy.getReal();

        return neighbors[step]->rank;
    }
}

int GridTask::neighborGetLeft(const task_t *t, const int step)
{
    int x = t->x;
    int y = t->y;

    if (y == 0)
    {
        return GRID_TASK_BORDER;
    }

    task_t *left = &this->tasks[x][y - 1];

    if (!left)
    {
        return GRID_TASK_BORDER;
    }
    else
    {
        if (step > left->redundancy.getRealSize())
        {
            std::cerr << "<" << __FUNCTION__ << ">"
                  << "Can't get neighbor" << std::endl;

#ifdef MPI_SUPPORT
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
#else
            exit(EXIT_FAILURE);
#endif /* MPI_SUPPORT */

        }

        std::vector<task_t *> neighbors = left->redundancy.getReal();

        return neighbors[step]->rank;
    }
}


/*****************************************************************************/
/* Redundancy ranks setter                                                   */
/*****************************************************************************/
void GridTask::redundancyTaskSetSimple(const int row, const int col)
{
    int commsize = this->rows_per_proc * this->cols_per_proc;
    int addRank  = this->tasks[row][col].rank;

    for (int i = 0; i < 3; i++)
    {
        if (addRank - 4 >= 0)
        {
            addRank -= 4;
        }
        else
        {
            addRank -= 4;
            addRank = commsize + addRank;
        }
        this->tasks[row][col].redundancy.addReal(this->taskGet(addRank));
    }
}

void GridTask::redundancyTaskSetBalanced(const int row,
                                         const int col,
                                         const int addRank)
{
    const int x_times = this->cols_per_proc / this->proc_per_node;
    const int y_times = this->rows_per_proc / this->proc_per_node;

    int ri = row + this->proc_per_node;
    for (int x = 0; x < x_times; x++)
    {
        int rj = col + this->proc_per_node;

        ri = checkOverflow(ri, this->cols_per_proc);

        for (int y = 0; y < y_times; y++)
        {
            rj = checkOverflow(rj, this->rows_per_proc);

            /*
             * Do not include self
             */
            task_t *self = this->tasks[ri][rj].redundancy.getSelfTask();
            if (self->rank != addRank)
            {
                task_t *addTask = this->taskGet(addRank);
                this->tasks[ri][rj].redundancy.addReal(addTask);
                this->tasks[ri][rj].redundancy.addRedundancy(addTask);
            }

            rj += this->proc_per_node;
        }

        ri += this->proc_per_node;
    }
}

void GridTask::redundancyTaskBalancedRedored(const int row,
                                             const int col,
                                             const int reorderRank)
{
    static int invokeCnt = -1;

    const int commsize = this->rows_per_proc * this->cols_per_proc;

    invokeCnt++;

    if (reorderRank == (commsize / 2))
    {
        invokeCnt = -this->proc_per_node;
    }

    if (reorderRank + 1 > (commsize / 2))
    {
        this->tasks[row][col].redundancy.realReorder();
        this->tasks[row][col].redundancy.redundancyReorder();
    }

    if (reorderRank + 1 > (commsize / 2))
    {
        if (invokeCnt >= 0 && invokeCnt < this->proc_per_node)
        {
/*
            std::cout << "2st half Swap rank "
                      << reorderRank
                      << " cnt " << invokeCnt << std::endl;
*/
            this->tasks[row][col].redundancy.realSwapLast();
            this->tasks[row][col].redundancy.redundancySwapLast();

            if (invokeCnt + 2 > this->proc_per_node)
            {
                invokeCnt = -this->proc_per_node - 1;
            }
        }
    }
    else
    {
        if (invokeCnt >= 0 && invokeCnt < this->proc_per_node)
        {
/*
            std::cout << "1st half Swap rank "
                      << reorderRank
                      << " cnt " << invokeCnt << std::endl;
*/
            this->tasks[row][col].redundancy.realSwapLast();
            this->tasks[row][col].redundancy.redundancySwapLast();

            if (invokeCnt + 2 > this->proc_per_node)
            {
                invokeCnt = -this->proc_per_node - 1;
            }
        }
    }
/*
    std::cout << "invokeCnt " << invokeCnt << std::endl;
*/
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

    std::cerr << "<" << __FUNCTION__ << ">"
              << "Can't find task by rank "
              << rank
              << std::endl;

#ifdef MPI_SUPPORT
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
#else
        exit(EXIT_FAILURE);
#endif /* MPI_SUPPORT */

    return NULL; // Never reached
}

int GridTask::realTaskGet(task_t *my_task, task_t **tasks)
{
    if (!my_task)
    {
        std::cerr << "<" << __FUNCTION__ << ">"
                  << " Bad task pointer"
                  << std::endl;

#ifdef MPI_SUPPORT
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
#else
        exit(EXIT_FAILURE);
#endif /* MPI_SUPPORT */

    }

    int size = my_task->redundancy.getRealSize();
    std::vector<task_t *> t = my_task->redundancy.getReal();
    for (int i = 0; i < size; i++)
    {
        tasks[i] = t[i];
    }

    return size;
}

int GridTask::replaceTaskGet(task_t *my_task, task_t **tasks)
{
    if (!my_task)
    {
        std::cerr << "<" << __FUNCTION__ << ">"
                  << " Bad task pointer"
                  << std::endl;

#ifdef MPI_SUPPORT
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
#else
        exit(EXIT_FAILURE);
#endif /* MPI_SUPPORT */

    }

    int size = my_task->redundancy.getReplaceSize();
    std::vector<task_t *> t = my_task->redundancy.getReplace();
    for (int i = 0; i < size; i++)
    {
        tasks[i] = t[i];
    }

    return size;
}

void GridTask::show()
{
    for (int i = 0; i < this->cols_per_proc; i++)
    {
        for (int j = 0; j < this->rows_per_proc; j++)
        {
            std::cout << "Task(" << i << ", " << j << ")"
                      << " rank (mpi)          : "
                      << this->tasks[i][j].rank
                      << std::endl;

            std::cout << "Task(" << i << ", " << j << ")"
                      << " top rank (mpi)      : ";
            if (this->tasks[i][j].top == GRID_TASK_BORDER)
            {
                std::cout << "MPI_PROC_NULL";
            }
            else
            {
                std::cout << this->tasks[i][j].top;
            }
            std::cout << std::endl;

            std::cout << "Task(" << i << ", " << j << ")"
                      << " bottom rank (mpi)   : ";
            if (this->tasks[i][j].bottom == GRID_TASK_BORDER)
            {
                std::cout << "MPI_PROC_NULL";
            }
            else
            {
                std::cout << this->tasks[i][j].bottom;
            }
            std::cout << std::endl;

            std::cout << "Task(" << i << ", " << j << ")"
                      << " left rank (mpi)     : ";
            if (this->tasks[i][j].left == GRID_TASK_BORDER)
            {
                std::cout << "MPI_PROC_NULL";
            }
            else
            {
                std::cout << this->tasks[i][j].left;
            }
            std::cout << std::endl;

            std::cout << "Task(" << i << ", " << j << ")"
                      << " right rank (mpi)    : ";
            if (this->tasks[i][j].right == GRID_TASK_BORDER)
            {
                std::cout << "MPI_PROC_NULL";
            }
            else
            {
                std::cout << this->tasks[i][j].right;
            }
            std::cout << std::endl;

            std::cout << "Task(" << i << ", " << j << ")"
                      << " real tasks          : "
                      << this->tasks[i][j].redundancy.getRealSize()
                      << std::endl;

            std::cout << "Task(" << i << ", " << j << ")"
                      << " red. tasks          : "
                      << this->tasks[i][j].redundancy.getRedundancySize()
                      << std::endl;

            std::cout << "Task(" << i << ", " << j << ")"
                      << " rep. tasks          : "
                      << this->tasks[i][j].redundancy.getReplaceSize()
                      << std::endl;

            std::cout << "Task(" << i << ", " << j << ")"
                      <<  " red. task (vector)  : [ ";
            this->tasks[i][j].redundancy.printRedundancyRank();
            std::cout << "]" << std::endl;

            std::cout << "Task(" << i << ", " << j << ")"
                      <<  " real task (vector)  : [ ";
            this->tasks[i][j].redundancy.printRealRank();
            std::cout << "]" << std::endl;

            std::cout << "Task(" << i << ", " << j << ")"
                      <<  " real task* (vector) : [ ";
            this->tasks[i][j].redundancy.printRealRankDetail();
            std::cout << "]" << std::endl;

            std::cout << "Task(" << i << ", " << j << ")"
                      <<  " rep. task* (vector) : [ ";
            this->tasks[i][j].redundancy.printReplaceRankDetail();
            std::cout << "]" << std::endl;
            std::cout << std::endl;
        }
        std::cout << "---------------------------------------" << std::endl;
    }
    std::cout << std::endl;
}

void GridTask::repair()
{
    for (int i = 0; i < this->cols_per_proc; i++)
    {
        for (int j = 0; j < this->rows_per_proc; j++)
        {
            // Step 1 : Find dead processor
            if (this->tasks[i][j].rank == GRID_TASK_DEAD_PROC)
            {
                // Step 2 : Get pointer to dead task
                task_t *deadTask = &this->tasks[i][j];

                // Step 3 : Find replace in redundancy tasks
                deadTask->rank = deadTask->redundancy.repair();

                if (deadTask->rank == GRID_TASK_BORDER) // null processor
                {
                    std::cerr << "<" << __FUNCTION__ << ">"
                              << " Repair failed"
                              << std::endl;

#ifdef MPI_SUPPORT
                    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
#else
                    exit(EXIT_FAILURE);
#endif /* MPI_SUPPORT */

                }

                task_t *replaceTask = this->taskGet(deadTask->rank);
                replaceTask->redundancy.addReplace(deadTask);
            }
        }
    }

    // Step 4 : Reset neighbors rank
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

void GridTask::kill(const int rank)
{
    for (int i = 0; i < this->cols_per_proc; i++)
    {
        for (int j = 0; j < this->rows_per_proc; j++)
        {
            if (this->tasks[i][j].rank == rank)
            {
                this->tasks[i][j].rank = GRID_TASK_DEAD_PROC;
            }
        }
    }
}
