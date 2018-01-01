#ifndef _GRID_H_
#define _GRID_H_

#ifdef MPI_SUPPORT
#include "mpi.h"
#endif /* MPI_SUPPORT */

#ifdef ULFM_SUPPORT
#include <mpi-ext.h>
#endif /* ULFM_SUPPORT */

#include "Task.h"

#include <vector>

/**
 * High-level abstraction of the two-dimensional decomposition
 * of the computational domain for grid computations.
 * Represents a grid of tasks and the processes assigned to them.
 */
class Grid {
public:
    /**
     * Main constructor
     * @input: columns, rows,
     *         cols by x, rows by y,
     *         processes by x, processes by y
     */
    Grid(int cols, int rows, int nx, int ny, int px, int py);

    /**
     * Destructor
     */
    ~Grid();

    /**
     * Get task by MPI-rank
     * @input: MPI-rank
     * @return: pointer to Task
     */
    Task* getTask(int rank);

    /**
     * Mark task as 'dead'
     * @input: MPI-rank
     */
    void kill(int rank);

    /**
     * Repair grid-tasks
     */
    void repair();

    /**
     * Show grid-tasks
     */
    void print();

private:
    /**
     * Amount of columns
     */
    const int cols_;

    /**
     * Amount of rows
     */
    const int rows_;

    /**
     * Amount of cells by x (per task)
     */
    const int nx_;

    /**
     * Amount of cells by y (per task)
     */
    const int ny_;

    /**
     * Amount of processes by x
     */
    const int px_;

    /**
     * Amount of processes by y
     */
    const int py_;

    /**
     *Amount of alive processes
     */
    int alive_;

    /**
     * Actual grid-tasks
     */
    std::vector<std::vector<Task> > tasks_;

    /**
     * Set neighbors for task
     * @input: reference to task, coordinates 'i' and 'j'
     */
    void setNeighbors_(Task& task, int i, int j);

    /**
     * Set mpi-rank for task
     * @input: MPI-rank
     */
    void setMpiRank_(Task& task, int rank);

    /**
     * Set MPI-tags for task
     * @input: reference to task,
     *         reference to tag,
     *         redundancy layer
     */
    void setTags_(Task& task, int& tag, int layer);

    /**
     * Compute next cordiantes for task (where task will be migrate)
     * @input: coordinates 'i' and 'j'
     */
    void computeNextCoordinates_(int& i, int& j);

    /**
     * Link MPI-ranks, migrating tasks with task
     * @input: pointer to task, coordinates 'i' and 'j'
     */
    void linkRanksTasks_(Task* task, int i, int j);
};

#endif /* _GRID_H_ */
