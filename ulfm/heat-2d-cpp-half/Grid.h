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
 * @brief High-level abstraction of the 2D calculation grid.
 * @details Represents a GRID of TASKS and the MPI-PROCESSES assigned to them.
 * @author Markov V.A.
 * @date 2018
 * @version 0.1
 * @warning Work In Progress
 */
class Grid {
public:
    /**
     * @brief Main constructor
     * @param columns an integer argument number of columns in calculation Grid
     * @param rows an integer argument, number of rows in calculation Grid
     * @param nx an integer argument, number of columns in one Grid cell
     * @param ny an integer argument, number of rows in one Grid cell
     * @param px an integer argument, number of processes by X axis
     * @param py an integer argument, number of processes by Y axis
     * @return Grid object
     * @details This method allocates Grid, allocates Task
     * (sets neighbors, MPI-ranks, MPI-tags, redundancy layers)
     */
    Grid(int cols, int rows, int nx, int ny, int px, int py);

    /** 
      * @brief Destructor
      * @param nothing
      * @return nothing
      */
    ~Grid();

    /**
     * @brief Gets Task object by MPI-rank
     * @param rank an integer argument, MPI-rank
     * @return pointer to Task object if found
     * @throw std::string If Task object not found
     */
    Task* getTask(int rank);

    /**
     * @brief Sets Task object status - 'DEAD_TASK'
     * @param rank an integer argument, MPI-rank
     * @return nothing
     * @throw std:string If reached the limit of reducibility
     * @details This method also reduces alive processes counter
     * and could shift left ranks
     * @code
     * Grid::shiftLeftMpiRank_(int rank);
     * @endcode
     */
    void kill(int rank);

    /**
     * @brief Repairs grid-tasks
     * @param nothing
     * @return nothing
     * @details This method finds dead Task object and invokes
     * @code
     * Task::repair()
     * @endcode
     * for it
     * @throw std::string If reached repair limit
     */
    void repair();

    /**
     * @brief Shows information about each Task in Grid
     * @param nothing
     * @return nothing
     */
    void print();

private:
    const int cols_;                        /**< Number of columns in calculation Grid */
    const int rows_;                        /**< Number of rows in calculation Grid */
    const int nx_;                          /**< Number of columns in one Grid cell */
    const int ny_;                          /**< Number of rows in one Grid cell */
    const int px_;                          /**< Number of MPI processes by X axis */
    const int py_;                          /**< Number of MPI processes by Y axis */
    int alive_;                             /**< Number of ALIVE MPI processes */
    std::vector<std::vector<Task> > tasks_; /**< Actual 2D Grid-Task */

    /**
     * @brief Sets neighbors for Task object
     * @param task a reference argument, Task object
     * @param i an integer argument, coordinate in Grid
     * @param j an integer argument, coordinate in Grid
     * @details Sets Up Down Left Right neighbor
     */
    void setNeighbors_(Task& task, int i, int j);

    /**
     * @brief Sets MPI-tags for Task object
     * @param task a reference argument, Task object
     * @param tag a reference argument, MPI-tag
     * @param layer an integer argument, redundancy layer
     * @details For each redundancy migration layer sets
     * unique MPI-tags betwen Task
     */
    void setTags_(Task& task, int& tag, int layer);

    /**
     * @brief Computes next coordinates for Task object
     * (where task will be migrate)
     * @param i a reference argument, grid coordinate
     * @param j a reference argument, grid coordinate
     * @details Redundancy migration
     */
    void computeNextCoordinates_(int& i, int& j);

    /**
     * @brief Links redundancy MPI-ranks, migrating Task objects
     * with Task object
     * @param task a poiner argument Task object
     * @param i an integer argument, Grid coordinate
     * @param j an integer argument, Grid coordinate
     * @throw std::string If can't gets self MPI-rank
     * @throw std::string If can't gets redundancy MPI-rank
     */
    void linkRanksTasks_(Task* task, int i, int j);

    /**
     * @brief Shifts left MPI-ranks
     * @param rank an integer argument, dead MPI-rank 
     * @return nothing
     * @details Updates all Tasks in Grid.
     * Shifts left (MPI-rank - 1) all MPI-ranks if DEAD MPI-rank > ALIVE MPI-rank
     */
    void shiftLeftMpiRank_(int rank);
};

#endif /* _GRID_H_ */
