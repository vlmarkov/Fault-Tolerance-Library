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
     * @param columns an integer argument number of columns in calculation grid
     * @param rows an integer argument, number of rows in calculation grid
     * @param nx an integer argument, number of columns in one grid cell
     * @param ny an integer argument, number of rows in one grid cell
     * @param px an integer argument, number of processes by X axis
     * @param py an integer argument, number of processes by Y axis
     * @return grid object
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
     * @return pointer to Task object
     * @throw std::string if Task object not found
     */
    Task* getTask(int rank);

    /**
     * @brief Marks Task object as 'dead-task'
     * @param rank an integer argument, MPI-rank
     * @return nothing
     */
    void kill(int rank);

    /**
     * @brief Repairs grid-tasks
     * @param nothing
     * @return nothing
     */
    void repair();

    /**
     * @brief Shows grid-tasks
     * @param nothing
     * @return nothing
     */
    void print();

private:
    
    const int cols_; /*!< Number of columns in calculation grid */
    const int rows_; /*!< Number of rows in calculation grid */
    const int nx_;   /*!< Number of columns in one grid cell */
    const int ny_;   /*!< Number of rows in one grid cell */
    const int px_;   /*!< Number of processes by X axis */
    const int py_;   /*!< Number of processes by Y axis */
    int alive_;      /*!< Number of ALIVE processes */
    std::vector<std::vector<Task> > tasks_; /*!< Actual calculation grid */

    /**
     * @brief Sets neighbors for Task object
     * @param task a reference argument, Task object
     * @param i an integer argument, coordinate in calculation grid
     * @param j an integer argument, coordinate in calculation grid
     */
    void setNeighbors_(Task& task, int i, int j);

    /**
     * @brief Sets MPI-rank for Task object
     * @param task a reference argument, Task object
     * @param rank an integer argument, MPI-rank
     * @details Also this task set status - 'ALIVE_TASK'
     */
    void setMpiRank_(Task& task, int rank);

    /**
     * @brief Sets MPI-tags for Task object
     * @param task a reference argument, Task object
     * @param tag a reference argument, MPI-tag
     * @param layer an integer argument, redundancy layer
     */
    void setTags_(Task& task, int& tag, int layer);

    /**
     * @brief Computes next coordinates for Task object (where task will be migrate)
     * @param i a reference argument, grid coordinate
     * @param j a reference argument, grid coordinate
     * @details Redundancy migration
     */
    void computeNextCoordinates_(int& i, int& j);

    /**
     * @brief Links redundancy MPI-ranks, migrating Task objects with Task object
     * @param task a poiner argument Task object
     * @param i an integer argument, grid coordinate
     * @param j an integer argument, grid coordinate
     * @throw std::string if can't gets self MPI-rank
     * @throw std::string if can't gets redundancy MPI-rank
     */
    void linkRanksTasks_(Task* task, int i, int j);

    /**
     * @brief Shifts left MPI-ranks
     * @param rank an integer argument, dead MPI-rank 
     * @return nothing
     * @throw nothing
     * @details Updates all Tasks in Grid. 
     * Shifts left (MPI-rank - 1) all MPI-ranks if DEAD MPI-rank > ALIVE MPI-rank
     */
    void shiftLeftMpiRank_(int rank);
};

#endif /* _GRID_H_ */
