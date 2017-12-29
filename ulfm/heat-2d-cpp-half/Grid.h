#ifndef _GRID_H_
#define _GRID_H_

#ifdef MPI_SUPPORT
#include "mpi.h"
#endif /* MPI_SUPPORT */

#include "Task.h"

#include <vector>

class Grid
{
    public:
        Grid(int cols, int rows, int nx, int ny, int px, int py);
        ~Grid();

        void setNeighbors(Task& task, int i, int j);
        void setMpiRank(Task& task, int rank);
        void setTags(Task& task, int tag);

        Task* getTask(int rank);

        void linkRanksTasks(Task* task, int i, int j);

        void repair();
        void kill(int rank);

        void print();

    private:
        /*const*/ int cols_; // Amount of cells by x
        /*const*/ int rows_; // Amount of cells by y
        /*const*/ int nx_; // Amount of cells by x (per task)
        /*const*/ int ny_; // Amount of cells by y (per task)
        /*const*/ int px_; // Amount of processes by x
        /*const*/ int py_; // Amount of processes by y

        std::vector<std::vector<Task> > tasks_;

        void computeNextCoordinates_(int& i, int& j);
};

#endif /* _GRID_H_ */
