#ifndef _TASK_H_
#define _TASK_H_

#ifdef MPI_SUPPORT
#include "mpi.h"
#endif /* MPI_SUPPORT */

#include <vector>

class Task
{
    public:
        Task();
        Task(int i, int j, int nx, int ny);
        Task(const Task &obj);
        ~Task();

        Task& operator=(const Task& rhs);

        void setUpNeighbor(Task* up);
        void setLeftNeighbor(Task* left);
        void setDownNeighbor(Task* down);
        void setRightNeighbor(Task* right);
        
        void setMpiRank(int rank);
        int* getMpiRankPtr();

        void addTag(int tag);
        void addRrank(int* rank);
        void addRtask(Task* task); 

        double* getLocalGrid();
        double* getLocalNewGrid();

        void setLocalGrid(double* ptr);
        void setLocalNewGrid(double* ptr);

        int getUpNeighbor();
        int getDownNeighbor();
        int getLeftNeighbor();
        int getRightNeighbor();

        int getUpTag();
        int getDownTag();
        int getLeftTag();
        int getRightTag();

        void swapLocalGrids();

        void print();

    private:
        /*const*/ int i_;
        /*const*/ int j_;
        /*const*/ int nx_;
        /*const*/ int ny_;

        /*const*/ Task* upNeighbor_;
        /*const*/ Task* downNeighbor_;
        /*const*/ Task* leftNeighbor_;
        /*const*/ Task* rightNeighbor_;

        /*const*/ double* grid_;
        /*const*/ double* newGrid_;

        /*const*/ std::vector<int> tags_;

        int mpiRank_;

        std::vector<int*> rRanks_;
        std::vector<Task*> rTasks_;
};

#endif /* _TASK_H_ */
