#ifndef _TASK_H_
#define _TASK_H_

#ifdef MPI_SUPPORT
#include "mpi.h"
#endif /* MPI_SUPPORT */

#include <vector>

enum
{
    DEAD_TASK    = -1,
    UNKNOWN_TASK = 0,
    ALIVE_TASK   = 1
};

class Task
{
    public:
        Task();
        Task(int i, int j, int nx, int ny, int repair);
        Task(const Task &obj);
        ~Task();

        Task& operator=(const Task& rhs);

        void setMpiRank(int rank);
        int getMpiRank();
        int* getMpiRankPtr();

        void setStatus(int status);
        int getStatus();

        void setUpNeighbor(Task* up);
        void setLeftNeighbor(Task* left);
        void setDownNeighbor(Task* down);
        void setRightNeighbor(Task* right);
        Task* getUpNeighbor();
        Task* getDownNeighbor();
        Task* getLeftNeighbor();
        Task* getRightNeighbor();
        int getUpNeighborRank(int layer);
        int getDownNeighborRank(int layer);
        int getLeftNeighborRank(int layer);
        int getRightNeighborRank(int layer);

        void setLocalGrid(double* ptr);
        void setLocalNewGrid(double* ptr);
        double* getLocalGrid();
        double* getLocalNewGrid();

        void addRrank(int* rank);
        void addRtask(Task* task);

        void addUpTag(int tag);
        void addDownTag(int tag);
        void addLeftTag(int tag);
        void addRightTag(int tag);
        int getUpTag(int layer);
        int getDownTag(int layer);
        int getLeftTag(int layer);
        int getRightTag(int layer);

        void swapLocalGrids();

        void repair();

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

        /*const*/ std::vector<int> upNeighborTags_;
        /*const*/ std::vector<int> downNeighborTags_;
        /*const*/ std::vector<int> leftNeighborTags_;
        /*const*/ std::vector<int> rightNeighborTags_;

        /*const*/ double* grid_;
        /*const*/ double* newGrid_;

        int status_;
        int repair_;
        int mpiRank_;

        std::vector<int*> rRanks_;
        std::vector<Task*> rTasks_;

        Task* getNextRepair_();
        int getNextRank_(int layer);
        int getNextUpTag_(int layer);
        int getNextDownTag_(int layer);
        int getNextLeftTag_(int layer);
        int getNextRightTag_(int layer);
        void reduceRepairAbility_();
};

#endif /* _TASK_H_ */
