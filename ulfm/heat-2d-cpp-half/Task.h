/** @file */
#ifndef _TASK_H_
#define _TASK_H_

#ifdef MPI_SUPPORT
#include "mpi.h"
#endif /* MPI_SUPPORT */

#ifdef ULFM_SUPPORT
#include <mpi-ext.h>
#endif /* ULFM_SUPPORT */

#include <vector>

/**
 * @brief Available task states
 */
enum
{
    DEAD_TASK    = -1, /**< This task is dead */
    UNKNOWN_TASK = 0,  /**< This task is in unkonw state */
    ALIVE_TASK   = 1   /**< This task is alive */
};

/**
 * @brief High-level abstraction of the task.
 * @details Each task has a rank attached to it.
 * Describes interaction with other tasks (neighbors).
 * Describes the redundancy mechanism.
 * @author Markov V.A.
 * @date 2018
 * @version 0.1
 * @warning Work In Progress
 */
class Task {
    public:
        /**
         * @brief Default constructor
         */
        Task();

        /**
         * @brief Main constructor
         * @param i an integer argument Grid coodinate
         * @param j an integer argument Grid coodinate
         * @param nx an integer argument, number of columns in one Grid cell
         * @param ny an integer argument, number of columns in one Grid cell
         * @param repair an integer argument, number of avialable repairs
         */
        Task(int i, int j, int nx, int ny, int repair);

        /**
         * @brief Copy constructor
         * @param rhs a refernce Task object
         */
        Task(const Task &rhs);

        /**
         * @brief Destructor
         * @details Free memory for fields
         * @code
         * Task::grid_
         * Task::newGrid_
         * @endcode
         */
        ~Task();

        /**
         * @brief Assign operator
         * @param rhs a refernce Task object
         */
        Task& operator=(const Task& rhs);

        /**
         * @brief Sets Task MPI-rank
         * @param rank an integer argument new MPI-rank
         * @details Also sets Task::status - 'ALIVE_TASK'
         */
        void setMpiRank(int rank);

        /**
         * @brief Gets Task MPI-rank
         * @return integer MPI-rank
         */
        int getMpiRank();

        /**
         * @brief Gets Task MPI-rank
         * @return pointer to MPI-rank
         */
        int* getMpiRankPtr();

        /**
         * @brief Sets Task status
         * @param status an integer argument new Task status
         */
        void setStatus(int status);

        /**
         * @brief Gets Task status
         */
        int getStatus();

        /**
         * @brief Sets up neighbor
         * @param up a pointer argument new Task neighbor
         */
        void setUpNeighbor(Task* up);

        /**
         * @brief left up neighbor
         * @param left a pointer argument new Task neighbor
         */
        void setLeftNeighbor(Task* left);

        /**
         * @brief Sets down neighbor
         * @param down a pointer argument new Task neighbor
         */
        void setDownNeighbor(Task* down);

        /**
         * @brief Sets right neighbor
         * @param right a pointer argument new Task neighbor
         */
        void setRightNeighbor(Task* right);

        /**
         * @brief Gets up neighbor
         * @return pointer to Task object
         */
        Task* getUpNeighbor();

        /**
         * @brief Gets down neighbor
         * @return pointer to Task object
         */
        Task* getDownNeighbor();

        /**
         * @brief Gets left neighbor
         * @return pointer to Task object
         */
        Task* getLeftNeighbor();

        /**
         * @brief Gets right neighbor
         * @return pointer to Task object
         */
        Task* getRightNeighbor();

        /**
         * @brief Gets up neighbor MPI-rank
         * @param layer an integer argument redundancy layer
         * @return integer MPI-rank
         */
        int getUpNeighborRank(int layer);

        /**
         * @brief Gets down neighbor MPI-rank
         * @param layer an integer argument redundancy layer
         * @return integer MPI-rank
         */
        int getDownNeighborRank(int layer);

        /**
         * @brief Gets left neighbor MPI-rank
         * @param layer an integer argument redundancy layer
         * @return integer MPI-rank
         */
        int getLeftNeighborRank(int layer);

        /**
         * @brief Gets right neighbor MPI-rank
         * @param layer an integer argument redundancy layer
         * @return integer MPI-rank
         */
        int getRightNeighborRank(int layer);

        /**
         * @brief Sets local grid
         * @param ptr a pointer argument
         */
        void setLocalGrid(double* ptr);

        /**
         * @brief Sets new local grid
         * @param ptr a pointer argument
         */
        void setLocalNewGrid(double* ptr);

        /**
         * @brief Gets local grid
         * @param layer an integer argument redundancy layer
         * @return pointer to local grid
         */
        double* getLocalGrid(int layer);

        /**
         * @brief Gets new local grid
         * @param layer an integer argument redundancy layer
         * @return pointer to new local grid
         */
        double* getLocalNewGrid(int layer);

        /**
         * @brief Adds redundancy MPI-rank
         * @param rank a pointer argument new redundancy MPI-rank
         */
        void addRrank(int* rank);

        /**
         * @brief Adds redundancy Task
         * @param task a pointer argument new redundancy Task
         */
        void addRtask(Task* task);

        /**
         * @brief Adds Up MPI-tag
         * @param tag an integer argument new MPI-tag
         */
        void addUpTag(int tag);

        /**
         * @brief Adds Down MPI-tag
         * @param tag an integer argument new MPI-tag
         */
        void addDownTag(int tag);

        /**
         * @brief Adds Left MPI-tag
         * @param tag an integer argument new MPI-tag
         */
        void addLeftTag(int tag);

        /**
         * @brief Adds Right MPI-tag
         * @param tag an integer argument new MPI-tag
         */
        void addRightTag(int tag);

        /**
         * @brief Gets Up MPI-tag
         * @param layer an integer argument redundancy layer
         * @return integer MPI-tag
         */
        int getUpTag(int layer);

        /**
         * @brief Gets Down MPI-tag
         * @param layer an integer argument redundancy layer
         * @return integer MPI-tag
         */
        int getDownTag(int layer);

        /**
         * @brief Gets Left MPI-tag
         * @param layer an integer argument redundancy layer
         * @return integer MPI-tag
         */
        int getLeftTag(int layer);

        /**
         * @brief Gets Right MPI-tag
         * @param layer an integer argument redundancy layer
         * @return integer MPI-tag
         */
        int getRightTag(int layer);

        /**
         * @brief Swaps Task::grid_ and Task::newGrid_ fields
         * @param layer an integer argument redundancy layer
         */
        void swapLocalGrids(int layer);

        /**
         * @brief Gets number of redundancy layers
         */
        int getLayersNumber();

        /**
         * @brief Gets replacements vector
         */
        std::vector<Task*> getReplacements();

        /**
         * @brief Gets X coordinate
         * @param layer an integer argument redundancy layer
         * @throw std::string If can't get X
         */
        int getX(int layer);

        /**
         * @brief Gets Y coordinate
         * @param layer an integer argument redundancy layer
         * @throw std::string If can't get Y
         */
        int getY(int layer);

        /**
         * @brief Repairs task
         */
        void repair();

        /**
         * @brief Shows whole infomation about Task object
         */
        void print();

        /**
         * @brief Show whole infomation about Task object
         * by redundancy layers
         */
        void printByLayers();

    private:
        const int i_;                        /**< Coordinate in Grid */
        const int j_;                        /**< Coordinate in Grid */
        const int nx_;                       /**< Number of columns in one Grid cell */
        const int ny_;                       /**< Number of rows in one Grid cell */

        Task* upNeighbor_;                   /**< Up neighbor */
        Task* downNeighbor_;                 /**< Down neighbor */
        Task* leftNeighbor_;                 /**< Left neighbor */
        Task* rightNeighbor_;                /**< Right neighbor */

        std::vector<int> upNeighborTags_;    /**< Up neighbor MPI-tag */
        std::vector<int> downNeighborTags_;  /**< Down neighbor MPI-tag */
        std::vector<int> leftNeighborTags_;  /**< Left neighbor MPI-tag */
        std::vector<int> rightNeighborTags_; /**< Right neighbor MPI-tag */

        double* grid_;                       /**< Local sub-grid */
        double* newGrid_;                    /**< Local sub-grid */

        int mpiRank_;                        /**< Attached MPI-rank to Task */
        int status_;                         /**< Task status */
        int repair_;                         /**< Repair counter */

        std::vector<int*> rRanks_;           /**< Redundancy MPI-ranks vector */
        std::vector<Task*> rTasks_;          /**< Redundancy Task vector */
        std::vector<Task*> replacements_;    /**< Replacements Task vector */

        /**
         * @brief Adds new replacment Task
         * @param task a pointer argument Task object
         */
        void addReplacement_(Task* task);

        /**
         * @brief Gets next replacment Task
         * @return pointer to Task object
         */
        Task* getNextReplacement_();

        /**
         * @brief Gets next right MPI-rank
         * @param layer an integer argument redundancy layer
         * @return integer MPI-rank
         */
        int getNextRank_(int layer);

        /**
         * @brief Gets next up MPI-tag
         * @param layer an integer argument redundancy layer
         * @return integer MPI-tag
         */
        int getNextUpTag_(int layer);

        /**
         * @brief Gets next down MPI-tag
         * @param layer an integer argument redundancy layer
         * @return integer MPI-tag
         */
        int getNextDownTag_(int layer);

        /**
         * @brief Gets next left MPI-tag
         * @param layer an integer argument redundancy layer
         * @return integer MPI-tag
         */
        int getNextLeftTag_(int layer);

        /**
         * @brief Gets next right MPI-tag
         * @param layer an integer argument redundancy layer
         * @return integer MPI-tag
         */
        int getNextRightTag_(int layer);

        /**
         * @brief Reduces repair ability of Task instance
         */
        void reduceRepairAbility_();
};

#endif /* _TASK_H_ */
