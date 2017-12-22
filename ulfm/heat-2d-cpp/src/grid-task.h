#ifndef _GRID_TASK_H_
#define _GRID_TASK_H_

#ifdef MPI_SUPPORT
#include "mpi.h"
#endif /* MPI_SUPPORT */

#include <vector>

typedef enum
{
    GRID_TASK_DEAD_PROC          = -999,
    GRID_TASK_DATA_REDUNDANCY    = -998,
    GRID_TASK_COMPUTE_REDUNDANCY = -997,
#ifdef MPI_SUPPORT
    GRID_TASK_BORDER             = MPI_PROC_NULL
#else
    GRID_TASK_BORDER             = -996
#endif /* MPI_SUPPORT */
} grid_task_e;

typedef struct task task_t;

class Redundancy
{
    public:
        Redundancy();
        ~Redundancy();

        task_t *getSelfTask();

        void addReal(task_t *t);
        void addReplace(task_t *t);
        void addRedundancy(task_t *t);

        void realReorder();
        void redundancyReorder();

        void redundancySort();

        void realSwapLast();
        void redundancySwapLast();

        std::vector<task_t *> getReal();
        std::vector<task_t *> getReplace();
        std::vector<task_t *> getRedundancy();

        int getRealSize();
        int getReplaceSize();
        int getRedundancySize();

        void printRealRank();
        void printRealRankDetail();
        void printRedundancyRank();
        void printReplaceRankDetail();

        int repair(int i, int j);

    private:
        std::vector<task_t *> real;
        std::vector<task_t *> replace;
        std::vector<task_t *> redundancy;
};

struct task
{
    int     x;                   // Grid coordinates
    int     y;                   // Grid coordinates
    int     rank;                // Assigned proc rank
    int     top;                 // Neighbor top
    int     bottom;              // Neighbor down
    int     left;                // Neighbor left
    int     right;               // Neighbor right

    double *local_grid;          // Local-grid
    double *local_newgrid;       // Local-grid (new)

    Redundancy redundancy;

};

class GridTask
{
    public:
        GridTask(const int cols,
                 const int rows,
                 const int nx,
                 const int ny,
                 const int commsize);

        ~GridTask();

        /*********************************************************************/
        /* Main intialize method                                             */
        /*********************************************************************/
        void init(grid_task_e mode);
        /*********************************************************************/
        /* Main kill process's task method                                   */
        /*********************************************************************/
        void kill(const int rank);
        /*********************************************************************/
        /* Main repair method                                                */
        /*********************************************************************/
        void repair();
        /*********************************************************************/
        /* Neighbor setters                                                  */
        /*********************************************************************/
        void neighborTopSet(const int x, const int y);
        void neighborBottomSet(const int x, const int y);
        void neighborLeftSet(const int x, const int y);
        void neighborRightSet(const int x, const int y);
        /*********************************************************************/
        /* Neighbor getter                                                   */
        /*********************************************************************/
        int neighborGet(const int rank, const int step);

        int neighborGetTop(const task_t *t, const int step);
        int neighborGetBottom(const task_t *t, const int step);
        int neighborGetLeft(const task_t *t, const int step);
        int neighborGetRight(const task_t *t, const int step);

        int neighborGetTopTag(const task_t *t, const int step);
        int neighborGetBottomTag(const task_t *t, const int step);
        int neighborGetLeftTag(const task_t *t, const int step);
        int neighborGetRightTag(const task_t *t, const int step);
        /*********************************************************************/
        /* Task getters                                                      */
        /*********************************************************************/
        task_t *taskGet(const int rank);
        /*********************************************************************/
        /* Redundancy task setters                                           */
        /*********************************************************************/
        void redundancyTaskSetSimple(const int row, const int col);
        void redundancyTaskSetBalanced(const int row,
                                       const int col,
                                       const int addRank);
        void redundancyTaskBalancedRedored(const int row,
                                           const int col,
                                           const int reorderRank);
        void redundancyTaskSort(const int row, const int col);
        /*********************************************************************/
        /* Raal task setter                                                  */
        /*********************************************************************/
        int realTaskGet(task_t *my_task, task_t **tasks);
        int replaceTaskGet(task_t *my_task, task_t **tasks);
        /*********************************************************************/
        /* Show/print method                                                 */
        /*********************************************************************/
        void show();
        void printGrid();
        void printGridRedundancy();

    private:
        grid_task_e mode;  // DATA or COMPUTE - redundancy

        int cols;          // Amount of cells by x (for whole grid)
        int rows;          // Amount of cells by y (for whole grid)
        int cols_per_proc; // Amount of processes by x (for whole grid)
        int rows_per_proc; // Amount of processes by y (for whole grid)
        int cols_per_task; // Amount of cells by x (per task)
        int rows_per_task; // Amount of cells by y (per task)
        int proc_per_node; // Amount of processes per node

        std::vector<std::vector<task_t> > tasks;

        void TwoDimDecomposition_(int &n, int &m, const int commsize);
};

#endif /* _GRID_TASK_H_ */
