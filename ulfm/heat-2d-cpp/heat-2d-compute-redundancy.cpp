/***************************************************************************************/
/* heat-2d.c: MPI implementation of Laplace equation solver by Jacobi iteration method.*/
/*                                                                                     */ 
/* 2D Laplace equation:                                                                */
/*   \Delta u = 0                                                                      */
/*   \frac{\partial^2 u}{\partial x^2} + \frac{\partial^2 u}{\partial y^2} = 0         */
/*                                                                                     */
/* Domain: x in [0, 1],  y in [0, 1]                                                   */
/* Boundary conditions:                                                                */
/*   u(x, 0) = sin(pi * x)                                                             */
/*   u(x, 1) = sin(pi * x) * exp(-pi)                                                  */
/*   u(0, y) = u(1, y) = 0                                                             */
/* Initial value for interior points is 0                                              */
/* Analytical solution:                                                                */
/*   u(x, y) = sin(pi * x) * exp(-pi * y)                                              */
/*                                                                                     */
/* Parallel implementation:                                                            */
/* 2D domain decomposition of grid [0..rows - 1][0..cols -1]                           */
/* Each process is assigned a subgrid [rows / nprocs][cols / nprocs]                   */
/*                                                                                     */
/* Input parameters: rows, cols, EPS                                                   */
/*                                                                                     */
/* Usage: mpiexec -np <p> ./heat-2d <rows> <cols>                                      */
/*                                                                                     */
/* (C) Mikhail Kurnosov, 2015                                                          */
/***************************************************************************************/

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <signal.h>
#include <unistd.h>
#include <inttypes.h>

#include <sys/stat.h>
#include <sys/types.h>

#include <mpi.h>
#include <mpi-ext.h>

#include "src/utils.h"
#include "src/grid-task.h"

#include <algorithm>


#define EPS       0.001
#define PI        3.14159265358979323846
#define NELEMS(x) (sizeof((x)) / sizeof((x)[0]))
#define IND(i, j) ((i) * (nx + 2) + (j))

FILE *fp;

typedef struct
{
    int          count;
    int          nx;
    int          ny;
    MPI_Datatype row;
    MPI_Datatype col;
    MPI_Comm     comm;
    task_t      *task;
    task_t      *neighbors;
    MPI_Request *reqs;
} halo_cookie_t;

static void updateInteriorPoints(double    *newGrid,
                                 double    *oldGrid,
                                 const int  ny,
                                 const int  nx)
{
    for (int i = 1; i <= ny; i++)
    {
        for (int j = 1; j <= nx; j++)
        {
            newGrid[IND(i, j)] = (oldGrid[IND(i - 1, j)] +
                                  oldGrid[IND(i + 1, j)] +
                                  oldGrid[IND(i, j - 1)] + 
                                  oldGrid[IND(i, j + 1)]) * 0.25;
        }
    }
}

static double checkTerminationCondition(double    *newGrid,
                                        double    *oldGrid,
                                        const int  ny,
                                        const int  nx)
{
    double maxdiff = 0;

    for (int i = 1; i <= ny; i++)
    {
        for (int j = 1; j <= nx; j++)
        {
            int ind = IND(i, j);
            maxdiff = fmax(maxdiff, fabs(oldGrid[ind] - newGrid[ind]));
        }
    }

    return maxdiff;
}

static int haloExchange(halo_cookie_t *cookie,
                        GridTask      &gridTask,
                        int            myRank,
                        int            cnt,
                        int            step,
                        char           mode)
{
    fprintf(fp, "<%s> Invoke\n", __FUNCTION__);
    fflush(fp);

    const int tag = 0;

    task_t *task      = cookie->task;
    task_t *neighbors = cookie->neighbors;
    double *buf       = task->local_grid;
    int count         = cookie->count;
    int nx            = cookie->nx;
    int ny            = cookie->ny;
    MPI_Datatype row  = cookie->row;
    MPI_Datatype col  = cookie->col;
    MPI_Comm comm     = cookie->comm;
    MPI_Request *reqs = cookie->reqs;

    /*
     * Pay attention to MPI-ranks who send or recv messages
     */
/*
    int top    = (neighbors->top == GRID_TASK_BORDER)   ? GRID_TASK_BORDER : task->top;
    int bottom = (neighbors->bottom == GRID_TASK_BORDER)? GRID_TASK_BORDER : task->bottom;
    int left   = (neighbors->left == GRID_TASK_BORDER)  ? GRID_TASK_BORDER : task->left;
    int right  = (neighbors->right == GRID_TASK_BORDER) ? GRID_TASK_BORDER : task->right;
*/
/*
    int top    = neighbors->top;
    int bottom = neighbors->bottom;
    int left   = neighbors->left;
    int right  = neighbors->right;
*/

    int top    = gridTask.neighborGetTop(task, step);
    int bottom = gridTask.neighborGetBottom(task, step);
    int left   = gridTask.neighborGetLeft(task, step);
    int right  = gridTask.neighborGetRight(task, step);

    int tagT   = gridTask.neighborGetTopTag(task, step);
    int tagB   = gridTask.neighborGetBottomTag(task, step);
    int tagL   = gridTask.neighborGetLeftTag(task, step);
    int tagR   = gridTask.neighborGetRightTag(task, step);
/*
    fprintf("<%s> <*%c> <%d> MPI-rank(%02d), Task(%02d)(%02d, %02d), "
           "Neighbors: U(%02d) D(%02d) L(%02d) R(%02d)\n",
             __FUNCTION__, mode,  step, myRank, task->rank, task->x, task->y,
            task->top, task->bottom, task->left, task->right);
*/
    fprintf(fp, "<%s> <%c> <%d> MPI-rank(%02d), Task(%02d)(%02d, %02d), "
           "Neighbors: U(%02d(%02d)) D(%02d(%02d)) L(%02d(%02d)) R(%02d(%02d))\n",
             __FUNCTION__, mode,  step, myRank, task->rank, task->x, task->y,
            top, tagT, bottom, tagB, left, tagL, right, tagR);
/*
    // Non-blocking recv from TOP neighbor
    MPI_Irecv(&buf[IND(0, 1)],      count, row, top,    tag, comm, &reqs[cnt++]);
    // Non-blocking recv from BOTTOM neighbor
    MPI_Irecv(&buf[IND(ny + 1, 1)], count, row, bottom, tag, comm, &reqs[cnt++]);
    // Non-blocking recv from LEFT neighbor
    MPI_Irecv(&buf[IND(1, 0)],      count, col, left,   tag, comm, &reqs[cnt++]);
    // Non-blocking recv from RIGHT neighbor
    MPI_Irecv(&buf[IND(1, nx + 1)], count, col, right,  tag, comm, &reqs[cnt++]);

    // Non-blocking send to TOP neighbor
    MPI_Isend(&buf[IND(1, 1)],      count, row, top,    tag, comm, &reqs[cnt++]);
    // Non-blocking send to BOTTOM neighbor
    MPI_Isend(&buf[IND(ny, 1)],     count, row, bottom, tag, comm, &reqs[cnt++]);
    // Non-blocking send to LEFT neighbor
    MPI_Isend(&buf[IND(1, 1)],      count, col, left,   tag, comm, &reqs[cnt++]);
    // Non-blocking send to RIGHT neighbor
    MPI_Isend(&buf[IND(1, nx)],     count, col, right,  tag, comm, &reqs[cnt++]);
*/

    // Non-blocking recv from TOP neighbor
    MPI_Irecv(&buf[IND(0, 1)],      count, row, top,    tagT, comm, &reqs[cnt++]);
    // Non-blocking recv from BOTTOM neighbor
    MPI_Irecv(&buf[IND(ny + 1, 1)], count, row, bottom, tagB, comm, &reqs[cnt++]);
    // Non-blocking recv from LEFT neighbor
    MPI_Irecv(&buf[IND(1, 0)],      count, col, left,   tagR, comm, &reqs[cnt++]);
    // Non-blocking recv from RIGHT neighbor
    MPI_Irecv(&buf[IND(1, nx + 1)], count, col, right,  tagL, comm, &reqs[cnt++]);

    // Non-blocking send to TOP neighbor
    MPI_Isend(&buf[IND(1, 1)],      count, row, top,    tagT, comm, &reqs[cnt++]);
    // Non-blocking send to BOTTOM neighbor
    MPI_Isend(&buf[IND(ny, 1)],     count, row, bottom, tagB, comm, &reqs[cnt++]);
    // Non-blocking send to LEFT neighbor
    MPI_Isend(&buf[IND(1, 1)],      count, col, left,   tagL, comm, &reqs[cnt++]);
    // Non-blocking send to RIGHT neighbor
    MPI_Isend(&buf[IND(1, nx)],     count, col, right,  tagR, comm, &reqs[cnt++]);


    fprintf(fp, "<%s> Done\n", __FUNCTION__);
    fflush(fp);

    return cnt;
}

static void errorHandler(MPI_Comm *pcomm, int err, GridTask &gridTask)
{
    fprintf(fp, "<%s> Invoke\n", __FUNCTION__);
    fflush(fp);

    MPI_Comm comm = *pcomm;

    char errStr[MPI_MAX_ERROR_STRING] = { 0 };

    int rank   = 0;
    int size   = 0;
    int nf     = 0;
    int len    = 0;
    int eclass = 0;

    MPI_Group groupC, groupF;

    int *ranksGc = NULL;
    int *ranksGf = NULL;

    MPI_Error_class(err, &eclass);

    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    MPIX_Comm_failure_ack(comm);
    MPIX_Comm_failure_get_acked(comm, &groupF);

    MPI_Group_size(groupF, &nf);
    MPI_Error_string(err, errStr, &len);

    if (MPI_ERR_PROC_FAILED  != eclass)
    {
        std::cerr << "Rank " << rank
                  << " got missmatch error "
                  << errStr << std::endl;


        return;
        //MPI_Abort(comm, err);
    }

    ranksGf = (int*)malloc(nf * sizeof(int));
    ranksGc = (int*)malloc(nf * sizeof(int));

    MPI_Comm_group(comm, &groupC);

    for (int i = 0; i < nf; i++)
    {
        ranksGf[i] = i;
    }

    MPI_Group_translate_ranks(groupF,
        nf, ranksGf, groupC, ranksGc);

    std::cout << "Rank " << rank << "/" << size - 1
            << ": Notified of error " << errStr << ". "
            << nf << " found dead: [ ";

    fprintf(fp, "Rank %04d/%04d Notified of error %s.%d found dead: [",
       rank, size - 1, errStr, nf);
    fflush(fp);

    for (int i = 0; i < nf; i++)
    {
        std::cout << ranksGc[i] << " ";
        fprintf(fp, "%04d ", ranksGc[i]);
        fflush(fp);
        gridTask.kill(ranksGc[i]);
        gridTask.repair();
    }

    fprintf(fp, "]\n");
    fflush(fp);
    std::cout << "]" << std::endl;

    free(ranksGc);
    free(ranksGf);

    fprintf(fp, "<%s> Done\n", __FUNCTION__);
    fflush(fp);
}

static MPI_Comm repairCommunicator(MPI_Comm comm, int err)
{
    fprintf(fp, "<%s> Invoke\n", __FUNCTION__);
    fflush(fp);

    int allSucceeded = (err == MPI_SUCCESS);
    MPIX_Comm_agree(comm, &allSucceeded);
    if (!allSucceeded)
    {
        MPI_Comm commSpare;
        MPIX_Comm_revoke(comm);
        MPIX_Comm_shrink(comm, &commSpare); // Shrink the communicator
        MPI_Comm_free(&comm);               // Release the communicator
        comm = commSpare;                   // Assign shrink
    }
    else
    {
        fprintf(fp, "<%s><%d> Can't repair communicator\n", __FUNCTION__, __LINE__);
        fflush(fp);

        std::cerr << "Can't repair communicator" << std::endl;
        MPI_Abort(comm, err);
    }

    fprintf(fp, "<%s> Done\n", __FUNCTION__);
    fflush(fp);

    return comm;
}

int main(int argc, char *argv[])
{
    int rank = 0;
    int commsize = 0;

    MPI_Init(&argc, &argv);


    double ttotal = -MPI_Wtime();

    /*
     * Create a new error handler for MPI_COMM_WORLD
     * This overrides the default MPI_ERRORS_ARE_FATAL
     * so that ranks in this communicator will not
     * automatically abort if a failure occurs.
     */
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Comm_set_errhandler(comm, MPI_ERRORS_RETURN);

    MPI_Comm_size(comm, &commsize);
    MPI_Comm_rank(comm, &rank);

    /*
     * Create 2D grid of processes: commsize = px * py
     */
    MPI_Comm cartcomm;
    int dims[2]     = {0, 0};
    int periodic[2] = {0, 0};

    MPI_Dims_create(commsize, 2, dims);
    int px = dims[0];
    int py = dims[1];

    if (px < 2 || py < 2)
    {
        fprintf(stderr, "Invalid number of processes %d: px %d and py %d"
                        "must be greater than 1\n", commsize, px, py);
        MPI_Abort(comm, EXIT_FAILURE);
    }

    MPI_Cart_create(comm, 2, dims, periodic, 0, &cartcomm);
    int coords[2];
    MPI_Cart_coords(cartcomm, rank, 2, coords);
    int rankx = coords[0];
    int ranky = coords[1];

    int rows, cols;

    /* 
     * Broadcast command line arguments
     */
    if (rank == 0)
    {
        rows = (argc > 1) ? atoi(argv[1]) : py * 100;
        cols = (argc > 2) ? atoi(argv[2]) : px * 100;

        if (rows < py)
        {
            fprintf(stderr, "Number of rows %d less then number of py processes %d\n", rows, py);
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        }

        if (cols < px)
        {
            fprintf(stderr, "Number of cols %d less then number of px processes %d\n", cols, px);
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        }

        int args[2] = {rows, cols};
        MPI_Bcast(&args, NELEMS(args), MPI_INT, 0, MPI_COMM_WORLD);
    }
    else
    {
        int args[2];
        MPI_Bcast(&args, NELEMS(args), MPI_INT, 0, MPI_COMM_WORLD);
        rows = args[0];
        cols = args[1];
    }

    /* 
     * Allocate memory for local 2D subgrids with halo cells
     * [0..ny + 1][0..nx + 1]
     */
    int ny = getBlockSize(rows, ranky, py);
    int nx = getBlockSize(cols, rankx, px);

    GridTask gridTask(cols, rows, nx, ny, commsize);

    /*
     * Each proccessor will be have 1 + n (task to compute)
     * - n by defaul is 3 tasks
     */
    gridTask.init(GRID_TASK_COMPUTE_REDUNDANCY);

    task_t *myTask = gridTask.taskGet(rank);

    task_t *realTask[commsize];
    task_t *replaceTask[commsize];

    memset(realTask, 0, sizeof(task_t *) * commsize);
    memset(replaceTask, 0, sizeof(task_t *) * commsize);

    int realCounter    = gridTask.realTaskGet(myTask, realTask);
    int replaceCounter = gridTask.replaceTaskGet(myTask, replaceTask);

    for (int i = 0; i < realCounter; ++i)
    {
        double *grid    = realTask[i]->local_grid;
        double *newgrid = realTask[i]->local_newgrid;
        rankx           = realTask[i]->x;
        ranky           = realTask[i]->y;

        /*
         * Fill boundary points: 
         *   - left and right borders are zero filled
         *   - top border: u(x, 0) = sin(pi * x)
         *   - bottom border: u(x, 1) = sin(pi * x) * exp(-pi)
         */
        double dx = 1.0 / (cols - 1.0); 
        int    sj = getSumOfPrevBlocks(cols, rankx, px);

        if (ranky == 0)
        {
            // Initialize top border: u(x, 0) = sin(pi * x)
            for (int j = 1; j <= nx; j++)
            {
                // Translate col index to x coord in [0, 1]
                double x     = dx * (sj + j - 1);
                int ind      = IND(0, j);
                newgrid[ind] = grid[ind] = sin(PI * x);
            }
        }

        if (ranky == py - 1)
        {
            // Initialize bottom border: u(x, 1) = sin(pi * x) * exp(-pi)
            for (int j = 1; j <= nx; j++)
            {
                // Translate col index to x coord in [0, 1]
                double x     = dx * (sj + j - 1);
                int ind      = IND(ny + 1, j);
                newgrid[ind] = grid[ind] = sin(PI * x) * exp(-PI);
            }
        }

        // Do sync
        realTask[i]->local_grid    = grid;
        realTask[i]->local_newgrid = newgrid;
    }

    // Restore ranks for statistic
    rankx = realTask[0]->x;
    ranky = realTask[0]->y;

    /*
     * Left and right borders type
     */
    MPI_Datatype col;
    MPI_Type_vector(ny, 1, nx + 2, MPI_DOUBLE, &col);
    MPI_Type_commit(&col);

    /*
     * Top and bottom borders type
     */
    MPI_Datatype row;
    MPI_Type_contiguous(nx, MPI_DOUBLE, &row);
    MPI_Type_commit(&row);

    double thalo   = 0;
    double treduce = 0;

    int niters = 0;
    double maxdiff = 0.0;
    int rc;

    if (rank == 0)
    {
        //gridTask.show();
    }



    char filename[256] = { 0 };
    mkdir("logs", 0777);
    sprintf(filename, "logs/log-%d-rank.txt", rank);

    fp = fopen(filename, "w+");
    fprintf(fp, "Start calculation\n");

    while (1)
    {
        niters++;
restart_step:

        if (niters > 1000)
        {
            std::cerr << "Overflow accepted iterations" << std::endl;
            fprintf(fp, "Overflow accepted iterations\n");
            break;
        }

        maxdiff = 0.0;

        // Loop by for process's tasks
        for (int i = 0; i < realCounter; i++)
        {
            // Step 1: Update interior points
            updateInteriorPoints(realTask[i]->local_newgrid,
                                 realTask[i]->local_grid,
                                 ny, nx);

            // Step 2: Check termination condition
            double tmp = checkTerminationCondition(
                realTask[i]->local_newgrid,
                realTask[i]->local_grid,
                ny, nx);

            if (tmp > maxdiff)
            {
                maxdiff = tmp;
            }

            /*
             * Step 3: Swap grids
             * (after termination local_grid will contain result)
             */
            std::swap(realTask[i]->local_grid,
                      realTask[i]->local_newgrid);
        }

        /*
         * Killing test
         */

        if (niters == 3 && rank == 15)
        {
            std::cerr << "Triggered SIGKILL for " << rank << " rank" << std::endl;
            fprintf(fp, "<%s><%d>Triggered SIGKILL\n", __FUNCTION__, __LINE__);
            fclose(fp);
            raise(SIGKILL);
        }

        if (niters == 3 && rank == 14)
        {
            std::cerr << "Triggered SIGKILL for " << rank << " rank" << std::endl;
            fprintf(fp, "<%s><%d>Triggered SIGKILL\n", __FUNCTION__, __LINE__);
            fclose(fp);
            raise(SIGKILL);
        }

        if (niters == 3 && rank == 13)
        {
            std::cerr << "Triggered SIGKILL for " << rank << " rank" << std::endl;
            fprintf(fp, "<%s><%d>Triggered SIGKILL\n", __FUNCTION__, __LINE__);
            fclose(fp);
            raise(SIGKILL);
        }

        if (niters == 3 && rank == 12)
        {
            std::cerr << "Triggered SIGKILL for " << rank << " rank" << std::endl;
            fprintf(fp, "<%s><%d>Triggered SIGKILL\n", __FUNCTION__, __LINE__);
            fclose(fp);
            raise(SIGKILL);
        }

        if (niters == 5 && rank == 11)
        {
            std::cerr << "Triggered SIGKILL for " << rank << " rank" << std::endl;
            fprintf(fp, "<%s><%d>Triggered SIGKILL\n", __FUNCTION__, __LINE__);
            fclose(fp);
            raise(SIGKILL);
        }

        if (niters == 5 && rank == 10)
        {
            std::cerr << "Triggered SIGKILL for " << rank << " rank" << std::endl;
            fprintf(fp, "<%s><%d>Triggered SIGKILL\n", __FUNCTION__, __LINE__);
            fclose(fp);
            raise(SIGKILL);
        }

        if (niters == 5 && rank == 9)
        {
            std::cerr << "Triggered SIGKILL for " << rank << " rank" << std::endl;
            fprintf(fp, "<%s><%d>Triggered SIGKILL\n", __FUNCTION__, __LINE__);
            fclose(fp);
            raise(SIGKILL);
        }

        if (niters == 5 && rank == 8)
        {
            std::cerr << "Triggered SIGKILL for " << rank << " rank" << std::endl;
            fprintf(fp, "<%s><%d>Triggered SIGKILL\n", __FUNCTION__, __LINE__);
            fclose(fp);
            raise(SIGKILL);
        }


        if (niters == 7 && rank == 7)
        {
            std::cerr << "Triggered SIGKILL for " << rank << " rank" << std::endl;
            fprintf(fp, "<%s><%d>Triggered SIGKILL\n", __FUNCTION__, __LINE__);
            fclose(fp);
            raise(SIGKILL);
        }

        if (niters == 7 && rank == 6)
        {
            std::cerr << "Triggered SIGKILL for " << rank << " rank" << std::endl;
            fprintf(fp, "<%s><%d>Triggered SIGKILL\n", __FUNCTION__, __LINE__);
            fclose(fp);
            raise(SIGKILL);
        }

        if (niters == 7 && rank == 3)
        {
            std::cerr << "Triggered SIGKILL for " << rank << " rank" << std::endl;
            fprintf(fp, "<%s><%d>Triggered SIGKILL\n", __FUNCTION__, __LINE__);
            fclose(fp);
            raise(SIGKILL);
        }

        if (niters == 7 && rank == 2)
        {
            std::cerr << "Triggered SIGKILL for " << rank << " rank" << std::endl;
            fprintf(fp, "<%s><%d>Triggered SIGKILL\n", __FUNCTION__, __LINE__);
            fclose(fp);
            raise(SIGKILL);
        }

        // Step 4: All reduce (may fail)
        treduce -= MPI_Wtime();
        rc = MPI_Allreduce(MPI_IN_PLACE, &maxdiff, 1, MPI_DOUBLE, MPI_MAX, comm);
        treduce += MPI_Wtime();

        /*
         * Collective operation may fail so
         * ULFM provide repair operation
         */
        if (MPI_ERR_PROC_FAILED == rc)
        {
            fprintf(fp, "<%s><%d> Error\n", __FUNCTION__, __LINE__);
            fflush(fp);
            // Handle error
            errorHandler(&comm, rc, gridTask);

            fprintf(fp, "<%s><%d>\n", __FUNCTION__, __LINE__);
            fflush(fp);

            // Repair communicator
            comm = repairCommunicator(comm, rc);

            fprintf(fp, "<%s><%d>\n", __FUNCTION__, __LINE__);
            fflush(fp);

            // Repair grid-task
            //gridTask.repair();

            fprintf(fp, "<%s><%d>\n", __FUNCTION__, __LINE__);
            fflush(fp);

            realCounter    = gridTask.realTaskGet(myTask, realTask);
            replaceCounter = gridTask.replaceTaskGet(myTask, replaceTask);

            fprintf(fp, "<%s><%d>\n", __FUNCTION__, __LINE__);
            fflush(fp);

            // Swap grids (after repair) loop
            for (int i = 0; i < realCounter; i++)
            {
                std::swap(realTask[i]->local_newgrid,
                          realTask[i]->local_grid);
            }

            fprintf(fp, "<%s><%d>\n", __FUNCTION__, __LINE__);
            fflush(fp);

            if (rank == 5)
            {
                //gridTask.show();
            }

            goto restart_step;
        }

        if (maxdiff < EPS)
        {
            break;
        }

        thalo -= MPI_Wtime();

        if (replaceCounter > 0)
        {
            fprintf(fp, "<%s><%d> Starting replaceExchange (%d)(%d)\n",
                __FUNCTION__, __LINE__, realCounter, replaceCounter);
            fflush(fp);

            int sizeRR = (8 * realCounter) * replaceCounter;

            int exchange = 0;
            MPI_Request reqs[sizeRR];
            for (int i = 0; i < replaceCounter; i++)
            {
                task_t *replaces[commsize];
                memset(replaces, 0, sizeof(task_t *) * commsize);

                gridTask.realTaskGet(replaceTask[i], replaces);

                for (int j = 0; j < realCounter; j++)
                {
                    halo_cookie_t halo_cookie;

                    halo_cookie.task      = replaces[j];
                    halo_cookie.neighbors = replaceTask[i];
                    halo_cookie.count     = 1;
                    halo_cookie.row       = row;
                    halo_cookie.col       = col;
                    halo_cookie.comm      = comm;
                    halo_cookie.reqs      = reqs;
                    halo_cookie.nx        = nx;
                    halo_cookie.ny        = ny;

                    /*
                     * Step 5: Halo exchange:
                     * T = 4 * (a + b * (rows / py)) + 4 * (a + b * (cols / px))
                     */
                    exchange = haloExchange(&halo_cookie, gridTask, rank, exchange, j, 'R');
                }
            }

            fprintf(fp, "<%s><%d> MPI_Waitall replaceExchange\n", __FUNCTION__, __LINE__);
            fflush(fp);
            rc = MPI_Waitall(sizeRR, reqs, MPI_STATUS_IGNORE);
            fprintf(fp, "<%s><%d> Ending replaceExchange\n", __FUNCTION__, __LINE__);
            fflush(fp);
        }

        fprintf(fp, "<%s><%d> Starting haloExchange (%d)\n",
            __FUNCTION__, __LINE__, realCounter);
        fflush(fp);

        int exchange = 0;
        MPI_Request reqs[8 * realCounter];
        memset(&reqs, 0, sizeof(MPI_Request) * 8 * realCounter);
        for (int i = 0; i < realCounter; i++)
        {
            halo_cookie_t halo_cookie;

            halo_cookie.task      = realTask[i];
            halo_cookie.neighbors = realTask[0];
            halo_cookie.count     = 1;
            halo_cookie.row       = row;
            halo_cookie.col       = col;
            halo_cookie.comm      = comm;
            halo_cookie.reqs      = reqs;
            halo_cookie.nx        = nx;
            halo_cookie.ny        = ny;

            /*
             * Step 5: Halo exchange:
             * T = 4 * (a + b * (rows / py)) + 4 * (a + b * (cols / px))
             */
            exchange = haloExchange(&halo_cookie, gridTask, rank, exchange, i, ' ');
        }

        fprintf(fp, "<%s><%d>MPI_Waitall haloExchange\n", __FUNCTION__, __LINE__);
        fflush(fp);

        // Step 6: Wait all
        rc = MPI_Waitall(8 * realCounter, reqs, MPI_STATUS_IGNORE);
        thalo += MPI_Wtime();

        fprintf(fp, "<%s><%d>Ending haloExchange\n", __FUNCTION__, __LINE__);
        fflush(fp);

        /*
         * Collective operation may fail so
         * ULFM provide repair operation
         */
        if (MPI_ERR_PROC_FAILED == rc)
        {
            fprintf(fp, "<%s><%d>Not handled error\n", __FUNCTION__, __LINE__);
            fflush(fp);
            goto exit; // TODO
        }
    }

exit:

    fprintf(fp, "<%s><%d> Ending calculation\n", __FUNCTION__, __LINE__);
    fflush(fp);

    MPI_Type_free(&row);
    MPI_Type_free(&col);

    ttotal += MPI_Wtime();

    if (rank == 0)
    {
        printf("# Heat 2D (mpi): grid: rows %d, cols %d, procs %d (px %d, py %d)\n",
               rows, cols, commsize, px, py);
    }
    else
    {
        sleep(1); // Small delay to print
    }

    int namelen;
    char procname[MPI_MAX_PROCESSOR_NAME];

    MPI_Get_processor_name(procname, &namelen);

    printf("# P %4d (%2d, %2d) on %s: grid ny %d nx %d, total %.6f,"
           " mpi %.6f (%.2f) = allred %.6f (%.2f) + halo %.6f (%.2f)\n", 
           rank, rankx, ranky, procname, ny, nx, ttotal, treduce + thalo, 
           (treduce + thalo) / ttotal, treduce, treduce / (treduce + thalo), 
           thalo, thalo / (treduce + thalo)); 

    double prof[3] = { ttotal, treduce, thalo };

    if (rank == 0)
    {
        MPI_Reduce(MPI_IN_PLACE, prof, NELEMS(prof), MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

        sleep(1); // Small delay to print last

        printf("# procs %d : grid %d %d : niters %d : total time %.6f :"
               " mpi time %.6f : allred %.6f : halo %.6f\n", 
               commsize, rows, cols, niters, prof[0],
               prof[1] + prof[2], prof[1], prof[2]);
    }
    else
    {
        MPI_Reduce(prof, NULL, NELEMS(prof), MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    }

    fclose(fp);

    MPI_Finalize();
    return 0;
}
