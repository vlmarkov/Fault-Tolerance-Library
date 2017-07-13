/**
 ** nbody-mpi.c
 **
 ** Calculates the solution to a 3D n-body problem in parallel using MPI.
 **
 ** The number of particles to simulate and the number of iterations, are given
 ** as constants at the top of this file. Each particle is assigned a mass of
 ** 1.0 and a random position by the "master" processor. The master processor
 ** then evenly divides particles amongst the available processors. During each
 ** iteration, processors calculate the new positions of particles assigned to
 ** them. Final particle positions are relayed back to the master processor once
 ** all iterations have been completed.
 **
 ** Forces are calculated within each iteration with a ring pipeline where the
 ** number of pipeline stages is equal to the number of available processors.
 ** During a given stage, each processor is calculating the force exerted by one
 ** set of particles upon the particles on that processor. Consider an example
 ** with 2500 particles and 4 processors. Each processor (0, 1, 2, 3) is given
 ** a set (A, B, C, D) of 625 particles. The force calculations done during the
 ** pipeline stages would look like:
 **
 **           +-------------------------------------------------------------+
 **           |                                                             |
 **           +- Processor 0 <- Processor 1 <- Processor 2 <- Processor 3 <-+
 **
 **  Stage 0      A on A          B on B         C on C         D on D
 **  Stage 1      B on A          C on B         D on C         A on D
 **  Stage 2      C on A          D on B         A on C         B on D
 **  Stage 3      D on A          A on B         B on C         C on D
 **
 ** A processor receives one message from the processor on its "right" and sends
 ** one message to the processor on its "left" during each stage. All processors
 ** know the total force acting on each of their particles after the pipeline
 ** stages are completed. New particle positions can then be calculated on each
 ** processor independently and in parallel. Variable time steps are used to
 ** keep things from getting out of control. The new time step is synchronized
 ** across all processors before beginning the next iteration.
 **
 ** This application was designed for testing MPI performance tools. Thus it was
 ** written to contain both collective and point-to-point operations. And in the
 ** grand tradition of physics simulation codes everywhere, it is written as a
 ** single, giant, main() function. It does NO error checking and produces NO
 ** output. It just calculates!
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <mpi.h>

#ifdef ULFM_SUPPORT
#include <mpi-ext.h> // ULFM support
#endif /* ULFM_SUPPORT */

#include "../../src/ulcp_lib/ulcp.h"


static const int MasterProcessor = 0;     // Number of the master processor

static const int NumParticles   = 8192;   // Number of particles
static const int NumIterations  = 100;    // Number of iterations
static const double Infinity    = 1.0E96; // Infinitely large value
static const double MinTimeStep = 1.0E-6; // Minimum allowed time step


typedef struct {
    double mass;
    double x, y, z;
} particle_t;


double       ttotal;   // Time measure
int          iteration;// ...
particle_t * local;    // Array containing our local particles
int        * num;      // Number of particles on each processor

/*****************************************************************************/
/* User defiend callbacks                                                    */
/* Additional functions helps save and get data                              */
/*****************************************************************************/
inline static void user_save_callback(int phase)
{
    MPI_File local_snapshot;

    ulcp_open_file(&local_snapshot, phase);

    int delta_idx = 1;
    int particles = NumParticles;
    int iterations= NumIterations;

    ulcp_snapshot_save_compressed(local_snapshot, &particles, 1, MPI_INT, delta_idx++);
    ulcp_snapshot_save_compressed(local_snapshot, &iterations, 1, MPI_INT, delta_idx++);
    ulcp_snapshot_save_compressed(local_snapshot, &ttotal, 1, MPI_DOUBLE, delta_idx++);
    ulcp_snapshot_save_compressed(local_snapshot, &iteration, 1, MPI_INT, delta_idx++);

    ulcp_snapshot_save_compressed_complex(local_snapshot, local,
                                          num[0], sizeof (particle_t),
                                          delta_idx++);


    ulcp_close_file(&local_snapshot);
}



inline static int checkpoint_get(particle_t * grid,
                                 int          size,
                                 int          rank,
                                 double     * ttotal,
                                 int        * niters)
{
    int readNumParticles, readNumIterations;

    char rank_str[256]             = { 0 };
    char last_checkpoint_path[256] = { 0 };
    char last_checkpoint[256]      = { 0 };

    int phase = ulcp_get_snapshot(last_checkpoint);

    sprintf(rank_str, "/%d/", rank);

    strcpy(last_checkpoint_path, ULCP_SNAPSHOT_DIR_NAME);
    strcat(last_checkpoint_path, rank_str);
    strcat(last_checkpoint_path, last_checkpoint);

    FILE *file = ulcp_open_snapshot(last_checkpoint_path, "rb");
    if (!file)
    {
        fprintf(stderr, "[ULCP] Can't read %s\n", last_checkpoint_path);
        exit(1);
    }

    // copy the file into the buffer:
    fread(&readNumParticles, sizeof(int), 1, file);
    fread(&readNumIterations, sizeof(int), 1, file);

    if ((readNumParticles  != NumParticles) &&
        (readNumIterations != NumIterations))
    {
        fprintf(stdout, "[ULCP] Snapshot size not match\n");
        fclose(file);
        exit(1);
    } 
    else
    {
        fprintf(stdout, "[ULCP] Snapshot size match\n");
    }

    fread(ttotal, sizeof(double), 1,    file);
    fread(niters, sizeof(int),    1,    file);
    fread(grid,   sizeof(double), size, file);

    fclose(file);
    return phase;
}

/*****************************************************************************/
/* Main function                                                             */
/*****************************************************************************/
int main(int argc, char* argv[])
{
    MPI_Datatype type;      // MPI data type for communicating particle data
    int num_processors;     // Number of processors being used
    int processor;          // My processor number
    int* offset;            // Offset to start of each processor's particles
    int buffer_size;        // Number of particles in pipeline data buffers

    // Allways FIRST
    MPI_Init(&argc, &argv);


    /*************************************************************************/
    /* Initialize checkpoint library                                         */
    /*************************************************************************/
    ulcp_init(2, argc, argv);

    ulcp_init_checkpoit(&&phase_one);
    ulcp_init_checkpoit(&&phase_two);


    // Create the MPI data type for communicating particle data
    MPI_Type_contiguous(4, MPI_DOUBLE, &type);
    MPI_Type_commit(&type);   

    // Determine the number of procesors being used and our processor number
    MPI_Comm_size(MPI_COMM_WORLD, &num_processors);
    MPI_Comm_rank(MPI_COMM_WORLD, &processor);

    // Time measure
    ttotal = MPI_Wtime();

    // Determine how the particles are allocated to the processors
    {
        num = (int*)malloc(num_processors * sizeof(int));
        offset = (int*)malloc(num_processors * sizeof(int));

        for (int p = 0; p < num_processors; p++)
        {

            num[p] = (NumParticles / num_processors) + ((p < (NumParticles % num_processors)) ? 1 : 0);
        }
        buffer_size = num[0];

        offset[0] = 0;
        for (int p = 1; p < num_processors; p++)
        {
            offset[p] = offset[p - 1] + num[p - 1];
        }
    }

    // Distribute the initial particle state
    {
        particle_t* particles = NULL;

        if (processor == MasterProcessor)
        {
            particles = (particle_t*)malloc(NumParticles * sizeof(particle_t));
            // Randomize the particles
            for (int i = 0; i < NumParticles; i++)
            {
                particles[i].mass = 1.0;
                particles[i].x = drand48();
                particles[i].y = drand48();
                particles[i].z = drand48();
            }
        }

        local = (particle_t*)malloc(num[processor] * sizeof(particle_t));

        MPI_Scatterv(particles,
            num,
            offset,
            type,
            local,
            num[processor],
            type,
            MasterProcessor,
            MPI_COMM_WORLD);
        
        if (processor == MasterProcessor)
        {
            free(particles);
        }
    }

    // Actual Simulation
    {
        if (ulcp_is_recovery_mode())
        {
            int checkpoint = checkpoint_get(local,
                                            num[processor],
                                            processor,
                                            &ttotal,
                                            &iteration);

            // Jumping to checkpoint
            ulpc_goto_checkpoint(checkpoint);
        }

        //ulcp_action_t action = {
        //    .mode = ULCP_SET_MODE_COMPLEX,
        //    .user_type = sizeof(particle_t),
        //    .size = num[processor]
        //};

        //ulcp_snapshot_set_diff(&action);
    
        ulcp_save_data(&&phase_one, user_save_callback);
        ulpc_set_checkpoint(phase_one);


        particle_t* buf_send = (particle_t*)malloc(buffer_size * sizeof(particle_t));
        particle_t* buf_recv = (particle_t*)malloc(buffer_size * sizeof(particle_t));

        double* tfx = (double*)malloc(num[processor] * sizeof(double));
        double* tfy = (double*)malloc(num[processor] * sizeof(double));
        double* tfz = (double*)malloc(num[processor] * sizeof(double));
        double* ox = (double*)malloc(num[processor] * sizeof(double));
        double* oy = (double*)malloc(num[processor] * sizeof(double));
        double* oz = (double*)malloc(num[processor] * sizeof(double));

        // Set the "old" position for each particle to its current position
        for (int i = 0; i < num[processor]; i++)
        {
            ox[i] = local[i].x;
            oy[i] = local[i].y;
            oz[i] = local[i].z;
        }

        // Time steps
        for (/*int iteration = 0*/; iteration < NumIterations; iteration++)
        {
            double f_max = -Infinity;

            // Show current iteration number
            if (processor == 1)
            {
                #ifdef DEBUG_ITER
                fprintf(stdout, "Iteration %d of %d...\n", iteration + 1, NumIterations);
                fflush(stdout);
                #endif // DEBUG_ITER
            }

            // Zero the total force for each particle
            for (int i = 0; i < num[processor]; i++)
            {
                tfx[i] = 0.0;
                tfy[i] = 0.0;
                tfz[i] = 0.0;
            }

            // Force computation pipeline
            for (int stage = 0; stage < num_processors; stage++)
            {
                MPI_Request request[2];
                MPI_Status status[2];

                // Prime the pipeline with our local data for stage zero
                if (stage == 0)
                {
                    memcpy(buf_send, local, num[processor] * sizeof(particle_t));
                }

                // Issue the send/receive pair for this pipeline stage
                if (stage < (num_processors - 1))
                {
                    MPI_Isend(buf_send, buffer_size, type, 
                        (processor - 1 + num_processors) % num_processors,
                        0, MPI_COMM_WORLD, &request[0]);

                    MPI_Irecv(buf_recv, buffer_size, type,
                        (processor + 1 + num_processors) % num_processors,
                        0, MPI_COMM_WORLD, &request[1]);
                }

                // Compute forces
                for (int i = 0; i < num[processor]; i++)
                {
                    double r_min = +Infinity;
                    double fx = 0.0, fy = 0.0, fz = 0.0, f = 0.0;

                    for (int j = 0; j < num[(processor + stage) % num_processors]; j++)
                    {
                        double rx = local[i].x - buf_send[j].x;
                        double ry = local[i].y - buf_send[j].y;
                        double rz = local[i].z - buf_send[j].z;
                        double r = (rx * rx) + (ry * ry) + (rz * rz);

                        if (r > 0.0)
                        {
                            if (r < r_min)
                                r_min = r;

                            fx -= buf_send[j].mass * (rx / r);
                            fy -= buf_send[j].mass * (ry / r);
                            fz -= buf_send[j].mass * (rz / r);
                        }
                    }

                    tfx[i] += fx;
                    tfy[i] += fy;
                    tfz[i] += fz;

                    // Rough estimate of 1/m|df/dx|
                    f = sqrt((fx * fx) + (fy * fy) + (fz * fz)) / r_min;
                    if (f > f_max)
                    {
                        f_max = f;
                    }
                }

                

                // Complete the send/receive pair for this pipeline stage
                if (stage < (num_processors - 1))
                {
                    MPI_Waitall(2, request, status);
                    memcpy(buf_send, buf_recv, buffer_size * sizeof(particle_t));
                }
            }

            // Compute new positions using a simple leapfrog time integration.
            // Use a variable step version to simplify time-step control.
            // Integration is (a0 * x^+) + (a1 * x) + (a2 * x^-) = f / m
            // Stability criteria is roughly 2.0 / sqrt(1/m|df/dx|) >= dt

            {
                static double dt_old = 0.001, dt_now = 0.001;
                double dt_est, dt_new;
                double a0 = +2.0 / (dt_now * (dt_old + dt_now));
                double a1 = -2.0 / (dt_old * dt_now);
                double a2 = +2.0 / (dt_old * (dt_old + dt_now));

                for (int i = 0; i < num[processor]; i++)
                {
                    double x = local[i].x;
                    double y = local[i].y;
                    double z = local[i].z;

                    local[i].x = (tfx[i] - (a1 * x) - (a2 * ox[i])) / a0;
                    local[i].y = (tfy[i] - (a1 * y) - (a2 * oy[i])) / a0;
                    local[i].z = (tfz[i] - (a1 * z) - (a2 * oz[i])) / a0;

                    ox[i] = x; oy[i] = y; oz[i] = z;
                }

                dt_est = 1.0 / sqrt(f_max);
                if (dt_est < MinTimeStep)
                {
                    dt_est = MinTimeStep;
                }

                MPI_Allreduce(&dt_est, &dt_new, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

                if (dt_new < dt_now)
                {
                    dt_old = dt_now;
                    dt_now = dt_new;
                }
                else if(dt_new > (4.0 * dt_now))
                {
                    dt_old = dt_now;
                    dt_now *= 2.0;
                }
            }

            #ifdef ULFM_TEST
            if(ulcp_get_comm_rank() == (ulcp_get_comm_size() - 1))
            {
                fprintf(stdout, "[ULCP] Rank %d: committing suicide\n", ulcp_get_comm_rank());
                raise(SIGKILL);
            }
            #endif // ULFM_TEST
            /*
             * 1. Manualy kill some process
             * 2. Get error by some collective operation, e.g. MPI_Barrier()
             * 3. Will be invoke ulcp_verbose_errhandler()
             */
            // TODO SAVE
            if (iteration % 50 == 0)
            {
                ulcp_save_data(&&phase_one, user_save_callback);
            }
        }

        free(buf_send);
        free(buf_recv);
        free(tfx);
        free(tfy);
        free(tfz);
        free(ox);
        free(oy);
        free(oz);
    }

    ulcp_save_data(&&phase_two, user_save_callback);
    ulpc_set_checkpoint(phase_two);

    // Gather the final particle state
    {
        particle_t* particles = NULL;

        if (processor == MasterProcessor)
        {
            particles = (particle_t*)malloc(NumParticles * sizeof(particle_t));
        }

        MPI_Gatherv(local, num[processor], type,
            particles, num, offset, type,
            MasterProcessor, MPI_COMM_WORLD);
    
        free(local);

        if (processor == MasterProcessor)
        {
            free(particles);
        }
    }

    ttotal = MPI_Wtime() - ttotal;
    fprintf(stdout, "Rank: %d, elapsed time %f sec\n", processor, ttotal);
    fprintf(stdout, "Particles per processor %d, Iterations %d",
        NumParticles/num_processors, NumIterations);

    // Free the particle distribution arrays
    free(num);
    free(offset);

    // Finalize MPI
    ulcp_finalize();
    MPI_Finalize();
    return 0;
}
