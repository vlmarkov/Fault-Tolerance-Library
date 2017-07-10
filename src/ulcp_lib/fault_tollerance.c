#include "utils.h"
#include "fault_tollerance.h"

#include <stdio.h>
#include <stdlib.h>

#ifdef ULFM_SUPPORT
void ulcp_verbose_errhandler(MPI_Comm *comm, int *err, ...)
{
    int len;
    char errstr[MPI_MAX_ERROR_STRING];

    MPI_Error_string(*err, errstr, &len);

    fprintf(stderr, "[ULCP] Rank %d/%d: Notified of error %s\n", 
        ulcp_get_comm_rank(), ulcp_get_comm_size(), errstr);

    exit(1);
}
#endif /* ULFM_SUPPORT */
