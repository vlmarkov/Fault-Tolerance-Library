#ifndef _ULCP_FT_H_
#define _ULCP_FT_H_

#include <mpi.h>

#ifdef ULFM_SUPPORT
void ulcp_verbose_errhandler(MPI_Comm *comm, int *err, ...);
#endif /* ULFM_SUPPORT */

#endif /* _ULCP_FT_H_ */
