#ifndef __SL_TUNE_H__
#define __SL_TUNE_H__

/* do reduce+bcast instead of allreduce on jugene */
#ifdef JUGENE
# define MPI_PARTITION_RADIX_REDUCEBCAST_THRESHOLD  0
#endif
#endif /* __SL_TUNE_H__ */
