#ifndef __SL_TUNE_H__
#define __SL_TUNE_H__

#define MSEG_BORDER_UPDATE_REDUCTION

#define MSEG_INFO

/*#define MSS_ROOT*/


/* do reduce+bcast instead of allreduce on jugene */
#ifdef JUGENE
# define GLOBAL_REDUCEBCAST_THRESHOLD  0
#endif
#endif /* __SL_TUNE_H__ */
