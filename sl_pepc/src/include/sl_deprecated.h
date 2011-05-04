/*
 *  SL - Sorting Library, v0.1, (michael.hofmann@informatik.tu-chemnitz.de)
 *  
 *  file: src/include/sl_deprecated.h
 *  timestamp: 2011-03-06 21:59:31 +0100
 *  
 */


#ifndef __SL_DEPRECATED_H__
#define __SL_DEPRECATED_H__


/* sl_macro MPI_BINNING_REDUCEBCAST_THRESHOLD */
#ifdef MB_REDUCEBCAST_THRESHOLD
# define MPI_BINNING_REDUCEBCAST_THRESHOLD  MB_REDUCEBCAST_THRESHOLD
#endif


#endif /* __SL_DEPRECATED_H__ */
