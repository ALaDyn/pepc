/*************************************************************************
>
>  different utils for directly accessing POSIX functions
>
*************************************************************************/

#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <sys/stat.h>
#include <limits.h>
#include <stdio.h>
#include <errno.h>

#ifndef PATH_MAX 
  #define MYPATH_MAX 255 
#else 
  #define MYPATH_MAX PATH_MAX 
#endif 


void create_directory_c(char dirname[])
{
  char cwd[MYPATH_MAX];
  char fullpath[MYPATH_MAX];

  getcwd(cwd, MYPATH_MAX);
  sprintf(fullpath, "%s/%s", cwd, dirname);

  if (0 != mkdir(fullpath, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH))
  {
    // ignore the error if the directory already existed
    if (EEXIST != errno) printf("Error while creating directory %s: %d\n", fullpath, errno);
  }
}


