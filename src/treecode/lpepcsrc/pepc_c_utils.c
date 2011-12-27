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

void create_directory_c(char dirname[])
{
  char cwd[PATH_MAX];
  char fullpath[PATH_MAX];

  printf("%s\n", dirname);

  getcwd(cwd, PATH_MAX);
  sprintf(fullpath, "%s/%s", cwd, dirname);

  if (0 != mkdir(fullpath, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH))
  {
    // ignore the error if the directory already existed
    if (EEXIST != errno) printf("Error while creating directory %s: %d\n", fullpath, errno);
  }
}


