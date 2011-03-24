/****************************************************************************
 *
 *               dummy implementation for libhpc
 *
 *               C. Pospiech ** ACTC
 *
 *               September 2010
 *
 ****************************************************************************/

#include <libhpc.h>
#include <stdio.h>
#include <string.h>

#ifdef __cplusplus
extern "C" {
#endif

  /*
   * This count is incremented by the error handler
   * (after getting the above mutex lock)
   */
  int hpm_error_count = 0;

  void _hpmInit_( int my_ID, const char* progName, int sz_progName ) {

    char out_text[300];
    const char *out_text_begin = "Called hpmInit with progName = ";
    int out_text_len = strlen(out_text_begin);

    strcpy(out_text, out_text_begin);
    strncpy(out_text + out_text_len, progName, sz_progName);
    out_text_len += sz_progName;
    out_text[out_text_len] = '\0';
    printf("%s with ID %d\n", out_text, my_ID);
  }

  void _hpm_start_( int inst_ID, int parent_ID, int line, const char* file,
		    const char* label, int sz_file, int sz_label,
		    int count_context) {

    char out_text[300];
    const char *out_text_begin = "Called hpmStart with label = ";
    int out_text_len = strlen(out_text_begin);

    strcpy(out_text, out_text_begin);
    strncpy(out_text + out_text_len, label, sz_label);
    out_text_len += sz_label;
    out_text[out_text_len] = '\0';
    printf("%s with ID %d\n", out_text, inst_ID);
  }

  void _hpm_stop_( int inst_ID, int line, int count_context) {

    printf("Called hpmStop with ID %d\n", inst_ID);
  }


  void _hpm_terminate_( int my_ID ) {

    printf("Called hpmTerminate with ID %d\n", my_ID);
  }

#ifdef __cplusplus
}
#endif

