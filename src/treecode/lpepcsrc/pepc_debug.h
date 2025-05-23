! CPP + Fortran do not allow for multiline replacements. We use ';' to combine several statements
! in a single line. To keep the line length to a minimum, we use ';;' which will be replaced by
! '\n' after preprocessing.

#define DEBUG_HEADER(file) \
    write(file,'("[PE ", I6.6, ", file: ",a, ", line ", I0, "] ")', advance='no') debug_my_rank, __FILE__, __LINE__


#define DEBUG_INFO(format, ...) \
        DEBUG_HEADER(debug_ipefile);; \
        write(debug_ipefile, format) __VA_ARGS__;; \
        flush(debug_ipefile);;

#define DEBUG_WARNING_ALL(format, ...) \
        DEBUG_HEADER(debug_ipefile);; \
        write(debug_ipefile, format) __VA_ARGS__;; \
        flush(debug_ipefile);; \
        DEBUG_HEADER(debug_stdout);; \
        write(debug_stdout, format) __VA_ARGS__;;

#define DEBUG_WARNING(format, ...) \
        DEBUG_HEADER(debug_ipefile);; \
        write(debug_ipefile, format) __VA_ARGS__;; \
        flush(debug_ipefile);; \
        if (debug_my_rank == 0) DEBUG_HEADER(debug_stdout);; \
        if (debug_my_rank == 0) write(debug_stdout, format) __VA_ARGS__;;

#define DEBUG_ERROR(format, ...) \
        DEBUG_HEADER(debug_ipefile);; \
        write(debug_ipefile, format) __VA_ARGS__;; \
        flush(debug_ipefile);; \
        DEBUG_HEADER(debug_stdout);; \
        write(debug_stdout, format) __VA_ARGS__;; \
        call debug_mpi_abort();;

#define DEBUG_ERROR_NO_HEADER(format, ...) \
        write(debug_ipefile, format) __VA_ARGS__;; \
        flush(debug_ipefile);; \
        write(debug_stdout, format) __VA_ARGS__;; \
        call debug_mpi_abort();;

#define DEBUG_ERROR_NO_DIAGFILE(format, ...) \
        write(debug_stdout, format) __VA_ARGS__;; \
        call debug_mpi_abort();;

#define DEBUG_DATA(format, ...) \
        write(debug_ipefile, format) __VA_ARGS__;; \
        flush(debug_ipefile);;

#define DEBUG_STRINGIFY_HELPER(s) #s

#define DEBUG_STRINGIFY(s) DEBUG_STRINGIFY_HELPER(s)

#define ERROR_ON_FAIL_HELPER(ret, sep, msg) \
        if (0 /= ret) then;; \
          DEBUG_ERROR('("ERROR_ON_FAIL: ",a," = ",I0,a,a)', DEBUG_STRINGIFY(ret), ret, sep, msg);; \
        end if;;

#define ERROR_ON_FAIL(ret) ERROR_ON_FAIL_HELPER(ret, "", "")

#define ERROR_ON_FAIL_MSG(ret, msg) ERROR_ON_FAIL_HELPER(ret, ": ", msg)

#ifdef NDEBUG

#define DEBUG_ASSERT(cond) ! Assertion: cond

#define DEBUG_ASSERT_MSG(cond, fmt, ...) ! Assertion: cond

#define DEBUG_PURE pure

#define DEBUG_ELEMENTAL elemental

#else

#define DEBUG_ASSERT_MSG(cond, fmt, ...) \
        if (.not. (cond)) then;; \
          DEBUG_HEADER(debug_ipefile);; \
          write(debug_ipefile, '("Assertion failed: ", a, " ")', advance='no') DEBUG_STRINGIFY(cond);; \
          write(debug_ipefile, fmt) __VA_ARGS__;; \
          flush(debug_ipefile);; \
          DEBUG_HEADER(debug_stdout);; \
          write(debug_stdout, '("Assertion failed: ", a, " ")', advance='no') DEBUG_STRINGIFY(cond);; \
          write(debug_stdout, fmt) __VA_ARGS__;; \
          call debug_mpi_abort();; \
        end if;;

#define DEBUG_ASSERT(cond) \
        if (.not. (cond)) then;; \
          DEBUG_ERROR('("Assertion failed: ", a)', DEBUG_STRINGIFY(cond));; \
        end if;;

#define DEBUG_PURE
#define DEBUG_ELEMENTAL

#endif

#ifdef MAC_NEEDS_PARTICLE
  #define IF_MAC_NEEDS_PARTICLE(param) param,
#else
  #define IF_MAC_NEEDS_PARTICLE(param)
#endif
