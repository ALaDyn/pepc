
#define DEBUG_HEADER(file)  \
    write(file,'("[PE ", I6.6, ", file: ",a, ", line ", I0, "] ")', advance='no') debug_my_rank, __FILE__, __LINE__


#define DEBUG_INFO(format, ...)                    \
        call debug_ipefile_open();                  \
          DEBUG_HEADER(debug_ipefile);              \
          write(debug_ipefile, format) __VA_ARGS__; \
        call debug_ipefile_close();

#define DEBUG_WARNING_ALL(format, ...)             \
        call debug_ipefile_open();                  \
          DEBUG_HEADER(debug_ipefile);              \
          write(debug_ipefile, format) __VA_ARGS__; \
        call debug_ipefile_close();                 \
        DEBUG_HEADER(debug_stdout);                 \
        write(debug_stdout, format) __VA_ARGS__;

#define DEBUG_WARNING(format, ...)                 \
        call debug_ipefile_open();                  \
          DEBUG_HEADER(debug_ipefile);              \
          write(debug_ipefile, format) __VA_ARGS__; \
        call debug_ipefile_close();                 \
        if (debug_my_rank == 0) DEBUG_HEADER(debug_stdout);              \
        if (debug_my_rank == 0) write(debug_stdout, format) __VA_ARGS__;

#define DEBUG_ERROR(format, ...)                   \
        call debug_ipefile_open();                  \
          DEBUG_HEADER(debug_ipefile);              \
          write(debug_ipefile, format) __VA_ARGS__; \
        call debug_ipefile_close();                 \
        DEBUG_HEADER(debug_stdout);                 \
        write(debug_stdout, format) __VA_ARGS__;    \
        call debug_mpi_abort();

#define DEBUG_ERROR_NO_HEADER(format, ...)         \
        call debug_ipefile_open();                  \
          write(debug_ipefile, format) __VA_ARGS__; \
        call debug_ipefile_close();                 \
        write(debug_stdout, format) __VA_ARGS__;    \
        call debug_mpi_abort();

#define DEBUG_ERROR_NO_DIAGFILE(format, ...)         \
        write(debug_stdout, format) __VA_ARGS__;    \
        call debug_mpi_abort();

#define DEBUG_DATA(format, ...)                    \
        call debug_ipefile_open();                  \
          write(debug_ipefile, format) __VA_ARGS__; \
        call debug_ipefile_close();

#define DEBUG_STRINGIFY_HELPER(s) #s

#define DEBUG_STRINGIFY(s) DEBUG_STRINGIFY_HELPER(s)

#ifdef NDEBUG
#define DEBUG_ASSERT(cond) ! DEBUG_ASSERT(cond)
#define DEBUG_ASSERT_MSG(cond, msg) ! DEBUG_ASSERT_MSG(cond, msg)
#else
#define DEBUG_ASSERT(cond) \
        if (.not. (cond)) then; \
          DEBUG_ERROR('("Assertion failed: ", a)', DEBUG_STRINGIFY(cond)); \
        end if;

#define DEBUG_ASSERT_MSG(cond, msg) \
        if (.not. (cond)) then; \
          DEBUG_ERROR('("Assertion failed: ", a, ": ", a)', DEBUG_STRINGIFY(cond), msg); \
        end if;
#endif
