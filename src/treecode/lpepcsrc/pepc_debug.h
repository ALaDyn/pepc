
#define DEBUG_HEADER(file) write(*,*) 'HDR '

#define DEBUG_INFO(format, ...) write(*,*) 'IFO ', __FILE__ ,__LINE__

#define DEBUG_WARNING_ALL(format, ...) write(*,*) 'WRN ', __FILE__ ,__LINE__

#define DEBUG_WARNING(format, ...) write(*,*) 'WRN ', __FILE__ ,__LINE__

#define DEBUG_ERROR(format, ...) write(*,*) 'ERR ', __FILE__ ,__LINE__

#define DEBUG_ERROR_NO_HEADER(format, ...) write(*,*) 'ERR '

#define DEBUG_ERROR_NO_DIAGFILE(format, ...) write(*,*) 'ERR ', __FILE__ ,__LINE__

#define DEBUG_DATA(format, ...) write(*,*) 'DTA ', __FILE__ ,__LINE__

#define DEBUG_STRINGIFY_HELPER(s) #s

#define DEBUG_STRINGIFY(s) DEBUG_STRINGIFY_HELPER(s)

#define ERROR_ON_FAIL_HELPER(ret, sep, msg) if (0/=ret) write(*,*) 'ERR ', __FILE__ ,__LINE__

#define ERROR_ON_FAIL(ret) ERROR_ON_FAIL_HELPER(ret, "", "")

#define ERROR_ON_FAIL_MSG(ret, msg) ERROR_ON_FAIL_HELPER(ret, ": ", msg)

#ifdef NDEBUG

#define DEBUG_ASSERT(cond) ! Assertion: cond

#define DEBUG_ASSERT_MSG(cond, fmt, ...) ! Assertion: cond

#else

#define DEBUG_ASSERT_MSG(cond, fmt, ...) if(.not.(cond)) write(*,*) 'ASS1 ', __FILE__ ,__LINE__

#define DEBUG_ASSERT(cond) if(.not.(cond)) write(*,*) 'ASS2 ', __FILE__ ,__LINE__

#endif

