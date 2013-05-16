
#define DEBUG_HEADER(file) write(*,*) 'DBG'

#define DEBUG_INFO(format, ...) write(*,*) 'DBG'

#define DEBUG_WARNING_ALL(format, ...) write(*,*) 'DBG'

#define DEBUG_WARNING(format, ...) write(*,*) 'DBG'

#define DEBUG_ERROR(format, ...) write(*,*) 'DBG'

#define DEBUG_ERROR_NO_HEADER(format, ...) write(*,*) 'DBG'

#define DEBUG_ERROR_NO_DIAGFILE(format, ...) write(*,*) 'DBG'

#define DEBUG_DATA(format, ...) write(*,*) 'DBG'

#define DEBUG_STRINGIFY_HELPER(s) #s

#define DEBUG_STRINGIFY(s) DEBUG_STRINGIFY_HELPER(s)

#define DEBUG_ERROR_ON_FAIL_HELPER(ret, sep, msg) if (0/=ret) write(*,*) 'DBG'

#define DEBUG_ERROR_ON_FAIL(ret) DEBUG_ERROR_ON_FAIL_HELPER(ret, "", "")

#define DEBUG_ERROR_ON_FAIL_MSG(ret, msg) DEBUG_ERROR_ON_FAIL_HELPER(ret, ": ", msg)

#ifdef NDEBUG

#define DEBUG_ASSERT(cond) ! Assertion: cond

#define DEBUG_ASSERT_MSG(cond, fmt, ...) ! Assertion: cond

#else

#define DEBUG_ASSERT_MSG(cond, fmt, ...) write(*,*) 'DBG'

#define DEBUG_ASSERT(cond) write(*,*) 'DBG'

#endif

