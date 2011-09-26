
/*************************************************************************
>
>  pthreads interface library
>
*************************************************************************/

#include <pthread.h>
#include <stdlib.h>
#include <stdio.h>
#include <sched.h>
#include "fortran2c_types.h"
#include <unistd.h>
#include <sys/syscall.h>
#include <sys/types.h>

pthread_t *__restrict__ my_threads;
pthread_rwlock_t *my_rwlocks;
pthread_attr_t thread_attr;
int maxnumthreads = 0;
int maxnumlocks = 0;

#define CHECKRES do {if (iret != 0) return iret;} while(0);


//////////////// PThreads //////////////////////

FINT_TYPE_C pthreads_init(FINT_TYPE_C numthreads)
{
    return pthreads_init_(numthreads);
}
FINT_TYPE_C pthreads_init_(FINT_TYPE_C numthreads)
{
    int iret = 0;

    maxnumthreads = numthreads;
    my_threads    = (pthread_t*)malloc(((unsigned int)maxnumthreads)*sizeof(pthread_t));

    iret = pthread_attr_init(&thread_attr);
    CHECKRES;

    iret = pthread_attr_setdetachstate(&thread_attr, PTHREAD_CREATE_JOINABLE);
    CHECKRES;

    iret = pthread_attr_setscope(&thread_attr, PTHREAD_SCOPE_SYSTEM);
    CHECKRES;

    return 0;
}


FINT_TYPE_C pthreads_uninit()
{
    return pthreads_uninit_();
}
FINT_TYPE_C pthreads_uninit_()
{
    int iret = 0;
    free(my_threads);

    iret = pthread_attr_destroy(&thread_attr);
    CHECKRES;

    return 0;
}


FINT_TYPE_C pthreads_createthread(FINT_TYPE_C id, void *(*start_routine) (void *), void *arg, FINT_TYPE_C relative_priority)
{
    return pthreads_createthread_(id, start_routine, arg, relative_priority);
}
FINT_TYPE_C pthreads_createthread_(FINT_TYPE_C id, void *(*start_routine) (void *), void *arg, FINT_TYPE_C relative_priority)
{
    return pthread_create(&(my_threads[id-1]), &thread_attr, start_routine, arg);
}


FINT_TYPE_C pthreads_jointhread(FINT_TYPE_C id)
{
    return pthreads_jointhread_(id);
}
FINT_TYPE_C pthreads_jointhread_(FINT_TYPE_C id)
{
    void *retval; // for convenience we do not pass it to fortran
    return pthread_join(my_threads[id-1], &retval);
}


FINT_TYPE_C pthreads_exitthread()
{
    return pthreads_exitthread_();
}
FINT_TYPE_C pthreads_exitthread_()
{
    pthread_exit(NULL);

    return 0;
}


FINT_TYPE_C pthreads_sched_yield()
{
    return pthreads_sched_yield_();
}
FINT_TYPE_C pthreads_sched_yield_()
{
    return sched_yield();
}

///////////////// RWLocks //////////////////////

FINT_TYPE_C rwlocks_init(FINT_TYPE_C numlocks)
{
    return rwlocks_init_(numlocks);
}
FINT_TYPE_C rwlocks_init_(FINT_TYPE_C numlocks)
{
    int iret = 0;
    int i = 0;

    maxnumlocks = numlocks;
    my_rwlocks    = (pthread_rwlock_t*)malloc(((unsigned int)maxnumlocks)*sizeof(pthread_rwlock_t));

    for (i=0;i<maxnumlocks;i++)
    {
      iret = pthread_rwlock_init(&my_rwlocks[i], NULL);
      CHECKRES;
    }

    return 0;
}


FINT_TYPE_C rwlocks_uninit()
{
    return rwlocks_uninit_();
}
FINT_TYPE_C rwlocks_uninit_()
{
    int iret = 0;
    int i = 0;


    for (i=0;i<maxnumlocks;i++)
    {
      iret = pthread_rwlock_destroy(&my_rwlocks[i]);
      CHECKRES;
    }

    free(my_rwlocks);

    return 0;
}


FINT_TYPE_C rwlocks_rdlock(FINT_TYPE_C id)
{
    return rwlocks_rdlock_(id);
}
FINT_TYPE_C rwlocks_rdlock_(FINT_TYPE_C id)
{
    return pthread_rwlock_rdlock(&my_rwlocks[id-1]);
}


FINT_TYPE_C rwlocks_wrlock(FINT_TYPE_C id)
{
    return rwlocks_wrlock_(id);
}
FINT_TYPE_C rwlocks_wrlock_(FINT_TYPE_C id)
{
    return pthread_rwlock_wrlock(&my_rwlocks[id-1]);
}


FINT_TYPE_C rwlocks_unlock(FINT_TYPE_C id)
{
    return rwlocks_unlock_(id);
}
FINT_TYPE_C rwlocks_unlock_(FINT_TYPE_C id)
{
    return pthread_rwlock_unlock(&my_rwlocks[id-1]);
}

FINT_TYPE_C get_my_tid()
{
    return get_my_tid_();
}
FINT_TYPE_C get_my_tid_()
{
    return (FINT_TYPE_C)syscall(SYS_gettid);
}

FINT_TYPE_C get_my_pid()
{
    return get_my_pid_();
}
FINT_TYPE_C get_my_pid_()
{
    return (FINT_TYPE_C)getpid();
}


