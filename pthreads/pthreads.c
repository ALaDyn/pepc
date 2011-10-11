
/*************************************************************************
>
>  pthreads interface library
>
*************************************************************************/

#include <pthread.h>
#include <stdlib.h>
#include <stdio.h>
#include <sched.h>
#include <unistd.h>
#include <sys/syscall.h>
#include <sys/types.h>

pthread_t *__restrict__ my_threads;
pthread_rwlock_t *my_rwlocks;
void* *my_thread_args; // array of void pointers for making packup copies of the pthread argument pointers
pthread_attr_t thread_attr;
int maxnumthreads = 0;
int maxnumlocks   = 0;

#define CHECKRES do {if (iret != 0) return iret;} while(0);


//////////////// PThreads //////////////////////

int pthreads_init(int numthreads)
{
    int iret = 0;

    maxnumthreads  = numthreads;
    my_threads     = (pthread_t*)malloc(((unsigned int)maxnumthreads)*sizeof(pthread_t));
    my_thread_args =      (void*)malloc(((unsigned int)maxnumthreads)*sizeof(void*));

    iret = pthread_attr_init(&thread_attr);
    CHECKRES;

    iret = pthread_attr_setdetachstate(&thread_attr, PTHREAD_CREATE_JOINABLE);
    CHECKRES;

    iret = pthread_attr_setscope(&thread_attr, PTHREAD_SCOPE_SYSTEM);
    CHECKRES;

    return 0;
}


int pthreads_uninit()
{
    int iret = 0;
    free(my_thread_args);
    free(my_threads);

    iret = pthread_attr_destroy(&thread_attr);
    CHECKRES;

    return 0;
}


int pthreads_createthread(int id, void *(*start_routine) (void *), void *arg, int relative_priority)
{
    // prepare a copy of the argument pointer to prevent it from being inaccessible when the thread actually starts
    my_thread_args[id] = arg;
printf("%p\n", my_thread_args[id]);
    return pthread_create(&(my_threads[id-1]), &thread_attr, start_routine, my_thread_args[id]);
}


int pthreads_jointhread(int id)
{
    void *retval; // for convenience we do not pass it to fortran
    return pthread_join(my_threads[id-1], &retval);
}


int pthreads_exitthread()
{
    pthread_exit(NULL);
    return 0;
}


int pthreads_sched_yield()
{
    return sched_yield();
}

///////////////// RWLocks //////////////////////

int rwlocks_init(int numlocks)
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


int rwlocks_uninit()
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


int rwlocks_rdlock(int id)
{
    return pthread_rwlock_rdlock(&my_rwlocks[id-1]);
}


int rwlocks_wrlock(int id)
{
    return pthread_rwlock_wrlock(&my_rwlocks[id-1]);
}


int rwlocks_unlock(int id)
{
    return pthread_rwlock_unlock(&my_rwlocks[id-1]);
}

///////////////// Utils //////////////////////

int get_my_tid()
{
    return (int)syscall(SYS_gettid);
}

int get_my_pid()
{
    return (int)getpid();
}


