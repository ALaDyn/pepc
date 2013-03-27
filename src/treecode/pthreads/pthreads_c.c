/*
* This file is part of PEPC - The Pretty Efficient Parallel Coulomb Solver.
* 
* Copyright (C) 2002-2013 Juelich Supercomputing Centre, 
*                         Forschungszentrum Juelich GmbH,
*                         Germany
* 
* PEPC is free software: you can redistribute it and/or modify
* it under the terms of the GNU Lesser General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
* 
* PEPC is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU Lesser General Public License for more details.
* 
* You should have received a copy of the GNU Lesser General Public License
* along with PEPC.  If not, see <http://www.gnu.org/licenses/>.
*/


/*************************************************************************
>
>  pthreads interface library
>
*************************************************************************/

#include <pthread.h>
#include <errno.h>
#include <stdlib.h>
#include <stdio.h>
#include <sched.h>
#include <unistd.h>
#include <sys/syscall.h>
#include <sys/types.h>

pthread_rwlock_t *my_rwlocks;
pthread_attr_t thread_attr;
int maxnumthreads  = 0;
int maxnumlocks    = 0;

int RWLOCKS_BUSY = EBUSY;

#define CHECKRES do {if (iret != 0) return iret;} while(0);


//////////////// PThreads //////////////////////

int pthreads_init()
{
    int iret = 0;

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

    iret = pthread_attr_destroy(&thread_attr);
    CHECKRES;

    return 0;
}


pthread_t* pthreads_alloc_thread()
{
  return (pthread_t*)malloc(sizeof(pthread_t));
}


void pthreads_free_thread(pthread_t *storage)
{
  free(storage);
}


int pthreads_createthread(pthread_t *thread, void *(*start_routine) (void *), void *arg, int relative_priority)
{
    // prepare a copy of the argument pointer to prevent it from being inaccessible when the thread actually starts
    // my_thread_args[id] = arg;
    //return pthread_create(thread, &thread_attr, start_routine, my_thread_args[id]);
    return pthread_create(thread, &thread_attr, start_routine, arg);
}


int pthreads_jointhread(pthread_t *thread)
{
    void *retval; // for convenience we do not pass it to fortran
    return pthread_join(*thread, &retval);
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


int rwlocks_tryrdlock(int id)
{
    return pthread_rwlock_tryrdlock(&my_rwlocks[id-1]);
}


int rwlocks_wrlock(int id)
{
    return pthread_rwlock_wrlock(&my_rwlocks[id-1]);
}


int rwlocks_trywrlock(int id)
{
    return pthread_rwlock_trywrlock(&my_rwlocks[id-1]);
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


