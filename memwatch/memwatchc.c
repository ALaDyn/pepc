
/*************************************************************************
>
>  Memory bookkeeping routines: Keeps track of every malloc() and free()
> (as well as dynamic allocation in fortran as long as it is handled via
>  malloc and companions)
>  Global sum of allocated memory is stored in OverallMemUsage and can
>  be queried via internal_getmemusage().
>  The routines care for correct deallocation by keeping a linked list
>  of every allocated block.
>
>  modified example from
>    [ http://www.gnu.org/s/libc/manual/html_node/Hooks-for-Malloc.html ]
>  additional information from
>    [ http://www.gnu.org/s/libc/manual/html_node/Statistics-of-Malloc.html ]
>
>
*************************************************************************/

// uncomment to use fancy method of hooking into malloc routines. Otherwise, this
// encapsulates a simple call to mallinfo()
//#define USE_MALLOC_HOOKS

// uncomment for extensive debug information during runtime
//#define DEBUG

/* Prototypes for __malloc_hook, __free_hook */
#include <malloc.h>
#include <stdio.h>
#include "fortran2c_types.h"

#ifndef __MALLOC_PMT
  #warning "Memwatch: malloc-hooks not supported by libc. Switching to use of mallinfo() instead."
  #warning "Memwatch: results may be slightly inaccurate (deviations up to +/-stacksize may occur)."
  #undef USE_MALLOC_HOOKS
#endif

#define MB ((unsigned long long int)1024*1024)
#define NEGOTIATE_LIMIT_MB 1
#define MAX_SIZE ((unsigned long long int)((size_t)-1))

#define MIN(X,Y) ((X) < (Y) ? (X) : (Y))
#define MAX(X,Y) ((X) > (Y) ? (X) : (Y))


//> Macro to retrieve the current stack pointer by simply allocating a dummy variable
#define GET_STACK_PTR(res)                      \
         do {                                    \
               int stackDummy_ = 0;              \
               res = ((size_t)&(stackDummy_));   \
         } while (0)


// we need to store the initial stack pointer address
static size_t InitialStackPtr = 0;
static void GetInitialStackPtr(void) { GET_STACK_PTR(InitialStackPtr); }


/* Override initializing hook from the C library. */
#ifndef USE_MALLOC_HOOKS
void (*__malloc_initialize_hook) (void) = GetInitialStackPtr;
#else
/* Prototypes for our hooks.  */
static void my_init_hook (void);
static void my_init_hook_first (void);
static void *my_malloc_hook (size_t, __const __malloc_ptr_t);
static void *my_memalign_hook (size_t, size_t, __const __malloc_ptr_t);
static void *my_realloc_hook (void*, size_t, __const __malloc_ptr_t);
static void my_free_hook (void*, __const __malloc_ptr_t);
/* Backups of the old ones */
void *(*old_malloc_hook) __MALLOC_PMT ((size_t __size, __const __malloc_ptr_t));
void (*old_free_hook) __MALLOC_PMT ((void *__ptr, __const __malloc_ptr_t));
void *(*old_realloc_hook) __MALLOC_PMT ((void *__ptr, size_t __size, __const __malloc_ptr_t));
void *(*old_memalign_hook) __MALLOC_PMT ((size_t __alignment, size_t __size, __const __malloc_ptr_t));
/* Re-Override initializing hook from the C library. */
void (*__malloc_initialize_hook) (void) = my_init_hook_first;


struct list_node {
    void *mem_ptr;
    __const __malloc_ptr_t call_ptr;
    size_t size;
    struct list_node *next;
};


static size_t OverallMemUsage = 0;
struct list_node *liststart = NULL;



void add_to_list(void* ptr, size_t bytes, __const __malloc_ptr_t caller, char* reason)
{
    struct list_node* new_node = (struct list_node*)malloc(sizeof(struct list_node));
    new_node->mem_ptr  = ptr;
    new_node->call_ptr = caller;
    new_node->size     = bytes;
    new_node->next     = liststart;

    liststart          = new_node;

    OverallMemUsage   += bytes;

#ifdef DEBUG
    printf ("MALLOC_DEBUG: OverallMemUsage +%15u bytes = %15u for address %p by call from %p [%s]\n",
                   (unsigned int)new_node->size, (unsigned int)OverallMemUsage, ptr, caller, reason);
#endif
}



struct list_node* find_in_list(void* ptr, struct list_node** parent)
{

    if (liststart)
    {
       if (liststart->mem_ptr == ptr)
       {
           *parent = NULL;
           return liststart;
       }

       struct list_node* nextnode = liststart->next;
       struct list_node* lastnode = liststart;

       while (nextnode)
       {
           if (nextnode->mem_ptr == ptr)
           {
               *parent = lastnode;
               return nextnode;
           }

           lastnode = nextnode;
           nextnode = lastnode->next;
       }

    }
}



void update_list(void* ptr, size_t bytes, void* ptrnew, __const __malloc_ptr_t caller, char* reason)
{
    struct list_node* my_node   = NULL;
    struct list_node* my_parent = NULL;

    my_node = find_in_list(ptr, &my_parent);

    if (my_node)
    {
        OverallMemUsage = (OverallMemUsage - my_node->size) + bytes;

#ifdef DEBUG
   printf ("MALLOC_DEBUG: OverallMemUsage  %+15lld bytes = %15u for address %p (old address: %p) by call from %p [%s]\n",
                  (long long int)my_node->size-(long long int)bytes, (unsigned int)OverallMemUsage, ptrnew, ptr, caller, reason);
#endif

        my_node->call_ptr = caller;
        my_node->mem_ptr  = ptrnew;
        my_node->size     = bytes;
    }
#ifdef DEBUG
    else
      printf ("MALLOC_DEBUG: no appropriate block for address %p by call from %p [%s]\n", ptr, caller, reason);
#endif

}



void remove_from_list(void* ptr, __const __malloc_ptr_t caller, char* reason)
{
    struct list_node* my_node   = NULL;
    struct list_node* my_parent = NULL;

    if (!ptr) return; // who the hell tries to free a NULL-pointer?

    my_node = find_in_list(ptr, &my_parent);

    if (my_node)
    {
       if (my_parent)
           my_parent->next = my_node->next;
       else
           liststart = my_node->next;

#ifdef DEBUG
           printf ("MALLOC_DEBUG: OverallMemUsage -%15u bytes = %15u for address %p by call from %p [%s]\n",
                          (unsigned int)my_node->size, (unsigned int)OverallMemUsage, ptr, caller, reason);
#endif
           OverallMemUsage -= my_node->size;

           free(my_node);
    }
#ifdef DEBUG
    else
      printf ("MALLOC_DEBUG: Free could not find appropriate block for address %p by call from %p [%s]\n", ptr, caller, reason);
#endif
}


static void my_init_hook_first (void)
{
   GetInitialStackPtr();
   my_init_hook();
}


static void my_init_hook (void)
{
  /* Save underlying hooks */
  old_malloc_hook   = __malloc_hook;
  old_free_hook     = __free_hook;
  old_realloc_hook  = __realloc_hook;
  old_memalign_hook = __memalign_hook;
  /* Restore our own hooks */
  __malloc_hook   = my_malloc_hook;
  __free_hook     = my_free_hook;
  __realloc_hook  = my_realloc_hook;
  __memalign_hook = my_memalign_hook;
}



static void my_uninit_hook (void)
{
  /* Restore all old hooks */
  __malloc_hook   = old_malloc_hook;
  __free_hook     = old_free_hook;
  __realloc_hook  = old_realloc_hook;
  __memalign_hook = old_memalign_hook;
}



static void * my_malloc_hook (size_t size, __const __malloc_ptr_t caller)
{
  void *result = NULL;
  /* Restore all old hooks */
  my_uninit_hook();

  /* Call original malloc recursively */
  result = malloc (size);

  if (result)
  {
    /* internally, we use malloc, so hooks must be uninstalled here */
    add_to_list(result, size, caller, "malloc");
  }
#ifdef DEBUG
  else
  {
      printf ("MALLOC_DEBUG: allocating %u bytes of memory failed [malloc]\n", (unsigned int)size);
  }
#endif

  /* reinstall our hooks */
  my_init_hook();

  return result;
}



static void *my_memalign_hook (size_t alignment, size_t size, __const __malloc_ptr_t caller)
{
  void *result = NULL;

  /* Restore all old hooks */
  my_uninit_hook();

  /* Call original malloc recursively */
  result = memalign(alignment, size);

  if (result)
  {
    /* internally, we use malloc, so hooks must be uninstalled here */
    add_to_list(result, size, caller, "memalign");
  }
#ifdef DEBUG
  else
  {
    printf ("MALLOC_DEBUG: allocating %u bytes of memory failed [memalign]\n", (unsigned int)size);
  }
#endif

  /* reinstall our hooks */
  my_init_hook();

  return result;
}



static void * my_realloc_hook (void *ptr, size_t size, __const __malloc_ptr_t caller)
{
  void *result = NULL;
  /* Restore all old hooks */
  my_uninit_hook();

  /* Call recursively */
  result = realloc (ptr, size);

  /* internally, we use memory updating etc., so hooks must be uninstalled here */
  if (ptr)
  {
    if (size == 0)
      remove_from_list(ptr, caller, "realloc");
    else
      update_list(ptr, size, result, caller, "realloc");
  }
  else
    add_to_list(result, size, caller, "realloc");

  /* reinstall our hooks */
  my_init_hook();

  return result;
}



static void my_free_hook (void *ptr, __const __malloc_ptr_t caller)
{
  /* Restore all old hooks */
  my_uninit_hook();

  /* Call recursively */
  free (ptr);

  /* internally, we use free, so hooks must be uninstalled here */
  remove_from_list(ptr, caller, "free");

  /* reinstall our hooks */
  my_init_hook();
}

#endif



/* two versions for compatibility without caring for underscore-stuff :-) */
long long int internal_getmemusage()
{
    size_t stack_used = 0;

#ifdef USE_MALLOC_HOOKS
    size_t stack_ptr;
    GET_STACK_PTR(stack_ptr);

    if (stack_ptr < InitialStackPtr)
        stack_used = InitialStackPtr - stack_ptr;
    else
        stack_used = stack_ptr - InitialStackPtr;
#ifdef DEBUG
    printf("== stack extent  [%10u] --> [%10u]    --------> occupied bytes on stack: %10u\n", InitialStackPtr, stack_ptr, stack_used);
#endif
#endif

#ifndef USE_MALLOC_HOOKS
  struct mallinfo info = mallinfo();

#ifdef DEBUG
  printf("================= MALLOC_DEBUG: mallinfo-data ===================\n", info.arena);
  printf("non-mmapped space allocated from system:      %10u\n", info.arena);
  printf("number of free chunks:                        %10u\n", info.ordblks);
  printf("number of fastbin blocks:                     %10u\n", info.smblks);
  printf("number of mmapped regions:                    %10u\n", info.hblks);
  printf("space in mmapped regions:                     %10u\n", info.hblkhd);
  printf("maximum total allocated space:                %10u\n", info.usmblks);
  printf("space available in freed fastbin blocks:      %10u\n", info.fsmblks);
  printf("total allocated space:                        %10u\n", info.uordblks);
  printf("total free space:                             %10u\n", info.fordblks);
  printf("top-most, releasable (via malloc_trim) space: %10u\n", info.keepcost);
#endif

  return (long long int)(info.hblkhd + info.arena + stack_used);
#else
  return (long long int)(OverallMemUsage + stack_used);
#endif
}

long long int internal_getmemusage_()
{
  return internal_getmemusage();
}



FINT_TYPE_C internal_negotiatefreemem(FINT_TYPE_C* startval)
{
  FINT_TYPE_C lastsuccess;
  FINT_TYPE_C actval, newval, megabytes;
  void* testfield = NULL;
  unsigned long long int bytes;

#ifdef USE_MALLOC_HOOKS
  /* Restore all old hooks */
  my_uninit_hook();
#endif

  megabytes = MAX(*startval, 2*NEGOTIATE_LIMIT_MB);
  lastsuccess = -1;
  newval = megabytes / 2;
  actval = megabytes;

  while ((newval >= NEGOTIATE_LIMIT_MB) && (actval <= megabytes))
  {
      bytes = (unsigned long long int)actval*MB;
      if (bytes > MAX_SIZE)
      {
#ifdef DEBUG
          printf("internal_negotiatefreemem: requested larger block of memory (%15llu bytes) than malloc() can handle (%15llu bytes)\n", bytes, MAX_SIZE);
#endif
          bytes  = MAX_SIZE; // limit value to maximum allowed number of bytes to be allocated in current architecture
          actval = bytes/MB;
          newval = actval / 2;
      }

      testfield = malloc((size_t)bytes);

#ifdef DEBUG
    printf("internal_negotiatefreemem: actval = %12d MB,   bytes = %15llu,   testfield = %p\n", actval, bytes, testfield);
#endif

    if (testfield)
    {
      lastsuccess = actval;
      free(testfield);
      testfield = NULL;
      actval = actval + newval;
    }
    else
      actval = actval - newval;

    newval = newval / 2;
  }

#ifdef USE_MALLOC_HOOKS
  /* reinstall our hooks */
  my_init_hook();
#endif

  return lastsuccess;
}


FINT_TYPE_C internal_negotiatefreemem_(FINT_TYPE_C* startval)
{
    return internal_negotiatefreemem(startval);
}
































