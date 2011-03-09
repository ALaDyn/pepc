


#include <signal.h>
#include <unistd.h>
#include <stdio.h>
#include <execinfo.h>
#include <stdlib.h>


void dump_stack()
{
    void *bt_pointers[32];
    size_t size, i;
    char **bt_names;

    size = backtrace (bt_pointers, 32);
    bt_names = backtrace_symbols (bt_pointers, size);

    for (i = 0; i < size; i++)
    {
       printf ("%s\n", bt_names[i]);
    }

    free (bt_names);

    printf("\n\nYou can inspect the given adresses with the addr2line utlility:\n");
    printf("     addr2line -f -s -C -e programfile 0x1234567\n\n\n");
}


void signal_handler(int sig)
{
    // Process is asked to terminate
    printf("\n\nProcess has been requested to terminate by signal %d - Trying to dump stack\n\n", sig);
    dump_stack();

    exit(1);
}


void init_signal_handler(void)
{
  signal(SIGINT,  signal_handler);
  signal(SIGKILL,  signal_handler);
  signal(SIGQUIT, signal_handler);
  signal(SIGSEGV, signal_handler);
  signal(SIGTERM, signal_handler);
}

void init_signal_handler_(void)
{
    init_signal_handler();
}
