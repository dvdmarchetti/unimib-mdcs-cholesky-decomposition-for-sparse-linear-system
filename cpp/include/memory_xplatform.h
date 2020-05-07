#ifndef MEMORY_XPLATFORM_H_
#define MEMORY_XPLATFORM_H_

#include <iostream>
#include <set>

#if defined(_WIN32) || defined(__CYGWIN__)
    #define OS_WIN
#elif defined(__linux__) || defined(unix) || defined(__unix__) || defined(__unix)
    #define OS_NIX
#else
    #error Unknown environment!
#endif

#ifdef OS_WIN
#include <windows.h>
#include <psapi.h>
#endif

#ifdef OS_NIX
#include <sys/types.h>
#include <sys/sysinfo.h>
#include <unistd.h>
#endif

namespace memory
{
    /**
     * Get total available virtual memory.
     */
    unsigned long long total_virtual();

    /**
     * Get total available physical memory.
     */
    unsigned long long total_physical();

    /**
     * Get the remaining available virtual memory.
     */
    unsigned long long current_virtual();

    /**
     * Get the remaning available physical memory.
     */
    unsigned long long current_physical();

    /**
     * Get the virtual memory used by the current process.
     */
    unsigned long long process_current_virtual();

    /**
     * Get the physical memory used by the current process.
     */
    unsigned long long process_current_physical();
};

#endif