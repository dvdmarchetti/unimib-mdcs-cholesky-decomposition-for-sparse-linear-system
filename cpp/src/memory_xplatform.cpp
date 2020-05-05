#include <memory_xplatform.h>

/**
 * Get total available virtual memory.
 */
DWORDLONG memory::total_virtual()
{
#ifdef OS_WIN
    MEMORYSTATUSEX memInfo;
    memInfo.dwLength = sizeof(MEMORYSTATUSEX);
    GlobalMemoryStatusEx(&memInfo);

    return memInfo.ullTotalPageFile;
#endif

#ifdef OS_NIX
    struct sysinfo memInfo;
    sysinfo(&memInfo);

    long long totalVirtual = memInfo.totalram;
    totalVirtual + memInfo.totalswap;
    totalVirtual *= memInfo.mem_unit;

    return totalVirtual;
#endif
}

/**
 * Get total available physical memory.
 */
DWORDLONG memory::total_physical()
{
#ifdef OS_WIN
    MEMORYSTATUSEX memInfo;
    memInfo.dwLength = sizeof(MEMORYSTATUSEX);
    GlobalMemoryStatusEx(&memInfo);

    return memInfo.ullTotalPhys;
#endif

#ifdef OS_NIX
    struct sysinfo memInfo;
    sysinfo(&memInfo);

    long long totalPhysical = memInfo.totalram;
    totalPhysical *= memInfo.mem_unit;

    return totalPhysical;
#endif
}

/**
 * Get the remaining available virtual memory.
 */
DWORDLONG memory::current_virtual()
{
#ifdef OS_WIN
    MEMORYSTATUSEX memInfo;
    memInfo.dwLength = sizeof(MEMORYSTATUSEX);
    GlobalMemoryStatusEx(&memInfo);

    return memInfo.ullTotalPageFile - memInfo.ullAvailPageFile;
#endif

#ifdef OS_NIX
    struct sysinfo memInfo;
    sysinfo(&memInfo);

    long long currentVirtual = memInfo.totalram - memInfo.freeram;
    currentVirtual += memInfo.totalswap - memInfo.freeswap;
    currentVirtual *= memInfo.mem_unit;

    return currentVirtual;
#endif
}

/**
 * Get the remaning available physical memory.
 */
DWORDLONG memory::current_physical()
{
#ifdef OS_WIN
    MEMORYSTATUSEX memInfo;
    memInfo.dwLength = sizeof(MEMORYSTATUSEX);
    GlobalMemoryStatusEx(&memInfo);

    return memInfo.ullTotalPhys - memInfo.ullAvailPhys;
#endif

#ifdef OS_NIX
    struct sysinfo memInfo;
    sysinfo(&memInfo);

    long long currentPhysical = memInfo.totalram - memInfo.freeram;
    currentPhysical *= memInfo.mem_unit;

    return currentPhysical;
#endif
}

/**
 * Get the virtual memory used by the current process.
 */
SIZE_T memory::process_current_virtual()
{
#ifdef OS_WIN
    PROCESS_MEMORY_COUNTERS_EX pmc;
    GetProcessMemoryInfo(GetCurrentProcess(), (PROCESS_MEMORY_COUNTERS*)&pmc, sizeof(pmc));

    return pmc.PrivateUsage;
#endif

#ifdef OS_NIX
    //
#endif
}

/**
 * Get the physical memory used by the current process.
 */
SIZE_T memory::process_current_physical()
{
#ifdef OS_WIN
    PROCESS_MEMORY_COUNTERS_EX pmc;
    GetProcessMemoryInfo(GetCurrentProcess(), (PROCESS_MEMORY_COUNTERS*)&pmc, sizeof(pmc));

    return pmc.WorkingSetSize;
#endif

#ifdef OS_NIX
    //
#endif
}