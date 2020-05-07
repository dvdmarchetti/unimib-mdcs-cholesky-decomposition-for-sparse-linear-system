#include <memory_xplatform.h>

/**
 * Get total available virtual memory.
 */
unsigned long long memory::total_virtual() {
#ifdef OS_WIN
    MEMORYSTATUSEX memInfo;
    memInfo.dwLength = sizeof(MEMORYSTATUSEX);
    GlobalMemoryStatusEx(&memInfo);

    return memInfo.ullTotalPageFile;
#elif defined(OS_NIX)
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
unsigned long long memory::total_physical() {
#ifdef OS_WIN
    MEMORYSTATUSEX memInfo;
    memInfo.dwLength = sizeof(MEMORYSTATUSEX);
    GlobalMemoryStatusEx(&memInfo);

    return memInfo.ullTotalPhys;
#elif defined(OS_NIX)
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
unsigned long long memory::current_virtual() {
#ifdef OS_WIN
    MEMORYSTATUSEX memInfo;
    memInfo.dwLength = sizeof(MEMORYSTATUSEX);
    GlobalMemoryStatusEx(&memInfo);

    return memInfo.ullTotalPageFile - memInfo.ullAvailPageFile;
#elif defined(OS_NIX)
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
unsigned long long memory::current_physical() {
#ifdef OS_WIN
    MEMORYSTATUSEX memInfo;
    memInfo.dwLength = sizeof(MEMORYSTATUSEX);
    GlobalMemoryStatusEx(&memInfo);

    return memInfo.ullTotalPhys - memInfo.ullAvailPhys;
#elif defined(OS_NIX)
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
unsigned long long memory::process_current_virtual() {
#ifdef OS_WIN
    PROCESS_MEMORY_COUNTERS_EX pmc;
    GetProcessMemoryInfo(GetCurrentProcess(), (PROCESS_MEMORY_COUNTERS*)&pmc, sizeof(pmc));

    return pmc.PrivateUsage;
#elif defined(OS_NIX)
    return 0L;
#endif
}

/**
 * Get the physical memory used by the current process.
 */
unsigned long long memory::process_current_physical() {
#ifdef OS_WIN
    PROCESS_MEMORY_COUNTERS_EX pmc;
    GetProcessMemoryInfo(GetCurrentProcess(), (PROCESS_MEMORY_COUNTERS*)&pmc, sizeof(pmc));

    return pmc.WorkingSetSize;
#elif defined(OS_NIX)
    long rss = 0L;
    FILE* fp = NULL;

    if ((fp = fopen("/proc/self/statm", "r")) == NULL) {
        return (size_t)0L;
    }

    if (fscanf(fp, "%*s%ld", &rss) != 1) {
        fclose(fp);
        return (size_t) 0L;
    }

    fclose(fp);
    return (size_t) rss * (size_t) sysconf(_SC_PAGESIZE);
#endif
}