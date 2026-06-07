/* Lightweight malloc interposer for memory profiling (LD_PRELOAD).
 *
 * Tracks live bytes / peak bytes / live allocation count, and a log2-size
 * histogram of LIVE allocations snapshotted AT THE PEAK. Reentrancy-safe:
 * uses malloc_usable_size() for accounting (no per-ptr map), a tiny static
 * bootstrap arena while dlsym resolves the real allocators, and fixed-size
 * atomic arrays for the histogram. Reports to stderr on exit.
 *
 * Goal: show *where* the memory peak comes from — how many buffers and of
 * what sizes are alive simultaneously at the high-water mark (the convolve
 * MemoryBank batch buffers).
 */
#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <dlfcn.h>
#include <stdatomic.h>
#include <malloc.h>
#include <unistd.h>
#include <sys/mman.h>

static void *(*real_malloc)(size_t)            = NULL;
static void  (*real_free)(void *)              = NULL;
static void *(*real_realloc)(void *, size_t)   = NULL;
static void *(*real_calloc)(size_t, size_t)    = NULL;
static int   (*real_posix_memalign)(void **, size_t, size_t) = NULL;
static void *(*real_aligned_alloc)(size_t, size_t) = NULL;
static void *(*real_memalign)(size_t, size_t)      = NULL;
static void *(*real_mmap)(void *, size_t, int, int, int, off_t) = NULL;
static int   (*real_munmap)(void *, size_t)        = NULL;
/* bytes attributed to explicit anonymous mmap (separate from malloc-family) */
static atomic_llong mmap_live = 0, mmap_peak = 0;

static atomic_llong cur_bytes  = 0;
static atomic_llong peak_bytes = 0;
static atomic_llong live_cnt   = 0;
/* histogram of LIVE allocation counts by floor(log2(usable_size)) */
static atomic_llong hist[64];
/* snapshot of (hist, cur, live) captured at each new peak */
static long long peak_hist[64];
static long long peak_live = 0;
static atomic_llong peak_snap_bytes = 0;
static atomic_int snap_lock = 0;

/* bootstrap arena for allocations issued during dlsym() init */
static char boot_arena[1 << 20];
static atomic_size_t boot_off = 0;
static int in_init = 0;

static int is_boot(void *p) {
    return (char *)p >= boot_arena && (char *)p < boot_arena + sizeof(boot_arena);
}
static void *boot_alloc(size_t n) {
    n = (n + 15) & ~((size_t)15);
    size_t off = atomic_fetch_add(&boot_off, n);
    if (off + n > sizeof(boot_arena)) { fprintf(stderr, "memtrace: boot arena overflow\n"); _exit(1); }
    return boot_arena + off;
}

static void init_real(void) {
    if (real_malloc) return;
    in_init = 1;
    real_malloc  = dlsym(RTLD_NEXT, "malloc");
    real_free    = dlsym(RTLD_NEXT, "free");
    real_realloc = dlsym(RTLD_NEXT, "realloc");
    real_calloc  = dlsym(RTLD_NEXT, "calloc");
    real_posix_memalign = dlsym(RTLD_NEXT, "posix_memalign");
    real_aligned_alloc  = dlsym(RTLD_NEXT, "aligned_alloc");
    real_memalign       = dlsym(RTLD_NEXT, "memalign");
    real_mmap           = dlsym(RTLD_NEXT, "mmap");
    real_munmap         = dlsym(RTLD_NEXT, "munmap");
    in_init = 0;
}

static inline int log2_bucket(size_t s) {
    int b = 0; while (s > 1) { s >>= 1; b++; } return b < 64 ? b : 63;
}

static atomic_llong last_snap = 0;
static void account_add(size_t usable) {
    long long c = atomic_fetch_add(&cur_bytes, (long long)usable) + (long long)usable;
    atomic_fetch_add(&live_cnt, 1);
    atomic_fetch_add(&hist[log2_bucket(usable)], 1);
    long long pk = atomic_load(&peak_bytes);
    if (c > pk) atomic_store(&peak_bytes, c);
    /* snapshot the live-histogram each +256MB of growth past the last snapshot,
       so the final snapshot sits within 256MB of the true peak */
    long long ls = atomic_load(&last_snap);
    if (c - ls > (256LL << 20) && !atomic_exchange(&snap_lock, 1)) {
        if (c - atomic_load(&last_snap) > (256LL << 20)) {
            for (int i = 0; i < 64; i++) peak_hist[i] = atomic_load(&hist[i]);
            peak_live = atomic_load(&live_cnt);
            atomic_store(&peak_snap_bytes, c);
            atomic_store(&last_snap, c);
        }
        atomic_store(&snap_lock, 0);
    }
}
static void account_sub(size_t usable) {
    atomic_fetch_sub(&cur_bytes, (long long)usable);
    atomic_fetch_sub(&live_cnt, 1);
    atomic_fetch_sub(&hist[log2_bucket(usable)], 1);
}

void *malloc(size_t size) {
    if (!real_malloc) { init_real(); if (in_init) return boot_alloc(size); }
    void *p = real_malloc(size);
    if (p) account_add(malloc_usable_size(p));
    return p;
}
void free(void *p) {
    if (!p || is_boot(p)) return;
    if (!real_free) init_real();
    account_sub(malloc_usable_size(p));
    real_free(p);
}
void *calloc(size_t n, size_t s) {
    if (!real_calloc) { init_real(); if (in_init) { void *p = boot_alloc(n*s); memset(p,0,n*s); return p; } }
    void *p = real_calloc(n, s);
    if (p) account_add(malloc_usable_size(p));
    return p;
}
void *realloc(void *old, size_t size) {
    if (!real_realloc) init_real();
    size_t oldu = (old && !is_boot(old)) ? malloc_usable_size(old) : 0;
    if (old && is_boot(old)) { void *p = malloc(size); if (p && oldu==0) memcpy(p, old, size); return p; }
    void *p = real_realloc(old, size);
    if (p) { if (oldu) account_sub(oldu); account_add(malloc_usable_size(p)); }
    return p;
}
int posix_memalign(void **out, size_t al, size_t size) {
    if (!real_posix_memalign) init_real();
    int rc = real_posix_memalign(out, al, size);
    if (rc == 0 && *out) account_add(malloc_usable_size(*out));
    return rc;
}
void *aligned_alloc(size_t al, size_t size) {
    if (!real_aligned_alloc) init_real();
    void *p = real_aligned_alloc(al, size);
    if (p) account_add(malloc_usable_size(p));
    return p;
}
void *memalign(size_t al, size_t size) {
    if (!real_memalign) init_real();
    void *p = real_memalign(al, size);
    if (p) account_add(malloc_usable_size(p));
    return p;
}
/* explicit anonymous mmap (big scratch some libs allocate directly) */
void *mmap(void *addr, size_t len, int prot, int flags, int fd, off_t off) {
    if (!real_mmap) init_real();
    void *p = real_mmap(addr, len, prot, flags, fd, off);
    if (p != MAP_FAILED && (flags & MAP_ANONYMOUS) && len >= (1u<<20)) {
        long long c = atomic_fetch_add(&mmap_live, (long long)len) + (long long)len;
        long long pk = atomic_load(&mmap_peak);
        if (c > pk) atomic_store(&mmap_peak, c);
    }
    return p;
}
int munmap(void *addr, size_t len) {
    if (!real_munmap) init_real();
    if (len >= (1u<<20)) atomic_fetch_sub(&mmap_live, (long long)len);
    return real_munmap(addr, len);
}

__attribute__((destructor))
static void report(void) {
    fprintf(stderr, "\n===== MEMTRACE peak report =====\n");
    fprintf(stderr, "peak_bytes = %lld (%.1f GB)\n", (long long)peak_bytes, peak_bytes / 1073741824.0);
    fprintf(stderr, "snapshot taken at cur=%.1f GB, live_allocs_in_snapshot = %lld\n",
            peak_snap_bytes / 1073741824.0, peak_live);
    fprintf(stderr, "LIVE-allocation size histogram at last snapshot (count x size-range -> total):\n");
    long long total = 0;
    for (int i = 0; i < 64; i++) {
        if (peak_hist[i] > 0) {
            double lo = (double)(1ULL << i);
            long long approx = (long long)(peak_hist[i] * lo * 1.5); /* mid-bucket estimate */
            total += approx;
            fprintf(stderr, "  2^%-2d (%8.1f KB) : %8lld allocs  ~%.1f GB\n",
                    i, lo/1024.0, peak_hist[i], approx/1073741824.0);
        }
    }
    fprintf(stderr, "  (sum of mid-bucket estimates ~%.1f GB)\n", total/1073741824.0);
    fprintf(stderr, "explicit-anon-mmap peak = %.1f GB (separate from malloc-family above)\n",
            mmap_peak / 1073741824.0);
    fprintf(stderr, "================================\n");
}
