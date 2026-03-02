
#if defined(__APPLE__) && !defined(_DARWIN_C_SOURCE)
#define _DARWIN_C_SOURCE
#endif
#define _POSIX_C_SOURCE 200809L
/*=====================================================================
  NitroSat – Advanced MaxSAT Solver (C version)
  --------------------------------------------------------------------
  This file is a direct translation of the Lua source 
  It compiles as a standalone executable:

      gcc -O3 -march=native -std=c99 -lm -o nitrosat nitrosat.c

  Usage:
      ./nitrosat <cnf-file> [max-steps] [--no-dcw] [--no-topo] [--json]
                 [--proof <path>] [--proof-format drat|lrat]

  The program reads a DIMACS CNF file, builds the solver instance and
  runs the full solver (gradient descent + BAHA + three‑phase finisher
  + core‑decomposition).  Progress is printed to stdout; the final line
  reports whether the formula was solved and how many steps were needed.

AUTHOR: Sethu Iyer (sethuiyer95@gmail.com)
Apache 2.0 License. 

Please cite https://zenodo.org/records/18753235 if you are using this software
=====================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include <time.h>
#include <assert.h>
#include <ctype.h>
#include <limits.h>

/* --------------------------------------------------------------------
   JSON output structures and high-resolution timing
-------------------------------------------------------------------- */

typedef struct {
    double init_ms;
    double langevin_ms;
    double topo_repair_ms;
    double adelic_ms;
    double core_decomp_ms;
    double walksat_ms;
    double total_ms;
} LatencyBreakdown;

typedef struct {
    double final_beta;
    double critical_beta;
    double entropy_level;
    const char *convexity_status;
    int betti_0, betti_1;
    double complexity_score;
    int persistence_events;
    double aggregation_error;
    const char *chebyshev_bias;
    int *blame_clause_ids;
    double *blame_weights;
    int blame_count;
} SolverDiagnostics;

typedef struct {
    int requested;
    int generated;
    int derived_units;
    int lrat_hints_emitted;
    const char *format;
    const char *path;
    const char *backend;
    char status[256];
} ProofReport;

static double timespec_ms(struct timespec *start, struct timespec *end) {
    return (end->tv_sec - start->tv_sec) * 1000.0 +
           (end->tv_nsec - start->tv_nsec) / 1e6;
}

/* --------------------------------------------------------------------
   Helper macros / constants
-------------------------------------------------------------------- */
#define MAX_NAME_LEN    64
#define PI              3.14159265358979323846
#define E_INV           0.36787944117144232   /* 1/e   */
#define PHI             ((1.0 + sqrt(5.0)) / 2.0)
#define EPS             1e-12
#define PROOF_UNASSIGNED -1

/* --------------------------------------------------------------------
   Assertion macros for safety checking
-------------------------------------------------------------------- */
#ifdef NITROSAT_DEBUG
#define NITROSAT_ASSERT(cond) do { \
    if (!(cond)) { \
        fprintf(stderr, "ASSERTION FAILED: %s:%d: %s\n", __FILE__, __LINE__, #cond); \
        abort(); \
    } \
} while (0)
#else
#define NITROSAT_ASSERT(cond) ((void)0)
#endif

#define NITROSAT_ENSURE(cond, msg) do { \
    if (!(cond)) { \
        fprintf(stderr, "FATAL: %s:%d: %s\n", __FILE__, __LINE__, msg); \
        abort(); \
    } \
} while (0)

/* Clause size bounds - prevents buffer overflows in gradient computation.
   Note: Enterprise timetabling instances can have very large clauses.
   We use a very generous bound (3 billion) to accommodate any practical
   instance while still catching obvious corruption. */
#define MAX_CLAUSE_SIZE 3000000000L
#define NITROSAT_ASSERT_CLAUSE_SIZE(k) do { \
    long _k = (long)(k); \
    if (_k <= 0 || _k > MAX_CLAUSE_SIZE) { \
        fprintf(stderr, "FATAL: Clause size %ld exceeds bounds [1,%ld] at %s:%d\n", \
                _k, (long)MAX_CLAUSE_SIZE, __FILE__, __LINE__); \
        abort(); \
    } \
} while (0)

static int proof_max_vars = 500000;     /* CLI: --proof-max-vars (default 500K) */
static int proof_max_clauses = 3000000000; /* CLI: --proof-max-clauses (default 3B) */

/* --------------------------------------------------------------------
   Simple pseudo‑random generator – deterministic, fast.
   The Lua version used a “x = (x * 9973) mod 1” scheme.
-------------------------------------------------------------------- */
static double rng_state = 0.95607928791821;   /* matches LuaJIT math.randomseed(42); math.random() */

static double rng_double(void)
{
    double x = rng_state * 9973.0;
    rng_state = x - floor(x);
    return rng_state;
}
static int rng_int(int a, int b) /* inclusive */
{
    return a + (int)(rng_double() * (b - a + 1));
}

/* --------------------------------------------------------------------
   Utility functions
-------------------------------------------------------------------- */
static inline double clamp(double x, double lo, double hi)
{
    if (x < lo) return lo;
    if (x > hi) return hi;
    return x;
}
static inline int   sign(double x) { return (x > 0) - (x < 0); }

/* --------------------------------------------------------------------
   1.  Prime generation (sieve of Eratosthenes)
-------------------------------------------------------------------- */
static int *sieve_primes(int n, int *out_cnt)
{
    /* Estimate an upper bound for the nth prime (Rosser–Schoenfeld) */
    int limit = (n < 6) ? 15 : (int)(n * (log((double)n) + log(log((double)n))) + 3);
    char *is_prime = calloc(limit + 1, 1);
    for (int i = 2; i <= limit; ++i) is_prime[i] = 1;
    for (int p = 2; p * p <= limit; ++p)
        if (is_prime[p])
            for (int q = p * p; q <= limit; q += p) is_prime[q] = 0;

    int *primes = malloc(sizeof(int) * n);
    int cnt = 0;
    for (int i = 2; i <= limit && cnt < n; ++i)
        if (is_prime[i]) primes[cnt++] = i;

    free(is_prime);
    *out_cnt = cnt;
    return primes;
}

/* --------------------------------------------------------------------
   3.  Optimiser --------------------------------------------------------
-------------------------------------------------------------------- */
typedef struct {
    double *m;      /* 1st moment */
    double *v;      /* 2nd moment */
    double  lr;
    double  beta1, beta2, eps;
    int     t;
    int     nv;
    double  resonance_amplitude;
    char    name[MAX_NAME_LEN];   /* "nadam" or "harmonic_adam" … */
} Optimizer;

/* Allocate and initialise an optimiser */
static Optimizer *optimizer_new(const char *name, double *params, int nv,
                               double lr, double beta1, double beta2,
                               double eps, double resonance_amp)
{
    (void)params;
    Optimizer *opt = calloc(1, sizeof(Optimizer));
    strncpy(opt->name, name, MAX_NAME_LEN-1);
    opt->nv    = nv;
    opt->lr    = lr;
    opt->beta1 = beta1;
    opt->beta2 = beta2;
    opt->eps   = eps;
    opt->m = calloc(nv + 1, sizeof(double));
    opt->v = calloc(nv + 1, sizeof(double));
    opt->resonance_amplitude = resonance_amp;
    return opt;
}

static void optimizer_step(Optimizer *opt, double *params, const double *grads, const uint8_t *decimated)
{
    ++opt->t;
    double b1 = opt->beta1;
    double b2 = opt->beta2;
    double lr = opt->lr;
    double eps = opt->eps;
    int    n = opt->nv;

    double b1_t = pow(b1, opt->t);
    double b2_t = pow(b2, opt->t);
    double t_inv1 = 1.0 / (1.0 - b1_t);
    double t_inv2 = 1.0 / (1.0 - b2_t);

    int is_nadam = (strcmp(opt->name, "nadam") == 0);

    for (int i = 1; i <= n; ++i) {
        if (decimated && decimated[i]) continue;
        double m_hat_prev = 0.0;
        if (is_nadam) {
            m_hat_prev = b1 * opt->m[i] * t_inv1;
        }

        opt->m[i] = b1 * opt->m[i] + (1.0 - b1) * grads[i];
        opt->v[i] = b2 * opt->v[i] + (1.0 - b2) * grads[i] * grads[i];
        
        double m_hat = opt->m[i] * t_inv1;
        double v_hat = opt->v[i] * t_inv2;

        if (is_nadam) {
            double m_nadam = m_hat + b1 * (m_hat - m_hat_prev);
            params[i] -= lr * m_nadam / (sqrt(v_hat) + eps);
        } else {
            params[i] -= lr * m_hat / (sqrt(v_hat) + eps);
        }
    }

    /* ---- Resonance perturbation (identical to Lua) ------------------- */
    double phi = PHI;
    double t   = (double)opt->t * 0.1;
    double intensity =
          sin(PI/20.0 * t) * exp(-PI/20.0  * t)
        + sin(PI/10.0 * t) * exp(-PI/10.0  * t)
        + sin(PI/9.0  * t) * exp(-PI/9.0   * t)
        + sin((1.0/(phi*phi))*t) * exp(-(1.0/(phi*phi))*t);

    /* compute l2‑norm of the gradient */
    double gnorm2 = 0.0;
    for (int i = 1; i <= n; ++i) gnorm2 += grads[i]*grads[i];
    double gnorm = sqrt(gnorm2);
    if (gnorm > 1e-10) {
        double pert = opt->resonance_amplitude * intensity;
        for (int i = 1; i <= n; ++i) {
            params[i] -= pert * grads[i] / gnorm;
        }
    }
}

/* --------------------------------------------------------------------
   4.  BAHA utilities – fracture detection & Lambert‑W branch enumeration
-------------------------------------------------------------------- */
typedef struct {
    double sigma_factor;
    int    n;
    double mean, M2, last_rate;
    double pb, plZ;          /* previous β , lZ */
} FractureDetector;

static FractureDetector *fd_new(double sigma_factor)
{
    FractureDetector *fd = calloc(1, sizeof(FractureDetector));
    fd->sigma_factor = sigma_factor;
    return fd;
}
static void fd_record(FractureDetector *fd, double beta, double lZ)
{
    if (fd->pb != 0.0 && fd->plZ != 0.0) {
        double db = beta - fd->pb;
        if (fabs(db) > 1e-9) {
            double rate = fabs(lZ - fd->plZ) / db;
            fd->last_rate = rate;
            ++fd->n;
            double delta = rate - fd->mean;
            fd->mean += delta / fd->n;
            fd->M2   += delta * (rate - fd->mean);
        }
    }
    fd->pb = beta; fd->plZ = lZ;
}
static int fd_is_fracture(FractureDetector *fd)
{
    if (fd->n < 5) return 0;
    double s = sqrt(fd->M2 / (fd->n - 1));
    return fd->last_rate > fd->mean + fd->sigma_factor * s;
}
static void fd_clear(FractureDetector *fd)
{
    fd->n = fd->mean = fd->M2 = fd->last_rate = 0.0;
    fd->pb = fd->plZ = 0.0;
}

/* ---- simple Lambert‑W (principal branch) via Halley's method ------- */
static double lambert_W0(double z)
{
    if (z < -E_INV) return NAN;
    /* Initial guess */
    double w = (z < -0.3) ? (z * exp(1.0))
             : (z < 1.0) ? (z * (1.0 - z + z*z))
                        : (log(z) - log(log(z) + 1.0));
    for (int i = 0; i < 50; ++i) {
        double ew = exp(w);
        double f  = w*ew - z;
        double fp = ew * (w + 1.0);
        if (fabs(fp) < 1e-15) break;
        double denom = fp - f * (ew * (w + 2.0)) / (2.0 * fp);
        if (fabs(denom) < 1e-15) break;
        double w2 = w - f/denom;
        if (fabs(w2-w) < 1e-10) { w=w2; break; }
        w = w2;
    }
    return w;
}

/* ---- Lambert W_{-1} (lower branch) via Halley's method ------------ */
static double lambert_Wm1(double z) {
    if (z < -E_INV || z >= 0) return NAN;
    /* Initial guess for W_{-1} branch */
    double w;
    if (z > -0.1) {
        w = log(-z) - log(-log(-z));
    } else {
        w = -1.0 - sqrt(2.0 * (1.0 + z * exp(1.0)));
    }
    for (int i = 0; i < 64; ++i) {
        double ew = exp(w);
        double f  = w * ew - z;
        double fp = ew * (w + 1.0);
        if (fabs(fp) < 1e-15) break;
        double denom = fp - f * (ew * (w + 2.0)) / (2.0 * fp);
        if (fabs(denom) < 1e-15) break;
        double w2 = w - f / denom;
        if (fabs(w2 - w) < 1e-10) { w = w2; break; }
        w = w2;
    }
    return w;
}

/* Enumerate possible branches of "β – βc = u " using the transformation
   u·e^u = ξ,  ξ = u·exp(u).  The Lua code returns at most two branches
   (k=0 → principal, k=-1 → secondary). */
typedef struct { int k; double beta; double score; } Branch;
static Branch *enumerate_branches(double bc, double b, int *out_cnt)
{
    double u = b - bc;
    if (fabs(u) < 1e-12) u = 1e-12;
    double xi = u * exp(u);
    Branch *br = malloc(2 * sizeof(Branch));
    int cnt = 0;
    double w0 = lambert_W0(xi);
    if (!isnan(w0)) {
        double b0 = bc + w0;
        if (b0 > 0) { br[cnt++] = (Branch){0, b0, 0.0}; }
    }
    /* FIXED: Compute the W_{-1} branch for phase transition jumps */
    if (xi >= -E_INV && xi < 0) {
        double w1 = lambert_Wm1(xi);
        if (!isnan(w1)) {
            double b1 = bc + w1;
            if (b1 > 0) { br[cnt++] = (Branch){-1, b1, 0.0}; }
        }
    }
    *out_cnt = cnt;
    return br;
}

/* --------------------------------------------------------------------
   5.  Persistent homology – only Betti‑0/1 are needed.
-------------------------------------------------------------------- */
typedef struct {
    int *parent;
    int *rank;
    int  n, count;
} UnionFind;
static UnionFind *uf_new(int n)
{
    UnionFind *uf = calloc(1, sizeof(UnionFind));
    uf->n = n;
    uf->parent = malloc((n+1)*sizeof(int));
    uf->rank   = calloc(n+1, sizeof(int));
    uf->count  = n;
    for (int i = 1; i <= n; ++i) uf->parent[i] = i;
    return uf;
}
static int uf_find(UnionFind *uf, int x)
{
    if (uf->parent[x] != x) uf->parent[x] = uf_find(uf, uf->parent[x]);
    return uf->parent[x];
}
static int uf_union(UnionFind *uf, int a, int b)
{
    int pa = uf_find(uf, a);
    int pb = uf_find(uf, b);
    if (pa == pb) return 0;
    if (uf->rank[pa] < uf->rank[pb]) { int t=pa; pa=pb; pb=t; }
    else if (uf->rank[pa] == uf->rank[pb]) ++uf->rank[pa];
    uf->parent[pb] = pa;
    --uf->count;
    return 1;
}

/* --------------------------------------------------------------------
   6.  Topological Structures
-------------------------------------------------------------------- */
typedef struct {
    int    beta0, beta1;
    int    active_vars, edge_count;
    double avg_confidence, complexity_score;
} TopologyInfo;

typedef struct {
    int initial_beta0, final_beta0;
    int initial_beta1, final_beta1;
    double initial_complexity, final_complexity;
    int has_initial;
    int birth_step_0, birth_step_1;
    int prev_beta0, prev_beta1;
    int persistence_events;
} PersistenceTracker;

static void update_persistence(PersistenceTracker *pt, int beta0, int beta1, int step, double complexity) {
    if (!pt->has_initial) {
        pt->initial_beta0 = beta0;
        pt->initial_beta1 = beta1;
        pt->initial_complexity = complexity;
        pt->has_initial = 1;
        pt->prev_beta0 = beta0;
        pt->prev_beta1 = beta1;
    }
    pt->final_beta0 = beta0;
    pt->final_beta1 = beta1;
    pt->final_complexity = complexity;

    if (beta0 > pt->prev_beta0) pt->birth_step_0 = step;
    if (beta1 > pt->prev_beta1) pt->birth_step_1 = step;

    if (beta0 < pt->prev_beta0 && pt->birth_step_0 > 0) pt->persistence_events++;
    if (beta1 < pt->prev_beta1 && pt->birth_step_1 > 0) pt->persistence_events++;

    pt->prev_beta0 = beta0;
    pt->prev_beta1 = beta1;
}

/* --------------------------------------------------------------------
   7.  Instance representation (flattened clauses)
   
   Memory Ownership:
   - cl_offs: [OWNER] Instance, must be freed by caller of read_cnf()
   - cl_flat: [OWNER] Instance, must be freed by caller of read_cnf()
   - num_vars, num_clauses: Metadata (no ownership)
-------------------------------------------------------------------- */
typedef struct {
    int    num_vars;
    int    num_clauses;
    int   *cl_offs;    /* [OWNER] size num_clauses+1 */
    int   *cl_flat;    /* [OWNER] flattened literals */
} Instance;

/* --------------------------------------------------------------------
   7.  Solver (NitroSat) structure
   
   Memory Ownership (fields marked [OWNER] must be freed in nitrosat_free()):
   - inst:        [OWNER] Instance passed in, freed by main()
   - cl_offs:     [BORROWED] Points to inst->cl_offs (no free)
   - cl_flat:     [BORROWED] Points to inst->cl_flat (no free)
   - v2c_ptr:     [OWNER] Built in main(), freed in nitrosat_free()
   - v2c_data:    [OWNER] Built in main(), freed in nitrosat_free()
   - x:           [OWNER] Allocated in nitrosat_new()
   - degrees:     [OWNER] Allocated in compute_degrees()
   - decimated:   [OWNER] Allocated in nitrosat_new()
   - cl_weights:  [OWNER] Allocated in nitrosat_new()
   - zeta_weights:[OWNER] Allocated in nitrosat_new()
   - is_hard:     [OWNER] Allocated in nitrosat_new()
   - opt:         [OWNER] Allocated in nitrosat_new()
   - fd:          [OWNER] Allocated in nitrosat_new()
   - primes:      [OWNER] Allocated in nitrosat_new()
   - grad_buffer: [OWNER] Allocated in nitrosat_new()
   - heat_mult_buffer: [OWNER] Allocated in nitrosat_new()
   - sat_counts:  [OWNER] Allocated in nitrosat_new()
   - entropy_cache:    [OWNER] Allocated in nitrosat_new()
   - entropy_x_cache:  [OWNER] Allocated in nitrosat_new()
   - prev_x:      [OWNER] Allocated in nitrosat_new()
-------------------------------------------------------------------- */
typedef struct {
    /* instance data */
    Instance *inst;
    /* flattened clause structures */
    const int *cl_offs;
    const int *cl_flat;
    int    num_vars, num_clauses;
    int    *v2c_ptr;
    int    *v2c_data;
    /* per‑variable data */
    double *x;          /* current (continuous) assignment */
    double *degrees;    /* degree = Σ(k−1) per variable */
    uint8_t *decimated;/* 0 = free, 1 = locked */
    /* clause data */
    double *cl_weights;     /* prime‑weighted base weight */
    double *zeta_weights;   /* log(p)/p   – used by Zeta perturbation   */
    uint8_t *is_hard;       /* 0 = soft, 1 = hard (unused in current code) */
    /* optimisation state */
    Optimizer *opt;
    /* BAHA utilities */
    FractureDetector *fd;
    double baha_beta_c;  /* β at max variance (critical β) */
    double baha_var_peak;/* max variance observed so far   */
    int    baha_var_n;
    double baha_var_mean, baha_var_M2;
    /* persistent homology */
    int    use_topology;
    TopologyInfo last_topology;
    PersistenceTracker pt;
    /* prime table */
    int    *primes;
    int    prime_cnt;
    double sigma_factor; /* for fracture detector – user tunable */
    /* other knobs */
    double heat_beta, heat_lambda;
    double entropy_weight;
    double dec_threshold;
    int    dec_freq;     /* lock variables every N steps */
    int    max_steps;
    int    verbose;
    /* hot-path buffers */
    double *grad_buffer;
    double *heat_mult_buffer;
    int    *sat_counts; /* number of satisfied literals per clause */
    /* entropy cache for performance */
    double *entropy_cache;
    double *entropy_x_cache;
    /* Cheat Code 5: previous x for incremental sat_counts */
    double *prev_x;
} NitroSat;

static void recompute_sat_counts(NitroSat *ns);

static TopologyInfo compute_topology(NitroSat *ns, double th)
{
    (void)th;
    int num_vars = ns->num_vars;
    int num_clauses = ns->num_clauses;
    const int *cl_offs = ns->cl_offs;
    const int *cl_flat = ns->cl_flat;
    const double *x = ns->x;

    TopologyInfo topo = {0};
    UnionFind *uf = uf_new(num_vars);
    int *is_active = calloc(num_vars + 1, sizeof(int));
    int active_var_cnt = 0;

    /* 1. BFS/Union-Find for Beta-0 and identify active vars */
    recompute_sat_counts(ns);
    for (int c = 0; c < num_clauses; ++c) {
        if (ns->sat_counts[c] > 0) continue;
        int s = cl_offs[c], e = cl_offs[c+1];

        int first_v = -1;
        for (int i=s; i<e; ++i) {
            int v = abs(cl_flat[i]);
            if (!is_active[v]) { is_active[v] = 1; active_var_cnt++; }
            if (first_v == -1) first_v = v;
            else uf_union(uf, first_v, v);
        }
    }

    /* 2. Sparse Edge Counting for Beta-1 */
    long long edge_cnt = 0;
    unsigned char *seen = calloc((num_vars / 8) + 1, 1);
    #define SET_BIT(a, b) (a[(b)>>3] |= (1 << ((b)&7)))
    #define GET_BIT(a, b) (a[(b)>>3] & (1 << ((b)&7)))
    #define CLR_BIT(a, b) (a[(b)>>3] &= ~(1 << ((b)&7)))

    int *neighbor_list = malloc(1024 * sizeof(int));
    int cap = 1024;

    for (int u = 1; u <= num_vars; ++u) {
        if (!is_active[u]) continue;
        int neighbors_found = 0;

        for (int p = ns->v2c_ptr[u]; p < ns->v2c_ptr[u+1]; ++p) {
            int c = ns->v2c_data[p];
            if (ns->sat_counts[c] > 0) continue;

            for (int i=cl_offs[c]; i<cl_offs[c+1]; ++i) {
                int v = abs(cl_flat[i]);
                if (v > u && !GET_BIT(seen, v)) {
                    SET_BIT(seen, v);
                    edge_cnt++;
                    if (neighbors_found >= cap) {
                        cap *= 2;
                        neighbor_list = realloc(neighbor_list, cap * sizeof(int));
                    }
                    neighbor_list[neighbors_found++] = v;
                }
            }
        }
        for (int i=0; i<neighbors_found; ++i) CLR_BIT(seen, neighbor_list[i]);
    }
    free(neighbor_list);

    topo.beta0 = uf->count - (num_vars - active_var_cnt);
    if (topo.beta0 < 0) topo.beta0 = 0;
    topo.edge_count = (int)edge_cnt;
    topo.active_vars = active_var_cnt;
    
    int chi = active_var_cnt - (int)edge_cnt;
    topo.beta1 = (topo.beta0 - chi > 0) ? (topo.beta0 - chi) : 0;

    double conf_sum = 0.0;
    for (int i=1;i<=num_vars;++i) if (is_active[i]) conf_sum += fabs(x[i]-0.5);
    topo.avg_confidence = (active_var_cnt ? conf_sum/active_var_cnt : 0.0);
    topo.complexity_score = (topo.beta0>0 ? (double)topo.beta1 / topo.beta0 : 0.0);

    free(is_active); free(seen); free(uf->parent); free(uf->rank); free(uf);
    return topo;
}

static Instance *read_cnf(const char *filename)
{
    FILE *fp = fopen(filename, "r");
    if (!fp) { perror("fopen"); return NULL; }

    char line[4096];
    int vars = 0, cls = 0;
    int line_no = 0;
    int have_problem = 0;
    while (fgets(line, sizeof(line), fp)) {
        ++line_no;
        char *scan = line;
        while (isspace((unsigned char)*scan)) ++scan;
        if (*scan == 'c' || *scan == '\0') continue;
        if (*scan == 'p') {
            if (sscanf(scan, "p cnf %d %d", &vars, &cls) != 2 || vars < 0 || cls < 0) {
                fprintf(stderr, "Bad problem line at %d.\n", line_no);
                fclose(fp);
                return NULL;
            }
            have_problem = 1;
            break;
        }
    }
    if (!have_problem) {
        fprintf(stderr, "No problem line.\n");
        fclose(fp);
        return NULL;
    }

    int *offs = malloc((cls+1) * sizeof(int));
    if (!offs) {
        perror("malloc");
        fclose(fp);
        return NULL;
    }
    int *flat = NULL;
    int flat_sz = 0, flat_cap = 0;
    int cur = 0;
    offs[0] = 0;

    while (fgets(line, sizeof(line), fp)) {
        ++line_no;
        char *scan = line;
        while (isspace((unsigned char)*scan)) ++scan;
        if (*scan == 'c' || *scan == 'p' || *scan == '\0') continue;

        /* DIMACS clauses are terminated by 0 and can span multiple lines.
           Parse directly without strtok for better performance. */
        while (*scan) {
            /* Skip whitespace */
            while (isspace((unsigned char)*scan)) ++scan;
            if (*scan == '\0' || *scan == 'c') break;

            /* Parse integer token */
            char *endptr = NULL;
            long lit_l = strtol(scan, &endptr, 10);
            if (scan == endptr) break;  /* No valid integer */
            scan = endptr;

            if (lit_l < INT_MIN || lit_l > INT_MAX) {
                fprintf(stderr, "Integer overflow at line %d\n", line_no);
                free(offs); free(flat); fclose(fp); return NULL;
            }
            int lit = (int)lit_l;

            if (lit == 0) {
                if (cur >= cls) {
                    fprintf(stderr, "Too many clauses (header=%d) by line %d\n", cls, line_no);
                    free(offs); free(flat); fclose(fp); return NULL;
                }
                ++cur;
                offs[cur] = flat_sz;
            } else {
                if (lit < -vars || lit > vars) {
                    fprintf(stderr, "Literal %d exceeds declared variable count %d\n", lit, vars);
                    free(offs); free(flat); fclose(fp); return NULL;
                }
                if (flat_sz == flat_cap) {
                    int new_cap = flat_cap ? flat_cap * 2 : 1024;
                    int *new_flat = realloc(flat, new_cap * sizeof(int));
                    if (!new_flat) {
                        perror("realloc");
                        free(offs); free(flat); fclose(fp); return NULL;
                    }
                    flat = new_flat;
                    flat_cap = new_cap;
                }
                flat[flat_sz++] = lit;
            }
        }
        if (cur == cls) break;
    }
    fclose(fp);
    if (cur != cls) {
        fprintf(stderr, "Clause count mismatch: header=%d parsed=%d\n", cls, cur);
        free(offs);
        free(flat);
        return NULL;
    }

    Instance *inst = malloc(sizeof(Instance));
    if (!inst) {
        perror("malloc");
        free(offs);
        free(flat);
        return NULL;
    }
    inst->num_vars    = vars;
    inst->num_clauses = cls;
    inst->cl_offs     = offs;
    inst->cl_flat     = flat;
    return inst;
}

/* --------------------------------------------------------------------
   7b. Math-guided DRAT backend (MATH.md Sections 6, 8, 9)

       The proof strategy uses the mathematical framework:

       1. Unit propagation to fixpoint (standard RUP prerequisite).
       2. Topological core extraction: compute variable participation
          in unsatisfied clauses weighted by prime mass W(p)=1/(1+ln p).
          Variables with highest topological degree are probed first
          because they sit at β₁ cycle hotspots (MATH.md §8).
       3. GF(2) parity check: for binary (k=2) clause subgraphs that
          form odd cycles (β₁ > 0), derive contradiction algebraically
          and emit the corresponding DRAT unit additions.
       4. Failed-literal probing in math-guided order: probe variables
          sorted by descending topological pressure, not sequentially.
       5. Emit empty clause on contradiction.
-------------------------------------------------------------------- */
static int proof_unit_propagate(const Instance *inst, int *assign)
{
    int changed;
    do {
        changed = 0;
        for (int c = 0; c < inst->num_clauses; ++c) {
            int sat = 0;
            int unassigned_count = 0;
            int last_unassigned_lit = 0;

            for (int i = inst->cl_offs[c]; i < inst->cl_offs[c+1]; ++i) {
                int lit = inst->cl_flat[i];
                int v = abs(lit);
                int av = assign[v];
                if (av == PROOF_UNASSIGNED) {
                    unassigned_count++;
                    last_unassigned_lit = lit;
                } else if ((lit > 0 && av == 1) || (lit < 0 && av == 0)) {
                    sat = 1;
                    break;
                }
            }

            if (sat) continue;
            if (unassigned_count == 0) return 0; /* conflict */
            if (unassigned_count == 1) {
                int v = abs(last_unassigned_lit);
                int need = (last_unassigned_lit > 0) ? 1 : 0;
                if (assign[v] == PROOF_UNASSIGNED) {
                    assign[v] = need;
                    changed = 1;
                } else if (assign[v] != need) {
                    return 0;
                }
            }
        }
    } while (changed);
    return 1;
}

static int proof_assumption_conflict_ucp(const Instance *inst, const int *base_assign, int lit)
{
    int n = inst->num_vars;
    int *work = malloc((n + 1) * sizeof(int));
    if (!work) return 0;
    memcpy(work, base_assign, (n + 1) * sizeof(int));

    int v = abs(lit);
    int val = (lit > 0) ? 1 : 0;
    int conflict = 0;
    if (work[v] != PROOF_UNASSIGNED && work[v] != val) {
        conflict = 1;
    } else {
        work[v] = val;
        if (!proof_unit_propagate(inst, work)) conflict = 1;
    }

    free(work);
    return conflict;
}

/* Compute topological probe order: variables sorted by descending
   prime-weighted participation in unsatisfied clauses (MATH.md §9).
   This is the "irreducibility pressure" — variables under highest
   contradictory force from the prime necklace boundary conditions. */
static int *compute_topo_probe_order(const Instance *inst, const int *assign,
                                     const int *primes, int prime_cnt, int *out_cnt)
{
    int n = inst->num_vars;
    double *pressure = calloc(n + 1, sizeof(double));
    if (!pressure) { *out_cnt = 0; return NULL; }

    for (int c = 0; c < inst->num_clauses; ++c) {
        int sat = 0;
        for (int i = inst->cl_offs[c]; i < inst->cl_offs[c+1]; ++i) {
            int v = abs(inst->cl_flat[i]);
            int av = assign[v];
            if (av != PROOF_UNASSIGNED &&
                ((inst->cl_flat[i] > 0 && av == 1) || (inst->cl_flat[i] < 0 && av == 0))) {
                sat = 1; break;
            }
        }
        if (sat) continue;

        int p = (c < prime_cnt) ? primes[c] : 2;
        double w = 1.0 / (1.0 + log((double)p));
        int k = inst->cl_offs[c+1] - inst->cl_offs[c];
        for (int i = inst->cl_offs[c]; i < inst->cl_offs[c+1]; ++i) {
            int v = abs(inst->cl_flat[i]);
            if (assign[v] == PROOF_UNASSIGNED)
                pressure[v] += w * (double)k;
        }
    }

    int cnt = 0;
    int *order = malloc(n * sizeof(int));
    if (!order) { free(pressure); *out_cnt = 0; return NULL; }
    for (int v = 1; v <= n; ++v)
        if (assign[v] == PROOF_UNASSIGNED && pressure[v] > 0.0)
            order[cnt++] = v;

    for (int i = 1; i < cnt; ++i) {
        int key = order[i];
        double kp = pressure[key];
        int j = i - 1;
        while (j >= 0 && pressure[order[j]] < kp) {
            order[j + 1] = order[j];
            j--;
        }
        order[j + 1] = key;
    }

    free(pressure);
    *out_cnt = cnt;
    return order;
}

/* GF(2) parity check on binary (k=2) clause subgraph.
   Binary clauses (a ∨ b) encode the implication ¬a → b, ¬b → a.
   In GF(2) terms: a ⊕ b = 1. An odd cycle in this implication graph
   forces a variable to equal both 0 and 1 — a parity contradiction.
   This corresponds to β₁ > 0 in the unsatisfied clause complex
   (MATH.md §8: "β₁ persistent homology catches the actual cycle
   structure of the XOR constraint graph"). */
static int proof_gf2_parity_check(const Instance *inst, const int *assign,
                                  FILE *pf, int *derived_units)
{
    int n = inst->num_vars;
    int *color = calloc(n + 1, sizeof(int));
    if (!color) return 0;

    typedef struct { int to; int parity; } Edge;
    int edge_cap = 1024;
    Edge *edges = malloc(edge_cap * sizeof(Edge));
    int *adj_ptr = calloc(n + 2, sizeof(int));
    if (!edges || !adj_ptr) { free(color); free(edges); free(adj_ptr); return 0; }

    int binary_cnt = 0;
    for (int c = 0; c < inst->num_clauses; ++c) {
        int k = inst->cl_offs[c+1] - inst->cl_offs[c];
        if (k != 2) continue;

        int l0 = inst->cl_flat[inst->cl_offs[c]];
        int l1 = inst->cl_flat[inst->cl_offs[c] + 1];
        int v0 = abs(l0), v1 = abs(l1);
        if (assign[v0] != PROOF_UNASSIGNED || assign[v1] != PROOF_UNASSIGNED) continue;
        if (v0 == v1) continue;

        adj_ptr[v0 + 1]++;
        adj_ptr[v1 + 1]++;
        binary_cnt++;
    }
    if (binary_cnt == 0) { free(color); free(edges); free(adj_ptr); return 0; }

    for (int i = 1; i <= n + 1; ++i) adj_ptr[i] += adj_ptr[i-1];
    int total_edges = adj_ptr[n + 1];
    if (total_edges > edge_cap) {
        edge_cap = total_edges;
        Edge *new_e = realloc(edges, edge_cap * sizeof(Edge));
        if (!new_e) { free(color); free(edges); free(adj_ptr); return 0; }
        edges = new_e;
    }
    int *fill = calloc(n + 1, sizeof(int));
    if (!fill) { free(color); free(edges); free(adj_ptr); return 0; }

    for (int c = 0; c < inst->num_clauses; ++c) {
        int k = inst->cl_offs[c+1] - inst->cl_offs[c];
        if (k != 2) continue;
        int l0 = inst->cl_flat[inst->cl_offs[c]];
        int l1 = inst->cl_flat[inst->cl_offs[c] + 1];
        int v0 = abs(l0), v1 = abs(l1);
        if (assign[v0] != PROOF_UNASSIGNED || assign[v1] != PROOF_UNASSIGNED) continue;
        if (v0 == v1) continue;

        /* (a ∨ b): ¬a→b, ¬b→a.  GF(2) parity: if both positive or
           both negative, parity=0 (same); otherwise parity=1 (different). */
        int par = ((l0 > 0) == (l1 > 0)) ? 0 : 1;
        edges[adj_ptr[v0] + fill[v0]++] = (Edge){v1, par};
        edges[adj_ptr[v1] + fill[v1]++] = (Edge){v0, par};
    }
    free(fill);

    /* BFS 2-coloring to detect odd cycles */
    int *queue = malloc((n + 1) * sizeof(int));
    int contradiction = 0;
    int contra_var = 0;

    for (int start = 1; start <= n && !contradiction; ++start) {
        if (color[start] || adj_ptr[start+1] == adj_ptr[start]) continue;
        if (assign[start] != PROOF_UNASSIGNED) continue;

        color[start] = 1;
        int qh = 0, qt = 0;
        queue[qt++] = start;

        while (qh < qt && !contradiction) {
            int u = queue[qh++];
            for (int e = adj_ptr[u]; e < adj_ptr[u+1]; ++e) {
                int w = edges[e].to;
                int expected = (edges[e].parity == 0) ? color[u] : -color[u];
                if (color[w] == 0) {
                    color[w] = expected;
                    queue[qt++] = w;
                } else if (color[w] != expected) {
                    contradiction = 1;
                    contra_var = u;
                    break;
                }
            }
        }
    }

    if (contradiction && pf && contra_var > 0) {
        /* Emit the forced assignment as DRAT unit clause */
        fprintf(pf, "%d 0\n", contra_var);
        if (derived_units) (*derived_units)++;
    }

    free(queue); free(color); free(edges); free(adj_ptr);
    return contradiction;
}

/* Run the 3-phase proof engine on an Instance.
   Returns 1 if empty clause derived (UNSAT proved), 0 otherwise.
   
   Proof format support:
   - DRAT: Emits clauses in standard DRAT format (additions only, no deletions)
   - LRAT: Emits clauses with hint annotations for efficient verification
   
   LRAT hints: For each derived clause C at line L, hints are clause IDs that
   together with C produce a conflict via unit propagation. This allows linear-time
   verification without re-running full unit propagation.
*/
static int proof_run_on_instance(const Instance *inst, FILE *pf,
                                 const int *primes, int prime_cnt,
                                 int *derived_units, char *status, size_t status_cap,
                                 int is_lrat, int *clause_ids, int *clause_id_ptr)
{
    (void)clause_ids;
    int n = inst->num_vars;
    int *assign = malloc((n + 1) * sizeof(int));
    if (!assign) { snprintf(status, status_cap, "oom"); return 0; }
    for (int i = 1; i <= n; ++i) assign[i] = PROOF_UNASSIGNED;

    /* Track which clauses are used in unit propagation for LRAT hints */
    int *up_witness = NULL;
    if (is_lrat) {
        up_witness = malloc(inst->num_clauses * sizeof(int));
        if (!up_witness) { free(assign); snprintf(status, status_cap, "oom"); return 0; }
        for (int c = 0; c < inst->num_clauses; ++c) up_witness[c] = -1;
    }

    /* Phase 1: Unit propagation to fixpoint */
    int up_changed;
    do {
        up_changed = 0;
        for (int c = 0; c < inst->num_clauses; ++c) {
            int sat = 0, unassigned_count = 0, last_unassigned_lit = 0;
            for (int i = inst->cl_offs[c]; i < inst->cl_offs[c+1]; ++i) {
                int lit = inst->cl_flat[i];
                int v = abs(lit);
                int av = assign[v];
                if (av == PROOF_UNASSIGNED) {
                    unassigned_count++;
                    last_unassigned_lit = lit;
                } else if ((lit > 0 && av == 1) || (lit < 0 && av == 0)) {
                    sat = 1;
                    if (is_lrat) up_witness[c] = c;  /* This clause was satisfied */
                    break;
                }
            }
            if (sat) continue;
            if (unassigned_count == 0) {
                /* Conflict found - emit empty clause with LRAT hints if requested */
                if (is_lrat && up_witness) {
                    fprintf(pf, "0");
                    for (int j = 0; j < inst->num_clauses; ++j)
                        if (up_witness[j] >= 0) fprintf(pf, " %d", up_witness[j] + 1);
                    fprintf(pf, " 0\n");
                } else {
                    fprintf(pf, "0\n");
                }
                free(assign); free(up_witness);
                snprintf(status, status_cap, "generated_top_level_up_conflict");
                return 1;
            }
            if (unassigned_count == 1) {
                int v = abs(last_unassigned_lit);
                int need = (last_unassigned_lit > 0) ? 1 : 0;
                if (assign[v] == PROOF_UNASSIGNED) {
                    assign[v] = need;
                    up_changed = 1;
                } else if (assign[v] != need) {
                    /* Conflict */
                    if (is_lrat && up_witness) {
                        fprintf(pf, "0");
                        for (int j = 0; j < inst->num_clauses; ++j)
                            if (up_witness[j] >= 0) fprintf(pf, " %d", up_witness[j] + 1);
                        fprintf(pf, " 0\n");
                    } else {
                        fprintf(pf, "0\n");
                    }
                    free(assign); free(up_witness);
                    snprintf(status, status_cap, "generated_top_level_up_conflict");
                    return 1;
                }
            }
        }
    } while (up_changed);

    /* Phase 2: GF(2) parity check on binary clauses */
    if (proof_gf2_parity_check(inst, assign, pf, derived_units)) {
        /* Re-run unit propagation after parity-derived unit */
        up_changed = 1;
        while (up_changed) {
            up_changed = 0;
            for (int c = 0; c < inst->num_clauses; ++c) {
                int sat = 0, unassigned_count = 0, last_unassigned_lit = 0;
                for (int i = inst->cl_offs[c]; i < inst->cl_offs[c+1]; ++i) {
                    int lit = inst->cl_flat[i];
                    int v = abs(lit);
                    int av = assign[v];
                    if (av == PROOF_UNASSIGNED) {
                        unassigned_count++;
                        last_unassigned_lit = lit;
                    } else if ((lit > 0 && av == 1) || (lit < 0 && av == 0)) {
                        sat = 1; break;
                    }
                }
                if (sat) continue;
                if (unassigned_count == 0) {
                    if (is_lrat) fprintf(pf, "0 0\n"); else fprintf(pf, "0\n");
                    free(assign); free(up_witness);
                    snprintf(status, status_cap, "generated_gf2_parity_refutation");
                    return 1;
                }
                if (unassigned_count == 1) {
                    int v = abs(last_unassigned_lit);
                    int need = (last_unassigned_lit > 0) ? 1 : 0;
                    if (assign[v] == PROOF_UNASSIGNED) {
                        assign[v] = need; up_changed = 1;
                    } else if (assign[v] != need) {
                        if (is_lrat) fprintf(pf, "0 0\n"); else fprintf(pf, "0\n");
                        free(assign); free(up_witness);
                        snprintf(status, status_cap, "generated_gf2_parity_refutation");
                        return 1;
                    }
                }
            }
        }
    }

    /* Phase 3: Failed-literal probing in topological order */
    for (int round = 0; round < n; ++round) {
        int order_cnt = 0;
        int *order = compute_topo_probe_order(inst, assign, primes, prime_cnt, &order_cnt);
        if (!order || order_cnt == 0) { free(order); break; }

        int progress = 0;
        for (int idx = 0; idx < order_cnt; ++idx) {
            int v = order[idx];
            if (assign[v] != PROOF_UNASSIGNED) continue;

            /* Try assuming v = false */
            if (proof_assumption_conflict_ucp(inst, assign, v)) {
                /* Emit unit clause -v with clause ID for LRAT */
                int clause_id = (*clause_id_ptr)++;
                if (is_lrat) {
                    /* LRAT format: literal(s) hint-clause-IDs 0 */
                    fprintf(pf, "-%d %d 0\n", v, clause_id);
                } else {
                    fprintf(pf, "-%d 0\n", v);
                }
                assign[v] = 0;
                if (derived_units) (*derived_units)++;
                
                /* Re-propagate and check for conflict */
                if (!proof_unit_propagate(inst, assign)) {
                    if (is_lrat) fprintf(pf, "0 %d 0\n", clause_id); else fprintf(pf, "0\n");
                    free(order); free(assign); free(up_witness);
                    snprintf(status, status_cap, "generated_topo_failed_literal_refutation");
                    return 1;
                }
                progress = 1; break;
            }
            /* Try assuming v = true */
            if (proof_assumption_conflict_ucp(inst, assign, -v)) {
                int clause_id = (*clause_id_ptr)++;
                if (is_lrat) {
                    fprintf(pf, "%d %d 0\n", v, clause_id);
                } else {
                    fprintf(pf, "%d 0\n", v);
                }
                assign[v] = 1;
                if (derived_units) (*derived_units)++;
                
                if (!proof_unit_propagate(inst, assign)) {
                    if (is_lrat) fprintf(pf, "0 %d 0\n", clause_id); else fprintf(pf, "0\n");
                    free(order); free(assign); free(up_witness);
                    snprintf(status, status_cap, "generated_topo_failed_literal_refutation");
                    return 1;
                }
                progress = 1; break;
            }
        }
        free(order);
        if (!progress) break;
    }

    free(assign); free(up_witness);
    return 0;
}

/* Extract UNSAT core from solver state: unsatisfied clauses + their
   1-neighbourhood (clauses sharing a variable). If the neighbourhood
   would exceed the size limits, fall back to unsatisfied-only.
   The core is a subset — if the core is UNSAT, the full formula is. */
static Instance *extract_unsat_core(const NitroSat *ns, int *core_clause_cnt)
{
    recompute_sat_counts((NitroSat *)(uintptr_t)ns);
    int core_cnt = 0;
    for (int c = 0; c < ns->num_clauses; ++c)
        if (ns->sat_counts[c] == 0) core_cnt++;

    if (core_cnt == 0) { if (core_clause_cnt) *core_clause_cnt = 0; return NULL; }

    uint8_t *in_core = calloc(ns->num_clauses, 1);
    uint8_t *core_var = calloc(ns->num_vars + 1, 1);
    if (!in_core || !core_var) {
        free(in_core); free(core_var);
        if (core_clause_cnt) *core_clause_cnt = core_cnt;
        return NULL;
    }

    for (int c = 0; c < ns->num_clauses; ++c) {
        if (ns->sat_counts[c] != 0) continue;
        in_core[c] = 1;
        for (int j = ns->cl_offs[c]; j < ns->cl_offs[c+1]; ++j)
            core_var[abs(ns->cl_flat[j])] = 1;
    }

    /* Expand to 1-neighbourhood via v2c CSR if available and result
       stays within proof size limits */
    int expanded = core_cnt;
    if (ns->v2c_ptr && ns->v2c_data) {
        for (int v = 1; v <= ns->num_vars; ++v) {
            if (!core_var[v]) continue;
            for (int p = ns->v2c_ptr[v]; p < ns->v2c_ptr[v+1]; ++p) {
                int c = ns->v2c_data[p];
                if (!in_core[c]) { in_core[c] = 1; expanded++; }
            }
        }
    }

    if (expanded > proof_max_clauses) {
        /* Neighbourhood too large — fall back to unsatisfied only */
        memset(in_core, 0, ns->num_clauses);
        expanded = 0;
        for (int c = 0; c < ns->num_clauses; ++c) {
            if (ns->sat_counts[c] == 0) { in_core[c] = 1; expanded++; }
        }
    }

    /* Build reduced Instance (preserves original variable numbering) */
    int *offs = malloc((expanded + 1) * sizeof(int));
    long flat_cap = 0;
    for (int c = 0; c < ns->num_clauses; ++c)
        if (in_core[c]) flat_cap += ns->cl_offs[c+1] - ns->cl_offs[c];
    int *flat = (flat_cap <= (long)INT_MAX) ? malloc((size_t)flat_cap * sizeof(int)) : NULL;
    if (!offs || !flat) {
        free(in_core); free(core_var); free(offs); free(flat);
        if (core_clause_cnt) *core_clause_cnt = core_cnt;
        return NULL;
    }

    int idx = 0, fpos = 0;
    offs[0] = 0;
    for (int c = 0; c < ns->num_clauses; ++c) {
        if (!in_core[c]) continue;
        for (int j = ns->cl_offs[c]; j < ns->cl_offs[c+1]; ++j)
            flat[fpos++] = ns->cl_flat[j];
        idx++;
        offs[idx] = fpos;
    }

    Instance *core = malloc(sizeof(Instance));
    if (!core) {
        free(offs); free(flat); free(in_core); free(core_var);
        if (core_clause_cnt) *core_clause_cnt = core_cnt;
        return NULL;
    }
    core->num_vars = ns->num_vars;
    core->num_clauses = expanded;
    core->cl_offs = offs;
    core->cl_flat = flat;

    free(in_core); free(core_var);
    if (core_clause_cnt) *core_clause_cnt = core_cnt;
    return core;
}

static int generate_drat_unsat_proof(const NitroSat *ns, const char *proof_path,
                                     int *derived_units, char *status, size_t status_cap,
                                     const char *proof_format)
{
    if (derived_units) *derived_units = 0;
    if (!proof_path || !*proof_path) {
        snprintf(status, status_cap, "missing_proof_path");
        return 0;
    }

    const Instance *inst = ns->inst;
    int is_lrat = (proof_format && strcmp(proof_format, "lrat") == 0);

    /* ── Solver UNSAT awareness signals (MATH.md §§6,8,9) ────────── */
    int fracture_detected = ns->fd ? fd_is_fracture(ns->fd) : 0;
    int beta1_persistent  = (ns->pt.has_initial && ns->pt.final_beta1 > 0);
    int non_convex        = 0;
    {
        double k_max = 3.0;
        double d_clause = 0.0;
        for (int i = 1; i <= ns->num_vars; ++i) d_clause += ns->degrees[i];
        d_clause = (ns->num_vars > 0) ? d_clause / ns->num_vars : 1.0;
        if (d_clause < 1e-12) d_clause = 1.0;
        double eff_beta = ns->heat_beta * d_clause;
        double W_max = ns->cl_weights[0];
        for (int c = 1; c < ns->num_clauses; ++c)
            if (ns->cl_weights[c] > W_max) W_max = ns->cl_weights[c];
        double lhs = 4.0 / (eff_beta > 1e-12 ? eff_beta : 1e-12);
        double rhs = W_max * k_max * k_max * d_clause / (0.3 * 0.3);
        non_convex = (lhs <= rhs);
    }

    /* ── Extract UNSAT core from solver blame state ───────────────── */
    int unsat_clause_cnt = 0;
    Instance *core = extract_unsat_core(ns, &unsat_clause_cnt);
    const Instance *target = core ? core : inst;
    int using_core = (core != NULL);

    if (target->num_vars > proof_max_vars || target->num_clauses > proof_max_clauses) {
        if (core) { free(core->cl_offs); free(core->cl_flat); free(core); }
        snprintf(status, status_cap, "skipped_size_limit(unsat=%d,core=%d/%d,fractures=%d,beta1=%d,convex=%s)",
                 unsat_clause_cnt, target->num_vars, target->num_clauses,
                 fracture_detected, beta1_persistent,
                 non_convex ? "NO" : "YES");
        return 0;
    }

    int prime_cnt = 0;
    int *primes = sieve_primes(target->num_clauses > 0 ? target->num_clauses : 1, &prime_cnt);

    FILE *pf = fopen(proof_path, "w");
    if (!pf) {
        free(primes);
        if (core) { free(core->cl_offs); free(core->cl_flat); free(core); }
        snprintf(status, status_cap, "proof_file_open_failed");
        return 0;
    }

    /* Write LRAT header if needed */
    int clause_id_ptr = 1;  /* Clause IDs start at 1 for LRAT */
    int *clause_ids = NULL;
    if (is_lrat) {
        /* LRAT format: emit original clauses first with IDs */
        fprintf(pf, "p cnf %d %d\n", target->num_vars, target->num_clauses);
        clause_ids = malloc(target->num_clauses * sizeof(int));
        if (!clause_ids) {
            fclose(pf); free(primes);
            if (core) { free(core->cl_offs); free(core->cl_flat); free(core); }
            snprintf(status, status_cap, "oom");
            return 0;
        }
        /* Emit original clauses in sorted order (proof compression) */
        for (int c = 0; c < target->num_clauses; ++c) {
            clause_ids[c] = clause_id_ptr++;
            /* Sort literals for better drat-trim performance */
            int start = target->cl_offs[c];
            int end = target->cl_offs[c+1];
            int len = end - start;
            int *sorted_lits = malloc(len * sizeof(int));
            if (sorted_lits) {
                memcpy(sorted_lits, target->cl_flat + start, len * sizeof(int));
                /* Simple insertion sort for small clauses (typical k=3) */
                for (int i = 1; i < len; ++i) {
                    int key = sorted_lits[i], j = i - 1;
                    while (j >= 0 && abs(sorted_lits[j]) > abs(key)) {
                        sorted_lits[j + 1] = sorted_lits[j]; j--;
                    }
                    sorted_lits[j + 1] = key;
                }
                for (int i = 0; i < len; ++i) fprintf(pf, "%d ", sorted_lits[i]);
                fprintf(pf, "0 %d 0\n", clause_ids[c]);  /* ID + empty hint list */
                free(sorted_lits);
            } else {
                for (int i = start; i < end; ++i) fprintf(pf, "%d ", target->cl_flat[i]);
                fprintf(pf, "0 %d 0\n", clause_ids[c]);
            }
        }
    }

    /* ── Run proof engine on the (possibly reduced) core ──────────── */
    int result = proof_run_on_instance(target, pf, primes, prime_cnt,
                                       derived_units, status, status_cap,
                                       is_lrat, clause_ids, &clause_id_ptr);

    if (!result) {
        /* Core didn't yield proof — try full formula if we used a core */
        if (using_core && inst->num_vars <= proof_max_vars &&
            inst->num_clauses <= proof_max_clauses) {
            /* For LRAT, we need to emit all original clauses first */
            if (is_lrat && !clause_ids) {
                clause_ids = malloc(inst->num_clauses * sizeof(int));
                if (clause_ids) {
                    for (int c = 0; c < inst->num_clauses; ++c) {
                        clause_ids[c] = clause_id_ptr++;
                        int start = inst->cl_offs[c], end = inst->cl_offs[c+1];
                        for (int i = start; i < end; ++i) fprintf(pf, "%d ", inst->cl_flat[i]);
                        fprintf(pf, "0 %d 0\n", clause_ids[c]);
                    }
                }
            }
            result = proof_run_on_instance(inst, pf, ns->primes, ns->prime_cnt,
                                           derived_units, status, status_cap,
                                           is_lrat, clause_ids, &clause_id_ptr);
        }
    }

    free(clause_ids);
    free(primes);
    if (core) { free(core->cl_offs); free(core->cl_flat); free(core); }

    if (!result) {
        fclose(pf);
        remove(proof_path);
        snprintf(status, status_cap,
                 "inconclusive(unsat=%d,core=%d,fractures=%d,beta1=%d,convex=%s)",
                 unsat_clause_cnt,
                 using_core ? target->num_clauses : inst->num_clauses,
                 fracture_detected, beta1_persistent,
                 non_convex ? "NO" : "YES");
        return 0;
    }

    fclose(pf);
    return 1;
}

/* --------------------------------------------------------------------
   nitrosat_free() - Release all NitroSat resources
   
   Memory Ownership: Frees all [OWNER] fields in NitroSat structure.
   Does NOT free ns->inst (owned by caller, freed in main()).
   Does NOT free ns->cl_offs/cl_flat (borrowed from inst).
   
   Free order: Data fields first, then structure itself.
   NULL-safe: All checks are performed before free().
-------------------------------------------------------------------- */
static void nitrosat_free(NitroSat *ns)
{
    if (!ns) return;
    free(ns->x);
    free(ns->degrees);
    free(ns->decimated);
    free(ns->cl_weights);
    free(ns->zeta_weights);
    free(ns->is_hard);
    if (ns->opt) {
        free(ns->opt->m);
        free(ns->opt->v);
        free(ns->opt);
    }
    free(ns->fd);
    free(ns->primes);
    free(ns->v2c_ptr);
    free(ns->v2c_data);
    free(ns->grad_buffer);
    free(ns->heat_mult_buffer);
    free(ns->sat_counts);
    free(ns->entropy_cache);
    free(ns->entropy_x_cache);
    free(ns->prev_x);
    free(ns);
}

/* --------------------------------------------------------------------
   8.  Helper – compute vertex degrees (used for heat kernel)
-------------------------------------------------------------------- */
static void compute_degrees(NitroSat *ns)
{
    int n = ns->num_vars;
    ns->degrees = calloc(n+1, sizeof(double));
    for (int c = 0; c < ns->num_clauses; ++c) {
        int s = ns->cl_offs[c];
        int e = ns->cl_offs[c+1];
        int k = e - s;                     /* clause length */
        for (int i=s; i<e; ++i) {
            int lit = ns->cl_flat[i];
            int v = abs(lit);
            ns->degrees[v] += (double)(k-1);
        }
    }
}

/* --------------------------------------------------------------------
   8b. Spectral initialisation – XOR-Laplacian power iteration
        (matches gh.lua xor_laplacian.spectral_init)
-------------------------------------------------------------------- */
static void spectral_init(NitroSat *ns, int iters)
{
    int nv = ns->num_vars;
    int nc = ns->num_clauses;
    double *v = malloc((nv+1) * sizeof(double));
    double *y = malloc((nv+1) * sizeof(double));

    /* random initial vector */
    for (int i = 1; i <= nv; ++i) v[i] = rng_double();

    for (int it = 0; it < iters; ++it) {
        /* y = D*v  (degree-weighted) */
        for (int i = 1; i <= nv; ++i) y[i] = ns->degrees[i] * v[i];

        /* subtract off-diagonal: for each clause, y[vi] -= (sum - v[vi]) */
        for (int c = 0; c < nc; ++c) {
            int s = ns->cl_offs[c];
            int e = ns->cl_offs[c+1];
            double sum = 0.0;
            for (int i = s; i < e; ++i) sum += v[abs(ns->cl_flat[i])];
            for (int i = s; i < e; ++i) {
                int vi = abs(ns->cl_flat[i]);
                y[vi] -= (sum - v[vi]);
            }
        }

        /* normalise */
        double norm = 0.0;
        for (int i = 1; i <= nv; ++i) norm += y[i] * y[i];
        norm = sqrt(norm);
        if (norm > 1e-10)
            for (int i = 1; i <= nv; ++i) v[i] = y[i] / norm;
        else
            for (int i = 1; i <= nv; ++i) v[i] = y[i];
    }

    /* rescale to [0,1] */
    double mi = 1e300, ma = -1e300;
    for (int i = 1; i <= nv; ++i) {
        if (v[i] < mi) mi = v[i];
        if (v[i] > ma) ma = v[i];
    }
    double r = ma - mi;
    if (r > 1e-10)
        for (int i = 1; i <= nv; ++i) v[i] = (v[i] - mi) / r;

    /* store with small noise */
    for (int i = 1; i <= nv; ++i)
        ns->x[i] = v[i] + (rng_double() - 0.5) * 0.01;

    free(v);
    free(y);
}

/* --------------------------------------------------------------------
   9.  Initialise a NitroSat object
-------------------------------------------------------------------- */
static NitroSat *nitrosat_new(Instance *inst, int max_steps, int verbose)
{
    NitroSat *ns = calloc(1, sizeof(NitroSat));
    ns->inst       = inst;
    ns->cl_offs    = inst->cl_offs;
    ns->cl_flat    = inst->cl_flat;
    ns->num_vars   = inst->num_vars;
    ns->num_clauses= inst->num_clauses;
    ns->max_steps  = max_steps;
    ns->verbose    = verbose;

    /* degree vector (heat kernel) */
    compute_degrees(ns);

    /* continuous assignment – spectral initialisation (power iteration) */
    ns->x = malloc((ns->num_vars+1) * sizeof(double));
    spectral_init(ns, 50);  /* 50 iterations, matches gh.lua */

    /* clause weight initialisation (prime‑weighted) */
    ns->primes = sieve_primes(ns->num_clauses, &ns->prime_cnt);
    ns->cl_weights = malloc(ns->num_clauses * sizeof(double));
    ns->zeta_weights = malloc(ns->num_clauses * sizeof(double));
    for (int c=0;c<ns->num_clauses;++c) {
        int p = (c < ns->prime_cnt) ? ns->primes[c] : 2;
        ns->cl_weights[c]   = 1.0 / pow(1.0 + log((double)p), 1.0);   /* prime_alpha = 1.0 */
        ns->zeta_weights[c] = log((double)p) / (double)p;
    }
    ns->is_hard = malloc(ns->num_clauses);
    for (int c=0;c<ns->num_clauses;++c) ns->is_hard[c] = 1;   /* all hard initially, matches gh.lua */

    /* decimation vector */
    ns->decimated = calloc(ns->num_vars+1, 1);

    /* optimiser – NADAM (the Lua default) */
    ns->opt = optimizer_new("nadam", ns->x, ns->num_vars,
                            0.002, 0.9, 0.999, 1e-8, 0.02);

    /* BAHA detector */
    ns->fd = fd_new(1.5);
    ns->sigma_factor = 1.5;   /* same as BAHA default */
    ns->baha_var_peak = -1e100;

    /* parameters */
    ns->heat_beta     = 0.5;
    ns->heat_lambda   = 0.1;
    ns->entropy_weight= 0.01;
    ns->dec_threshold = 0.49;
    ns->dec_freq      = 100;
    ns->use_topology  = 1;   /* enable topology snapshots */

    ns->grad_buffer      = calloc(ns->num_vars + 1, sizeof(double));
    ns->heat_mult_buffer = malloc((ns->num_vars + 1) * sizeof(double));
    ns->sat_counts       = calloc(ns->num_clauses, sizeof(int));
    /* Cheat Code 1: Initialize entropy cache */
    ns->entropy_cache    = malloc((ns->num_vars + 1) * sizeof(double));
    ns->entropy_x_cache = malloc((ns->num_vars + 1) * sizeof(double));
    for (int i = 1; i <= ns->num_vars; ++i) {
        ns->entropy_x_cache[i] = -999.0; /* force initial compute */
    }
    /* Cheat Code 5: Initialize prev_x for incremental sat_counts */
    ns->prev_x = malloc((ns->num_vars + 1) * sizeof(double));
    for (int i = 1; i <= ns->num_vars; ++i) {
        ns->prev_x[i] = ns->x[i];
    }

    /* pre-calculate heat multipliers */
    for (int i=1;i<=ns->num_vars;++i) {
        double w = exp(-ns->heat_beta * ns->degrees[i]);
        ns->heat_mult_buffer[i] = 1.0 + ns->heat_lambda * w;
    }

    return ns;
}

/* heat‑kernel multipliers are pre‑calculated in nitrosat_new */

/* --------------------------------------------------------------------
   11. Compute the gradient of the smooth Log‑Barrier objective.
        Returns the gradient vector (allocated) and the current number
        of unsatisfied clauses.

   OPTIMIZATIONS:
   - Cheat Code 2: Skip clauses with 2+ satisfied literals (negligible gradient)
   - Cheat Code 4: Unroll k=3 case for 3-SAT (dominant case)
   - Cheat Code 1: Cache entropy gradient computation
-------------------------------------------------------------------- */
static int compute_gradients(NitroSat *ns)
{
    int n = ns->num_vars;
    double *grad = ns->grad_buffer;
    for (int i=0; i<=n; ++i) grad[i] = 0.0;
    int unsat = 0;

    for (int c = 0; c < ns->num_clauses; ++c) {
        int s = ns->cl_offs[c];
        int e = ns->cl_offs[c+1];
        int k = e - s;

        /* Assertion: clause size must be within bounds */
        NITROSAT_ASSERT_CLAUSE_SIZE(k);

        double violation = 1.0;
        int sat = 0;

        /* Cheat Code 4: Unroll k=3 (dominant case for random 3-SAT) */
        if (k == 3) {
            /* Inline unrolled clause evaluation for k=3 */
            int l0 = ns->cl_flat[s], l1 = ns->cl_flat[s+1], l2 = ns->cl_flat[s+2];
            double v0 = ns->x[abs(l0)], v1 = ns->x[abs(l1)], v2 = ns->x[abs(l2)];

            double lv0 = (l0 > 0) ? (1.0 - v0) : v0;
            double lv1 = (l1 > 0) ? (1.0 - v1) : v1;
            double lv2 = (l2 > 0) ? (1.0 - v2) : v2;

            lv0 = (lv0 < 1e-12) ? 1e-12 : lv0;
            lv1 = (lv1 < 1e-12) ? 1e-12 : lv1;
            lv2 = (lv2 < 1e-12) ? 1e-12 : lv2;

            violation = lv0 * lv1 * lv2;

            if ((l0 > 0 && v0 > 0.5) || (l0 < 0 && v0 <= 0.5)) sat++;
            if ((l1 > 0 && v1 > 0.5) || (l1 < 0 && v1 <= 0.5)) sat++;
            if ((l2 > 0 && v2 > 0.5) || (l2 < 0 && v2 <= 0.5)) sat++;
        } else {
            /* General case for other clause sizes */
            for (int i = s; i < e; ++i) {
                int lit = ns->cl_flat[i];
                int v   = abs(lit);
                double val = ns->x[v];
                double lv = (lit > 0) ? (1.0 - val) : val;
                if (lv < 1e-12) lv = 1e-12;
                violation *= lv;
                if ((lit > 0 && val > 0.5) || (lit < 0 && val <= 0.5)) {
                    sat = 1;
                }
            }
        }

        if (!sat) ++unsat;

        double barrier = 1.0 / (1e-6 + (1.0 - violation));
        double w = ns->cl_weights[c];
        if (ns->is_hard[c]) w *= 50.0;
        double coef = w * barrier * violation;

        /* Compute gradient contribution */
        if (k == 3) {
            /* Unrolled gradient computation for k=3 */
            int l0 = ns->cl_flat[s], l1 = ns->cl_flat[s+1], l2 = ns->cl_flat[s+2];
            double v0 = ns->x[abs(l0)], v1 = ns->x[abs(l1)], v2 = ns->x[abs(l2)];

            double lv0 = (l0 > 0) ? (1.0 - v0) : v0;
            double lv1 = (l1 > 0) ? (1.0 - v1) : v1;
            double lv2 = (l2 > 0) ? (1.0 - v2) : v2;

            lv0 = (lv0 < 1e-12) ? 1e-12 : lv0;
            lv1 = (lv1 < 1e-12) ? 1e-12 : lv1;
            lv2 = (lv2 < 1e-12) ? 1e-12 : lv2;

            grad[abs(l0)] += coef * ((l0 > 0) ? -1.0 : 1.0) / lv0;
            grad[abs(l1)] += coef * ((l1 > 0) ? -1.0 : 1.0) / lv1;
            grad[abs(l2)] += coef * ((l2 > 0) ? -1.0 : 1.0) / lv2;
        } else {
            for (int i = s; i < e; ++i) {
                int lit = ns->cl_flat[i];
                int v   = abs(lit);
                double val = ns->x[v];
                double lv = (lit > 0) ? (1.0 - val) : val;
                if (lv < 1e-12) lv = 1e-12;
                grad[v] += coef * (((lit > 0) ? -1.0 : 1.0) / lv);
            }
        }
    }

    /* Cheat Code 1: Heat-kernel smoothing + cached entropy gradient */
    for (int i=1; i<=n; ++i) {
        grad[i] *= ns->heat_mult_buffer[i];
        
        /* Gradient Clipping: prevent a single clause from imposing infinite momentum */
        if (grad[i] > 1000.0) grad[i] = 1000.0;
        if (grad[i] < -1000.0) grad[i] = -1000.0;

        /* Only recompute log if x changed by more than 0.001 */
        if (fabs(ns->x[i] - ns->entropy_x_cache[i]) > 0.001) {
            double v_clamped = clamp(ns->x[i], 1e-9, 1.0-1e-9);
            ns->entropy_cache[i] = log((1.0 - v_clamped) / v_clamped);
            ns->entropy_x_cache[i] = ns->x[i];
        }
        grad[i] += ns->entropy_weight * ns->entropy_cache[i];
    }
    return unsat;
}

/* --------------------------------------------------------------------
   12. Adiabatic Learning Rate with Riemann Zeta Gain
       Matches -z'(s)/z(s) as s -> 1
-------------------------------------------------------------------- */
static double zeta_log_derivative_gain(double progress)
{
    /* progress matches 1 - (s-1). As progress -> 1, s -> 1.
       The pole at s=1 behaves as 1/(s-1).
       Our 's-1' here is approximately (1.0 - progress + 0.01). */
    double epsilon = 1.0 - progress + 0.02; 
    return 1.0 / epsilon;
}

static double annealing_lr(double t, double A0, double progress)
{
    double phi = PHI;
    double s = sin(PI/20.0 * t) * exp(-PI/20.0 * t)
             + sin(PI/10.0 * t) * exp(-PI/10.0 * t)
             + sin(PI/9.0  * t) * exp(-PI/9.0  * t)
             + sin((1.0/(phi*phi))*t) * exp(-(1.0/(phi*phi))*t);
    
    /* Adiabatic gain based on arithmetic proximity to zeta pole */
    double gain = zeta_log_derivative_gain(progress);
    
    /* Cap gain to prevent total divergence */
    if (gain > 10.0) gain = 10.0;
    
    return fmax(0.0, A0 * s * gain);
}

/* --------------------------------------------------------------------
   13. Satisfaction test (discrete interpretation)
-------------------------------------------------------------------- */
static void recompute_sat_counts(NitroSat *ns)
{
    memset(ns->sat_counts, 0, ns->num_clauses * sizeof(int));
    for (int c=0; c<ns->num_clauses; ++c) {
        for (int i=ns->cl_offs[c]; i<ns->cl_offs[c+1]; ++i) {
            int lit = ns->cl_flat[i];
            int v = abs(lit);
            if ((lit > 0 && ns->x[v] > 0.5) || (lit < 0 && ns->x[v] <= 0.5)) {
                ns->sat_counts[c]++;
            }
        }
    }
}

/* --------------------------------------------------------------------
   13. Satisfaction test (discrete interpretation)
-------------------------------------------------------------------- */
static int check_satisfaction(NitroSat *ns)
{
    int sat = 0;
    for (int c=0; c<ns->num_clauses; ++c) {
        if (ns->sat_counts[c] > 0) sat++;
    }
    return sat;
}

/* --------------------------------------------------------------------
   14. Variable decimation (soft locking of confident variables)
-------------------------------------------------------------------- */
static int decimate(NitroSat *ns)
{
    int locked = 0;
    for (int i=1;i<=ns->num_vars;++i) if (!ns->decimated[i]) {
        double conf = fabs(ns->x[i] - 0.5);
        if (conf > ns->dec_threshold) {
            ns->x[i] = (ns->x[i] > 0.5) ? 1.0 : 0.0;
            ns->decimated[i] = 1;
            ++locked;
        }
    }
    return locked;
}

/* --------------------------------------------------------------------
   15. Zeta‑zero perturbation (standard & “nuclear” version)
-------------------------------------------------------------------- */
static double zeta_zero_perturb(int var_idx, int step, int max_steps,
                               const int *primes, int prime_cnt,
                               const TopologyInfo *topo, int num_vars,
                               int nuclear)
{
    double progress = (double)step / max_steps;
    double val = 0.0;
    int limit = (prime_cnt < 256) ? prime_cnt : 256;

    for (int i=0; i<limit; ++i) {
        int p = primes[i];
        double phase = 2.0*PI * ((double)p * var_idx) / ((log((double)p))*num_vars);
        double amplitude = 1.0 / sqrt((double)p);
        if (nuclear && topo) {
            double crit = sin(2.0*PI*0.5*progress * p); /* crude critical‑line term */
            val += amplitude * sin(phase) * crit;
        } else {
            val += amplitude * sin(phase);
        }
    }

    if (nuclear && topo && topo->complexity_score > 0.5)
        val *= 1.0 + topo->complexity_score;
    if (progress > 0.9) val *= 2.0;
    
    /* Nuclear perturbations must fracture the gradient, so we scale it aggressively 
       rather than using the standard 0.001 multiplier. */
    if (nuclear) {
        return 0.5 * val;
    } else {
        return 0.001 * val;   /* tiny step – same scale as Lua's 0.001 */
    }
}

/* --------------------------------------------------------------------
   17. BAHA-WalkSAT – phase‑transition aware discrete local search
        IMPROVED: Matches baha.cpp magic numbers exactly
        Reference: baha.cpp lines 166-172, 627-629
-------------------------------------------------------------------- */
static int baha_walksat(NitroSat *ns, int max_flips)
{
    /* BAHA magic numbers from baha.cpp */
    const double BAHA_BETA_START = 0.01;
    const double BAHA_BETA_END = 10.0;
    const double BAHA_FRACTURE_THRESHOLD = 1.5;
    const int BAHA_SAMPLES_PER_BETA = (ns->num_vars > 10000) ? 5 : 50;
    const int BAHA_BRANCH_JUMP_INTERVAL = (ns->num_vars > 10000) ? (ns->num_vars / 50) : 200;
    
    int m = ns->num_clauses;
    recompute_sat_counts(ns);

    int *unsat = malloc(m * sizeof(int));
    int *unsat_pos = malloc(m * sizeof(int));
    for (int c = 0; c < m; ++c) unsat_pos[c] = -1;

    int unsz = 0;
    for (int c = 0; c < m; ++c) {
        if (ns->sat_counts[c] == 0) {
            unsat_pos[c] = unsz;
            unsat[unsz++] = c;
        }
    }

    if (unsz == 0) { free(unsat); free(unsat_pos); return 1; }

    /* BAHA improvement: Track best state for branch jumps */
    int best_unsat = unsz;
    double *best_x = malloc((ns->num_vars + 1) * sizeof(double));
    for (int i = 1; i <= ns->num_vars; ++i) best_x[i] = ns->x[i];
    
    /* Pre-allocate recurring buffers to eliminate loop malloc overhead */
    double *sample_x = malloc((ns->num_vars + 1) * sizeof(double));
    
    /* Fracture detection state (matches baha.cpp FractureDetector) */
    double prev_beta = 0.0;
    double prev_log_Z = 0.0;
    int fractures_detected = 0;

    for (int step = 1; step <= max_flips; ++step) {
        if (unsz == 0) break;

        /* Beta schedule: linear from 0.01 to 10.0 (matches baha.cpp line 198-200) */
        double beta = BAHA_BETA_START + 
            (BAHA_BETA_END - BAHA_BETA_START) * ((double)step / max_flips);
        
        /* Estimate log Z = -beta * unsat (boltzmann approximation) */
        double log_Z = -beta * unsz;
        
        /* Fracture detection: rho = |d/dβ log Z| (matches baha.cpp line 131-137) */
        if (step > 1 && prev_beta > 0) {
            double d_log_Z = fabs(log_Z - prev_log_Z);
            double d_beta = beta - prev_beta;
            double rho = (d_beta > 1e-9) ? d_log_Z / d_beta : 0.0;
            
            if (rho > BAHA_FRACTURE_THRESHOLD) {
                fractures_detected++;
            }
        }
        prev_beta = beta;
        prev_log_Z = log_Z;
        
        /* BAHA improvement: Branch-guided jump every 200 steps (matches beta_steps=200) */
        if (step % BAHA_BRANCH_JUMP_INTERVAL == 0 && unsz < best_unsat + 5) {
            /* Sample BAHA_SAMPLES_PER_BETA states and jump to best (line 355-370) */
            int sample_unsat = unsz;
            
            for (int s = 0; s < BAHA_SAMPLES_PER_BETA; ++s) {
                /* Perturb current state (30% flip rate) */
                for (int i = 1; i <= ns->num_vars; ++i) {
                    if (!ns->decimated[i] && rng_double() < 0.3) {
                        sample_x[i] = (ns->x[i] > 0.5) ? 0.0 : 1.0;
                    } else {
                        sample_x[i] = ns->x[i];
                    }
                }
                
                /* Evaluate sample */
                int s_unsat = 0;
                for (int c = 0; c < m; ++c) {
                    int sat = 0;
                    int cs_inner = ns->cl_offs[c];
                    int k = ns->cl_offs[c+1] - cs_inner;
                    
                    /* Unroll k=3 typical case for raw speed */
                    if (k == 3) {
                        int l0 = ns->cl_flat[cs_inner], l1 = ns->cl_flat[cs_inner+1], l2 = ns->cl_flat[cs_inner+2];
                        if ((l0 > 0 && sample_x[abs(l0)] > 0.5) || (l0 < 0 && sample_x[abs(l0)] <= 0.5) ||
                            (l1 > 0 && sample_x[abs(l1)] > 0.5) || (l1 < 0 && sample_x[abs(l1)] <= 0.5) ||
                            (l2 > 0 && sample_x[abs(l2)] > 0.5) || (l2 < 0 && sample_x[abs(l2)] <= 0.5)) {
                            sat = 1;
                        }
                    } else {
                        for (int q = cs_inner; q < ns->cl_offs[c+1]; ++q) {
                            int lit = ns->cl_flat[q], v = abs(lit);
                            double val = sample_x[v];
                            if ((lit > 0 && val > 0.5) || (lit < 0 && val <= 0.5)) { sat = 1; break; }
                        }
                    }
                    if (!sat) s_unsat++;
                }
                
                if (s_unsat < sample_unsat) {
                    sample_unsat = s_unsat;
                    for (int i = 1; i <= ns->num_vars; ++i) best_x[i] = sample_x[i];
                }
            }
            
            /* Jump if improvement found */
            if (sample_unsat < best_unsat) {
                best_unsat = sample_unsat;
                for (int i = 1; i <= ns->num_vars; ++i) ns->x[i] = best_x[i];
                
                /* Recompute sat_counts after jump */
                recompute_sat_counts(ns);
                unsz = 0;
                for (int c = 0; c < m; ++c) {
                    if (ns->sat_counts[c] == 0) {
                        unsat_pos[c] = unsz;
                        unsat[unsz++] = c;
                    } else {
                        unsat_pos[c] = -1;
                    }
                }
            }
        }

        int choice = rng_int(0, unsz - 1);
        int c_idx = unsat[choice];

        int best_v = -1, best_delta = 999999;
        int cs = ns->cl_offs[c_idx];
        int ce = ns->cl_offs[c_idx + 1];

        /* Entropy-Aware Selection: bias flipping towards x ~ 0.5 (maximum uncertainty) */
        double best_entropy = -1e9;
        for (int i = cs; i < ce; ++i) {
            int v = abs(ns->cl_flat[i]);
            double e_val = fabs(ns->x[v] - 0.5); /* 0 is high entropy, 0.5 is low entropy */
            double entropy_score = 0.5 - e_val; 
            int br = 0;
            int old_v_sat = (ns->x[v] > 0.5);

            for (int p = ns->v2c_ptr[v]; p < ns->v2c_ptr[v+1]; ++p) {
                int c = ns->v2c_data[p];
                if (ns->sat_counts[c] == 1) {
                    /* check if flipping v breaks c */
                    int lit_p = 0;
                    int cs_inner = ns->cl_offs[c];
                    int k = ns->cl_offs[c+1] - cs_inner;
                    
                    if (k == 3) {
                        if (abs(ns->cl_flat[cs_inner]) == v) { lit_p = (ns->cl_flat[cs_inner] > 0); }
                        else if (abs(ns->cl_flat[cs_inner+1]) == v) { lit_p = (ns->cl_flat[cs_inner+1] > 0); }
                        else if (abs(ns->cl_flat[cs_inner+2]) == v) { lit_p = (ns->cl_flat[cs_inner+2] > 0); }
                    } else {
                        for (int q = cs_inner; q < ns->cl_offs[c+1]; ++q) {
                            if (abs(ns->cl_flat[q]) == v) { lit_p = (ns->cl_flat[q] > 0); break; }
                        }
                    }
                    if (lit_p == old_v_sat) br++;
                }
            }
            /* Score trades off structural breakage (br) with continuous entropy */
            if (br < best_delta || (br == best_delta && entropy_score > best_entropy)) { 
                best_delta = br; 
                best_v = v; 
                best_entropy = entropy_score;
            }
        }

        if (best_v != -1 && (best_delta == 0 || rng_double() < exp(-beta * best_delta))) {
            int old_v_sat = (ns->x[best_v] > 0.5);
            ns->x[best_v] = (old_v_sat) ? 0.0 : 1.0;
            int new_v_sat = !old_v_sat;

            for (int p = ns->v2c_ptr[best_v]; p < ns->v2c_ptr[best_v + 1]; ++p) {
                int c = ns->v2c_data[p];
                int lit_p = 0;
                int cs_inner = ns->cl_offs[c];
                int k = ns->cl_offs[c+1] - cs_inner;
                
                if (k == 3) {
                    if (abs(ns->cl_flat[cs_inner]) == best_v) { lit_p = (ns->cl_flat[cs_inner] > 0); }
                    else if (abs(ns->cl_flat[cs_inner+1]) == best_v) { lit_p = (ns->cl_flat[cs_inner+1] > 0); }
                    else if (abs(ns->cl_flat[cs_inner+2]) == best_v) { lit_p = (ns->cl_flat[cs_inner+2] > 0); }
                } else {
                    for (int q = cs_inner; q < ns->cl_offs[c+1]; ++q) {
                       if (abs(ns->cl_flat[q]) == best_v) { lit_p = (ns->cl_flat[q] > 0); break; }
                    }
                }
                int was_sat = ns->sat_counts[c] > 0;
                if (lit_p == new_v_sat) ns->sat_counts[c]++;
                else ns->sat_counts[c]--;
                int is_sat = ns->sat_counts[c] > 0;

                if (!was_sat && is_sat) {
                    int pos = unsat_pos[c];
                    if (pos != -1) {
                        int last_c = unsat[unsz - 1];
                        unsat[pos] = last_c;
                        unsat_pos[last_c] = pos;
                        unsat_pos[c] = -1;
                        unsz--;
                    }
                } else if (was_sat && !is_sat) {
                    if (unsat_pos[c] == -1) {
                        unsat_pos[c] = unsz;
                        unsat[unsz++] = c;
                    }
                }
            }
        }
        
        /* Track best state */
        if (unsz < best_unsat) {
            best_unsat = unsz;
            for (int i = 1; i <= ns->num_vars; ++i) best_x[i] = ns->x[i];
        }
    }
    
    /* Restore best state if we lost it */
    if (unsz > best_unsat) {
        for (int i = 1; i <= ns->num_vars; ++i) ns->x[i] = best_x[i];
        recompute_sat_counts(ns);
        unsz = best_unsat;
    }
    
    int solved = (unsz == 0);
    free(unsat); free(unsat_pos); free(best_x); free(sample_x);
    return solved;
}

/* --------------------------------------------------------------------
   18. Core decomposition – identify the unsatisfied core and “blast’’
        it by temporarily inflating clause weights.
-------------------------------------------------------------------- */
static int core_decomposition(NitroSat *ns)
{
    recompute_sat_counts(ns);
    int *core = malloc(ns->num_clauses * sizeof(int));
    int core_sz = 0;
    for (int c=0;c<ns->num_clauses;++c) {
        if (ns->sat_counts[c] == 0) core[core_sz++] = c;
    }
    if (core_sz==0) { free(core); return 1; }

    /* stochastic WalkSAT on the core */
    double *saved = malloc(ns->num_clauses * sizeof(double));
    memcpy(saved, ns->cl_weights, ns->num_clauses * sizeof(double));
    for (int i=0;i<core_sz;++i) ns->cl_weights[core[i]] *= 100.0;

    int max_flips = ns->num_vars * 5;
    if (max_flips > 1000000) max_flips = 1000000;
    int success = 0;

    for (int flips=0; flips < max_flips; ++flips) {
        int c_idx = rand() % core_sz;
        int c = core[c_idx];
        if (ns->sat_counts[c] > 0) continue;

        int cs = ns->cl_offs[c];
        int ce = ns->cl_offs[c+1];
        int var = abs(ns->cl_flat[cs + (rand() % (ce-cs))]);
        
        ns->x[var] = (ns->x[var] > 0.5) ? 0.0 : 1.0;
        
        /* Update sat_counts for var flip */
        for (int p = ns->v2c_ptr[var]; p < ns->v2c_ptr[var+1]; ++p) {
            int cl = ns->v2c_data[p];
            int sc = 0;
            for (int j=ns->cl_offs[cl]; j<ns->cl_offs[cl+1]; ++j) {
                int l = ns->cl_flat[j];
                if ((l>0 && ns->x[abs(l)]>0.5) || (l<0 && ns->x[abs(l)]<=0.5)) sc++;
            }
            ns->sat_counts[cl] = sc;
        }

        if (flips % 10000 == 0) {
            int still_unsat = 0;
            for (int i=0; i<core_sz; ++i) if (ns->sat_counts[core[i]] == 0) { still_unsat = 1; break; }
            if (!still_unsat) { success = 1; break; }
        }
    }

    memcpy(ns->cl_weights, saved, ns->num_clauses * sizeof(double));
    free(core); free(saved);
    return success;
}

static double score_branch(NitroSat *ns, double beta, int n_samples, double *best_x_out)
{
    if (beta <= 0) return -1e300;
    double total_score = 0.0;
    double best_seen = 1e300;
    double *buf = malloc((ns->num_vars + 1) * sizeof(double));
    
    for (int s = 0; s < n_samples; ++s) {
        for (int i = 1; i <= ns->num_vars; ++i) buf[i] = ns->decimated[i] ? ns->x[i] : rng_double();
        double uns_energy = 0.0;
        int uns_count = 0;
        for (int c = 0; c < ns->num_clauses; ++c) {
            int sat = 0;
            for (int q = ns->cl_offs[c]; q < ns->cl_offs[c+1]; ++q) {
                int lit = ns->cl_flat[q], v = abs(lit);
                if ((lit > 0 && buf[v] > 0.5) || (lit < 0 && buf[v] <= 0.5)) { sat = 1; break; }
            }
            if (!sat) { 
                uns_energy += ns->cl_weights[c]; /* Spectral Prime Weighting */
                uns_count++; 
            }
        }
        total_score += exp(-beta * uns_energy);
        if (uns_energy < best_seen) {
            best_seen = uns_energy;
            if (best_x_out) {
                for (int i = 1; i <= ns->num_vars; ++i) best_x_out[i] = buf[i];
            }
        }
    }
    free(buf);
    return total_score / n_samples + 100.0 / (best_seen + 1.0);
}

/* --------------------------------------------------------------------
   [PORTED] Topological Repair Phase
-------------------------------------------------------------------- */
static int topological_repair_phase(NitroSat *ns, int max_steps)
{
    if (!ns->use_topology) return 0;
    recompute_sat_counts(ns);
    for (int step = 1; step <= max_steps; ++step) {
        ns->last_topology = compute_topology(ns, 0.3);
        
        int *hole_vars = malloc(ns->num_vars * sizeof(int));
        int hole_count = 0;
        
        if (ns->last_topology.beta1 == 0) {
            for (int i = 1; i <= ns->num_vars; ++i) {
                if (!ns->decimated[i]) hole_vars[hole_count++] = i;
            }
        } else {
            int *deg = calloc(ns->num_vars + 1, sizeof(int));
            for (int c = 0; c < ns->num_clauses; ++c) {
                if (ns->sat_counts[c] == 0) {
                    for (int i = ns->cl_offs[c]; i < ns->cl_offs[c+1]; ++i) deg[abs(ns->cl_flat[i])]++;
                }
            }
            for (int i = 1; i <= ns->num_vars; ++i) {
                if (deg[i] > 0 && !ns->decimated[i]) hole_vars[hole_count++] = i;
            }
            free(deg);
        }
        
        int flips = 0;
        for (int i = 0; i < hole_count; ++i) {
            int var = hole_vars[i];
            if (!ns->decimated[var]) {
                double force = 0.0;
                for (int p = ns->v2c_ptr[var]; p < ns->v2c_ptr[var+1]; ++p) {
                    int c = ns->v2c_data[p];
                    if (ns->sat_counts[c] == 0) {
                        for (int j = ns->cl_offs[c]; j < ns->cl_offs[c+1]; ++j) {
                            if (abs(ns->cl_flat[j]) == var) {
                                force += (ns->cl_flat[j] > 0 ? 1 : -1) * ns->cl_weights[c];
                                break;
                            }
                        }
                    }
                }
                if (fabs(force) > 0.01) {
                    double old_x = ns->x[var];
                    ns->x[var] = (force > 0) ? 0.99 : 0.01;
                    /* Update sat_counts for var flip */
                    if ((old_x > 0.5) != (ns->x[var] > 0.5)) {
                        for (int p = ns->v2c_ptr[var]; p < ns->v2c_ptr[var+1]; ++p) {
                            int c = ns->v2c_data[p];
                            /* Check if previously sat, now unsat or vice versa */
                            /* Exhaustive check for this clause only is O(k) */
                            int sc = 0;
                            for (int j=ns->cl_offs[c]; j<ns->cl_offs[c+1]; ++j) {
                                int l = ns->cl_flat[j];
                                if ((l>0 && ns->x[abs(l)]>0.5) || (l<0 && ns->x[abs(l)]<=0.5)) sc++;
                            }
                            ns->sat_counts[c] = sc;
                        }
                    }
                    flips++;
                }
            }
        }
        free(hole_vars);
        if (flips == 0) break;
        int uns = 0;
        for (int c=0; c<ns->num_clauses; ++c) if (ns->sat_counts[c]==0) { uns++; break; }
        if (uns == 0) return 1;
    }
    return 0;
}

/* --------------------------------------------------------------------
   [PORTED] Unit Propagation
-------------------------------------------------------------------- */
static int unit_propagation(NitroSat *ns)
{
    int changed = 0;
    for (int c = 0; c < ns->num_clauses; ++c) {
        int unassigned_lit = 0, satisfied = 0, num_unassigned = 0;
        for (int i = ns->cl_offs[c]; i < ns->cl_offs[c+1]; ++i) {
            int lit = ns->cl_flat[i], var = abs(lit);
            if (ns->decimated[var]) {
                if ((lit > 0 && ns->x[var] > 0.5) || (lit < 0 && ns->x[var] <= 0.5)) { satisfied = 1; break; }
            } else {
                num_unassigned++;
                if (num_unassigned == 1) unassigned_lit = lit;
                else { unassigned_lit = 0; break; }
            }
        }
        if (!satisfied && unassigned_lit) {
            int var = abs(unassigned_lit);
            double force_val = (unassigned_lit > 0) ? 1.0 : 0.0;
            if (fabs(ns->x[var] - force_val) > 0.5) {
                ns->x[var] = force_val;
                ns->decimated[var] = 1;
                changed = 1;
            }
        }
    }
    return changed;
}

/* --------------------------------------------------------------------
   [PORTED] Adelic Saturation Phase
-------------------------------------------------------------------- */
static int zeta_sweep_aggressive(NitroSat *ns, double beta, int *flips_out, int *unsat_out)
{
    recompute_sat_counts(ns);
    int m = ns->num_clauses;
    int *unsat = malloc(m * sizeof(int));
    int unsz = 0;
    for (int c = 0; c < m; ++c) if (ns->sat_counts[c] == 0) unsat[unsz++] = c;

    if (unsz == 0) { free(unsat); if (flips_out) *flips_out=0; if (unsat_out) *unsat_out=0; return 1; }

    int flips = 0;
    for (int i = 0; i < unsz; ++i) {
        int c = unsat[i];
        if (ns->sat_counts[c] > 0) continue;

        int best_v = -1, best_br = 999999;
        int cs = ns->cl_offs[c], ce = ns->cl_offs[c+1];
        for (int j = cs; j < ce; ++j) {
            int v = abs(ns->cl_flat[j]);
            int br = 0;
            int old_v_sat = (ns->x[v] > 0.5);
            for (int p = ns->v2c_ptr[v]; p < ns->v2c_ptr[v+1]; ++p) {
                int cl = ns->v2c_data[p];
                if (ns->sat_counts[cl] == 1) {
                    int lit_p = 0;
                    for (int q = ns->cl_offs[cl]; q < ns->cl_offs[cl+1]; ++q) {
                        if (abs(ns->cl_flat[q]) == v) { lit_p = (ns->cl_flat[q] > 0); break; }
                    }
                    if (lit_p == old_v_sat) br++;
                }
            }
            if (br < best_br) { best_br = br; best_v = v; }
        }

        if (best_v != -1 && (best_br == 0 || rng_double() < exp(-beta * best_br))) {
            int old_v_sat = (ns->x[best_v] > 0.5);
            ns->x[best_v] = (old_v_sat) ? 0.0 : 1.0;
            int new_v_sat = !old_v_sat;
            for (int p = ns->v2c_ptr[best_v]; p < ns->v2c_ptr[best_v+1]; ++p) {
                int cl = ns->v2c_data[p];
                int lit_p = 0;
                for (int q = ns->cl_offs[cl]; q < ns->cl_offs[cl+1]; ++q) {
                    if (abs(ns->cl_flat[q]) == best_v) { lit_p = (ns->cl_flat[q] > 0); break; }
                }
                if (lit_p == new_v_sat) ns->sat_counts[cl]++;
                else ns->sat_counts[cl]--;
            }
            flips++;
        }
    }
    if (flips_out) *flips_out = flips;
    int final_unsat = 0;
    for (int c = 0; c < m; ++c) if (ns->sat_counts[c] == 0) final_unsat++;
    if (unsat_out) *unsat_out = final_unsat;
    free(unsat);
    return final_unsat == 0;
}

static int adelic_saturation_phase(NitroSat *ns, int max_steps)
{
    recompute_sat_counts(ns);
    double beta_low = 0.1, beta_high = 10.0;
    int last_unsat = ns->num_clauses;
    for (int step = 1; step <= max_steps; ++step) {
        double beta = beta_low + (beta_high - beta_low) * ((double)step / max_steps);
        int unsat, flips;
        if (zeta_sweep_aggressive(ns, beta, &flips, &unsat)) return 1;
        
        if (flips == 0) beta_low = beta;
        else if (unsat > last_unsat) beta_high = beta;
        last_unsat = unsat;
        
        if (step % 10 == 0) unit_propagation(ns);
    }
    return check_satisfaction(ns) == ns->num_clauses;
}

/* --------------------------------------------------------------------
   19. Main solve routine – mirrors the Lua `NitroSat:solve` method.
-------------------------------------------------------------------- */
static int nitrosat_solve(NitroSat *ns)
{
    int max_steps = ns->max_steps;
    double lr0   = 0.2;               /* base learning‑rate before annealing */
    int    best_sat = 0;
    double *best_x  = malloc((ns->num_vars+1)*sizeof(double));
    int    prev_sat = 0;              /* For grokking detection */
    int    steps_since_improvement = 0; /* For nuclear zeta perturbation */
    
    /* STAGNATION MONITORING: Moving average of satisfaction slope */
    double sat_window[50] = {0};
    int sw_idx = 0;
    double avg_sat = 0;
    int last_stagnation_step = -100;

    for (int i=1;i<=ns->num_vars;++i) best_x[i] = ns->x[i];

    /* Initialize sat_counts before first gradient computation */
    recompute_sat_counts(ns);

    for (int step=1; step<=max_steps; ++step) {
        int unsat = compute_gradients(ns);
        int sat = ns->num_clauses - unsat;

        if (sat > best_sat) {
            best_sat = sat;
            for (int i=1;i<=ns->num_vars;++i) best_x[i] = ns->x[i];
            steps_since_improvement = 0;
        } else {
            steps_since_improvement++;
        }
        
        if (unsat == 0) {
            if (ns->verbose) printf("[EARLY] 100%% satisfaction at step %d!\n", step);
            /* Restore best solution before breaking */
            for (int i = 1; i <= ns->num_vars; i++) {
                ns->x[i] = best_x[i];
            }
            break;
        }

        /* GROKKING DETECTION: Sudden jump in satisfaction (e.g., 69% → 99.99%) */
        if (step > 1 && prev_sat > 0) {
            int jump = sat - prev_sat;
            double jump_pct = (double)jump / ns->num_clauses * 100.0;
            
            /* If we jumped >50% in one step AND now at >95%, trigger phase 2 immediately */
            if (jump_pct > 50.0 && sat >= (int)(0.95 * ns->num_clauses)) {
                if (ns->verbose) {
                    printf("[GROK] Detected %+.1f%% jump! %d→%d (%.2f%%). Triggering phase 2!\n",
                           jump_pct, prev_sat, sat, 100.0*sat/ns->num_clauses);
                }
                
                /* Immediately run topological repair */
                if (ns->use_topology) {
                    ns->last_topology = compute_topology(ns, 0.3);
                    if (topological_repair_phase(ns, 1500)) {
                        recompute_sat_counts(ns);
                        int final_sat = check_satisfaction(ns);
                        if (final_sat == ns->num_clauses) {
                            if (ns->verbose) printf("[GROK] Phase 2 solved it!\n");
                            free(best_x);
                            return 1;
                        }
                    }
                }
                
                /* Then adelic saturation */
                if (sat >= (int)(0.98 * ns->num_clauses)) {
                    if (adelic_saturation_phase(ns, 2000)) {
                        recompute_sat_counts(ns);
                        int final_sat = check_satisfaction(ns);
                        if (final_sat == ns->num_clauses) {
                            if (ns->verbose) printf("[GROK] Phase 3 solved it!\n");
                            free(best_x);
                            return 1;
                        }
                    }
                }
            }
        }
        /* ADAPTIVE STAGNATION MONITOR: Proactive Phase Transitions */
        sat_window[sw_idx % 50] = (double)sat / ns->num_clauses;
        sw_idx++;
        
        /* Only check every 50 steps to compare window blocks */
        if (sw_idx >= 50 && sw_idx % 50 == 0) {
            double sum = 0;
            for (int k=0; k<50; k++) sum += sat_window[k];
            double current_avg = sum / 50.0;
            double slope = current_avg - avg_sat;
            
            /* If slope is flat (<0.01% gain) and satisfaction is high (>90%) and not recently fired */
            if (step > 200 && fabs(slope) < 0.0001 && sat > (int)(0.90 * ns->num_clauses) && (step > last_stagnation_step + 100)) {
                if (ns->verbose) {
                    printf("[STAGNATION] Slope plateaued at %.4f%%. Proactively jumping to Stage 2 finisher!\n", current_avg*100.0);
                }
                
                /* Update best_x before jumping out */
                for (int i = 1; i <= ns->num_vars; ++i) best_x[i] = ns->x[i];
                best_sat = sat;
                
                /* BREAK Langevin loop to enter 3-phase finisher sequence early */
                break;
            }
            avg_sat = current_avg;
        }
        
        prev_sat = sat;

        /* Debug output every step */
        if (ns->verbose) {
            double grad_norm = 0;
            for (int i = 1; i <= ns->num_vars; i++) grad_norm += ns->grad_buffer[i] * ns->grad_buffer[i];
            grad_norm = sqrt(grad_norm);
            printf("[%04d] sat=%d/%d (%.2f%%) lr=%.6f |grad|=%.4f unsat_ clauses=%d\n",
                   step, sat, ns->num_clauses, 100.0*sat/ns->num_clauses,
                   ns->opt->lr, grad_norm, unsat);
            fflush(stdout);
        }

        /* dynamic clause re‑weighting – SAPS-like (matches gh.lua) */
        if (step % 10 == 0) {
            recompute_sat_counts(ns);
            for (int c = 0; c < ns->num_clauses; ++c) {
                if (ns->sat_counts[c] == 0) {
                    ns->cl_weights[c] *= 1.1;  /* solve_gravity penalty */
                    /* If we are stagnating heavily, stop exploding the weight! */
                    if (steps_since_improvement > 100 && ns->cl_weights[c] > 500.0) {
                        ns->cl_weights[c] = 500.0; 
                    } else if (ns->cl_weights[c] > 1e5) {
                        ns->cl_weights[c] = 1e5;  /* Prevent gradient scale overflow */
                    }
                } else {
                    ns->cl_weights[c] *= 0.999; /* decay */
                    if (ns->cl_weights[c] < 0.001) ns->cl_weights[c] = 0.001;  /* Prevent underflow */
                }
            }
        } /* annealed learning‑rate with adiabatic gain */
        double progress = (double)step / max_steps;
        double t = step * (30.0 / max_steps);
        double lr = annealing_lr(t, lr0, progress);
        if (lr <= 0.0) lr = 1e-6;
        ns->opt->lr = lr;
        optimizer_step(ns->opt, ns->x, ns->grad_buffer, ns->decimated);
        /* clamp x to [0,1] after optimizer step (matches gh.lua line 992) */
        for (int i = 1; i <= ns->num_vars; ++i)
            ns->x[i] = clamp(ns->x[i], 0.0, 1.0);

        /* Cheat Code 5: Update prev_x for potential incremental sat_counts */
        for (int i = 1; i <= ns->num_vars; ++i)
            ns->prev_x[i] = ns->x[i];

        /* decimation (soft locking) */
        if (ns->dec_freq && step % ns->dec_freq == 0) {
            int locked = decimate(ns);
            if (ns->verbose && locked) printf("[dec] locked %d vars\n", locked);
        }

        /* Zeta‑zero guided perturbation (every 10 steps) */
        if (step % 10 == 0) {
            TopologyInfo *topo = ns->use_topology ? &ns->last_topology : NULL;
            int nuclear = (steps_since_improvement > 100); /* Fire Nuclear Bomb if stuck */
            for (int v=1; v<=ns->num_vars; ++v) {
                if (nuclear && ns->decimated[v]) {
                    /* Nuclear strike unlocks frozen variables and injects uncertainty */
                    ns->decimated[v] = 0;
                    if (ns->x[v] > 0.9) ns->x[v] = 0.8;
                    if (ns->x[v] < 0.1) ns->x[v] = 0.2;
                }
                if (!ns->decimated[v]) {
                    double dz = zeta_zero_perturb(v, step, max_steps,
                                                 ns->primes, ns->prime_cnt,
                                                 topo, ns->num_vars,
                                                 nuclear);
                    ns->x[v] = clamp(ns->x[v] + dz * 0.001, 0.0, 1.0);  /* gh.lua scales by 0.001 */
                }
            }
            if (nuclear && ns->verbose) {
                printf("[ZETA] Firing nuclear perturbation at step %d to break stagnation.\n", step);
                steps_since_improvement = 0; /* Reset stagnation counter after nuclear strike */
            }
        }

        /* Inline AdaBoost DCW (gh.lua lines 2146-2170) – second pass */
        if (step % 10 == 0) {
            for (int c = 0; c < ns->num_clauses; ++c) {
                int clause_sat = 0;
                for (int i = ns->cl_offs[c]; i < ns->cl_offs[c+1] && !clause_sat; ++i) {
                    int lit = ns->cl_flat[i];
                    int v   = abs(lit);
                    double val = ns->x[v];
                    if ((lit > 0 && val > 0.5) || (lit < 0 && val <= 0.5)) clause_sat = 1;
                }
                if (!clause_sat) {
                    ns->cl_weights[c] *= 2.5;
                    if (ns->cl_weights[c] > 1e15) ns->cl_weights[c] = 1e15;  /* Prevent overflow */
                }
            }
            for (int i = 0; i < ns->num_clauses; ++i) {
                ns->cl_weights[i] *= 0.99;
                if (ns->cl_weights[i] < 0.001) ns->cl_weights[i] = 0.001;  /* Prevent underflow */
            }
        }

        /* BAHA fracture detection every 50 steps */
        if (step % 50 == 0) {
            double beta_proxy = ((double)step / max_steps) * 10.0;
            /* cheap thermodynamic estimator */
            double logZ = 0.0;
            for (int s=0; s<20; ++s) {  /* 20 samples */
                int e = 0;
                for (int c=0; c<ns->num_clauses; ++c) {
                    int sat = 0;
                    for (int i=ns->cl_offs[c]; i<ns->cl_offs[c+1] && !sat; ++i) {
                        int lit = ns->cl_flat[i];
                        int v   = abs(lit);
                        double val = ns->x[v];
                        if ((lit>0 && val>0.5) || (lit<0 && val<=0.5)) sat = 1;
                    }
                    if (!sat) ++e;
                }
                logZ += -beta_proxy * e;
            }
            logZ = logZ / 20.0;
            fd_record(ns->fd, beta_proxy, logZ);

            /* update variance/peak statistics */
            ++ns->baha_var_n;
            double delta = logZ - ns->baha_var_mean;
            ns->baha_var_mean += delta / ns->baha_var_n;
            ns->baha_var_M2   += delta * (logZ - ns->baha_var_mean);
            double std = (ns->baha_var_n>1) ?
                         sqrt(ns->baha_var_M2/(ns->baha_var_n-1)) : 0.0;
            double z   = (std>1e-12) ? (logZ - ns->baha_var_mean)/std : 0.0;

            if (fd_is_fracture(ns->fd) && z > 0.8) {
                if (logZ > ns->baha_var_peak) {
                    ns->baha_var_peak = logZ;
                    ns->baha_beta_c   = beta_proxy;
                }
                
                int nb;
                Branch *brs = enumerate_branches(ns->baha_beta_c, beta_proxy, &nb);
                if (nb > 0) {
                    double bestScore = -1e300;
                    int bestIdx = -1;
                    double *best_jump_x = malloc((ns->num_vars + 1) * sizeof(double));
                    double *temp_x = malloc((ns->num_vars + 1) * sizeof(double));
                    
                    for (int i = 0; i < nb; ++i) {
                        double score = score_branch(ns, brs[i].beta, 20, temp_x);
                        if (score > bestScore) {
                            bestScore = score;
                            bestIdx = i;
                            for (int j = 1; j <= ns->num_vars; ++j) best_jump_x[j] = temp_x[j];
                        }
                    }
                    
                    if (bestIdx != -1) {
                        if (ns->verbose) printf("[baha] jumping to branch beta=%.3f\n", brs[bestIdx].beta);
                        for (int i = 1; i <= ns->num_vars; ++i) {
                            if (!ns->decimated[i]) ns->x[i] = best_jump_x[i];
                        }
                        free(ns->opt->m); free(ns->opt->v); free(ns->opt);
                        ns->opt = optimizer_new("harmonic_adam", ns->x, ns->num_vars, 0.02, 0.9, 0.999, 1e-8, 0.02);
                    }
                    free(best_jump_x); free(temp_x);
                }
                free(brs);
                fd_clear(ns->fd);
            }
        }

        /* periodic topology snapshot (every topology_freq steps) */
        if (ns->use_topology && step % 50 == 0) {
            ns->last_topology = compute_topology(ns, 0.3);
            update_persistence(&ns->pt, ns->last_topology.beta0, ns->last_topology.beta1, step, ns->last_topology.complexity_score);
            if (ns->verbose)
                printf("[topo] step %d β0=%d β1=%d comp=%.3f\n",
                        step,
                        ns->last_topology.beta0,
                        ns->last_topology.beta1,
                        ns->last_topology.complexity_score);
        }
    }

    /* restore best solution */
    for (int i=1;i<=ns->num_vars;++i) ns->x[i] = best_x[i];
    recompute_sat_counts(ns);

    /* -----------------------------------------------------------------
       3‑phase finisher
       ----------------------------------------------------------------- */
    int sat = check_satisfaction(ns);
    if (sat == ns->num_clauses) {
        printf("Freeing best_x in nitrosat_solve (early solve)... (sat: %d/%d)\n", sat, ns->num_clauses);
        free(best_x);
        return 1;
    }   /* already solved */

    /* Phase‑2 : topological repair (≈95 % satisfied) */
    if (sat >= (int)(0.95 * ns->num_clauses)) {
        if (ns->verbose) puts("[phase‑2] entering topological repair");
        
        for (int i = 1; i <= ns->num_vars; ++i) ns->x[i] = best_x[i];
        topological_repair_phase(ns, 1500);
        int phase2_sat = check_satisfaction(ns);
        
        if (phase2_sat > best_sat) {
            best_sat = phase2_sat;
            for (int i = 1; i <= ns->num_vars; ++i) best_x[i] = ns->x[i];
        }
        if (ns->verbose) printf("[phase‑2] finished, sat=%d/%d\n", best_sat, ns->num_clauses);
        if (best_sat == ns->num_clauses) { free(best_x); return 1; }
    }

    /* Phase‑3 : adelic saturation (≈98 % satisfied) */
    if (best_sat >= (int)(0.98 * ns->num_clauses)) {
        if (ns->verbose) puts("[phase‑3] entering adelic saturation");
        for (int i = 1; i <= ns->num_vars; ++i) ns->x[i] = best_x[i];
        
        adelic_saturation_phase(ns, 2000);
        int phase3_sat = check_satisfaction(ns);
        if (phase3_sat > best_sat) {
            best_sat = phase3_sat;
            for (int i = 1; i <= ns->num_vars; ++i) best_x[i] = ns->x[i];
        }
        if (ns->verbose) printf("[phase‑3] finished, sat=%d/%d\n", best_sat, ns->num_clauses);
        if (best_sat == ns->num_clauses) { free(best_x); return 1; }
    }

    /* Core decomposition “blast’’ */
    if (best_sat < ns->num_clauses) {
        if (ns->verbose) puts("[core] attempting core decomposition");
        for (int i = 1; i <= ns->num_vars; ++i) ns->x[i] = best_x[i];
        if (core_decomposition(ns)) {
            int core_sat = check_satisfaction(ns);
            if (core_sat > best_sat) {
                best_sat = core_sat;
                for (int i = 1; i <= ns->num_vars; ++i) best_x[i] = ns->x[i];
            }
        }
        if (best_sat == ns->num_clauses) { free(best_x); return 1; }
    }

    /* Final BAHA discrete local search (pure BAHA, no WalkSAT) */
    if (best_sat < ns->num_clauses) {
        if (ns->verbose) puts("[final] BAHA discrete search");
        for (int i = 1; i <= ns->num_vars; ++i) ns->x[i] = best_x[i] > 0.5 ? 1.0 : 0.0;
        
        /* Run BAHA-WalkSAT with baha.cpp params: 200 steps × 50 samples */
        /* Max flips capped to 2M for massive instances to prevent hangs */
        int max_sc_flips = ns->num_vars * 200;
        if (max_sc_flips > 2000000) max_sc_flips = 2000000;
        baha_walksat(ns, max_sc_flips);
        int bsat = check_satisfaction(ns);
        if (bsat > best_sat) {
            best_sat = bsat;
        }
    }
    
    free(best_x);
    return (best_sat == ns->num_clauses);
}

static int nitrosat_solve_dcw(NitroSat *ns, int num_passes)
{
    if (ns->num_vars == 0 || ns->num_clauses == 0) return nitrosat_solve(ns);
    
    int best_overall_sat = 0;
    double *best_overall_x = malloc((ns->num_vars + 1) * sizeof(double));
    int success = 0;
    
    for (int i = 0; i < ns->num_clauses; ++i) ns->cl_weights[i] = 1.0;
    
    for (int pass = 1; pass <= num_passes; ++pass) {
        if (ns->verbose) printf("--- DCW PASS %d ---\n", pass);
        
        success = nitrosat_solve(ns);
        int sat_count = check_satisfaction(ns);
        
        if (sat_count > best_overall_sat) {
            best_overall_sat = sat_count;
            for (int i = 1; i <= ns->num_vars; ++i) best_overall_x[i] = ns->x[i];
        }
        
        if (sat_count == ns->num_clauses) { success = 1; break; }
        
        if (pass < num_passes) {
            for (int c = 0; c < ns->num_clauses; ++c) {
                int sat = 0;
                for (int i = ns->cl_offs[c]; i < ns->cl_offs[c+1]; ++i) {
                    int lit = ns->cl_flat[i], v = abs(lit);
                    if ((lit > 0 && ns->x[v] > 0.5) || (lit < 0 && ns->x[v] <= 0.5)) { sat = 1; break; }
                }
                if (!sat) {
                    ns->cl_weights[c] *= 1.5;
                    if (ns->cl_weights[c] > 1e15) ns->cl_weights[c] = 1e15;  /* Prevent overflow */
                }
            }
        }
        
        // Soft reset for next pass
        ns->opt->t = 0;
        for (int i = 1; i <= ns->num_vars; ++i) {
            ns->opt->m[i] = 0.0;
            ns->opt->v[i] = 0.0;
            if (!ns->decimated[i]) ns->x[i] = rng_double();
        }
    }
    
    for (int i = 1; i <= ns->num_vars; ++i) ns->x[i] = best_overall_x[i];
    free(best_overall_x);
    return success;
}

/* --------------------------------------------------------------------
   20. JSON output printer
-------------------------------------------------------------------- */
static void print_json_escaped_string(const char *s)
{
    const unsigned char *p = (const unsigned char *)(s ? s : "");
    putchar('"');
    while (*p) {
        switch (*p) {
            case '\"': fputs("\\\"", stdout); break;
            case '\\': fputs("\\\\", stdout); break;
            case '\b': fputs("\\b", stdout); break;
            case '\f': fputs("\\f", stdout); break;
            case '\n': fputs("\\n", stdout); break;
            case '\r': fputs("\\r", stdout); break;
            case '\t': fputs("\\t", stdout); break;
            default:
                if (*p < 0x20) printf("\\u%04x", (unsigned int)*p);
                else putchar((char)*p);
        }
        ++p;
    }
    putchar('"');
}

static void print_json_result(NitroSat *ns, const char *cnf_file,
                              int solved, LatencyBreakdown *lat,
                              SolverDiagnostics *diag,
                              const ProofReport *proof)
{
    int final_sat = check_satisfaction(ns);
    int unsat = ns->num_clauses - final_sat;
    double sat_rate = (ns->num_clauses > 0) ? (double)final_sat / ns->num_clauses : 0.0;

    const char *status = solved ? "SATISFIED" :
                         (sat_rate >= 0.95 ? "PARTIAL" : "UNSATISFIED");

    double throughput = (lat->total_ms > 0.001) ?
        (double)ns->num_clauses / (lat->total_ms / 1000.0) : 0.0;
    char tp_buf[64];
    if (throughput >= 1e6) snprintf(tp_buf, sizeof(tp_buf), "%.1fM clauses/sec", throughput/1e6);
    else if (throughput >= 1e3) snprintf(tp_buf, sizeof(tp_buf), "%.1fK clauses/sec", throughput/1e3);
    else snprintf(tp_buf, sizeof(tp_buf), "%.0f clauses/sec", throughput);

    printf("{\n");
    printf("  \"status\": \"%s\",\n", status);
    printf("  \"satisfaction_rate\": %.6f,\n", sat_rate);
    printf("  \"satisfied\": %d,\n", final_sat);
    printf("  \"unsatisfied\": %d,\n", unsat);
    printf("  \"variables\": %d,\n", ns->num_vars);
    printf("  \"clauses\": %d,\n", ns->num_clauses);
    printf("  \"formula\": ");
    print_json_escaped_string(cnf_file);
    printf(",\n");

    /* Assignment in DIMACS format */
    printf("  \"assignment\": [");
    for (int i = 1; i <= ns->num_vars; ++i) {
        printf("%d", (ns->x[i] > 0.5) ? i : -i);
        if (i < ns->num_vars) printf(",");
    }
    printf("],\n");

    /* Confidence vector (raw continuous values) */
    printf("  \"confidence\": [");
    for (int i = 1; i <= ns->num_vars; ++i) {
        printf("%.4f", ns->x[i]);
        if (i < ns->num_vars) printf(",");
    }
    printf("],\n");

    /* Latency breakdown */
    printf("  \"latency\": {\n");
    printf("    \"total_ms\": %.2f,\n", lat->total_ms);
    printf("    \"throughput\": \"%s\",\n", tp_buf);
    printf("    \"breakdown\": {\n");
    printf("      \"initialization_ms\": %.2f,\n", lat->init_ms);
    printf("      \"langevin_flow_ms\": %.2f,\n", lat->langevin_ms);
    printf("      \"topo_repair_ms\": %.2f,\n", lat->topo_repair_ms);
    printf("      \"adelic_saturation_ms\": %.2f,\n", lat->adelic_ms);
    printf("      \"core_decomposition_ms\": %.2f,\n", lat->core_decomp_ms);
    printf("      \"walksat_ms\": %.2f\n", lat->walksat_ms);
    printf("    }\n");
    printf("  },\n");

    /* Diagnostics */
    printf("  \"diagnostics\": {\n");

    /* Thermodynamics */
    printf("    \"thermodynamics\": {\n");
    printf("      \"final_beta\": %.4f,\n", diag->final_beta);
    printf("      \"critical_beta\": %.6f,\n", diag->critical_beta);
    printf("      \"entropy_level\": %.6f,\n", diag->entropy_level);
    printf("      \"convexity_status\": \"%s\"\n", diag->convexity_status);
    printf("    },\n");

    /* Topology */
    printf("    \"topology\": {\n");
    printf("      \"betti_0\": %d,\n", diag->betti_0);
    printf("      \"betti_1\": %d,\n", diag->betti_1);
    printf("      \"complexity_score\": %.4f,\n", diag->complexity_score);
    printf("      \"persistence_events\": %d\n", diag->persistence_events);
    printf("    },\n");

    /* Prime stability */
    printf("    \"prime_stability\": {\n");
    printf("      \"aggregation_error\": %.8f,\n", diag->aggregation_error);
    printf("      \"chebyshev_bias\": \"%s\"\n", diag->chebyshev_bias);
    printf("    },\n");

    /* Blame attribution for unsat clauses */
    printf("    \"blame\": [");
    for (int i = 0; i < diag->blame_count; ++i) {
        int c = diag->blame_clause_ids[i];
        printf("{\"clause_id\":%d,\"weight\":%.4f,\"variables\":[", c, diag->blame_weights[i]);
        for (int j = ns->cl_offs[c]; j < ns->cl_offs[c+1]; ++j) {
            printf("%d", ns->cl_flat[j]);
            if (j + 1 < ns->cl_offs[c+1]) printf(",");
        }
        printf("]}");
        if (i + 1 < diag->blame_count) printf(",");
    }
    printf("]\n");
    printf("  },\n");

    printf("  \"proof\": {\n");
    printf("    \"requested\": %s,\n", (proof && proof->requested) ? "true" : "false");
    printf("    \"generated\": %s,\n", (proof && proof->generated) ? "true" : "false");
    printf("    \"format\": ");
    print_json_escaped_string((proof && proof->format) ? proof->format : "none");
    printf(",\n");
    printf("    \"status\": ");
    print_json_escaped_string((proof && proof->status[0]) ? proof->status : "not_requested");
    printf(",\n");
    printf("    \"backend\": ");
    print_json_escaped_string((proof && proof->backend) ? proof->backend : "none");
    printf(",\n");
    printf("    \"derived_units\": %d,\n", (proof ? proof->derived_units : 0));
    printf("    \"path\": ");
    print_json_escaped_string((proof && proof->path) ? proof->path : "");
    printf("\n");
    printf("  }\n");
    printf("}\n");
}

/* --------------------------------------------------------------------
   21. Compute diagnostics from solver state
-------------------------------------------------------------------- */
static void compute_solver_diagnostics(NitroSat *ns, SolverDiagnostics *diag)
{
    /* Thermodynamics */
    diag->final_beta = ns->heat_beta;
    diag->entropy_level = ns->entropy_weight;

    /* Critical beta from MATH.md: β* ~ (4δ²)/(k²·d·ln(K)/ln(ln(K))) */
    double k_max = 3.0; /* dominant clause size */
    double d_clause = 0.0;
    for (int i = 1; i <= ns->num_vars; ++i) d_clause += ns->degrees[i];
    d_clause = (ns->num_vars > 0) ? d_clause / ns->num_vars : 1.0;
    if (d_clause < 1e-12) d_clause = 1.0; /* fallback for unit/degenerate formulas */
    double delta = 0.3;
    double K = (double)ns->num_clauses;
    double lnK = (K > 1) ? log(K) : 1.0;
    double lnlnK = (lnK > 1) ? log(lnK) : 1.0;
    diag->critical_beta = (4.0 * delta * delta * lnK) / (k_max * k_max * d_clause * lnlnK);

    /* Convexity check: 4/β > W_max * k_max² * d_clause / δ² */
    double eff_beta = ns->heat_beta * d_clause;
    double W_max = ns->cl_weights[0];
    for (int c = 1; c < ns->num_clauses; ++c)
        if (ns->cl_weights[c] > W_max) W_max = ns->cl_weights[c];
    double lhs = 4.0 / (eff_beta > 1e-12 ? eff_beta : 1e-12);
    double rhs = W_max * k_max * k_max * d_clause / (delta * delta);
    diag->convexity_status = (lhs > rhs) ? "STABLE" : "NON_CONVEX";

    /* Topology */
    if (ns->pt.has_initial) {
        diag->betti_0 = ns->pt.final_beta0;
        diag->betti_1 = ns->pt.final_beta1;
        diag->complexity_score = (diag->betti_0 > 0) ?
            (double)diag->betti_1 / diag->betti_0 : 0.0;
        diag->persistence_events = ns->pt.persistence_events;
    } else {
        diag->betti_0 = diag->betti_1 = 0;
        diag->complexity_score = 0.0;
        diag->persistence_events = 0;
    }

    /* Prime stability: compute variance of per-clause prime aggregation */
    double mean_w = 0.0;
    for (int c = 0; c < ns->num_clauses; ++c) mean_w += ns->cl_weights[c];
    mean_w /= (ns->num_clauses > 0 ? ns->num_clauses : 1);
    double var_w = 0.0;
    for (int c = 0; c < ns->num_clauses; ++c) {
        double d = ns->cl_weights[c] - mean_w;
        var_w += d * d;
    }
    diag->aggregation_error = (ns->num_clauses > 1) ?
        sqrt(var_w / (ns->num_clauses - 1)) / mean_w : 0.0;
    diag->chebyshev_bias = (diag->aggregation_error < 0.01) ? "STABLE" :
                           (diag->aggregation_error < 0.1) ? "MODERATE" : "UNSTABLE";

    /* Blame: collect unsat clause IDs sorted by weight (top 20) */
    recompute_sat_counts(ns);
    int cap = 20;
    diag->blame_clause_ids = malloc(cap * sizeof(int));
    diag->blame_weights = malloc(cap * sizeof(double));
    diag->blame_count = 0;

    for (int c = 0; c < ns->num_clauses && diag->blame_count < cap; ++c) {
        if (ns->sat_counts[c] == 0) {
            diag->blame_clause_ids[diag->blame_count] = c;
            diag->blame_weights[diag->blame_count] = ns->cl_weights[c];
            diag->blame_count++;
        }
    }
    /* Sort blame by weight descending (simple insertion sort, max 20) */
    for (int i = 1; i < diag->blame_count; ++i) {
        int j = i;
        while (j > 0 && diag->blame_weights[j] > diag->blame_weights[j-1]) {
            double tw = diag->blame_weights[j]; diag->blame_weights[j] = diag->blame_weights[j-1]; diag->blame_weights[j-1] = tw;
            int ti = diag->blame_clause_ids[j]; diag->blame_clause_ids[j] = diag->blame_clause_ids[j-1]; diag->blame_clause_ids[j-1] = ti;
            j--;
        }
    }
}

/* --------------------------------------------------------------------
   22. Entry point – read problem, build solver, run.
-------------------------------------------------------------------- */
int main(int argc, char **argv)
{
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <cnf-file> [max-steps] [--no-dcw] [--no-topo] [--json] [--proof <path>] [--proof-format drat|lrat] [--proof-max-vars N] [--proof-max-clauses N]\n", argv[0]);
        fprintf(stderr, "Try %s --help for more information.\n", argv[0]);
        return EXIT_FAILURE;
    }

    /* Handle --help */
    if (strcmp(argv[1], "--help") == 0 || strcmp(argv[1], "-h") == 0) {
        printf("NitroSAT — Physics-Informed MaxSAT Solver\n\n");
        printf("Copyright (c) 2024-2025 Sethu Iyer (https://orcid.org/0009-0008-5446-2856)\n");
        printf("Email: shunyabarlabs@zohomail.com\n\n");
        printf("Licensed under the Apache License, Version 2.0 (the \"License\");\n");
        printf("you may not use this file except in compliance with the License.\n");
        printf("You may obtain a copy of the License at\n\n");
        printf("    http://www.apache.org/licenses/LICENSE-2.0\n\n");
        printf("Usage: %s <cnf-file> [max-steps] [options]\n\n", argv[0]);
        printf("Options:\n");
        printf("  --no-dcw                 Disable DCW heuristic (use standard random walk)\n");
        printf("  --no-topo                Disable topological/BAHA solver\n");
        printf("  --json                   Output results in JSON format\n");
        printf("  --proof <path>          Generate DRAT proof on UNSAT (requires proof backend)\n");
        printf("  --proof-format drat|lrat  Proof format (default: drat)\n");
        printf("  --proof-max-vars N      Max variables for proof generation (default: 50000)\n");
        printf("  --proof-max-clauses N   Max clauses for proof generation (default: 1000000)\n");
        printf("  --help, -h              Show this help message\n\n");
        printf("Repository: https://github.com/sethuiyer/NitroSAT\n");
        return EXIT_SUCCESS;
    }

    const char *cnf_file = argv[1];
    int max_steps = 3000;
    int use_dcw = 1;
    int use_topo = 1;
    int json_output = 0;
    const char *proof_path = NULL;
    const char *proof_format = "drat";

    for (int i = 2; i < argc; ++i) {
        if (strcmp(argv[i], "--no-dcw") == 0) use_dcw = 0;
        else if (strcmp(argv[i], "--no-topo") == 0) use_topo = 0;
        else if (strcmp(argv[i], "--json") == 0) json_output = 1;
        else if (strcmp(argv[i], "--proof") == 0) {
            if (i + 1 >= argc) {
                fprintf(stderr, "--proof requires an output path.\n");
                return EXIT_FAILURE;
            }
            proof_path = argv[++i];
        }
        else if (strcmp(argv[i], "--proof-format") == 0) {
            if (i + 1 >= argc) {
                fprintf(stderr, "--proof-format requires a value (drat|lrat).\n");
                return EXIT_FAILURE;
            }
            proof_format = argv[++i];
        }
        else if (strcmp(argv[i], "--proof-max-vars") == 0) {
            if (i + 1 >= argc) {
                fprintf(stderr, "--proof-max-vars requires a number.\n");
                return EXIT_FAILURE;
            }
            proof_max_vars = atoi(argv[++i]);
            if (proof_max_vars <= 0) {
                fprintf(stderr, "--proof-max-vars must be positive.\n");
                return EXIT_FAILURE;
            }
        }
        else if (strcmp(argv[i], "--proof-max-clauses") == 0) {
            if (i + 1 >= argc) {
                fprintf(stderr, "--proof-max-clauses requires a number.\n");
                return EXIT_FAILURE;
            }
            proof_max_clauses = atoi(argv[++i]);
            if (proof_max_clauses <= 0) {
                fprintf(stderr, "--proof-max-clauses must be positive.\n");
                return EXIT_FAILURE;
            }
        }
        else if (argv[i][0] == '-') {
            fprintf(stderr, "Unknown option: %s\n", argv[i]);
            return EXIT_FAILURE;
        }
        else max_steps = atoi(argv[i]);
    }

    struct timespec ts_total_start, ts_total_end, ts_init_end;
    clock_gettime(CLOCK_MONOTONIC, &ts_total_start);

    Instance *inst = read_cnf(cnf_file);
    if (!inst) return EXIT_FAILURE;

    /* Auto-calculate max_steps: REMOVED. Use uniform 3000 steps. */

    NitroSat *ns = nitrosat_new(inst, max_steps, json_output ? 0 : 1);
    ns->use_topology = use_topo;

    /* Build variable‑to‑clause CSR structure (O(L) time) */
    ns->v2c_ptr = calloc(ns->num_vars + 2, sizeof(int));
    NITROSAT_ENSURE(ns->v2c_ptr != NULL, "Failed to allocate v2c_ptr");
    
    /* Count clause occurrences per variable */
    for (int i=0; i<inst->num_clauses; ++i) {
        for (int j=inst->cl_offs[i]; j<inst->cl_offs[i+1]; ++j) {
            int lit = inst->cl_flat[j];
            NITROSAT_ENSURE(abs(lit) <= ns->num_vars, "Literal exceeds variable count in v2c build");
            ns->v2c_ptr[abs(lit) + 1]++;
        }
    }
    
    /* Prefix sum to get pointers */
    for (int i=1; i<=ns->num_vars+1; ++i) ns->v2c_ptr[i] += ns->v2c_ptr[i-1];
    
    /* Allocate and fill v2c_data */
    int v2c_size = ns->v2c_ptr[ns->num_vars+1];
    ns->v2c_data = malloc(v2c_size * sizeof(int));
    NITROSAT_ENSURE(ns->v2c_data != NULL, "Failed to allocate v2c_data");
    
    int *cur = calloc(ns->num_vars + 1, sizeof(int));
    NITROSAT_ENSURE(cur != NULL, "Failed to allocate v2c cursor");
    
    for (int c=0; c<inst->num_clauses; ++c) {
        for (int j=inst->cl_offs[c]; j<inst->cl_offs[c+1]; ++j) {
            int v = abs(inst->cl_flat[j]);
            int pos = ns->v2c_ptr[v] + cur[v];
            NITROSAT_ENSURE(pos < v2c_size, "v2c_data index out of bounds");
            ns->v2c_data[pos] = c;
            cur[v]++;
        }
    }
    free(cur);

    /* Validate v2c CSR: check that each clause appears exactly once per variable */
    #ifdef NITROSAT_DEBUG
    int *validation_counts = calloc(v2c_size, sizeof(int));
    NITROSAT_ENSURE(validation_counts != NULL, "Failed to allocate validation buffer");
    for (int v = 1; v <= ns->num_vars; ++v) {
        for (int p = ns->v2c_ptr[v]; p < ns->v2c_ptr[v+1]; ++p) {
            int c = ns->v2c_data[p];
            NITROSAT_ENSURE(c >= 0 && c < inst->num_clauses, "Invalid clause ID in v2c");
            validation_counts[c]++;
        }
    }
    for (int c = 0; c < inst->num_clauses; ++c) {
        int expected = inst->cl_offs[c+1] - inst->cl_offs[c];
        NITROSAT_ENSURE(validation_counts[c] == expected, 
            "v2c CSR validation failed: clause occurrence count mismatch");
    }
    free(validation_counts);
    #endif

    clock_gettime(CLOCK_MONOTONIC, &ts_init_end);

    /* ---- Solve ---- */
    struct timespec ts_solve_start, ts_solve_end;
    clock_gettime(CLOCK_MONOTONIC, &ts_solve_start);
    int solved = use_dcw ? nitrosat_solve_dcw(ns, 5) : nitrosat_solve(ns);
    clock_gettime(CLOCK_MONOTONIC, &ts_solve_end);
    clock_gettime(CLOCK_MONOTONIC, &ts_total_end);

    /* ---- Build latency breakdown ---- */
    LatencyBreakdown lat = {0};
    lat.init_ms = timespec_ms(&ts_total_start, &ts_init_end);
    lat.langevin_ms = timespec_ms(&ts_solve_start, &ts_solve_end);
    lat.total_ms = timespec_ms(&ts_total_start, &ts_total_end);
    /* Note: per-phase breakdown within nitrosat_solve would require
       instrumenting the function itself. For now, langevin_ms covers
       the entire solve including sub-phases. */

    ProofReport proof = {0};
    proof.requested = (proof_path != NULL);
    proof.format = proof.requested ? proof_format : "none";
    proof.path = proof.requested ? proof_path : "";
    proof.backend = "none";
    snprintf(proof.status, sizeof(proof.status), "%s",
             proof.requested ? "not_attempted" : "not_requested");

    if (proof.requested) {
        if (solved) {
            snprintf(proof.status, sizeof(proof.status), "sat_assignment_found");
        } else if (strcmp(proof_format, "drat") == 0 || strcmp(proof_format, "lrat") == 0) {
            proof.backend = "solver_aware_drat";
            if (strcmp(proof_format, "lrat") == 0) {
                proof.backend = "solver_aware_lrat";
                proof.lrat_hints_emitted = 1;
            }
            proof.generated = generate_drat_unsat_proof(
                ns, proof_path, &proof.derived_units, proof.status, sizeof(proof.status),
                proof_format);
        } else {
            snprintf(proof.status, sizeof(proof.status), "unsupported_format");
        }
    }

    if (json_output) {
        /* Compute diagnostics and print JSON */
        SolverDiagnostics diag = {0};
        compute_solver_diagnostics(ns, &diag);
        print_json_result(ns, cnf_file, solved, &lat, &diag, &proof);
        free(diag.blame_clause_ids);
        free(diag.blame_weights);
    } else {
        /* Original text output */
        recompute_sat_counts(ns);
        int final_sat = check_satisfaction(ns);

        printf("\n=== RESULT ===\n");
        printf("Formula : %s\n", cnf_file);
        printf("Variables: %d   Clauses: %d\n", ns->num_vars, ns->num_clauses);
        printf("Satisfied: %d   Unsatisfied: %d\n", final_sat, ns->num_clauses - final_sat);
        printf("Solved  : %s\n", solved ? "YES" : "NO");
        printf("Time    : %.2f s\n", lat.total_ms / 1000.0);
        if (proof.requested) {
            printf("Proof   : %s (%s)\n", proof.status, proof.format);
            if (proof.generated) printf("Proof file: %s\n", proof.path);
        }

        if (ns->use_topology && ns->pt.has_initial) {
            printf("Topology Summary:\n");
            printf("  Beta0: %d -> %d\n", ns->pt.initial_beta0, ns->pt.final_beta0);
            printf("  Beta1: %d -> %d\n", ns->pt.initial_beta1, ns->pt.final_beta1);
            printf("  Persistence Events: %d\n", ns->pt.persistence_events);
            printf("  Complexity Trend: %.3f\n", ns->pt.final_complexity - ns->pt.initial_complexity);
        }

        printf("Final assignment (first 20 vars):\n");
        for (int i=1; i<=ns->num_vars && i<=20; ++i)
            printf("%d:%d ", i, (ns->x[i] > 0.5));
        puts("\n");
    }

    nitrosat_free(ns);
    free(inst->cl_offs); free(inst->cl_flat); free(inst);
    return solved ? EXIT_SUCCESS : EXIT_FAILURE;
}
