
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
                 [--proof <path>] [--proof-format drat]

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
#define GAMMA           0.5772156649015328606   /* Euler-Mascheroni constant */
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
static long proof_max_clauses = 3000000000L; /* CLI: --proof-max-clauses (default 3B) */

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
    if (b <= a) return a;
    return a + (int)(rng_double() * (b - a + 1));
}

/* --------------------------------------------------------------------
   Utility
*/
#define PI  3.14159265358979323846
/* PHI already defined as ((1.0 + sqrt(5.0)) / 2.0) */
#define EPS 1e-12

/* Bit-Parallel Adelic Evaluation Macros */
#define SET_BIT(a, b) (a[(b)>>3] |= (1 << ((b)&7)))
#define GET_BIT(a, b) (a[(b)>>3] & (1 << ((b)&7)))
#define CLR_BIT(a, b) (a[(b)>>3] &= ~(1 << ((b)&7)))

static inline double clamp(double x, double lo, double hi) {
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
    if (out_cnt) *out_cnt = 0;
    if (n <= 0) return NULL;

    /* Estimate an upper bound for the nth prime (Rosser–Schoenfeld) */
    int limit = (n < 6) ? 15 : (int)(n * (log((double)n) + log(log((double)n))) + 3);
    char *is_prime = calloc(limit + 1, 1);
    if (!is_prime) return NULL;
    for (int i = 2; i <= limit; ++i) is_prime[i] = 1;
    for (int p = 2; p * p <= limit; ++p)
        if (is_prime[p])
            for (int q = p * p; q <= limit; q += p) is_prime[q] = 0;

    int *primes = malloc(sizeof(int) * n);
    if (!primes) {
        free(is_prime);
        return NULL;
    }
    int cnt = 0;
    for (int i = 2; i <= limit && cnt < n; ++i)
        if (is_prime[i]) primes[cnt++] = i;

    free(is_prime);
    if (out_cnt) *out_cnt = cnt;
    return primes;
}

/* --------------------------------------------------------------------
   3.  Optimiser --------------------------------------------------------
-------------------------------------------------------------------- */
typedef struct {
    double *m_standard; /* Exact signed momentum for 1D reals */
    double *m_phase;    /* Phase (angular) momentum */
    double *v_amp;      /* Amplitude variance */
    double *v_phase;    /* Phase velocity variance */
    double *prev_phase; /* Previous phase state */
    double  lr;
    double  beta1, beta2, eps, kappa;
    int     t;
    int     nv;
    double  resonance_amplitude;
    char    name[MAX_NAME_LEN];
} Optimizer;

static Optimizer *optimizer_new(const char *name, double *params, int nv,
                               double lr, double beta1, double beta2,
                               double eps, double resonance_amp)
{
    (void)params;
    Optimizer *opt = calloc(1, sizeof(Optimizer));
    NITROSAT_ENSURE(opt != NULL, "Failed to allocate optimizer");
    strncpy(opt->name, name, MAX_NAME_LEN-1);
    opt->nv    = nv;
    opt->lr    = lr;
    opt->beta1 = beta1;
    opt->beta2 = beta2;
    opt->eps   = eps;
    opt->kappa = 0.1;
    
    opt->m_standard = calloc(nv + 1, sizeof(double));
    opt->m_phase = calloc(nv + 1, sizeof(double));
    opt->v_amp = calloc(nv + 1, sizeof(double));
    opt->v_phase = calloc(nv + 1, sizeof(double));
    opt->prev_phase = calloc(nv + 1, sizeof(double));
    NITROSAT_ENSURE(opt->m_standard != NULL && opt->m_phase != NULL &&
                    opt->v_amp != NULL && opt->v_phase != NULL &&
                    opt->prev_phase != NULL,
                    "Failed to allocate optimizer buffers");
    
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

    double kappa = opt->kappa;

    for (int i = 1; i <= n; ++i) {
        if (decimated && decimated[i]) continue;
        
        double grad = grads[i];
        
        /* 1D WAdam: Decompose into Amplitude and Phase-Velocity momentum
           Amplitude momentum tracks the integrated force. 
           Phase momentum tracks the average directional bias (Wasserstein drift). */
        opt->m_standard[i] = b1 * opt->m_standard[i] + (1.0 - b1) * grad;
        opt->v_amp[i] = b2 * opt->v_amp[i] + (1.0 - b2) * grad * grad;
        
        double phase_cur = (grad < 0) ? PI : 0.0;
        double phase_diff = phase_cur - opt->prev_phase[i];
        while (phase_diff > PI) phase_diff -= 2.0 * PI;
        while (phase_diff <= -PI) phase_diff += 2.0 * PI;
        
        opt->m_phase[i] = b1 * opt->m_phase[i] + (1.0 - b1) * phase_cur;
        opt->v_phase[i] = b2 * opt->v_phase[i] + (1.0 - b2) * phase_diff * phase_diff;
        
        /* Bias correction */
        double m_hat = opt->m_standard[i] * t_inv1;
        double v_amp_hat = opt->v_amp[i] * t_inv2;
        double m_phase_hat = opt->m_phase[i] * t_inv1;
        double v_phase_hat = opt->v_phase[i] * t_inv2;
        
        /* Directional Damping (cos(m_phase)) - The 1D Projective WFR Step
           If a variable is oscillating between 0 and PI, m_phase_hat -> PI/2, and cos -> 0.
           The force cancels exactly at high resonance, preventing limit cycles. */
        double m_wfr = fabs(m_hat) * cos(m_phase_hat);
        if (m_hat < 0) m_wfr = -fabs(m_wfr); /* Preserve original direction or flip? User spec exp(1j*phase) */
        
        /* Preconditioning (Variance Interpolation) */
        double v_wfr = kappa * v_amp_hat + (1.0 - kappa) * v_phase_hat * (v_amp_hat / (PI * PI + 1e-12)); 
        
        double h_inv = 1.0 / (sqrt(v_wfr) + eps);
        params[i] -= lr * m_wfr * h_inv;
        
        opt->prev_phase[i] = phase_cur;
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
    fd->n = 0;
    fd->mean = fd->M2 = fd->last_rate = 0.0;
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
    NITROSAT_ENSURE(uf != NULL, "Failed to allocate union-find");
    uf->n = n;
    uf->parent = malloc((n+1)*sizeof(int));
    uf->rank   = calloc(n+1, sizeof(int));
    NITROSAT_ENSURE(uf->parent != NULL && uf->rank != NULL, "Failed to allocate union-find buffers");
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
    int    progress;   /* print milestone events to stderr (always on) */
    /* hot-path buffers */
    double *grad_buffer;
    double *heat_mult_buffer;
    int    *sat_counts; /* number of satisfied literals per clause */
    /* entropy cache for performance */
    double *entropy_cache;
    double *entropy_x_cache;
    /* Cheat Code 5: previous x for incremental sat_counts */
    double *prev_x;
    int    grok_detected;

    /* Incremental Unsat Tracking (For O(Unsat) performance) */
    int    *unsat_list;
    int    *unsat_pos;
    int    unsat_sz;

    /* Topology Buffers (Shared across compute_topology calls) */
    UnionFind *uf_buffer;
    uint8_t   *bitset_buffer;
    int       *active_vars_buffer;
} NitroSat;

static void recompute_sat_counts(NitroSat *ns);

static int clause_sat_count(const NitroSat *ns, int c)
{
    int sc = 0;
    for (int j = ns->cl_offs[c]; j < ns->cl_offs[c+1]; ++j) {
        int lit = ns->cl_flat[j];
        int v = abs(lit);
        if ((lit > 0 && ns->x[v] > 0.5) || (lit < 0 && ns->x[v] <= 0.5)) sc++;
    }
    return sc;
}

static int literal_satisfied_by_value(int lit, int value)
{
    return (lit > 0) ? value : !value;
}

static int clause_sat_delta_for_flip(const NitroSat *ns, int c, int var,
                                     int old_value, int new_value)
{
    int delta = 0;
    for (int j = ns->cl_offs[c]; j < ns->cl_offs[c+1]; ++j) {
        int lit = ns->cl_flat[j];
        if (abs(lit) != var) continue;
        delta += literal_satisfied_by_value(lit, new_value) -
                 literal_satisfied_by_value(lit, old_value);
    }
    return delta;
}

static void update_local_unsat_status(int c, int was_sat, int is_sat,
                                      int *unsat_list, int *unsat_pos, int *unsz)
{
    if (!was_sat && is_sat) {
        int pos = unsat_pos[c];
        if (pos != -1) {
            --(*unsz);
            unsat_list[pos] = unsat_list[*unsz];
            unsat_pos[unsat_list[pos]] = pos;
            unsat_pos[c] = -1;
        }
    } else if (was_sat && !is_sat) {
        if (unsat_pos[c] == -1) {
            unsat_pos[c] = *unsz;
            unsat_list[(*unsz)++] = c;
        }
    }
}

static TopologyInfo compute_topology(NitroSat *ns, double th)
{
    (void)th;
    int num_vars = ns->num_vars;
    const int *cl_offs = ns->cl_offs;
    const int *cl_flat = ns->cl_flat;

    TopologyInfo topo = {0};
    UnionFind *uf = ns->uf_buffer;
    uf->count = num_vars;
    /* Only reset UF for active variables (O(Unsat) reset) */
    for (int i = 0; i < ns->unsat_sz; ++i) {
        int c = ns->unsat_list[i];
        for (int j = cl_offs[c]; j < cl_offs[c+1]; ++j) {
            int v = abs(cl_flat[j]);
            uf->parent[v] = v; uf->rank[v] = 0;
        }
    }
    
    int *is_active = ns->active_vars_buffer; /* Used as a work list here */
    uint8_t *mark = ns->bitset_buffer;
    memset(mark, 0, (num_vars / 8) + 1);
    int active_cnt = 0;

    /* 1. O(Unsat) Beta-0 Identification */
    for (int i = 0; i < ns->unsat_sz; ++i) {
        int c = ns->unsat_list[i];
        int s = cl_offs[c], e = cl_offs[c+1];
        int first_v = -1;
        for (int j=s; j<e; ++j) {
            int v = abs(cl_flat[j]);
            if (!GET_BIT(mark, v)) {
                SET_BIT(mark, v);
                is_active[active_cnt++] = v;
            }
            if (first_v == -1) first_v = v;
            else uf_union(uf, first_v, v);
        }
    }

    /* 2. O(Unsat * K^2) Beta-1 Identification (Optimized Edge Count) */
    long long edge_cnt = 0;
    for (int i = 0; i < active_cnt; ++i) {
        int u = is_active[i];
        for (int k = ns->v2c_ptr[u]; k < ns->v2c_ptr[u+1]; ++k) {
            int c = ns->v2c_data[k];
            if (ns->sat_counts[c] > 0) continue;
            for (int m = cl_offs[c]; m < cl_offs[c+1]; ++m) {
                int v = abs(cl_flat[m]);
                if (v > u) edge_cnt++; 
            }
        }
    }
    /* Correction: The above counts edges multiply. Correcting for k-clique effects in Beta-1 proxy */
    edge_cnt /= 2; 

    topo.active_vars = active_cnt;
    topo.beta0 = 0;
    /* Count components among active vars only */
    for (int i=0; i<active_cnt; ++i) {
        int v = is_active[i];
        if (uf->parent[v] == v) topo.beta0++;
    }

    int chi = active_cnt - (int)edge_cnt + ns->unsat_sz;
    topo.beta1 = (topo.beta0 - chi > 0) ? (topo.beta0 - chi) : 0;
    topo.complexity_score = (topo.beta0 > 0 ? (double)topo.beta1 / topo.beta0 : 0.0);

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

    /* BFS 2-coloring to detect parity inconsistencies.  This signal is useful
       diagnostically, but a binary CNF clause is not generally an XOR equation.
       Do not emit a DRAT/LRAT unit from this check. */
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

    (void)pf;
    (void)derived_units;
    (void)contra_var;

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

    /* Phase 2: parity diagnostics only.  CNF binary clauses are implications,
       not XOR equations, so this must not generate proof clauses. */
    (void)proof_gf2_parity_check(inst, assign, NULL, NULL);

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
        double W_max = 0.0;
        for (int c = 0; c < ns->num_clauses; ++c)
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

    if (target->num_vars > proof_max_vars || (long)target->num_clauses > proof_max_clauses) {
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
            (long)inst->num_clauses <= proof_max_clauses) {
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
        free(ns->opt->m_standard);
        free(ns->opt->m_phase);
        free(ns->opt->v_amp);
        free(ns->opt->v_phase);
        free(ns->opt->prev_phase);
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
    free(ns->unsat_list);
    free(ns->unsat_pos);
    if (ns->uf_buffer) {
        free(ns->uf_buffer->parent);
        free(ns->uf_buffer->rank);
        free(ns->uf_buffer);
    }
    free(ns->bitset_buffer);
    free(ns->active_vars_buffer);
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
    if (nv <= 0) return;
    double *v = malloc((nv+1) * sizeof(double));
    double *y = malloc((nv+1) * sizeof(double));
    NITROSAT_ENSURE(v != NULL && y != NULL, "Failed to allocate spectral initialization buffers");

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
    NITROSAT_ENSURE(ns != NULL, "Failed to allocate solver");
    ns->inst       = inst;
    ns->cl_offs    = inst->cl_offs;
    ns->cl_flat    = inst->cl_flat;
    ns->num_vars   = inst->num_vars;
    ns->num_clauses= inst->num_clauses;
    ns->max_steps  = max_steps;
    ns->verbose    = verbose;
    ns->progress   = 1;

    /* degree vector (heat kernel) */
    compute_degrees(ns);

    /* continuous assignment – spectral initialisation (power iteration) */
    ns->x = malloc((ns->num_vars+1) * sizeof(double));
    NITROSAT_ENSURE(ns->x != NULL, "Failed to allocate assignment vector");
    spectral_init(ns, 50);  /* 50 iterations, matches gh.lua */

    /* clause weight initialisation (prime‑weighted) */
    ns->primes = sieve_primes(ns->num_clauses, &ns->prime_cnt);
    ns->cl_weights = malloc(ns->num_clauses * sizeof(double));
    ns->zeta_weights = malloc(ns->num_clauses * sizeof(double));
    NITROSAT_ENSURE(ns->num_clauses == 0 || (ns->cl_weights != NULL && ns->zeta_weights != NULL),
                    "Failed to allocate clause weights");
    for (int c=0;c<ns->num_clauses;++c) {
        int p = (c < ns->prime_cnt) ? ns->primes[c] : 2;
        ns->cl_weights[c]   = 1.0 / pow(1.0 + log((double)p), 1.0);   /* prime_alpha = 1.0 */
        ns->zeta_weights[c] = log((double)p) / (double)p;
    }
    ns->is_hard = malloc(ns->num_clauses);
    NITROSAT_ENSURE(ns->num_clauses == 0 || ns->is_hard != NULL, "Failed to allocate hard-clause flags");
    for (int c=0;c<ns->num_clauses;++c) ns->is_hard[c] = 1;   /* all hard initially, matches gh.lua */

    /* decimation vector */
    ns->decimated = calloc(ns->num_vars+1, 1);
    NITROSAT_ENSURE(ns->decimated != NULL, "Failed to allocate decimation vector");
    ns->grok_detected = 0;

    /* Incremental Unsat Tracking Initialization */
    ns->unsat_list = malloc(ns->num_clauses * sizeof(int));
    ns->unsat_pos = malloc(ns->num_clauses * sizeof(int));
    NITROSAT_ENSURE(ns->num_clauses == 0 || (ns->unsat_list != NULL && ns->unsat_pos != NULL),
                    "Failed to allocate unsat tracking buffers");
    for (int i=0; i<ns->num_clauses; ++i) ns->unsat_pos[i] = -1;
    ns->unsat_sz = 0;

    /* Topology Buffers */
    ns->uf_buffer = uf_new(ns->num_vars);
    ns->bitset_buffer = malloc((ns->num_vars / 8) + 8);
    ns->active_vars_buffer = malloc((ns->num_vars + 1) * sizeof(int));
    NITROSAT_ENSURE(ns->uf_buffer != NULL && ns->bitset_buffer != NULL && ns->active_vars_buffer != NULL,
                    "Failed to allocate topology buffers");

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
    NITROSAT_ENSURE(ns->grad_buffer != NULL && ns->heat_mult_buffer != NULL &&
                    (ns->num_clauses == 0 || ns->sat_counts != NULL) &&
                    ns->entropy_cache != NULL && ns->entropy_x_cache != NULL,
                    "Failed to allocate hot-path buffers");
    for (int i = 1; i <= ns->num_vars; ++i) {
        ns->entropy_x_cache[i] = -999.0; /* force initial compute */
    }
    /* Cheat Code 5: Initialize prev_x for incremental sat_counts */
    ns->prev_x = malloc((ns->num_vars + 1) * sizeof(double));
    NITROSAT_ENSURE(ns->prev_x != NULL, "Failed to allocate previous assignment buffer");
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
        if (k == 0) {
            ++unsat;
            continue;
        }

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

        /* Correct Inverse Metric Tensor (Inverted Poincaré Disk) 
           Metric g = 1/(1-z^2)^2, so Inverse Metric g^-1 = (1-z^2)^2
           This prevents binary collapse (x=0,1) and allows full momentum at x=0.5. */
        double z = 2.0 * ns->x[i] - 1.0;
        double metric_inv = (1.0 - z * z);
        grad[i] *= (metric_inv * metric_inv);
    }
    return unsat;
}

/* --------------------------------------------------------------------
   12. Riemann Zeta Function and Derivatives
       Used for the self-consistent annealing schedule
-------------------------------------------------------------------- */

/* Riemann zeta function ζ(s) for s > 1 using Euler-Maclaurin summation
   This is the generating function of the prime-weighted clause system */
static double riemann_zeta(double s)
{
    if (s <= 1.0) {
        /* Handle pole at s=1 - return large value to prevent division by zero */
        return 1e300;
    }

    /* Use Euler-Maclaurin for better convergence */
    double n = 1.0;
    double sum = 0.0;
    double term;

    /* Direct summation for first 1000 terms */
    for (int i = 0; i < 1000; i++) {
        term = pow(n, -s);
        sum += term;
        n += 1.0;
        if (term < 1e-15) break;
    }

    /* Add Euler-Maclaurin correction for remainder */
    double rm = pow(n, -s) / (s - 1.0);
    rm += 0.5 * pow(n, -s);
    rm += pow(n, -s-1) / 12.0;

    return sum + rm;
}

/* Derivative of Riemann zeta ζ'(s) = d/ds ζ(s) */
static double riemann_zeta_prime(double s)
{
    if (s <= 1.0) {
        return -1e300;  /* Pole behavior */
    }

    double n = 1.0;
    double sum = 0.0;
    double term;

    /* Sum -ln(n) * n^(-s) */
    for (int i = 0; i < 1000; i++) {
        term = -log(n) * pow(n, -s);
        sum += term;
        n += 1.0;
        if (fabs(term) < 1e-15) break;
    }

    return sum;
}

/* --------------------------------------------------------------------
   12b. Learning Rate Annealing - Oscillating Schedule (V2)
       A(t) = A_0 * max(0, I(t))
       where I(t) is a sum of four critically damped oscillators
-------------------------------------------------------------------- */
static double annealing_lr_v2(double t, double A0)
{
    double phi = PHI;
    /* Four incommensurate frequencies: π/20, π/10, π/9, 1/φ² */
    double I = sin(PI/20.0 * t) * exp(-PI/20.0 * t)
             + sin(PI/10.0 * t) * exp(-PI/10.0 * t)
             + sin(PI/9.0  * t) * exp(-PI/9.0  * t)
             + sin((t / (phi * phi))) * exp(-(t / (phi * phi)));

    /* Clamp to non-negative (max(0, I(t))) */
    if (I < 0.0) I = 0.0;

    return A0 * I;
}

/* --------------------------------------------------------------------
   12c. Temperature from Riemann Zeta Dynamics (V2)
       Self-consistent annealing: β_eff uses ζ(s) to sweep from s=1.5 → 1+
       This anneals along the natural spectral evolution of the energy landscape
-------------------------------------------------------------------- */
static double compute_beta_eff_v2(double beta, int step, int max_steps, double s_cumulative, double *s_out)
{
    double tau = max_steps / 4.0;
    double epsilon = 0.5;
    double s = 1.0 + epsilon * exp(-step / tau) + s_cumulative;
    if (s > 1.5) s = 1.5;
    *s_out = s;

    /* Pole-aware approximation for fast execution: ζ(s) ≈ 1/(s-1) + γ */
    double zeta_s = 1.0 / (s - 1.0) + GAMMA;
    double zeta_prime_s = -1.0 / ((s-1.0) * (s-1.0));
    
    double zeta_dynamics = (zeta_s + 0.1 * zeta_prime_s) / fmax(fabs(zeta_s), 0.1);
    return beta * (1.0 + zeta_dynamics);
}

/* --------------------------------------------------------------------
   13. Satisfaction test (discrete interpretation)
-------------------------------------------------------------------- */
static void recompute_sat_counts(NitroSat *ns)
{
    ns->unsat_sz = 0;
    for (int c = 0; c < ns->num_clauses; ++c) {
        int sc = clause_sat_count(ns, c);
        ns->sat_counts[c] = sc;
        if (sc == 0) {
            ns->unsat_pos[c] = ns->unsat_sz;
            ns->unsat_list[ns->unsat_sz++] = c;
        } else {
            ns->unsat_pos[c] = -1;
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
   17. BAHA-WalkSAT – Adelic Bit-Parallel Discrete Local Search
        "Secret Math Technique": Bit-Parallel Adelic Evaluation (BPAE)
        We evaluate 64 parallel truth assignments in one pass over the clauses.
-------------------------------------------------------------------- */
static int baha_walksat(NitroSat *ns, int max_flips)
{
    const int num_v = ns->num_vars;
    const int num_c = ns->num_clauses;
    /* Calculate jump interval for roughly 20 global re-evaluations */
    int BP_JUMP_INTERVAL = max_flips / 20;
    if (BP_JUMP_INTERVAL < 500) BP_JUMP_INTERVAL = 500;

    recompute_sat_counts(ns);
    
    /* Best state tracking */
    int *best_assignment = malloc((num_v + 1) * sizeof(int));
    for (int i=1; i<=num_v; ++i) best_assignment[i] = (ns->x[i] > 0.5);
    int current_unsat = 0;
    for (int c=0; c<num_c; ++c) if (ns->sat_counts[c] == 0) current_unsat++;
    int best_unsat = current_unsat;

    /* Bit-parallel state (64 parallel worlds) */
    uint64_t *state64 = malloc((num_v + 1) * sizeof(uint64_t));
    
    /* Incremental SAT list management */
    int *unsat_list = malloc(num_c * sizeof(int));
    int *unsat_pos = malloc(num_c * sizeof(int));
    for (int c=0; c<num_c; ++c) unsat_pos[c] = -1;
    int unsz = 0;
    for (int c=0; c<num_c; ++c) {
        if (ns->sat_counts[c] == 0) {
            unsat_pos[c] = unsz;
            unsat_list[unsz++] = c;
        }
    }

    for (int step = 1; step <= max_flips; ++step) {
        if (unsz == 0) break;

        /* ADELIC JUMP: Bit-Parallel Evaluation of 64 trajectories (Miles Faster!) */
        if (step % BP_JUMP_INTERVAL == 0) {
            for (int i=1; i<=num_v; ++i) {
                uint64_t bits = 0;
                int current = best_assignment[i];
                for (int b=0; b<64; ++b) {
                    if (rng_double() < 0.2) bits |= ((uint64_t)(1-current) << b);
                    else bits |= ((uint64_t)current << b);
                }
                state64[i] = bits;
            }

            uint32_t scores[64] = {0};
            for (int c=0; c<num_c; ++c) {
                uint64_t sat_mask = 0;
                for (int j=ns->cl_offs[c]; j<ns->cl_offs[c+1]; ++j) {
                    int lit = ns->cl_flat[j];
                    uint64_t v_bits = state64[abs(lit)];
                    sat_mask |= (lit > 0) ? v_bits : ~v_bits;
                }
                uint64_t unsat_mask = ~sat_mask;
                /* Efficient bit-counting: only iterate set bits (worlds where clause is unsat) */
                while (unsat_mask) {
                    int b = __builtin_ctzll(unsat_mask);
                    scores[b]++;
                    unsat_mask &= (unsat_mask - 1);
                }
            }

            int found_jump = 0;
            for (int b=0; b<64; ++b) {
                if (scores[b] < (uint32_t)best_unsat) {
                    best_unsat = scores[b];
                    for (int i=1; i<=num_v; ++i) best_assignment[i] = (int)((state64[i] >> b) & 1);
                    found_jump = 1;
                }
            }
            if (found_jump) {
                for (int i=1; i<=num_v; ++i) ns->x[i] = (double)best_assignment[i];
                recompute_sat_counts(ns);
                unsz = 0;
                for (int c=0; c<num_c; ++c) {
                    if (ns->sat_counts[c] == 0) {
                        unsat_pos[c] = unsz;
                        unsat_list[unsz++] = c;
                    } else unsat_pos[c] = -1;
                }
            }
        }

        /* LOCAL SEARCH: Greedy with Adelic Gradient Tie-Breaking */
        int target_c = unsat_list[rng_int(0, unsz - 1)];
        int c_start = ns->cl_offs[target_c];
        int c_size = ns->cl_offs[target_c+1] - c_start;
        
        int best_var = -1, best_break = 1e9;
        for (int k=0; k<c_size; ++k) {
            int v = abs(ns->cl_flat[c_start + k]);
            if (ns->decimated[v]) continue;
            int brk = 0;
            for (int p = ns->v2c_ptr[v]; p < ns->v2c_ptr[v+1]; ++p) {
                int cl = ns->v2c_data[p];
                if (ns->sat_counts[cl] == 1) {
                    int sc = 0;
                    for (int j=ns->cl_offs[cl]; j<ns->cl_offs[cl+1]; ++j) {
                        int l = ns->cl_flat[j];
                        if ((l>0 && ns->x[abs(l)]>0.5) || (l<0 && ns->x[abs(l)]<=0.5)) sc++;
                    }
                    if (sc == 1) {
                        /* Check if v is THE literal that satisfies it */
                        int lit_p = 0;
                        for (int j=ns->cl_offs[cl]; j<ns->cl_offs[cl+1]; ++j) {
                            if (abs(ns->cl_flat[j]) == v) {
                                int l = ns->cl_flat[j];
                                if ((l>0 && ns->x[v]>0.5) || (l<0 && ns->x[v]<=0.5)) lit_p = 1;
                                break;
                            }
                        }
                        if (lit_p) brk++;
                    }
                }
            }
            /* Score biasing using Adelic continuous gradient */
            double force = ns->grad_buffer[v] * (ns->x[v] > 0.5 ? 1.0 : -1.0);
            if (brk < best_break || (brk == best_break && force > 0)) {
                best_break = brk; best_var = v;
            }
        }

        if (best_var != -1 && (best_break == 0 || rng_double() < 0.3)) {
            int old_at_v = (ns->x[best_var] > 0.5);
            int new_at_v = 1 - old_at_v;
            ns->x[best_var] = (double)new_at_v;
            int last_cl = -1;
            for (int p = ns->v2c_ptr[best_var]; p < ns->v2c_ptr[best_var+1]; ++p) {
                int cl = ns->v2c_data[p];
                if (cl == last_cl) continue;
                last_cl = cl;
                int was_sat = (ns->sat_counts[cl] > 0);
                ns->sat_counts[cl] += clause_sat_delta_for_flip(ns, cl, best_var, old_at_v, new_at_v);
                int is_sat = (ns->sat_counts[cl] > 0);
                update_local_unsat_status(cl, was_sat, is_sat, unsat_list, unsat_pos, &unsz);
            }
        }
        if (unsz < best_unsat) {
            best_unsat = unsz;
            for (int i=1; i<=num_v; ++i) best_assignment[i] = (ns->x[i] > 0.5);
        }
    }
    for (int i=1; i<=num_v; ++i) ns->x[i] = (double)best_assignment[i];
    recompute_sat_counts(ns);
    free(state64); free(best_assignment); free(unsat_list); free(unsat_pos);
    return (best_unsat == 0);
}

/* --------------------------------------------------------------------
   18. Core Decomposition (Blasting Phase)
-------------------------------------------------------------------- */
static int core_decomposition(NitroSat *ns)
{
    recompute_sat_counts(ns);
    int *core = malloc(ns->unsat_sz * sizeof(int));
    int core_sz = 0;
    for (int i=0; i<ns->unsat_sz; ++i) core[core_sz++] = ns->unsat_list[i];
    if (core_sz == 0) { free(core); return 1; }

    double *saved = malloc(ns->num_clauses * sizeof(double));
    memcpy(saved, ns->cl_weights, ns->num_clauses * sizeof(double));
    for (int i=0; i<core_sz; ++i) ns->cl_weights[core[i]] *= 10.0;

    int max_flips = ns->num_vars * 10;
    if (max_flips > 1000000) max_flips = 1000000;
    int success = 0;

    for (int flips=0; flips < max_flips; ++flips) {
        int target_c = core[rand() % core_sz];
        if (ns->sat_counts[target_c] > 0) continue;
        int cs = ns->cl_offs[target_c], ce = ns->cl_offs[target_c + 1];
        if (ce <= cs) continue;
        int var = abs(ns->cl_flat[cs + (rand() % (ce - cs))]);
        if (ns->decimated[var]) continue;

        ns->x[var] = (ns->x[var] > 0.5) ? 0.0 : 1.0;
        for (int p = ns->v2c_ptr[var]; p < ns->v2c_ptr[var+1]; ++p) {
            int cl = ns->v2c_data[p];
            int sc = 0;
            for (int j=ns->cl_offs[cl]; j < ns->cl_offs[cl+1]; ++j) {
                int l = ns->cl_flat[j];
                if ((l > 0 && ns->x[abs(l)] > 0.5) || (l < 0 && ns->x[abs(l)] <= 0.5)) sc++;
            }
            ns->sat_counts[cl] = sc;
        }
    }
    recompute_sat_counts(ns);
    int uns = 0;
    for (int i=0; i<core_sz; ++i) if (ns->sat_counts[core[i]] == 0) uns++;
    success = (uns == 0);
    memcpy(ns->cl_weights, saved, ns->num_clauses * sizeof(double));
    free(core); free(saved);
    return success;
}

static double score_branch(NitroSat *ns, double beta, int n_samples, double *best_x_out)
{
    if (beta <= 0) return -1e300;
    double total_score = 0.0, best_seen = 1e300;
    double *buf = malloc((ns->num_vars + 1) * sizeof(double));
    for (int s = 0; s < n_samples; ++s) {
        for (int i = 1; i <= ns->num_vars; ++i) buf[i] = ns->decimated[i] ? ns->x[i] : (rng_double() > 0.5);
        double energy = 0.0;
        for (int c = 0; c < ns->num_clauses; ++c) {
            int sat = 0;
            for (int q = ns->cl_offs[c]; q < ns->cl_offs[c+1]; ++q) {
                int lit = ns->cl_flat[q], v = abs(lit);
                if ((lit > 0 && buf[v] > 0.5) || (lit < 0 && buf[v] <= 0.5)) { sat = 1; break; }
            }
            if (!sat) energy += ns->cl_weights[c];
        }
        total_score += exp(-beta * energy);
        if (energy < best_seen) {
            best_seen = energy;
            if (best_x_out) for (int i = 1; i <= ns->num_vars; ++i) best_x_out[i] = buf[i];
        }
    }
    free(buf);
    return total_score / n_samples + 1.0 / (best_seen + 1.0);
}

/* --------------------------------------------------------------------
   Phase 2: Discrete Greedy Repair (WalkSAT-style, O(1) per flip)
   
   The old continuous-diffusion approach caused catastrophic Jacobi
   oscillations on large instances. This version operates purely in
   discrete space: pick random unsat clause → find min-break variable
   → flip it → incrementally update sat_counts. Provably stable.
-------------------------------------------------------------------- */
static int topological_repair_phase(NitroSat *ns, int max_steps)
{
    if (!ns->use_topology) return 0;

    /* Skip topological repair for dense random 5-SAT: the clause-variable
       graph is essentially a dense random graph with little exploitable
       topology.  Heat kernel propagation becomes expensive and ineffective.
       Density > 4.0 means avg clause size * clause ratio exceeds WalkSAT
       crossover territory — bail out early. */
    double density = (double)ns->num_clauses / ns->num_vars;
    if (density > 4.0) {
        if (ns->verbose) fprintf(stderr, "[topo] skipping: dense formula (density=%.1f > 4.0)\n", density);
        return 0;
    }

    /* Snap continuous x to discrete 0/1 for stable repair */
    for (int i = 1; i <= ns->num_vars; ++i)
        ns->x[i] = (ns->x[i] > 0.5) ? 1.0 : 0.0;

    recompute_sat_counts(ns);
    if (ns->unsat_sz == 0) return 1;

    /* Build incremental unsat list for O(1) random access */
    int *ulist = malloc(ns->num_clauses * sizeof(int));
    int *upos  = malloc(ns->num_clauses * sizeof(int));
    int usz = 0;
    for (int c = 0; c < ns->num_clauses; ++c) {
        if (ns->sat_counts[c] == 0) {
            upos[c] = usz;
            ulist[usz++] = c;
        } else {
            upos[c] = -1;
        }
    }

    int best_unsat = usz;
    double *best_x = malloc((ns->num_vars + 1) * sizeof(double));
    for (int i = 1; i <= ns->num_vars; ++i) best_x[i] = ns->x[i];

    int max_flips = max_steps * ns->num_vars;
    /* Scale flip budget by clause density: moderate density benefits from
       more attempts; very high density is already caught by the early skip. */
    if (density > 2.0) max_flips = (int)(max_flips * density / 2.0);
    if (max_flips > 10000000) max_flips = 10000000;
    if (max_flips < 200000) max_flips = 200000;

    /* Precompute heat kernel normalization */
    double mean_deg = 0;
    for (int i = 1; i <= ns->num_vars; ++i) mean_deg += ns->degrees[i];
    mean_deg /= ns->num_vars;
    double hk_beta = ns->heat_beta / (mean_deg + 1.0);

    for (int flip = 0; flip < max_flips && usz > 0; ++flip) {

        /* Pick a random unsatisfied clause */
        int ci = rng_int(0, usz - 1);
        int target_c = ulist[ci];

        /* Find the variable with minimum break score */
        int cs = ns->cl_offs[target_c];
        int ce = ns->cl_offs[target_c + 1];
        int best_var = -1, best_break = 1 << 30;

        for (int j = cs; j < ce; ++j) {
            int v = abs(ns->cl_flat[j]);
            if (ns->decimated[v]) continue;

            /* Count how many currently-satisfied clauses this flip would break */
            int brk = 0;
            for (int p = ns->v2c_ptr[v]; p < ns->v2c_ptr[v+1]; ++p) {
                int cl = ns->v2c_data[p];
                if (ns->sat_counts[cl] == 1) {
                    /* Check if v is the sole satisfier */
                    int pol = 0;
                    for (int q = ns->cl_offs[cl]; q < ns->cl_offs[cl+1]; ++q) {
                        if (abs(ns->cl_flat[q]) == v) {
                            pol = (ns->cl_flat[q] > 0);
                            break;
                        }
                    }
                    int v_satisfies = (pol == (int)(ns->x[v] > 0.5));
                    if (v_satisfies) brk++;
                }
            }

            /* Use Adelic gradient as tie-breaker */
            if (brk < best_break || (brk == best_break && best_var != -1 &&
                fabs(ns->grad_buffer[v]) > fabs(ns->grad_buffer[best_var]))) {
                best_break = brk;
                best_var = v;
            }
        }

        if (best_var == -1) continue;

        /* Only flip if zero-break, or with noise probability */
        if (best_break > 0 && rng_double() > 0.4) continue;

        /* Flip the variable */
        int old_val = (int)(ns->x[best_var] > 0.5);
        ns->x[best_var] = (double)(1 - old_val);
        int new_val = 1 - old_val;

        /* Incrementally update sat_counts and unsat list */
        int last_cl = -1;
        for (int p = ns->v2c_ptr[best_var]; p < ns->v2c_ptr[best_var+1]; ++p) {
            int cl = ns->v2c_data[p];
            if (cl == last_cl) continue;
            last_cl = cl;
            int was_sat = (ns->sat_counts[cl] > 0);
            ns->sat_counts[cl] += clause_sat_delta_for_flip(ns, cl, best_var, old_val, new_val);
            int is_sat = (ns->sat_counts[cl] > 0);
            update_local_unsat_status(cl, was_sat, is_sat, ulist, upos, &usz);
        }

        /* ── HEAT KERNEL PROPAGATION ──
           After flipping best_var, propagate a heat wave through the
           local neighborhood. For each clause still unsatisfied after
           the flip, nudge neighboring variables toward resolution,
           weighted by exp(-β_eff · degree). If they cross 0.5, cascade. */
        for (int p = ns->v2c_ptr[best_var]; p < ns->v2c_ptr[best_var+1]; ++p) {
            int cl = ns->v2c_data[p];
            if (ns->sat_counts[cl] > 0) continue; /* only propagate from still-unsat clauses */

            for (int q = ns->cl_offs[cl]; q < ns->cl_offs[cl+1]; ++q) {
                int neighbor = abs(ns->cl_flat[q]);
                if (neighbor == best_var || ns->decimated[neighbor]) continue;

                /* Heat kernel weight: normalized by mean degree */
                double hk_weight = exp(-hk_beta * ns->degrees[neighbor]);

                /* Direction: which way does this clause want the neighbor? */
                double desired = (ns->cl_flat[q] > 0) ? 1.0 : 0.0;
                double nudge = hk_weight * 0.3 * (desired - ns->x[neighbor]);

                double old_x = ns->x[neighbor];
                ns->x[neighbor] += nudge;
                /* Clamp to [0,1] */
                if (ns->x[neighbor] > 1.0) ns->x[neighbor] = 1.0;
                if (ns->x[neighbor] < 0.0) ns->x[neighbor] = 0.0;

                /* Cascade: if crossing 0.5 threshold, treat as discrete flip */
                int was = (old_x > 0.5);
                int now = (ns->x[neighbor] > 0.5);
                if (was != now) {
                    /* Snap to discrete */
                    ns->x[neighbor] = now ? 1.0 : 0.0;
                    /* Incremental SAT update for cascade flip */
                    int last_cl2 = -1;
                    for (int r = ns->v2c_ptr[neighbor]; r < ns->v2c_ptr[neighbor+1]; ++r) {
                        int cl2 = ns->v2c_data[r];
                        if (cl2 == last_cl2) continue;
                        last_cl2 = cl2;
                        int was_sat2 = (ns->sat_counts[cl2] > 0);
                        ns->sat_counts[cl2] += clause_sat_delta_for_flip(ns, cl2, neighbor, was, now);
                        int is_sat2 = (ns->sat_counts[cl2] > 0);
                        update_local_unsat_status(cl2, was_sat2, is_sat2, ulist, upos, &usz);
                    }
                }
            }
        }

        /* Track best */
        if (usz < best_unsat) {
            best_unsat = usz;
            for (int i = 1; i <= ns->num_vars; ++i) best_x[i] = ns->x[i];
        }

        /* Progress reporting */
        if (ns->verbose && (flip % 50000 == 0)) {
            printf("[phase‑2] flip %d: unsat=%d best=%d\n", flip, usz, best_unsat);
            fflush(stdout);
        }
    }

    /* Restore best assignment found */
    for (int i = 1; i <= ns->num_vars; ++i) ns->x[i] = best_x[i];
    recompute_sat_counts(ns);

    if (ns->verbose) {
        printf("[phase‑2] finished, sat=%d/%d\n",
               ns->num_clauses - ns->unsat_sz, ns->num_clauses);
    }

    free(ulist); free(upos); free(best_x);
    return (ns->unsat_sz == 0);
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
            int new_v_sat = 1 - old_v_sat;
            ns->x[best_v] = (double)new_v_sat;
            int last_cl = -1;
            for (int p = ns->v2c_ptr[best_v]; p < ns->v2c_ptr[best_v+1]; ++p) {
                int cl = ns->v2c_data[p];
                if (cl == last_cl) continue;
                last_cl = cl;
                ns->sat_counts[cl] += clause_sat_delta_for_flip(ns, cl, best_v, old_v_sat, new_v_sat);
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
    
    double s_cumulative = 0.0;
    int total_triggers = 0;
    double last_rho = 1.0;
    int last_trigger_step = -1000;
    
    #define CORR_WINDOW 100
    double heat_buf[CORR_WINDOW] = {0};
    double prime_buf[CORR_WINDOW] = {0};
    int buf_n = 0;

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

        /* STAGE 1 HANDOFF: Detect if we've captured the global basin. 
           Wait at least 50 steps to ensure the Adelic Gradient has geometric signal. */
        int is_grok = (step > 50 && prev_sat > 0 && (sat - prev_sat) > (int)(0.5 * ns->num_clauses));
        if (is_grok || (step > 50 && sat >= (int)(0.95 * ns->num_clauses))) {
            if (ns->verbose || ns->progress) {
                if (is_grok) fprintf(stderr, "[GROK] Detected rapid transition! %d%% jump at step %d. Handoff to polishing.\n", 
                                     (int)(100.0*(sat - prev_sat)/ns->num_clauses), step);
                else fprintf(stderr, "[%02d%%] Handoff to polishing at step %d.\n", (int)(100.0*sat/ns->num_clauses), step);
            }
            ns->grok_detected = 1;
            break;
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
            
            /* If slope is flat (<0.01% gain) and satisfaction is high (>99%) and not recently fired */
            if (step > 200 && fabs(slope) < 0.0001 && sat > (int)(0.99 * ns->num_clauses) && (step > last_stagnation_step + 100)) {
                if (ns->verbose) {
                    printf("[STAGNATION] Slope plateaued at %.4f%%. Proactively jumping to Stage 2 finisher!\n", current_avg*100.0);
                } else if (ns->progress) {
                    fprintf(stderr, "  [stagnation] %.2f%% plateau → jumping to finisher\n", current_avg*100.0);
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

        /* NitroSATv2: New oscillating LR schedule and zeta temperature dynamics */
        double t = step * (30.0 / max_steps);
        double lr = annealing_lr_v2(t, lr0);  /* V2: A(t) = A_0 * max(0, I(t)) */
        if (lr <= 0.0) lr = 1e-6;
        ns->opt->lr = lr;

        /* V2: Compute effective temperature from Riemann zeta dynamics
           This anneals along the natural spectral evolution of the energy landscape */
        double base_beta = 1.0;  /* Base inverse temperature */
        double s_val = 0.0;
        double beta_eff = compute_beta_eff_v2(base_beta, step, max_steps, s_cumulative, &s_val);

        /* Use beta_eff for entropy weight modulation (cooling toward zero) */
        ns->entropy_weight = 0.01 * beta_eff;

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

        /* Zeta Thermostat: Correlation order parameter tracking & Glass transition decoupling detection */
        double e_heat_val = 0.0;
        double p_prime_val = 0.0;
        double inv_zeta2 = 6.0 / (PI * PI);
        
        /* Efficient order parameter update (every 10 steps to save cycles) */
        if (step % 10 == 0) {
            for (int i = 1; i <= ns->num_vars; ++i) {
                e_heat_val += exp(-ns->heat_beta * ns->degrees[i]);
            }
            for (int c = 0; c < ns->num_clauses; ++c) {
                double vc = 1.0;
                for (int i = ns->cl_offs[c]; i < ns->cl_offs[c+1]; ++i) {
                    int lit = ns->cl_flat[i];
                    double val = ns->x[abs(lit)];
                    vc *= (lit > 0) ? (1.0 - val) : val;
                }
                p_prime_val -= log(inv_zeta2 + 1e-6 * vc);
            }
            int wpos = buf_n % CORR_WINDOW;
            heat_buf[wpos] = e_heat_val;
            prime_buf[wpos] = p_prime_val;
            buf_n++;
        }
        int w_filled = (buf_n < CORR_WINDOW) ? buf_n : CORR_WINDOW;
        
        /* Betti-1 Targeted Harmonic Perturbation (every 100 steps) */
        if (step % 100 == 0) {
            /* Compute Pearson */
            if (w_filled >= 10) {
                double mx = 0, my = 0;
                for (int i = 0; i < w_filled; ++i) { mx += heat_buf[i]; my += prime_buf[i]; }
                mx /= w_filled; my /= w_filled;
                double cov = 0, vx = 0, vy = 0;
                for (int i = 0; i < w_filled; ++i) {
                    double dx = heat_buf[i] - mx;
                    double dy = prime_buf[i] - my;
                    cov += dx * dy;
                    vx += dx * dx;
                    vy += dy * dy;
                }
                double denom = sqrt(vx * vy);
                last_rho = (denom < 1e-12) ? 1.0 : (cov / denom);
            }
            
            if (last_rho < 0.95 && (step - last_trigger_step) > 200) {
                last_trigger_step = step;
                total_triggers++;
                
                double bump = (1.0 - last_rho) * 0.05;
                s_cumulative += bump;
                if (s_cumulative > 0.5) s_cumulative = 0.5; /* never exceed max */

                if (ns->verbose) {
                    printf("[ZETASAT] step=%4d  TRIGGER #%d  rho=%.4f  s=%.4f  beta_eff=%.4f\n",
                           step, total_triggers, last_rho, s_val, beta_eff);
                }

                double *frust_w = calloc(ns->num_vars + 1, sizeof(double));
                for (int c = 0; c < ns->num_clauses; ++c) {
                    double vc = 1.0;
                    for (int i = ns->cl_offs[c]; i < ns->cl_offs[c+1]; ++i) {
                        int lit = ns->cl_flat[i];
                        int v = abs(lit);
                        double val = ns->x[v];
                        double lv = (lit > 0) ? (1.0 - val) : val;
                        vc *= lv;
                    }
                    if (vc > 0.5) { /* Frustration threshold */
                        for (int i = ns->cl_offs[c]; i < ns->cl_offs[c+1]; ++i) {
                            int v = abs(ns->cl_flat[i]);
                            frust_w[v] += vc;
                        }
                    }
                }
                
                double alpha = 0.02;
                double amplitude = 0.05;
                double decay = exp(-alpha * step);
                for (int j = 1; j <= ns->num_vars; ++j) {
                    double wj = frust_w[j] > 0.0 ? frust_w[j] : 1.0;
                    double w1 = PI * j / (20.0 * ns->num_vars);
                    double w2 = PI * j / (10.0 * ns->num_vars);
                    double w3 = j / (PHI * PHI * ns->num_vars);
                    
                    double p = amplitude * wj * decay * (
                        sin(w1 * step) * exp(-w1 * step) +
                        0.5 * sin(w2 * step) * exp(-w2 * step) +
                        0.25 * sin(w3 * step) * exp(-w3 * step)
                    );
                    if (!ns->decimated[j]) {
                        ns->x[j] = clamp(ns->x[j] + p, 0.0, 1.0);
                    }
                }
                free(frust_w);
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
        fprintf(stderr, "Freeing best_x in nitrosat_solve (early solve)... (sat: %d/%d)\n", sat, ns->num_clauses);
        free(best_x);
        return 1;
    }   /* already solved */

    /* Phase‑2 : topological repair (≈95 % satisfied) */
    if (sat >= (int)(0.95 * ns->num_clauses)) {
        if (ns->verbose) puts("[phase‑2] entering topological repair");
        else if (ns->progress) fprintf(stderr, "  [phase-2] topological repair... ");
        
        for (int i = 1; i <= ns->num_vars; ++i) ns->x[i] = best_x[i];
        topological_repair_phase(ns, 1500);
        int phase2_sat = check_satisfaction(ns);
        
        if (phase2_sat > best_sat) {
            best_sat = phase2_sat;
            for (int i = 1; i <= ns->num_vars; ++i) best_x[i] = ns->x[i];
        }
        if (ns->verbose) printf("[phase‑2] finished, sat=%d/%d\n", best_sat, ns->num_clauses);
        else if (ns->progress) fprintf(stderr, "done → %d/%d (%.2f%%)\n", best_sat, ns->num_clauses, 100.0*best_sat/ns->num_clauses);
        if (best_sat == ns->num_clauses) { free(best_x); return 1; }
    }

    /* Phase‑3 : adelic saturation (≈98 % satisfied) */
    if (best_sat >= (int)(0.98 * ns->num_clauses)) {
        if (ns->verbose) puts("[phase‑3] entering adelic saturation");
        else if (ns->progress) fprintf(stderr, "  [phase-3] adelic saturation... ");
        for (int i = 1; i <= ns->num_vars; ++i) ns->x[i] = best_x[i];
        
        adelic_saturation_phase(ns, 2000);
        int phase3_sat = check_satisfaction(ns);
        if (phase3_sat > best_sat) {
            best_sat = phase3_sat;
            for (int i = 1; i <= ns->num_vars; ++i) best_x[i] = ns->x[i];
        }
        if (ns->verbose) printf("[phase‑3] finished, sat=%d/%d\n", best_sat, ns->num_clauses);
        else if (ns->progress) fprintf(stderr, "done → %d/%d (%.2f%%)\n", best_sat, ns->num_clauses, 100.0*best_sat/ns->num_clauses);
        if (best_sat == ns->num_clauses) { free(best_x); return 1; }
    }

    /* Core decomposition “blast’’ */
    if (best_sat < ns->num_clauses) {
        if (ns->verbose) puts("[core] attempting core decomposition");
        else if (ns->progress) fprintf(stderr, "  [core]    decomposition... ");
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
        else if (ns->progress) fprintf(stderr, "  [final]  BAHA discrete search... ");
        for (int i = 1; i <= ns->num_vars; ++i) ns->x[i] = best_x[i] > 0.5 ? 1.0 : 0.0;
        
        /* Run BAHA-WalkSAT with baha.cpp params: 200 steps × 50 samples */
        /* Max flips capped to 2M for massive instances to prevent hangs */
        int max_sc_flips = ns->num_vars * 200;
        if (max_sc_flips > 2000000) max_sc_flips = 2000000;
        baha_walksat(ns, max_sc_flips);
        int bsat = check_satisfaction(ns);
        if (bsat > best_sat) {
            best_sat = bsat;
            for (int i = 1; i <= ns->num_vars; ++i) best_x[i] = ns->x[i];
        }
    }

    /* Always restore best known assignment before returning */
    for (int i = 1; i <= ns->num_vars; ++i) ns->x[i] = best_x[i];
    recompute_sat_counts(ns);
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
        else if (ns->progress) fprintf(stderr, "[pass %d/%d] ", pass, num_passes);
        
        success = nitrosat_solve(ns);
        recompute_sat_counts(ns);
        int sat_count = check_satisfaction(ns);
        if (ns->progress && !ns->verbose)
            fprintf(stderr, "sat=%d/%d (%.2f%%)\n", sat_count, ns->num_clauses, 100.0*sat_count/ns->num_clauses);
        
        if (sat_count > best_overall_sat) {
            best_overall_sat = sat_count;
            for (int i = 1; i <= ns->num_vars; ++i) best_overall_x[i] = ns->x[i];
        }
        
        if (sat_count == ns->num_clauses) { success = 1; break; }
        
        /* For truly massive instances (>50k vars), one grok event is often the limit.
           For small/medium instances, keep hunting via DCW passes. */
        if (ns->grok_detected && sat_count > (int)(0.99 * ns->num_clauses) && ns->num_vars > 50000) {
            if (ns->verbose) printf("[GROK] High satisfaction on large instance. Stopping.\n");
            break;
        }
        
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
        ns->grok_detected = 0;
        for (int i = 1; i <= ns->num_vars; ++i) {
            ns->opt->m_standard[i] = 0.0;
            ns->opt->m_phase[i] = 0.0;
            ns->opt->v_amp[i] = 0.0;
            ns->opt->v_phase[i] = 0.0;
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
    recompute_sat_counts(ns);
    int final_sat = check_satisfaction(ns);
    int unsat = ns->num_clauses - final_sat;
    double sat_rate = (ns->num_clauses > 0) ? (double)final_sat / ns->num_clauses : 1.0;

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

    /* Assignment as 0/1 binary array (index = variable number - 1) */
    printf("  \"assignment\": [");
    for (int i = 1; i <= ns->num_vars; ++i) {
        printf("%d", (ns->x[i] > 0.5) ? 1 : 0);
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
    double W_max = 0.0;
    for (int c = 0; c < ns->num_clauses; ++c)
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
    diag->aggregation_error = (ns->num_clauses > 1 && fabs(mean_w) > 1e-12) ?
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
        fprintf(stderr, "Usage: %s <cnf-file> [max-steps] [--no-dcw] [--no-topo] [--cinematic] [--proof <path>] [--proof-format drat] [--proof-max-vars N] [--proof-max-clauses N]\n", argv[0]);
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
        printf("  --cinematic              Output results in verbose text format (legacy)\n");
        printf("  --proof <path>          Generate DRAT proof on UNSAT (requires proof backend)\n");
        printf("  --proof-format drat    Proof format (default: drat)\n");
        printf("  --proof-max-vars N      Max variables for proof generation (default: 500000)\n");
        printf("  --proof-max-clauses N   Max clauses for proof generation (default: 3000000000)\n");
        printf("  --help, -h              Show this help message\n\n");
        printf("Repository: https://github.com/sethuiyer/NitroSAT\n");
        return EXIT_SUCCESS;
    }

    const char *cnf_file = argv[1];
    int max_steps = 3000;
    int use_dcw = 1;
    int use_topo = 1;
    int json_output = 1;
    const char *proof_path = NULL;
    const char *proof_format = "drat";

    for (int i = 2; i < argc; ++i) {
        if (strcmp(argv[i], "--no-dcw") == 0) use_dcw = 0;
        else if (strcmp(argv[i], "--no-topo") == 0) use_topo = 0;
        else if (strcmp(argv[i], "--cinematic") == 0) json_output = 0;
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
                fprintf(stderr, "--proof-format requires a value (drat).\n");
                return EXIT_FAILURE;
            }
            proof_format = argv[++i];
        }
        else if (strcmp(argv[i], "--proof-max-vars") == 0) {
            if (i + 1 >= argc) {
                fprintf(stderr, "--proof-max-vars requires a number.\n");
                return EXIT_FAILURE;
            }
            char *endptr = NULL;
            long parsed = strtol(argv[++i], &endptr, 10);
            if (endptr == argv[i] || *endptr != '\0' || parsed > INT_MAX) {
                fprintf(stderr, "--proof-max-vars requires a valid integer.\n");
                return EXIT_FAILURE;
            }
            proof_max_vars = (int)parsed;
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
            char *endptr = NULL;
            proof_max_clauses = strtol(argv[++i], &endptr, 10);
            if (endptr == argv[i] || *endptr != '\0') {
                fprintf(stderr, "--proof-max-clauses requires a valid integer.\n");
                return EXIT_FAILURE;
            }
            if (proof_max_clauses <= 0) {
                fprintf(stderr, "--proof-max-clauses must be positive.\n");
                return EXIT_FAILURE;
            }
        }
        else if (argv[i][0] == '-') {
            fprintf(stderr, "Unknown option: %s\n", argv[i]);
            return EXIT_FAILURE;
        }
        else {
            char *endptr = NULL;
            long parsed = strtol(argv[i], &endptr, 10);
            if (endptr == argv[i] || *endptr != '\0' || parsed <= 0 || parsed > INT_MAX) {
                fprintf(stderr, "max-steps must be a positive integer.\n");
                return EXIT_FAILURE;
            }
            max_steps = (int)parsed;
        }
    }

    if (strcmp(proof_format, "drat") != 0) {
        fprintf(stderr, "Unsupported proof format: %s (supported: drat)\n", proof_format);
        return EXIT_FAILURE;
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
    if (json_output)
        fprintf(stderr, "NitroSAT  %s  | %d vars  %d clauses\n", cnf_file, ns->num_vars, ns->num_clauses);
    else {
        /* Cinematic mode header */
        fprintf(stderr, "NitroSAT - Powerful linear-time MaxSAT approximator that consistently hits 99.5%%+ satisfaction across diverse categories including Graph Coloring, Clique, and Ramsey instances. Built by Sethu Iyer (https://sethuiyer.github.io)\n");
        fprintf(stderr, "NitroSAT  %s  | %d vars  %d clauses\n", cnf_file, ns->num_vars, ns->num_clauses);
    }
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
        } else if (strcmp(proof_format, "drat") == 0) {
            proof.backend = "solver_aware_drat";
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

        /* Save solution to .cnf.sol file */
        char sol_file[1024];
        snprintf(sol_file, sizeof(sol_file), "%s.sol", cnf_file);
        FILE *sol_fp = fopen(sol_file, "w");
        if (sol_fp) {
            for (int i = 1; i <= ns->num_vars; ++i) {
                fprintf(sol_fp, "%d ", (ns->x[i] > 0.5) ? i : -i);
            }
            fprintf(sol_fp, "0\n");
            fclose(sol_fp);
            fprintf(stderr, "Solution saved to %s\n", sol_file);
        } else {
            fprintf(stderr, "Warning: Could not save solution to %s\n", sol_file);
        }
    }

    nitrosat_free(ns);
    free(inst->cl_offs); free(inst->cl_flat); free(inst);
    return solved ? EXIT_SUCCESS : EXIT_FAILURE;
}
