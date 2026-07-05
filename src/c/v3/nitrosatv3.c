#define _POSIX_C_SOURCE 200809L

/*
 * NitroSAT V3 - bounded-memory, disk-streaming CNF solver.
 *
 * V2 intentionally remains unchanged.  V3 stores clauses as an appendable
 * binary stream and never constructs an in-memory clause table or a global
 * variable-to-clause index.  The resident set is O(variables + batch literals
 * + active-cache literals), independent of the total number of clauses.
 */

#include <ctype.h>
#include <errno.h>
#include <inttypes.h>
#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <time.h>
#include <unistd.h>

#define STORE_VERSION 2u
#define STORE_HEADER_SIZE 40u
#define STORE_FLAG_WEIGHTED 1u
#define CLAUSE_FLAG_HARD 1u
#define CLAUSE_FLAG_WEIGHTED 2u
#define DEFAULT_BATCH_CLAUSES 65536u
#define DEFAULT_ACTIVE_CLAUSES 65536u
#define DEFAULT_ACTIVE_LITERALS 1048576u
#define DEFAULT_EPOCHS 20u
#define DEFAULT_FINISHER_PASSES 256u
#define DEFAULT_FINISHER_MAX_CLAUSES UINT64_C(100000)
#define DEFAULT_FINISHER_BATCH_FLIPS 64u
#define DEFAULT_LR 0.015
#define MAX_CLAUSE_LITERALS UINT32_MAX

static const unsigned char STORE_MAGIC[8] = {'N','S','A','T','V','3','\r','\n'};

typedef struct {
    uint32_t num_vars;
    uint32_t flags;
    uint64_t num_clauses;
    uint64_t num_literals;
} StoreMeta;

typedef struct {
    FILE *fp;
    StoreMeta meta;
    uint64_t clauses_read;
    int32_t *lits;
    size_t lit_cap;
} ClauseReader;

typedef struct {
    uint32_t *offs;
    int32_t *lits;
    uint64_t *weights;
    uint8_t *hard;
    uint32_t clauses;
    uint32_t clause_cap;
    uint32_t literals;
    uint32_t literal_cap;
} ActiveCache;

typedef struct {
    double *x;
    double *m;
    double *v;
    double *grad;
    double *degrees;
    double mean_degree;
    uint32_t *marks;
    uint32_t *touched;
    uint32_t mark_generation;
    uint32_t touched_count;
    uint32_t num_vars;
    uint64_t rng;
    uint64_t batches;
} SolverState;

typedef struct {
    uint64_t satisfied;
    uint64_t hard_unsatisfied;
    uint64_t soft_unsatisfied;
    uint64_t soft_cost;
    uint64_t total_soft_weight;
    int cost_overflow;
} Verification;

typedef struct {
    const char *base_cnf;
    const char **add_files;
    int add_count;
    const char *store_path;
    const char *solution_path;
    int load_store;
    int keep_store;
    int store_only;
    uint32_t epochs;
    uint32_t finisher_passes;
    uint32_t finisher_batch_flips;
    uint64_t finisher_max_clauses;
    uint32_t batch_clauses;
    uint32_t active_clauses;
    uint32_t active_literals;
    double learning_rate;
    uint64_t seed;
    int json;
} Options;

static void put_u32(unsigned char *p, uint32_t v)
{
    p[0] = (unsigned char)v; p[1] = (unsigned char)(v >> 8);
    p[2] = (unsigned char)(v >> 16); p[3] = (unsigned char)(v >> 24);
}

static void put_u64(unsigned char *p, uint64_t v)
{
    for (unsigned i = 0; i < 8; ++i) p[i] = (unsigned char)(v >> (i * 8));
}

static uint32_t get_u32(const unsigned char *p)
{
    return (uint32_t)p[0] | ((uint32_t)p[1] << 8) |
           ((uint32_t)p[2] << 16) | ((uint32_t)p[3] << 24);
}

static uint64_t get_u64(const unsigned char *p)
{
    uint64_t v = 0;
    for (unsigned i = 0; i < 8; ++i) v |= (uint64_t)p[i] << (i * 8);
    return v;
}

static int write_header(FILE *fp, const StoreMeta *m)
{
    unsigned char h[STORE_HEADER_SIZE] = {0};
    memcpy(h, STORE_MAGIC, sizeof(STORE_MAGIC));
    put_u32(h + 8, STORE_VERSION);
    put_u32(h + 12, STORE_HEADER_SIZE);
    put_u32(h + 16, m->num_vars);
    put_u32(h + 20, m->flags);
    put_u64(h + 24, m->num_clauses);
    put_u64(h + 32, m->num_literals);
    if (fseeko(fp, 0, SEEK_SET) != 0) return 0;
    return fwrite(h, 1, sizeof(h), fp) == sizeof(h);
}

static int read_header(FILE *fp, StoreMeta *m)
{
    unsigned char h[STORE_HEADER_SIZE];
    if (fseeko(fp, 0, SEEK_SET) != 0 || fread(h, 1, sizeof(h), fp) != sizeof(h)) return 0;
    if (memcmp(h, STORE_MAGIC, sizeof(STORE_MAGIC)) != 0 ||
        get_u32(h + 8) != STORE_VERSION || get_u32(h + 12) != STORE_HEADER_SIZE) return 0;
    m->num_vars = get_u32(h + 16);
    m->flags = get_u32(h + 20);
    m->num_clauses = get_u64(h + 24);
    m->num_literals = get_u64(h + 32);
    return 1;
}

static int write_clause(FILE *fp, const int32_t *lits, uint32_t count,
                        uint64_t weight, uint32_t flags)
{
    unsigned char b[8];
    put_u32(b, count);
    put_u32(b + 4, flags);
    if (fwrite(b, 1, sizeof(b), fp) != sizeof(b)) return 0;
    put_u64(b, weight);
    if (fwrite(b, 1, sizeof(b), fp) != sizeof(b)) return 0;
    for (uint32_t i = 0; i < count; ++i) {
        put_u32(b, (uint32_t)lits[i]);
        if (fwrite(b, 1, 4, fp) != 4) return 0;
    }
    return 1;
}

static int grow_literals(int32_t **lits, size_t *cap, size_t need)
{
    if (need <= *cap) return 1;
    size_t n = *cap ? *cap : 16;
    while (n < need) {
        if (n > SIZE_MAX / 2) { errno = ENOMEM; return 0; }
        n *= 2;
    }
    if (n > SIZE_MAX / sizeof(**lits)) { errno = ENOMEM; return 0; }
    int32_t *p = realloc(*lits, n * sizeof(*p));
    if (!p) return 0;
    *lits = p; *cap = n;
    return 1;
}

/* Token scanner handles arbitrary line length, comments, and clauses spanning lines. */
static int next_dimacs_token(FILE *fp, int64_t *value, int *is_problem,
                             uint32_t *problem_vars, uint64_t *problem_clauses,
                             int *problem_weighted, uint64_t *problem_top,
                             int *problem_has_top, int *is_hard_marker)
{
    int ch;
    *is_problem = 0;
    *is_hard_marker = 0;
    for (;;) {
        ch = fgetc(fp);
        if (ch == EOF) return feof(fp) ? 0 : -1;
        if (isspace((unsigned char)ch)) continue;
        if (ch == 'c') {
            while ((ch = fgetc(fp)) != EOF && ch != '\n') {}
            continue;
        }
        if (ch == 'p') {
            char line[512], kind[16], extra;
            unsigned long long clauses, vars, top = 0;
            if (!fgets(line, sizeof(line), fp)) return -1;
            int fields = sscanf(line, "%15s %llu %llu %llu %c",
                                kind, &vars, &clauses, &top, &extra);
            int weighted = !strcmp(kind, "wcnf");
            if (vars > INT32_MAX ||
                (!weighted && (strcmp(kind, "cnf") != 0 || fields != 3)) ||
                (weighted && fields != 3 && fields != 4) ||
                (weighted && fields == 4 && top == 0))
                return -1;
            *problem_vars = (uint32_t)vars;
            *problem_clauses = (uint64_t)clauses;
            *problem_weighted = weighted;
            *problem_has_top = weighted && fields == 4;
            *problem_top = (uint64_t)top;
            *is_problem = 1;
            return 1;
        }
        ungetc(ch, fp);
        char token[64];
        size_t n = 0;
        while ((ch = fgetc(fp)) != EOF && !isspace((unsigned char)ch)) {
            if (n + 1 >= sizeof(token)) return -1;
            token[n++] = (char)ch;
        }
        token[n] = '\0';
        if (!strcmp(token, "h")) {
            *is_hard_marker = 1;
            *value = 0;
            return 1;
        }
        char *end = NULL;
        errno = 0;
        long long parsed = strtoll(token, &end, 10);
        if (errno || end == token || *end != '\0') return -1;
        *value = (int64_t)parsed;
        return 1;
    }
}

static int append_dimacs(FILE *store, const char *path, StoreMeta *total)
{
    FILE *in = fopen(path, "rb");
    if (!in) { fprintf(stderr, "cannot open %s: %s\n", path, strerror(errno)); return 0; }
    uint32_t declared_vars = 0;
    uint64_t declared_clauses = 0, parsed_clauses = 0;
    uint64_t top = 0, clause_weight_value = 1;
    int have_header = 0, weighted = 0, has_top = 0, have_clause_weight = 0;
    int clause_explicit_hard = 0;
    int32_t *clause = NULL;
    size_t clause_cap = 0, clause_len = 0;
    int ok = 1;

    for (;;) {
        int64_t tok = 0;
        int is_problem = 0;
        uint32_t pv = 0;
        uint64_t pc = 0;
        int pw = 0, pht = 0, hard_marker = 0;
        uint64_t pt = 0;
        int status = next_dimacs_token(in, &tok, &is_problem, &pv, &pc,
                                       &pw, &pt, &pht, &hard_marker);
        if (status == 0) break;
        if (status < 0) { fprintf(stderr, "%s: invalid DIMACS token\n", path); ok = 0; break; }
        if (is_problem) {
            if (have_header || parsed_clauses || clause_len || have_clause_weight) {
                fprintf(stderr, "%s: misplaced or duplicate problem line\n", path); ok = 0; break;
            }
            declared_vars = pv; declared_clauses = pc; weighted = pw;
            top = pt; has_top = pht; have_header = 1;
            if (weighted) total->flags |= STORE_FLAG_WEIGHTED;
            continue;
        }
        if (!have_header) { fprintf(stderr, "%s: clause before problem line\n", path); ok = 0; break; }
        if (weighted && !have_clause_weight) {
            if (hard_marker) {
                clause_weight_value = 1;
                clause_explicit_hard = 1;
                have_clause_weight = 1;
                continue;
            }
            if (tok <= 0) {
                fprintf(stderr, "%s: WCNF clause weight must be positive\n", path); ok = 0; break;
            }
            clause_weight_value = (uint64_t)tok;
            have_clause_weight = 1;
            continue;
        }
        if (hard_marker) {
            fprintf(stderr, "%s: hard marker must start a WCNF clause\n", path); ok = 0; break;
        }
        if (tok == 0) {
            uint32_t clause_flags = weighted ? CLAUSE_FLAG_WEIGHTED : CLAUSE_FLAG_HARD;
            if (weighted && (clause_explicit_hard ||
                             (has_top && clause_weight_value >= top)))
                clause_flags |= CLAUSE_FLAG_HARD;
            if (clause_len > MAX_CLAUSE_LITERALS ||
                !write_clause(store, clause, (uint32_t)clause_len,
                              clause_weight_value, clause_flags)) {
                fprintf(stderr, "%s: failed writing clause: %s\n", path, strerror(errno)); ok = 0; break;
            }
            parsed_clauses++;
            if (UINT64_MAX - total->num_literals < clause_len) { ok = 0; errno = EOVERFLOW; break; }
            total->num_literals += (uint64_t)clause_len;
            clause_len = 0;
            clause_weight_value = 1;
            have_clause_weight = 0;
            clause_explicit_hard = 0;
        } else {
            if (tok < INT32_MIN || tok > INT32_MAX || tok == INT32_MIN ||
                (uint64_t)llabs((long long)tok) > declared_vars) {
                fprintf(stderr, "%s: literal %" PRId64 " exceeds declared variable range\n", path, tok);
                ok = 0; break;
            }
            if (!grow_literals(&clause, &clause_cap, clause_len + 1)) { ok = 0; break; }
            clause[clause_len++] = (int32_t)tok;
        }
    }
    if (ok && (clause_len || have_clause_weight)) {
        fprintf(stderr, "%s: unterminated final clause\n", path); ok = 0;
    }
    if (ok && parsed_clauses != declared_clauses) {
        fprintf(stderr, "%s: clause count mismatch (header=%" PRIu64 ", parsed=%" PRIu64 ")\n",
                path, declared_clauses, parsed_clauses); ok = 0;
    }
    if (ok) {
        if (UINT64_MAX - total->num_clauses < parsed_clauses) { ok = 0; errno = EOVERFLOW; }
        else total->num_clauses += parsed_clauses;
        if (declared_vars > total->num_vars) total->num_vars = declared_vars;
    }
    if (!ok && errno) fprintf(stderr, "%s: %s\n", path, strerror(errno));
    free(clause);
    fclose(in);
    return ok;
}

static int create_store(const Options *opt, StoreMeta *meta)
{
    FILE *fp = fopen(opt->store_path, "w+b");
    if (!fp) { fprintf(stderr, "cannot create %s: %s\n", opt->store_path, strerror(errno)); return 0; }
    memset(meta, 0, sizeof(*meta));
    int ok = write_header(fp, meta) && fseeko(fp, STORE_HEADER_SIZE, SEEK_SET) == 0;
    if (ok) ok = append_dimacs(fp, opt->base_cnf, meta);
    for (int i = 0; ok && i < opt->add_count; ++i) ok = append_dimacs(fp, opt->add_files[i], meta);
    if (ok) ok = write_header(fp, meta) && fflush(fp) == 0;
    if (fclose(fp) != 0) ok = 0;
    if (!ok) unlink(opt->store_path);
    return ok;
}

static int load_or_extend_store(const Options *opt, StoreMeta *meta)
{
    FILE *fp = fopen(opt->store_path, opt->add_count ? "r+b" : "rb");
    if (!fp) { fprintf(stderr, "cannot open store %s: %s\n", opt->store_path, strerror(errno)); return 0; }
    StoreMeta original = {0};
    int ok = read_header(fp, &original);
    off_t original_size = 0;
    if (ok && fseeko(fp, 0, SEEK_END) == 0) original_size = ftello(fp);
    else ok = 0;
    *meta = original;
    for (int i = 0; ok && i < opt->add_count; ++i) ok = append_dimacs(fp, opt->add_files[i], meta);
    if (ok && opt->add_count) ok = write_header(fp, meta) && fflush(fp) == 0;
    if (!ok && opt->add_count) {
        /* Roll back both bytes and metadata if any increment is invalid. */
        clearerr(fp);
        /* Flush pending append bytes before truncating, otherwise fclose could
           replay the stdio buffer beyond the restored end of file. */
        int rollback_ok = fflush(fp) == 0 && ftruncate(fileno(fp), original_size) == 0;
        if (!rollback_ok || !write_header(fp, &original) || fflush(fp) != 0)
            fprintf(stderr, "warning: failed to roll back store %s\n", opt->store_path);
        *meta = original;
    }
    if (fclose(fp) != 0) ok = 0;
    return ok;
}

static int reader_open(ClauseReader *r, const char *path)
{
    memset(r, 0, sizeof(*r));
    r->fp = fopen(path, "rb");
    if (!r->fp) return 0;
    if (!read_header(r->fp, &r->meta) || fseeko(r->fp, STORE_HEADER_SIZE, SEEK_SET) != 0) {
        fclose(r->fp); r->fp = NULL; return 0;
    }
    return 1;
}

static int reader_rewind(ClauseReader *r)
{
    r->clauses_read = 0;
    clearerr(r->fp);
    return fseeko(r->fp, STORE_HEADER_SIZE, SEEK_SET) == 0;
}

static int reader_next(ClauseReader *r, const int32_t **lits, uint32_t *count,
                       uint64_t *weight, uint32_t *flags)
{
    if (r->clauses_read >= r->meta.num_clauses) return 0;
    unsigned char b[8];
    if (fread(b, 1, 8, r->fp) != 8) return -1;
    uint32_t n = get_u32(b);
    uint32_t clause_flags = get_u32(b + 4);
    if (clause_flags & ~(CLAUSE_FLAG_HARD | CLAUSE_FLAG_WEIGHTED)) return -1;
    if (fread(b, 1, 8, r->fp) != 8) return -1;
    uint64_t clause_weight_value = get_u64(b);
    if (!clause_weight_value) return -1;
    if (!grow_literals(&r->lits, &r->lit_cap, n)) return -1;
    for (uint32_t i = 0; i < n; ++i) {
        if (fread(b, 1, 4, r->fp) != 4) return -1;
        r->lits[i] = (int32_t)get_u32(b);
        uint32_t v = (uint32_t)(r->lits[i] < 0 ? -(int64_t)r->lits[i] : r->lits[i]);
        if (!v || v > r->meta.num_vars) return -1;
    }
    r->clauses_read++;
    *lits = r->lits; *count = n;
    *weight = clause_weight_value; *flags = clause_flags;
    return 1;
}

static void reader_close(ClauseReader *r)
{
    if (r->fp) fclose(r->fp);
    free(r->lits);
    memset(r, 0, sizeof(*r));
}

static int validate_store_file(const char *path, const StoreMeta *expected)
{
    ClauseReader r = {0};
    if (!reader_open(&r, path)) return 0;
    uint64_t literals = 0;
    const int32_t *lits;
    uint32_t count;
    uint64_t weight;
    uint32_t flags;
    int status;
    while ((status = reader_next(&r, &lits, &count, &weight, &flags)) > 0) {
        (void)lits;
        (void)weight;
        (void)flags;
        if (UINT64_MAX - literals < count) { status = -1; break; }
        literals += count;
    }
    int trailing = fgetc(r.fp);
    int ok = status == 0 && r.clauses_read == expected->num_clauses &&
             literals == expected->num_literals && trailing == EOF && feof(r.fp);
    reader_close(&r);
    return ok;
}

static uint64_t rng_next(SolverState *s)
{
    uint64_t x = s->rng;
    x ^= x >> 12; x ^= x << 25; x ^= x >> 27;
    s->rng = x;
    return x * UINT64_C(2685821657736338717);
}

static double rng_unit(SolverState *s)
{
    return (double)(rng_next(s) >> 11) * (1.0 / 9007199254740992.0);
}

static int solver_init(SolverState *s, uint32_t vars, uint64_t seed)
{
    memset(s, 0, sizeof(*s));
    s->num_vars = vars;
    s->rng = seed ? seed : UINT64_C(0x9e3779b97f4a7c15);
    size_t n = (size_t)vars + 1;
    if (n > SIZE_MAX / sizeof(double)) return 0;
    s->x = malloc(n * sizeof(double));
    s->m = calloc(n, sizeof(double));
    s->v = calloc(n, sizeof(double));
    s->grad = calloc(n, sizeof(double));
    s->degrees = calloc(n, sizeof(double));
    s->marks = calloc(n, sizeof(uint32_t));
    s->touched = malloc(n * sizeof(uint32_t));
    if (!s->x || !s->m || !s->v || !s->grad || !s->degrees ||
        !s->marks || !s->touched) return 0;
    for (uint32_t i = 1; i <= vars; ++i) s->x[i] = 0.45 + 0.1 * rng_unit(s);
    s->mark_generation = 1;
    return 1;
}

static void solver_free(SolverState *s)
{
    free(s->x); free(s->m); free(s->v); free(s->grad); free(s->degrees);
    free(s->marks); free(s->touched);
    memset(s, 0, sizeof(*s));
}

static int compute_stream_degrees(ClauseReader *reader, SolverState *solver)
{
    if (!reader_rewind(reader)) return 0;
    const int32_t *lits;
    uint32_t n, flags;
    uint64_t weight;
    int status;
    while ((status = reader_next(reader, &lits, &n, &weight, &flags)) > 0) {
        (void)weight; (void)flags;
        double contribution = n ? (double)(n - 1) : 0.0;
        for (uint32_t i = 0; i < n; ++i) {
            uint32_t v = (uint32_t)(lits[i] < 0 ? -(int64_t)lits[i] : lits[i]);
            solver->degrees[v] += contribution;
        }
    }
    if (status < 0 || reader->clauses_read != reader->meta.num_clauses) return 0;
    double sum = 0.0;
    for (uint32_t v = 1; v <= solver->num_vars; ++v) sum += solver->degrees[v];
    solver->mean_degree = solver->num_vars ? sum / solver->num_vars : 0.0;
    return 1;
}

static int cache_init(ActiveCache *c, uint32_t clauses, uint32_t literals)
{
    memset(c, 0, sizeof(*c));
    c->clause_cap = clauses; c->literal_cap = literals;
    c->offs = malloc(((size_t)clauses + 1) * sizeof(uint32_t));
    c->lits = malloc((size_t)(literals ? literals : 1) * sizeof(int32_t));
    c->weights = malloc((size_t)clauses * sizeof(uint64_t));
    c->hard = malloc((size_t)clauses * sizeof(uint8_t));
    if (!c->offs || !c->lits || !c->weights || !c->hard) return 0;
    c->offs[0] = 0;
    return 1;
}

static void cache_reset(ActiveCache *c)
{
    c->clauses = 0; c->literals = 0; c->offs[0] = 0;
}

static void cache_add(ActiveCache *c, const int32_t *lits, uint32_t n,
                      uint64_t weight, int hard)
{
    if (c->clauses >= c->clause_cap || n > c->literal_cap - c->literals) return;
    memcpy(c->lits + c->literals, lits, (size_t)n * sizeof(*lits));
    c->weights[c->clauses] = weight;
    c->hard[c->clauses] = (uint8_t)hard;
    c->literals += n;
    c->offs[++c->clauses] = c->literals;
}

static void cache_free(ActiveCache *c)
{
    free(c->offs); free(c->lits); free(c->weights); free(c->hard);
    memset(c, 0, sizeof(*c));
}

static int literal_is_satisfied(int32_t lit, double x)
{
    return lit > 0 ? x > 0.5 : x <= 0.5;
}

static double clause_weight(uint64_t clause_id)
{
    /* Procedural approximation to V2's nth-prime weighting; no 1B-entry sieve. */
    double n = (double)clause_id + 2.0;
    double p = n < 6.0 ? (2.0 + n) : n * (log(n) + log(log(n)) - 1.0);
    return 50.0 / (1.0 + log(fmax(p, 2.0)));
}

static void touch_gradient(SolverState *s, uint32_t var, double contribution)
{
    if (s->marks[var] != s->mark_generation) {
        s->marks[var] = s->mark_generation;
        s->grad[var] = 0.0;
        s->touched[s->touched_count++] = var;
    }
    s->grad[var] += contribution;
}

static void apply_batch(SolverState *s, double lr)
{
    if (!s->touched_count) return;
    s->batches++;
    double b1 = 0.9, b2 = 0.999;
    double b1c = 1.0 - pow(b1, (double)s->batches);
    double b2c = 1.0 - pow(b2, (double)s->batches);
    for (uint32_t i = 0; i < s->touched_count; ++i) {
        uint32_t var = s->touched[i];
        double heat = exp(-0.5 * s->degrees[var] / (s->mean_degree + 1.0));
        double g = s->grad[var] * (1.0 + 0.1 * heat);
        if (g > 1000.0) g = 1000.0;
        if (g < -1000.0) g = -1000.0;
        s->m[var] = b1 * s->m[var] + (1.0 - b1) * g;
        s->v[var] = b2 * s->v[var] + (1.0 - b2) * g * g;
        s->x[var] -= lr * (s->m[var] / b1c) / (sqrt(s->v[var] / b2c) + 1e-8);
        if (s->x[var] < 1e-6) s->x[var] = 1e-6;
        if (s->x[var] > 1.0 - 1e-6) s->x[var] = 1.0 - 1e-6;
    }
    s->touched_count = 0;
    if (++s->mark_generation == 0) {
        memset(s->marks, 0, ((size_t)s->num_vars + 1) * sizeof(uint32_t));
        s->mark_generation = 1;
    }
}

static int process_epoch(ClauseReader *r, SolverState *s, ActiveCache *cache,
                         uint32_t batch_clauses, double lr, uint64_t *unsat_out)
{
    if (!reader_rewind(r)) return 0;
    cache_reset(cache);
    uint64_t unsat = 0, clause_id = 0;
    uint32_t in_batch = 0;
    const int32_t *lits;
    uint32_t n;
    uint64_t stored_weight;
    uint32_t clause_flags;
    int status;
    while ((status = reader_next(r, &lits, &n, &stored_weight, &clause_flags)) > 0) {
        int satisfied = 0;
        double violation = 1.0;
        for (uint32_t i = 0; i < n; ++i) {
            uint32_t var = (uint32_t)(lits[i] < 0 ? -(int64_t)lits[i] : lits[i]);
            double x = s->x[var];
            if (literal_is_satisfied(lits[i], x)) satisfied = 1;
            double lv = lits[i] > 0 ? 1.0 - x : x;
            violation *= fmax(lv, 1e-12);
        }
        int hard = (clause_flags & CLAUSE_FLAG_HARD) != 0;
        if (!satisfied) { unsat++; cache_add(cache, lits, n, stored_weight, hard); }
        double base_weight = (clause_flags & CLAUSE_FLAG_WEIGHTED)
                           ? (hard ? 1e6 : fmin((double)stored_weight, 1e5))
                           : clause_weight(clause_id);
        double coef = base_weight * violation;
        for (uint32_t i = 0; i < n; ++i) {
            uint32_t var = (uint32_t)(lits[i] < 0 ? -(int64_t)lits[i] : lits[i]);
            double lv = lits[i] > 0 ? 1.0 - s->x[var] : s->x[var];
            touch_gradient(s, var, coef * (lits[i] > 0 ? -1.0 : 1.0) / fmax(lv, 1e-12));
        }
        clause_id++;
        if (++in_batch == batch_clauses) { apply_batch(s, lr); in_batch = 0; }
    }
    apply_batch(s, lr);
    if (status < 0 || r->clauses_read != r->meta.num_clauses) return 0;
    *unsat_out = unsat;
    return 1;
}

static void repair_active_cache(ActiveCache *c, SolverState *s)
{
    /* Bounded stochastic local repair. Recheck because earlier flips overlap. */
    if (!c->clauses) return;
    uint64_t rounds = (uint64_t)c->clauses * 2;
    for (uint64_t r = 0; r < rounds; ++r) {
        uint32_t ci = (uint32_t)(rng_next(s) % c->clauses);
        uint32_t alternative = (uint32_t)(rng_next(s) % c->clauses);
        if ((c->hard[alternative] && !c->hard[ci]) ||
            (c->hard[alternative] == c->hard[ci] && c->weights[alternative] > c->weights[ci]))
            ci = alternative;
        uint32_t begin = c->offs[ci], end = c->offs[ci + 1];
        int sat = 0;
        for (uint32_t j = begin; j < end; ++j) {
            int32_t lit = c->lits[j];
            uint32_t v = (uint32_t)(lit < 0 ? -(int64_t)lit : lit);
            if (literal_is_satisfied(lit, s->x[v])) { sat = 1; break; }
        }
        if (!sat && end > begin) {
            int32_t lit = c->lits[begin + (uint32_t)(rng_next(s) % (end - begin))];
            uint32_t v = (uint32_t)(lit < 0 ? -(int64_t)lit : lit);
            s->x[v] = lit > 0 ? 0.999 : 0.001;
        }
    }
}

static void add_cost(uint64_t *total, uint64_t value, int *overflow)
{
    if (UINT64_MAX - *total < value) { *total = UINT64_MAX; *overflow = 1; }
    else *total += value;
}

static int verify_assignment(ClauseReader *r, SolverState *s, Verification *result)
{
    if (!reader_rewind(r)) return 0;
    memset(result, 0, sizeof(*result));
    const int32_t *lits;
    uint32_t n;
    uint64_t weight;
    uint32_t flags;
    int status;
    while ((status = reader_next(r, &lits, &n, &weight, &flags)) > 0) {
        int clause_sat = 0;
        for (uint32_t i = 0; i < n; ++i) {
            uint32_t v = (uint32_t)(lits[i] < 0 ? -(int64_t)lits[i] : lits[i]);
            if (literal_is_satisfied(lits[i], s->x[v])) { clause_sat = 1; break; }
        }
        result->satisfied += (uint64_t)clause_sat;
        if (!(flags & CLAUSE_FLAG_HARD)) {
            add_cost(&result->total_soft_weight, weight, &result->cost_overflow);
            if (!clause_sat) {
                result->soft_unsatisfied++;
                add_cost(&result->soft_cost, weight, &result->cost_overflow);
            }
        } else if (!clause_sat) result->hard_unsatisfied++;
    }
    if (status < 0 || r->clauses_read != r->meta.num_clauses) return 0;
    return 1;
}

static int verification_better(const Verification *candidate, const Verification *best)
{
    if (candidate->hard_unsatisfied != best->hard_unsatisfied)
        return candidate->hard_unsatisfied < best->hard_unsatisfied;
    if (candidate->soft_cost != best->soft_cost)
        return candidate->soft_cost < best->soft_cost;
    return candidate->satisfied > best->satisfied;
}

static void retain_best(const Verification *candidate, Verification *best,
                        int *have_best, double *best_x, const SolverState *solver)
{
    if (!*have_best || verification_better(candidate, best)) {
        *best = *candidate;
        memcpy(best_x, solver->x, ((size_t)solver->num_vars + 1) * sizeof(double));
        *have_best = 1;
    }
}

static double optimization_weight(uint64_t clause_id, uint64_t stored_weight,
                                  uint32_t clause_flags)
{
    if (!(clause_flags & CLAUSE_FLAG_WEIGHTED)) return clause_weight(clause_id);
    return (clause_flags & CLAUSE_FLAG_HARD)
         ? 1e6 : fmin((double)stored_weight, 1e5);
}

static int streaming_make_break_finisher(ClauseReader *reader, SolverState *solver,
                                         uint32_t max_passes, uint32_t batch_flips,
                                         Verification *best,
                                         int *have_best, double *best_x,
                                         uint32_t *passes_run)
{
    size_t vars_size = (size_t)solver->num_vars + 1;
    double *make = calloc(vars_size, sizeof(double));
    double *brk = calloc(vars_size, sizeof(double));
    double *heat_force = calloc(vars_size, sizeof(double));
    uint32_t *batch_marks = calloc(vars_size, sizeof(uint32_t));
    uint32_t *candidates = malloc((size_t)batch_flips * sizeof(uint32_t));
    double *candidate_scores = malloc((size_t)batch_flips * sizeof(double));
    int32_t *escape_clause = NULL;
    size_t escape_cap = 0;
    uint32_t previous_var = 0;
    if (!make || !brk || !heat_force || !batch_marks || !candidates || !candidate_scores) {
        free(make); free(brk); free(heat_force); free(batch_marks);
        free(candidates); free(candidate_scores);
        return 0;
    }

    for (uint32_t pass = 0; pass < max_passes; ++pass) {
        memset(make, 0, vars_size * sizeof(double));
        memset(brk, 0, vars_size * sizeof(double));
        memset(heat_force, 0, vars_size * sizeof(double));
        if (!reader_rewind(reader)) {
            free(make); free(brk); free(heat_force); free(batch_marks); free(candidates);
            free(candidate_scores); free(escape_clause); return 0;
        }
        uint64_t clause_id = 0;
        uint64_t escape_seen = 0;
        uint32_t escape_n = 0;
        int escape_is_hard = 0;
        const int32_t *lits;
        uint32_t n, flags;
        uint64_t stored_weight;
        int status;
        while ((status = reader_next(reader, &lits, &n, &stored_weight, &flags)) > 0) {
            uint32_t sat_count = 0, sole_var = 0;
            for (uint32_t i = 0; i < n; ++i) {
                uint32_t var = (uint32_t)(lits[i] < 0 ? -(int64_t)lits[i] : lits[i]);
                if (literal_is_satisfied(lits[i], solver->x[var])) {
                    sat_count++;
                    sole_var = var;
                }
            }
            double weight = optimization_weight(clause_id, stored_weight, flags);
            if (sat_count == 0) {
                for (uint32_t i = 0; i < n; ++i) {
                    uint32_t var = (uint32_t)(lits[i] < 0 ? -(int64_t)lits[i] : lits[i]);
                    make[var] += weight;
                    double hk = exp(-0.5 * solver->degrees[var] / (solver->mean_degree + 1.0));
                    heat_force[var] += weight * hk;
                }
                int hard = (flags & CLAUSE_FLAG_HARD) != 0;
                if (n && hard && !escape_is_hard) {
                    escape_seen = 0;
                    escape_n = 0;
                    escape_is_hard = 1;
                }
                if (n && hard == escape_is_hard) {
                    escape_seen++;
                    int replace = !escape_n || rng_next(solver) % escape_seen == 0;
                    if (!replace) { clause_id++; continue; }
                    if (!grow_literals(&escape_clause, &escape_cap, n)) {
                        free(make); free(brk); free(heat_force); free(batch_marks); free(candidates);
                        free(candidate_scores); free(escape_clause); return 0;
                    }
                    memcpy(escape_clause, lits, (size_t)n * sizeof(*lits));
                    escape_n = n;
                    escape_is_hard = hard;
                }
            } else if (sat_count == 1) brk[sole_var] += weight;
            clause_id++;
        }
        if (status < 0 || reader->clauses_read != reader->meta.num_clauses) {
            free(make); free(brk); free(heat_force); free(batch_marks); free(candidates);
            free(candidate_scores); free(escape_clause); return 0;
        }
        if (!escape_n) break;

        uint32_t selected = 0, candidate_count = 0;
        double best_gain = -INFINITY;
        for (uint32_t v = 1; v <= solver->num_vars; ++v) {
            if (v == previous_var) continue;
            double gain = make[v] - brk[v];
            if (gain > best_gain) { best_gain = gain; selected = v; }
            if (gain > 1e-12) {
                double score = gain + 0.1 * heat_force[v];
                if (candidate_count < batch_flips) {
                    candidates[candidate_count] = v;
                    candidate_scores[candidate_count++] = score;
                } else {
                    uint32_t min_i = 0;
                    for (uint32_t i = 1; i < candidate_count; ++i)
                        if (candidate_scores[i] < candidate_scores[min_i]) min_i = i;
                    if (score > candidate_scores[min_i]) {
                        candidates[min_i] = v; candidate_scores[min_i] = score;
                    }
                }
            }
        }
        /* At a local minimum, use the least-damaging variable from a currently
           unsatisfied clause. This is a bounded WalkSAT escape informed by the
           global streamed break scores. */
        if (!(best_gain > 1e-12)) {
            selected = 0;
            best_gain = -INFINITY;
            if (rng_unit(solver) < 0.15) {
                for (uint32_t tries = 0; tries < escape_n * 2; ++tries) {
                    int32_t lit = escape_clause[rng_next(solver) % escape_n];
                    uint32_t v = (uint32_t)(lit < 0 ? -(int64_t)lit : lit);
                    if (v != previous_var || escape_n == 1) { selected = v; break; }
                }
            } else {
                for (uint32_t i = 0; i < escape_n; ++i) {
                    uint32_t v = (uint32_t)(escape_clause[i] < 0
                                          ? -(int64_t)escape_clause[i] : escape_clause[i]);
                    if (v == previous_var && escape_n > 1) continue;
                    double gain = make[v] - brk[v];
                    if (gain > best_gain ||
                        (gain == best_gain && (rng_next(solver) & 1))) {
                        best_gain = gain; selected = v;
                    }
                }
            }
        }
        if (!selected) break;
        if (best_gain > 1e-12 && candidate_count) {
            memset(batch_marks, 0, vars_size * sizeof(uint32_t));
            for (uint32_t i = 0; i < candidate_count; ++i) batch_marks[candidates[i]] = i + 1;
            /* Conflict screen: a clause may contain at most one simultaneous
               flip. Keep its highest heat/gain candidate and suppress peers. */
            if (candidate_count > 1) {
                if (!reader_rewind(reader)) {
                    free(make); free(brk); free(heat_force); free(batch_marks);
                    free(candidates); free(candidate_scores); free(escape_clause); return 0;
                }
                while ((status = reader_next(reader, &lits, &n, &stored_weight, &flags)) > 0) {
                    uint32_t keep = 0;
                    for (uint32_t j = 0; j < n; ++j) {
                        uint32_t v = (uint32_t)(lits[j] < 0 ? -(int64_t)lits[j] : lits[j]);
                        uint32_t mark = batch_marks[v];
                        if (mark && (!keep || candidate_scores[mark-1] > candidate_scores[keep-1]))
                            keep = mark;
                    }
                    if (keep) for (uint32_t j = 0; j < n; ++j) {
                        uint32_t v = (uint32_t)(lits[j] < 0 ? -(int64_t)lits[j] : lits[j]);
                        if (batch_marks[v] && batch_marks[v] != keep) batch_marks[v] = 0;
                    }
                }
                if (status < 0) {
                    free(make); free(brk); free(heat_force); free(batch_marks);
                    free(candidates); free(candidate_scores); free(escape_clause); return 0;
                }
            }
            for (uint32_t i = 0; i < candidate_count; ++i) {
                uint32_t v = candidates[i];
                if (batch_marks[v]) solver->x[v] = solver->x[v] > 0.5 ? 0.001 : 0.999;
            }
            previous_var = 0;
        } else {
            previous_var = selected;
            solver->x[selected] = solver->x[selected] > 0.5 ? 0.001 : 0.999;
        }
        Verification current;
        if (!verify_assignment(reader, solver, &current)) {
            free(make); free(brk); free(heat_force); free(batch_marks); free(candidates);
            free(candidate_scores); free(escape_clause); return 0;
        }
        retain_best(&current, best, have_best, best_x, solver);
        *passes_run = pass + 1;
        if (current.satisfied == reader->meta.num_clauses) break;
    }
    free(make); free(brk); free(heat_force); free(batch_marks); free(candidates);
    free(candidate_scores); free(escape_clause);
    return 1;
}

static double elapsed_ms(const struct timespec *a, const struct timespec *b)
{
    return (double)(b->tv_sec - a->tv_sec) * 1000.0 +
           (double)(b->tv_nsec - a->tv_nsec) / 1e6;
}

static int parse_u32(const char *s, uint32_t *out)
{
    char *end = NULL; errno = 0;
    unsigned long long v = strtoull(s, &end, 10);
    if (errno || end == s || *end || v > UINT32_MAX) return 0;
    *out = (uint32_t)v; return 1;
}

static int parse_u64(const char *s, uint64_t *out)
{
    char *end = NULL; errno = 0;
    unsigned long long v = strtoull(s, &end, 10);
    if (errno || end == s || *end) return 0;
    *out = (uint64_t)v; return 1;
}

static void usage(const char *argv0)
{
    fprintf(stderr,
        "Usage: %s BASE.cnf [options]\n"
        "       %s --load-store FILE [options]\n"
        "  --add FILE              append an incremental DIMACS constraint set (repeatable)\n"
        "  --store FILE            binary clause-store path (default: temporary)\n"
        "  --load-store FILE       solve/extend an existing V3 store\n"
        "  --solution FILE         write the final Boolean assignment\n"
        "  --keep-store            retain the binary store after solving\n"
        "  --store-only            build/validate the store without solving\n"
        "  --epochs N              streaming passes (default: %u)\n"
        "  --finisher-passes N     global streamed make/break flips (default: %u)\n"
        "  --finisher-batch-flips N  positive-gain flips per pass (default: %u)\n"
        "  --finisher-max-clauses N  clause ceiling for finisher (default: %" PRIu64 ")\n"
        "  --batch-clauses N       optimizer batch size (default: %u)\n"
        "  --active-clauses N      local-repair clause bound (default: %u)\n"
        "  --active-literals N     local-repair literal bound (default: %u)\n"
        "  --learning-rate X       optimizer learning rate (default: %.3f)\n"
        "  --seed N                deterministic random seed\n"
        "  --cinematic             text output instead of JSON\n",
        argv0, argv0, DEFAULT_EPOCHS, DEFAULT_FINISHER_PASSES,
        DEFAULT_FINISHER_BATCH_FLIPS,
        DEFAULT_FINISHER_MAX_CLAUSES, DEFAULT_BATCH_CLAUSES,
        DEFAULT_ACTIVE_CLAUSES, DEFAULT_ACTIVE_LITERALS, DEFAULT_LR);
}

static int parse_options(int argc, char **argv, Options *o)
{
    memset(o, 0, sizeof(*o));
    o->epochs = DEFAULT_EPOCHS; o->finisher_passes = DEFAULT_FINISHER_PASSES;
    o->finisher_batch_flips = DEFAULT_FINISHER_BATCH_FLIPS;
    o->finisher_max_clauses = DEFAULT_FINISHER_MAX_CLAUSES;
    o->batch_clauses = DEFAULT_BATCH_CLAUSES;
    o->active_clauses = DEFAULT_ACTIVE_CLAUSES; o->active_literals = DEFAULT_ACTIVE_LITERALS;
    o->learning_rate = DEFAULT_LR; o->seed = 42; o->json = 1;
    if (argc < 2) return 0;
    o->add_files = calloc((size_t)argc, sizeof(char *));
    if (!o->add_files) return 0;
    int first_option;
    if (!strcmp(argv[1], "--load-store")) {
        if (argc < 3) return 0;
        o->load_store = 1; o->store_path = argv[2]; o->keep_store = 1; first_option = 3;
    } else if (argv[1][0] != '-') {
        o->base_cnf = argv[1]; first_option = 2;
    } else return 0;
    for (int i = first_option; i < argc; ++i) {
        if (!strcmp(argv[i], "--add") && i + 1 < argc) o->add_files[o->add_count++] = argv[++i];
        else if (!strcmp(argv[i], "--store") && i + 1 < argc) o->store_path = argv[++i];
        else if (!strcmp(argv[i], "--solution") && i + 1 < argc) o->solution_path = argv[++i];
        else if (!strcmp(argv[i], "--keep-store")) o->keep_store = 1;
        else if (!strcmp(argv[i], "--store-only")) o->store_only = 1;
        else if (!strcmp(argv[i], "--cinematic")) o->json = 0;
        else if (!strcmp(argv[i], "--epochs") && i + 1 < argc && parse_u32(argv[i+1], &o->epochs)) i++;
        else if (!strcmp(argv[i], "--finisher-passes") && i + 1 < argc && parse_u32(argv[i+1], &o->finisher_passes)) i++;
        else if (!strcmp(argv[i], "--finisher-batch-flips") && i + 1 < argc && parse_u32(argv[i+1], &o->finisher_batch_flips)) i++;
        else if (!strcmp(argv[i], "--finisher-max-clauses") && i + 1 < argc && parse_u64(argv[i+1], &o->finisher_max_clauses)) i++;
        else if (!strcmp(argv[i], "--batch-clauses") && i + 1 < argc && parse_u32(argv[i+1], &o->batch_clauses)) i++;
        else if (!strcmp(argv[i], "--active-clauses") && i + 1 < argc && parse_u32(argv[i+1], &o->active_clauses)) i++;
        else if (!strcmp(argv[i], "--active-literals") && i + 1 < argc && parse_u32(argv[i+1], &o->active_literals)) i++;
        else if (!strcmp(argv[i], "--seed") && i + 1 < argc && parse_u64(argv[i+1], &o->seed)) i++;
        else if (!strcmp(argv[i], "--learning-rate") && i + 1 < argc) {
            char *end = NULL; errno = 0; o->learning_rate = strtod(argv[++i], &end);
            if (errno || end == argv[i] || *end || !(o->learning_rate > 0.0)) return 0;
        } else return 0;
    }
    return o->batch_clauses > 0 && o->active_clauses > 0 &&
           o->active_literals > 0 && o->finisher_batch_flips > 0;
}

static int write_solution(const char *path, const SolverState *solver)
{
    if (!path) return 1;
    FILE *fp = fopen(path, "wb");
    if (!fp) return 0;
    for (uint32_t i = 1; i <= solver->num_vars; ++i) {
        int64_t lit = solver->x[i] > 0.5 ? (int64_t)i : -(int64_t)i;
        if (fprintf(fp, "%" PRId64 " ", lit) < 0) { fclose(fp); return 0; }
    }
    int ok = fprintf(fp, "0\n") >= 0;
    if (fclose(fp) != 0) ok = 0;
    return ok;
}

int main(int argc, char **argv)
{
    Options opt;
    if (!parse_options(argc, argv, &opt)) { usage(argv[0]); free(opt.add_files); return 2; }
    char temp_path[] = "/tmp/nitrosatv3-XXXXXX";
    if (!opt.store_path) {
        int fd = mkstemp(temp_path);
        if (fd < 0) { perror("mkstemp"); free(opt.add_files); return 2; }
        close(fd); opt.store_path = temp_path;
    } else if (!opt.load_store) {
        /* An explicitly named store is a useful artifact unless removal is requested by omission. */
        opt.keep_store = 1;
    }

    struct timespec start, built, finished;
    clock_gettime(CLOCK_MONOTONIC, &start);
    StoreMeta meta;
    int store_ok = opt.load_store ? load_or_extend_store(&opt, &meta) : create_store(&opt, &meta);
    if (!store_ok) { free(opt.add_files); return 2; }
    clock_gettime(CLOCK_MONOTONIC, &built);

    if (opt.store_only) {
        if (opt.load_store && !validate_store_file(opt.store_path, &meta)) {
            fprintf(stderr, "store validation failed: %s\n", opt.store_path);
            free(opt.add_files);
            return 2;
        }
        printf("{\"store\":\"%s\",\"format\":\"%s\",\"variables\":%u,"
               "\"clauses\":%" PRIu64 ",\"literals\":%" PRIu64
               ",\"increment_files\":%d}\n",
               opt.store_path, (meta.flags & STORE_FLAG_WEIGHTED) ? "wcnf" : "cnf",
               meta.num_vars, meta.num_clauses, meta.num_literals, opt.add_count);
        free(opt.add_files);
        return 0;
    }

    ClauseReader reader = {0};
    SolverState solver = {0};
    ActiveCache cache = {0};
    if (!reader_open(&reader, opt.store_path) || !solver_init(&solver, meta.num_vars, opt.seed) ||
        !compute_stream_degrees(&reader, &solver) ||
        !cache_init(&cache, opt.active_clauses, opt.active_literals)) {
        fprintf(stderr, "V3 initialization failed: %s\n", strerror(errno));
        reader_close(&reader); solver_free(&solver); cache_free(&cache);
        if (!opt.keep_store) unlink(opt.store_path);
        free(opt.add_files); return 2;
    }

    uint64_t unsat = meta.num_clauses;
    uint32_t epochs_run = 0;
    int io_ok = 1;
    int weighted_store = (meta.flags & STORE_FLAG_WEIGHTED) != 0;
    Verification best_verification = {0};
    int have_best = 0;
    double *best_x = malloc(((size_t)meta.num_vars + 1) * sizeof(double));
    if (!best_x) io_ok = 0;
    for (uint32_t e = 0; io_ok && e < opt.epochs; ++e) {
        if (!process_epoch(&reader, &solver, &cache, opt.batch_clauses,
                           opt.learning_rate / sqrt((double)e + 1.0), &unsat)) {
            io_ok = 0; break;
        }
        epochs_run = e + 1;
        if (!opt.json) fprintf(stderr, "epoch %u: unsatisfied=%" PRIu64 "/%" PRIu64 " cached=%u\n",
                               epochs_run, unsat, meta.num_clauses, cache.clauses);
        if (weighted_store) {
            Verification before_repair;
            if (!verify_assignment(&reader, &solver, &before_repair)) { io_ok = 0; break; }
            retain_best(&before_repair, &best_verification, &have_best, best_x, &solver);
            if (before_repair.satisfied == meta.num_clauses) break;
        }
        if (unsat) repair_active_cache(&cache, &solver);
        Verification current;
        if (!verify_assignment(&reader, &solver, &current)) { io_ok = 0; break; }
        retain_best(&current, &best_verification, &have_best, best_x, &solver);
        if (current.satisfied == meta.num_clauses) break;
    }

    Verification verification = {0};
    uint32_t finisher_passes_run = 0;
    if (io_ok && have_best && best_verification.satisfied < meta.num_clauses &&
        opt.finisher_passes && meta.num_clauses <= opt.finisher_max_clauses) {
        memcpy(solver.x, best_x, ((size_t)meta.num_vars + 1) * sizeof(double));
        io_ok = streaming_make_break_finisher(&reader, &solver, opt.finisher_passes,
                                              opt.finisher_batch_flips,
                                              &best_verification, &have_best, best_x,
                                              &finisher_passes_run);
    }
    if (io_ok && have_best) {
        memcpy(solver.x, best_x, ((size_t)meta.num_vars + 1) * sizeof(double));
        verification = best_verification;
    } else if (io_ok) io_ok = verify_assignment(&reader, &solver, &verification);
    if (io_ok && !write_solution(opt.solution_path, &solver)) {
        fprintf(stderr, "cannot write solution %s: %s\n", opt.solution_path, strerror(errno));
        io_ok = 0;
    }
    clock_gettime(CLOCK_MONOTONIC, &finished);
    int solved = io_ok && verification.satisfied == meta.num_clauses;
    int feasible = io_ok && verification.hard_unsatisfied == 0;
    double resident_mb = (((double)meta.num_vars + 1.0) *
                          (5.0 * sizeof(double) + 3.0 * sizeof(uint32_t)) +
                          ((double)opt.active_clauses + 1.0) * sizeof(uint32_t) +
                          (double)opt.active_clauses * (sizeof(uint64_t) + sizeof(uint8_t)) +
                          ((double)meta.num_vars + 1.0) * sizeof(double) +
                          3.0 * ((double)meta.num_vars + 1.0) * sizeof(double) +
                          (double)opt.finisher_batch_flips *
                          (sizeof(uint32_t) + sizeof(double)) +
                          (double)opt.active_literals * sizeof(int32_t)) / (1024.0 * 1024.0);
    if (opt.json) {
        printf("{\"solver\":\"NitroSAT V3\",\"format\":\"%s\",\"solved\":%s,"
               "\"feasible\":%s,\"variables\":%u,"
               "\"clauses\":%" PRIu64 ",\"literals\":%" PRIu64 ","
               "\"satisfied\":%" PRIu64 ",\"unsatisfied\":%" PRIu64 ","
               "\"hard_unsatisfied\":%" PRIu64 ",\"soft_unsatisfied\":%" PRIu64 ","
               "\"soft_cost\":%" PRIu64 ",\"total_soft_weight\":%" PRIu64 ","
               "\"cost_overflow\":%s,"
               "\"epochs\":%u,\"finisher_passes\":%u,\"increment_files\":%d,"
               "\"bounded_memory_mb\":%.2f,"
               "\"store_build_ms\":%.3f,\"solve_ms\":%.3f}\n",
               (meta.flags & STORE_FLAG_WEIGHTED) ? "wcnf" : "cnf",
               solved ? "true" : "false", feasible ? "true" : "false",
               meta.num_vars, meta.num_clauses, meta.num_literals,
               verification.satisfied, meta.num_clauses - verification.satisfied,
               verification.hard_unsatisfied, verification.soft_unsatisfied,
               verification.soft_cost, verification.total_soft_weight,
               verification.cost_overflow ? "true" : "false",
               epochs_run, finisher_passes_run, opt.add_count, resident_mb,
               elapsed_ms(&start, &built), elapsed_ms(&built, &finished));
    } else {
        printf("NitroSAT V3 (%s): %s\nVariables: %u  Clauses: %" PRIu64 "  Literals: %" PRIu64
               "\nSatisfied: %" PRIu64 "  Unsatisfied: %" PRIu64
               "\nHard unsatisfied: %" PRIu64 "  Soft cost: %" PRIu64
               "\nBounded solver allocation: %.2f MiB\n",
               (meta.flags & STORE_FLAG_WEIGHTED) ? "WCNF" : "CNF",
               solved ? "all clauses satisfied" :
               (weighted_store && feasible ? "weighted feasible assignment" :
                                               "best assignment (not proven UNSAT)"),
               meta.num_vars, meta.num_clauses, meta.num_literals, verification.satisfied,
               meta.num_clauses - verification.satisfied, verification.hard_unsatisfied,
               verification.soft_cost, resident_mb);
    }

    free(best_x);
    reader_close(&reader); solver_free(&solver); cache_free(&cache);
    if (!opt.keep_store) unlink(opt.store_path);
    free(opt.add_files);
    if (!io_ok) { fprintf(stderr, "clause store read failed\n"); return 2; }
    return (weighted_store ? feasible : solved) ? 0 : 1;
}
