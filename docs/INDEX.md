# NitroSAT: Pedagogical Slides & Code

This document pairs the cinematic educational slides with the core C implementation found in [nitrosat.c](nitrosat.c).

---

## Slide 0: Performance Overview
![Slide 0](slides/0.png)

**Concept**: NitroSAT is a high-performance linear-time MaxSAT approximator achieving 99.5%+ satisfaction on hard structured instances.

**Code Reference**: [nitrosat.c:1745-1753](nitrosat.c#L1745)
```c
    int final_sat = check_satisfaction(ns);
    printf("Satisfied: %d   Unsatisfied: %d\n", final_sat, ns->num_clauses - final_sat);
    printf("Solved  : %s\n", solved ? "YES" : "NO");
    printf("Time    : %.2f s\n", end - start);
```

---

## Slide 1: Boolean vs Energy Landscape
![Slide 1](slides/1.png)

**Concept**: Mapping discrete Boolean satisfaction into a continuous energy functional.

**Code Reference**: [nitrosat.c:384](nitrosat.c#L384)
```c
    double *x;          /* current (continuous) assignment [0, 1] */
```

---

## Slide 2: Continuous Relaxation
![Slide 2](slides/2.png)

**Concept**: Using Adelic geometry to relax variables over the reals.

**Code Reference**: [nitrosat.c:1448-1449](nitrosat.c#L1448)
```c
        for (int i = 1; i <= ns->num_vars; ++i)
            ns->x[i] = clamp(ns->x[i], 0.0, 1.0);
```

---

## Slide 3: Spectral Initialization
![Slide 3](slides/3.png)

**Concept**: XOR-Laplacian power iteration to find balanced starting points.

**Code Reference**: [nitrosat.c:603](nitrosat.c#L603)
```c
static void spectral_init(NitroSat *ns, int iters) {
    /* ... power iteration logic ... */
    for (int it = 0; it < iters; ++it) {
        /* y = D*v subtract off-diagonal ... */
    }
}
```

---

## Slide 4: Prime Weighting
![Slide 4](slides/4.png)

**Concept**: Number-theoretic weighting using primes to prioritize clause satisfaction.

**Code Reference**: [nitrosat.c:682-686](nitrosat.c#L682)
```c
    for (int c=0;c<ns->num_clauses;++c) {
        int p = (c < ns->prime_cnt) ? ns->primes[c] : 2;
        ns->cl_weights[c]   = 1.0 / pow(1.0 + log((double)p), 1.0);
        ns->zeta_weights[c] = log((double)p) / (double)p;
    }
```

---

## Slide 5: Log-Zeta Barrier
![Slide 5](slides/5.png)

**Concept**: The energy functional uses log-barrier terms to smooth the landscape.

**Code Reference**: [nitrosat.c:760-763](nitrosat.c#L760)
```c
        double barrier = 1.0 / (1e-6 + (1.0 - violation));
        double w = ns->cl_weights[c];
        if (ns->is_hard[c]) w *= 50.0;
        double coef = w * barrier * violation;
```

---

## Slide 6: Entropy & Momentum
![Slide 6](slides/6.png)

**Concept**: NADAM optimization paired with Adelic momentum and entropy regularization.

**Code Reference**: [nitrosat.c:165-170](nitrosat.c#L165)
```c
        if (is_nadam) {
            double m_nadam = m_hat + b1 * (m_hat - m_hat_prev);
            params[i] -= lr * m_nadam / (sqrt(v_hat) + eps);
        }
```

---

## Slide 7: BAHA & Lambert-W
![Slide 7](slides/7.png)

**Concept**: Branch-Aware Holonomy Annealing for escaping local minima via branch jumps.

**Code Reference**: [nitrosat.c:263](nitrosat.c#L263)
```c
static Branch *enumerate_branches(double bc, double b, int *out_cnt) {
    double u = b - bc;
    double xi = u * exp(u);
    double w0 = lambert_W0(xi);
    /* ... jumps to bc + w0 ... */
}
```

---

## Slide 8: Fracture Detection
![Slide 8](slides/8.png)

**Concept**: Detecting phase transitions in complexity to trigger jumps.

**Code Reference**: [nitrosat.c:225](nitrosat.c#L225)
```c
static int fd_is_fracture(FractureDetector *fd) {
    if (fd->n < 5) return 0;
    double s = sqrt(fd->M2 / (fd->n - 1));
    return fd->last_rate > fd->mean + fd->sigma_factor * s;
}
```

---

## Slide 9: Persistent Homology
![Slide 9](slides/9.png)

**Concept**: Analyzing the "holes" in the satisfaction manifold.

**Code Reference**: [nitrosat.c:485-491](nitrosat.c#L485)
```c
    topo.beta0 = uf->count - (num_vars - active_var_cnt);
    int chi = active_var_cnt - (int)edge_cnt;
    topo.beta1 = (topo.beta0 - chi > 0) ? (topo.beta0 - chi) : 0;
```

---

## Slide 10: Topological Repair
![Slide 10](slides/10.png)

**Concept**: Targeted local repair of topological defects (holes).

**Code Reference**: [nitrosat.c:1207](nitrosat.c#L1207)
```c
static int topological_repair_phase(NitroSat *ns, int max_steps) {
    /* Identifies and repairs manifold holes ... */
}
```

---

## Slide 11: Adelic Saturation
![Slide 11](slides/11.png)

**Concept**: Final multi-scale repairing phase using harmonic saturation.

**Code Reference**: [nitrosat.c:1372](nitrosat.c#L1372)
```c
static int adelic_saturation_phase(NitroSat *ns, int max_steps) {
    /* Aggressive saturation repair ... */
}
```

---

## Slide 12: Performance Results
![Slide 12](slides/12.png)

**Concept**: Verified O(M) scaling on benchmark suites.

**Code Reference**: [nitrosat.c:1425-1428](nitrosat.c#L1425)
```c
            printf("[%04d] sat=%d/%d (%.2f%%) lr=%.6f |grad|=%.4f\n",
                   step, sat, ns->num_clauses, 100.0*sat/ns->num_clauses,
                   ns->opt->lr, grad_norm);
```

---

## Slide 13: Summary & Conclusion
![Slide 13](slides/13.png)

**Concept**: Efficient physics-inspired optimization for complex MaxSAT problems.

**Code Reference**: [nitrosat.c:1764](nitrosat.c#L1764)
```c
    nitrosat_free(ns);
    return solved ? EXIT_SUCCESS : EXIT_FAILURE;
```
