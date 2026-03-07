--[[
NitroSat - Advanced MaxSAT Solver
---------------------------------
High-performance MaxSAT solver with Heat Kernel + BAHA optimization.
Developed with cutting-edge mathematical techniques for industrial applications.

Copyright 2026 Sethu Iyer (sethuiyer95@gmail.com)
Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

Key Features:
1.  Heat Kernel Smoothing: Diffuses local minima.
2.  BAHA: Fracture Detection & Branch Jumping.
]]

local math = require("math")
local os = require("os")
local ffi = require("ffi")
local has_json, json = pcall(require, "cjson")

-- FFI Helpers for massive scale
local function ffi_double(n) return ffi.new("double[?]", n + 1) end
local function ffi_int(n) return ffi.new("int32_t[?]", n + 1) end
local function ffi_uint8(n) return ffi.new("uint8_t[?]", n + 1) end

local function emit(event, fields)
    if not has_json then return end
    fields = fields or {}
    fields.event = event
    fields.component = "nitrosat"
    print(json.encode(fields))
end

-- ============================================================================
-- [INLINED] Optimizer Library
-- ============================================================================
local optimizer = {}
do
    local M = {}
    M.BaseOptimizer = {}
    M.BaseOptimizer.__index = M.BaseOptimizer
    function M.BaseOptimizer.new(params, config)
        local self = setmetatable({}, M.BaseOptimizer)
        self.params = params
        self.config = config or {}
        self.t = 0
        return self
    end
    
    M.HarmonicAdamOptimizer = {}
    setmetatable(M.HarmonicAdamOptimizer, {__index = M.BaseOptimizer})
    M.HarmonicAdamOptimizer.__index = M.HarmonicAdamOptimizer
    function M.HarmonicAdamOptimizer.new(params, config, nv)
        local self = setmetatable({}, M.HarmonicAdamOptimizer); self.params = params; self.config = config or {}; self.t = 0; self.nv = nv
        self.m = ffi_double(nv); self.v = ffi_double(nv)
        for i = 1, nv do self.m[i] = 0.0; self.v[i] = 0.0 end
        return self
    end
    function M.HarmonicAdamOptimizer:step(grads)
        self.t = self.t + 1
        local b1, b2, lr, eps, A0 = self.config.beta1 or 0.9, self.config.beta2 or 0.999, self.config.lr or 0.001, self.config.eps or 1e-8, self.config.resonance_amplitude or 0.02
        local n = self.nv
        local t_inv1, t_inv2 = 1/(1 - b1^self.t), 1/(1 - b2^self.t)
        for i = 1, n do
            self.m[i] = b1 * self.m[i] + (1 - b1) * grads[i]
            self.v[i] = b2 * self.v[i] + (1 - b2) * grads[i] * grads[i]
            self.params[i] = self.params[i] - lr * (self.m[i] * t_inv1) / (math.sqrt(self.v[i] * t_inv2) + eps)
        end
        local phi = (1 + math.sqrt(5)) / 2; local t_s = self.t * 0.1
        local intensity = math.sin(math.pi/20 * t_s) * math.exp(-math.pi/20 * t_s) + math.sin(math.pi/10 * t_s) * math.exp(-math.pi/10 * t_s) + math.sin(math.pi/9 * t_s) * math.exp(-math.pi/9 * t_s) + math.sin(1/(phi*phi) * t_s) * math.exp(-1/(phi*phi) * t_s)
        local g_norm = 0; for i = 1, n do g_norm = g_norm + grads[i] * grads[i] end; g_norm = math.sqrt(g_norm)
        if g_norm > 1e-10 then
            local pert = A0 * intensity
            for i = 1, n do self.params[i] = self.params[i] - pert * grads[i] / g_norm end
        end
        return self.params
    end

    -- NadamOptimizer: m_hat_prev lookahead (the one that got 36/316)
    M.NadamOptimizer = {}
    setmetatable(M.NadamOptimizer, {__index = M.BaseOptimizer})
    M.NadamOptimizer.__index = M.NadamOptimizer
    function M.NadamOptimizer.new(params, config, nv)
        local self = setmetatable({}, M.NadamOptimizer)
        self.params = params
        self.config = config or {}
        self.t = 0
        self.nv = nv
        self.m = ffi_double(nv)
        self.v = ffi_double(nv)
        for i = 1, nv do self.m[i] = 0.0; self.v[i] = 0.0 end
        return self
    end
    function M.NadamOptimizer:step(grads)
        self.t = self.t + 1
        local b1 = self.config.beta1 or 0.9
        local b2 = self.config.beta2 or 0.999
        local lr = self.config.lr or 0.001
        local eps = self.config.eps or 1e-8
        local A0 = self.config.resonance_amplitude or 0.02
        local n = self.nv

        local b1_t = b1 ^ self.t
        local b2_t = b2 ^ self.t
        local t_inv1 = 1 / (1 - b1_t)
        local t_inv2 = 1 / (1 - b2_t)

        -- m_hat_prev lookahead
        local m_hat_prev = ffi_double(n)
        for i = 1, n do
            m_hat_prev[i] = b1 * self.m[i] * t_inv1
        end

        -- Update momentum and velocity
        for i = 1, n do
            self.m[i] = b1 * self.m[i] + (1 - b1) * grads[i]
            self.v[i] = b2 * self.v[i] + (1 - b2) * grads[i] * grads[i]
        end

        -- Bias-corrected estimates
        local m_hat = ffi_double(n)
        local v_hat = ffi_double(n)
        for i = 1, n do
            m_hat[i] = self.m[i] * t_inv1
            v_hat[i] = self.v[i] * t_inv2
        end

        -- Nesterov momentum correction
        for i = 1, n do
            local m_nadam = m_hat[i] + b1 * (m_hat[i] - m_hat_prev[i])
            local denom = math.sqrt(v_hat[i]) + eps
            self.params[i] = self.params[i] - lr * m_nadam / denom
        end

        -- Resonance perturbation
        local phi = (1 + math.sqrt(5)) / 2
        local t_s = self.t * 0.1
        local intensity = math.sin(math.pi/20 * t_s) * math.exp(-math.pi/20 * t_s)
                     + math.sin(math.pi/10 * t_s) * math.exp(-math.pi/10 * t_s)
                     + math.sin(math.pi/9 * t_s) * math.exp(-math.pi/9 * t_s)
                     + math.sin(1/(phi*phi) * t_s) * math.exp(-1/(phi*phi) * t_s)

        local g_norm = 0
        for i = 1, n do g_norm = g_norm + grads[i] * grads[i] end
        g_norm = math.sqrt(g_norm)

        if g_norm > 1e-10 then
            local pert = A0 * intensity
            for i = 1, n do
                self.params[i] = self.params[i] - pert * grads[i] / g_norm
            end
        end
        return self.params
    end

    M.HarmonicNadamOptimizer = {}
    setmetatable(M.HarmonicNadamOptimizer, {__index = M.BaseOptimizer})
    M.HarmonicNadamOptimizer.__index = M.HarmonicNadamOptimizer
    function M.HarmonicNadamOptimizer.new(params, config, nv)
        local self = setmetatable({}, M.HarmonicNadamOptimizer); self.params = params; self.config = config or {}; self.t = 0; self.nv = nv
        self.m = ffi_double(nv); self.v = ffi_double(nv)
        for i = 1, nv do self.m[i] = 0.0; self.v[i] = 0.0 end
        return self
    end
    function M.HarmonicNadamOptimizer:step(grads)
        self.t = self.t + 1
        local b1, b2, lr, eps, A0 = self.config.beta1 or 0.9, self.config.beta2 or 0.999, self.config.lr or 0.002, self.config.eps or 1e-8, self.config.resonance_amplitude or 0.02
        local n = self.nv
        local t = self.t
        local b1_t = b1^t
        local b1_t_next = b1^(t+1)
        local b2_t = b2^t
        
        for i = 1, n do
            local g = grads[i]
            self.m[i] = b1 * self.m[i] + (1 - b1) * g
            self.v[i] = b2 * self.v[i] + (1 - b2) * g * g
            
            local m_hat = (b1 * self.m[i] / (1 - b1_t_next)) + ((1 - b1) * g / (1 - b1_t))
            local v_hat = self.v[i] / (1 - b2_t)
            
            self.params[i] = self.params[i] - lr * m_hat / (math.sqrt(v_hat) + eps)
        end
        
        local phi = (1 + math.sqrt(5)) / 2; local t_s = self.t * 0.1
        local intensity = math.sin(math.pi/20 * t_s) * math.exp(-math.pi/20 * t_s) + math.sin(math.pi/10 * t_s) * math.exp(-math.pi/10 * t_s) + math.sin(math.pi/9 * t_s) * math.exp(-math.pi/9 * t_s) + math.sin(1/(phi*phi) * t_s) * math.exp(-1/(phi*phi) * t_s)
        local g_norm = 0; for i = 1, n do g_norm = g_norm + grads[i] * grads[i] end; g_norm = math.sqrt(g_norm)
        if g_norm > 1e-10 then
            local pert = A0 * intensity
            for i = 1, n do self.params[i] = self.params[i] - pert * grads[i] / g_norm end
        end
        return self.params
    end

    function M.create_optimizer(name, params, config, nv)
        if name == "nadam" then return M.NadamOptimizer.new(params, config, nv) end
        if name == "harmonic_adam" then return M.HarmonicAdamOptimizer.new(params, config, nv) end
        if name == "harmonic_nadam" then return M.HarmonicNadamOptimizer.new(params, config, nv) end
        error("Optimizer not found: " .. name)
    end
    optimizer = M
end


-- ============================================================================
-- [INLINED] BAHA Utils Library
-- ============================================================================
local baha_utils = {}
do
    local M = {}
    local E_INV = 0.36787944117144232
    function M.log_sum_exp(log_terms, nsamp)
        local mt = -math.huge; for i=1,nsamp do if log_terms[i] > mt then mt = log_terms[i] end end
        if mt == -math.huge then return -math.huge end
        local sum = 0.0; for i=1,nsamp do sum = sum + math.exp(log_terms[i] - mt) end
        return mt + math.log(sum)
    end
    local function halley(z, w)
        for i = 1, 50 do
            local ew = math.exp(w); local f = w * ew - z; local fp = ew * (w + 1.0)
            if math.abs(fp) < 1e-15 then break end
            local denom = fp - f * (ew * (w + 2.0)) / (2.0 * fp)
            if math.abs(denom) < 1e-15 then break end
            local w_new = w - f / denom; if math.abs(w_new - w) < 1e-10 then return w_new end; w = w_new
        end
        return w
    end
    M.LambertW = {
        W0 = function(z)
            if z < -E_INV then return 0/0 end
            local w = (z < -0.3) and (z * math.exp(1.0)) or ((z < 1.0) and (z * (1.0 - z + z * z)) or (math.log(z) - math.log(math.log(z) + 1.0)))
            return halley(z, w)
        end,
        Wm1 = function(z)
            if z < -E_INV or z >= 0.0 then return 0/0 end
            return halley(z, math.log(-z) - math.log(-math.log(-z)))
        end
    }
    local FD = {}; FD.__index = FD
    function M.FractureDetector_new(sf)
        return setmetatable({sigma_factor = sf or 2.0, n=0, mean=0.0, M2=0.0, last_rate=0.0}, FD)
    end
    function FD:record(beta, lZ)
        if self.pb and self.plZ then
            local db = beta - self.pb
            if math.abs(db) > 1e-9 then
                local rate = math.abs(lZ - self.plZ) / db; self.last_rate = rate; self.n = self.n + 1
                local d = rate - self.mean; self.mean = self.mean + d / self.n; self.M2 = self.M2 + d * (rate - self.mean)
            end
        end
        self.pb, self.plZ = beta, lZ
    end
    function FD:is_fracture()
        if self.n < 5 then return false end
        local s = math.sqrt(self.M2 / (self.n - 1))
        return self.last_rate > self.mean + self.sigma_factor * s
    end
    function FD:clear() self.n, self.mean, self.M2, self.last_rate, self.pb, self.plZ = 0, 0.0, 0.0, 0.0, nil, nil end
    function M.enumerate_branches(bc, b)
        local br = {}; local u = b - bc; if math.abs(u) < 1e-10 then u = 1e-10 end; local xi = u * math.exp(u)
        local w0 = M.LambertW.W0(xi); if w0 == w0 then local b0 = bc + w0; if b0 > 0 then table.insert(br, {k=0, beta=b0, score=0.0}) end end
        if xi >= -E_INV and xi < 0 then local wm1 = M.LambertW.Wm1(xi); if wm1 == wm1 then local bm1 = bc + wm1; if bm1 > 0 then table.insert(br, {k=-1, beta=bm1, score=0.0}) end end end
        return br
    end
    baha_utils = M
end

-- ============================================================================
local xor_laplacian = {}
do
    local M = {}
    function M.matvec(num_clauses, c_flat, c_offsets, degs, x, nv)
        local y = ffi_double(nv)
        for i=1,nv do y[i] = degs[i]*x[i] end
        for c=0,num_clauses-1 do
            local start = c_offsets[c]
            local stop = c_offsets[c+1]
            local sum = 0.0
            for i=start,stop-1 do sum = sum + x[math.abs(c_flat[i])] end
            for i=start,stop-1 do
                local vi = math.abs(c_flat[i])
                y[vi] = y[vi] - (sum - x[vi])
            end
        end
        return y
    end
    function M.spectral_init(num_clauses, c_flat, c_offsets, degs, nv, iters, rng_fn)
        local v = ffi_double(nv)
        for i=1,nv do v[i] = rng_fn() end
        for it = 1, (iters or 20) do
            v = M.matvec(num_clauses, c_flat, c_offsets, degs, v, nv)
            local n = 0.0; for i=1,nv do n = n + v[i]*v[i] end; n = math.sqrt(n)
            if n > 1e-10 then for i=1,nv do v[i] = v[i]/n end end
        end
        local mi, ma = math.huge, -math.huge; for i=1,nv do if v[i]<mi then mi=v[i] end; if v[i]>ma then ma=v[i] end end
        local r = ma - mi; if r>1e-10 then for i=1,nv do v[i] = (v[i]-mi)/r end end; return v
    end
    xor_laplacian = M
end

-- ============================================================================
-- [NEW] Persistent Homology Module for Traceability
-- ============================================================================
local persistent_homology = {}
do
    local M = {}

    -- Union-Find (Disjoint Set Union) for connected components
    local UnionFind = {}
    UnionFind.__index = UnionFind
    function UnionFind.new(n)
        local self = setmetatable({}, UnionFind)
        self.parent = ffi_int(n + 1)
        self.rank = ffi_int(n + 1)
        self.count = n
        for i = 1, n do
            self.parent[i] = i
            self.rank[i] = 0
        end
        return self
    end
    function UnionFind:find(x)
        if self.parent[x] ~= x then
            self.parent[x] = self:find(self.parent[x])  -- Path compression
        end
        return self.parent[x]
    end
    function UnionFind:union(x, y)
        local px, py = self:find(x), self:find(y)
        if px == py then return false end
        -- Union by rank
        if self.rank[px] < self.rank[py] then
            px, py = py, px
        elseif self.rank[px] == self.rank[py] then
            self.rank[px] = self.rank[px] + 1
        end
        self.parent[py] = px
        self.count = self.count - 1
        return true
    end
    function UnionFind:num_components()
        return self.count
    end

    M.UnionFind = UnionFind

    -- Compute Betti numbers from variable similarity graph
    -- beta0: number of connected components
    -- beta1: number of independent cycles (holes)
    function M.compute_betti_numbers(num_clauses, c_flat, c_offsets, x, num_vars, threshold)
        threshold = threshold or 0.3
        local adj = {}
        local degrees = ffi_double(num_vars)

        for c = 0, num_clauses - 1 do
            local start = c_offsets[c]
            local stop = c_offsets[c+1]
            local is_sat = false
            
            for i = start, stop - 1 do
                local lit = c_flat[i]
                local var = math.abs(lit)
                if (lit > 0 and x[var] > 0.5) or (lit < 0 and x[var] <= 0.5) then
                    is_sat = true
                    break
                end
            end

            if not is_sat then
                for i = start, stop - 1 do
                    local vi = math.abs(c_flat[i])
                    for j = i + 1, stop - 1 do
                        local vj = math.abs(c_flat[j])
                        adj[vi] = adj[vi] or {}
                        adj[vj] = adj[vj] or {}
                        local w = (adj[vi][vj] or 0) + 1
                        adj[vi][vj] = w
                        adj[vj][vi] = w
                        degrees[vi] = degrees[vi] + 1
                        degrees[vj] = degrees[vj] + 1
                    end
                end
            end
        end

        local edges = {}
        for i = 1, num_vars do
            if adj[i] then
                for j, w in pairs(adj[i]) do
                    if i < j then
                        local similarity = w / math.max(degrees[i], degrees[j], 1)
                        table.insert(edges, {i, j, similarity})
                    end
                end
            end
        end
        table.sort(edges, function(a, b) return a[3] > b[3] end)

        local uf = UnionFind.new(num_vars)
        local edge_count = 0
        for _, edge in ipairs(edges) do
            if edge[3] >= threshold then
                uf:union(edge[1], edge[2])
                edge_count = edge_count + 1
            end
        end

        local beta0 = uf:num_components()
        local active_vars = 0
        for i = 1, num_vars do if degrees[i] > 0 then active_vars = active_vars + 1 end end
        local chi = active_vars - edge_count
        local beta1 = math.max(0, beta0 - chi)
        local solution_density = 0
        local active_count = 0
        for i = 1, num_vars do
            if degrees[i] > 0 then
                solution_density = solution_density + math.abs(x[i] - 0.5)
                active_count = active_count + 1
            end
        end
        local avg_confidence = (active_count > 0) and (solution_density / active_count) or 0
        return {beta0=beta0, beta1=beta1, active_vars=active_vars, edge_count=edge_count, avg_confidence=avg_confidence, complexity_score=beta1/math.max(beta0,1)}
    end

    -- Track persistence of topological features across steps
    function M.create_persistence_tracker(num_vars)
        return {
            beta0_history = {},
            beta1_history = {},
            complexity_history = {},
            birth_step = {},
            persistence_pairs = {},  -- {birth, death, dimension}
            current_dim = 0
        }
    end

    function M.update_persistence(tracker, beta0, beta1, step, complexity)
        table.insert(tracker.beta0_history, beta0)
        table.insert(tracker.beta1_history, beta1)
        if complexity then
            table.insert(tracker.complexity_history, complexity)
        end

        -- Track birth/death of features
        local prev_beta0 = tracker.beta0_history[#tracker.beta0_history - 1] or beta0
        local prev_beta1 = tracker.beta1_history[#tracker.beta1_history - 1] or beta1

        -- Feature born
        if beta0 > prev_beta0 then
            tracker.birth_step[0] = step
        end
        if beta1 > prev_beta1 then
            tracker.birth_step[1] = step
        end

        -- Feature died (persistence = death - birth)
        if beta0 < prev_beta0 and tracker.birth_step[0] then
            local persistence = step - tracker.birth_step[0]
            table.insert(tracker.persistence_pairs, {0, persistence, step})
        end
        if beta1 < prev_beta1 and tracker.birth_step[1] then
            local persistence = step - tracker.birth_step[1]
            table.insert(tracker.persistence_pairs, {1, persistence, step})
        end
    end

    persistent_homology = M
end

local function clamp(x, min_val, max_val)
    if x < min_val then return min_val end
    if x > max_val then return max_val end
    return x
end

local function math_sign(x)
    if x > 0 then return 1 elseif x < 0 then return -1 else return 0 end
end

-- Sieve of Eratosthenes (Ported from Pro)
local function sieve_primes(n)
    local limit = math.floor(n * math.log(n) * 1.2) + 100
    if n <= 10 then limit = 100 end
    
    local sieve = {}
    for i = 0, limit do sieve[i] = true end
    sieve[0], sieve[1] = false, false
    
    for i = 2, math.floor(math.sqrt(limit)) do
        if sieve[i] then
            for j = i*i, limit, i do
                sieve[j] = false
            end
        end
    end
    
    local primes = {}
    local count = 0
    for i = 2, limit do
        if sieve[i] then
            primes[count] = i
            count = count + 1
            if count >= n then break end
        end
    end
    
    return primes
end

-- Zeta zero guided perturbation (Ported from Pro)
-- Uses critical line Re(s) = 1/2 and prime harmonics
local function zeta_zero_perturbation(var_idx, step, max_steps, primes, prime_count)
    local progress = step / max_steps
    local critical_line = 0.5  -- Re(s) = 1/2 (Riemann Hypothesis!)

    local zero_guidance = 0
    for i = 0, prime_count - 1 do  -- Adaptive prime count for symmetry breaking
        local p = primes[i]
        if p then
            local zeta_map = math.sin(2 * math.pi * (p / var_idx)) * progress
            zero_guidance = zero_guidance + zeta_map * (1 / math.sqrt(i + 1))
        end
    end

    local critical_osc = math.sin(2 * math.pi * critical_line * progress) * 0.1

    return zero_guidance + critical_osc
end

-- ============================================================================
-- NUCLEAR ZETA MODE: Adelic force saturation with critical line resonance
-- ============================================================================
local function zeta_zero_perturbation_nuclear(var_idx, step, max_steps, primes, topology, num_vars, prime_count)
    local progress = step / max_steps
    local critical_line = 0.5  -- Re(s) = 1/2

    -- Base: Riemann-Siegel theta approximation
    local theta = 0
    for i = 0, math.min(prime_count - 1, #primes) do
        local p = primes[i]
        if p then
            -- Zeta zero alignment: sin(2π · p · var_idx / log(p))
            local phase = 2 * math.pi * (p * var_idx) / (math.log(p) * num_vars)
            local amplitude = 1 / math.sqrt(p)  -- Prime reciprocity law

            -- Critical line oscillation
            local critical_osc = math.sin(2 * math.pi * critical_line * progress * p)

            theta = theta + amplitude * math.sin(phase) * critical_osc
        end
    end

    -- Topology-guided amplification: boost in high-complexity regions
    if topology and topology.complexity_score and topology.complexity_score > 0.5 then
        theta = theta * (1 + topology.complexity_score)
    end

    -- Phase transition resonance near convergence
    if progress > 0.9 then
        theta = theta * 2.0  -- Amplify near endgame
    end

    return theta
end

-- ============================================================================
-- SHUNYABAR SPECTRAL-ARITHMETIC CORE
-- ============================================================================

-- Logarithmic derivative of Riemann Zeta function (von Mangoldt force)
-- -ζ'(s)/ζ(s) = Σ Λ(n) n^{-s}
local function zeta_log_derivative(s, primes, num_terms)
    local force = 0.0
    for i = 0, num_terms - 1 do
        local p = primes[i]
        if not p then break end
        -- Λ(p^k) = log p. For square-free/prime approximation:
        force = force + math.log(p) * math.pow(p, -s)
    end
    return force
end

-- Spectral-Arithmetic Partition Function Gradient
-- ∇ log Z(s) = ∇ log ζ(s) + ∇ log Tr(e^{-s L})
-- At s=1 (Arithmetic Pole), the ζ(s) term dominates the symmetry breaking.
local function compute_shunyabar_force(x, s, primes, l_weights, num_vars)
    local forces = ffi_double(num_vars)
    for i = 1, num_vars do forces[i] = 0.0 end
    
    -- Arithmetic component: von Mangoldt gain at s
    local gain = zeta_log_derivative(s, primes, 100)
    
    -- Spectral component: Geometric instability from Laplacian L
    -- F_spectral = -s * ∇ Tr(e^{-s L}) / Tr(e^{-s L})
    -- We approximate this using the prime-weighted signatures
    for i = 1, num_vars do
        forces[i] = gain * (l_weights[i] or 1.0) * (x[i] - 0.5)
    end
    
    return forces
end

-- ==============================================================================
-- NitroSat Solver
-- ==============================================================================

local NitroSat = {}
NitroSat.__index = NitroSat

function NitroSat.new(instance, opts)
    opts = opts or {}
    local t_init_start = os.clock()

    local self = setmetatable({}, NitroSat)
    self.num_vars = instance.num_vars
    self.num_clauses = #instance.clauses
    -- Flatten clauses into FFI arrays for massive scale (avoids Lua heap limits)
    local total_literals = 0
    for c = 1, self.num_clauses do
        total_literals = total_literals + #instance.clauses[c]
    end
    
    self.clauses_flat = ffi_int(total_literals)
    self.clauses_offsets = ffi_int(self.num_clauses + 1)
    
    local curr_idx = 0
    for c = 1, self.num_clauses do
        self.clauses_offsets[c-1] = curr_idx -- 0-indexed for FFI
        local clause = instance.clauses[c]
        for i = 1, #clause do
            self.clauses_flat[curr_idx] = clause[i]
            curr_idx = curr_idx + 1
        end
    end
    self.clauses_offsets[self.num_clauses] = curr_idx
    
    -- We can now drop the Lua tables to free heap memory
    self.clauses = nil 
    self.trace_id = opts.trace_id

    emit("solver_init", {
        duration_ms = (os.clock() - t_init_start) * 1000,
        num_vars = self.num_vars,
        num_clauses = self.num_clauses,
        trace_id = self.trace_id
    })
    
    self.seed = opts.seed or 42
    
    -- Heat kernel parameters
    self.heat_beta = opts.heat_beta or 0.5
    self.heat_lambda = opts.heat_lambda or 0.1
    self.use_heat_kernel = true -- Always enabled in NitroSat

    -- Prime-weighted log-barrier parameter (paper Eq. 11)
    self.prime_alpha = opts.prime_alpha or 1.0

    -- BAHA Configuration
    self.use_baha = true -- Always enabled in NitroSat
    self.dec_freq = opts.dec_freq or 100 
    self.decimated = ffi_uint8(self.num_vars)
    for i = 1, self.num_vars do self.decimated[i] = 0 end

    math.randomseed(self.seed)
    self._rng_state = math.random()
    if self._rng_state == 0 then self._rng_state = 0.1 end
    
    local sl = os.clock()
    self.degrees = ffi_double(self.num_vars)
    -- Re-calculate degrees using flattened clauses
    for c_idx = 0, self.num_clauses - 1 do
        local start_idx = self.clauses_offsets[c_idx]
        local end_idx = self.clauses_offsets[c_idx+1]
        local k = end_idx - start_idx -- Length of clause
        for i = start_idx, end_idx - 1 do
            local lit = self.clauses_flat[i]
            self.degrees[math.abs(lit)] = self.degrees[math.abs(lit)] + (k - 1)
        end
    end
    emit("degrees_built", {duration_ms = (os.clock() - sl) * 1000})
    
    -- Hard Constraint Logic
    -- For now, treat all clauses as hard unless marked otherwise.
    -- In SAT, usually unit clauses or specific structures differ.
    -- Pro uses detailed logic, here we'll just implement the mechanism:
    -- Hard constraints get 50x weight and priority in separation.
    -- Zeta Zero Magic (The Final Piece)
    self.zeta_guidance = opts.zeta_guidance ~= false -- Default to TRUE
    self.primes = sieve_primes(self.num_clauses)

    -- Adaptive prime count: 3 * log(n) * log(α+1), clamped [10, 256]
    local alpha = (self.num_vars > 0) and (self.num_clauses / self.num_vars) or 1
    self.prime_count = math.max(10, math.min(256,
        math.floor(3 * math.log(self.num_vars > 0 and self.num_vars or 1) * math.log(alpha + 1))))

    self.is_hard_constraint = ffi_uint8(self.num_clauses)
    self.clause_weights = ffi_double(self.num_clauses)
    for i = 0, self.num_clauses - 1 do
        self.is_hard_constraint[i] = 1 -- Default to Hard (1 for true)
        -- Prime-weighted base: w_c = (1 + log p_c)^{-alpha} (paper Eq. 11)
        local p = self.primes[i] or 2
        self.clause_weights[i] = 1.0 / math.pow(1 + math.log(p), self.prime_alpha)
    end

    -- [NEW] Persistent Homology Configuration
    self.use_topology = opts.use_topology ~= false -- Default to TRUE
    self.topology_freq = opts.topology_freq or 50  -- Compute every N steps
    self.topology_threshold = opts.topology_threshold or 0.3
    self.persistence_tracker = persistent_homology.create_persistence_tracker(self.num_vars)
    self.last_topology = nil
    
    -- Pre-calculate Zeta weights for high-precision sweep
    -- w_c = log(p_c) / p_c (Logarithmic derivative components)
    self.zeta_weights = ffi_double(self.num_clauses)
    for i = 0, self.num_clauses - 1 do
        local p = self.primes[i]
        if p then self.zeta_weights[i+1] = math.log(p) / p end
    end
    
    -- Initialization
    if opts.x_init then
        self.x = ffi_double(self.num_vars)
        for i = 1, #opts.x_init do self.x[i] = opts.x_init[i] end
    else
        -- Spectral Initialization (Heat Kernel pre-conditioning)
        local rng = function() return self:_rng() end
        local v_init = xor_laplacian.spectral_init(self.num_clauses, self.clauses_flat, self.clauses_offsets, self.degrees, self.num_vars, 50, rng)
        self.x = ffi_double(self.num_vars)
        for i = 1, self.num_vars do
            self.x[i] = v_init[i] + (self:_rng() - 0.5) * 0.01
        end
    end
    
    self.m = ffi_double(self.num_vars)
    self.v = ffi_double(self.num_vars)
    for i = 1, self.num_vars do self.m[i] = 0.0; self.v[i] = 0.0 end
    
    self.baha_threshold = opts.baha_threshold or 1.5
    self.baha_detector = baha_utils.FractureDetector_new(self.baha_threshold)
    self.baha_history_size = 50 
    self.last_energies = {}

    -- BAHA-WalkSAT configuration
    self.bw_beta_start    = opts.bw_beta_start    or 1.0
    self.bw_beta_end      = opts.bw_beta_end      or 20.0
    self.bw_fracture_freq = opts.bw_fracture_freq or 50
    self.bw_threshold     = opts.bw_threshold     or 1.5
    self.bw_n_samples     = opts.bw_n_samples     or 20
    self.bw_var_spike_z   = opts.bw_var_spike_z   or 0.8
    self.bw_var_jump_gain = opts.bw_var_jump_gain or 0.25
    self.baha_var_spike_z = opts.baha_var_spike_z or 0.8
    
    self.entropy_weight = opts.entropy_weight or 0.01
    self.dec_threshold = opts.dec_threshold or 0.4

    -- Nuclear Zeta Mode
    self.nuclear_zeta = opts.nuclear_zeta or false
    self.nuclear_threshold = opts.nuclear_threshold or 0.995

    -- Three-Phase Finisher Configuration
    self.phase2_start = opts.phase2_start or 0.95
    self.phase3_start = opts.phase3_start or 0.98
    
    self.optimizer = optimizer.create_optimizer('nadam', self.x, {
        lr = 0.002,
        beta1 = 0.9,
        beta2 = 0.999,
        resonance_amplitude = 0.02
    }, self.num_vars)
    
    -- Var map (CSR format for memory efficiency)
    local v_counts = ffi_int(self.num_vars + 1)
    for i=0,self.num_vars do v_counts[i] = 0 end
    
    for i=0, self.clauses_offsets[self.num_clauses]-1 do
        local v = math.abs(self.clauses_flat[i])
        v_counts[v] = v_counts[v] + 1
    end
    
    self.v2c_ptr = ffi_int(self.num_vars + 2)
    self.v2c_ptr[1] = 0
    for i=1,self.num_vars do self.v2c_ptr[i+1] = self.v2c_ptr[i] + v_counts[i] end
    
    self.v2c_data = ffi_int(self.clauses_offsets[self.num_clauses])
    local v_cur = ffi_int(self.num_vars + 1)
    for i=0,self.num_vars do v_cur[i] = 0 end
    
    for c=0, self.num_clauses - 1 do
        local start = self.clauses_offsets[c]
        local stop = self.clauses_offsets[c+1]
        for i=start, stop-1 do
            local v = math.abs(self.clauses_flat[i])
            self.v2c_data[self.v2c_ptr[v] + v_cur[v]] = c
            v_cur[v] = v_cur[v] + 1
        end
    end
    
    return self
end

-- Fast deterministic PRNG: x = (x * 9973) mod 1
function NitroSat:_rng()
    local x = self._rng_state * 9973
    x = x - math.floor(x)
    self._rng_state = x
    return x
end

function NitroSat:_rng_int(a, b)
    return a + math.floor(self:_rng() * (b - a + 1))
end

-- Compute heat kernel multipliers
function NitroSat:compute_heat_multipliers()
    local multipliers = ffi_double(self.num_vars)
    -- Implicit Laplacian for heat kernel
    local x_sums = {} -- Optimization: compute vertex-clique-sums first
    -- This is actually harder implicitly. For now, let's use a simplified
    -- heat kernel that only considers degrees if we are at this scale.
    for i = 1, self.num_vars do
        local heat_weight = math.exp(-self.heat_beta * (self.degrees[i] or 0))
        multipliers[i] = 1 + self.heat_lambda * heat_weight
    end
    return multipliers
end

-- [NEW] Compute topological features (Betti numbers)
function NitroSat:compute_topology()
    if not self.use_topology then return nil end

    local topology = persistent_homology.compute_betti_numbers(
        self.num_clauses,
        self.clauses_flat,
        self.clauses_offsets,
        self.x,
        self.num_vars,
        self.topology_threshold
    )

    -- Update persistence tracker
    persistent_homology.update_persistence(
        self.persistence_tracker,
        topology.beta0,
        topology.beta1,
        self.current_step or 0,
        topology.complexity_score
    )

    self.last_topology = topology
    return topology
end

-- [NEW] Zeta Sweep Finisher: High-precision stochastic repair
function NitroSat:zeta_sweep(beta)
    beta = beta or 1.0
    local dF = ffi_double(self.num_vars)
    for i=1,self.num_vars do dF[i] = 0.0 end
    
    local unsat_clauses = {}
    local unsat_count = 0
    
    -- Identify unsat clauses efficiently
    for c = 0, self.num_clauses - 1 do
        local is_sat = false
        local start = self.clauses_offsets[c]
        local stop = self.clauses_offsets[c+1]
        for i = start, stop - 1 do
            local lit = self.clauses_flat[i]
            local var = math.abs(lit)
            if (lit > 0 and self.x[var] > 0.5) or (lit < 0 and self.x[var] <= 0.5) then
                is_sat = true
                break
            end
        end
        if not is_sat then
            unsat_count = unsat_count + 1
            unsat_clauses[unsat_count] = c
        end
    end
    
    if unsat_count == 0 then return true, 0 end
    
    -- Apply Adelic prime-forces to variables in unsat clauses
    for idx = 1, unsat_count do
        local c = unsat_clauses[idx]
        local w = self.zeta_weights[c+1] or 0.1
        local start = self.clauses_offsets[c]
        local stop = self.clauses_offsets[c+1]
        for i = start, stop - 1 do
            local lit = self.clauses_flat[i]
            local var = math.abs(lit)
            local force = (lit > 0 and 1 or -1) * w
            dF[var] = dF[var] + force
        end
    end
    
    local flips = 0
    for i = 1, self.num_vars do
        if math.abs(dF[i]) > 1e-9 then
            -- Flip probability proportional to Zeta force intensity
            local p = 1 / (1 + math.exp(-beta * math.abs(dF[i])))
            if self:_rng() < p then
                self.x[i] = (dF[i] > 0) and 0.99 or 0.01
                flips = flips + 1
            end
        end
    end
    
    return false, unsat_count, flips
end

-- Compute Gradients (Standard + Heat Kernel Smoothing)
function NitroSat:compute_gradients()
    local grads = ffi_double(self.num_vars)
    for i = 1, self.num_vars do grads[i] = 0.0 end
    local unsat_count = 0
    
    for c = 0, self.num_clauses - 1 do
        local start = self.clauses_offsets[c]
        local stop = self.clauses_offsets[c+1]
        
        local violation = 1.0
        local is_discrete_sat = false
        
        for i = start, stop - 1 do
            local lit = self.clauses_flat[i]
            local var = math.abs(lit)
            local val = self.x[var]
            
            if (lit > 0 and val > 0.5) or (lit < 0 and val <= 0.5) then is_discrete_sat = true end
            
            local lit_val = lit > 0 and (1.0 - val) or val
            violation = violation * lit_val
        end
        
        if not is_discrete_sat then unsat_count = unsat_count + 1 end
        
        -- Log-barrier: gradient of -w_c * log(eps + s_C) where s_C = 1 - v_C
        local slack = 1.0 - violation
        local barrier = 1.0 / (1e-6 + slack)
        
        local coef = self.clause_weights[c]
        if self.is_hard_constraint[c] ~= 0 then
            coef = coef * 50.0
        end
        
        for i = start, stop - 1 do
            local lit = self.clauses_flat[i]
            local var = math.abs(lit)
            local val = self.x[var]
            local lit_val = lit > 0 and (1.0 - val) or val
            if math.abs(lit_val) < 1e-9 then lit_val = 1e-9 end
            
            local term_grad = (violation / lit_val) * (lit > 0 and -1 or 1)
            grads[var] = grads[var] + coef * barrier * term_grad
        end
    end
    
    -- HEAT KERNEL SMOOTHING
    local multipliers = self:compute_heat_multipliers()
    for i = 1, self.num_vars do
        grads[i] = grads[i] * multipliers[i]
        
        -- Entropy Regularization Gradient
        -- H(x) = -sum(x log x + (1-x) log(1-x))
        -- dH/dx = -log(x) + log(1-x) = log((1-x)/x)
        -- We want to MAXIMIZE entropy, so we ADD gradient of H.
        -- grad += weight * log((1-x)/x)
        local val = self.x[i]
        val = clamp(val, 1e-9, 1.0 - 1e-9)
        local entropy_grad = math.log((1.0 - val) / val)
        grads[i] = grads[i] + self.entropy_weight * entropy_grad
    end
    
    return grads, unsat_count
end

function NitroSat:adam_step(grads, lr)
    if lr then
        self.optimizer.config.lr = lr
    end
    self.x = self.optimizer:step(grads)
    for i = 1, self.num_vars do self.x[i] = clamp(self.x[i], 0.0, 1.0) end
end

-- Oscillating annealing schedule (Ported from Pro)
function NitroSat:get_annealing_lr(t, A0)
    A0 = A0 or 0.1
    local phi = (1 + math.sqrt(5)) / 2
    
    local term_solid  = math.sin(math.pi/20 * t) * math.exp(-math.pi/20 * t)
    local term_liquid = math.sin(math.pi/10 * t) * math.exp(-math.pi/10 * t)
    local term_gas    = math.sin(math.pi/9 * t)  * math.exp(-math.pi/9 * t)
    local term_golden = math.sin(1/(phi*phi) * t) * math.exp(-1/(phi*phi) * t)
    
    local intensity = term_solid + term_liquid + term_gas + term_golden
    return math.max(0.0, A0 * intensity)
end

function NitroSat:check_satisfaction()
    local sat_count = 0
    for c = 0, self.num_clauses - 1 do
        local satisfied = false
        local start = self.clauses_offsets[c]
        local stop = self.clauses_offsets[c+1]
        for i = start, stop - 1 do
            local lit = self.clauses_flat[i]
            local var = math.abs(lit)
            local val = self.x[var]
            if (lit > 0 and val > 0.5) or (lit < 0 and val <= 0.5) then
                satisfied = true
                break
            end
        end
        if satisfied then sat_count = sat_count + 1 end
    end
    return sat_count
end

-- BAHA Core Functions (Ported from C++)
-- Optimized BAHA: pre-allocate sample buffer, localize math, avoid table churn
local _abs = math.abs
local _exp = math.exp
local _log = math.log

function NitroSat:_fill_sample(buf)
    -- Fill buffer with random state, respecting decimation
    local x, dec = self.x, self.decimated
    for j = 1, self.num_vars do
        buf[j] = (dec[j] ~= 0) and x[j] or self:_rng()
    end
end

function NitroSat:calculate_energy(x)
    local unsat = 0
    for c = 0, self.num_clauses - 1 do
        local sat = false
        local start = self.clauses_offsets[c]
        local stop = self.clauses_offsets[c+1]
        for i = start, stop - 1 do
            local lit = self.clauses_flat[i]
            local v = math.abs(lit)
            if (lit > 0 and x[v] > 0.5) or (lit < 0 and x[v] <= 0.5) then
                sat = true
                break
            end
        end
        if not sat then unsat = unsat + 1 end
    end
    return unsat
end

function NitroSat:estimate_thermo_stats(beta, n_samples)
    local log_terms = ffi_double(n_samples)
    local energies = ffi_double(n_samples)
    local buf = self._sample_buf
    if not buf then buf = ffi_double(self.num_vars); self._sample_buf = buf end

    for i = 1, n_samples do
        self:_fill_sample(buf)
        local e = self:calculate_energy(buf)
        energies[i] = e
        log_terms[i] = -beta * e
    end

    local logZ = baha_utils.log_sum_exp(log_terms, n_samples)
    local meanE = 0.0
    local meanE2 = 0.0
    local meanE3 = 0.0
    for i = 1, n_samples do
        local w = _exp(log_terms[i] - logZ)
        local e = energies[i]
        meanE = meanE + w * e
        meanE2 = meanE2 + w * e * e
        meanE3 = meanE3 + w * e * e * e
    end
    local varE = meanE2 - meanE * meanE
    if varE < 0 then varE = 0 end

    -- Third central moment proxy for asymmetry near branch splitting.
    local m3 = meanE3 - 3.0 * meanE * meanE2 + 2.0 * meanE * meanE * meanE
    local skew = 0.0
    if varE > 1e-12 then
        skew = m3 / math.pow(varE, 1.5)
    end

    return {
        logZ = logZ,
        meanE = meanE,
        varE = varE,
        skew = skew
    }
end

function NitroSat:estimate_log_Z_sampled(beta, n_samples)
    return self:estimate_thermo_stats(beta, n_samples).logZ
end

function NitroSat:score_branch(beta, n_samples)
    if beta <= 0 then return -math.huge end
    local total_score = 0.0
    local best_seen = math.huge
    local buf = self._sample_buf
    if not buf then buf = ffi_double(self.num_vars); self._sample_buf = buf end
    for i = 1, n_samples do
        self:_fill_sample(buf)
        local energy = self:calculate_energy(buf)
        total_score = total_score + _exp(-beta * energy)
        if energy < best_seen then best_seen = energy end
    end
    return total_score / n_samples + 100.0 / (best_seen + 1.0)
end

function NitroSat:jump_to_branch(beta, n_samples)
    local best_x = nil
    local best_energy = math.huge
    local buf = self._sample_buf
    if not buf then buf = ffi_double(self.num_vars); self._sample_buf = buf end
    for i = 1, n_samples do
        self:_fill_sample(buf)
        local energy = self:calculate_energy(buf)
        if energy < best_energy then
            best_energy = energy
            -- Copy buf to a new best (only when improving)
            best_x = ffi_double(self.num_vars)
            for j = 1, self.num_vars do best_x[j] = buf[j] end
        end
    end
    if best_x then
        local dec = self.decimated
        for i = 1, self.num_vars do
            if dec[i] == 0 then self.x[i] = best_x[i] end
        end
    end
    return best_energy
end

-- Decimation (Soft Lock)
function NitroSat:decimate()
    local locked_count = 0
    for i = 1, self.num_vars do
        if self.decimated[i] == 0 then
            local confidence = math.abs(self.x[i] - 0.5)
            if confidence > self.dec_threshold then
                self.x[i] = self.x[i] > 0.5 and 1.0 or 0.0
                self.decimated[i] = 1
                locked_count = locked_count + 1
            end
        end
    end
    return locked_count
end

-- Local Search Repair (WalkSAT) - Ported from Pro
function NitroSat:local_search_repair(max_flips)
    local clause_sat_counts = ffi_int(self.num_clauses)
    local unsat_list, unsat_pos = {}, {}

    for c = 0, self.num_clauses - 1 do
        local count = 0
        local start, stop = self.clauses_offsets[c], self.clauses_offsets[c+1]
        for i = start, stop - 1 do
            local lit = self.clauses_flat[i]
            local var = math.abs(lit)
            if (lit > 0 and self.x[var] > 0.5) or (lit < 0 and self.x[var] <= 0.5) then
                count = count + 1
            end
        end
        clause_sat_counts[c] = count
        if count == 0 then
            table.insert(unsat_list, c)
            unsat_pos[c] = #unsat_list
        end
    end

    if #unsat_list == 0 then return true end

    for step = 1, max_flips do
        if #unsat_list == 0 then return true end
        local clause_idx = unsat_list[self:_rng_int(1, #unsat_list)]
        local start, stop = self.clauses_offsets[clause_idx], self.clauses_offsets[clause_idx+1]

        local flip_var = -1
        if self:_rng() < 0.3 then 
             flip_var = math.abs(self.clauses_flat[self:_rng_int(start, stop - 1)])
        else
            local min_breaks, candidates = 999999, {}
            for i = start, stop - 1 do
                local var = math.abs(self.clauses_flat[i])
                local breaks = 0
                local v_start, v_stop = self.v2c_ptr[var], self.v2c_ptr[var+1]
                for k_idx = v_start, v_stop - 1 do
                    local k = self.v2c_data[k_idx]
                    if clause_sat_counts[k] == 1 then
                        local k_start, k_stop = self.clauses_offsets[k], self.clauses_offsets[k+1]
                        for j = k_start, k_stop - 1 do
                            local k_lit = self.clauses_flat[j]
                            if math.abs(k_lit) == var then
                                if (k_lit > 0 and self.x[var] > 0.5) or (k_lit < 0 and self.x[var] <= 0.5) then
                                    breaks = breaks + 1
                                end; break
                            end
                        end
                    end
                end
                if breaks < min_breaks then min_breaks = breaks; candidates = {var}
                elseif breaks == min_breaks then table.insert(candidates, var) end
            end
            if #candidates > 0 then flip_var = candidates[self:_rng_int(1, #candidates)] end
        end

        if flip_var ~= -1 then
            local new_val = (self.x[flip_var] > 0.5) and 0.0 or 1.0
            self.x[flip_var] = new_val
            local v_start, v_stop = self.v2c_ptr[flip_var], self.v2c_ptr[flip_var+1]
            for k_idx = v_start, v_stop - 1 do
                local k = self.v2c_data[k_idx]
                local old_count = clause_sat_counts[k]
                local delta = 0
                local k_start, k_stop = self.clauses_offsets[k], self.clauses_offsets[k+1]
                for j = k_start, k_stop - 1 do
                    local lit = self.clauses_flat[j]
                    if math.abs(lit) == flip_var then
                        delta = ((lit > 0 and new_val > 0.5) or (lit < 0 and new_val <= 0.5)) and 1 or -1
                        break
                    end
                end
                local new_count = old_count + delta
                clause_sat_counts[k] = new_count
                if old_count == 0 and new_count > 0 then
                    local pos = unsat_pos[k]
                    if pos then
                        local last_c = unsat_list[#unsat_list]
                        unsat_list[pos] = last_c; unsat_pos[last_c] = pos
                        table.remove(unsat_list); unsat_pos[k] = nil
                    end
                elseif old_count > 0 and new_count == 0 then
                    if not unsat_pos[k] then
                        table.insert(unsat_list, k)
                        unsat_pos[k] = #unsat_list
                    end
                end
            end
        end
    end
    return #unsat_list == 0
end

-- Standard WalkSAT (Refinement - Minimal version)
function NitroSat:refine(max_flips)
    return self:local_search_repair(max_flips)
end

-- ============================================================================
-- BAHA-WalkSAT helpers
-- ============================================================================

-- Compute the break/make counts for flipping a variable.
-- breaks = clauses that go from satisfied to unsatisfied
-- makes  = clauses that go from unsatisfied to satisfied
function NitroSat:compute_flip_delta(flip_var, clause_sat_counts)
    local breaks, makes = 0, 0
    local cur_true = self.x[flip_var] > 0.5
    local v_start, v_stop = self.v2c_ptr[flip_var], self.v2c_ptr[flip_var + 1]
    for k_idx = v_start, v_stop - 1 do
        local c = self.v2c_data[k_idx]
        local cnt = clause_sat_counts[c]
        -- Determine whether this variable currently satisfies clause c
        local c_start, c_stop = self.clauses_offsets[c], self.clauses_offsets[c + 1]
        local var_satisfies = false
        for j = c_start, c_stop - 1 do
            local lit = self.clauses_flat[j]
            if math.abs(lit) == flip_var then
                var_satisfies = (lit > 0 and cur_true) or (lit < 0 and not cur_true)
                break
            end
        end
        if var_satisfies then
            -- Flipping removes this literal's satisfaction
            if cnt == 1 then breaks = breaks + 1 end
        else
            -- Flipping adds satisfaction to this clause
            if cnt == 0 then makes = makes + 1 end
        end
    end
    return breaks, makes
end

-- Execute a flip and incrementally update clause_sat_counts, unsat_list, unsat_pos.
function NitroSat:execute_flip(flip_var, clause_sat_counts, unsat_list, unsat_pos)
    local new_val = (self.x[flip_var] > 0.5) and 0.0 or 1.0
    self.x[flip_var] = new_val
    local v_start, v_stop = self.v2c_ptr[flip_var], self.v2c_ptr[flip_var + 1]
    for k_idx = v_start, v_stop - 1 do
        local k = self.v2c_data[k_idx]
        local old_count = clause_sat_counts[k]
        local delta = 0
        local k_start, k_stop = self.clauses_offsets[k], self.clauses_offsets[k + 1]
        for j = k_start, k_stop - 1 do
            local lit = self.clauses_flat[j]
            if math.abs(lit) == flip_var then
                delta = ((lit > 0 and new_val > 0.5) or (lit < 0 and new_val <= 0.5)) and 1 or -1
                break
            end
        end
        local new_count = old_count + delta
        clause_sat_counts[k] = new_count
        if old_count == 0 and new_count > 0 then
            local pos = unsat_pos[k]
            if pos then
                local last_c = unsat_list[#unsat_list]
                unsat_list[pos] = last_c; unsat_pos[last_c] = pos
                table.remove(unsat_list); unsat_pos[k] = nil
            end
        elseif old_count > 0 and new_count == 0 then
            if not unsat_pos[k] then
                table.insert(unsat_list, k)
                unsat_pos[k] = #unsat_list
            end
        end
    end
end

-- ============================================================================
-- BAHA-WalkSAT: Phase-transition-aware discrete local search
-- Combines Metropolis-Hastings acceptance with BAHA fracture detection
-- and coordinated multi-flip branch jumping.
-- ============================================================================
function NitroSat:baha_walksat(max_flips, opts)
    opts = opts or {}
    local beta_start     = opts.beta_start     or self.bw_beta_start     or 1.0
    local beta_end       = opts.beta_end       or self.bw_beta_end       or 20.0
    local fracture_freq  = opts.fracture_freq  or self.bw_fracture_freq  or 50
    local threshold      = opts.threshold      or self.bw_threshold      or 1.5
    local n_samples      = opts.n_samples      or self.bw_n_samples      or 20
    local _abs, _exp, _floor, _sqrt = math.abs, math.exp, math.floor, math.sqrt

    -- Initialize clause satisfaction counts and unsatisfied-clause list
    local clause_sat_counts = ffi_int(self.num_clauses)
    local unsat_list, unsat_pos = {}, {}
    for c = 0, self.num_clauses - 1 do
        local count = 0
        local start, stop = self.clauses_offsets[c], self.clauses_offsets[c + 1]
        for i = start, stop - 1 do
            local lit = self.clauses_flat[i]
            local var = _abs(lit)
            if (lit > 0 and self.x[var] > 0.5) or (lit < 0 and self.x[var] <= 0.5) then
                count = count + 1
            end
        end
        clause_sat_counts[c] = count
        if count == 0 then
            table.insert(unsat_list, c)
            unsat_pos[c] = #unsat_list
        end
    end
    if #unsat_list == 0 then return true end

    -- Fresh fracture detector for the discrete phase
    local detector = baha_utils.FractureDetector_new(threshold)
    local bw_beta_c = beta_start
    local bw_var_peak = -math.huge
    local bw_var_n, bw_var_mean, bw_var_M2 = 0, 0.0, 0.0

    for step = 1, max_flips do
        if #unsat_list == 0 then return true end

        -- Anneal beta linearly over the walk
        local beta = beta_start + (beta_end - beta_start) * (step / max_flips)

        -- Pick a random unsatisfied clause
        local clause_idx = unsat_list[self:_rng_int(1, #unsat_list)]
        local c_start, c_stop = self.clauses_offsets[clause_idx], self.clauses_offsets[clause_idx + 1]

        -- Greedy: pick variable in clause with minimum break count
        local best_var, best_breaks, best_makes = -1, 999999, 0
        local candidates = {}
        for i = c_start, c_stop - 1 do
            local var = _abs(self.clauses_flat[i])
            local br, mk = self:compute_flip_delta(var, clause_sat_counts)
            local delta = br - mk
            if delta < (best_breaks - best_makes) then
                best_breaks = br; best_makes = mk
                candidates = {var}
            elseif delta == (best_breaks - best_makes) then
                table.insert(candidates, var)
            end
        end
        if #candidates == 0 then goto continue end
        best_var = candidates[self:_rng_int(1, #candidates)]

        -- Recompute delta for chosen variable
        local br, mk = self:compute_flip_delta(best_var, clause_sat_counts)
        local delta_E = br - mk

        -- Metropolis-Hastings acceptance
        if delta_E <= 0 then
            self:execute_flip(best_var, clause_sat_counts, unsat_list, unsat_pos)
        else
            if self:_rng() < _exp(-beta * delta_E) then
                self:execute_flip(best_var, clause_sat_counts, unsat_list, unsat_pos)
            end
        end

        -- BAHA fracture detection (periodic)
        if step % fracture_freq == 0 then
            local thermo = self:estimate_thermo_stats(beta, n_samples)
            detector:record(beta, thermo.logZ)

            -- Track critical beta from heat-capacity proxy peak.
            if thermo.varE > bw_var_peak then
                bw_var_peak = thermo.varE
                bw_beta_c = beta
            end

            -- Rolling variance baseline for z-score spike detection.
            bw_var_n = bw_var_n + 1
            local delta_v = thermo.varE - bw_var_mean
            bw_var_mean = bw_var_mean + delta_v / bw_var_n
            bw_var_M2 = bw_var_M2 + delta_v * (thermo.varE - bw_var_mean)
            local var_std = (bw_var_n > 1) and math.sqrt(math.max(0.0, bw_var_M2 / (bw_var_n - 1))) or 0.0
            local var_z = (bw_var_n > 3 and var_std > 1e-12) and ((thermo.varE - bw_var_mean) / var_std) or 0.0

            local slope_fracture = detector:is_fracture()
            local variance_spike = var_z > self.bw_var_spike_z
            if slope_fracture and variance_spike then
                emit("baha_walksat_fracture", {
                    step = step,
                    beta = beta,
                    beta_c = bw_beta_c,
                    unsat = #unsat_list,
                    varE = thermo.varE,
                    var_z = var_z,
                    meanE = thermo.meanE
                })

                -- Rank variables by instability: count how many unsat clauses each appears in
                local var_unsat_score = {}
                for _, c in ipairs(unsat_list) do
                    local cs, ce = self.clauses_offsets[c], self.clauses_offsets[c + 1]
                    for j = cs, ce - 1 do
                        local v = _abs(self.clauses_flat[j])
                        var_unsat_score[v] = (var_unsat_score[v] or 0) + 1
                    end
                end

                -- Collect and sort by descending instability
                local ranked = {}
                for v, sc in pairs(var_unsat_score) do
                    table.insert(ranked, {var = v, score = sc})
                end
                table.sort(ranked, function(a, b) return a.score > b.score end)

                -- Select top-K unstable variables with variance-scaled jump size.
                local base_K = _floor(_sqrt(#unsat_list) + 0.5)
                local jump_scale = 1.0 + self.bw_var_jump_gain * math.max(0.0, var_z)
                local K = _floor(base_K * jump_scale + 0.5)
                if K < 2 then K = 2 end
                local max_K = _floor(self.num_vars / 4)
                if max_K < 2 then max_K = 2 end
                if K > max_K then K = max_K end
                if K > #ranked then K = #ranked end

                -- Enumerate branches around estimated critical beta.
                local branches = baha_utils.enumerate_branches(bw_beta_c, beta)
                local branch_beta = beta
                if #branches > 0 then
                    local best_score = -math.huge
                    for _, b in ipairs(branches) do
                        local s = self:score_branch(b.beta, n_samples)
                        if s > best_score then
                            best_score = s
                            branch_beta = b.beta
                        end
                    end
                end

                -- Save current state for potential revert
                local saved = ffi_double(self.num_vars)
                for i = 1, self.num_vars do saved[i] = self.x[i] end
                local energy_before = #unsat_list

                -- Coordinated multi-flip: flip top-K unstable variables
                local flipped_vars = {}
                for i = 1, K do
                    local fv = ranked[i].var
                    self:execute_flip(fv, clause_sat_counts, unsat_list, unsat_pos)
                    table.insert(flipped_vars, fv)
                end

                local energy_after = #unsat_list
                local accept_jump = false
                if energy_after <= energy_before then
                    accept_jump = true
                else
                    -- M-H acceptance at branch beta
                    local d = energy_after - energy_before
                    if self:_rng() < _exp(-branch_beta * d) then
                        accept_jump = true
                    end
                end

                if not accept_jump then
                    -- Revert: restore saved state and rebuild unsat structures
                    for i = 1, self.num_vars do self.x[i] = saved[i] end
                    unsat_list, unsat_pos = {}, {}
                    for c = 0, self.num_clauses - 1 do
                        local count = 0
                        local s, e = self.clauses_offsets[c], self.clauses_offsets[c + 1]
                        for j = s, e - 1 do
                            local lit = self.clauses_flat[j]
                            local var = _abs(lit)
                            if (lit > 0 and self.x[var] > 0.5) or (lit < 0 and self.x[var] <= 0.5) then
                                count = count + 1
                            end
                        end
                        clause_sat_counts[c] = count
                        if count == 0 then
                            table.insert(unsat_list, c)
                            unsat_pos[c] = #unsat_list
                        end
                    end
                end

                detector:clear()
            end
        end

        ::continue::
    end
    return #unsat_list == 0
end

-- ============================================================================
-- PHASE 2: TOPOLOGICAL REPAIR - Focus on cycle-breaking variables
-- ============================================================================
function NitroSat:identify_cycle_variables(topology)
    local cycle_vars = {}
    if not topology or not topology.beta1 or topology.beta1 == 0 then
        return cycle_vars
    end
    local adj = {}
    local degrees = ffi_double(self.num_vars + 1)
    for c = 0, self.num_clauses - 1 do
        local start = self.clauses_offsets[c]
        local stop = self.clauses_offsets[c+1]
        local is_sat = false
        for i = start, stop - 1 do
            local lit = self.clauses_flat[i]
            local var = math.abs(lit)
            if (lit > 0 and self.x[var] > 0.5) or (lit < 0 and self.x[var] <= 0.5) then
                is_sat = true; break
            end
        end
        if not is_sat then
            for i = start, stop - 1 do
                local vi = math.abs(self.clauses_flat[i])
                for j = i + 1, stop - 1 do
                    local vj = math.abs(self.clauses_flat[j])
                    adj[vi] = adj[vi] or {}; adj[vj] = adj[vj] or {}
                    local w = (adj[vi][vj] or 0) + 1
                    adj[vi][vj] = w; adj[vj][vi] = w
                    degrees[vi] = degrees[vi] + 1; degrees[vj] = degrees[vj] + 1
                end
            end
        end
    end
    local threshold = topology.complexity_score or 0.3
    for i = 1, self.num_vars do
        if degrees[i] > 0 and self.decimated[i] == 0 then
            local neighbors = adj[i] or {}
            local neighbor_count, total_weight = 0, 0
            for j, w in pairs(neighbors) do
                neighbor_count = neighbor_count + 1
                total_weight = total_weight + w / math.max(degrees[j], 1)
            end
            if neighbor_count > 0 and (total_weight / neighbor_count) >= threshold then
                table.insert(cycle_vars, i)
            end
        end
    end
    return cycle_vars
end

function NitroSat:compute_clause_force(var)
    local force = 0
    local start_ptr = self.v2c_ptr[var]
    local stop_ptr = self.v2c_ptr[var+1]
    for k_idx = start_ptr, stop_ptr - 1 do
        local c = self.v2c_data[k_idx]
        local satisfied = false
        local start = self.clauses_offsets[c]
        local stop = self.clauses_offsets[c+1]
        for i = start, stop - 1 do
            local lit = self.clauses_flat[i]
            local v = math.abs(lit)
            if (lit > 0 and self.x[v] > 0.5) or (lit < 0 and self.x[v] <= 0.5) then
                satisfied = true; break
            end
        end
        if not satisfied then
            for i = start, stop - 1 do
                local lit = self.clauses_flat[i]
                if math.abs(lit) == var then
                    force = force + (lit > 0 and 1 or -1) * (self.clause_weights[c] or 1.0)
                    break
                end
            end
        end
    end
    return force
end

function NitroSat:topological_repair_phase(max_steps)
    emit("phase2_topological_repair_start", {max_steps = max_steps})

    for step = 1, max_steps do
        -- Compute persistent homology
        local topo = self:compute_topology()

        -- Identify "holes" in solution space (beta1 features)
        local hole_vars = self:identify_cycle_variables(topo)

        if #hole_vars == 0 then
            -- No cycles found, switch to local search
            hole_vars = {}
            for i = 1, self.num_vars do
                if self.decimated[i] == 0 then table.insert(hole_vars, i) end
            end
        end

        -- Focused perturbation on cycle-breaking variables
        for _, var in ipairs(hole_vars) do
            if self.decimated[var] == 0 then
                -- Strong flip toward satisfaction
                local force = self:compute_clause_force(var)
                self.x[var] = clamp(self.x[var] + 0.1 * math_sign(force), 0, 1)
            end
        end

        -- Local search on hole region
        local repaired = self:restricted_walksat(hole_vars, self.num_vars * 10)

        local sat = self:check_satisfaction()
        if sat == self.num_clauses then
            emit("phase2_topological_repair_complete", {success = true, steps = step})
            return true
        end
    end

    emit("phase2_topological_repair_complete", {success = false, steps = max_steps})
    return false
end

function NitroSat:restricted_walksat(var_list, max_flips)
    if #var_list == 0 then return true end
    local clause_sat_counts = ffi_int(self.num_clauses)
    local unsat_list, unsat_pos = {}, {}
    for c = 0, self.num_clauses - 1 do
        local count = 0
        local start, stop = self.clauses_offsets[c], self.clauses_offsets[c+1]
        for i = start, stop - 1 do
            local lit = self.clauses_flat[i]
            local var = math.abs(lit)
            if (lit > 0 and self.x[var] > 0.5) or (lit < 0 and self.x[var] <= 0.5) then
                count = count + 1
            end
        end
        clause_sat_counts[c] = count
        if count == 0 then table.insert(unsat_list, c); unsat_pos[c] = #unsat_list end
    end
    if #unsat_list == 0 then return true end
    local var_set = {}
    for _, v in ipairs(var_list) do var_set[v] = true end
    for step = 1, max_flips do
        if #unsat_list == 0 then return true end
        local clause_idx = unsat_list[self:_rng_int(1, #unsat_list)]
        local start, stop = self.clauses_offsets[clause_idx], self.clauses_offsets[clause_idx+1]
        local flip_var = -1
        if self:_rng() < 0.3 then
            local candidates = {}
            for i = start, stop - 1 do
                local v = math.abs(self.clauses_flat[i])
                if var_set[v] then table.insert(candidates, v) end
            end
            if #candidates > 0 then flip_var = candidates[self:_rng_int(1, #candidates)] end
        else
            local min_breaks, candidates = 999999, {}
            for i = start, stop - 1 do
                local var = math.abs(self.clauses_flat[i])
                if not var_set[var] then goto continue end
                local breaks = 0
                local v_start, v_stop = self.v2c_ptr[var], self.v2c_ptr[var+1]
                for k_idx = v_start, v_stop - 1 do
                    local k = self.v2c_data[k_idx]
                    if clause_sat_counts[k] == 1 then
                        local k_start, k_stop = self.clauses_offsets[k], self.clauses_offsets[k+1]
                        for j = k_start, k_stop - 1 do
                            local k_lit = self.clauses_flat[j]
                            if math.abs(k_lit) == var then
                                if (k_lit > 0 and self.x[var] > 0.5) or (k_lit < 0 and self.x[var] <= 0.5) then
                                    breaks = breaks + 1
                                end; break
                            end
                        end
                    end
                end
                if breaks < min_breaks then min_breaks = breaks; candidates = {var}
                elseif breaks == min_breaks then table.insert(candidates, var) end
                ::continue::
            end
            if #candidates > 0 then flip_var = candidates[self:_rng_int(1, #candidates)] end
        end
        if flip_var == -1 then
            for i = start, stop - 1 do
                local v = math.abs(self.clauses_flat[i])
                if var_set[v] then flip_var = v; break end
            end
        end
        if flip_var ~= -1 then
            local new_val = (self.x[flip_var] > 0.5) and 0.0 or 1.0
            self.x[flip_var] = new_val
            local v_start, v_stop = self.v2c_ptr[flip_var], self.v2c_ptr[flip_var+1]
            for k_idx = v_start, v_stop - 1 do
                local k = self.v2c_data[k_idx]
                local old_count = clause_sat_counts[k]
                local delta = 0
                local k_start, k_stop = self.clauses_offsets[k], self.clauses_offsets[k+1]
                for j = k_start, k_stop - 1 do
                    local lit = self.clauses_flat[j]
                    if math.abs(lit) == flip_var then
                        delta = ((lit > 0 and new_val > 0.5) or (lit < 0 and new_val <= 0.5)) and 1 or -1
                        break
                    end
                end
                local new_count = old_count + delta
                clause_sat_counts[k] = new_count
                if old_count == 0 and new_count > 0 then
                    local pos = unsat_pos[k]
                    if pos then
                        local last_c = unsat_list[#unsat_list]
                        unsat_list[pos] = last_c; unsat_pos[last_c] = pos
                        table.remove(unsat_list); unsat_pos[k] = nil
                    end
                elseif old_count > 0 and new_count == 0 then
                    if not unsat_pos[k] then
                        table.insert(unsat_list, k)
                        unsat_pos[k] = #unsat_list
                    end
                end
            end
        end
    end
    return #unsat_list == 0
end

-- ============================================================================
-- PHASE 3: ADELIC SATURATION - Binary search on temperature
-- ============================================================================
function NitroSat:zeta_sweep_aggressive(beta)
    local dF = ffi_double(self.num_vars)
    for i = 1, self.num_vars do dF[i] = 0.0 end
    local unsat_clauses, unsat_count = {}, 0
    for c = 0, self.num_clauses - 1 do
        local is_sat = false
        local start, stop = self.clauses_offsets[c], self.clauses_offsets[c+1]
        for i = start, stop - 1 do
            local lit = self.clauses_flat[i]
            local var = math.abs(lit)
            if (lit > 0 and self.x[var] > 0.5) or (lit < 0 and self.x[var] <= 0.5) then
                is_sat = true; break
            end
        end
        if not is_sat then
            unsat_count = unsat_count + 1
            unsat_clauses[unsat_count] = c
        end
    end
    if unsat_count == 0 then return true, 0, 0 end
    for idx = 1, unsat_count do
        local c = unsat_clauses[idx]
        local w = (self.zeta_weights[c+1] or 0.1) * (self.clause_weights[c] or 1.0) * 2.0
        local start, stop = self.clauses_offsets[c], self.clauses_offsets[c+1]
        for i = start, stop - 1 do
            local lit = self.clauses_flat[i]
            local var = math.abs(lit)
            dF[var] = dF[var] + (lit > 0 and 1 or -1) * w
        end
    end
    local flips = 0
    for i = 1, self.num_vars do
        if math.abs(dF[i]) > 1e-9 then
            if self:_rng() < (1 / (1 + math.exp(-beta * math.abs(dF[i])))) then
                self.x[i] = (dF[i] > 0) and 0.99 or 0.01
                flips = flips + 1
            end
        end
    end
    return false, unsat_count, flips
end

function NitroSat:unit_propagation()
    local changed = false
    for c = 0, self.num_clauses - 1 do
        local start, stop = self.clauses_offsets[c], self.clauses_offsets[c+1]
        local unassigned_lit, satisfied, num_unassigned = nil, false, 0
        for i = start, stop - 1 do
            local lit = self.clauses_flat[i]
            local var = math.abs(lit)
            local val = self.x[var]
            if self.decimated[var] == 1 then
                if (lit > 0 and val > 0.5) or (lit < 0 and val <= 0.5) then satisfied = true; break end
            else
                num_unassigned = num_unassigned + 1
                if num_unassigned == 1 then unassigned_lit = lit
                else unassigned_lit = nil; break end
            end
        end
        if not satisfied and unassigned_lit then
            local var = math.abs(unassigned_lit)
            local force_val = (unassigned_lit > 0) and 1.0 or 0.0
            if math.abs(self.x[var] - force_val) > 0.5 then
                self.x[var] = force_val; self.decimated[var] = 1; changed = true
            end
        end
    end
    return changed
end

function NitroSat:adelic_saturation_phase(max_steps)
    emit("phase3_adelic_saturation_start", {max_steps = max_steps})
    local beta_low, beta_high = 0.1, 10.0
    local last_unsat = self.num_clauses
    for step = 1, max_steps do
        local beta = beta_low + (beta_high - beta_low) * (step / max_steps)
        local done, unsat, flips = self:zeta_sweep_aggressive(beta)
        if done then return true end
        if flips == 0 then beta_low = beta
        elseif unsat > last_unsat then beta_high = beta end
        last_unsat = unsat
        if step % 10 == 0 then self:unit_propagation() end
    end
    return false
end

-- ============================================================================
-- CLAUSE-CORE DECOMPOSITION: Identify and blast the unsat core
-- ============================================================================
function NitroSat:identify_and_blast_core()
    emit("core_decomposition_start", {})
    local core_clauses = {}
    for c = 0, self.num_clauses - 1 do
        local satisfied = false
        local start, stop = self.clauses_offsets[c], self.clauses_offsets[c+1]
        for i = start, stop - 1 do
            local lit = self.clauses_flat[i]
            local var = math.abs(lit)
            if (lit > 0 and self.x[var] > 0.5) or (lit < 0 and self.x[var] <= 0.5) then
                satisfied = true; break
            end
        end
        if not satisfied then table.insert(core_clauses, c) end
    end
    if #core_clauses == 0 then return true end
    emit("core_identified", {core_size = #core_clauses})
    local original_weights = {}
    for _, c in ipairs(core_clauses) do
        original_weights[c] = self.clause_weights[c]
        self.clause_weights[c] = self.clause_weights[c] * 100
    end
    local core_vars = {}, {}
    local var_in_core = {}
    for _, c in ipairs(core_clauses) do
        local start, stop = self.clauses_offsets[c], self.clauses_offsets[c+1]
        for i = start, stop - 1 do
            local v = math.abs(self.clauses_flat[i])
            if not var_in_core[v] then
                var_in_core[v] = true; table.insert(core_vars, v)
            end
        end
    end
    emit("core_blast_start", {core_vars = #core_vars})
    local success = self:restricted_walksat(core_vars, self.num_vars * 20)
    for c, w in pairs(original_weights) do self.clause_weights[c] = w end
    return success
end

-- Solver Loop
function NitroSat:solve(opts)
    opts = opts or {}
    local max_steps = opts.max_steps or 3000
    local lr = opts.lr or 0.2 -- Pro uses 0.2 base
    local use_annealing = opts.use_annealing ~= false
    local time_scale = 30.0 / max_steps
    local verbose = opts.verbose or false
    local start_time = os.clock()
    
    emit("solve_start", {max_steps = max_steps, trace_id = self.trace_id})
    
    local final_step = max_steps
    local success = false
    local best_sat_count = 0
    local best_x = ffi_double(self.num_vars)
    
    local stagnant_count = 0
    local last_sat_rate = 0.0
    local baha_beta_c = 1.0
    local baha_var_peak = -math.huge
    local baha_var_n, baha_var_mean, baha_var_M2 = 0, 0.0, 0.0
    
    for step = 1, max_steps do
        final_step = step
        self.current_step = step  -- Track for persistence

        local grads, unsat_count = self:compute_gradients()
        local current_sat = self.num_clauses - unsat_count
        
        -- Track best
        if current_sat > best_sat_count then
            best_sat_count = current_sat
            for i = 1, self.num_vars do best_x[i] = self.x[i] end
        end
        
        if unsat_count == 0 then
            success = true
            break
        end
        
        -- Dynamic Clause Weighting (SAPS-like)
        -- Periodically flush soft constraints to encourage solving them
        if step % 10 == 0 then
            for c = 0, self.num_clauses - 1 do
                -- Skip hard constraints (fixed high priority)
                -- Actually, if we mark ALL as hard (current default), this loop does nothing.
                -- We need to mark some as soft to see effect.
                -- But Pro defaults to soft (1.0) unless specified.
                -- Let's change init to default SOFT for this to work, or just apply to all.
                -- Pro skips hard constraints.
                -- If is_hard_constraint is ALL true (current NitroSat default), this block is useless.
                -- I'll modify init to make them soft by default for standard CNF, 
                -- or just apply to all for now to test "adaptive re-weighting".
                -- Pro: "if not self.is_hard_constraint[c] then"
                -- I'll relax the check for now since we don't have IsHard input.
                
                -- Check satisfaction
                 local satisfied = false
                 local start, stop = self.clauses_offsets[c], self.clauses_offsets[c+1]
                 for i = start, stop - 1 do
                     local lit = self.clauses_flat[i]
                     local var = math.abs(lit)
                     local val = self.x[var]
                     if (lit > 0 and val > 0.5) or (lit < 0 and val <= 0.5) then
                         satisfied = true; break
                     end
                 end
                 
                 if not satisfied then
                     self.clause_weights[c] = self.clause_weights[c] + 0.01
                 end
                 self.clause_weights[c] = self.clause_weights[c] * 0.9999
            end
        end

        -- Annealing schedule
        local current_lr
        if use_annealing then
            local t = step * time_scale
            current_lr = self:get_annealing_lr(t, lr)
        else
            current_lr = lr
        end

        self:adam_step(grads, current_lr)
        
        -- BAHA: Fracture Detection & Branch Jumping
        -- Use step-based beta proxy
        local beta_proxy = (step / max_steps) * 10.0
        
        -- Use REAL independent sampling for Z estimation (The "Power")
        -- We use a small sample size (e.g. 50) to keep it fast enough for Lua
        -- C++ uses 100, checking every step is too slow.
        -- Let's check every 10 steps to balance speed/accuracy
        if self.use_baha and step % 50 == 0 then
            -- Satisfaction gate: stop jumping when we're close to solving
            local sat_rate = current_sat / self.num_clauses
            if sat_rate < 0.98 then
            local thermo = self:estimate_thermo_stats(beta_proxy, 20)
            self.baha_detector:record(beta_proxy, thermo.logZ)

            if thermo.varE > baha_var_peak then
                baha_var_peak = thermo.varE
                baha_beta_c = beta_proxy
            end
            baha_var_n = baha_var_n + 1
            local delta_v = thermo.varE - baha_var_mean
            baha_var_mean = baha_var_mean + delta_v / baha_var_n
            baha_var_M2 = baha_var_M2 + delta_v * (thermo.varE - baha_var_mean)
            local var_std = (baha_var_n > 1) and math.sqrt(math.max(0.0, baha_var_M2 / (baha_var_n - 1))) or 0.0
            local var_z = (baha_var_n > 3 and var_std > 1e-12) and ((thermo.varE - baha_var_mean) / var_std) or 0.0

            local slope_fracture = self.baha_detector:is_fracture()
            local variance_spike = var_z > self.baha_var_spike_z
            if slope_fracture and variance_spike then
                -- [NEW] Get topology for smarter branch jumping
                local topology = self:compute_topology()
                local topology_info = {}
                if topology then
                    topology_info = {
                        beta0 = topology.beta0,
                        beta1 = topology.beta1,
                        complexity_score = topology.complexity_score
                    }
                    if verbose then print(string.format("  Topology: beta0=%d beta1=%d complexity=%.3f",
                        topology.beta0, topology.beta1, topology.complexity_score)) end
                end

                if verbose then print(string.format("⚡ NITROSAT FRACTURE DETECTED at step %d (beta=%.2f) ⚡", step, beta_proxy)) end
                emit("baha_fracture", {
                    step = step,
                    beta = beta_proxy,
                    beta_c = baha_beta_c,
                    trace_id = self.trace_id,
                    topology = topology_info,
                    varE = thermo.varE,
                    var_z = var_z,
                    meanE = thermo.meanE
                })
                
                -- Enumerate branches using Lambert W
                local branches = baha_utils.enumerate_branches(baha_beta_c, beta_proxy, 5)
                
                if #branches > 0 then
                    -- Score branches
                    local best_branch = nil
                    local best_score = -math.huge
                    
                    for _, b in ipairs(branches) do
                        local s = self:score_branch(b.beta, 20) -- Branch scoring sample count reduced to 20
                        if s > best_score then -- Corrected variable name from 'score' to 's'
                            best_score = s
                            best_branch = b
                        end
                    end
                    
                    if best_branch then
                        if verbose then 
                            print(string.format("  -> Jumping to branch k=%d, beta=%.2f (Score: %.2f)", 
                                best_branch.k, best_branch.beta, best_score)) 
                        end
                        
                        -- Execute Jump: Reset state to best sample in new branch
                        self:jump_to_branch(best_branch.beta, 20)
                        
                        -- Reset Optimizer (Clear momentum)
                        self.optimizer = optimizer.create_optimizer('harmonic_adam', self.x, {
                            lr = 0.02,
                            beta1 = 0.9,
                            beta2 = 0.999
                        }, self.num_vars)
                        
                        self.baha_detector:clear()
                        -- Recalculate satisfactions for next step
                        local _, u_count = self:compute_gradients()
                        local c_sat = self.num_clauses - u_count
                        if c_sat > best_sat_count then
                            best_sat_count = c_sat
                            for i=1, self.num_vars do best_x[i] = self.x[i] end
                        end
                    end
                end
            end
            end -- sat_rate gate
        end

        -- [NEW] Persistent Homology Topology Snapshots
        if self.use_topology and step % self.topology_freq == 0 then
            local topology = self:compute_topology()
            if topology then
                emit("topology_snapshot", {
                    step = step,
                    beta0 = topology.beta0,
                    beta1 = topology.beta1,
                    active_vars = topology.active_vars,
                    edge_count = topology.edge_count,
                    complexity_score = topology.complexity_score,
                    avg_confidence = topology.avg_confidence,
                    unsat_count = unsat_count,
                    sat_rate = current_sat / self.num_clauses,
                    trace_id = self.trace_id
                })

                -- Detect topological phase transitions
                if self.last_topology then
                    local beta0_change = math.abs(topology.beta0 - self.last_topology.beta0)
                    local beta1_change = math.abs(topology.beta1 - self.last_topology.beta1)

                    if beta0_change > 5 or beta1_change > 3 then
                        emit("topology_phase_transition", {
                            step = step,
                            beta0_change = beta0_change,
                            beta1_change = beta1_change,
                            trace_id = self.trace_id
                        })
                    end
                end
            end
        end

        if self.dec_freq and step > 0 and step % self.dec_freq == 0 then
            local locked = self:decimate()
            if verbose and locked > 0 then print(string.format("Decimated (locked) %d variables", locked)) end
        end

        -- Zeta zero guided perturbation (The Final Piece)
        -- Nuclear Zeta Mode: Use enhanced perturbation when enabled
        if self.zeta_guidance and step % 10 == 0 then
            local current_topo = nil
            if self.use_topology and self.last_topology then
                current_topo = self.last_topology
            end

            for i = 1, self.num_vars do
                if self.decimated[i] == 0 then
                    local zeta_perturb
                    if self.nuclear_zeta then
                        -- Nuclear Zeta: Adelic force saturation with critical line resonance
                        zeta_perturb = zeta_zero_perturbation_nuclear(i, step, max_steps, self.primes, current_topo, self.num_vars, self.prime_count) * 0.001
                    else
                        -- Standard: Gentle prime perturbation
                        zeta_perturb = zeta_zero_perturbation(i, step, max_steps, self.primes, self.prime_count) * 0.001
                    end
                    self.x[i] = clamp(self.x[i] + zeta_perturb, 0.0, 1.0)
                end
            end
        end
        
        -- Inline DCW: AdaBoost-style clause gravity
        -- Increase weight for persistently unsatisfied clauses (like AdaFactor)
        if step % 10 == 0 then
            for c = 0, self.num_clauses - 1 do
                local satisfied = false
                local start, stop = self.clauses_offsets[c], self.clauses_offsets[c+1]
                for i = start, stop - 1 do
                    local lit = self.clauses_flat[i]
                    local var = math.abs(lit)
                    local val = self.x[var]
                    if (lit > 0 and val > 0.5) or (lit < 0 and val <= 0.5) then
                        satisfied = true
                        break
                    end
                end
                if not satisfied then
                    -- Gravity increases for unsatisfied clauses
                    self.clause_weights[c] = self.clause_weights[c] * 2.5
                end
            end
            -- Decay all weights slightly to prevent explosion
            for i = 0, self.num_clauses - 1 do
                self.clause_weights[i] = self.clause_weights[i] * 0.99
            end
        end

        -- Stagnation Check
        if step % 50 == 0 then
             local sat_now = self:check_satisfaction()
             local sat_rate = sat_now / self.num_clauses

             if sat_rate - last_sat_rate < 0.001 then
                 stagnant_count = stagnant_count + 1
             else
                 stagnant_count = 0
             end
             last_sat_rate = sat_rate
        end

        if step % 100 == 0 and verbose then
           local sat = self:check_satisfaction()
           print(string.format("Step %d: %d/%d (%.2f%%)", step, sat, self.num_clauses, sat/self.num_clauses*100))
        end
    end
    
    -- Restore best state
    for i = 1, self.num_vars do self.x[i] = best_x[i] end
    if best_sat_count == self.num_clauses then success = true end
    
    -- Final check (restore best if needed, but for now just return stats)
    -- If we ended on a worse state, strictly speaking for "solve" we might want the current state,
    -- but for MaxSAT we return the best stats.
    
    -- ========================================================================
    -- Zeta Sweep: High-precision finisher (Adelic Repair Phase)
    -- ========================================================================
    if best_sat_count < self.num_clauses then
        emit("zeta_sweep_start", {best_sat = best_sat_count, total = self.num_clauses})
        for i = 1, 50 do
            local done, unsat, flips = self:zeta_sweep(1.0 + i*0.1)
            if done then
                best_sat_count = self.num_clauses
                success = true
                break
            end
            if unsat < (self.num_clauses - best_sat_count) then
                best_sat_count = self.num_clauses - unsat
            end
            if flips == 0 then break end
        end
    end

    -- ========================================================================
    -- THREE-PHASE FINISHER PROTOCOL
    -- Phase 1 already completed above (gradient convergence)
    -- Phase 2: Topological Repair (when at 95%+ satisfaction)
    -- Phase 3: Adelic Saturation (when at 98%+ satisfaction)
    -- ========================================================================

    local sat_rate = best_sat_count / self.num_clauses

    -- PHASE 2: Topological Repair
    -- Activates when standard gradients plateau near completion
    if best_sat_count < self.num_clauses and sat_rate >= self.phase2_start then
        emit("three_phase_finisher_start", {phase = 2, sat_rate = sat_rate})

        -- Restore best state
        for i = 1, self.num_vars do self.x[i] = best_x[i] end

        local phase2_success = self:topological_repair_phase(1500)
        local phase2_sat = self:check_satisfaction()

        if phase2_sat > best_sat_count then
            best_sat_count = phase2_sat
            for i = 1, self.num_vars do best_x[i] = self.x[i] end
        end

        if best_sat_count == self.num_clauses then
            success = true
        end
        emit("phase2_complete", {sat = best_sat_count, success = phase2_success})
    end

    -- PHASE 3: Adelic Saturation
    -- Activates when Phase 2 plateaus - Nuclear option with binary search
    if best_sat_count < self.num_clauses and sat_rate >= self.phase3_start then
        emit("three_phase_finisher_phase3_start", {sat_rate = best_sat_count / self.num_clauses})

        -- Restore best state
        for i = 1, self.num_vars do self.x[i] = best_x[i] end

        local phase3_success = self:adelic_saturation_phase(2000)
        local phase3_sat = self:check_satisfaction()

        if phase3_sat > best_sat_count then
            best_sat_count = phase3_sat
            for i = 1, self.num_vars do best_x[i] = self.x[i] end
        end

        if best_sat_count == self.num_clauses then
            success = true
        end
        emit("phase3_complete", {sat = best_sat_count, success = phase3_success})
    end

    -- ========================================================================
    -- CLAUSE-CORE DECOMPOSITION: Blast the unsat core
    -- ========================================================================
    if best_sat_count < self.num_clauses then
        -- Restore best state
        for i = 1, self.num_vars do self.x[i] = best_x[i] end

        local core_success = self:identify_and_blast_core()
        local core_sat = self:check_satisfaction()

        if core_sat > best_sat_count then
            best_sat_count = core_sat
            for i = 1, self.num_vars do best_x[i] = self.x[i] end
        end

        if best_sat_count == self.num_clauses then
            success = true
        end
    end

    -- ========================================================================
    -- BAHA-WalkSAT finisher (phase-transition-aware discrete local search)
    -- ========================================================================
    if best_sat_count < self.num_clauses and self.num_vars > 0 then
        emit("baha_walksat_start", {best_sat = best_sat_count, total = self.num_clauses})
        -- Restore best state and snap to discrete
        for i = 1, self.num_vars do
            self.x[i] = best_x[i] > 0.5 and 1.0 or 0.0
        end

        local repaired = self:baha_walksat(self.num_vars * 100)
        local final_sat = self:check_satisfaction()

        if final_sat > best_sat_count then
            best_sat_count = final_sat
            for i = 1, self.num_vars do best_x[i] = self.x[i] end
        end
        if best_sat_count == self.num_clauses then
            success = true
        end
        emit("baha_walksat_complete", {satisfied = best_sat_count, total = self.num_clauses})
    end

    -- [NEW] Final topology summary
    local topology_summary = nil
    if self.use_topology and self.persistence_tracker then
        local tracker = self.persistence_tracker
        topology_summary = {
            initial_beta0 = tracker.beta0_history[1] or 0,
            final_beta0 = tracker.beta0_history[#tracker.beta0_history] or 0,
            initial_beta1 = tracker.beta1_history[1] or 0,
            final_beta1 = tracker.beta1_history[#tracker.beta1_history] or 0,
            persistence_events = #tracker.persistence_pairs,
            complexity_trend = (#tracker.complexity_history > 0) and
                (tracker.complexity_history[#tracker.complexity_history] - tracker.complexity_history[1]) or 0
        }
    end

    emit("solve_complete", {
        success = success,
        duration_ms = (os.clock() - start_time) * 1000,
        final_step = final_step,
        best_sat = best_sat_count,
        trace_id = self.trace_id,
        topology = topology_summary
    })
    
    return success, final_step, best_sat_count
end

-- ==============================================================================
-- Dynamic Clause Weighting (DCW): 5-pass SAPS-style refinement
-- ==============================================================================
function NitroSat:solve_dcw(opts)
    opts = opts or {}
    local num_passes = opts.dcw_passes or 5

    -- Skip DCW for edge cases
    if self.num_vars == 0 or self.num_clauses == 0 then
        return self:solve(opts)
    end

    emit("dcw_start", {num_passes = num_passes, trace_id = self.trace_id})

    local best_overall_sat = 0
    local best_overall_x = ffi_double(self.num_vars)
    local overall_success = false

    -- Initialize clause weights to 1.0
    for i = 1, self.num_clauses do
        self.clause_weights[i] = 1.0
    end

    for pass = 1, num_passes do
        emit("dcw_pass_start", {pass = pass, num_passes = num_passes})

        -- Run main solver
        local success, steps, sat_count = self:solve(opts)

        if sat_count > best_overall_sat then
            best_overall_sat = sat_count
            for i = 1, self.num_vars do
                best_overall_x[i] = self.x[i]
            end
        end

        if sat_count == self.num_clauses then
            overall_success = true
            break
        end

        -- Increase weights for unsatisfied clauses
        if pass < num_passes then
            emit("dcw_reweight_start", {pass = pass})
            local unsat_count = 0
            for c = 0, self.num_clauses - 1 do
                local satisfied = false
                local start, stop = self.clauses_offsets[c], self.clauses_offsets[c+1]
                for i = start, stop - 1 do
                    local lit = self.clauses_flat[i]
                    local var = math.abs(lit)
                    if (lit > 0 and self.x[var] > 0.5) or (lit < 0 and self.x[var] <= 0.5) then
                        satisfied = true; break
                    end
                end
                if not satisfied then
                    self.clause_weights[c] = self.clause_weights[c] * 1.5
                    unsat_count = unsat_count + 1
                end
            end
            emit("dcw_reweight_done", {pass = pass, unsat_increased = unsat_count})

            -- Decay all weights slightly to prevent explosion
            for i = 1, self.num_clauses do
                self.clause_weights[i] = self.clause_weights[i] * 0.95
                if self.clause_weights[i] < 0.1 then self.clause_weights[i] = 0.1 end
            end

            -- Reset optimizer for fresh start
            self.optimizer = optimizer.create_optimizer('nadam', self.x, {
                lr = 0.002,
                beta1 = 0.9,
                beta2 = 0.999,
                resonance_amplitude = 0.02
            }, self.num_vars)
        end
    end

    -- Restore best state
    for i = 1, self.num_vars do self.x[i] = best_overall_x[i] end

    emit("dcw_complete", {
        success = overall_success,
        best_sat = best_overall_sat,
        trace_id = self.trace_id
    })

    return overall_success, num_passes, best_overall_sat
end

-- ==============================================================================
-- SHUNYABAR SOLVER MODE
-- ==============================================================================
function NitroSat:solve_shunyabar(opts)
    opts = opts or {}
    local max_steps = opts.max_steps or 3500
    local verbose = opts.verbose ~= false
    local start_time = os.clock()
    
    emit("shunyabar_start", {
        max_steps = max_steps,
        num_vars = self.num_vars,
        num_clauses = self.num_clauses,
        trace_id = self.trace_id
    })

    -- RG Sweep Parameter Range: s in [0.5, 1.0]
    local s_initial = 0.5
    local s_final = 1.0
    
    local best_sat_count = 0
    local final_step = 0
    local success = false
    
    for step = 1, max_steps do
        -- Quasi-static RG Sweep: Progressively move toward the arithmetic pole
        local t = step / max_steps
        local s = s_initial + (s_final - s_initial) * t
        
        -- Compute standard gradients (Geometric structure)
        local grads, unsat_count = self:compute_gradients()
        local current_sat = self.num_clauses - unsat_count
        
        -- Compute ShunyaBar Spectral-Arithmetic Force (Arithmetic structure)
        -- Using Laplacian-weighted arithmetic gain
        local shunyabar_forces = compute_shunyabar_force(self.x, s, self.primes, self.degrees, self.num_vars)
        
        -- Update gradients: combine geometric descent with arithmetic instability
        -- The arithmetic force destroys illegal regions by making them unstable.
        for i = 1, self.num_vars do
            grads[i] = grads[i] + shunyabar_forces[i]
        end
        
        -- Step using Adams with Lookahead (KMS state stabilization)
        self:adam_step(grads, opts.lr or 0.002)
        
        -- Symmetry Breaking: Detect 1-RSB transition as s -> 1
        -- As s approaches 1, the arithmetic pole forces clusters to shatter
        if s > 0.9 and step % 10 == 0 then
            local pole_gain = zeta_log_derivative(s, self.primes, 200)
            if pole_gain > 10.0 then -- High gain indicates proximity to singularity
                for i = 1, self.num_vars do
                    if self.decimated[i] == 0 then
                        -- Harmonic perturbation to shatter isomorphic variable clusters
                        local pert = zeta_zero_perturbation_nuclear(i, step, max_steps, self.primes, self.last_topology, self.num_vars, self.prime_count)
                        self.x[i] = clamp(self.x[i] + pert * 0.005, 0.0, 1.0)
                    end
                end
            end
        end
        
        -- Check for 100% satisfaction (ShunyaBar terminates immediately)
        local sat_now = self:check_satisfaction()
        if sat_now > best_sat_count then
            best_sat_count = sat_now
            if sat_now == self.num_clauses then
                success = true
                final_step = step
                break
            end
        end
        
        if step % 100 == 0 and verbose then
            print(string.format("ShunyaBar Step %d (s=%.3f): %d/%d", step, s, sat_now, self.num_clauses))
        end
    end
    
    emit("shunyabar_complete", {
        success = success,
        duration_ms = (os.clock() - start_time) * 1000,
        final_step = final_step,
        best_sat = best_sat_count,
        trace_id = self.trace_id
    })
    
    return success, final_step, best_sat_count
end

return { NitroSat = NitroSat }


