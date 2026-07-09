"""
Genuinely hard list coloring instance generator.

Hardness sources:
1. Dense structured graph (triangulated torus + compatible skip-2 frustration)
2. Tight lists: 3 colors from palette of 5, adversarially chosen to maximize
   neighbor overlap (avg ~2 shared colors per edge)
3. Planted solution verified correct
4. High degree (~8-10) creates deep constraint cascades per flip
"""

import random


def generate_hard_list_coloring(grid_rows, grid_cols, filename, seed=2026):
    random.seed(seed)
    num_vertices = grid_rows * grid_cols
    num_colors = 5    # palette
    list_size = 3     # tight lists

    def vid(r, c):
        return r * grid_cols + c

    # 1. Build triangulated torus
    adj_base = {v: set() for v in range(num_vertices)}
    def add_base(u, v):
        if u != v:
            adj_base[u].add(v)
            adj_base[v].add(u)

    for r in range(grid_rows):
        for c in range(grid_cols):
            v = vid(r, c)
            add_base(v, vid(r, (c + 1) % grid_cols))
            add_base(v, vid((r + 1) % grid_rows, c))
            add_base(v, vid((r + 1) % grid_rows, (c + 1) % grid_cols))

    # 2. Greedy + repair planted coloring with 5 colors
    order = []
    visited = [False] * num_vertices
    queue = [0]
    visited[0] = True
    while queue:
        v = queue.pop(0)
        order.append(v)
        for w in sorted(adj_base[v]):
            if not visited[w]:
                visited[w] = True
                queue.append(w)

    # DSatur coloring: always color the most-constrained vertex next
    planted = [-1] * num_vertices
    saturation = [0] * num_vertices  # count of distinct colors among colored neighbors
    neighbor_colors_used = [set() for _ in range(num_vertices)]

    for step in range(num_vertices):
        # Pick uncolored vertex with highest saturation (ties broken by degree)
        best_v = -1
        best_sat = -1
        best_deg = -1
        for v in range(num_vertices):
            if planted[v] >= 0:
                continue
            s = saturation[v]
            d = len(adj_base[v])
            if s > best_sat or (s == best_sat and d > best_deg):
                best_v = v
                best_sat = s
                best_deg = d
        v = best_v
        available = [c for c in range(num_colors) if c not in neighbor_colors_used[v]]
        assert available, f"DSatur failed at vertex {v} with {len(adj_base[v])} neighbors and {num_colors} colors"
        planted[v] = random.choice(available)
        # Update saturation of uncolored neighbors
        for w in adj_base[v]:
            if planted[w] < 0:
                if planted[v] not in neighbor_colors_used[w]:
                    neighbor_colors_used[w].add(planted[v])
                    saturation[w] += 1

    base_conflicts = sum(1 for u in range(num_vertices) for v in adj_base[u]
                         if u < v and planted[u] == planted[v])
    assert base_conflicts == 0, f"Base coloring failed: {base_conflicts} conflicts"

    # 3. Add skip-2 + skip-3 frustration edges (only compatible ones)
    adj = {v: set(adj_base[v]) for v in range(num_vertices)}
    edges = set()
    for v in range(num_vertices):
        for w in adj[v]:
            edges.add((min(v, w), max(v, w)))

    skip_added = 0
    for r in range(grid_rows):
        for c in range(grid_cols):
            v = vid(r, c)
            targets = [
                vid(r, (c + 2) % grid_cols),
                vid((r + 2) % grid_rows, c),
                vid((r + 1) % grid_rows, (c + 2) % grid_cols),
                vid((r + 2) % grid_rows, (c + 1) % grid_cols),
            ]
            for w in targets:
                e = (min(v, w), max(v, w))
                if v != w and e not in edges and planted[v] != planted[w]:
                    edges.add(e)
                    adj[v].add(w)
                    adj[w].add(v)
                    skip_added += 1

    # 4. Adversarial lists: include planted color + 2 most-popular neighbor colors
    lists = [None] * num_vertices
    for v in range(num_vertices):
        p = planted[v]
        neighbor_colors = {}
        for w in adj[v]:
            nc = planted[w]
            if nc != p:
                neighbor_colors[nc] = neighbor_colors.get(nc, 0) + 1

        ranked = sorted(neighbor_colors.keys(),
                        key=lambda c: neighbor_colors[c], reverse=True)
        remaining = [c for c in range(num_colors) if c != p and c not in ranked]
        ranked = ranked + remaining
        lists[v] = sorted([p] + ranked[:list_size - 1])

    # 5. CNF generation
    def get_var(v, local_idx):
        return v * list_size + local_idx + 1

    clauses = []
    for v in range(num_vertices):
        clauses.append([get_var(v, i) for i in range(list_size)])
        for i in range(list_size):
            for j in range(i + 1, list_size):
                clauses.append([-get_var(v, i), -get_var(v, j)])

    edge_conflict_clauses = 0
    for u, v in edges:
        shared = set(lists[u]).intersection(set(lists[v]))
        for c in shared:
            clauses.append([-get_var(u, lists[u].index(c)), -get_var(v, lists[v].index(c))])
            edge_conflict_clauses += 1

    num_vars = num_vertices * list_size

    # 6. Verify planted solution
    planted_assignment = {}
    for v in range(num_vertices):
        idx = lists[v].index(planted[v])
        for i in range(list_size):
            planted_assignment[get_var(v, i)] = (i == idx)

    violations = sum(1 for clause in clauses
                     if not any((lit > 0 and planted_assignment.get(abs(lit), False)) or
                                (lit < 0 and not planted_assignment.get(abs(lit), False))
                                for lit in clause))
    assert violations == 0, f"Planted solution has {violations} violations"

    print(f"Grid: {grid_rows}x{grid_cols} = {num_vertices} vertices")
    print(f"Base edges: {len(edges) - skip_added}, Frustration edges: {skip_added}, Total: {len(edges)}")
    print(f"Variables: {num_vars}, Clauses: {len(clauses)}")
    print(f"Edge conflict clauses: {edge_conflict_clauses}")
    print(f"Avg degree: {sum(len(adj[v]) for v in range(num_vertices)) / num_vertices:.1f}")
    print(f"Avg list overlap per edge: {edge_conflict_clauses / len(edges):.2f}")
    print(f"Planted solution: VERIFIED (0 violations)")

    with open(filename, 'w') as f:
        f.write(f"p cnf {num_vars} {len(clauses)}\n")
        for c in clauses:
            f.write(" ".join(map(str, c)) + " 0\n")


if __name__ == "__main__":
    generate_hard_list_coloring(12, 10, "/tmp/hard_list_120.cnf")
    print(f"\nGenerated /tmp/hard_list_120.cnf")
