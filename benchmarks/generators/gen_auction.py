"""
Hard benchmark: Combinatorial Auction Winner Determination as WCNF

This is a classic NP-hard optimization problem used in MaxSAT evaluations
(e.g., spectrum auctions, procurement). It creates extremely dense conflict
structures where bids sharing items form large cliques.

Problem:
- I items for sale
- B bids, each requesting a subset of items with a monetary value
- Hard: each item allocated to at most one winning bid
- Soft: maximize total auction revenue (accept high-value bids)

This is equivalent to Maximum Weight Independent Set on the bid conflict
graph — one of the hardest problems for local search solvers due to
dense clique interactions and non-local value propagation.
"""

import random
import itertools


def generate_auction(num_items, num_bids, max_bundle_size, filename, seed=2026):
    random.seed(seed)
    top = 10000000

    # Generate bids: each bid requests a random subset of items
    bids = []
    for b in range(num_bids):
        bundle_size = random.randint(1, max_bundle_size)
        items = sorted(random.sample(range(num_items), bundle_size))
        # Value correlated with bundle size (synergy) + random premium
        base_value = bundle_size * 10
        synergy = bundle_size * bundle_size  # superadditive synergy
        noise = random.randint(-5, 15)
        value = max(1, base_value + synergy + noise)
        bids.append((items, value))

    # Variable b+1 = bid b is accepted
    num_vars = num_bids
    hard_clauses = []
    soft_clauses = []

    # HARD: each item allocated to at most one bid
    # For each item, collect all bids that contain it
    item_to_bids = {i: [] for i in range(num_items)}
    for b, (items, value) in enumerate(bids):
        for item in items:
            item_to_bids[item].append(b)

    conflict_clauses = 0
    for item in range(num_items):
        bidders = item_to_bids[item]
        # At-most-one: for each pair of bids on this item, they conflict
        for b1, b2 in itertools.combinations(bidders, 2):
            hard_clauses.append([-(b1 + 1), -(b2 + 1)])
            conflict_clauses += 1

    # SOFT: maximize revenue by accepting bids
    for b, (items, value) in enumerate(bids):
        soft_clauses.append((value, [b + 1]))

    num_clauses = len(hard_clauses) + len(soft_clauses)
    total_soft_weight = sum(w for w, _ in soft_clauses)

    # Compute conflict graph density
    adj = {b: set() for b in range(num_bids)}
    for item in range(num_items):
        bidders = item_to_bids[item]
        for b1, b2 in itertools.combinations(bidders, 2):
            adj[b1].add(b2)
            adj[b2].add(b1)
    avg_conflicts = sum(len(adj[b]) for b in range(num_bids)) / num_bids

    print(f"Combinatorial Auction Winner Determination Benchmark")
    print(f"  Items: {num_items}, Bids: {num_bids}, Max bundle: {max_bundle_size}")
    print(f"  Variables: {num_vars}")
    print(f"  Hard clauses (item conflicts): {len(hard_clauses)}")
    print(f"  Soft clauses (bid values): {len(soft_clauses)}")
    print(f"  Total clauses: {num_clauses}")
    print(f"  Total potential revenue: {total_soft_weight}")
    print(f"  Avg bid conflicts: {avg_conflicts:.1f}")
    print(f"  Constraint density: {num_clauses / num_vars:.1f}")

    # Analyze item contention
    max_contention = max(len(item_to_bids[i]) for i in range(num_items))
    avg_contention = sum(len(item_to_bids[i]) for i in range(num_items)) / num_items
    print(f"  Max item contention: {max_contention} bids")
    print(f"  Avg item contention: {avg_contention:.1f} bids")

    with open(filename, 'w') as f:
        f.write(f"p wcnf {num_vars} {num_clauses} {top}\n")
        for c in hard_clauses:
            f.write(f"{top} " + " ".join(map(str, c)) + " 0\n")
        for w, c in soft_clauses:
            f.write(f"{w} " + " ".join(map(str, c)) + " 0\n")


if __name__ == "__main__":
    # 30 items, 200 bids competing, bundles up to size 6
    # This creates dense clique interactions and high contention
    generate_auction(
        num_items=30, num_bids=200, max_bundle_size=6,
        filename="/tmp/auction_200.wcnf"
    )
    print(f"\nGenerated /tmp/auction_200.wcnf")
