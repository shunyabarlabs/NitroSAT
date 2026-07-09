from eval_auction import generate_bids, greedy_baseline, nitrosat_revenue

num_items = 30
num_bids = 200
max_bundle_size = 6
bids = generate_bids(num_items, num_bids, max_bundle_size)

greedy_accepted, greedy_rev = greedy_baseline(bids, num_items)
nitro_accepted, nitro_rev, _ = nitrosat_revenue("/tmp/auction_200.wcnf", "src/c/v3/nitrosatv3", bids)

def analyze(name, accepted_bids):
    items_used = set()
    total_val = 0
    valid = True
    print(f"\n--- {name} ---")
    for b in accepted_bids:
        items, val = bids[b]
        overlap = items_used.intersection(set(items))
        if overlap:
            print(f"  CONFLICT! Bid {b} overlaps on {overlap}")
            valid = False
        items_used.update(items)
        total_val += val
        print(f"  Bid {b}: val {val}, items {items}")
    
    print(f"  Valid: {valid}")
    print(f"  Revenue: {total_val}")
    print(f"  Total items allocated: {len(items_used)} / {num_items}")

analyze("Greedy", greedy_accepted)
analyze("NitroSAT", nitro_accepted)
