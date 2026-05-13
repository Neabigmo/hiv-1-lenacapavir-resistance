import pandas as pd
from pathlib import Path
from upsetplot import UpSet, from_contents
import matplotlib.pyplot as plt

DATA_DIR = Path(r'H:\2026try\4.20\JMV\data\processed\revision_v2')
df = pd.read_csv(DATA_DIR / 'harmonized_phenotype_data.csv')

resistance_mutations = ['N57H','M66I','Q67H','Q67K','Q67H+K70R','Q67H+N74D',
                        'K70N','K70R','N74D','L56V','M66I+A105T','M66I+T107A',
                        'M66I+N74D','K70N+N74K','M66I+N74D+A105T']
df = df[df['Mutation'].isin(resistance_mutations)]

print("df shape:", df.shape)
print("df Mutation values:", df['Mutation'].tolist())

# Build from_contents dict exactly as the script does
mutation_items = {}
for _, row in df.iterrows():
    mut = row['Mutation']
    item = (row['Subtype'], row['context_tier'])
    if mut not in mutation_items:
        mutation_items[mut] = []
    mutation_items[mut].append(item)

print("\nmutation_items:", dict(mutation_items))

# Test 1: from_contents with tuples
print("\n--- Test 1: from_contents with tuples ---")
try:
    upset_data = from_contents(mutation_items)
    print("OK: from_contents succeeded")
    print("upset_data type:", type(upset_data))
    print("upset_data:\n", upset_data)
except Exception as e:
    print(f"FAILED: {type(e).__name__}: {e}")

# Test 2: Try with string items instead
print("\n--- Test 2: from_contents with string items ---")
mutation_items_str = {}
for _, row in df.iterrows():
    mut = row['Mutation']
    item = f"{row['Subtype']}_{row['context_tier']}"
    if mut not in mutation_items_str:
        mutation_items_str[mut] = []
    mutation_items_str[mut].append(item)

try:
    upset_data2 = from_contents(mutation_items_str)
    print("OK: from_contents with strings succeeded")

    fig, ax = plt.subplots(figsize=(8, 6))
    upset_obj = UpSet(upset_data2, subset_size='count',
                      intersection_plot_elements=6,
                      sort_by='cardinality', sort_category_by='cardinality')
    upset_obj.plot(ax=ax)
    plt.savefig(r'H:\2026try\4.20\JMV\manuscript\figures\test_upset.png', dpi=150, bbox_inches='tight')
    print("Saved test_upset.png successfully")
except Exception as e:
    print(f"FAILED: {type(e).__name__}: {e}")