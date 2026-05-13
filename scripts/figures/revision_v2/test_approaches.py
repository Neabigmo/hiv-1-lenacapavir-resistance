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

# Approach A: from_contents with deduplicated (unique) items — loses count info
print("=== Approach A: from_contents (unique items) ===")
mutation_items = {}
for _, row in df.iterrows():
    mut = row['Mutation']
    item = (row['Subtype'], row['context_tier'])
    if mut not in mutation_items:
        mutation_items[mut] = set()
    mutation_items[mut].add(item)
mutation_items = {k: list(v) for k, v in mutation_items.items()}
print("mutation_items:", dict(mutation_items))

try:
    upset_data = from_contents(mutation_items)
    print("from_contents OK")
    fig, ax = plt.subplots(figsize=(10, 6))
    upset_obj = UpSet(upset_data, subset_size='count',
                      intersection_plot_elements=6,
                      sort_by='cardinality', sort_category_by='cardinality')
    upset_obj.plot(ax=ax)
    plt.savefig(r'H:\2026try\4.20\JMV\manuscript\figures\test_upset_a.png', dpi=150, bbox_inches='tight')
    print("Saved test_upset_a.png")
    plt.close()
except Exception as e:
    print(f"FAILED: {type(e).__name__}: {e}")

# Approach B: from_contents with unique IDs to preserve duplicates
print("\n=== Approach B: from_contents with unique IDs ===")
mutation_items2 = {}
for idx, (_, row) in enumerate(df.iterrows()):
    mut = row['Mutation']
    item = f"{row['Subtype']}_{row['context_tier']}_{idx}"  # unique ID per observation
    if mut not in mutation_items2:
        mutation_items2[mut] = []
    mutation_items2[mut].append(item)
print("mutation_items2:", dict(mutation_items2))

try:
    upset_data2 = from_contents(mutation_items2)
    print("from_contents OK")
    fig, ax = plt.subplots(figsize=(10, 6))
    upset_obj2 = UpSet(upset_data2, subset_size='count',
                       intersection_plot_elements=6,
                       sort_by='cardinality', sort_category_by='cardinality')
    upset_obj2.plot(ax=ax)
    plt.savefig(r'H:\2026try\4.20\JMV\manuscript\figures\test_upset_b.png', dpi=150, bbox_inches='tight')
    print("Saved test_upset_b.png")
    plt.close()
except Exception as e:
    print(f"FAILED: {type(e).__name__}: {e}")

# Approach C: from_contents using string items (Subtype+Context combined)
print("\n=== Approach C: from_contents with string items ===")
mutation_items3 = {}
for _, row in df.iterrows():
    mut = row['Mutation']
    item = f"{row['Subtype']}_{row['context_tier']}"
    if mut not in mutation_items3:
        mutation_items3[mut] = set()
    mutation_items3[mut].add(item)
mutation_items3 = {k: list(v) for k, v in mutation_items3.items()}
print("mutation_items3:", dict(mutation_items3))

try:
    upset_data3 = from_contents(mutation_items3)
    print("from_contents OK")
    fig, ax = plt.subplots(figsize=(10, 6))
    upset_obj3 = UpSet(upset_data3, subset_size='count',
                       intersection_plot_elements=6,
                       sort_by='cardinality', sort_category_by='cardinality')
    upset_obj3.plot(ax=ax)
    plt.savefig(r'H:\2026try\4.20\JMV\manuscript\figures\test_upset_c.png', dpi=150, bbox_inches='tight')
    print("Saved test_upset_c.png")
    plt.close()
except Exception as e:
    print(f"FAILED: {type(e).__name__}: {e}")