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

# Build unique items per mutation (sets to deduplicate)
mutation_items = {}
for _, row in df.iterrows():
    mut = row['Mutation']
    item = (row['Subtype'], row['context_tier'])
    if mut not in mutation_items:
        mutation_items[mut] = set()
    mutation_items[mut].add(item)
mutation_items = {k: list(v) for k, v in mutation_items.items()}
print("mutation_items:", dict(mutation_items))
print()

# Test 1: Create UpSet and plot
print("--- Test 1: Create UpSet, plot to new figure ---")
try:
    upset_data = from_contents(mutation_items)
    print("from_contents OK, type:", type(upset_data))

    fig_test, ax_test = plt.subplots(figsize=(8, 6))
    upset_obj = UpSet(upset_data, subset_size='count',
                      intersection_plot_elements=6,
                      sort_by='cardinality')
    print("UpSet object created OK")
    upset_obj.plot(ax=ax_test)
    print("plot() OK")
    fig_test.savefig(r'H:\2026try\4.20\JMV\manuscript\figures\test_upset_plot.png', dpi=150, bbox_inches='tight')
    print("Saved test_upset_plot.png")
    plt.close()
except Exception as e:
    import traceback
    print(f"FAILED: {type(e).__name__}: {e}")
    traceback.print_exc()

# Test 2: Plot onto a subplot of a larger figure (exactly like the script does)
print("\n--- Test 2: Plot to subplot of existing figure ---")
try:
    upset_data2 = from_contents(mutation_items)
    fig2 = plt.figure(figsize=(20, 15))
    # Simulate exactly what the script does: fig.add_subplot(2, 2, 1)
    ax_upset = fig2.add_subplot(2, 2, 1)

    upset_obj2 = UpSet(upset_data2, subset_size='count',
                        intersection_plot_elements=6,
                        sort_by='cardinality')
    print("UpSet object created OK")
    upset_obj2.plot(ax=ax_upset)
    print("plot() OK")
    fig2.savefig(r'H:\2026try\4.20\JMV\manuscript\figures\test_upset_subplot.png', dpi=150, bbox_inches='tight')
    print("Saved test_upset_subplot.png")
    plt.close()
except Exception as e:
    import traceback
    print(f"FAILED: {type(e).__name__}: {e}")
    traceback.print_exc()

# Test 3: Exact copy of the panel A code
print("\n--- Test 3: Exact panel A logic ---")
try:
    df2 = pd.read_csv(DATA_DIR / 'harmonized_phenotype_data.csv')
    df2 = df2[df2['Mutation'].isin(resistance_mutations)]

    fig3 = plt.figure(figsize=(20, 15))
    ax_upset3 = fig3.add_subplot(2, 2, 1)
    ax_upset3.set_title('A. Test Title', fontsize=11, weight='bold', loc='left')

    mutation_items3 = {}
    for _, row in df2.iterrows():
        mut = row['Mutation']
        item = (row['Subtype'], row['context_tier'])
        if mut not in mutation_items3:
            mutation_items3[mut] = set()
        mutation_items3[mut].add(item)
    mutation_items3 = {k: list(v) for k, v in mutation_items3.items()}

    upset_data3 = from_contents(mutation_items3)
    upset_obj3 = UpSet(upset_data3, subset_size='count',
                       intersection_plot_elements=6,
                       sort_by='cardinality')
    print("UpSet object created OK")
    upset_obj3.plot(ax=ax_upset3)
    print("plot() OK")
    fig3.savefig(r'H:\2026try\4.20\JMV\manuscript\figures\test_upset_exact.png', dpi=150, bbox_inches='tight')
    print("Saved test_upset_exact.png")
    plt.close()
except Exception as e:
    import traceback
    print(f"FAILED: {type(e).__name__}: {e}")
    traceback.print_exc()