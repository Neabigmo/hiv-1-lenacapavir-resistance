import pandas as pd
from pathlib import Path
from upsetplot import UpSet, from_contents
import upsetplot

print("upsetplot version:", upsetplot.__version__)

DATA_DIR = Path(r'H:\2026try\4.20\JMV\data\processed\revision_v2')
df = pd.read_csv(DATA_DIR / 'harmonized_phenotype_data.csv')

resistance_mutations = ['N57H', 'M66I', 'Q67H', 'Q67K', 'Q67H+K70R', 'Q67H+N74D',
                        'K70N', 'K70R', 'N74D', 'L56V', 'M66I+A105T', 'M66I+T107A',
                        'M66I+N74D', 'K70N+N74K', 'M66I+N74D+A105T']
df = df[df['Mutation'].isin(resistance_mutations)]

print("\nDataFrame shape:", df.shape)
print("Columns:", list(df.columns[:10]), "...")
print("\nMutation unique:", df['Mutation'].unique())
print("Subtype unique:", df['Subtype'].unique())
print("context_tier unique:", df['context_tier'].unique())

# Test 1: Series with MultiIndex
grouped = df.groupby(['Mutation', 'Subtype', 'context_tier']).size()
print("\n--- Test 1: Series with MultiIndex ---")
print("type:", type(grouped))
print("index type:", type(grouped.index))
print("index class:", grouped.index.__class__.__name__)
print("\nValue:\n", grouped.head(10))

try:
    upset_obj = UpSet(grouped, subset_size='count', intersection_plot_elements=6)
    print("UpSet created successfully!")
except Exception as e:
    print(f"UpSet failed: {e}")

# Test 2: from_contents format
print("\n--- Test 2: from_contents format ---")
from_contents_data = {}
for (mut, subtype, ctx), count in grouped.items():
    key = f"{mut}|{subtype}|{ctx}"
    from_contents_data[key] = list(range(count))
print("from_contents data keys:", list(from_contents_data.keys())[:5])

try:
    upset_obj2 = UpSet(from_contents(from_contents_data), subset_size='count', intersection_plot_elements=6)
    print("UpSet from_contents created successfully!")
except Exception as e:
    print(f"UpSet from_contents failed: {e}")

# Test 3: Binary incidence matrix
print("\n--- Test 3: Binary incidence matrix ---")
binary_df = df[['Mutation', 'Subtype', 'context_tier']].drop_duplicates()
binary_df['present'] = 1
pivot = binary_df.pivot_table(index='Mutation', columns=['Subtype', 'context_tier'], values='present', fill_value=0)
print("pivot shape:", pivot.shape)
print("pivot index type:", type(pivot.index))
print("\nPivot head:\n", pivot.head())

try:
    upset_obj3 = UpSet(pivot, subset_size='count', intersection_plot_elements=6)
    print("UpSet pivot created successfully!")
except Exception as e:
    print(f"UpSet pivot failed: {type(e).__name__}: {e}")