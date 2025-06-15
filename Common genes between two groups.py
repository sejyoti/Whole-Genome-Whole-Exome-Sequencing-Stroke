# ✅ STEP 1: Install necessary library
!pip install matplotlib-venn

# ✅ STEP 2: Import libraries
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn2

# ✅ STEP 3: Load annotated files
rec = pd.read_csv("X.csv")
nonrec = pd.read_csv("Y.csv")

# ✅ STEP 4: Extract unique gene names
genes_rec = set(rec['Gene.refGene'].dropna().unique())
genes_nonrec = set(nonrec['Gene.refGene'].dropna().unique())

# ✅ STEP 5: Compute intersections
common_genes = genes_rec & genes_nonrec
recurrent_only = genes_rec - genes_nonrec
nonrecurrent_only = genes_nonrec - genes_rec

# ✅ STEP 6: Save common genes to CSV
pd.DataFrame(sorted(common_genes), columns=["Common_Stroke_Genes"]).to_csv("common_genes.csv", index=False)

# ✅ STEP 7: Beautiful Venn diagram
plt.figure(figsize=(8, 6))
venn = venn2([genes_rec, genes_nonrec],
             set_labels=("Recurrent Event", "Non-Recurrent"),
             set_colors=("#FF6F61", "#6B5B95"),
             alpha=0.7)

# Custom font and title
plt.title("🧬 Venn Diagram", fontsize=14, fontweight='bold')
plt.grid(False)
plt.show()
