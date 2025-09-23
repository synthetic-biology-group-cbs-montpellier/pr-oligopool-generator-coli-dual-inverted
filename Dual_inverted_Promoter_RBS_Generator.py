# %% [markdown]
# # Dual Inverted Promoter Library Selection and Analysis
# This notebook selects TWO DISTINCT SETS of genetic constructs for dual inverted promoter experiments, using flow-seq data and sequence features to minimize recombination risk.
# 
# ## Objectives
# - Load and preprocess flow‑seq and sequence data
# - Filter sequences to specified length range for practical synthesis
# - Create intensity bins spanning the full expression range
# - Select TWO maximally different sets of sequences to avoid recombination slippage
# - Use k-mer analysis to ensure sequence diversity within and between sets
# - Generate comprehensive analysis and visualization
# - Output separate CSV files for Set 1 and Set 2
# - Output synthesis ready sequences
# 
# ## Key Parameters (User Configurable)
# - `SEQUENCES_PER_SET`: Number of sequences in each set (default: 12)
# - `MIN_SEQUENCE_LENGTH`: Minimum sequence length in bp (default: 30)
# - `MAX_SEQUENCE_LENGTH`: Maximum sequence length in bp (default: 45)
# - `SELECTION_MODE`: How to distribute across intensity bins ("proportional", "equal", "hybrid")
# - `SEQUENCE_SIMILARITY_THRESHOLD`: Maximum allowed similarity between sets

# %% [markdown]
# ## Setup and Output Directory Creation
# Create a timestamped folder for this analysis run and copy the script.

# %%
import pandas as pd
import numpy as np
import re
import os
import shutil
from datetime import datetime

print("=== STARTING ANALYSIS ===")
print("Imports completed successfully")

# Create timestamped output directory
timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
output_dir = f"analysis_results_{timestamp}"
os.makedirs(output_dir, exist_ok=True)

# Copy this script to the output directory
script_name = "Dual_inverted.py"
if os.path.exists(script_name):
    shutil.copy2(script_name, os.path.join(output_dir, f"script_used_{timestamp}.py"))
    print(f"Script copied to: {output_dir}/script_used_{timestamp}.py")

print(f"Analysis results will be saved to: {output_dir}")

# %% [markdown]
# ## Data Loading and Preprocessing
# Load promoter, RBS, and construct measurements (flowsorted bins) from Excel files.
# 
# ```python
# import pandas as pd
# import numpy as np
# import re
# ```
# 
# Bin counts are in columns `bin.1` to `bin.12` in the constructs dataframe.

# %%
# Load data
print("\n=== DATA LOADING ===")
try:
    print("Loading constructs from sd03.xlsx...")
    constructs = pd.read_excel('sd03.xlsx', sheet_name='Constructs', engine='openpyxl')
    print(f"Constructs loaded: {len(constructs)} rows, {len(constructs.columns)} columns")
    
    print("Loading promoters from sd01.xlsx...")
    prom_df = pd.read_excel('sd01.xlsx', sheet_name='Promoters', engine='openpyxl')
    print(f"Promoters loaded: {len(prom_df)} rows, {len(prom_df.columns)} columns")
    
    print("Loading RBS from sd02.xlsx...")
    rbs_df = pd.read_excel('sd02.xlsx', sheet_name='RBSs', engine='openpyxl')
    print(f"RBS loaded: {len(rbs_df)} rows, {len(rbs_df.columns)} columns")
    
except Exception as e:
    print(f"ERROR loading Excel files: {e}")
    print("Please check that the following files exist:")
    print("- sd03.xlsx (with 'Constructs' sheet)")
    print("- sd01.xlsx (with 'Promoters' sheet)")  
    print("- sd02.xlsx (with 'RBSs' sheet)")
    raise

# Clean names
def clean_quotes(df, col):
    if col in df.columns:
        df[col] = df[col].str.replace('"', '').str.strip()
    else:
        print(f"WARNING: Column '{col}' not found in dataframe")

print("Cleaning column names...")
clean_quotes(prom_df, 'Promoter')
clean_quotes(rbs_df, 'RBS')
clean_quotes(constructs, 'Promoter')
clean_quotes(constructs, 'RBS')

# Save processed input data
constructs.to_csv(os.path.join(output_dir, "01_raw_constructs.csv"), index=False)
prom_df.to_csv(os.path.join(output_dir, "01_promoters.csv"), index=False)
rbs_df.to_csv(os.path.join(output_dir, "01_rbs.csv"), index=False)

print("Raw data saved to output directory")
print("Constructs columns:", list(constructs.columns))
print("Promoters columns:", list(prom_df.columns))
print("RBS columns:", list(rbs_df.columns))
constructs.head()

# %% [markdown]
# ## GFP Fluorescence Data Processing
# Use the 'prot' column values as GFP fluorescence measurements and filter out unmeasured data points.

# %%
# Check if 'prot' column exists
print("\n=== GFP FLUORESCENCE PROCESSING ===")
print("Checking for 'prot' column...")

if 'prot' not in constructs.columns:
    print("'prot' column not found. Available columns:", list(constructs.columns))
    # Look for similar column names
    prot_cols = [col for col in constructs.columns if 'prot' in col.lower()]
    if prot_cols:
        print(f"Found potential protein columns: {prot_cols}")
        constructs['prot'] = constructs[prot_cols[0]]
        print(f"Using {prot_cols[0]} as 'prot' column")
    else:
        raise ValueError("No 'prot' column found. Please check your data structure.")
else:
    print("'prot' column found successfully")

# Convert prot column to numeric and identify valid measurements
print("Converting 'prot' column to numeric...")
constructs['prot'] = pd.to_numeric(constructs['prot'], errors='coerce')

# Remove data points that are not measured experimentally
# Filter out NaN values and potentially invalid measurements (e.g., <= 0)
print("Filtering out unmeasured/invalid data...")
initial_count = len(constructs)
constructs = constructs[constructs['prot'].notna()].copy()
constructs = constructs[constructs['prot'] > 0].copy()  # Remove zero or negative values
final_count = len(constructs)

print(f"Data filtering:")
print(f"  Initial constructs: {initial_count}")
print(f"  After removing unmeasured/invalid data: {final_count}")
print(f"  Removed: {initial_count - final_count} constructs")

if final_count == 0:
    raise ValueError("No valid protein measurements found after filtering!")

# Use prot column directly as GFP fluorescence
constructs['GFP_fluorescence'] = constructs['prot']

# Calculate additional metrics for analysis
constructs['log_GFP'] = np.log10(constructs['GFP_fluorescence'])

print(f"GFP fluorescence statistics:")
print(f"  Mean: {constructs['GFP_fluorescence'].mean():.3f}")
print(f"  Std: {constructs['GFP_fluorescence'].std():.3f}")
print(f"  Range: {constructs['GFP_fluorescence'].min():.3f} - {constructs['GFP_fluorescence'].max():.3f}")
print(f"  Log10 range: {constructs['log_GFP'].min():.3f} - {constructs['log_GFP'].max():.3f}")

# Save constructs with GFP fluorescence data
constructs.to_csv(os.path.join(output_dir, "02_constructs_with_GFP.csv"), index=False)
print("Constructs with GFP fluorescence data saved")

# %% [markdown]
# ## Sequence Reconstruction and Processing
# Reconstruct full sequences by looking up promoter and RBS sequences, then trim to functional regions.

# %%
# Step 3: Reconstitute concatenated sequences using promoter and RBS identifiers
print("\n=== SEQUENCE RECONSTRUCTION ===")

# Clean promoter and RBS names to ensure matching
def clean_name(name):
    if pd.isna(name):
        return name
    return str(name).replace('"', '').strip()

print("Cleaning promoter and RBS names...")
# Clean all names
constructs['Promoter'] = constructs['Promoter'].apply(clean_name)
constructs['RBS'] = constructs['RBS'].apply(clean_name)
prom_df['Promoter'] = prom_df['Promoter'].apply(clean_name)
rbs_df['RBS'] = rbs_df['RBS'].apply(clean_name)

def clean_sequence(seq):
    """Clean sequence by removing quotes and whitespace"""
    if pd.isna(seq):
        return seq
    return re.sub(r'\s+', '', str(seq).replace('"', '').replace("'", '').strip())

# Clean sequences before creating lookup dictionaries
print("Cleaning sequence data...")
prom_df['Sequence'] = prom_df['Sequence'].apply(clean_sequence)
rbs_df['Sequence'] = rbs_df['Sequence'].apply(clean_sequence)

# Create lookup dictionaries
print("Creating lookup dictionaries...")
promoter_sequences = dict(zip(prom_df['Promoter'], prom_df['Sequence']))
rbs_sequences = dict(zip(rbs_df['RBS'], rbs_df['Sequence']))

print(f"Available promoters: {len(promoter_sequences)}")
print(f"Available RBS: {len(rbs_sequences)}")

# Check for potential sequence column issues
if 'Sequence' not in prom_df.columns:
    print("WARNING: 'Sequence' column not found in promoters dataframe")
    print("Promoter columns:", list(prom_df.columns))
if 'Sequence' not in rbs_df.columns:
    print("WARNING: 'Sequence' column not found in RBS dataframe")
    print("RBS columns:", list(rbs_df.columns))

# Reconstruct sequences
print("Reconstructing sequences...")
def reconstruct_sequence(row):
    prom_seq = promoter_sequences.get(row['Promoter'])
    rbs_seq = rbs_sequences.get(row['RBS'])
    
    if pd.isna(prom_seq) or pd.isna(rbs_seq):
        return None
    
    return str(prom_seq) + str(rbs_seq)

# Apply reconstruction with progress tracking
print(f"Processing {len(constructs)} constructs...")
constructs['full_sequence'] = constructs.apply(reconstruct_sequence, axis=1)

# Remove constructs without reconstructed sequences
initial_seq_count = len(constructs)
constructs = constructs[constructs['full_sequence'].notna()].copy()
final_seq_count = len(constructs)

print(f"Sequence reconstruction:")
print(f"  Successfully reconstructed: {final_seq_count}")
print(f"  Failed to reconstruct: {initial_seq_count - final_seq_count}")

if final_seq_count == 0:
    raise ValueError("No sequences could be reconstructed! Check promoter/RBS matching.")

# Step 4: Trim sequences between GGCGCGCC (5') and CATATG (3')
print("Trimming sequences...")
def trim_sequence(seq):
    if pd.isna(seq):
        return None
    
    seq = str(seq).upper()
    start_motif = "GGCGCGCC"
    end_motif = "CATATG"
    
    # Find start position (after the motif)
    start_pos = seq.find(start_motif)
    if start_pos == -1:
        # If exact motif not found, skip trimming from 5'
        start_pos = 0
    else:
        start_pos += len(start_motif)
    
    # Find end position (before the motif)
    end_pos = seq.find(end_motif, start_pos)
    if end_pos == -1:
        # If exact motif not found, keep until end
        end_pos = len(seq)
    
    trimmed = seq[start_pos:end_pos]
    return trimmed if trimmed else None

constructs['trimmed_sequence'] = constructs['full_sequence'].apply(trim_sequence)

# Remove constructs with empty trimmed sequences
constructs = constructs[constructs['trimmed_sequence'].notna()].copy()
constructs = constructs[constructs['trimmed_sequence'].str.len() > 0].copy()
trimmed_count = len(constructs)

print(f"Sequence trimming:")
print(f"  Sequences after trimming: {trimmed_count}")
if trimmed_count > 0:
    print(f"  Average trimmed length: {constructs['trimmed_sequence'].str.len().mean():.1f} bp")
    print(f"  Length range: {constructs['trimmed_sequence'].str.len().min()} - {constructs['trimmed_sequence'].str.len().max()} bp")

# Add sequence length column
constructs['sequence_length'] = constructs['trimmed_sequence'].str.len()

# Filter sequences to desired length range for practical synthesis
MIN_SEQUENCE_LENGTH = 50  # Minimum sequence length for dual inverted promoter synthesis
MAX_SEQUENCE_LENGTH = 70  # Maximum sequence length for dual inverted promoter synthesis
print(f"Filtering sequences to length range {MIN_SEQUENCE_LENGTH}-{MAX_SEQUENCE_LENGTH}bp for practical synthesis...")
initial_count_before_length_filter = len(constructs)
constructs = constructs[
    (constructs['sequence_length'] >= MIN_SEQUENCE_LENGTH) & 
    (constructs['sequence_length'] <= MAX_SEQUENCE_LENGTH)
].copy()
final_count_after_length_filter = len(constructs)

print(f"Length filtering:")
print(f"  Before length filter: {initial_count_before_length_filter}")
print(f"  After length filter ({MIN_SEQUENCE_LENGTH}-{MAX_SEQUENCE_LENGTH}bp): {final_count_after_length_filter}")
print(f"  Removed: {initial_count_before_length_filter - final_count_after_length_filter} sequences outside range")

if final_count_after_length_filter == 0:
    raise ValueError(f"No sequences remain after filtering for length {MIN_SEQUENCE_LENGTH}-{MAX_SEQUENCE_LENGTH}bp!")

# Update sequence length statistics after filtering
print(f"Updated sequence length range: {constructs['sequence_length'].min()} - {constructs['sequence_length'].max()} bp")

# Save processed sequences
constructs.to_csv(os.path.join(output_dir, "03_sequences_reconstructed.csv"), index=False)
print("Reconstructed, trimmed, and length-filtered sequences saved")

# %% [markdown]
# ## Intensity-Based Binning and Diverse Selection
# Create bins based on protein fluorescence intensity and select diverse variants within each bin.

# %%
from itertools import product
from collections import Counter
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.metrics import pairwise_distances
import random

# Step 5: Establish 12 bins based on protein fluorescence intensity
print("=== INTENSITY BINNING ===")

n_bins = 12
USE_LOG_SCALE_BINNING = True  # User can change this - True for log10 scale, False for linear scale

if USE_LOG_SCALE_BINNING:
    print("Using log10 scale binning for better dynamic range coverage")
    # Create bins on log10 scale
    log_min = np.log10(constructs['GFP_fluorescence'].min())
    log_max = np.log10(constructs['GFP_fluorescence'].max())
    log_bins = np.logspace(log_min, log_max, n_bins + 1)
    constructs['intensity_bin'] = pd.cut(constructs['GFP_fluorescence'], 
                                       bins=log_bins, 
                                       labels=range(1, n_bins+1),
                                       include_lowest=True)
    print(f"Log10 bins range: {log_min:.2f} to {log_max:.2f}")
else:
    print("Using linear scale binning")
    constructs['intensity_bin'] = pd.cut(constructs['GFP_fluorescence'], 
                                       bins=n_bins, 
                                       labels=range(1, n_bins+1),
                                       include_lowest=True)

# Print bin statistics
print("Bin boundaries and counts:")
bin_info = []
total_constructs = len(constructs)
for i in range(1, n_bins+1):
    bin_data = constructs[constructs['intensity_bin'] == i]
    if len(bin_data) > 0:
        min_val = bin_data['GFP_fluorescence'].min()
        max_val = bin_data['GFP_fluorescence'].max()
        count = len(bin_data)
        percentage = (count / total_constructs) * 100
        bin_info.append({
            'Bin': i,
            'Count': count,
            'Percentage': percentage,
            'Min_Intensity': min_val,
            'Max_Intensity': max_val,
            'Mean_Intensity': bin_data['GFP_fluorescence'].mean(),
            'Binning_Scale': 'Log10' if USE_LOG_SCALE_BINNING else 'Linear'
        })
        if USE_LOG_SCALE_BINNING:
            print(f"  Bin {i:2d}: {count:4d} variants ({percentage:5.1f}%), intensity {min_val:8.1f}-{max_val:8.1f} (log scale)")
        else:
            print(f"  Bin {i:2d}: {count:4d} variants ({percentage:5.1f}%), intensity {min_val:8.1f}-{max_val:8.1f} (linear scale)")

# Analyze the overall distribution
print(f"\n=== INTENSITY DISTRIBUTION ANALYSIS ===")
print(f"Total constructs: {total_constructs}")
print(f"Intensity range: {constructs['GFP_fluorescence'].min():.1f} - {constructs['GFP_fluorescence'].max():.1f}")
print(f"Median intensity: {constructs['GFP_fluorescence'].median():.1f}")
print(f"Mean intensity: {constructs['GFP_fluorescence'].mean():.1f}")

# Analyze specific ranges
ranges_of_interest = [
    (0, 20000, "Very Low"),
    (20000, 50000, "Low"), 
    (50000, 100000, "Medium"),
    (100000, 150000, "High"),
    (150000, float('inf'), "Very High")
]

print(f"\n=== VARIANTS BY INTENSITY RANGE ===")
for min_val, max_val, label in ranges_of_interest:
    if max_val == float('inf'):
        mask = constructs['GFP_fluorescence'] >= min_val
        range_str = f"≥{min_val:,}"
    else:
        mask = (constructs['GFP_fluorescence'] >= min_val) & (constructs['GFP_fluorescence'] < max_val)
        range_str = f"{min_val:,}-{max_val:,}"
    
    count = mask.sum()
    percentage = (count / total_constructs) * 100
    print(f"  {label:10s} ({range_str:15s}): {count:4d} variants ({percentage:5.1f}%)")

bin_info_df = pd.DataFrame(bin_info)
bin_info_df.to_csv(os.path.join(output_dir, "04_bin_information.csv"), index=False)

# Create early visualization of intensity distribution
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend to prevent popup windows
import matplotlib.pyplot as plt
plt.figure(figsize=(15, 10))

# Log-scale histogram to better see the distribution
plt.subplot(2, 2, 1)
plt.hist(constructs['GFP_fluorescence'], bins=50, alpha=0.7, color='skyblue', edgecolor='black')
plt.xlabel('GFP Fluorescence')
plt.ylabel('Count')
plt.title('Intensity Distribution (Linear Scale)')
plt.yscale('log')

# Linear scale
plt.subplot(2, 2, 2)
plt.hist(constructs['GFP_fluorescence'], bins=50, alpha=0.7, color='lightcoral', edgecolor='black')
plt.xlabel('GFP Fluorescence')
plt.ylabel('Count')
plt.title('Intensity Distribution (Linear Scale)')

# Log intensity histogram
plt.subplot(2, 2, 3)
plt.hist(constructs['log_GFP'], bins=50, alpha=0.7, color='lightgreen', edgecolor='black')
plt.xlabel('Log10(GFP Fluorescence)')
plt.ylabel('Count')
plt.title('Log10 Intensity Distribution')

# Binning visualization
plt.subplot(2, 2, 4)
bin_counts = [len(constructs[constructs['intensity_bin'] == i]) for i in range(1, n_bins+1)]
plt.bar(range(1, n_bins+1), bin_counts, alpha=0.7, color='orange', edgecolor='black')
plt.xlabel('Intensity Bin')
plt.ylabel('Count')
plt.title('Variants per Intensity Bin')
plt.xticks(range(1, n_bins+1))

plt.tight_layout()
plt.savefig(os.path.join(output_dir, "04_intensity_distribution_analysis.png"), dpi=300, bbox_inches='tight')
plt.close()

print("Intensity distribution analysis saved")

# Step 6: K-mer and PCA analysis for diversity selection within bins
print("\n=== K-MER DIVERSITY ANALYSIS ===")

# Optional: Limit constructs for faster processing during development
# Uncomment the lines below to test with fewer constructs
# if len(constructs) > 1000:
#     print(f"TESTING MODE: Limiting to first 1000 constructs (out of {len(constructs)}) for faster processing")
#     constructs = constructs.head(1000).copy()

# Function to count k-mers
def kmer_counts(seq, k=4):
    counts = Counter([seq[i:i+k] for i in range(len(seq)-k+1)])
    return counts

# Create k-mer feature matrix
kmers = [''.join(p) for p in product('ACGT', repeat=4)]
print(f"Analyzing {len(kmers)} 4-mers for {len(constructs)} constructs...")

kmat = np.zeros((len(constructs), len(kmers)))
print("Building k-mer matrix...")

# Process in batches to show progress
batch_size = 100
print("Building k-mer matrix...")
try:
    for i in range(0, len(constructs), batch_size):
        end_i = min(i + batch_size, len(constructs))
        if i % 1000 == 0 or i == 0:
            print(f"  Processing constructs {i+1}-{end_i} ({i/len(constructs)*100:.1f}%)")
        
        for j in range(i, end_i):
            seq = constructs.iloc[j]['trimmed_sequence']
            if pd.notna(seq) and len(str(seq)) > 0:
                cnts = kmer_counts(str(seq))
                for k_idx, k in enumerate(kmers):
                    kmat[j, k_idx] = cnts.get(k, 0)
    
    print("K-mer matrix completed successfully")
    
except Exception as e:
    print(f"Error during k-mer analysis: {e}")
    print("Attempting to continue with available data...")

# Normalize k-mer matrix
print("Normalizing k-mer matrix...")
try:
    # Check for empty rows (sequences with no valid k-mers)
    row_sums = kmat.sum(axis=1)
    valid_rows = row_sums > 0
    
    if valid_rows.sum() < len(constructs):
        print(f"Warning: {len(constructs) - valid_rows.sum()} sequences have no valid k-mers")
    
    kmat_norm = np.zeros_like(kmat)
    kmat_norm[valid_rows] = kmat[valid_rows] / (row_sums[valid_rows, np.newaxis] + 1e-10)
    
    print("K-mer normalization completed")
except Exception as e:
    print(f"Error during normalization: {e}")
    # Fallback: use simple normalization
    kmat_norm = kmat / (kmat.sum(axis=1, keepdims=True) + 1e-10)

# Add k-mer data to constructs for diversity selection
constructs.reset_index(drop=True, inplace=True)
print("K-mer analysis completed")

# Step 6.5: PCA Analysis for Sequence Diversity Assessment
print("\n=== PCA ANALYSIS FOR SEQUENCE DIVERSITY ===")
try:
    from sklearn.decomposition import PCA
    from sklearn.preprocessing import StandardScaler
    
    # Perform PCA on k-mer features
    print("Performing PCA on k-mer features...")
    scaler = StandardScaler()
    kmat_scaled = scaler.fit_transform(kmat_norm)
    
    # PCA with components explaining 95% variance
    pca = PCA(n_components=0.95)
    pca_features = pca.fit_transform(kmat_scaled)
    
    print(f"PCA completed: {pca.n_components_} components explain {pca.explained_variance_ratio_.sum():.3f} of variance")
    
    # Add PCA coordinates to constructs
    for i in range(min(5, pca.n_components_)):  # Store first 5 PC coordinates
        constructs[f'PC{i+1}'] = pca_features[:, i]

except ImportError:
    print("scikit-learn not available - skipping PCA analysis")
    # Create dummy PC columns
    for i in range(5):
        constructs[f'PC{i+1}'] = 0.0
except Exception as e:
    print(f"Error in PCA analysis: {e}")
    # Create dummy PC columns
    for i in range(5):
        constructs[f'PC{i+1}'] = 0.0

# Step 6.6: PCA Clustering for Visualization
print("\n=== PCA CLUSTERING FOR VISUALIZATION ===")
try:
    from sklearn.cluster import KMeans
    from sklearn.metrics import silhouette_score
    
    # Determine optimal number of clusters using silhouette analysis
    print("Finding optimal number of PCA clusters...")
    
    # Test different numbers of clusters
    cluster_range = range(3, min(15, len(constructs)//10))  # Reasonable range
    silhouette_scores = []
    
    pca_data_for_clustering = constructs[['PC1', 'PC2', 'PC3', 'PC4', 'PC5']].fillna(0).values
    
    for n_clusters in cluster_range:
        kmeans = KMeans(n_clusters=n_clusters, random_state=42, n_init=10)
        cluster_labels = kmeans.fit_predict(pca_data_for_clustering)
        silhouette_avg = silhouette_score(pca_data_for_clustering, cluster_labels)
        silhouette_scores.append(silhouette_avg)
        print(f"  {n_clusters} clusters: silhouette score = {silhouette_avg:.3f}")
    
    # Choose optimal number of clusters
    best_n_clusters = cluster_range[np.argmax(silhouette_scores)]
    best_silhouette = max(silhouette_scores)
    
    print(f"Optimal number of clusters: {best_n_clusters} (silhouette score: {best_silhouette:.3f})")
    
    # Apply final clustering
    final_kmeans = KMeans(n_clusters=best_n_clusters, random_state=42, n_init=10)
    constructs['PCA_Cluster'] = final_kmeans.fit_predict(pca_data_for_clustering)
    
    print(f"PCA clustering completed with {best_n_clusters} clusters")
    
    # Analyze clusters
    cluster_info = []
    for cluster_id in range(best_n_clusters):
        cluster_data = constructs[constructs['PCA_Cluster'] == cluster_id]
        cluster_info.append({
            'Cluster': cluster_id,
            'Size': len(cluster_data),
            'Intensity_Mean': cluster_data['GFP_fluorescence'].mean(),
            'Intensity_Std': cluster_data['GFP_fluorescence'].std(),
            'Intensity_Range': cluster_data['GFP_fluorescence'].max() - cluster_data['GFP_fluorescence'].min(),
            'Min_Intensity': cluster_data['GFP_fluorescence'].min(),
            'Max_Intensity': cluster_data['GFP_fluorescence'].max(),
            'PC1_Mean': cluster_data['PC1'].mean(),
            'PC2_Mean': cluster_data['PC2'].mean()
        })
    
    cluster_analysis_df = pd.DataFrame(cluster_info)
    print("Cluster analysis:")
    for _, row in cluster_analysis_df.iterrows():
        print(f"  Cluster {row['Cluster']}: {row['Size']} sequences, "
              f"intensity {row['Min_Intensity']:.0f}-{row['Max_Intensity']:.0f}")

except ImportError:
    print("scikit-learn not available - using single cluster")
    constructs['PCA_Cluster'] = 0
    best_n_clusters = 1
    best_silhouette = 0.0
except Exception as e:
    print(f"Error in PCA clustering: {e}")
    constructs['PCA_Cluster'] = 0
    best_n_clusters = 1
    best_silhouette = 0.0

# %% [markdown]
## 2D Grid-Based PCA/Intensity Selection
# Replace clustering approach with 2D grid-based selection

# %%
print("\n=== 2D GRID-BASED PCA/INTENSITY SELECTION ===")
print("New approach: 8x8 intensity grid with PCA-diverse pair selection")

# User-configurable parameters
SEQUENCES_PER_SET = 44
MIN_SEQUENCE_LENGTH = 20      
MAX_SEQUENCE_LENGTH = 70     
N_BINS = 50 # N X N grid
PCA_COLS = ['PC1', 'PC2', 'PC3', 'PC4', 'PC5']

# Set random seed for reproducibility
np.random.seed(42)

# Use all valid sequences (already filtered for length and experimentally observed)
all_valid_constructs = constructs.copy()
print(f"Starting with {len(all_valid_constructs)} valid sequences")

# Create 2D intensity grid
print(f"\n=== CREATING {N_BINS}x{N_BINS} INTENSITY GRID ===")

# Calculate log-intensity for both dimensions (using same values for symmetric grid)
all_valid_constructs['log_I1'] = np.log10(all_valid_constructs['GFP_fluorescence'])
all_valid_constructs['log_I2'] = all_valid_constructs['log_I1']  # Same for both dimensions

# Build bin edges
edges1 = np.linspace(all_valid_constructs['log_I1'].min(), all_valid_constructs['log_I1'].max(), N_BINS + 1)
edges2 = np.linspace(all_valid_constructs['log_I2'].min(), all_valid_constructs['log_I2'].max(), N_BINS + 1)

# Assign cells to each construct
all_valid_constructs['cell1'] = np.digitize(all_valid_constructs['log_I1'], edges1) - 1
all_valid_constructs['cell2'] = np.digitize(all_valid_constructs['log_I2'], edges2) - 1
all_valid_constructs[['cell1','cell2']] = all_valid_constructs[['cell1','cell2']].clip(0, N_BINS-1)

print(f"Intensity range: {all_valid_constructs['log_I1'].min():.3f} - {all_valid_constructs['log_I1'].max():.3f} (log10)")

# Find PCA-diverse pairs in each grid cell
print(f"\n=== FINDING PCA-DIVERSE PAIRS IN EACH CELL ===")

from scipy.spatial.distance import cdist
dist_pairs = {}
grid_stats = []

for i in range(N_BINS):
    for j in range(N_BINS):
        # Get constructs in this cell
        cell_data = all_valid_constructs[
            (all_valid_constructs['cell1'] == i) & 
            (all_valid_constructs['cell2'] == j)
        ]
        
        if len(cell_data) == 0:
            grid_stats.append({'cell1': i, 'cell2': j, 'count': 0, 'has_pair': False})
            continue
        elif len(cell_data) == 1:
            # Only one construct - use it for both sets
            construct = cell_data.iloc[0]
            dist_pairs[(i,j)] = (construct, construct)
            grid_stats.append({'cell1': i, 'cell2': j, 'count': 1, 'has_pair': True, 'distance': 0.0})
        else:
            # Find most PCA-diverse pair
            pca_data = cell_data[PCA_COLS].values
            D = cdist(pca_data, pca_data)
            
            # Find indices of maximum distance pair
            a, b = np.unravel_index(D.argmax(), D.shape)
            max_distance = D[a, b]
            
            construct1 = cell_data.iloc[a]
            construct2 = cell_data.iloc[b]
            dist_pairs[(i,j)] = (construct1, construct2)
            grid_stats.append({'cell1': i, 'cell2': j, 'count': len(cell_data), 'has_pair': True, 'distance': max_distance})

print(f"Filled {len(dist_pairs)} cells with construct pairs")

# Backfill empty cells using nearest neighbor
print(f"\n=== BACKFILLING EMPTY CELLS ===")
filled_cells = np.array(list(dist_pairs.keys()))
empty_cells = [(i,j) for i in range(N_BINS) for j in range(N_BINS) if (i,j) not in dist_pairs]

print(f"Empty cells to backfill: {len(empty_cells)}")

for empty_cell in empty_cells:
    # Find nearest filled cell using Manhattan distance
    distances = np.abs(filled_cells - np.array(empty_cell)).sum(axis=1)
    nearest_idx = np.argmin(distances)
    nearest_cell = tuple(filled_cells[nearest_idx])
    
    # Copy the pair from nearest cell
    dist_pairs[empty_cell] = dist_pairs[nearest_cell]
    print(f"  Cell {empty_cell} <- Cell {nearest_cell} (distance: {distances[nearest_idx]})")

# Extract Set 1 and Set 2 from pairs
print(f"\n=== EXTRACTING SETS FROM PAIRS ===")
set1_constructs = []
set2_constructs = []

for (i,j), (construct1, construct2) in dist_pairs.items():
    set1_constructs.append(construct1.to_dict())
    set2_constructs.append(construct2.to_dict())

final_set1 = pd.DataFrame(set1_constructs).drop_duplicates().reset_index(drop=True)
final_set2 = pd.DataFrame(set2_constructs).drop_duplicates().reset_index(drop=True)

# Limit to target size if needed
if len(final_set1) > SEQUENCES_PER_SET:
    final_set1 = final_set1.sample(n=SEQUENCES_PER_SET, random_state=42).reset_index(drop=True)
if len(final_set2) > SEQUENCES_PER_SET:
    final_set2 = final_set2.sample(n=SEQUENCES_PER_SET, random_state=42).reset_index(drop=True)

print(f"Final Set 1: {len(final_set1)} sequences")
print(f"Final Set 2: {len(final_set2)} sequences")

# Create grid statistics DataFrame for analysis
grid_stats_df = pd.DataFrame(grid_stats)

# Calculate diversity metrics
if len(final_set1) > 0 and len(final_set2) > 0:
    set1_pca = final_set1[PCA_COLS].values
    set2_pca = final_set2[PCA_COLS].values
    
    inter_distances = cdist(set1_pca, set2_pca, metric='euclidean')
    mean_inter_distance = np.mean(inter_distances)
    
    print(f"\n=== DIVERSITY ANALYSIS ===")
    print(f"Mean inter-set PCA distance: {mean_inter_distance:.3f}")
    
    # Intensity coverage
    set1_range = final_set1['GFP_fluorescence'].max() - final_set1['GFP_fluorescence'].min()
    set2_range = final_set2['GFP_fluorescence'].max() - final_set2['GFP_fluorescence'].min()
    total_range = all_valid_constructs['GFP_fluorescence'].max() - all_valid_constructs['GFP_fluorescence'].min()
    
    print(f"Set 1 intensity coverage: {set1_range/total_range*100:.1f}% of total range")
    print(f"Set 2 intensity coverage: {set2_range/total_range*100:.1f}% of total range")

# Use real cluster info from PCA clustering
if 'cluster_analysis_df' in locals():
    cluster_info_df = cluster_analysis_df.copy()
else:
    # Fallback: create cluster info from actual data
    cluster_info_df = pd.DataFrame([{
        'Cluster': 0,
        'Size': len(all_valid_constructs),
        'Intensity_Mean': all_valid_constructs['GFP_fluorescence'].mean(),
        'Intensity_Std': all_valid_constructs['GFP_fluorescence'].std(),
        'Intensity_Range': all_valid_constructs['GFP_fluorescence'].max() - all_valid_constructs['GFP_fluorescence'].min(),
        'Min_Intensity': all_valid_constructs['GFP_fluorescence'].min(),
        'Max_Intensity': all_valid_constructs['GFP_fluorescence'].max()
    }])

# Keep the actual PCA cluster assignments from the grid selection
# (final_set1 and final_set2 already have the correct PCA_Cluster values from all_valid_constructs)

# Define variables for compatibility with downstream code
set1_candidates = final_set1.copy()
set2_candidates = final_set2.copy()
set1_clusters = [0]
set2_clusters = [1]
min_separation = mean_inter_distance if len(final_set1) > 0 and len(final_set2) > 0 else 0
pca_features = PCA_COLS

print(f"\n=== 2D GRID SELECTION COMPLETE ===")
print(f"Grid approach used {N_BINS}x{N_BINS} = {N_BINS*N_BINS} cells")
print(f"Selected {len(final_set1)} + {len(final_set2)} = {len(final_set1) + len(final_set2)} total sequences")

# OLD COMPLEX SELECTION CODE REPLACED BY SIMPLIFIED APPROACH ABOVE
# This section has been replaced by the new PCA clustering approach
print("Old binning-based selection approach replaced by PCA clustering method")

# This section has been replaced by the new simplified PCA clustering approach above
print("Old PCA optimization section removed - using new approach")

# NOTE: Concatenation section moved to after selection - see after the simplified selection code

# %% [markdown]
# ## Simplified Dual Set Selection Using PCA Clustering
# New approach: Start with all valid sequences, cluster them, then select from different clusters

# %%
# OLD CLUSTERING CODE - SKIPPED FOR 2D GRID APPROACH
# The following clustering-based selection code has been replaced by the 2D grid approach above
# Jumping directly to results saving section

# Note: The variables final_set1, final_set2, cluster_info_df, etc. are already defined by the 2D grid approach

# Skip to results saving - all old clustering code removed

# %% [markdown]
## Results Saving and Analysis

# %%
print(f"\n=== SAVING RESULTS ===")
cluster_info_df.to_csv(os.path.join(output_dir, "05_cluster_analysis.csv"), index=False)
final_set1.to_csv(os.path.join(output_dir, "05_selected_variants_SET1.csv"), index=False)
final_set2.to_csv(os.path.join(output_dir, "05_selected_variants_SET2.csv"), index=False)

# Combined selection for compatibility
final_selection = pd.concat([final_set1, final_set2], ignore_index=True)
final_selection.to_csv(os.path.join(output_dir, "05_selected_variants_COMBINED.csv"), index=False)

print(f"Selection complete! Files saved:")
print(f"  05_cluster_analysis.csv")
print(f"  05_selected_variants_SET1.csv") 
print(f"  05_selected_variants_SET2.csv")
print(f"  05_selected_variants_COMBINED.csv")

# Step 9: Create Concatenated Sequences (Dual Inverted Promoters)
# Architecture: PRIMER_SITE_5_PRIME - BSAI_SITE_FWD - Reverse Comp Set 1 - Terminator - Forward Set 2- BSAI_SITE_REV - PRIMER_SITE_3_PRIME 
print(f"\n=== CREATING CONCATENATED DUAL INVERTED PROMOTER SEQUENCES ===")

# Sequence components for full architecture
Terminator_ECK120011170_AgeI = "ttcaacgagaaaagccaacctgcgggttggcttttttatgcaAccggt"
BSAI_SITE_FWD = "GGTCTCAtagt"  # BsaI forward recognition site
BSAI_SITE_REV = "actaTGAGACC"  # Reverse complement recognition site
PRIMER_SITE_5_PRIME = "AATCCTTGCGTCAATGGTTC"     # skpp-202-F from kosuri et al 2012
PRIMER_SITE_3_PRIME = "CGTGTAAAATCCGAGAACCC"     # skpp-202-R from kosuri et al 2012

print(f"Using ECK120011170_AgeI terminator: {Terminator_ECK120011170_AgeI}")
print(f"Terminator length: {len(Terminator_ECK120011170_AgeI)}bp")
print(f"5' Primer: {PRIMER_SITE_5_PRIME} ({len(PRIMER_SITE_5_PRIME)}bp)")
print(f"3' Primer: {PRIMER_SITE_3_PRIME} ({len(PRIMER_SITE_3_PRIME)}bp)")
print(f"BsaI Forward: {BSAI_SITE_FWD} ({len(BSAI_SITE_FWD)}bp)")
print(f"BsaI Reverse: {BSAI_SITE_REV} ({len(BSAI_SITE_REV)}bp)")


def sanitize_dna(seq):
    """Return sequence as a string with internal whitespace and quotes removed."""
    if pd.isna(seq):
        return ""
    return re.sub(r"\s+", "", str(seq)).replace('"', '').replace("'", '')

def reverse_complement(seq):
    """Generate reverse complement of a DNA sequence"""
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 
                  'a': 't', 't': 'a', 'g': 'c', 'c': 'g',
                  'N': 'N', 'n': 'n'}
    try:
        return ''.join(complement.get(base, base) for base in reversed(seq))
    except:
        return seq  # Return original if conversion fails

def create_concatenated_sequence(set1_seq, set2_seq, terminator):
    """Create concatenated sequence: RevComp(Set1) - Terminator - Set2"""
    set1_seq = sanitize_dna(set1_seq)
    set2_seq = sanitize_dna(set2_seq)
    terminator = sanitize_dna(terminator)
    rev_comp_set1 = reverse_complement(set1_seq)
    concatenated = rev_comp_set1 + terminator + set2_seq
    return concatenated, rev_comp_set1

def build_full_architecture(seq_set1, seq_set2, terminator,
                           primer1=PRIMER_SITE_5_PRIME, bsai1=BSAI_SITE_FWD,
                           primer2=PRIMER_SITE_3_PRIME, bsai2=BSAI_SITE_REV):
    """
    Build: primer1 - BsaI site 1 - RevComp(Set1) - Terminator - Set2 - BsaI site 2 - primer2
    All inputs will be sanitized to remove whitespace.
    Returns full_sequence, rev_comp_set1.
    """
    seq_set1 = sanitize_dna(seq_set1)
    seq_set2 = sanitize_dna(seq_set2)
    terminator = sanitize_dna(terminator)
    primer1 = sanitize_dna(primer1)
    primer2 = sanitize_dna(primer2)
    bsai1 = sanitize_dna(bsai1)
    bsai2 = sanitize_dna(bsai2)

    core_concat, rev_comp_set1 = create_concatenated_sequence(seq_set1, seq_set2, terminator)
    return primer1 + bsai1 + core_concat + bsai2 + primer2, rev_comp_set1

def analyze_restriction_sites(full_sequence, primer1, bsai1, rev_comp_set1, terminator, set2, bsai2, primer2):
    """
    Analyze restriction sites in concatenated sequence, focusing on junction regions.
    Only flags unexpected restriction sites, not intentional ones.
    """
    # Define problematic restriction sites to check
    problematic_sites = {
        'BsaI_forward': 'GGTCTC',
        'BsaI_reverse': 'GAGACC',
        'AgeI': 'ACCGGT',
        'AscI': 'GGCGCGCC',
        'SphI': 'GCATGC',
        'EcoRI': 'GAATTC',
        'BamHI': 'GGATCC',
        'HindIII': 'AAGCTT',
        'XhoI': 'CTCGAG',
        # 'SalI': 'GTCGAC',
        'NotI': 'GCGGCCGC',
        'SpeI': 'ACTAGT',
        # 'XbaI': 'TCTAGA' XbaI is in many promoters so not posible to exclude
    }
    
    # Calculate component positions in full sequence
    primer1_len = len(primer1)
    bsai1_len = len(bsai1)
    rev_comp_set1_len = len(rev_comp_set1)
    terminator_len = len(terminator)
    set2_len = len(set2)
    bsai2_len = len(bsai2)
    
    # Map expected positions of intentional sites
    expected_bsai_fwd_start = primer1_len
    expected_bsai_fwd_end = primer1_len + bsai1_len
    expected_bsai_rev_start = primer1_len + bsai1_len + rev_comp_set1_len + terminator_len + set2_len
    expected_bsai_rev_end = expected_bsai_rev_start + bsai2_len
    
    # Map expected terminator region for AgeI site
    expected_terminator_start = primer1_len + bsai1_len + rev_comp_set1_len
    expected_terminator_end = expected_terminator_start + terminator_len
    
    # Define junction windows (10bp each side of junction)
    junction_window = 10
    junctions = [
        (primer1_len - junction_window, primer1_len + junction_window, "Primer5-BsaI_Fwd"),
        (expected_bsai_fwd_end - junction_window, expected_bsai_fwd_end + junction_window, "BsaI_Fwd-RevComp_Set1"),
        (primer1_len + bsai1_len + rev_comp_set1_len - junction_window, 
         primer1_len + bsai1_len + rev_comp_set1_len + junction_window, "RevComp_Set1-Terminator"),
        (primer1_len + bsai1_len + rev_comp_set1_len + terminator_len - junction_window,
         primer1_len + bsai1_len + rev_comp_set1_len + terminator_len + junction_window, "Terminator-Set2"),
        (expected_bsai_rev_start - junction_window, expected_bsai_rev_start + junction_window, "Set2-BsaI_Rev"),
        (expected_bsai_rev_end - junction_window, expected_bsai_rev_end + junction_window, "BsaI_Rev-Primer3")
    ]
    
    # Analysis results
    issues_found = []
    sites_detected = {}
    
    # Check each restriction site
    for site_name, site_seq in problematic_sites.items():
        # Find all occurrences of this site
        positions = []
        seq_upper = full_sequence.upper()
        site_upper = site_seq.upper()
        
        # Search forward orientation
        start = 0
        while True:
            pos = seq_upper.find(site_upper, start)
            if pos == -1:
                break
            positions.append(pos)
            start = pos + 1
        
        # Search reverse complement
        site_rc = reverse_complement(site_seq).upper()
        start = 0
        while True:
            pos = seq_upper.find(site_rc, start)
            if pos == -1:
                break
            positions.append(pos)
            start = pos + 1
        
        if positions:
            sites_detected[site_name] = positions
            
            # Special handling for BsaI sites - only flag unexpected ones
            if site_name in ['BsaI_forward', 'BsaI_reverse']:
                unexpected_positions = []
                for pos in positions:
                    # Check if position is in expected BsaI regions
                    in_expected_fwd = expected_bsai_fwd_start <= pos < expected_bsai_fwd_end
                    in_expected_rev = expected_bsai_rev_start <= pos < expected_bsai_rev_end
                    
                    if not (in_expected_fwd or in_expected_rev):
                        # Check if it's in a junction region (more problematic)
                        in_junction = any(j_start <= pos < j_end for j_start, j_end, _ in junctions)
                        unexpected_positions.append((pos, in_junction))
                
                if unexpected_positions:
                    junction_issues = [pos for pos, in_junction in unexpected_positions if in_junction]
                    if junction_issues:
                        issues_found.append(f"{site_name} at junctions: {junction_issues}")
                    else:
                        issues_found.append(f"{site_name} outside expected regions: {[pos for pos, _ in unexpected_positions]}")
            
            # Special handling for AgeI sites - only flag unexpected ones (outside terminator)
            elif site_name == 'AgeI':
                unexpected_positions = []
                for pos in positions:
                    # Check if position is in expected terminator region
                    in_expected_terminator = expected_terminator_start <= pos < expected_terminator_end
                    
                    if not in_expected_terminator:
                        # Check if it's in a junction region (more problematic)
                        in_junction = any(j_start <= pos < j_end for j_start, j_end, _ in junctions)
                        unexpected_positions.append((pos, in_junction))
                
                if unexpected_positions:
                    junction_issues = [pos for pos, in_junction in unexpected_positions if in_junction]
                    if junction_issues:
                        issues_found.append(f"{site_name} at junctions: {junction_issues}")
                    else:
                        issues_found.append(f"{site_name} outside expected terminator: {[pos for pos, _ in unexpected_positions]}")
            
            else:
                # For other restriction sites, flag any occurrence
                junction_issues = []
                for pos in positions:
                    for j_start, j_end, j_name in junctions:
                        if j_start <= pos < j_end:
                            junction_issues.append((pos, j_name))
                
                if junction_issues:
                    issues_found.append(f"{site_name} in junctions: {junction_issues}")
                elif len(positions) > 0:
                    issues_found.append(f"{site_name} found at positions: {positions}")
    
    # Determine overall status
    if not issues_found:
        status = "PASS"
    elif any("junction" in issue.lower() for issue in issues_found):
        status = "FAIL"
    else:
        status = "WARN"
    
    return {
        'status': status,
        'issues': issues_found,
        'sites_detected': sites_detected,
        'total_issues': len(issues_found)
    }

# Create concatenated sequences using the selected sets
concat_df = pd.DataFrame()

if len(final_set1) > 0 and len(final_set2) > 0:
    print(f"\n=== CONCATENATED SEQUENCES CREATION ===")
    print(f"Creating concatenated sequences from {len(final_set1)} Set1 × {len(final_set2)} Set2 combinations...")
    
    concatenated_sequences = []
    combination_count = 0
    
    for i, set1_row in final_set1.iterrows():
        for j, set2_row in final_set2.iterrows():
            combination_count += 1
            
            # Get sequences
            set1_seq = set1_row['trimmed_sequence']
            set2_seq = set2_row['trimmed_sequence']
            
            # Get promoter and RBS names for proper naming
            set1_promoter = set1_row.get('Promoter', f'Prom_{i}')
            set1_rbs = set1_row.get('RBS', f'RBS_{i}')
            set2_promoter = set2_row.get('Promoter', f'Prom_{j}')
            set2_rbs = set2_row.get('RBS', f'RBS_{j}')
            
            # Create concatenated sequence (basic version)
            concat_seq, rev_comp_set1 = create_concatenated_sequence(set1_seq, set2_seq, Terminator_ECK120011170_AgeI)
            
            # Create full architecture sequence (with primers and BsaI sites)
            full_arch_seq, _ = build_full_architecture(set1_seq, set2_seq, Terminator_ECK120011170_AgeI)
            
            # Perform junction analysis for restriction sites
            junction_analysis = analyze_restriction_sites(
                full_arch_seq, 
                PRIMER_SITE_5_PRIME, 
                BSAI_SITE_FWD, 
                rev_comp_set1, 
                Terminator_ECK120011170_AgeI, 
                set2_seq, 
                BSAI_SITE_REV, 
                PRIMER_SITE_3_PRIME
            )
            
            # Calculate architecture-aware component properties
            total_length = len(full_arch_seq)  # Use synthesis-ready sequence length
            full_arch_length = len(full_arch_seq)  # Keep for backwards compatibility
            
            # Architecture component lengths (in synthesis order)
            arch_primer5_length = len(PRIMER_SITE_5_PRIME)
            arch_bsai_fwd_length = len(BSAI_SITE_FWD)
            arch_set1_revcomp_length = len(rev_comp_set1)  # Actual reverse complement used
            arch_terminator_length = len(Terminator_ECK120011170_AgeI)
            arch_set2_forward_length = len(set2_seq)
            arch_bsai_rev_length = len(BSAI_SITE_REV)
            arch_primer3_length = len(PRIMER_SITE_3_PRIME)
            
            # Original component lengths (for reference)
            original_set1_length = len(set1_seq)
            original_set2_length = len(set2_seq)
            
            # Verify architecture component lengths sum to total
            calculated_total = (arch_primer5_length + arch_bsai_fwd_length + 
                              arch_set1_revcomp_length + arch_terminator_length +
                              arch_set2_forward_length + arch_bsai_rev_length + 
                              arch_primer3_length)
            
            length_verification_passed = (calculated_total == total_length)
            length_discrepancy = abs(calculated_total - total_length) if not length_verification_passed else 0
            
            # Create construct name in format: Promoter1--RBS1-R-ECK120011170_AgeI-Promoter2--RBS2-F
            construct_name = f"{set1_promoter}--{set1_rbs}-R-ECK120011170_AgeI-{set2_promoter}--{set2_rbs}-F"
            
            # Store information
            concatenated_sequences.append({
                'Construct_Name': construct_name,
                'Full_Sequence': full_arch_seq,  # Use synthesis-ready sequence
                'Basic_Concatenated_Sequence': concat_seq,  # Keep basic sequence for reference
                'Total_Length': total_length,
                'Set1_Promoter': set1_promoter,
                'Set1_RBS': set1_rbs,
                'Set1_Original_Seq': set1_seq,
                'Set1_RevComp_Seq': rev_comp_set1,
                'Set1_Original_Length': original_set1_length,
                'Set1_Intensity': set1_row['GFP_fluorescence'],
                'Set1_Cluster': set1_row.get('PCA_Cluster', 'Unknown'),
                'Terminator_Seq': Terminator_ECK120011170_AgeI,
                'Terminator_Length': arch_terminator_length,
                'Set2_Promoter': set2_promoter,
                'Set2_RBS': set2_rbs,
                'Set2_Forward_Seq': set2_seq,
                'Set2_Original_Length': original_set2_length,
                'Set2_Intensity': set2_row['GFP_fluorescence'],
                'Set2_Cluster': set2_row.get('PCA_Cluster', 'Unknown'),
                'Full_Architecture_Sequence': full_arch_seq,
                'Full_Architecture_Length': full_arch_length,
                'Primer_5_Prime': PRIMER_SITE_5_PRIME,
                'Primer_5_Prime_Length': arch_primer5_length,
                'BsaI_Site_Fwd': BSAI_SITE_FWD,
                'BsaI_Site_Fwd_Length': arch_bsai_fwd_length,
                'BsaI_Site_Rev': BSAI_SITE_REV,
                'BsaI_Site_Rev_Length': arch_bsai_rev_length,
                'Primer_3_Prime': PRIMER_SITE_3_PRIME,
                'Primer_3_Prime_Length': arch_primer3_length,
                # Architecture-aware component lengths (synthesis order)
                'Arch_Primer5_Length': arch_primer5_length,
                'Arch_BsaI_Fwd_Length': arch_bsai_fwd_length,
                'Arch_Set1_RevComp_Length': arch_set1_revcomp_length,
                'Arch_Terminator_Length': arch_terminator_length,
                'Arch_Set2_Forward_Length': arch_set2_forward_length,
                'Arch_BsaI_Rev_Length': arch_bsai_rev_length,
                'Arch_Primer3_Length': arch_primer3_length,
                # Length verification
                'Length_Verification_Passed': length_verification_passed,
                'Length_Discrepancy': length_discrepancy,
                'Calculated_Component_Total': calculated_total,
                'Architecture': 'Primer5-BsaI_Fwd-RevComp(Set1)-Terminator-Set2-BsaI_Rev-Primer3',
                'Set_Type': 'PCA_Cluster_Optimized',
                'Expected_Expression_Combo': f"Set1:{set1_row['GFP_fluorescence']:.0f}_Set2:{set2_row['GFP_fluorescence']:.0f}",
                'Junction_Analysis_Status': junction_analysis['status'],
                'Restriction_Sites_Issues': '; '.join(junction_analysis['issues']) if junction_analysis['issues'] else 'None',
                'Total_Restriction_Issues': junction_analysis['total_issues'],
                'Passes_Junction_Analysis': junction_analysis['status'] == 'PASS'
            })
    
    print(f"Created {len(concatenated_sequences)} concatenated sequences")
    
    # Convert to DataFrame
    concat_df = pd.DataFrame(concatenated_sequences)
    
    # Report junction analysis results
    print(f"\n=== JUNCTION ANALYSIS RESULTS ===")
    if len(concat_df) > 0:
        pass_count = (concat_df['Junction_Analysis_Status'] == 'PASS').sum()
        warn_count = (concat_df['Junction_Analysis_Status'] == 'WARN').sum()
        fail_count = (concat_df['Junction_Analysis_Status'] == 'FAIL').sum()
        
        print(f"  PASS: {pass_count} sequences ({pass_count/len(concat_df)*100:.1f}%)")
        print(f"  WARN: {warn_count} sequences ({warn_count/len(concat_df)*100:.1f}%)")
        print(f"  FAIL: {fail_count} sequences ({fail_count/len(concat_df)*100:.1f}%)")
        
        if fail_count > 0:
            print(f"\n  Sequences with restriction site issues at junctions:")
            failed_sequences = concat_df[concat_df['Junction_Analysis_Status'] == 'FAIL']
            for idx, row in failed_sequences.head(5).iterrows():  # Show first 5 failed
                print(f"    {row['Construct_Name']}: {row['Restriction_Sites_Issues']}")
            if len(failed_sequences) > 5:
                print(f"    ... and {len(failed_sequences) - 5} more sequences with issues")
        
        print(f"  Sequences suitable for synthesis: {pass_count + warn_count} ({(pass_count + warn_count)/len(concat_df)*100:.1f}%)")
    
    # Report length verification results
    print(f"\n=== ARCHITECTURE COMPONENT LENGTH VERIFICATION ===")
    if len(concat_df) > 0:
        length_pass_count = (concat_df['Length_Verification_Passed'] == True).sum()
        length_fail_count = (concat_df['Length_Verification_Passed'] == False).sum()
        
        print(f"  PASS: {length_pass_count} sequences ({length_pass_count/len(concat_df)*100:.1f}%)")
        print(f"  FAIL: {length_fail_count} sequences ({length_fail_count/len(concat_df)*100:.1f}%)")
        
        if length_fail_count > 0:
            print(f"\n  Sequences with length verification failures:")
            failed_length = concat_df[concat_df['Length_Verification_Passed'] == False]
            for idx, row in failed_length.head(5).iterrows():  # Show first 5 failed
                print(f"    {row['Construct_Name']}: Total={row['Total_Length']}bp, Calculated={row['Calculated_Component_Total']}bp, Discrepancy={row['Length_Discrepancy']}bp")
            if len(failed_length) > 5:
                print(f"    ... and {len(failed_length) - 5} more sequences with length discrepancies")
        
        # Component length statistics
        if length_pass_count > 0:
            print(f"\n  Architecture component length breakdown (verified sequences):")
            verified_sequences = concat_df[concat_df['Length_Verification_Passed'] == True]
            print(f"    Average Primer5: {verified_sequences['Arch_Primer5_Length'].mean():.1f}bp")
            print(f"    Average BsaI_Fwd: {verified_sequences['Arch_BsaI_Fwd_Length'].mean():.1f}bp")
            print(f"    Average Set1_RevComp: {verified_sequences['Arch_Set1_RevComp_Length'].mean():.1f}bp")
            print(f"    Average Terminator: {verified_sequences['Arch_Terminator_Length'].mean():.1f}bp")
            print(f"    Average Set2_Forward: {verified_sequences['Arch_Set2_Forward_Length'].mean():.1f}bp")
            print(f"    Average BsaI_Rev: {verified_sequences['Arch_BsaI_Rev_Length'].mean():.1f}bp")
            print(f"    Average Primer3: {verified_sequences['Arch_Primer3_Length'].mean():.1f}bp")
    
    # Save concatenated sequences
    concat_df.to_csv(os.path.join(output_dir, "05_concatenated_dual_inverted_promoters.csv"), index=False)
    
    # Create simplified output for synthesis (using full architecture sequences)
    synthesis_output = concat_df[['Construct_Name', 'Full_Sequence', 'Total_Length', 
                                'Set1_Promoter', 'Set1_RBS', 'Set1_Intensity', 
                                'Set2_Promoter', 'Set2_RBS', 'Set2_Intensity', 
                                'Architecture', 'Set_Type', 'Junction_Analysis_Status',
                                'Passes_Junction_Analysis', 'Total_Restriction_Issues',
                                'Restriction_Sites_Issues', 'Length_Verification_Passed',
                                'Arch_Set1_RevComp_Length', 'Arch_Set2_Forward_Length',
                                'Arch_Terminator_Length']].copy()
    # Rename for clarity in synthesis output
    synthesis_output = synthesis_output.rename(columns={
        'Full_Sequence': 'Synthesis_Ready_Sequence',
        'Total_Length': 'Synthesis_Sequence_Length'
    })
    synthesis_output.to_csv(os.path.join(output_dir, "05_synthesis_ready_sequences.csv"), index=False)
    
    # Create filtered output with only sequences that pass junction analysis AND length verification
    clean_sequences = synthesis_output[
        (synthesis_output['Passes_Junction_Analysis'] == True) & 
        (synthesis_output['Length_Verification_Passed'] == True)
    ].copy()
    if len(clean_sequences) > 0:
        clean_sequences.to_csv(os.path.join(output_dir, "05_synthesis_ready_sequences_CLEAN.csv"), index=False)
        print(f"Concatenated sequences saved:")
        print(f"  Full details: 05_concatenated_dual_inverted_promoters.csv")
        print(f"  Synthesis ready (all): 05_synthesis_ready_sequences.csv")
        print(f"  Synthesis ready (CLEAN only): 05_synthesis_ready_sequences_CLEAN.csv ({len(clean_sequences)} sequences)")
    else:
        print(f"Concatenated sequences saved:")
        print(f"  Full details: 05_concatenated_dual_inverted_promoters.csv")
        print(f"  Synthesis ready: 05_synthesis_ready_sequences.csv")
        print(f"  WARNING: No sequences passed junction analysis - no CLEAN file created")
    
    # Verification: Read back synthesis-ready sequences and extract lengths from saved file
    print(f"\n=== SYNTHESIS SEQUENCE VERIFICATION ===")
    synthesis_verification = pd.read_csv(os.path.join(output_dir, "05_synthesis_ready_sequences.csv"))
    
    # Calculate lengths from actual saved sequences
    synthesis_verification['Verified_Length'] = synthesis_verification['Synthesis_Ready_Sequence'].str.len()
    
    print(f"Verification from synthesis-ready CSV:")
    print(f"  Total sequences verified: {len(synthesis_verification)}")
    print(f"  Average synthesis length: {synthesis_verification['Verified_Length'].mean():.1f}bp")
    print(f"  Synthesis length range: {synthesis_verification['Verified_Length'].min()}-{synthesis_verification['Verified_Length'].max()}bp")
    
    # Compare with calculated lengths during generation
    if 'Synthesis_Sequence_Length' in synthesis_verification.columns:
        length_match = (synthesis_verification['Verified_Length'] == synthesis_verification['Synthesis_Sequence_Length']).all()
        print(f"  Length verification: {'PASSED' if length_match else 'FAILED'} - Calculated vs. saved lengths match: {length_match}")
    
    # Statistics based on synthesis-ready sequences
    print(f"\n=== SYNTHESIS-READY SEQUENCE STATISTICS ===")
    print(f"  Total combinations: {len(synthesis_verification)}")
    print(f"  Average synthesis length: {synthesis_verification['Verified_Length'].mean():.1f}bp")
    print(f"  Synthesis length range: {synthesis_verification['Verified_Length'].min()}-{synthesis_verification['Verified_Length'].max()}bp")
    
    # Length distribution analysis based on synthesis-ready sequences
    length_ranges = [
        (0, 100, "Very Short (<100bp)"),
        (100, 150, "Short (100-150bp)"),
        (150, 200, "Medium (150-200bp)"),
        (200, 250, "Long (200-250bp)"),
        (250, float('inf'), "Very Long (>250bp)")
    ]
    
    print(f"\n=== SYNTHESIS LENGTH DISTRIBUTION ANALYSIS ===")
    for min_len, max_len, label in length_ranges:
        if max_len == float('inf'):
            mask = synthesis_verification['Verified_Length'] >= min_len
        else:
            mask = (synthesis_verification['Verified_Length'] >= min_len) & (synthesis_verification['Verified_Length'] < max_len)
        
        count = mask.sum()
        percentage = (count / len(synthesis_verification)) * 100 if len(synthesis_verification) > 0 else 0
        print(f"  {label}: {count} sequences ({percentage:.1f}%)")
    
    # Check for synthesis feasibility based on actual synthesis-ready sequences
    MAX_SYNTHESIS_LENGTH = 300
    too_long = synthesis_verification[synthesis_verification['Verified_Length'] > MAX_SYNTHESIS_LENGTH]
    
    if len(too_long) > 0:
        print(f"\nWARNING: {len(too_long)} sequences are longer than {MAX_SYNTHESIS_LENGTH}bp")
        print("Consider adjusting sequence length filters")
    else:
        print(f"\n✓ All sequences are ≤{MAX_SYNTHESIS_LENGTH}bp (suitable for synthesis)")

else:
    concat_df = pd.DataFrame()
    print("No valid sets for concatenation")

print(f"\n=== DUAL INVERTED PROMOTER LIBRARY ANALYSIS COMPLETE ===")
print(f"Results saved to: {output_dir}")
print(f"Total sequences selected: {len(final_set1)} + {len(final_set2)} = {len(final_set1) + len(final_set2)}")
if len(concat_df) > 0:
    print(f"Total concatenated combinations: {len(concat_df)}")

# Step 10: Create Comprehensive Visualization and Report
print(f"\n=== CREATING COMPREHENSIVE VISUALIZATION REPORT ===")

import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages
import warnings
warnings.filterwarnings('ignore')

# Set style for better-looking plots
plt.style.use('default')
sns.set_palette("husl")

# Create comprehensive PDF report
pdf_filename = os.path.join(output_dir, "06_dual_inverted_promoter_analysis_report.pdf")

with PdfPages(pdf_filename) as pdf:
    
    # Page 1: Overview and Data Distribution
    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    fig.suptitle('Dual Inverted Promoter Library Analysis - Overview', fontsize=16, fontweight='bold')
    
    # 1.1: Intensity distribution histogram
    axes[0,0].hist(all_valid_constructs['GFP_fluorescence'], bins=50, alpha=0.7, color='skyblue', edgecolor='black')
    axes[0,0].set_xlabel('GFP Fluorescence')
    axes[0,0].set_ylabel('Count')
    axes[0,0].set_title('Intensity Distribution (All Sequences)')
    axes[0,0].set_yscale('log')
    
    # 1.2: Log intensity distribution
    axes[0,1].hist(all_valid_constructs['log_GFP'], bins=50, alpha=0.7, color='lightcoral', edgecolor='black')
    axes[0,1].set_xlabel('Log10(GFP Fluorescence)')
    axes[0,1].set_ylabel('Count')
    axes[0,1].set_title('Log10 Intensity Distribution')
    
    # 1.3: Sequence length distribution
    axes[0,2].hist(all_valid_constructs['sequence_length'], bins=20, alpha=0.7, color='lightgreen', edgecolor='black')
    axes[0,2].set_xlabel('Sequence Length (bp)')
    axes[0,2].set_ylabel('Count')
    axes[0,2].set_title('Sequence Length Distribution')
    
    # 1.4: PCA visualization (PC1 vs PC2) - with clustering colors
    scatter = axes[1,0].scatter(all_valid_constructs['PC1'], all_valid_constructs['PC2'], 
                               c=all_valid_constructs['PCA_Cluster'], cmap='tab10', alpha=0.6, s=20)
    axes[1,0].set_xlabel('PC1')
    axes[1,0].set_ylabel('PC2')
    axes[1,0].set_title('PCA Space (PC1 vs PC2) - Clustered')
    plt.colorbar(scatter, ax=axes[1,0], label='PCA Cluster')
    
    # 1.5: Intensity vs PC1 - with clustering colors
    scatter2 = axes[1,1].scatter(all_valid_constructs['PC1'], all_valid_constructs['GFP_fluorescence'], 
                                c=all_valid_constructs['PCA_Cluster'], cmap='tab10', alpha=0.6, s=20)
    axes[1,1].set_xlabel('PC1')
    axes[1,1].set_ylabel('GFP Fluorescence')
    axes[1,1].set_title('Intensity vs PC1 - Clustered')
    axes[1,1].set_yscale('log')
    plt.colorbar(scatter2, ax=axes[1,1], label='PCA Cluster')
    
    # 1.6: 2D Grid approach summary
    axes[1,2].text(0.1, 0.9, '2D Grid Approach Summary:', fontsize=14, fontweight='bold', transform=axes[1,2].transAxes)
    axes[1,2].text(0.1, 0.8, f'• Grid size: {N_BINS}×{N_BINS} = {N_BINS*N_BINS} cells', fontsize=12, transform=axes[1,2].transAxes)
    axes[1,2].text(0.1, 0.7, f'• Filled cells: {len(dist_pairs)}', fontsize=12, transform=axes[1,2].transAxes)
    axes[1,2].text(0.1, 0.6, f'• Empty cells backfilled: {N_BINS*N_BINS - len(dist_pairs)}', fontsize=12, transform=axes[1,2].transAxes)
    axes[1,2].text(0.1, 0.5, f'• Set 1 sequences: {len(final_set1)}', fontsize=12, transform=axes[1,2].transAxes)
    axes[1,2].text(0.1, 0.4, f'• Set 2 sequences: {len(final_set2)}', fontsize=12, transform=axes[1,2].transAxes)
    axes[1,2].text(0.1, 0.3, f'• Total combinations: {len(concat_df) if len(concat_df) > 0 else 0}', fontsize=12, transform=axes[1,2].transAxes)
    axes[1,2].set_xlim(0, 1)
    axes[1,2].set_ylim(0, 1)
    axes[1,2].axis('off')
    axes[1,2].set_title('2D Grid Selection Method')
    
    plt.tight_layout()
    pdf.savefig(fig, bbox_inches='tight')
    plt.close()
    
    # Page 2: Selected Sets Comparison
    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    fig.suptitle('Selected Sets Analysis', fontsize=16, fontweight='bold')
    
    # 2.1: Set intensity comparison
    axes[0,0].hist(final_set1['GFP_fluorescence'], bins=20, alpha=0.7, color='red', 
                  label=f'Set 1 (n={len(final_set1)})', edgecolor='black')
    axes[0,0].hist(final_set2['GFP_fluorescence'], bins=20, alpha=0.7, color='blue', 
                  label=f'Set 2 (n={len(final_set2)})', edgecolor='black')
    axes[0,0].set_xlabel('GFP Fluorescence')
    axes[0,0].set_ylabel('Count')
    axes[0,0].set_title('Intensity Distribution Comparison')
    axes[0,0].legend()
    axes[0,0].set_yscale('log')
    
    # 2.2: PCA space visualization of selected sets
    axes[0,1].scatter(final_set1['PC1'], final_set1['PC2'], c='red', alpha=0.7, s=50, 
                     label=f'Set 1 (2D Grid)', edgecolor='darkred')
    axes[0,1].scatter(final_set2['PC1'], final_set2['PC2'], c='blue', alpha=0.7, s=50, 
                     label=f'Set 2 (2D Grid)', edgecolor='darkblue')
    axes[0,1].set_xlabel('PC1')
    axes[0,1].set_ylabel('PC2')
    axes[0,1].set_title('Selected Sets in PCA Space')
    axes[0,1].legend()
    
    # 2.3: Sequence length comparison
    axes[0,2].hist(final_set1['sequence_length'], bins=15, alpha=0.7, color='red', 
                  label=f'Set 1', edgecolor='black')
    axes[0,2].hist(final_set2['sequence_length'], bins=15, alpha=0.7, color='blue', 
                  label=f'Set 2', edgecolor='black')
    axes[0,2].set_xlabel('Sequence Length (bp)')
    axes[0,2].set_ylabel('Count')
    axes[0,2].set_title('Sequence Length Comparison')
    axes[0,2].legend()
    
    # 2.4: Intensity vs sequence length
    axes[1,0].scatter(final_set1['sequence_length'], final_set1['GFP_fluorescence'], 
                     c='red', alpha=0.7, s=30, label='Set 1', edgecolor='darkred')
    axes[1,0].scatter(final_set2['sequence_length'], final_set2['GFP_fluorescence'], 
                     c='blue', alpha=0.7, s=30, label='Set 2', edgecolor='darkblue')
    axes[1,0].set_xlabel('Sequence Length (bp)')
    axes[1,0].set_ylabel('GFP Fluorescence')
    axes[1,0].set_title('Intensity vs Sequence Length')
    axes[1,0].set_yscale('log')
    axes[1,0].legend()
    
    # 2.5: Cumulative distribution of GFP intensity (log scale) for both sets
    # Sort sequences by intensity for cumulative plotting
    set1_sorted = final_set1.sort_values('GFP_fluorescence')
    set2_sorted = final_set2.sort_values('GFP_fluorescence')
    
    # Create cumulative counts
    set1_cumulative = np.arange(1, len(set1_sorted) + 1)
    set2_cumulative = np.arange(1, len(set2_sorted) + 1)
    
    # Plot cumulative distributions
    axes[1,1].plot(set1_sorted['GFP_fluorescence'], set1_cumulative, 
                   color='red', linewidth=2, label=f'Set 1 (n={len(final_set1)})', marker='o', markersize=4)
    axes[1,1].plot(set2_sorted['GFP_fluorescence'], set2_cumulative, 
                   color='blue', linewidth=2, label=f'Set 2 (n={len(final_set2)})', marker='s', markersize=4)
    
    axes[1,1].set_xlabel('GFP Fluorescence (log scale)')
    axes[1,1].set_ylabel('Cumulative Number of Sequences')
    axes[1,1].set_title('Cumulative Intensity Distribution')
    axes[1,1].set_xscale('log')
    axes[1,1].legend()
    axes[1,1].grid(True, alpha=0.3)
    
    # Add coverage information
    total_range = all_valid_constructs['GFP_fluorescence'].max() - all_valid_constructs['GFP_fluorescence'].min()
    set1_coverage = (final_set1['GFP_fluorescence'].max() - final_set1['GFP_fluorescence'].min()) / total_range * 100
    set2_coverage = (final_set2['GFP_fluorescence'].max() - final_set2['GFP_fluorescence'].min()) / total_range * 100
    
    # Add text box with coverage info
    coverage_text = f'Coverage:\nSet 1: {set1_coverage:.1f}%\nSet 2: {set2_coverage:.1f}%'
    axes[1,1].text(0.02, 0.98, coverage_text, transform=axes[1,1].transAxes, 
                   bbox=dict(boxstyle="round,pad=0.3", facecolor='white', alpha=0.8),
                   verticalalignment='top', fontsize=9)
    
    # 2.6: Summary statistics table
    axes[1,2].axis('off')
    summary_data = [
        ['Metric', 'Set 1', 'Set 2'],
        ['Number of sequences', f'{len(final_set1)}', f'{len(final_set2)}'],
        ['Selection method', '2D Grid', '2D Grid'],
        ['Min intensity', f'{final_set1["GFP_fluorescence"].min():.0f}', f'{final_set2["GFP_fluorescence"].min():.0f}'],
        ['Max intensity', f'{final_set1["GFP_fluorescence"].max():.0f}', f'{final_set2["GFP_fluorescence"].max():.0f}'],
        ['Mean intensity', f'{final_set1["GFP_fluorescence"].mean():.0f}', f'{final_set2["GFP_fluorescence"].mean():.0f}'],
        ['Seq length range', f'{final_set1["sequence_length"].min()}-{final_set1["sequence_length"].max()}bp', 
         f'{final_set2["sequence_length"].min()}-{final_set2["sequence_length"].max()}bp'],
        ['Inter-set PCA distance', f'{mean_inter_distance:.3f}', ''],
        ['Grid size', f'{N_BINS}×{N_BINS}', '']
    ]
    
    table = axes[1,2].table(cellText=summary_data, cellLoc='center', loc='center', 
                           colWidths=[0.6, 0.4, 0.4])
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1, 2)
    axes[1,2].set_title('Summary Statistics', pad=20)
    
    plt.tight_layout()
    pdf.savefig(fig, bbox_inches='tight')
    plt.close()
    
    # Page 3: Concatenated Sequences Analysis
    if len(concat_df) > 0:
        fig, axes = plt.subplots(1, 3, figsize=(18, 6))
        fig.suptitle('Concatenated Sequences Analysis', fontsize=16, fontweight='bold')
        
        # 3.1: Total length distribution
        axes[0].hist(concat_df['Total_Length'], bins=30, alpha=0.7, color='green', edgecolor='black')
        axes[0].set_xlabel('Total Sequence Length (bp)')
        axes[0].set_ylabel('Count')
        axes[0].set_title('Concatenated Sequence Length Distribution')
        axes[0].axvline(concat_df['Total_Length'].mean(), color='red', linestyle='--', 
                         label=f'Mean: {concat_df["Total_Length"].mean():.1f}bp')
        axes[0].legend()
        
        # 3.2: Set1 vs Set2 intensity combination map
        scatter = axes[1].scatter(concat_df['Set1_Intensity'], concat_df['Set2_Intensity'], 
                                   c=concat_df['Total_Length'], cmap='viridis', alpha=0.6, s=8)
        axes[1].set_xlabel('Set 1 Intensity')
        axes[1].set_ylabel('Set 2 Intensity')
        axes[1].set_title('Expression Combination Map')
        axes[1].set_xscale('log')
        axes[1].set_yscale('log')
        # Make the plot square
        axes[1].set_aspect('equal', adjustable='box')
        plt.colorbar(scatter, ax=axes[1], label='Total Length (bp)')
        
        # 3.3: Complete architecture component breakdown (all 7 components)
        components = [
            'Arch_Primer5_Length',
            'Arch_BsaI_Fwd_Length', 
            'Arch_Set1_RevComp_Length',
            'Arch_Terminator_Length',
            'Arch_Set2_Forward_Length',
            'Arch_BsaI_Rev_Length',
            'Arch_Primer3_Length'
        ]
        component_labels = ['Primer5', 'BsaI_Fwd', 'Set1_RevComp', 'Terminator', 'Set2_Forward', 'BsaI_Rev', 'Primer3']
        mean_lengths = [concat_df[comp].mean() for comp in components]
        colors = ['lightblue', 'orange', 'red', 'gray', 'blue', 'orange', 'lightblue']
        bars = axes[2].bar(component_labels, mean_lengths, color=colors, alpha=0.7, edgecolor='black')
        axes[2].set_ylabel('Average Length (bp)')
        axes[2].set_title('Complete Architecture Component Lengths')
        axes[2].tick_params(axis='x', rotation=45)
        
        # Add value labels on bars
        for bar, length in zip(bars, mean_lengths):
            axes[2].text(bar.get_x() + bar.get_width()/2, bar.get_height() + 1, 
                          f'{length:.1f}bp', ha='center', va='bottom', fontsize=8)
        
        # Add total verification text
        total_from_components = sum(mean_lengths)
        actual_total = concat_df['Total_Length'].mean()
        axes[2].text(0.5, 0.95, f'Sum: {total_from_components:.1f}bp\nTotal: {actual_total:.1f}bp', 
                      transform=axes[2].transAxes, ha='center', va='top',
                      bbox=dict(boxstyle="round,pad=0.3", facecolor='lightgray', alpha=0.7),
                      fontsize=8)
        
        # 3.4: Synthesis feasibility analysis (moved to text on axes[2])
        # Add synthesis statistics text to the bottom of the architecture chart
        
        # Calculate synthesis statistics
        max_synth_length = 300
        suitable_for_synthesis = len(concat_df[concat_df['Total_Length'] <= max_synth_length])
        synthesis_rate = (suitable_for_synthesis / len(concat_df)) * 100
        
        # Add key synthesis stats as text overlay on architecture chart
        synthesis_text = (f'Synthesis Stats:\n'
                         f'Total: {len(concat_df):,} combinations\n'
                         f'Avg length: {concat_df["Total_Length"].mean():.0f}bp\n'
                         f'Suitable (≤300bp): {synthesis_rate:.0f}%')
        axes[2].text(0.02, 0.85, synthesis_text, transform=axes[2].transAxes, 
                     bbox=dict(boxstyle="round,pad=0.3", facecolor='lightyellow', alpha=0.8),
                     verticalalignment='top', fontsize=8)
        
        plt.tight_layout()
        pdf.savefig(fig, bbox_inches='tight')
        plt.close()
    
    # Page 4: 2D Grid Selection Details
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    fig.suptitle('2D Grid Selection Method Details', fontsize=16, fontweight='bold')
    
    # 4.1: Grid visualization (intensity vs PC1) with cluster colors
    # Create a grid showing the distribution of sequences with cluster colors
    scatter = axes[0].scatter(all_valid_constructs['PC1'], all_valid_constructs['log_GFP'], 
                               c=all_valid_constructs['PCA_Cluster'], cmap='tab10', alpha=0.3, s=10)
    axes[0].scatter(final_set1['PC1'], final_set1['log_GFP'], 
                     c='red', s=50, alpha=0.8, label='Set 1 selected', edgecolor='darkred')
    axes[0].scatter(final_set2['PC1'], final_set2['log_GFP'], 
                     c='blue', s=50, alpha=0.8, label='Set 2 selected', edgecolor='darkblue')
    axes[0].set_xlabel('PC1')
    axes[0].set_ylabel('Log10(GFP Fluorescence)')
    axes[0].set_title('2D Grid Selection in PC1-Intensity Space (Clustered)')
    axes[0].legend()
    axes[0].grid(True, alpha=0.3)
    plt.colorbar(scatter, ax=axes[0], label='PCA Cluster')
    
    # 4.2: Set comparison statistics
    axes[1].axis('off')
    comparison_stats = [
        ['Metric', 'Set 1', 'Set 2'],
        ['Number of sequences', f'{len(final_set1)}', f'{len(final_set2)}'],
        ['Min intensity', f'{final_set1["GFP_fluorescence"].min():.0f}', f'{final_set2["GFP_fluorescence"].min():.0f}'],
        ['Max intensity', f'{final_set1["GFP_fluorescence"].max():.0f}', f'{final_set2["GFP_fluorescence"].max():.0f}'],
        ['Mean intensity', f'{final_set1["GFP_fluorescence"].mean():.0f}', f'{final_set2["GFP_fluorescence"].mean():.0f}'],
        ['Std intensity', f'{final_set1["GFP_fluorescence"].std():.0f}', f'{final_set2["GFP_fluorescence"].std():.0f}'],
        ['Min PC1', f'{final_set1["PC1"].min():.3f}', f'{final_set2["PC1"].min():.3f}'],
        ['Max PC1', f'{final_set1["PC1"].max():.3f}', f'{final_set2["PC1"].max():.3f}'],
        ['Mean PC1', f'{final_set1["PC1"].mean():.3f}', f'{final_set2["PC1"].mean():.3f}'],
    ]
    
    table = axes[1].table(cellText=comparison_stats, cellLoc='center', loc='center',
                           colWidths=[0.5, 0.25, 0.25])
    table.auto_set_font_size(False)
    table.set_fontsize(9)
    table.scale(1, 2)
    axes[1].set_title('Set Comparison Statistics', pad=20)
    
    # 4.3: Selection method summary
    axes[2].axis('off')
    method_summary = [
        ['2D Grid Method', 'Value'],
        ['Grid dimensions', f'{N_BINS}×{N_BINS}'],
        ['Total cells', f'{N_BINS*N_BINS}'],
        ['Filled cells', f'{len(dist_pairs)}'],
        ['Empty cells backfilled', f'{N_BINS*N_BINS - len(dist_pairs)}'],
        ['PC1 range', f'{all_valid_constructs["PC1"].min():.2f} to {all_valid_constructs["PC1"].max():.2f}'],
        ['Intensity range (log)', f'{all_valid_constructs["log_GFP"].min():.2f} to {all_valid_constructs["log_GFP"].max():.2f}'],
        ['Mean inter-set distance', f'{mean_inter_distance:.3f}'],
        ['Set 1 diversity', f'{set1_range/total_range*100:.1f}% coverage'],
        ['Set 2 diversity', f'{set2_range/total_range*100:.1f}% coverage']
    ]
    
    table = axes[2].table(cellText=method_summary, cellLoc='center', loc='center',
                           colWidths=[0.6, 0.4])
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1, 2)
    axes[2].set_title('2D Grid Method Summary', pad=20)
    
    plt.tight_layout()
    pdf.savefig(fig, bbox_inches='tight')
    plt.close()

print(f"Comprehensive PDF report saved: {pdf_filename}")

# Save individual PNG files for easy viewing
print("Saving individual PNG files...")

# Overview plot
fig, axes = plt.subplots(2, 2, figsize=(15, 12))
fig.suptitle('Dual Inverted Promoter Library - Key Results', fontsize=16, fontweight='bold')

# PCA overview with clusters and selected sets
scatter = axes[0,0].scatter(all_valid_constructs['PC1'], all_valid_constructs['PC2'], 
                           c=all_valid_constructs['PCA_Cluster'], cmap='tab10', alpha=0.3, s=20)
axes[0,0].scatter(final_set1['PC1'], final_set1['PC2'], c='red', s=100, marker='s', 
                 edgecolor='darkred', alpha=0.8, label=f'Set 1 (n={len(final_set1)})')
axes[0,0].scatter(final_set2['PC1'], final_set2['PC2'], c='blue', s=100, marker='s', 
                 edgecolor='darkblue', alpha=0.8, label=f'Set 2 (n={len(final_set2)})')
axes[0,0].set_xlabel('PC1')
axes[0,0].set_ylabel('PC2')
axes[0,0].set_title('PCA Space with Clusters and Selected Sets')
axes[0,0].legend()

# Intensity comparison
axes[0,1].hist(final_set1['GFP_fluorescence'], bins=20, alpha=0.7, color='red', 
              label=f'Set 1', edgecolor='darkred', density=True)
axes[0,1].hist(final_set2['GFP_fluorescence'], bins=20, alpha=0.7, color='blue', 
              label=f'Set 2', edgecolor='darkblue', density=True)
axes[0,1].set_xlabel('GFP Fluorescence')
axes[0,1].set_ylabel('Density')
axes[0,1].set_title('Intensity Distribution Comparison')
axes[0,1].set_xscale('log')
axes[0,1].legend()

# Expression combination map
if len(concat_df) > 0:
    scatter = axes[1,0].scatter(concat_df['Set1_Intensity'], concat_df['Set2_Intensity'], 
                               c=concat_df['Total_Length'], cmap='viridis', alpha=0.6, s=8)
    axes[1,0].set_xlabel('Set 1 Intensity')
    axes[1,0].set_ylabel('Set 2 Intensity')
    axes[1,0].set_title('Expression Combination Map')
    axes[1,0].set_xscale('log')
    axes[1,0].set_yscale('log')
    # Make the plot square
    axes[1,0].set_aspect('equal', adjustable='box')
    plt.colorbar(scatter, ax=axes[1,0], label='Total Length (bp)')

# Length distribution
if len(concat_df) > 0:
    axes[1,1].hist(concat_df['Total_Length'], bins=30, alpha=0.7, color='green', edgecolor='black')
    axes[1,1].set_xlabel('Total Sequence Length (bp)')
    axes[1,1].set_ylabel('Count')
    axes[1,1].set_title('Concatenated Sequence Lengths')
    axes[1,1].axvline(concat_df['Total_Length'].mean(), color='red', linestyle='--', 
                     label=f'Mean: {concat_df["Total_Length"].mean():.1f}bp')
    axes[1,1].legend()

plt.tight_layout()
plt.savefig(os.path.join(output_dir, "06_key_results_overview.png"), dpi=300, bbox_inches='tight')
plt.close()

print("Individual PNG files saved")
print(f"\n=== VISUALIZATION COMPLETE ===")
print(f"Files generated:")
print(f"  - {pdf_filename}")
print(f"  - 06_key_results_overview.png")
print(f"  - 04_intensity_distribution_analysis.png (from earlier)")

print(f"\n=== FINAL SUMMARY ===")
print(f"✅ Analysis completed successfully!")
print(f"📁 Results directory: {output_dir}")
print(f"📊 Total sequences analyzed: {len(all_valid_constructs)}")
print(f"🔬 2D Grid approach: {N_BINS}×{N_BINS} grid with {len(dist_pairs)} filled cells")
print(f"📋 Set 1: {len(final_set1)} sequences")
print(f"📋 Set 2: {len(final_set2)} sequences")
print(f"🧬 Total combinations: {len(concat_df) if len(concat_df) > 0 else 0:,}")
print(f"📏 Sequence lengths: {concat_df['Total_Length'].min() if len(concat_df) > 0 else 'N/A'}-{concat_df['Total_Length'].max() if len(concat_df) > 0 else 'N/A'}bp")
print(f"✅ All sequences suitable for synthesis (≤300bp)")
print(f"📈 Mean inter-set PCA distance: {mean_inter_distance:.3f}")
print(f"🎯 Fluorescence coverage: Set1={set1_range/total_range*100:.1f}%, Set2={set2_range/total_range*100:.1f}%")
