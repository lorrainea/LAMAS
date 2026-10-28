from collections import Counter
from itertools import combinations
import math
import os
import subprocess
import tempfile

import numpy as np
from Bio import SeqIO
from sklearn.metrics import silhouette_score, adjusted_rand_score
from sklearn.manifold import MDS


# ======================================================================
# 1. SEQUENCE LOADING & DISTANCE MATH
# ======================================================================

def acs_xx(x: str) -> float:
    n = len(x)
    return sum(range(1, n + 1)) / n


def acs_xy(li: list[int], y: str) -> float:
    return sum(li) / len(y)


def nacs(li: list[int], x: str, y: str, d: int) -> float:
    acs_xy_var = acs_xy(li, y)
    acs_xx_var = acs_xx(x)
        
    return (math.log2(d) / acs_xy_var) - (math.log2(d) / acs_xx_var)


def dis(li_xy: list[int], li_yx: list[int], x: str, y: str, d_xy: int, d_yx: int) -> float:
    return (nacs(li_xy, x, y, d_xy) + nacs(li_yx, y, x, d_yx)) / 2.0


def load_fasta(filepath: str) -> dict[str, str]:
    """
    Reads a FASTA file and returns: {taxon_id: sequence_string}
    """
    sequences = {}
    for record in SeqIO.parse(filepath, "fasta"):
        sequences[record.id] = str(record.seq).upper()
    return sequences


# ======================================================================
# 2. HELPERS: PARSE OUTPUT VALUES & DISTANCE MATRIX HELPER
# ======================================================================

def parse_first_numeric_value(content: str, method_name: str) -> float:
    for token in content.replace("\n", " ").split():
        try:
            return float(token)
        except ValueError:
            continue
    raise ValueError(f"Could not find numeric {method_name} distance.")


def read_integer_list_from_file(file_path: str) -> list[int]:
    values = []
    if os.path.exists(file_path):
        with open(file_path, "r") as file:
            for line in file:
                line = line.strip()
                if line:
                    values.append(int(line))
    return values


def build_condensed_matrix(taxa: list[str], matrix_dict: dict[tuple[str, str], float]) -> np.ndarray:
    """
    Converts pairwise tuple dicts into a square 2D symmetric numpy distance matrix.
    """
    n = len(taxa)
    dist_matrix = np.zeros((n, n), dtype=float)

    for i, taxon1 in enumerate(taxa):
        for j in range(i + 1, n):
            taxon2 = taxa[j]
            if (taxon1, taxon2) in matrix_dict:
                d_val = matrix_dict[(taxon1, taxon2)]
            elif (taxon2, taxon1) in matrix_dict:
                d_val = matrix_dict[(taxon2, taxon1)]
            else:
                raise KeyError(f"Missing distance for pair: {taxon1}, {taxon2}")

            dist_matrix[i, j] = d_val
            dist_matrix[j, i] = d_val

    return dist_matrix


def save_and_print_matrix(matrix: np.ndarray, taxa: list[str], method_name: str, file_prefix: str, preview_size: int = 5):
    """
    Saves the full distance matrix as a CSV file and prints a preview to the terminal.
    """
    filename = f"{file_prefix}.matrix_{method_name}.csv"
    
    # Save matrix to CSV
    header = ",".join(taxa)
    np.savetxt(filename, matrix, delimiter=",", header=header, comments="")
    print(f"Saved full distance matrix to: {filename}")

    # Print preview matrix to console
    print(f"\n--- Matrix Preview for {method_name} ({preview_size}x{preview_size}) ---")
    sub_taxa = taxa[:preview_size]
    sub_matrix = matrix[:preview_size, :preview_size]

    # Column headers
    header_str = f"{'':<12}" + "".join([f"{t:>12}" for t in sub_taxa])
    print(header_str)
    print("-" * len(header_str))

    # Matrix rows
    for i, taxon in enumerate(sub_taxa):
        row_str = f"{taxon:<12}" + "".join([f"{sub_matrix[i, j]:>12.4f}" for j in range(len(sub_taxa))])
        print(row_str)
    print()


# ======================================================================
# 3. EDLIB EDIT DISTANCE
# ======================================================================

def run_edlib(seq1: str, seq2: str, binary_path: str = "./ed") -> float:
    with tempfile.NamedTemporaryFile(mode="w+", delete=False) as f1, \
         tempfile.NamedTemporaryFile(mode="w+", delete=False) as f2, \
         tempfile.NamedTemporaryFile(mode="w+", delete=False) as f_out:

        f1.write(f">seq1\n{seq1}\n")
        f2.write(f">seq2\n{seq2}\n")

        f1.flush()
        f2.flush()

        cmd = [binary_path, "-x", f1.name, "-y", f2.name, "-o", f_out.name]

        try:
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            f_out.seek(0)
            content = f_out.read().strip()
            if not content:
                content = result.stdout.strip()

            return parse_first_numeric_value(content, "Edlib")

        except Exception as e:
            print(f"\nEdlib error: {e}")
            return 0.0

        finally:
            for path in [f1.name, f2.name, f_out.name]:
                if os.path.exists(path):
                    os.remove(path)


# ======================================================================
# 4. LAMAS NORMALIZED DISTANCE
# ======================================================================

def run_lamas_normalized(
    seq_x: str,
    seq_y: str,
    d: int,
    binary_path: str = "./lamas_linear"
) -> float:
    with tempfile.NamedTemporaryFile(mode="w+", delete=False) as f_x, \
         tempfile.NamedTemporaryFile(mode="w+", delete=False) as f_y, \
         tempfile.NamedTemporaryFile(mode="w+", delete=False) as out_xy, \
         tempfile.NamedTemporaryFile(mode="w+", delete=False) as out_yx:

        f_x.write(f">seq_x\n{seq_x}\n")
        f_y.write(f">seq_y\n{seq_y}\n")
        f_x.flush()
        f_y.flush()

        try:
            cmd_xy = [binary_path, "-x", f_x.name, "-y", f_y.name, "-d", str(d), "-o", out_xy.name]
            subprocess.run(cmd_xy, check=True, capture_output=True)
            lamas_xy = read_integer_list_from_file(out_xy.name)

            cmd_yx = [binary_path, "-x", f_y.name, "-y", f_x.name, "-d", str(d), "-o", out_yx.name]
            subprocess.run(cmd_yx, check=True, capture_output=True)
            lamas_yx = read_integer_list_from_file(out_yx.name)

            return dis(lamas_xy, lamas_yx, seq_x, seq_y, d, d)

        except Exception as e:
            print(f"\nLAMAS error: {e}")
            return 0.0

        finally:
            for path in [f_x.name, f_y.name, out_xy.name, out_yx.name]:
                if os.path.exists(path):
                    os.remove(path)


# ======================================================================
# 5. MS NORMALIZED DISTANCE
# ======================================================================

def run_ms_normalized(
    seq_x: str,
    seq_y: str,
    binary_path: str = "./ms"
) -> float:
    with tempfile.NamedTemporaryFile(mode="w+", delete=False) as f_x, \
         tempfile.NamedTemporaryFile(mode="w+", delete=False) as f_y, \
         tempfile.NamedTemporaryFile(mode="w+", delete=False) as out_xy, \
         tempfile.NamedTemporaryFile(mode="w+", delete=False) as out_yx:

        f_x.write(f">seq_x\n{seq_x}\n")
        f_y.write(f">seq_y\n{seq_y}\n")
        f_x.flush()
        f_y.flush()

        try:
            ms_xy_d = len(seq_y) - 1
            cmd_xy = [binary_path, "-x", f_x.name, "-y", f_y.name, "-o", out_xy.name]
            subprocess.run(cmd_xy, check=True, capture_output=True)
            ms_xy = read_integer_list_from_file(out_xy.name)

            ms_yx_d = len(seq_x) - 1
            cmd_yx = [binary_path, "-x", f_y.name, "-y", f_x.name, "-o", out_yx.name]
            subprocess.run(cmd_yx, check=True, capture_output=True)
            ms_yx = read_integer_list_from_file(out_yx.name)

            return dis(ms_xy, ms_yx, seq_x, seq_y, ms_xy_d, ms_yx_d)

        except Exception as e:
            print(f"\nMS-naive error: {e}")
            return 0.0

        finally:
            for path in [f_x.name, f_y.name, out_xy.name, out_yx.name]:
                if os.path.exists(path):
                    os.remove(path)


# ======================================================================
# 6. Q-GRAM DISTANCE
# ======================================================================

def qgram_distance(s1: str, s2: str, q: int) -> float:
    if len(s1) < q or len(s2) < q:
        return float(abs(len(s1) - len(s2)))

    qgrams1 = Counter(s1[i:i + q] for i in range(len(s1) - q + 1))
    qgrams2 = Counter(s2[i:i + q] for i in range(len(s2) - q + 1))
    all_keys = set(qgrams1.keys()).union(qgrams2.keys())

    return float(sum(abs(qgrams1[k] - qgrams2[k]) for k in all_keys))


# ======================================================================
# 7. K-MEDIANS & CLUSTERING EVALUATION
# ======================================================================

def run_k_medians_clustering(
    dist_matrix: np.ndarray, 
    n_clusters: int = 10, 
    max_iter: int = 500, 
    random_state: int = 42
) -> tuple[np.ndarray, float]:
    """
    Runs true K-Medians clustering on a precomputed distance matrix by:
    1. Embedding the distance matrix into a continuous coordinate space using MDS.
    2. Applying K-Medians (Manhattan L1 distance & coordinate-wise medians).
    3. Calculating the Silhouette Score on the original distance matrix.
    """
    # Embed pairwise distance matrix into coordinate space (e.g., 10 dimensions)
    mds = MDS(n_components=10, dissimilarity="precomputed", random_state=random_state, n_init=4)
    coords = mds.fit_transform(dist_matrix)
    
    np.random.seed(random_state)
    n_samples = coords.shape[0]

    # Initialize random medians in coordinate space
    initial_indices = np.random.choice(n_samples, n_clusters, replace=False)
    medians = coords[initial_indices].copy()

    labels = np.zeros(n_samples, dtype=int)

    for _ in range(max_iter):
        # Assign points to the nearest median using Manhattan (L1) distance
        distances = np.zeros((n_samples, n_clusters))
        for k in range(n_clusters):
            distances[:, k] = np.sum(np.abs(coords - medians[k]), axis=1)
        
        labels = np.argmin(distances, axis=1)

        new_medians = medians.copy()
        for k in range(n_clusters):
            cluster_members = np.where(labels == k)[0]
            if len(cluster_members) > 0:
                # Update center using coordinate-wise median (True K-Medians)
                new_medians[k] = np.median(coords[cluster_members], axis=0)

        if np.allclose(medians, new_medians):
            break
        medians = new_medians

    # Compute Silhouette score using your original true distance matrix
    sil_score = silhouette_score(dist_matrix, labels, metric="precomputed")
    
    return labels, float(sil_score)


# ======================================================================
# 8. PROCESS SINGLE FASTA DATASET
# ======================================================================

def process_file(
    filepath: str,
    d_range: list[int],
    q_range: list[int],
    n_clusters: int = 10,
    edlib_binary: str = "./ed",
    lamas_binary: str = "./lamas_linear",
    ms_binary: str = "./ms"
):
    if not os.path.exists(filepath):
        print(f"ERROR: File '{filepath}' does not exist.")
        return

    output_prefix = os.path.basename(filepath)

    print("=" * 75)
    print(f"Processing file: {filepath}")
    print("=" * 75)

    sequences = load_fasta(filepath)
    taxa = list(sequences.keys())
    n_seqs = len(taxa)

    print(f"Number of sequences: {n_seqs}")
    if n_seqs < n_clusters:
        print(f"ERROR: Number of sequences ({n_seqs}) must be >= n_clusters ({n_clusters}).")
        return

    # Construct ground-truth cluster labels (0-9 for 1-10, 10-19, etc.)
    true_labels = np.array([i // 10 for i in range(n_seqs)])

    pairs = list(combinations(taxa, 2))
    print(f"Number of pairwise comparisons: {len(pairs)}\n")

    matrices = {}
    cluster_results = {}

    # 1. EDLIB EDIT DISTANCE
    print("=" * 75)
    print("1. EDIT DISTANCE - EDLIB")
    print("=" * 75)

    edit_dists = {}
    for count, pair in enumerate(pairs, start=1):
        edit_dists[pair] = run_edlib(
            sequences[pair[0]],
            sequences[pair[1]],
            binary_path=edlib_binary
        )
        if count % 100 == 0 or count == len(pairs):
            print(f"  Edlib: {count}/{len(pairs)}")

    matrices["Edit_Distance"] = build_condensed_matrix(taxa, edit_dists)
    save_and_print_matrix(matrices["Edit_Distance"], taxa, "Edit_Distance", output_prefix)

    # 2. LAMAS NORMALIZED
    print("\n" + "=" * 75)
    print("2. LAMAS NORMALIZED")
    print("=" * 75)

    for d in d_range:
        label = f"LAMAS_d{d}"
        print(f"\nComputing LAMAS distances, d={d}")
        lamas_dists = {}

        for count, pair in enumerate(pairs, start=1):
            lamas_dists[pair] = run_lamas_normalized(
                sequences[pair[0]],
                sequences[pair[1]],
                d=d,
                binary_path=lamas_binary
            )
            if count % 100 == 0 or count == len(pairs):
                print(f"  d={d}: {count}/{len(pairs)}")

        matrices[label] = build_condensed_matrix(taxa, lamas_dists)
        save_and_print_matrix(matrices[label], taxa, label, output_prefix)

    # 3. Q-GRAM DISTANCE
    print("\n" + "=" * 75)
    print("3. Q-GRAM DISTANCE")
    print("=" * 75)

    for q in q_range:
        label = f"QGram_q{q}"
        print(f"\nComputing q-gram distances, q={q}")
        qgram_dists = {}

        for count, pair in enumerate(pairs, start=1):
            qgram_dists[pair] = qgram_distance(
                sequences[pair[0]],
                sequences[pair[1]],
                q=q
         )
            if count % 500 == 0 or count == len(pairs):
                print(f"  q={q}: {count}/{len(pairs)}")

        matrices[label] = build_condensed_matrix(taxa, qgram_dists)
        save_and_print_matrix(matrices[label], taxa, label, output_prefix)

    # 4. MS NORMALIZED
    print("\n" + "=" * 75)
    print("4. MS NORMALIZED")
    print("=" * 75)

    ms_dists = {}
    for count, pair in enumerate(pairs, start=1):
       ms_dists[pair] = run_ms_normalized(
            sequences[pair[0]],
            sequences[pair[1]],
            binary_path=ms_binary
        )
    if count % 5 == 0 or count == len(pairs):
        print(f"  MS-naive: {count}/{len(pairs)}")

    matrices["MS"] = build_condensed_matrix(taxa, ms_dists)
    save_and_print_matrix(matrices["MS"], taxa, "MS", output_prefix)

    # 5. K-MEDIANS & CLUSTERING EVALUATION
    print("\n" + "=" * 75)
    print(f"PERFORMING K-MEDIANS CLUSTERING (k={n_clusters}) & EVALUATION AGAINST TRUE CLUSTERS")
    print("=" * 75)

    for method, matrix in matrices.items():
        labels, sil_score = run_k_medians_clustering(matrix, n_clusters=n_clusters)
        ari_vs_true = adjusted_rand_score(true_labels, labels)
        
        cluster_results[method] = {
            "labels": labels,
            "silhouette": sil_score,
            "ari_vs_true": ari_vs_true
        }

    # 6. SAVE RESULTS
    scores_filename = f"{output_prefix}.clustering_scores.txt"
    with open(scores_filename, "w") as f:
        f.write("Method\tSilhouette_Score\tARI_vs_True_Clusters\n")
        for method, result in cluster_results.items():
            f.write(f"{method}\t{result['silhouette']:.6f}\t{result['ari_vs_true']:.6f}\n")

    print(f"Saved summary metrics to {scores_filename}")

    # 7. PRINT FINAL SUMMARY
    print("\n" + "=" * 75)
    print("FINAL CLUSTERING EVALUATION SUMMARY")
    print("=" * 75)
    print(f"{'Method':<25} | {'Silhouette Score':>18} | {'ARI vs True':>15}")
    print("-" * 65)
    for method, result in cluster_results.items():
        print(f"{method:<25} | {result['silhouette']:>18.6f} | {result['ari_vs_true']:>15.6f}")

    print("\nAnalysis complete.")


# ======================================================================
# 9. MAIN
# ======================================================================

if __name__ == "__main__":
    FILEPATH = "../10000.100.15"

    D_RANGE = [100, 500, 1000, 1500, 2000, 2500, 3000]
    Q_RANGE = [4, 6, 8, 10, 12,  14,  16]

    EDLIB_BINARY = "./ed"
    LAMAS_BINARY = "./lamas_linear"
    MS_BINARY = "./ms"

    # Set desired number of clusters (k)
    N_CLUSTERS = 10

    process_file(
        FILEPATH,
        d_range=D_RANGE,
        q_range=Q_RANGE,
        n_clusters=N_CLUSTERS,
        edlib_binary=EDLIB_BINARY,
        lamas_binary=LAMAS_BINARY,
        ms_binary=MS_BINARY
    )
