#!/usr/bin/env python3
import os
import subprocess
import re
import math
import tempfile
import sys
from pathlib import Path

# List of algorithms to run (make sure their binaries are compiled and in PATH or current dir)
ALGORITHMS = {
    "qgram": "../qgram",
    "ed": "../ed",                     
    "lamas_linear": "../lamas_linear",
    "lamas": "../lamas",          
    "ms": "../ms",
    "lamas_improved": "../lamas_improved",                
}

def parse_time_output(stderr_output):
    """
    Parses the output of /usr/bin/time -v
    Extracts Elapsed Time and Maximum Resident Set Size (KB), returning time in milliseconds.
    """
    elapsed_seconds = 0.0
    max_rss_kb = 0

    # More flexible regex for GNU time -v elapsed time
    match_time = re.search(r'Elapsed[^\n]*?:\s+([0-9:.]+)', stderr_output)
    if match_time:
        t_str = match_time.group(1).strip()
        parts = t_str.split(':')
        try:
            if len(parts) == 3:
                elapsed_seconds = float(parts[0]) * 3600 + float(parts[1]) * 60 + float(parts[2])
            elif len(parts) == 2:
                elapsed_seconds = float(parts[0]) * 60 + float(parts[1])
            else:
                elapsed_seconds = float(parts[0])
        except ValueError:
            elapsed_seconds = 0.0

    # Convert seconds to milliseconds
    elapsed_ms = elapsed_seconds * 1000.0

    # Parse Maximum Resident Set Size: e.g., "Maximum resident set size (kbytes): 4212"
    match_mem = re.search(r'Maximum resident set size[^\n]*?:\s+(\d+)', stderr_output, re.IGNORECASE)
    if match_mem:
        try:
            max_rss_kb = int(match_mem.group(1))
        except ValueError:
            max_rss_kb = 0

    return elapsed_ms, max_rss_kb

def parse_seq_file(file_path):
    pairs = []
    current_x = None
    
    with open(file_path, 'r', encoding='utf-8', errors='ignore') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                current_x = line[1:].strip()
            elif line.startswith('<') and current_x is not None:
                current_y = line[1:].strip()
                pairs.append((current_x, current_y))
                current_x = None
    return pairs

def compute_stats(values):
    if not values:
        return 0.0, 0.0, 0.0, 0.0
    n = len(values)
    mean = sum(values) / n
    minimum = min(values)
    maximum = max(values)
    variance = sum((v - mean) ** 2 for v in values) / n
    stddev = math.sqrt(variance)
    return mean, minimum, maximum, stddev

def run_benchmark(directory_path):
    dir_path = Path(directory_path)
    if not dir_path.is_dir():
        print(f"[ERROR] {directory_path} is not a valid directory.")
        return

    # Verify binaries before starting
    print("[INFO] Verifying binaries...")
    for algo, binary in ALGORITHMS.items():
        if not os.path.exists(binary):
            print(f"[WARNING] Binary for '{algo}' not found at path: {binary}")
        else:
            print(f"[OK] Found {algo} binary at {binary}")

    seq_files = sorted(list(dir_path.glob("*.seq")))
    if not seq_files:
        print(f"[ERROR] No .seq files found in directory: {directory_path}")
        return

    results = {algo: {'time': [], 'mem': []} for algo in ALGORITHMS}
    seq_lengths = []
    total_pairs = 0

    print(f"\n[INFO] Found {len(seq_files)} .seq files. Starting step-by-step benchmark...\n")

    for file_idx, file_path in enumerate(seq_files, 1):
        print(f"--- Processing File {file_idx}/{len(seq_files)}: {file_path.name} ---")
        pairs = parse_seq_file(file_path)
        if not pairs:
            print(f"    [WARNING] No sequence pairs found in {file_path.name}")
            continue

        for pair_idx, (seq_x, seq_y) in enumerate(pairs, 1):
            len_x = len(seq_x)
            len_y = len(seq_y)
            
            total_pairs += 1
            
            # Record individual sequence lengths for overall statistics
            seq_lengths.append(len_x)
            seq_lengths.append(len_y)
            
            d_val = max(1, int(0.1 * min(len_x, len_y)))
            q_val = 4

            # Write sequences to temporary files
            with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.txt') as fx, \
                 tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.txt') as fy, \
                 tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.txt') as fout:
                
                fx.write(seq_x)
                fx_name = fx.name
                fy.write(seq_y)
                fy_name = fy.name
                fout_name = fout.name

            try:
                for algo, binary in ALGORITHMS.items():
                    if not os.path.exists(binary):
                        continue

                    if algo == 'qgram':
                        cmd = f"/usr/bin/time -v {binary} -x {fx_name} -y {fy_name} -o {fout_name} -q {q_val}"
                    elif algo == 'ed':
                        cmd = f"/usr/bin/time -v {binary} -x {fx_name} -y {fy_name} -o {fout_name}"
                    elif algo in ['lamas_linear', 'lamas', 'lamas_improved']:
                        cmd = f"/usr/bin/time -v {binary} -x {fx_name} -y {fy_name} -o {fout_name} -d {d_val}"
                    else: 
                        cmd = f"/usr/bin/time -v {binary} -x {fx_name} -y {fy_name} -o {fout_name}"

                    # Run process
                    process = subprocess.run(cmd, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE, text=True)
                    
                    if process.returncode != 0:
                        print(f"    [FAIL] Pair {pair_idx} | Algo: {algo} exited with code {process.returncode}")
                        print(f"        Stderr: {process.stderr.strip()}")
                        continue

                    elapsed, mem = parse_time_output(process.stderr)
                    results[algo]['time'].append(elapsed)
                    results[algo]['mem'].append(mem)
                    print(f"    [SUCCESS] Pair {pair_idx} | {algo:<14} -> Time: {elapsed:.2f} ms | Mem: {mem} KB")

            finally:
                for fname in [fx_name, fy_name, fout_name]:
                    if os.path.exists(fname):
                        os.remove(fname)

    print(f"\n[INFO] Finished processing {total_pairs} total sequence pairs across {len(seq_files)} files.\n")

    # Compute overall sequence length statistics
    avg_seq_len, min_seq_len, max_seq_len, std_seq_len = compute_stats(seq_lengths)

    # Print Overall Summary Report
    print("="*85)
    print(f"{'OVERALL BENCHMARK RESULTS SUMMARY':^85}")
    print("="*85)
    print(f"Dataset Overview:")
    print(f"  Total Sequence Pairs    : {total_pairs}")
    print(f"  Total Sequences         : {len(seq_lengths)}")
    print(f"  Average Sequence Length: {avg_seq_len:.2f} chars (Min: {min_seq_len}, Max: {max_seq_len}, Std Dev: {std_seq_len:.2f})")
    print("-" * 85)

    for algo in ALGORITHMS:
        times = results[algo]['time']
        mems = results[algo]['mem']
        
        if not times:
            print(f"\nAlgorithm: {algo.upper()} -> No data collected (binary missing or all runs failed).")
            continue

        t_mean, t_min, t_max, t_std = compute_stats(times)
        m_mean, m_min, m_max, m_std = compute_stats(mems)

        print(f"\nAlgorithm: {algo.upper()}")
        print(f"  --- Execution Time (milliseconds) ---")
        print(f"      Average : {t_mean:.2f} ms")
        print(f"      Min     : {t_min:.2f} ms")
        print(f"      Max     : {t_max:.2f} ms")
        print(f"      Std Dev : {t_std:.2f} ms")
        print(f"  --- Peak Memory (KB / Max RSS) ---")
        print(f"      Average : {m_mean:.2f} KB")
        print(f"      Min     : {m_min} KB")
        print(f"      Max     : {m_max} KB")
        print(f"      Std Dev : {m_std:.2f} KB")
    print("="*85)

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Usage: python benchmark_debug.py <directory_path>")
        sys.exit(1)
    
    dir_path = sys.argv[1]
    run_benchmark(dir_path)
