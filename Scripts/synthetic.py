#!/usr/bin/env python3

import math
import subprocess
import matplotlib.pyplot as plt
from scipy.stats import pearsonr, spearmanr, kendalltau


##############################################################################
# ACS / DIS FRAMEWORK
##############################################################################

def acs_xx(x):
    n = len(x)
    return sum(range(1, n + 1)) / n

def acs_xy(li, y):
    return sum(li) / len(y)
    
def nacs(li, x, y, d):
    acs_xy_var = acs_xy(li, y)
    acs_xx_var = acs_xx(x)

    return (math.log2(d) / acs_xy_var) - (math.log2(d) / acs_xx_var)

def dis(li_xy, li_yx, x, y, d_xy, d_yx):
    return ( nacs(li_xy, x, y, d_xy) + nacs(li_yx, y, x, d_yx)) / 2


##############################################################################
# IO HELPERS
##############################################################################

def write_file(name, content):
    with open(name, "w") as f:
        f.write(content)

##############################################################################
# ED
##############################################################################

def compute_ed(pairs):
    ed = []
    for x, y in pairs:
        write_file("x_temp", str(x))
        write_file("y_temp", str(y))

        subprocess.run(["../ed", "-x", "x_temp", "-y", "y_temp", "-o", "out_ed"], check=True)
        with open("out_ed") as f:
            ed.append( int(f.readline()))
            
    return ed


##############################################################################
# LAMAS
##############################################################################

def evaluate_lamas(pairs, ed, d_min=1000, d_max=3000):

    results = {}
  
    for d in range(d_min, d_max + 1, 500):

        vals = []

        for x, y in pairs:
            write_file("x_temp", str(x))
            write_file("y_temp", str(y))

            subprocess.run(["../lamas", "-x", "x_temp", "-y", "y_temp", "-d", str(d), "-o", "out_xy"], check=True)
            subprocess.run(["../lamas", "-x", "y_temp", "-y", "x_temp", "-d", str(d), "-o", "out_yx"], check=True)

            with open("out_xy") as f:
                li_xy = [int(line.strip()) for line in f if line.strip()]

            with open("out_yx") as f:
                li_yx = [int(line.strip()) for line in f if line.strip()]

            vals.append(dis(li_xy, li_yx, x, y, d, d))

        r = pearsonr(vals, ed)[0]
        s = spearmanr(vals, ed)[0]
        k = kendalltau(vals, ed)[0]

        results[d] = {"P": r, "S": s, "K": k,}

        print(f"LAMAS d={d} | P={r:.4f} S={s:.4f} K={k:.4f}")

    return results


##############################################################################
# MS (ACS BASED)
##############################################################################

def evaluate_ms(pairs, ed):

    vals = []

    for x, y in pairs:

        d1 = len(y) - 1
        d2 = len(x) - 1

        write_file("x_temp", str(x))
        write_file("y_temp", str(y))

        subprocess.run(["../lamas", "-x", "x_temp", "-y", "y_temp", "-d", str(d1), "-o", "out_xy"], check=True)
        subprocess.run(["../lamas", "-x", "y_temp", "-y", "x_temp", "-d", str(d2), "-o", "out_yx"], check=True)

        with open("out_xy") as f:
            li_xy = [int(line.strip()) for line in f if line.strip()]

        with open("out_yx") as f:
            li_yx = [int(line.strip()) for line in f if line.strip()]

        vals.append(dis(li_xy, li_yx, x, y, d1, d2))

    r = pearsonr(vals, ed)[0]
    s = spearmanr(vals, ed)[0]
    k = kendalltau(vals, ed)[0]

    print(f"MS | P={r:.4f} S={s:.4f} K={k:.4f}")

    return r, s, k

##############################################################################
# QGRAM
##############################################################################

def evaluate_qgram(pairs, ed, q_min=4, q_max=20):
    results = {}

    for q in range(q_min, q_max + 1):

        vals = []

        for x, y in pairs:
            write_file("x_temp", str(x))
            write_file("y_temp", str(y))

            subprocess.run(["../qgram", "-x", "x_temp", "-y", "y_temp", "-q", str(q), "-o", "out_q"], check=True)

            with open("out_q") as f:
                vals.append( int(f.readline()))

        r = pearsonr(vals, ed)[0]
        s = spearmanr(vals, ed)[0]
        k = kendalltau(vals, ed)[0]

        results[q] = {"P": r, "S": s, "K": k,}

        print(f"QGRAM q={q} | P={r:.4f} S={s:.4f} K={k:.4f}")

    return results

##############################################################################
# MAIN
##############################################################################
file_path = "../../100000/100000.15"

xs, ys = [], []

with open(file_path) as f:
    for line in f:
        line = line.strip()
        if line.startswith(">"):
            xs.append(line[1:])
        elif line.startswith("<"):
            ys.append(line[1:])

pairs = list(zip(xs[:100], ys[:100]))

ed = compute_ed(pairs)
lamas_results = evaluate_lamas(pairs, ed)
ms_result = evaluate_ms(pairs, ed)
qgram_results = evaluate_qgram(pairs, ed)

##############################################################################
# PLOTTING
##############################################################################

ms_pearson, ms_spearman, ms_kendall = ms_result

all_values = []


for metric in ["P", "S", "K"]:

    all_values.extend(lamas_results[d][metric] for d in lamas_results)
    all_values.extend(qgram_results[q][metric] for q in qgram_results)

all_values.extend([ms_pearson, ms_spearman, ms_kendall])

ymin = min(all_values)
ymax = max(all_values)

padding = 0.05 * (ymax - ymin)

ymin -= padding
ymax += padding

print(f"Y range = [{ymin:.4f}, {ymax:.4f}]")

fig, axes = plt.subplots(3, 1, figsize=(14, 12),)

metrics = [("P", ms_pearson, "Pearson's Correlation"), ("S", ms_spearman, "Spearman Rank Correlation"), ("K", ms_kendall, "Kendall Rank Correlation"),]

for ax, (metric, ms_value, title) in zip(axes, metrics):

    ##########################################################################
    # LAMAS
    ##########################################################################

    d_vals = sorted(lamas_results.keys())

    d_corr = [lamas_results[d][metric] for d in d_vals]

    ax.plot(d_vals, d_corr, color="tab:blue", linewidth=5, label="LAMAS",)

    ##########################################################################
    # MS
    ##########################################################################

    ax.axhline(ms_value, color="tab:green", linestyle="--", linewidth=5, label="MS",)

    ##########################################################################
    # AXIS SETTINGS
    ##########################################################################

    ax.set_xlim(64, 1024)
    ax.set_ylim(ymin, ymax)
    ax.set_ylabel("Correlation")
    ax.grid(True)
    ax.tick_params(axis="x", labelbottom=True, labelsize=18)
    ax.set_xlabel("d", fontsize=20)

    ##########################################################################
    # QGRAM AXIS
    ##########################################################################

    q_vals = sorted(qgram_results.keys())

    qax = ax.twiny()
    qax.set_xlim(4, 21)
    qax.plot(q_vals, [qgram_results[q][metric] for q in q_vals], color="tab:orange", marker="o", linewidth=5, label="Q-GRAM",)
    qax.set_xlabel("q",  fontsize=20)
    qax.tick_params(axis='both', labelsize=18)

    ##########################################################################
    # COMBINED LEGEND
    ##########################################################################

    lines1, labels1 = ax.get_legend_handles_labels()
    lines2, labels2 = qax.get_legend_handles_labels()

    ax.legend(lines1 + lines2, labels1 + labels2, loc="lower right", fontsize=20)
    ax.set_title(title, fontsize=20)

fig.subplots_adjust(hspace=1)
plt.tight_layout()

plt.savefig("100000.25.pdf", bbox_inches="tight")
plt.savefig("100000.25.png", dpi=300,bbox_inches="tight")

plt.show()

