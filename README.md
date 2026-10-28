# LAMAS: Location-Aware MAtching Statistics

This directory contains the source code, binaries, datasets, and benchmarking scripts for evaluating the LAMAS algorithm alongside baseline methods (`q-gram`, `ED`, `MS`).

---

### Installation
Navigate to your chosen LAMAS algorithm (`Improved`, `Linear` or `Simple`) directory, run the pre-installation script, and compile the source code:

```
cd <algorithm_directory>
./pre-install.sh
make
```

## Usage

Run the compiled executable using the following command structure:
```
./lamas -x <seq_x_filename> -y <seq_y_filename> -d <d_range> -o <output_filename>
```

### Parameters

| Flag | Argument | Description |
| :--- | :--- | :--- |
| **`-x`** | `<seq_x_filename>` | Path to the first input sequence file. |
| **`-y`** | `<seq_y_filename>` | Path to the second input sequence file. |
| **`-d`** | `<d_range>` | Range value for $d$. |
| **`-o`** | `<output_filename>` | Path to the designated output file. |

---

## Benchmarking & Baselines

This project includes built-in scripts to benchmark LAMAS against standard baseline methods:
* **`q-gram`**
* **`ED`** (Edit Distance)
* **`MS`** (Matching Statistics)
