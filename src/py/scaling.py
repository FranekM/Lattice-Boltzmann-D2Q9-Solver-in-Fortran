import os
import csv
import time
import subprocess
import matplotlib.pyplot as plt

# ---------------------------------
# konfiguracja
# ---------------------------------

EXE = "../../lbm.exe"
EXPERIMENT = "sliding-lid"

RESULTS_FILE = "benchmark.tsv"
PLOT_FILE = "benchmark_mlups.png"

Nx = 300
Ny = 300
Nt = 10000

results = []

# ---------------------------------
# benchmark
# ---------------------------------

for side in range(1, 11):
    nproc = side * side

    env = os.environ.copy()
    env["FOR_COARRAY_NUM_IMAGES"] = str(nproc)

    print(f"Running with {nproc} images...")

    start = time.perf_counter()

    subprocess.run(
        [EXE, EXPERIMENT],
        env = env,
        check = True,
    )

    elapsed = time.perf_counter() - start

    mlups = Nx * Ny * Nt / elapsed / 1e6

    results.append((nproc, elapsed, mlups))

    print(
        f"{nproc:3d} images   "
        f"{elapsed:8.2f} s    "
        f"{mlups:8.2f} MLUPS"
    )

# ---------------------------------
# zapis wyników
# ---------------------------------

with open(RESULTS_FILE, "w", newline="") as f:

    writer = csv.writer(f, delimiter="\t")

    writer.writerow([
        "images",
        "time_s",
        "mlups"
    ])

    writer.writerows(results)

print(f"Results saved to {RESULTS_FILE}")

# ---------------------------------
# wykres
# ---------------------------------

procs = [r[0] for r in results]
mlups = [r[2] for r in results]

plt.figure(figsize=(7,5))
plt.plot(procs, mlups, "o-", lw=2)

plt.xlabel("Number of images")
plt.ylabel("MLUPS")
plt.title("LBM parallel performance")
plt.grid(True)

plt.tight_layout()
plt.savefig(PLOT_FILE, dpi=300)

print(f"Plot saved to {PLOT_FILE}")