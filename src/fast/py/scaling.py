import pandas as pd
import matplotlib.pyplot as plt

parent_local = "../../../out/fast"
parent_bw_1 = "../../../out/fast_bw_1/fast"
parent_bw_2 = "../../../out/fast_bw_2/fast"

parent = parent_bw_2

df = pd.read_csv(
        f"{parent}/stats.tsv",
        sep = "\t",
        names = ["n", "i", "mlups", "total", "str", "sld", "cpy", "com", "col"]
)

# średnia po wszystkich obrazach
avg = (
    df.groupby("n", as_index=False)
      .agg({
          "mlups": "mean",
          "total": "mean",
          "str": "mean",
          "sld": "mean",
          "cpy": "mean",
          "com": "mean",
          "col": "mean",
      })
)

# MLUPS
plt.figure(figsize=(7, 4))
plt.plot(avg["n"], avg["mlups"], "o-")
plt.xlabel("Number of processes")
plt.ylabel("MLUPS")
plt.grid(True)
plt.tight_layout()
plt.savefig(f"{parent}/mlups_scaling.png", dpi=300)


# Relative execution time
relative = avg.copy()

for col in ["str", "sld", "cpy", "com", "col"]:
    relative[col] = 100.0 * relative[col] / relative["total"]

plt.figure(figsize=(8,5))

for col, label in [
    ("str", "Streaming"),
    # ("sld", "Slide"),
    # ("cpy", "Copy"),
    ("com", "Communication"),
    ("col", "Collision"),
]:
    plt.plot(relative["n"], relative[col], "o-", label=label)

plt.xlabel("Number of processes")
plt.ylabel("Share of total runtime [%]")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig(f"{parent}/runtime_share.png", dpi=300)


# Absolute times
plt.figure(figsize=(8,5))

for col, label in [
    ("str", "Streaming"),
    # ("sld", "Slide"),
    # ("cpy", "Copy"),
    ("com", "Communication"),
    ("col", "Collision"),
]:
    plt.plot(avg["n"], avg[col], "o-", label=label)

plt.xlabel("Number of processes")
plt.ylabel("Time [s]")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig(f"{parent}/runtime_absolute.png", dpi=300)