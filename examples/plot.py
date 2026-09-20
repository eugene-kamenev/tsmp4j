#!/usr/bin/env python3
# Reads the sample time series and matrix profile produced by GenerateSample.java
# and renders the README figure (examples/matrix-profile-sample.png).
#
#   python3 examples/plot.py
#
# Requires matplotlib and numpy.
import csv
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
DATA = HERE / "data"
OUT = HERE / "matrix-profile-sample.png"

WINDOW_SIZE = 30
MOTIF_LENGTH = 60
PERIOD = 2 * MOTIF_LENGTH
DISCORD_REPETITION = 4
DISCORD_START = DISCORD_REPETITION * PERIOD
DISCORD_END = DISCORD_START + MOTIF_LENGTH


def read_csv(path):
    with path.open(newline="") as fh:
        rows = list(csv.DictReader(fh))
    return rows


def main():
    series_rows = read_csv(DATA / "timeseries.csv")
    profile_rows = read_csv(DATA / "matrix_profile.csv")

    ts_x = np.array([int(r["index"]) for r in series_rows])
    ts_y = np.array([float(r["value"]) for r in series_rows])

    mp_x = np.array([int(r["index"]) for r in profile_rows])
    mp_raw = np.array([float(r["profile"]) for r in profile_rows])

    # Degenerate windows report +Infinity; keep them out of the plot but keep their x.
    finite = np.isfinite(mp_raw)
    mp_y = np.where(finite, mp_raw, np.nan)

    discord = int(np.argmax(np.where(finite, mp_raw, -np.inf)))
    discord_dist = mp_raw[discord]
    motif_min = float(np.nanmin(mp_y))

    fig, (ax_ts, ax_mp) = plt.subplots(
        2, 1, figsize=(12, 7), dpi=150, sharex=True,
        gridspec_kw={"height_ratios": [1.15, 1.0], "hspace": 0.08},
    )

    # --- Time series ---
    ax_ts.axvspan(DISCORD_START, DISCORD_END, color="#e8b4b4", alpha=0.45, lw=0)
    ax_ts.plot(ts_x, ts_y, color="#1f77b4", lw=1.0)
    ax_ts.set_ylabel("value")
    ax_ts.set_title("Time series — two repeating motifs (sine, sawtooth) and one discord", fontsize=11)
    ax_ts.annotate(
        "discord\n(unique shape)",
        xy=((DISCORD_START + DISCORD_END) / 2, ts_y[DISCORD_START]),
        xytext=(DISCORD_END + 70, 1.1),
        fontsize=9, color="#a33", ha="left", va="center",
        arrowprops=dict(arrowstyle="->", color="#a33", lw=1.0),
    )
    ax_ts.grid(True, ls=":", alpha=0.4)

    # --- Matrix profile ---
    ax_mp.axvspan(DISCORD_START, DISCORD_END, color="#e8b4b4", alpha=0.45, lw=0)
    ax_mp.plot(mp_x, mp_y, color="#2ca02c", lw=1.0)
    ax_mp.axhline(motif_min, color="#888", ls="--", lw=0.8)
    ax_mp.scatter([discord], [discord_dist], color="#d62728", zorder=5)
    ax_mp.annotate(
        f"discord  d={discord_dist:.2f}",
        xy=(discord, discord_dist), xytext=(discord + 90, discord_dist * 0.9),
        fontsize=9, color="#d62728",
        arrowprops=dict(arrowstyle="->", color="#d62728", lw=1.0),
    )
    ax_mp.annotate(
        "recurring motifs match well → low profile",
        xy=(120, motif_min), xytext=(150, motif_min + 0.6),
        fontsize=9, color="#2ca02c",
        arrowprops=dict(arrowstyle="->", color="#2ca02c", lw=1.0),
    )
    ax_mp.set_xlabel("subsequence index")
    ax_mp.set_ylabel("matrix profile distance")
    ax_mp.set_title(f"Matrix profile (window size = {WINDOW_SIZE})", fontsize=11)
    ax_mp.grid(True, ls=":", alpha=0.4)

    fig.tight_layout()
    fig.savefig(OUT)
    print(f"wrote {OUT}")
    print(
        f"series n={ts_x.size}, profile n={mp_x.size}, "
        f"non-finite windows={int((~finite).sum())}, "
        f"discord at {discord} (d={discord_dist:.3f}), motif min={motif_min:.3f}"
    )


if __name__ == "__main__":
    main()
