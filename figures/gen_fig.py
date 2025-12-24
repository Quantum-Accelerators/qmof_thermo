import logging
import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm
from PIL import Image
from scipy.stats import linregress
from sklearn.metrics import mean_absolute_error

title_fontsize = 9
LOGGER = logging.getLogger(__name__)


def get_energy_columns(calc_type):
    if calc_type == "Hull":
        return "dft_ehull", "energy_above_hull"
    elif calc_type == "Formation":
        return "dft_form", "form_energy_per_atom"
    else:
        return "dft_norm_energy", "ml_norm_energy"


def format_labels(calc_type, model_type, subset): 
    calc_label = "Total" if calc_type == "Total Energy" else calc_type
    model_label = "UMA-ODAC"
    subset_label = "MP Reference" if subset == "MP_ref" else subset
    energy_symbol = "E" if calc_type == "Total Energy" else r"\Delta E"
    return calc_label, model_label, subset_label, energy_symbol


def create_parity_plot(x, y, calc_type, model_type, subset, OUTPUT_PATH, label=""):
    mask = ~(np.isnan(x) | np.isnan(y))
    x, y = x[mask], y[mask]

    slope, intercept, r_value, p_value, std_err = linregress(x, y)
    r2 = r_value**2
    slope * x + intercept
    mae = mean_absolute_error(x, y)

    # set some bounds
    x_min, x_max = x.min() - 0.05, x.max() + 0.05
    y_min, y_max = y.min() - 0.05, y.max() + 0.05

    # graph specific bounds
    if calc_type == "Hull":
        if subset == "QMOF":
            x_min, x_max = -0.05, 0.805
            y_min, y_max = -0.05, 0.805
        else:
            y_min, y_max = -0.05, 2.505
            x_min, x_max = -0.05, 2.505
    elif calc_type == "Formation":
        if subset == "QMOF":
            y_min, y_max = -1.75, 0.85
            x_min, x_max = -1.75, 0.85
        else:
            y_min, y_max = -3.25, 2.65
    else:  # total energy
        if subset == "MP_ref":
            y_min, y_max = -15.50, 0
            x_min, x_max = y_min, y_max
        else:
            x_min, x_max = y_min, y_max

    fig, ax = plt.subplots(figsize=(3.25, 2.5))
    hexbin = ax.hexbin(
        x,
        y,
        gridsize=35,
        mincnt=1,
        norm=LogNorm(vmin=1, vmax=1e3),
        extent=[x_min, x_max, y_min, y_max],
        linewidths=0.2,
    )
    cbar = plt.colorbar(hexbin, ax=ax, pad=0.02)
    cbar.ax.tick_params(labelsize=8)

    location = "upper left"
    if calc_type == "Hull" and subset == "MP_ref":
        ax.plot([], [], " ", label=f"MAE = {mae:.3f} eV/atom")
        location = "upper right"
    else:
        ax.plot([-200, 200], [-200, 200], "b--", lw=2, alpha=0.8, zorder=10)
        x_fit = np.array([x_min, x_max])
        slope * x_fit + intercept
        ax.plot([], [], " ", label=f"MAE = {mae:.3f} eV/atom")

    calc_label, model_label, subset_label, energy_symbol = format_labels(
        calc_type, model_type, subset
    )

    ax.set_xlabel(
        rf"PBE-D3(BJ) ${energy_symbol}_{{\mathrm{{{str.lower(calc_label)}}}}}$ (eV/atom)",
        fontsize=title_fontsize,
    )
    ax.set_ylabel(
        rf"{model_label} ${energy_symbol}_{{\mathrm{{{str.lower(calc_label)}}}}}$ (eV/atom)",
        fontsize=title_fontsize,
    )

    if calc_type == "Formation":
        ax.set_xlabel(
            rf"PBE-D3(BJ) ${energy_symbol}_{{\mathrm{{form}}}}$ (eV/atom)",
            fontsize=title_fontsize,
        )
        ax.set_ylabel(
            rf"{model_label} ${energy_symbol}_{{\mathrm{{form}}}}$ (eV/atom)",
            fontsize=title_fontsize,
        )

    title_text = (
        f"Relaxed {subset_label} (Corrected) Structures"
        if model_type in ("uma_dft", "esen_dft")
        else f"Relaxed {subset_label} Structures"
    )
    ax.set_title(title_text, fontsize=8, pad=6)

    ax.tick_params(which="major", direction="in", length=10, width=1.25, labelsize=8)
    ax.tick_params(which="minor", direction="in", length=5, width=1.25)
    ax.minorticks_on()
    for spine in ax.spines.values():
        spine.set_linewidth(1.25)

    legend = ax.legend(
        loc="lower right",
        bbox_to_anchor=(1.00, 0.06),  # move a bit right (x=1.02) and keep same height
        fontsize=7,
        framealpha=0.9,
        handlelength=0,
        handletextpad=0,
    )
    legend.get_frame().set_edgecolor("black")
    ax.set_xlim(x_min, x_max)
    ax.set_ylim(y_min, y_max)

    # Simplified grid (optional: can disable)
    ax.grid(False)
    fig.patch.set_facecolor("white")

    plt.tight_layout(pad=0.2)
    ax.text(
        -0.295,  # x position
        1.1,  # y position
        label,  # the label
        transform=ax.transAxes,
        fontsize=11,
        fontweight="bold",
        va="top",
        ha="left",
    )
    output_file = OUTPUT_PATH.format(
        label.replace(" ", "_"), subset_label, model_label, calc_label
    )
    plt.savefig(output_file, dpi=600, bbox_inches="tight")
    plt.show()

    LOGGER.info(f"Total points plotted: {len(x):,}")
    LOGGER.info(f"R² = {r2:.3f}")
    LOGGER.info(f"MAE = {mae:.4f} eV/atom")
    LOGGER.info(f"Linear fit: y = {slope:.3f}x + {intercept:.4f}\n")


def gen_grid(OUTPUT_DIR: str | Path):
    imgs = [Image.open(Path(OUTPUT_DIR + p)) for p in sorted(os.listdir(OUTPUT_DIR))]

    w, h = imgs[0].size

    combined = Image.new("RGB", (2 * w, 2 * h), color="white")

    combined.paste(imgs[0], (0, 0))
    combined.paste(imgs[1], (w, 0))
    if len(imgs) > 2:
        combined.paste(imgs[2], (0, h))
        combined.paste(imgs[3], (w, h))
    combined.save(f"{OUTPUT_DIR}grid.jpg")
    LOGGER.info(f"Saved grid image to {OUTPUT_DIR}grid.jpg")
