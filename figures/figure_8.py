import logging
from pathlib import Path

import pandas as pd
from gen_fig import create_parity_plot, gen_grid, get_energy_columns

LOGGER = logging.getLogger(__name__)

DATA_PATH = "data/external/{}_12_25_results.csv"
OUTPUT_DIR = "figures/figure_8/"
OUTPUT_PATH = OUTPUT_DIR + "{}_{}_{}_{}.png"
Path(OUTPUT_DIR).mkdir(parents=True, exist_ok=True)


def main():
    logging.basicConfig(level=logging.INFO)

    graphs = [
        ["QMOF", "uma", "Hull", "A"],
        ["QMOF", "uma", "Formation", "B"],
        ["MP_ref", "uma", "Formation", "C"],
        ["QMOF", "uma_dft", "Hull", "D"],
    ]

    for graph_set in graphs:
        subset, model_type, calc_type, label = graph_set
        LOGGER.info(f"Processing: {subset}, {model_type}, {calc_type}")
        df = pd.read_csv(DATA_PATH.format(model_type))
        df = df[df["type"] == subset]
        x_col, y_col = get_energy_columns(calc_type)
        x, y = df[x_col].values, df[y_col].values
        create_parity_plot(x, y, calc_type, model_type, subset, OUTPUT_PATH, label)
    gen_grid(OUTPUT_DIR)


if __name__ == "__main__":
    main()
