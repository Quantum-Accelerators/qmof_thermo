from pathlib import Path

import pandas as pd
from gen_fig import create_parity_plot, gen_grid, get_energy_columns

DATA_PATH = "data/external/{}_12_25_results.csv"
OUTPUT_DIR = "figures/figure_s18/"
OUTPUT_PATH = OUTPUT_DIR + "{}_{}_{}_{}.png"
Path(OUTPUT_DIR).parent.mkdir(parents=True, exist_ok=True)

graphs = [["QMOF", "uma_dft", "Formation", "A"], ["QMOF", "esen_dft", "Formation", "B"]]

for graph_set in graphs:
    subset, model_type, calc_type, label = graph_set
    print(f"Processing: {subset}, {model_type}, {calc_type}")
    df = pd.read_csv(DATA_PATH.format(model_type))
    df = df[df["type"] == subset]
    x_col, y_col = get_energy_columns(calc_type)
    x, y = df[x_col].values, df[y_col].values
    create_parity_plot(x, y, calc_type, model_type, subset, OUTPUT_PATH, label)
gen_grid(OUTPUT_DIR)
