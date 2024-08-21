import pandas as pd
import pyproj

import matplotlib.pyplot as plt
from sklearn.cluster import DBSCAN


def read_data(in_file):
    df_input = pd.read_parquet(in_file)
    orig_varnames = set(df_input.columns)
    # These should use the same values as CRS_LONGLAT_int and CRS_PROJECT_int in shared_functions.r
    # see https://epsg.io/4326 and https://epsg.io/6350
    proj = pyproj.Transformer.from_crs(
        "EPSG:4326", "EPSG:6350", always_xy=True, allow_ballpark=False
    )
    df_input["x_coord"], df_input["y_coord"] = proj.transform(
        df_input["well_pad_lon"], df_input["well_pad_lat"], errcheck=True
    )
    return df_input, orig_varnames


def get_basin_limits(df):
    """the following steps are used to ensure the same grid size in subplots
    so the size of each basin can directly be compared
    """
    basin_lim = (
        df.groupby("basin")
        .agg(
            {
                "x_coord": ["max", "min"],
                "y_coord": ["max", "min"],
            }
        )
        .reset_index()
    )
    basin_lim.columns = basin_lim.columns.droplevel(0)
    basin_lim.columns = ["basin", "x_max", "x_min", "y_max", "y_min"]
    basin_lim["x_d"] = basin_lim["x_max"] - basin_lim["x_min"]
    basin_lim["y_d"] = basin_lim["y_max"] - basin_lim["y_min"]
    x_range = 1.1 * basin_lim["x_d"].max()
    y_range = 1.1 * basin_lim["y_d"].max()
    x_range, y_range = max(x_range, y_range), max(x_range, y_range)
    lim_parameters = {}
    for b in basin_lim.basin:
        temp_df = basin_lim.loc[basin_lim["basin"] == b]
        x_mid = (temp_df["x_max"].iloc[0] + temp_df["x_min"].iloc[0]) / 2
        y_mid = (temp_df["y_max"].iloc[0] + temp_df["y_min"].iloc[0]) / 2
        lim_parameters[b] = [
            x_mid - x_range / 2,
            x_mid + x_range / 2,
            y_mid - y_range / 2,
            y_mid + y_range / 2,
        ]
    return lim_parameters


def plot_before_clustering(df_input, filename):
    """plot all basins before clustering"""
    lim_parameters = get_basin_limits(df_input)
    fig, axes = plt.subplots(2, 2, figsize=(40, 40), dpi=300)
    for i, b in enumerate(df_input["basin"].unique()):
        df = df_input.loc[df_input["basin"] == b]
        ax = axes[i // 2, i % 2]
        for basin_label in df["basin"].unique():
            basin_data = df[df["basin"] == basin_label]
            ax.scatter(
                basin_data["x_coord"],
                basin_data["y_coord"],
                alpha=0.05,
                label=f"{basin_label}",
            )
        ax.set_xlim(lim_parameters[b][0:2])
        ax.set_ylim(lim_parameters[b][2:4])
        ax.set_title(f"{b}")
        ax.set_xlabel("X Coordinate")
        ax.set_ylabel("Y Coordinate")
        ax.legend()
    plt.savefig(filename, dpi=300)


def plot_after_clustering(df_output, filename):
    """plot all basins after clustering"""
    lim_parameters = get_basin_limits(df_output)
    fig, axes = plt.subplots(2, 2, figsize=(40, 40), dpi=300)
    for i, b in enumerate(df_output["basin"].unique()):
        df = df_output.loc[df_output["basin"] == b]
        ax = axes[i // 2, i % 2]
        for cluster_label in df["cluster"].unique():
            cluster_data = df[df["cluster"] == cluster_label]
            ax.scatter(
                cluster_data["x_coord"],
                cluster_data["y_coord"],
                alpha=0.05,
                label=f"Cluster {cluster_label}",
            )
        ax.set_xlim(lim_parameters[b][0:2])
        ax.set_ylim(lim_parameters[b][2:4])
        ax.set_title(f"{b}")
        ax.set_xlabel("X Coordinate")
        ax.set_ylabel("Y Coordinate")
        ax.legend()
    plt.savefig(filename, dpi=300)


def run_clustering(df_input):
    """run dbscan for clustering"""
    dbscan_param = {
        "San Joaquin": [10500, 150],
        "Other California": [10500, 15],
        "San Juan": [10500, 150],
    }

    df_output = pd.DataFrame()
    for b in df_input.basin.unique():
        temp_df = df_input.loc[df_input["basin"] == b].copy()
        dbscan = DBSCAN(eps=dbscan_param[b][0], min_samples=dbscan_param[b][1])
        temp_df["cluster"] = (
            dbscan.fit_predict(temp_df[["x_coord", "y_coord"]])
        ).astype(str)
        temp_df["cluster"] = b + "_" + temp_df["cluster"]
        df_output = pd.concat([df_output, temp_df])
    # df_output.loc[df_output["cluster"] == "-1", "cluster"] = "outlier"
    df_output = df_output.sort_values(by="cluster")

    return df_output


if __name__ == "__main__":
    in_file = snakemake.input["cleaned_matched_obs"]
    df_input, orig_varnames = read_data(in_file)
    df_output = run_clustering(df_input)
    out_file = snakemake.output[0]
    df_output.drop(columns=["x_coord", "y_coord"]).to_parquet(out_file)
    assert orig_varnames.issubset(set(df_output.columns))

    # Optional:
    # plot_before_clustering(df_input, snakemake.output["plot_well_maps"])
    # plot_after_clustering(df_output, snakemake.output["map_dbscan_clusters"])
