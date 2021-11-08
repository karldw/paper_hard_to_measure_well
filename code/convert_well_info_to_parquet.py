#!/usr/bin/env python3

import sys
import itertools
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
import pandas as pd
import pyarrow
import pyarrow.parquet as pq
import logging
from functools import partial

if sys.version_info < (3, 7):
    raise AssertionError("Need python version 3.7 or higher")

# Create a logger to output any messages we might have...
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
handler = logging.StreamHandler()
formatter = logging.Formatter("%(message)s")
handler.setFormatter(formatter)
logger.addHandler(handler)


def dtypes_to_schema(dtype_dict):
    schema = []
    type_conversion = {
        "Int64": pyarrow.int64(),
        "Int32": pyarrow.int32(),
        "int64": pyarrow.int64(),
        "int32": pyarrow.int32(),
        "int": pyarrow.int32(),
        "float64": pyarrow.float64(),
        "float32": pyarrow.float32(),
        "float": pyarrow.float32(),
        "str": pyarrow.string(),
        "O": pyarrow.string(),
        # Note that the arguments to pyarrow.dictionary changed between v0.13.0
        # and v0.14.1
        "category": pyarrow.dictionary(
            pyarrow.int32(), pyarrow.string(), ordered=False
        ),
        "date32": pyarrow.date32(),  # Note: date32 isn't a python type
    }
    for field_name, dtype in dtype_dict.items():
        type = type_conversion[dtype]
        schema.append(pyarrow.field(field_name, type, nullable=True))
    return pyarrow.schema(schema)


def extension_int_to_float(df, exclude=[]):
    """
    You can delete this function once this issue is resolved:
    https://issues.apache.org/jira/browse/ARROW-5379
    """
    extension_int_types = {
        pd.Int8Dtype(),
        pd.Int16Dtype(),
        pd.Int32Dtype(),
        pd.Int64Dtype(),
    }
    new_types = {
        col: "float64"
        for col in df.columns
        if df[col].dtypes in extension_int_types and col not in exclude
    }
    return df.astype(new_types)


def fix_monthly_dates(df):
    df = coerce_to_date(df, ["date"], drop_bad=True)
    # Note: these are sometimes NA
    df["month"] = df["date"].dt.month.astype("Int32")
    df["year"] = df["date"].dt.year.astype("Int32")
    del df["date"]
    return df


def read_csv(filepath_or_buffer, *args, **kwargs):
    """Thin wrapper around pd.read_csv"""
    logger.info(f"  Reading {filepath_or_buffer}")
    return pd.read_csv(filepath_or_buffer, *args, **kwargs)


def convert_monthly(csv_file, out_dir):
    csv_file = Path(csv_file)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    # Note: we're using python 3.7+, so dicts keep their order
    pd_dtypes = {
        "entity_id": "Int32",
        "api": "str",
        "api_list": "str",
        "date": "str",
        "oil": "str",  # will convert to float later
        "gas": "str",  # will convert to float later
        "water": "str",  # will convert to float later
        "well_count": "str",  # will convert to Int later
        "days": "str",  # will convert to Int later
        "daily_avg_oil": "str",  # will convert to float later
        "daily_avg_gas": "str",  # will convert to float later
        "daily_avg_water": "str",  # will convert to float later
        "reservoir": "str",
        "well_name": "str",
        "well_number": "str",
        "operator_alias": "str",
        "production_type": "str",
        "production_status": "str",
        "entity_type": "str",
    }
    # Need to read these as str, then convert to float, because of errors.
    float_cols = [
        "oil",
        "gas",
        "water",
        "daily_avg_oil",
        "daily_avg_gas",
        "daily_avg_water",
    ]
    int_cols = ["well_count", "days"]
    output_dtypes = pd_dtypes.copy()
    del output_dtypes["date"]
    output_dtypes["year"] = "int32"
    output_dtypes["month"] = "int32"
    for c in float_cols:
        output_dtypes[c] = "float"
    for c in int_cols:
        output_dtypes[c] = "int"
    pq_schema = dtypes_to_schema(output_dtypes)

    csv_iter = read_csv(
        csv_file,
        dtype=pd_dtypes,
        chunksize=3_000_000,
        index_col=False,
        names=list(pd_dtypes.keys()),
        on_bad_lines="warn",
        header=0,
        low_memory=True,
        na_values={"(N/A)"},
        keep_default_na=True,
    )
    partition_cols = ["year"]
    for chunk in csv_iter:
        chunk = (
            chunk.pipe(fix_monthly_dates)
            .pipe(coerce_to_float, cols=float_cols)
            .pipe(coerce_to_integer, cols=int_cols)
            .pipe(extension_int_to_float)
            .pipe(fill_na_for_partitioning, cols=partition_cols)
        )
        table = pyarrow.Table.from_pandas(chunk, preserve_index=False, schema=pq_schema)
        pq.write_to_dataset(
            table, root_path=str(out_dir), partition_cols=partition_cols, version="2.0"
        )
    logger.info(f"  Done with monthly production {csv_file}")


def fill_na_for_partitioning(df, cols):
    """
    Replace NaN or null values with non-NA placeholders.

    These placeholders vary by column type: "NA" for string or category,
    -9 for numeric.

    The goal here is to prevent write_to_dataset from dropping groups with
    missing-valued partition columns. This function can be cut out when this
    pyarrow issue is fixed:
    https://issues.apache.org/jira/projects/ARROW/issues/ARROW-7345
    """
    orig_order = df.columns
    assert orig_order.nunique() == df.shape[1]

    to_fill = df[cols]
    fill_str = to_fill.select_dtypes(include=["object", "category"]).fillna("NA")
    fill_num = to_fill.select_dtypes(include="number").fillna(-9)
    unfilled = to_fill.columns.difference(fill_str.columns).difference(fill_num.columns)
    if len(unfilled) > 0:
        raise NotImplementedError(f"Can't fill columns: {unfilled.tolist()}")
    # Everything else:
    remainder = df.drop(cols, axis=1)
    out = pd.concat([fill_str, fill_num, remainder], axis=1)[orig_order.tolist()]
    # Reset columns to original ordering:
    return out


def coerce_cols(df, cols, target_type, drop_bad=True, **kwargs):
    """Coerce `cols` to `target_type`, optionally dropping rows that fail.

    A very small number of rows fail to parse because they have things like
    "WELL" in the oil quantity column. Just drop these.
    But because of this, we have to read all of these column as str, then convert

    args: df a dataframe
    cols: columns, currently string, that will be coerced
    target_type: what type to coerce to?
    drop_bad: should we drop the whole row if the coersion fails? (Default True)
    kwargs: arguments passed on to pd.to_datetime or pd.to_numeric.
    """

    rename_dict = {c: c + "_str" for c in cols}
    if df.columns.isin(rename_dict.values()).any():
        raise ValueError("temp rename column already exists")
    df = df.rename(columns=rename_dict)

    if target_type == "datetime":
        conversion_fn = pd.to_datetime
    else:
        conversion_fn = partial(pd.to_numeric, downcast=target_type)

    for new_col, old_col in rename_dict.items():
        try:
            # Happy case first
            df[new_col] = conversion_fn(df[old_col], errors="raise", **kwargs)
        except ValueError:
            df[new_col] = conversion_fn(df[old_col], errors="coerce", **kwargs)
            newly_missing = df[new_col].isna() & df[old_col].notna()
            newly_missing_count = newly_missing.sum()
            logger.info(
                f"  {newly_missing_count} values failed to parse as {target_type}."
                + f" Here are a few: "
                + ", ".join(df.loc[newly_missing, old_col].iloc[:6])
            )
            if drop_bad:
                # Drop rows that failed to parse (very few of these)
                df = df.loc[~newly_missing]
        finally:
            del df[old_col]
    return df


coerce_to_date = partial(coerce_cols, target_type="datetime", format="%Y-%m-%d")
coerce_to_float = partial(coerce_cols, target_type="float")
coerce_to_integer = partial(coerce_cols, target_type="integer")


def convert_headers(csv_file, out_dir):
    csv_file = Path(csv_file)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    pd_dtypes = {
        "api": "str",  # API/UWI
        "operator_alias_legacy": "str",  # Operator Alias (Legacy)
        "operator_company_name": "str",  # Operator Company Name
        "operator_reported": "str",  # Operator (Reported)
        "operator_ticker": "str",  # Operator Ticker
        "well_name": "str",  # Well/Lease Name
        "well_number": "str",  # Well Number
        "entity_type": "str",  # Entity Type
        "county": "str",  # County/Parish
        "di_basin": "str",  # DI Basin
        "di_play": "str",  # DI Play
        "di_subplay": "str",  # DI Subplay
        "reservoir": "str",  # Reservoir
        "production_type": "str",  # Production Type
        "producing_status": "str",  # Producing Status
        "drill_type": "str",  # Drill Type
        "first_prod_date": "str",  # First Prod Date
        "last_prod_date": "str",  # Last Prod Date
        "cum_gas": "float",  # Cum Gas
        "cum_oil": "float",  # Cum Oil
        "cum_boe": "float",  # Cum BOE
        "cum_water": "float",  # Cum Water
        "cum_mmcfge": "float",  # Cum MMCFGE
        "cum_bcfge": "float",  # Cum BCFGE
        "daily_gas": "float",  # Daily Gas
        "daily_oil": "float",  # Daily Oil
        "first_month_oil": "float",  # First Month Oil
        "first_month_gas": "float",  # First Month Gas
        "first_month_water": "float",  # First Month Water
        "first_6_oil": "float",  # First 6 Oil
        "first_6_gas": "float",  # First 6 Gas
        "first_6_boe": "float",  # First 6 BOE
        "first_6_water": "float",  # First 6 Water
        "first_12_oil": "float",  # First 12 Oil
        "first_12_gas": "float",  # First 12 Gas
        "first_12_boe": "float",  # First 12 BOE
        "first_12_mmcfge": "float",  # First 12 MMCFGE
        "first_12_water": "float",  # First 12 Water
        "first_24_oil": "float",  # First 24 Oil
        "first_24_gas": "float",  # First 24 Gas
        "first_24_boe": "float",  # First 24 BOE
        "first_24_mmcfge": "float",  # First 24 MMCFGE
        "first_24_water": "float",  # First 24 Water
        "first_60_oil": "float",  # First 60 Oil
        "first_60_gas": "float",  # First 60 Gas
        "first_60_boe": "float",  # First 60 BOE
        "first_60_water": "float",  # First 60 Water
        "first_60_mmcfge": "float",  # First 60 MMCFGE
        "prac_ip_oil_daily": "float",  # Prac IP Oil Daily
        "prac_ip_gas_daily": "float",  # Prac IP Gas Daily
        "prac_ip_cfged": "float",  # Prac IP CFGED
        "prac_ip_boe": "float",  # Prac IP BOE
        "latest_oil": "float",  # Latest Oil
        "latest_gas": "float",  # Latest Gas
        "latest_water": "float",  # Latest Water
        "prior_12_oil": "float",  # Prior 12 Oil
        "prior_12_gas": "float",  # Prior 12 Gas
        "prior_12_water": "float",  # Prior 12 Water
        "last_test_date": "str",  # Last Test Date
        "last_flow_pressure": "float",  # Last Flow Pressure
        "last_whsip": "float",  # Last WHSIP
        "2nd_month_gor": "float",  # 2nd Month GOR
        "latest_gor": "float",  # Latest GOR
        "cum_gor": "float",  # Cum GOR
        "last_12_yield": "float",  # Last 12 Yield
        "2nd_month_yield": "float",  # 2nd Month Yield
        "latest_yield": "float",  # Latest Yield
        "peak_gas": "float",  # Peak Gas
        "peak_gas_month_no.": "Int32",  # Peak Gas Month No.
        "peak_oil": "float",  # Peak Oil
        "peak_oil_month_no.": "Int32",  # Peak Oil Month No.
        "peak_boe": "float",  # Peak BOE
        "peak_boe_month_no.": "Int32",  # Peak BOE Month No.
        "peak_mmcfge": "float",  # Peak MMCFGE
        "peak_mmcfge_month_no.": "Int32",  # Peak MMCFGE Month No.
        "upper_perforation": "float",  # Upper Perforation
        "lower_perforation": "float",  # Lower Perforation
        "gas_gravity": "float",  # Gas Gravity
        "oil_gravity": "float",  # Oil Gravity
        "completion_date": "str",  # Completion Date
        "well_count": "Int32",  # Well Count
        "max_active_wells": "Int32",  # Max Active Wells
        "months_produced": "Int32",  # Months Produced
        "gas_gatherer": "str",  # Gas Gatherer
        "oil_gatherer": "str",  # Oil Gatherer
        "lease_number": "str",  # Lease Number
        "spud_date": "str",  # Spud Date
        "measured_depth_td": "float",  # Measured Depth (TD)
        "true_vertical_depth": "float",  # True Vertical Depth
        "gross_perforated_interval": "float",  # Gross Perforated Interval
        "field": "str",  # Field
        "state": "str",  # State
        "district": "str",  # District
        "aapg_geologic_province": "str",  # AAPG Geologic Province
        "country": "str",  # Country
        "section": "str",  # Section
        "township": "str",  # Township
        "range": "str",  # Range
        "abstract": "str",  # Abstract
        "block": "str",  # Block
        "survey": "str",  # Survey
        "ocs_area": "str",  # OCS Area
        "pgc_area": "str",  # PGC Area
        "surface_latitude_wgs84": "float",  # Surface Latitude (WGS84)
        "surface_longitude_wgs84": "float",  # Surface Longitude (WGS84)
        "last_12_oil": "float",  # Last 12 Oil
        "last_12_gas": "float",  # Last 12 Gas
        "last_12_water": "float",  # Last 12 Water
        "entity_id": "int32",  # Entity ID
    }
    date_cols = [
        "first_prod_date",
        "last_prod_date",
        "last_test_date",
        "completion_date",
        "spud_date",
    ]

    output_dtypes = pd_dtypes.copy()
    for d in date_cols:
        output_dtypes[d] = "date32"
    output_dtypes["first_prod_year"] = "int32"
    pq_schema = dtypes_to_schema(output_dtypes)
    csv_iter = read_csv(
        csv_file,
        dtype=pd_dtypes,
        chunksize=1_000_000,
        index_col=False,
        names=list(pd_dtypes.keys()),
        on_bad_lines="warn",
        parse_dates=date_cols,
        header=0,
        low_memory=True,
    )
    partition_cols = ["first_prod_year"]
    for chunk in csv_iter:
        chunk["first_prod_year"] = chunk["first_prod_date"].dt.year
        chunk = chunk.pipe(extension_int_to_float, exclude=partition_cols).pipe(
            fill_na_for_partitioning, cols=partition_cols
        )
        table = pyarrow.Table.from_pandas(chunk, preserve_index=False, schema=pq_schema)
        pq.write_to_dataset(
            table, root_path=str(out_dir), partition_cols=partition_cols, version="2.0"
        )
    logger.info(f"  Done with headers: {csv_file}")


def make_uniform_bins(x, num):
    import numpy as np

    # Bins are evenly spaced in the range of unsigned 64bit ints
    bins = np.linspace(0, np.iinfo("uint64").max - 1024, num=num, dtype="uint64")
    # Hash (determanistic) of the values is a uint64
    hashes = pd.util.hash_pandas_object(x, index=False)
    assert hashes.dtypes == "uint64"
    return np.digitize(hashes, bins)


def regroup_parquet_fragments(directory, recursive=True):
    """
    Read all the parquet files in a directory and write a single file.
    If there are inner directories, `regroup_parquet_fragments` will be called
    recursively on those.

    Why:
    Because I'm calling write_to_dataset repeatedly for multiple chunks of data,
    I sometimes get multiple fragments for the same partition. That's a pain,
    so this function fixes that. After running it, there will be only one file
    (or zero if the directory was initially empty), named `file.parquet`.
    """
    assert directory.is_dir()
    pq_files = []
    dirs = []
    other = []
    for p in directory.iterdir():
        if p.is_file() and p.suffix == ".parquet":
            pq_files.append(p)
        elif p.is_dir():
            dirs.append(p)
        else:
            other.append(p)
    if len(other):
        raise ValueError("Unexpected files found: " + ", ".join(other))
    if recursive:
        [regroup_parquet_fragments(d, recursive=recursive) for d in dirs]

    outfile = directory.joinpath("file.parquet")
    assert not outfile.exists()
    if len(pq_files) == 0:
        return
    elif len(pq_files) == 1:
        pq_files[0].rename(outfile)
        return
    else:
        data = pq.ParquetDataset(pq_files, memory_map=True).read(
            use_pandas_metadata=True
        )
        pq.write_table(data, outfile, version="2.0")
        [f.unlink() for f in pq_files]


def main(snakemake):
    # Fill with partial so we can easily map with one argument.
    monthly = partial(convert_monthly, out_dir=Path(snakemake.output["monthly"]))
    headers = partial(convert_headers, out_dir=Path(snakemake.output["headers"]))

    with ProcessPoolExecutor(max_workers=snakemake.threads) as ex:
        ex.map(monthly, snakemake.input["monthly"])
        ex.map(headers, snakemake.input["headers"])
    regroup_parquet_fragments(Path(snakemake.output["headers"]))
    regroup_parquet_fragments(Path(snakemake.output["monthly"]))


if __name__ == "__main__":
    main(snakemake)
