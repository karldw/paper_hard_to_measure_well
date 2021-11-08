#!/usr/bin/env python3


import argparse
import linecache
import itertools
import io
import pandas as pd
from pyarrow.parquet import ParquetWriter
import pyarrow


def parse_command_line():
    """
    Parse command line arguments. See the -h option.
    :param argv: arguments on the command line must include caller file name.
    """
    parser = argparse.ArgumentParser(
        description="""
        Convert a CSV file to parquet. May be larger than memory.

        Note that because pandas -> parquet cannot yet handle extentiontypes,
        numpy ints can't handle NA values, and we're not reading all the data
        in one go, all integers will be converted to floats until this issue is
        resolved: https://issues.apache.org/jira/browse/ARROW-5379
        """
    )
    parser.add_argument("input_csv")
    parser.add_argument("output_parquet")
    parser.add_argument(
        "--dtypes",
        help="""Specify CSV's data types. Overrides Pandas' type inference.
        Example:
        --dtypes "{'a': 'float64', 'b': 'Int64'}"
        """,
    )
    arguments = parser.parse_args()
    return arguments


def get_csv_types(csv, **kwargs):
    """Infer column types from `csv`, only using the first 10000 rows.

    Pass argument dtype in kwargs overrides pandas dtype inference.
    Reads any type of file pd.read_csv can work with.
    """
    sample_data = pd.read_csv(csv, low_memory=False, nrows=10000, **kwargs)

    # Note: get schema before making nullable. This whole dance is mainly useful
    # to deal with https://issues.apache.org/jira/browse/ARROW-5379
    schema = dtypes_to_schema(sample_data.dtypes)
    dtypes = dict(make_dtypes_nullable(sample_data.dtypes))
    # If the csv was an open file, rather than a filename, we just changed the
    # position. Reset to the beginning of the file so pd.read_csv works.
    try:
        csv.seek(0)
    except AttributeError:
        pass
    return schema, dtypes


def make_dtypes_nullable(dtypes):
    dtypes[dtypes == "int64"] = "float64"
    dtypes[dtypes == "int32"] = "float32"
    return dtypes


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
        "object": pyarrow.string(),
        # Note that the arguments to pyarrow.dictionary changed between v0.13.0
        # and v0.14.1
        "category": pyarrow.dictionary(
            pyarrow.int32(), pyarrow.string(), ordered=False
        ),
        "date32": pyarrow.date32(),  # Note: date32 isn't a python type
    }
    for field_name, dtype in dtype_dict.items():
        type = type_conversion[str(dtype)]
        schema.append(pyarrow.field(field_name, type, nullable=True))
    return pyarrow.schema(schema)


def convert_csv_to_parquet(
    csv, parquet, csv_dtypes=None, chunksize=1_000_000, **kwargs
):
    """
    Read a CSV in chunks and write to a parquet file.

    kwargs are passed to pandas' read_csv. To specify input types, use
    `csv_dtypes` rather than `dtype`.

    Note that this creates a single parquet file. It may be a better idea to
    write out a partitioned parquet dataset.
    """
    schema, dtype = get_csv_types(csv, dtype=csv_dtypes, **kwargs)
    reader = pd.read_csv(csv, dtype=dtype, chunksize=chunksize, **kwargs)
    with ParquetWriter(parquet, schema) as pqwriter:
        for chunk in reader:
            table = pyarrow.Table.from_pandas(chunk, schema=schema)
            pqwriter.write_table(table)


if __name__ == "__main__":
    args = parse_command_line()
    convert_csv_to_parquet(args.input_csv, args.output_parquet, csv_dtypes=args.dtypes)
