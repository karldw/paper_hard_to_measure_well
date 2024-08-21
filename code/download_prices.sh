#!/bin/bash

# This script downloads the files below, but once they're downloaded, I view
# them as immutable data, not something to re-run with Snakemake.
# (The files below get updated over time)
# Here's what the snakemake rule would look like:
#
# rule download_prices:
#     output:
#         "data/fred/CPILFENS.csv",
#         "data/eia/NG_PRI_FUT_S1_D.xls",
#         "data/eia/NG_PROD_WELLS_S1_A.xls",
#         "data/eia/NG_PROD_OILWELLS_S1_A.xls",
#         "data/eia/NG_PROD_SUM_A_EPG0_VGM_MMCF_M.xls",
#     script:
#         "code/download_prices.sh"


set -eu -o pipefail

code_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null && pwd )"
# eia_dir="$code_dir/../data/eia"
fred_dir="$code_dir/../data/fred"

# cd "$eia_dir"
# curl --remote-name-all \
#     --url "https://www.eia.gov/dnav/ng/xls/NG_PRI_FUT_S1_D.xls" \
#     --url "https://www.eia.gov/dnav/ng/xls/NG_PROD_WELLS_S1_A.xls" \
#     --url "https://www.eia.gov/dnav/ng/xls/NG_PROD_OILWELLS_S1_A.xls" \
#     --url "https://www.eia.gov/dnav/ng/xls/NG_PROD_SUM_A_EPG0_VGM_MMCF_M.xls"

cd "$fred_dir"
curl --remote-name-all \
    --url "https://fred.stlouisfed.org/series/CPILFENS/downloaddata/CPILFENS.csv"
