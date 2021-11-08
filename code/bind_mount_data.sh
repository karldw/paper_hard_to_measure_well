#!/bin/bash


set -eu -o pipefail


code_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
base_dir="$code_dir/.."

mkdir -p "$base_dir/data" "$base_dir/scratch"

data_drive="/home/karl/data/Dropbox/large_files/methane_abatement/data"
scratch_drive="/home/karl/data/scratch/methane_abatement"

if [ ! -d "$data_drive" ]; then
    echo "expected data source files don't exist"
    exit 1
fi
if [ ! -d "$scratch_drive" ]; then
    echo "expected scratch directory don't exist"
    exit 1
fi

sudo mount --bind "$data_drive"    "$base_dir/data"
sudo mount --bind "$scratch_drive" "$base_dir/scratch"
