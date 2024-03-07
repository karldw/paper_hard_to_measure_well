from zipfile import ZipFile
from pathlib import Path


drillinginfo_zip = Path(snakemake.output["drillinginfo_zip"])
public_zip = Path(snakemake.output["public_zip"])


def write_drillinginfo_data(zip_filename):
    di_dir = Path("data/drillinginfo")
    assert di_dir.is_dir()

    to_write = [
        di_dir / "entity_production_headers" / f
        for f in (
            "CA_Kern_Production_Headers.CSV.gz",
            "CA_notKern_Production_Headers.CSV.gz",
            "CO_Production_Headers.CSV.gz",
            "NM_Production_Headers.CSV.gz",
        )
    ]
    to_write += [
        di_dir / "entity_production_monthly" / f
        for f in (
            "CA_Kern_Producing_Entity_Monthly_Production.CSV.gz",
            "CA_notKern_Producing_Entity_Monthly_Production.CSV.gz",
            "CO_Producing_Entity_Monthly_Production.CSV.gz",
            "NM_Producing_Entity_Monthly_Production.CSV.gz",
        )
    ]

    with ZipFile(zip_filename, "w") as z:
        [z.write(f) for f in to_write]


def write_public(zip_filename):
    # Public / sharable data
    data = Path("data")
    to_write = [
        data / "drillinginfo" / "docs" / f
        for f in (
            "data_dict_production.pdf",
            "data_dict_wells.pdf",
            "Drillinginfo Well dictionary.pdf",
            "reference_eur_calculation.pdf",
            "reference_recompletion.pdf",
        )
    ]
    to_write += (data / "bea").iterdir()
    to_write += (data / "eia").iterdir()
    to_write += [
        data / "epa" / f
        for f in ("EPA_GHG_report_2020_table_data.zip", "egrid2018_summary_tables.pdf")
    ]
    to_write += (data / "fred").iterdir()
    to_write += [
        data / "snl" / f
        for f in (
            "gas_prices_download_notes.txt",
            "na_gas_methodology.pdf",
            "snl_gas_price_series.xls",
            "snl_gas_volume_series.xls",
        )
    ]
    to_write += (data / "studies/alvarez_etal_2018").iterdir()
    to_write += (data / "studies/duren_etal_2019").iterdir()
    to_write += (data / "studies/frankenberg_etal_2016").iterdir()
    to_write += (data / "studies/lyon_etal_2016").iterdir()
    to_write += (data / "studies/omara_etal_2018").iterdir()
    to_write += (data / "studies/zavala-araiza_etal_2018").iterdir()

    # Code
    to_write += [
        Path(".") / f
        for f in (
            ".condarc",
            ".gitignore",
            "README.md",
            "Snakefile",
        )
    ]
    code_dir = Path("code")
    to_write += code_dir.glob("*.[rR]")
    to_write += code_dir.glob("*.py")
    to_write += code_dir.glob("*.sh")
    to_write += [
        code_dir / f
        for f in (
            "constants.json",
            "version_info_cmdstan.txt",
            "snl_nat_gas_basin_match.csv",
        )
    ]
    to_write += (code_dir / "envs").iterdir()
    to_write += (code_dir / "stan_models").glob("*.stan")  # exclude compiled models

    # Graphics (that aren't recreated)
    to_write += [
        Path("graphics") / f
        for f in (
            "game_tree_audit_covariates.tikz",
            "game_tree_audit_leak.tikz",
            "game_tree_audit_leak_censor.tikz",
            "marks_2018_fig4a.png",
            "natural_gas_leakage_percentages_marks_fig1.png",
            "ORCIDiD_iconvector.pdf",
            "frankenberg_etal_2016_fig1.jpg",
            "frankenberg_etal_2016_fig4-4.jpg",
        )
    ]

    # Output
    output_dir = Path("output")
    to_write += [
        output_dir / f
        for f in (
            "define_acronyms.tex",
            "methane_measurement_refs.bib",
            "paper.tex",
            "refs.bib",
            "software_cites_python.bib",
            "software_cites_r.bib",
            "tex_fragments/README.md",
        )
    ]
    to_write += (output_dir / "individual_figures").iterdir()

    with ZipFile(zip_filename, "w") as z:
        [z.write(f) for f in to_write]


if __name__ == "__main__":
    write_drillinginfo_data(snakemake.output["drillinginfo_zip"])
    write_public(snakemake.output["public_zip"])
