from dataclasses import dataclass
from pathlib import Path
import subprocess
import yaml
import gzip
import shutil

# Config
@dataclass
class Config:
    #Define what we're comparing
    species_1: str      
    species_2: str
    tissue: str

    #Define file ouput locations
    output_dir: Path    #Directory to store outputs
    temp_dir: Path      #Directory to store temporsry files and intermediate results --> May be deleted later

    #Input ATAC-seq file paths
    species_1_peak_file: Path
    species_2_peak_file: Path

    #HALPER mapping file paths
    species_1_to_species_2: Path | None
    species_2_to_species_1: Path | None

    #TSS file paths
    species_1_tss_file: Path
    species_2_tss_file: Path

    #Check configuration automatically after intialization
    def __post_init__(self):
        #Validate ATAC-seq files exist
        for p in [self.species_1_peak_file, self.species_2_peak_file]:
            ensure_exists(p, "Peak file")

        #Validate HALPER files exist (if provided)
        for h in [self.species_1_to_species_2, self.species_2_to_species_1]:
            if h is None:
                continue
            if not h.exists() and not (h.with_suffix(h.suffix + ".gz")).exists():
                raise FileNotFoundError(f"HALPER file missing: {h}")

        #Validate TSS files exist
        for t in [self.species_1_tss_file, self.species_2_tss_file]:
            ensure_exists(t, "TSS file")

        #Create output directory
        self.output_dir.mkdir(parents=True, exist_ok=True)



#HELPER FUNCTIONS

#ensure_exists() --> Checks that the file in the file path exists
def ensure_exists(path: Path, label: str):
    if not path.exists():
        raise FileNotFoundError(f"{label} does not exist: {path}")

#gunzip_keep() --> Unzip .gz file without deleting original gunzip file
def gunzip_keep(src: Path, dst: Path):
    with gzip.open(src, "rb") as f_in:
        with open(dst, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)

#ensure_unzipped --> Checks that the gunzip file has been properly unzipped by checking extension
def ensure_unzipped(path: Path) -> Path:
    if path.suffix != ".gz":
        return path

    unzipped = path.with_suffix("")
    if not unzipped.exists():
        gunzip_keep(path, unzipped)
        print(f"Unzipped {path} -> {unzipped}")
    return unzipped

#extract_bed3() --> Extracts the first 3 columns of a BED file
#Keeps only chrom, start, end
def extract_bed3(input_path: Path, output_path: Path):
    with open(output_path, "w") as out:
        subprocess.run(
            ["cut", "-f1-3", str(input_path)],
            stdout=out,
            check=True
        )


#load_config() --> Reads in a yaml file and intializes a config object
def load_config(config_path: Path, output_dir_key: str) -> Config:
    with open(config_path, "r") as f:
        cfg = yaml.safe_load(f)

    return Config(
        species_1=cfg["species_1"],
        species_2=cfg["species_2"],
        tissue=cfg["tissue"],

        output_dir=Path(cfg[output_dir_key]),
        temp_dir=Path(cfg["temp_dir"]),

        species_1_peak_file=Path(cfg["species_1_peak_file_cleaned"]),
        species_2_peak_file=Path(cfg["species_2_peak_file_cleaned"]),

        species_1_to_species_2=Path(cfg["species_1_to_species_2_cleaned"])
            if "species_1_to_species_2_cleaned" in cfg else None,
        species_2_to_species_1=Path(cfg["species_2_to_species_1_cleaned"])
            if "species_2_to_species_1_cleaned" in cfg else None,

        species_1_tss_file=Path(cfg["species_1_TSS_file"]),
        species_2_tss_file=Path(cfg["species_2_TSS_file"]),
    )


# Preprocessing
HALPER_KEYS = [
    "species_1_to_species_2",
    "species_2_to_species_1",
]

ATACSEQ_KEYS = [
    "species_1_peak_file",
    "species_2_peak_file",
]

"""
preprocess_config()
    Input: YAML file for the configuration

    Tasks:
        - Unzips HALPER files
        - Extracts first three columns from ATAC-SEQ + HALPER files
        - Writes new processed config

    Returns:
        Path to processed config
 """
def preprocess_config(config_path: Path) -> Path:

    with open(config_path, "r") as f:
        config = yaml.safe_load(f)

    cleaned_dir = Path(config["bedtool_preprocess_output_dir"])
    cleaned_dir.mkdir(parents=True, exist_ok=True)


    #Process HALPER files
    for key in HALPER_KEYS:
        if key not in config:
            continue

        halper_path = Path(config[key])

        #If HALPER files are unzipped/missing -> Try .gz extension
        if not halper_path.exists():
            gz_path = halper_path.with_name(halper_path.name + ".gz")
            if gz_path.exists():
                halper_path = gz_path
            else:
                raise FileNotFoundError(f"{halper_path} not found")

        #Ensure HALPER files are unzipped
        halper_path = ensure_unzipped(halper_path)

        # Extract first three columns of BED files
        cleaned_file = cleaned_dir / halper_path.name
        if not cleaned_file.exists():
            extract_bed3(halper_path, cleaned_file)

        config[key + "_cleaned"] = str(cleaned_file)


    # Process ATACSEQ files
    for key in ATACSEQ_KEYS:
        peak_path = Path(config[key])
        ensure_exists(peak_path, "Peak file")

        cleaned_file = cleaned_dir / peak_path.name
        if not cleaned_file.exists():
            extract_bed3(peak_path, cleaned_file)

        config[key + "_cleaned"] = str(cleaned_file)


    # Write new config
    processed_path = config_path.with_suffix(".processed.yaml")
    with open(processed_path, "w") as f:
        yaml.dump(config, f)

    print(f"Processed config written to: {processed_path}")

    return processed_path