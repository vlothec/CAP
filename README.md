# CAP

1. /modules/TRASH_2/src/TRASH.R
Inputs:
	Required: *assembly.fasta* <- external input
	Optional: *templates.fasta* <- external input
optional TRASH2 arguments:
	-p --cores_no           number of cores for parallel run, default: 1
	-m --max_rep_size       maximum repeat size, default: 1000
	-i --min_rep_size       minimum repeat size, default: 7
	-t --templates          fasta file with repeat templates and their names 
First run: modules/TRASH_2/src/TRASH.R -o [output directory] -f [fasta input directory] [other optional arguments and optional input template]
Output: *[fasta_file]_repeats_with_seq.csv* *_arrays.csv*

2. /bin/filter_trash.R
Inputs: 
	Required: *[fasta_file]_repeats_with_seq.csv*, *_arrays.csv*
Output: *[fasta_file]_repeats_filtered.csv* *[fasta_file]_arrays_filtered.csv* 

3. /bin/merge_classes.R
Inputs: 
	Required: *[fasta_file]_repeats_filtered.csv*, *_arrays_filtered.csv*
Output: *[fasta_file]_repeats_filtered_reclassed.csv*, *_arrays_filtered_reclassed.csv*, *[fasta_file]_genome_classes.csv*

4. /bin/parse_TEs.R optional module
Inputs: 
	Required: *[fasta_file]_EDTA.TEanno.split.gff3* <- external input
Output: *[fasta_file]_TEs_parsed.csv*

5. /bin/filter_TEs.R optional module, run if 4 complete
Inputs: 
	Required: *[fasta_file]_TEs_parsed.csv*
Output: *[fasta_file]_TEs_filtered.csv*

6. /bin/parse_genes.R optional module
Inputs: 
	Required: *[fasta_file]_helixer.gff* <- external input
Output: *[fasta_file]_genes_parsed.csv*

7. /bin/filter_genes.R optional module, run if 6 complete
Inputs: 
	Required: *[fasta_file]_genes_parsed.csv*
Output: *[fasta_file]_genes_filtered.csv*

8. /bin/score_centromeric_classes.R 
Inputs: 
	Required: *[fasta_file]_repeats_filtered_reclassed.csv*, 
			  *_arrays_filtered_reclassed.csv*, 
			  *[fasta_file]_genome_classes.csv*
	Optional: *[fasta_file]_TEs_filtered.csv*, 
			  *[fasta_file]_genes_filtered.csv*
Output: *[fasta_file]_centromeric_scores.csv*

9. /bin/predict_centromeric.py
Inputs: 
	Required: *[fasta_file]_centromeric_scores.csv*, *[fasta_file]_metadata.csv*
	Additional model input: /model/centromeric_model_v2.pkl <- provided in the repository 
Output: *[fasta_file]_predictions.csv*

10. /bin/GC.R
Inputs: 
	Required: *assembly.fasta*
Output: *[fasta_file]_GC.csv*

11. /bin/CTW.R
Inputs: 
	Required: *assembly.fasta*
Output: *[fasta_file]_CTW.csv*

12. /bin/CAP.R 
Inputs:
	Required: *[fasta_file]_predictions.csv*, 
			  *[fasta_file]_repeats_filtered_reclassed.csv*, 
			  *[fasta_file]_arrays_filtered_reclassed.csv*, 
			  *[fasta_file]_genome_classes.csv*,
			  *[fasta_file]_metadata.csv*
	Optional: *[fasta_file]_TEs_filtered.csv*, 
			  *[fasta_file]_genes_filtered.csv*
Output: *[fasta_file]_CAP_plot.png*, 
		*[fasta_file]_CAP_repeat_families.csv*, 
		*[fasta_file]_CAP_dotplot.png*, 
		*[fasta_file]_CAP_model.txt*
		
		
11. /bin/get_metadata.R
Inputs:
	Required: *assembly.fasta*
Output: *[fasta_file]_metadata.csv*



sudo apt update && sudo apt install -y \
    build-essential \
    gfortran \
    libblas-dev liblapack-dev \
    zlib1g-dev \
	libicu-dev


R installation:
install.packages("renv")
install.packages("yaml")

# CAP: Centromere Analysis Pipeline

**CAP** is a reproducible, containerized Nextflow pipeline for **centromere prediction and repeat analysis** in genomic assemblies. It integrates **TRASH2** for repeat detection, processes transposon and gene annotations, and generates publication-ready plots and tables.

---

## Outputs

For input `genome.fasta`, CAP produces:

| File | Description |
|------|-------------|
| `genome_CAP_plot.png` | Main centromere visualization |
| `genome_CAP_dotplot.png` | Repeat array dotplot |

---

## Installation

Choose **one** of the three methods below.

---

### Option 1: Docker

```bash
# Pull and run
nextflow run yourname/cap-pipeline -profile docker --assembly data/genome.fasta
```


---

### Option 2: Conda

```bash
git clone https://github.com/yourname/CAP.git
cd CAP

# Create environment
conda env create -f environment.yml
conda activate cap-pipeline

# Run
nextflow run . -profile conda --assembly data/genome.fasta
```


---

### Option 3: R + renv

```bash
git clone https://github.com/yourname/CAP.git
cd CAP

# Install exact R packages
R -e 'renv::restore()'

# Run
nextflow run . -profile renv --assembly data/genome.fasta
```

---

## Usage

```bash
nextflow run . [options] --assembly <fasta>
```

### Required
```bash
--assembly PATH      Assembly in FASTA format
```

### Optional
```bash
--te_gff PATH        EDTA TE annotation (.gff3)
--gene_gff PATH      Helixer gene annotation (.gff)
--templates PATH     Custom repeat templates (FASTA)
--cores INT          Number of cores [default: 4]
--outdir DIR         Output directory [default: ./results]
```

### Example
```bash
nextflow run . -profile docker \
  --assembly data/arabidopsis.fasta \
  --te_gff data/EDTA.TEanno.split.gff3 \
  --gene_gff data/helixer.gff \
  --cores 8 \
  --outdir results_arabidopsis
```

---

## Directory Structure

```
CAP/
├── main.nf                  # Nextflow workflow
├── nextflow.config          # Profiles and params
├── environment.yml          # Conda dependencies
├── Dockerfile               # Docker image
├── renv.lock                # R package versions
├── bin/                     # R and Python scripts
├── modules/TRASH_2/         # TRASH2 (submodule)
├── model/                   # ML model
└── test/sample.fasta        # Test data
```

---

## HPC Usage (SLURM)

Create `conf/hpc.config`:
```groovy
process {
  executor = 'slurm'
  queue = 'normal'
  cpus = 8
  memory = '32 GB'
}
```

Run:
```bash
nextflow run . -profile docker,hpc --assembly genome.fasta
```

---

## Citation

If you use CAP, please cite:

---

## License

---

## Contact

---

