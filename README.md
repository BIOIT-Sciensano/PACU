[![PyPI version](https://badge.fury.io/py/pacu_snp.svg)](https://badge.fury.io/py/pacu_snp)
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat-square)](http://bioconda.github.io/recipes/pacu_snp/README.html)

# PACU
PACU is a workflow for whole genome sequencing based phylogeny of Illumina and ONT R9/R10 data.

PACU stands for the Prokaryotic Awesome variant Calling Utility and is named after an omnivorous fish (that eats both 
Illumina and ONT reads).

#### PACU is also available on our public [Galaxy instance](https://galaxy.sciensano.be/) (registration required), or on [UseGalaxy.eu](https://usegalaxy.eu/root?tool_id=toolshed.g2.bx.psu.edu/repos/iuc/pacu_snp/pacu_snp/0.0.5+galaxy0) (no registration required).

----

## INSTALLATION

### Pixi installation (recommended)

PACU can be installed easily using [Pixi](https://github.com/prefix-dev/pixi). 

```bash
mkdir pixi_pacu
cd pixi_pacu/
pixi init -c conda-forge -c bioconda
pixi add pacu_snp
pixi run PACU --help
``` 

### Conda installation

PACU can be installed in a new environment using the following command:
```bash
conda create -n pacu_snp -c conda-forge -c bioconda pacu_snp
```

**Note:** `MEGA` is currently not available through Conda, it can be installed manually from the link below, or 
`IQ-TREE` can be used instead (default).

### Manual installation

The PACU workflow has the following dependencies:
- [BEDTools](https://github.com/arq5x/bedtools2) (tested with 2.31.1)
- [bcftools](https://github.com/samtools/bcftools) (tested with 1.24)
- [FigTree](http://tree.bio.ed.ac.uk/software/figtree/) (tested with 1.4.4)
- [samtools](https://github.com/samtools/samtools) (tested with 1.24)
- [Gubbins](https://github.com/nickjcroucher/gubbins) (tested with 3.4.3)
- [snp-dists](https://github.com/tseemann/snp-dists) (tested with 1.2.0)
- [IQ-TREE 2](https://github.com/iqtree/iqtree2) (tested with 2.4.0, IQ-TREE 3 is not supported yet)
- [MEGA](https://www.megasoftware.net/) (optional, only required when using `--use-mega`)

The mapping script has the following additional dependencies:
- [Trimmomatic](https://github.com/usadellab/Trimmomatic) (tested with 0.40)
- [SeqKit](https://github.com/shenwei356/seqkit) (tested with 2.14.0)
- [Bowtie2](https://github.com/BenLangmead/bowtie2) (tested with 2.5.5)
- [Minimap2](https://github.com/lh3/minimap2) (tested with 2.31)

The corresponding binaries should be in your PATH to run the workflow. 
Other versions of these tools may work, but have not been tested.

The required Python packages are listed in the `pyproject.toml` file and are installed automatically by `pip`.
A `requirements.txt` file with pinned versions (generated with `pip-compile`) is also available.
PACU requires Python 3.10 or newer, but has only been tested with Python 3.10 (Gubbins is currently not available for 
newer Python versions on Bioconda).

```bash
python3.10 -m venv pacu_env
. pacu_env/bin/activate
pip install pacu_snp
```

## USAGE

```
usage: PACU [-h] [--ilmn-in ILMN_IN] [--ont-in ONT_IN] --ref-fasta REF_FASTA [--ref-bed REF_BED]
            [--dir-working DIR_WORKING] --output OUTPUT [--output-html OUTPUT_HTML] [--use-mega] [--include-ref]
            [--min-snp-af MIN_SNP_AF] [--min-snp-qual MIN_SNP_QUAL] [--min-snp-depth MIN_SNP_DEPTH]
            [--min-snp-dist MIN_SNP_DIST] [--skip-gubbins] [--min-global-depth MIN_GLOBAL_DEPTH]
            [--min-mq-depth MIN_MQ_DEPTH] [--bcftools-filt-af1] [--image-width IMAGE_WIDTH]
            [--image-height IMAGE_HEIGHT] [--threads THREADS] [--version]

options:
  -h, --help            show this help message and exit
  --ilmn-in ILMN_IN     Directory with Illumina input BAM files
  --ont-in ONT_IN       Directory with ONT input BAM files
  --ref-fasta REF_FASTA
                        Reference FASTA file
  --ref-bed REF_BED     BED file with phage regions
  --dir-working DIR_WORKING
                        Working directory
  --output OUTPUT       Output directory
  --output-html OUTPUT_HTML
                        Output report name
  --use-mega            If set, MEGA is used for the construction of the phylogeny (instead of IQ-TREE)
  --include-ref         If set, the reference genome is included in the phylogeny
  --min-snp-af MIN_SNP_AF
                        Minimum allele frequency for variants
  --min-snp-qual MIN_SNP_QUAL
                        Minimum SNP quality
  --min-snp-depth MIN_SNP_DEPTH
                        Minimum SNP depth
  --min-snp-dist MIN_SNP_DIST
                        Minimum distance between SNPs
  --skip-gubbins        If set, gubbins is skipped
  --min-global-depth MIN_GLOBAL_DEPTH
                        Minimum depth for all samples to include positions in SNP analysis
  --min-mq-depth MIN_MQ_DEPTH
                        MQ cutoff for samtools depth
  --bcftools-filt-af1   If enabled, allele frequency filtering also considers the VAF value
  --image-width IMAGE_WIDTH
                        Image width
  --image-height IMAGE_HEIGHT
                        Image height
  --threads THREADS
  --version             Print version and exit
```

**Note:** The location of the temporary directory can be changed by setting the `TMPDIR` environment variable.

### Basic usage example

The PACU workflow requires BAM files as input with reads mapped to a reference genome. 
Illumina data can be provided using the `--ilmn-in` option, ONT data can be provided using the `--ont-in` option.

- At least four input BAM files are required (for bootstrapping). The sample names are derived from the BAM file names.
- The reference genome should consist of a single sequence, unless Gubbins is disabled using `--skip-gubbins`.
- Phage regions (or other regions to exclude) can be provided in BED format using the `--ref-bed` option.

```
PACU \
    --ilmn-in in/ilmn/ \
    --ont-in in/ont/ \
    --ref-fasta ref.fasta \
    --output output/ \
    --dir-working work/ \
    --threads 8
```

### Read mapping

A script is included to map reads to a reference genome in FASTA format for both ONT and Illumina data.
The resulting BAM files can be used as input for the SNP workflow. The `--trim` option can be used to perform read
trimming before mapping.

The workflow generates a novel index for the reference genome if one is not already available.

*Illumina data*
```
PACU_map \
    --ref-fasta genome.fasta \
    --read-type illumina \
    --fastq-illumina reads_1.fastq.gz reads_2.fastq.gz \
    --output mapped.bam \
    --threads 4
```
*ONT data*
```
PACU_map \
    --ref-fasta genome.fasta \
    --read-type ont \
    --fastq-ont reads_ont.fastq.gz \
    --output mapped.bam \
    --threads 4
```

## TESTING

A test dataset is available under `pacu/resources/testdata/bam`, these files contain *Escherichia coli* reads mapped to
a small part of the *E. coli* NC_002695.2 genome. This is a not a real dataset, and should only be used for testing.

The tests are included in the package and can be executed for an existing installation (requires `pytest`):
```bash
pytest --pyargs pacu.tests
```

The complete workflow can be tested using the following command (from the repository root):
```bash
pytest --log-cli-level=DEBUG pacu/tests/test_workflow.py
```

**Note:** The MEGA tests fail if MEGA is not installed.

### Development environment

A development environment with all dependencies (and PACU installed in editable mode) can be created from the
repository using [Pixi](https://github.com/prefix-dev/pixi). The exact versions are pinned in `pixi.lock`.

```bash
pixi install
pixi run test
```

The `environment.yml` file is generated from `pixi.toml` and can be used to create the same environment with Conda
(`conda env create -f environment.yml`, from the repository root).

## CONTACT
[Create an issue](https://github.com/BIOIT-Sciensano/PACU/issues) to report bugs, propose new functions or ask for help.

## CITATION
If you use this tool, please consider citing our [publication](https://pubmed.ncbi.nlm.nih.gov/38441926/).

-----

Copyright - 2024-2026 Bert Bogaerts <bert.bogaerts@sciensano.be>
