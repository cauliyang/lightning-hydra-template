# Analysis Commands

## Table of Contents

1. [Arriba](#arriba)

## [Arriba](https://arriba.readthedocs.io/en/latest/)

- run_arriba_ly.sh

__NOTE__: Add output name compared origin script

```bash
#!/bin/bash

# For more information about this demo workflow, visit: https://arriba.readthedocs.io/en/latest/workflow/

if [ $# -lt 9 -o $# -gt 10 ]; then
	echo "Usage: $(basename $0) STAR_genomeDir/ annotation.gtf assembly.fa blacklist.tsv known_fusions.tsv protein_domains.gff3 threads outname read1.fastq.gz [read2.fastq.gz]" 1>&2
	exit 1
fi

# tell bash to be verbose and to abort on error
set -o pipefail
set -x -e -u

# get arguments
STAR_INDEX_DIR="$1"
ANNOTATION_GTF="$2"
ASSEMBLY_FA="$3"
BLACKLIST_TSV="$4"
KNOWN_FUSIONS_TSV="$5"
TAGS_TSV="$KNOWN_FUSIONS_TSV" # different files can be used for filtering and tagging, but the provided one can be used for both
PROTEIN_DOMAINS_GFF3="$6"
THREADS="$7"
OUTNAME="$8" ## here is my change
READ1="$9"
READ2="${10-}"

# find installation directory of arriba
BASE_DIR=$(dirname "$0")

# align FastQ files (STAR >=2.7.10a recommended)
STAR \
	--runThreadN "$THREADS" \
	--genomeDir "$STAR_INDEX_DIR" --genomeLoad NoSharedMemory \
	--readFilesIn "$READ1" "$READ2" --readFilesCommand zcat \
	--outStd BAM_Unsorted --outSAMtype BAM Unsorted --outSAMunmapped Within --outBAMcompression 0 \
	--outFilterMultimapNmax 50 --peOverlapNbasesMin 10 --alignSplicedMateMapLminOverLmate 0.5 --alignSJstitchMismatchNmax 5 -1 5 5 \
	--chimSegmentMin 10 --chimOutType WithinBAM HardClip --chimJunctionOverhangMin 10 --chimScoreDropMax 30 --chimScoreJunctionNonGTAG 0 --chimScoreSeparation 1 --chimSegmentReadGapMax 3 --chimMultimapNmax 50 |

tee "$OUTNAME".Aligned.out.bam |

# call arriba
"$BASE_DIR/arriba" \
	-x /dev/stdin \
	-o "$OUTNAME".fusions.tsv -O "$OUTNAME".fusions.discarded.tsv \
	-a "$ASSEMBLY_FA" -g "$ANNOTATION_GTF" -b "$BLACKLIST_TSV" -k "$KNOWN_FUSIONS_TSV" -t "$TAGS_TSV" -p "$PROTEIN_DOMAINS_GFF3" \
#	-d structural_variants_from_WGS.tsv

# sorting and indexing is only required for visualization
if [[ $(samtools --version-only 2> /dev/null) =~ ^1\. ]]; then
	samtools sort -@ "$THREADS" -m $((40000/THREADS))M -T tmp -O bam "$OUTNAME".Aligned.out.bam > "$OUTNAME".Aligned.sortedByCoord.out.bam
	rm -f "$OUTNAME".Aligned.out.bam
	samtools index "$OUTNAME".Aligned.sortedByCoord.out.bam
else
	echo "samtools >= 1.0 required for sorting of alignments" 1>&2

fi
```

- snakefile.py

```snakemake
from pathlib import Path


## TODO: Change this data folder
# DataFolder = "/panfs/home/yang4414/li002252/data/rna_ccle/success/paired_data/"
# DataFolder = "/panfs/home/yang4414/li002252/data/rna_ccle/success/test/"
DataFolder ="/projects/b1171/ylk4626/rna_ccle/test_arriba/"


SAMPLES = [ file.stem[:-5] for file in Path(DataFolder).iterdir() if file.suffix == ".gz"]
SAMPLES = list(set(SAMPLES))

ARRIBA_FILES= Path("/home/ylk4626/miniconda3/envs/arriba/var/lib/arriba")
arriba=f"/home/ylk4626/miniconda3/envs/arriba/bin/run_arriba_ly.sh"

black_list=ARRIBA_FILES/"blacklist_hg38_GRCh38_v2.3.0.tsv.gz"
star_index= "/projects/b1171/ylk4626/project/arriba/STAR_index_GRCh38viral_RefSeq_hg38"
known_fusions=ARRIBA_FILES/"known_fusions_hg38_GRCh38_v2.3.0.tsv.gz"
protein_domains=ARRIBA_FILES/"protein_domains_hg38_GRCh38_v2.3.0.gff3"

gtf="/projects/b1171/ylk4626/project/arriba/RefSeq_hg38.gtf"
fa="/projects/b1171/ylk4626/project/arriba/GRCh38viral.fa"


cmd = (f"{arriba} {star_index}/ {gtf} {fa} "
       f"{black_list} {known_fusions} {protein_domains} 8 " )

dst ="/projects/b1171/ylk4626/project/arriba/"


rule all:
    input:
        expand(dst + "{sample}/{sample}.fusions.tsv", sample=SAMPLES),
        expand(dst + "{sample}/{sample}.fusions.discarded.tsv", sample=SAMPLES)
        # expand(dst +"{sample}/{sample}.Aligned.sortedByCoord.out.bam.bai", sample=SAMPLES)

rule arriba:
    input:
        f1=DataFolder + "{sample}.1.fq.gz",
        f2=DataFolder + "{sample}.2.fq.gz"
    params:
        outdir=dst + '{sample}'
    output:
        dst + "{sample}/{sample}.fusions.tsv",
        dst + "{sample}/{sample}.fusions.discarded.tsv"
    resources:
        mem_mb="37G"
    shell:
        "rm -rf {params.outdir} && mkdir {params.outdir} && cd {params.outdir} && " + cmd + "{wildcards.sample} {input.f1} {input.f2} "

```

- run the pipeline

https://kb.northwestern.edu/page.php?id=69247

__Result__: /panfs/home/yang4414/li002252/project/arriba
