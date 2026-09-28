# *Patiriella regularis* genome – code

Code for the genome assembly, decontamination and quality assessment of *Patiriella regularis*, the New Zealand common cushion star, known in te reo Māori as kapu parahua. Asteroidea: Valvatida: Asterinidae. Sequenced with PacBio Revio HiFi from a single individual, specimen S8, a uniformly blue animal. File names carry the identifier B1, which is the tube the extracted DNA was sent to the sequencing provider in.

> Citation to follow.

The sequencing data and assembly will be available from NCBI under the accessions given in the paper. This page contains the code only.

The analyses were run on the University of Otago SLURM cluster. In the code below, cluster-specific settings (job headers, module and conda environment loading, absolute paths) have been removed. Input and output files are given generic names, set as variables at the top of each block. Unless stated otherwise, programs were run with default settings.

Five helper scripts referenced below are in [`scripts/`](scripts/): `bin_bacteria.py`, `plot_bacteria.py`, `mito_atlas.py`, `karyotype.py` and `build_report.py`. All are plain Python with numpy and matplotlib as the only dependencies.

## Contents

- [Software](#software)
- [1. Read extraction and quality control](#1-read-extraction-and-quality-control)
  - [1.1 Checksum verification](#11-checksum-verification)
  - [1.2 BAM to FASTQ](#12-bam-to-fastq)
  - [1.3 Read quality control](#13-read-quality-control)
- [2. k-mer profiling](#2-k-mer-profiling)
  - [2.1 FastK](#21-fastk)
  - [2.2 GenomeScope2](#22-genomescope2)
  - [2.3 Smudgeplot](#23-smudgeplot)
- [3. Genome assembly](#3-genome-assembly)
- [4. Long-read scaffolding](#4-long-read-scaffolding)
- [5. Assembly quality assessment](#5-assembly-quality-assessment)
  - [5.1 Assembly statistics](#51-assembly-statistics)
  - [5.2 BUSCO](#52-busco)
  - [5.3 MerquryFK](#53-merquryfk)
  - [5.4 Read coverage](#54-read-coverage)
  - [5.5 Contig-level plots](#55-contig-level-plots)
- [6. Taxonomic screening and decontamination](#6-taxonomic-screening-and-decontamination)
  - [6.1 DIAMOND BLASTx](#61-diamond-blastx)
  - [6.2 BlobTools](#62-blobtools)
  - [6.3 Candidate identification](#63-candidate-identification)
  - [6.4 Kraken2 confirmation](#64-kraken2-confirmation)
  - [6.5 Removal](#65-removal)
- [7. Bacterial fraction](#7-bacterial-fraction)
  - [7.1 Binning by taxon](#71-binning-by-taxon)
  - [7.2 Completeness and composition](#72-completeness-and-composition)
- [8. Mitochondrial genome](#8-mitochondrial-genome)
  - [8.1 Identification](#81-identification)
  - [8.2 Reference selection](#82-reference-selection)
  - [8.3 Annotation](#83-annotation)
  - [8.4 tRNA annotation](#84-trna-annotation)
  - [8.5 Genome atlas](#85-genome-atlas)
- [9. Chromosome features](#9-chromosome-features)
- [10. Preparation for annotation](#10-preparation-for-annotation)
  - [10.1 Separating the mitochondrial genome](#101-separating-the-mitochondrial-genome)
  - [10.2 Renaming scaffolds](#102-renaming-scaffolds)
  - [10.3 Verification](#103-verification)
- [11. Assembly report](#11-assembly-report)
- [License](#license)

## Software

| Software | Version | Section |
|---|---|---|
| SAMtools | 1.21 | 1, 5, 6, 7 |
| pigz | – | 1 |
| SeqKit | – | 1, 5, 6, 7, 9 |
| NanoPlot | 1.47.2 | 1 |
| FastK / Histex | – | 2, 5 |
| GenomeScope | 2.0 | 2 |
| Smudgeplot | 0.5.3 | 2 |
| Hifiasm | 0.25.0-r726 | 3 |
| LongStitch (ntLink, ARKS/LINKS) | 1.0.5 | 4 |
| assembly-stats | – | 5 |
| BUSCO (metazoa_odb12, bacteria_odb12) | 6.1.0 | 5, 7 |
| BUSCO (metazoa_odb10) | 5.8.3 | 5 |
| MerquryFK | – | 5 |
| Minimap2 | 2.30-r1287 | 5, 6 |
| DIAMOND | 2.1.16 | 6 |
| BlobToolKit | – | 6 |
| Kraken2 (standard_16gb, 2025-04-02) | – | 6, 7 |
| MitoFinder | 1.4.1 | 8 |
| tRNAscan-SE | 2.0.13 | 8 |
| Python / NumPy / Matplotlib | 3.10 / – / – | 5, 7, 8, 9, 10 |

– version not recorded.

## 1. Read extraction and quality control

### 1.1 Checksum verification

The first transfer of the HiFi BAM to the cluster was corrupt: the file failed its md5 check and read extraction stopped partway through. The provider's drive and an intermediate copy both verified clean, so the file was transferred again. All subsequent steps used the verified copy, and a checksum gate was added to the extraction script so a corrupt input cannot propagate silently.

```bash
BAM=hifi_reads.bam

for f in "$BAM" "${BAM}.pbi"; do
    exp=$(awk '{print $1}' "${f}.md5.txt")
    act=$(md5sum "$f" | awk '{print $1}')
    [ "$exp" = "$act" ] && st=MATCH || st=MISMATCH
    printf '%-8s %s (%s bytes)\n' "$st" "$f" "$(stat -c %s "$f")"
done

# -u : unaligned input, do not require @SQ targets in the header
samtools quickcheck -u -v "$BAM"
```

The expected read count can be read from the PacBio `.pbi` index and compared against the number extracted.

```python
import gzip, struct, sys

with gzip.open(sys.argv[1], "rb") as f:   # hifi_reads.bam.pbi
    f.read(8)                             # magic + version
    struct.unpack("<H", f.read(2))        # flags
    print(struct.unpack("<I", f.read(4))[0])
```

### 1.2 BAM to FASTQ

```bash
BAM=hifi_reads.bam
READS=hifi_reads.fastq
THREADS=16

samtools fastq -@ "$THREADS" "$BAM" | pigz -p "$THREADS" > "${READS}.gz"
pigz -d -k -p "$THREADS" "${READS}.gz"     # FastK is faster on uncompressed input
```

### 1.3 Read quality control

```bash
BAM=hifi_reads.bam
READS=hifi_reads.fastq.gz
THREADS=16

seqkit stats -a -j "$THREADS" "$READS" > read_stats.txt

NanoPlot --ubam "$BAM" \
         -t "$THREADS" -o nanoplot -p preg_ \
         --N50 --loglength --tsv_stats --info_in_report \
         --format png --plots dot kde

# per-read length, predicted accuracy (rq) and pass count (np)
samtools view -@ "$THREADS" "$BAM" \
    | awk 'BEGIN {OFS = "\t"; print "len", "rq", "np"}
           { len = length($10); rq = "NA"; np = "NA"
             for (i = 12; i <= NF; i++) {
                 if ($i ~ /^rq:f:/) { split($i, a, ":"); rq = a[3] }
                 if ($i ~ /^np:i:/) { split($i, b, ":"); np = b[3] }
             }
             print len, rq, np }' \
    | pigz > read_metrics.tsv.gz
```

## 2. k-mer profiling

### 2.1 FastK

`-t1` retains singleton k-mers, so the error tail stays intact for GenomeScope's model fit.

```bash
READS=hifi_reads.fastq
K=31
THREADS=32
TMPDIR=/tmp/fastk

mkdir -p "$TMPDIR"
FastK -v -t1 -k"$K" -T"$THREADS" -M120 -P"$TMPDIR" -Npreg_k"$K" "$READS"
```

FastK writes `<prefix>.hist` itself. Histex output must therefore go to a different filename, or the shell truncates FastK's file before Histex reads it.

```bash
Histex -h1:10000 -G preg_k31 > preg_k31.genomescope.hist
```

### 2.2 GenomeScope2

Fitted under a diploid and a tetraploid model for comparison.

```bash
HIST=preg_k31.genomescope.hist
K=31

for P in 2 4; do
    genomescope2 -i "$HIST" -o "genomescope_p${P}" -k "$K" -p "$P" -n "preg_k${K}_p${P}"
done
```

### 2.3 Smudgeplot

`cutoff` derives the lower count threshold from the histogram. The coverage passed to `plot` is the `kcov` estimated by GenomeScope under the diploid model.

```bash
HIST=preg_k31.genomescope.hist
TABLE=preg_k31.ktab
COV=33
THREADS=16
TMPDIR=/tmp/smudge

mkdir -p "$TMPDIR"
L=$(smudgeplot cutoff "$HIST" L)

smudgeplot hetmers -L "$L" -t "$THREADS" -tmp "$TMPDIR" -o kmerpairs "$TABLE"
smudgeplot all -o preg -cov "$COV" --json_report kmerpairs.smu
```

## 3. Genome assembly

No polishing was applied. HiFi consensus is already above Q60, and polishing at that accuracy introduces about as many errors as it corrects; MerquryFK (section 5.3) measured QV 64.8 on the unpolished assembly.

```bash
READS=hifi_reads.fastq
PREFIX=preg
THREADS=48

hifiasm -o "$PREFIX" -t "$THREADS" --primary "$READS"

for gfa in "${PREFIX}.p_ctg.gfa" "${PREFIX}.a_ctg.gfa"; do
    awk '/^S/ {print ">"$2"\n"$3}' "$gfa" > "${gfa%.gfa}.fa"
done
```

## 4. Long-read scaffolding

The `ntLink-arks` target was used rather than `run`, so Tigmint-long did not break contigs. Gap filling was enabled, which fills joins with read sequence during ntLink rather than as a separate step. `G` is the haploid genome size from GenomeScope.

```bash
DRAFT=preg_primary          # expects preg_primary.fa
READS=preg_reads            # expects preg_reads.fq.gz
GENOME_SIZE=358000000
THREADS=32

longstitch ntLink-arks \
    draft="$DRAFT" \
    reads="$READS" \
    t="$THREADS" \
    G="$GENOME_SIZE" \
    longmap=hifi \
    gap_fill=True \
    rounds=2 \
    out_prefix=preg_longstitch
```

## 5. Assembly quality assessment

### 5.1 Assembly statistics

```bash
assembly-stats preg.p_ctg.fa preg.a_ctg.fa preg_longstitch.scaffolds.fa
seqkit stats -a preg_longstitch.scaffolds.fa
```

### 5.2 BUSCO

Run under BUSCO 6 with `odb12` for the comparison across stages, and under BUSCO 5 with `odb10` for comparability with older assemblies. The two use different gene sets (672 against 954) and are not directly comparable, so one should be reported consistently.

```bash
THREADS=32

for asm in contigs.fa scaffolds.fa decontaminated.fa; do
    busco -i "$asm" -o "$(basename "$asm" .fa)_metazoa" --out_path busco \
          -l metazoa_odb12 -m genome -c "$THREADS" -f
done

busco -i bacterial_scaffolds.fa -o bacterial --out_path busco \
      -l bacteria_odb12 -m genome -c "$THREADS" -f
```

### 5.3 MerquryFK

Run on the primary alone, on the primary and alternate together, and on the scaffolds. The primary alone recovers only one haplotype at each heterozygous site, so its k-mer completeness is expected to be low; the combined figure is the meaningful one.

```bash
READ_TABLE=preg_k31          # FastK table from section 2.1
THREADS=16

MerquryFK -v -T"$THREADS" -pdf "$READ_TABLE" primary.fa primary_only
MerquryFK -v -T"$THREADS" -pdf "$READ_TABLE" primary.fa alternate.fa both
MerquryFK -v -T"$THREADS" -pdf "$READ_TABLE" scaffolds.fa scaffolds
```

### 5.4 Read coverage

BlobTools requires a CSI index rather than BAI.

```bash
ASSEMBLY=scaffolds.fa
READS=hifi_reads.fastq
THREADS=32

minimap2 -ax map-hifi -t "$THREADS" --secondary=no "$ASSEMBLY" "$READS" \
    | samtools sort -@ 8 -o coverage.bam -
samtools index -c -@ 8 coverage.bam

samtools depth -a coverage.bam \
    | awk '{s[$1] += $3; n[$1]++} END {for (c in s) print c "\t" n[c] "\t" s[c]/n[c]}' \
    > per_scaffold_depth.tsv

samtools view -@ "$THREADS" coverage.bam \
    | awk '$3 != "*" {s[$3] += $5; n[$3]++} END {for (c in s) print c "\t" s[c]/n[c]}' \
    > per_scaffold_mapq.tsv
```

### 5.5 Contig-level plots

`contig_analysis.py` writes the contig size distribution, Nx curves, GC against length, an assembly comparison panel, and a trim-loss curve showing exactly how much sequence a given minimum-length filter would discard. `contig_classification.py` plots per-contig depth against length, coloured by mapping quality.

No length filter was applied to this assembly. The trim curve has no inflection point, so there is no natural cutoff separating junk from real sequence, and length is in any case a poor proxy for contamination: a 500 kb contaminant passes any sensible filter while a 16 kb mitochondrial genome fails it.

```bash
python3 scripts/contig_analysis.py primary.fa alternate.fa out_dir primary alternate
python3 scripts/contig_classification.py scaffolds.fa per_scaffold_depth.tsv per_scaffold_mapq.tsv out_dir
```

`contig_classification.py` calibrates against the length-weighted median depth of the assembly, not the haploid k-mer coverage from GenomeScope. Mapping all reads to a single collapsed assembly puts both haplotypes on the same sequence, so the modal depth is roughly twice `kcov`.

## 6. Taxonomic screening and decontamination

### 6.1 DIAMOND BLASTx

The output format is the one BlobTools expects: `qseqid staxids bitscore` followed by standard tabular fields.

```bash
ASSEMBLY=scaffolds.fa
DB=swissprot_tax.dmnd
THREADS=32

diamond blastx \
    --query "$ASSEMBLY" \
    --db "$DB" \
    --outfmt 6 qseqid staxids bitscore qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore \
    --sensitive \
    --max-target-seqs 1 \
    --evalue 1e-25 \
    --threads "$THREADS" \
    --out diamond.blastx.out
```

### 6.2 BlobTools

```bash
ASSEMBLY=scaffolds.fa
HITS=diamond.blastx.out
BAM=coverage.bam
BUSCO_TABLE=busco/scaffolds_metazoa/run_metazoa_odb12/full_table.tsv
TAXDUMP=taxdump

blobtools create --fasta "$ASSEMBLY" BlobDir
blobtools add --hits "$HITS" --taxrule bestsumorder --taxdump "$TAXDUMP" BlobDir
blobtools add --cov "$BAM" BlobDir
blobtools add --busco "$BUSCO_TABLE" BlobDir

blobtools view --plot --format png --view blob BlobDir
blobtools view --plot --format png --view snail BlobDir

blobtools filter --table blobtools_table.tsv \
                 --table-fields gc,length,cov,bestsumorder_phylum \
                 BlobDir
```

The phylum labels in the blob plot are strongly affected by database composition. SwissProt is dominated by human, mouse, *Drosophila* and *C. elegans* proteins and contains almost no echinoderm sequence, so conserved sea star genes best-hit vertebrate or arthropod orthologues and large host scaffolds carry labels such as Chordata. These labels were not used for filtering.

### 6.3 Candidate identification

Candidates were selected on coverage and mapping quality, not on length or taxonomy. Scaffolds below 10× mean depth with mapping quality above 40 form a coherent population: unique, confident mapping at a fraction of genome coverage is the signature of a separate organism sequenced at low depth.

```bash
CLASSIFICATION=contig_classification.tsv
awk -F'\t' 'NR > 1 && $3 < 10 && $4 > 40 {print $1}' "$CLASSIFICATION" > contaminant_ids.txt
```

GC alone was misleading here. The candidates sat at 25–50 % GC, mostly 35–38 %, below the host mean of 40.7 % rather than above it, which argues against the usual high-GC bacterial signature. Spirochaetes are AT-rich, which is why the GC signal pointed the wrong way. The DIAMOND hits settled it: gyrase A and B, aminoacyl-tRNA synthetases, elongation factor G, MutS and RNA polymerase, matching *Borrelia*, *Treponema*, *Geobacillus* and *Marinobacter* at 45–72 % identity.

A BlobTools filter on coverage and length was used as a cross-check. It flagged a superset of the manual selection; the extra scaffolds were assigned to Echinodermata or returned no protein hit, and were kept.

### 6.4 Kraken2 confirmation

```bash
DB=kraken2_standard_16gb
THREADS=16

kraken2 --db "$DB" --threads "$THREADS" \
        --report kraken2.report \
        --output kraken2.out \
        bacterial_scaffolds.fa
```

### 6.5 Removal

```bash
SCAFFOLDS=scaffolds.fa
IDS=contaminant_ids.txt

seqkit grep -f "$IDS" "$SCAFFOLDS" > bacterial_scaffolds.fa
seqkit grep -v -f "$IDS" "$SCAFFOLDS" > decontaminated.fa
seqkit stats -a "$SCAFFOLDS" decontaminated.fa bacterial_scaffolds.fa
```

To rebuild a BlobDir on the decontaminated assembly the coverage BAM must be subset as well. Reheadering alone is not sufficient, because the alignment records keep the original reference numbering; the BAM has to be rewritten through SAM so the reference IDs are renumbered.

```bash
samtools view -H coverage.bam | grep '^@SQ' | sed 's/.*SN:\([^\t]*\).*/\1/' > all_refs.txt
grep -v -F -f contaminant_ids.txt all_refs.txt > keep_refs.txt

samtools view -h -@ 8 coverage.bam $(tr '\n' ' ' < keep_refs.txt) \
    | grep -v -F -f contaminant_ids.txt \
    | samtools view -b -@ 8 -o coverage.decontam.bam -
samtools index -c -@ 8 coverage.decontam.bam
```

## 7. Bacterial fraction

The removed scaffolds were retained and assessed rather than discarded.

### 7.1 Binning by taxon

```bash
python3 scripts/bin_bacteria.py kraken2.out kraken2.report bacterial_scaffolds.fa bins_by_taxon
```

### 7.2 Completeness and composition

```bash
python3 scripts/plot_bacteria.py kraken2.report plots
```

Reassembly was not attempted. Only 441 reads (7.52 Mb) map to these scaffolds, about 2.7× coverage against the implied genome sizes. Autocycler's documented lower bound is 25×, with 50× as a working target, so no assembler can recover genomes from this depth. The scaffolds assign across many genera with roughly one scaffold each, and BUSCO against `bacteria_odb12` gives 44.0 % complete, 16.4 % fragmented and 39.7 % missing: fragments of a community, not recovered genomes.

## 8. Mitochondrial genome

### 8.1 Identification

One scaffold stood out on every metric: 16,388 bp at 659.6× mean depth (about 12× nuclear), mapping quality 60, GC 38.3 %, with a best DIAMOND hit to NADH-ubiquinone oxidoreductase chain 5 of *Patiria pectinifera* at E = 8.17e-224.

```bash
seqkit grep -f mito_id.txt scaffolds.fa > mitogenome.fa
seqkit fx2tab -nlg mitogenome.fa
```

Circularity was tested by comparing the terminal kilobases against each other. No overlap was found, so the sequence is linear as assembled.

```bash
seqkit subseq -r 1:1000    mitogenome.fa > mito_start.fa
seqkit subseq -r -1000:-1  mitogenome.fa > mito_end.fa
blastn -query mito_start.fa -subject mito_end.fa \
       -outfmt '6 qstart qend sstart send pident length' -evalue 1e-10
```

### 8.2 Reference selection

Annotation was run twice. The first reference, *Coscinasterias muricata* (PV535852.1, 16,178 bp), was chosen on length and New Zealand provenance before the sample was identified to species. It is Asteriidae, order Forcipulatida, whereas *Patiriella regularis* is Asterinidae, order Valvatida. The second reference, *Patiria pectinifera* (NC_001627.1, 16,260 bp), is an asterinid, the same genus as the best protein hit for the contig, and carries a complete curated feature table of 13 CDS, 22 tRNAs and 2 rRNAs.

Reference distance changed the result. Nothing else differed between the two runs.

| Feature | *Coscinasterias* (Forcipulatida) | *Patiria* (Valvatida) |
|---|---|---|
| Protein-coding genes | 12 of 13 | 13 of 13 |
| rRNA genes | 1 of 2 | 2 of 2 |
| tRNA genes | 19 of 22 | 19 of 22 |

ATP8 is 93 bp, the shortest and most divergent mitochondrial protein, and at that length a reference from another order carries too little similarity to place it. The candidate search below is restricted to Asterinidae; the table it prints should be read for a record whose title says "complete genome" rather than one flagged `UNVERIFIED:`, since MitoFinder transfers the reference's annotation and inherits its errors.

```bash
ACCESSION=NC_001627.1

curl -sL "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=nuccore&term=Asterinidae%5BOrganism%5D+AND+complete+genome%5BTitle%5D+AND+mitochondrion%5BTitle%5D&retmax=20&retmode=json"

curl -sL "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=${ACCESSION}&rettype=gb&retmode=text" \
    > reference_mito.gb
```

### 8.3 Annotation

Genetic code 9 is the echinoderm and flatworm mitochondrial code. Using the invertebrate or vertebrate code produces wrong translations.

```bash
MITO=mitogenome.fa
REF=reference_mito.gb
THREADS=8

mitofinder -j preg_mito \
           -a "$MITO" \
           -r "$REF" \
           -o 9 \
           -p "$THREADS" \
           --out-gb --new-genes \
           -t trnascan
```

The bioconda package of MitoFinder 1.4.1 does not ship `Mitofinder.config` or the `install.sh.ok` marker, and aborts without them. The config holds folder paths to the bundled tools, all of which the conda package places in the environment's `bin`:

```bash
ENV_BIN=/path/to/env/bin

touch "$ENV_BIN/install.sh.ok"
for key in megahitfolder blastfolder idbafolder metaspadesfolder \
           arwenfolder trnascanfolder mitfifolder; do
    echo "$key = $ENV_BIN/"
done > "$ENV_BIN/Mitofinder.config"
```

### 8.4 tRNA annotation

MitoFinder's bundled tRNAscan-SE call fails silently: it looks for the codon table at a relative path (`lib/gcode/`) rather than where the conda package installs it, reports "tRNA annotation run well", and writes no tRNAs. tRNAscan-SE was therefore run directly in organellar mode with the echinoderm codon table given explicitly.

```bash
MITO=mitogenome.fa
GCODE=/path/to/env/lib/tRNAscan-SE/gcode/gcode.echdmito

tRNAscan-SE -O -g "$GCODE" -o trnascan.out -f trnascan.ss "$MITO"
```

The codon table matters: without it the tRNA at the Trp position is called `Sup` (suppressor), because TCA reads as a stop under the standard code. Under the echinoderm code TCA is tryptophan, and the call corrects to Trp — an independent confirmation that code 9 is correct for this genome.

### 8.5 Genome atlas

`mito_atlas.py` draws a circular atlas following the usual convention, outside in: ruler, forward-strand genes, reverse-strand genes, GC content as deviation from the genome mean, and GC skew as a two-colour deviation from zero. It also writes a linear map and a feature table. The optional fifth argument merges tRNA calls from a tRNAscan-SE output file, which is needed because MitoFinder does not write them into its GFF.

```bash
python3 scripts/mito_atlas.py mitogenome.fa \
        preg_mito/*_Final_Results/*_mtDNA_contig.gff \
        plots \
        "Patiriella regularis mitochondrial genome" \
        trnascan.out
```

Against the asterinid reference the annotation recovered 13 of 13 protein-coding genes, 2 of 2 rRNAs and 19 of 22 tRNAs. Gene order matches the echinoderm arrangement, including the cluster of fourteen tRNAs between ND1 and COX1 that is diagnostic for the group.

Two limits remain. `rrnL` is called at 448 to 914, only 467 bp, where a metazoan 16S rRNA runs 1,300 to 1,600 bp; the unannotated 660 bp immediately following it is almost certainly the remainder of the gene, so that boundary is a partial call. The three absent tRNAs are Asp, Lys and the second Ser, the last being a recognised difficulty because its secondary structure is aberrant in echinoderms. Both should be curated by hand before submission.

## 9. Chromosome features

A cursory scan for telomeric repeat arrays and tandem satellite arrays, run to establish how much of the assembly is already chromosome-complete and what further scaffolding would stand to gain.

Telomeres are called by counting the canonical animal repeat `TTAGGG` and its reverse complement in the terminal 20 kb of each scaffold; twenty or more copies calls a telomere. Satellite arrays are found by walking each scaffold in 50 kb windows and recording the fraction of positions holding a distinct 21-mer. Ordinary sequence sits near 1.0 because almost every k-mer is unique, while a tandem array of period *n* contains only about *n* distinct k-mers however long it runs, so the fraction collapses toward zero.

```bash
ASSEMBLY=decontaminated.fa
MIN_LENGTH=1000000
MOTIF=TTAGGG

python3 scripts/karyotype.py "$ASSEMBLY" karyotype "$MIN_LENGTH" "$MOTIF"
```

The satellite threshold (`SAT_CUT` at the top of the script) was set to 0.05. The window distribution is strongly bimodal — median 0.985, fifth percentile 0.062, minimum 0.012 — and a cut at 0.5 flagged 13 % of windows, which in a genome that is 32 % repetitive is mostly ordinary repeat-rich sequence rather than satellite.

These satellite arrays are not confirmed centromeres. A centromere is defined by what binds it, not by its sequence, and cannot be identified from an assembly alone. In most eukaryotes centromeres sit within large satellite arrays, so these regions are reasonable candidates, but confirming one requires Hi-C contact data, CENP-A ChIP, or a related reference genome.

## 10. Preparation for annotation

Repeat and gene annotation are run on a renamed, mitochondrion-free copy of the decontaminated assembly. Renaming happens once, before any annotation, because every downstream coordinate file inherits the sequence identifiers and a later rename breaks them.

### 10.1 Separating the mitochondrial genome

The mitochondrial scaffold is removed before annotation. It uses genetic code 9 and BRAKER would apply the nuclear code to it, and its copy number is two orders of magnitude above the nuclear average, which distorts repeat modelling.

```bash
ASSEMBLY=decontaminated.fa
MITO_ID=mito_scaffold_id

samtools faidx "$ASSEMBLY"
samtools faidx "$ASSEMBLY" "$MITO_ID" > mitochondrion.fa
cut -f1 "${ASSEMBLY}.fai" | grep -v -x "$MITO_ID" > nuclear.ids
seqkit grep -f nuclear.ids "$ASSEMBLY" > nuclear_unsorted.fa
```

### 10.2 Renaming scaffolds

Scaffolds are sorted longest first and given sequential identifiers. Identifiers are kept under 50 characters and restricted to letters, digits and underscores, which is what NCBI accepts and what BRAKER, RepeatMasker and GFF tools handle without quoting problems. Descriptions are dropped, as trailing text after the first whitespace is silently truncated by some tools and retained by others.

```bash
IN=nuclear_unsorted.fa
PREFIX=Preg_scaf
OUT=Patiriella_regularis_B1_v1.nuclear.fa

seqkit sort -lr "$IN" \
  | seqkit replace -p '.*' -r "${PREFIX}{nr}" --nr-width 3 \
  | seqkit seq -w 60 > "$OUT"

seqkit fx2tab -nl "$OUT" | head
```

The mapping from old to new identifiers is written out and kept, because it is the only way to relate the annotation back to the assembly graph or to any earlier coordinate file.

```bash
paste <(seqkit sort -lr "$IN" | seqkit fx2tab -n) \
      <(seqkit fx2tab -n "$OUT") > scaffold_name_map.tsv
```

### 10.3 Verification

Renaming is verified by length, not by name: the sorted input and the renamed output are compared position by position, and any difference in length means the order changed.

```bash
paste <(seqkit sort -lr "$IN" | seqkit fx2tab -nl | awk '{print $1"\t"$NF}') \
      <(seqkit fx2tab -nl "$OUT" | awk '{print $1"\t"$NF}') \
  | awk '$2!=$4 {print "MISMATCH:", $0; n++} END{print (n?n:0), "mismatched lengths of", FNR}'
```

`0 mismatched lengths of N`, where N is the scaffold count, is the expected result.

## 11. Assembly report

`build_report.py` collects the outputs of every stage above into one self-contained HTML file: figures are base64-embedded, only system fonts are used, and nothing is loaded from a network, so the file can be moved or emailed on its own. Any missing input is marked in the report rather than causing a failure, and the script prints how many were not found.

```bash
python3 scripts/build_report.py [PROJECT_DIR] [OUTPUT.html]
```

MerquryFK writes PDFs; the script rasterises them with `pdftoppm` or ImageMagick if either is available.

## License

MIT License

Copyright (c) 2026 Marc A. Bailie and Nathan J. Kenny

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
