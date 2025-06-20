import sys
import os
import pandas as pd

# Load configuration
configfile: "./config.json"

# Global variables from config
OUTPUT_DIR = config["general"]["output_dir"]
TOOLS_DIR = config["paths"]["tools_dir"]
SCRIPTS = config["paths"]["scripts"]
REFERENCE_GENOME_DIR = config["paths"]["reference_genome_dir"] # For convenience
PREPARED_GTF_FILE = config["paths"]["gtf_file"] # This is the output of prepare_reference

# Sample handling
SAMPLE_FILE_PATH = config["general"]["samples_tsv_path"]
SAMPLE_FILE = pd.read_table(SAMPLE_FILE_PATH, sep="\s+", dtype=str).set_index("sample", drop=False)
SAMPLE_LIST = SAMPLE_FILE["sample"].values.tolist()

# Helper function for constructing script paths
def script_path(script_name):
    return os.path.join(TOOLS_DIR, SCRIPTS[script_name])

# Rule: all - Specifies all final output files expected from the workflow.
rule all:
    input:
        OUTPUT_DIR + "/Final_analysis/diffExpr.P1e-3_C1.matrix",
        OUTPUT_DIR + "/featureCount.cnt_for_tpm.tpm.tab",
        directory(expand(OUTPUT_DIR + "/{rundeg_dir}", rundeg_dir=config["rundeg-output"])),
        # OUTPUT_DIR + "/featureCount.cnt.fixed", # Intermediate file
        expand(OUTPUT_DIR + "/{sample}/{sample}Aligned.sortedByCoord.out.bam", sample=SAMPLE_LIST),
        expand(OUTPUT_DIR + "/{sample}/TRIM/{sample}_pass_2_val_2.fq.gz", sample=SAMPLE_LIST),
        expand(OUTPUT_DIR + "/{sample}/TRIM/{sample}_pass_1_val_1.fq.gz", sample=SAMPLE_LIST),
        # Also ensure reference is prepared as it's a prerequisite for mapping and readscount
        os.path.join(REFERENCE_GENOME_DIR, "Genome"), # STAR Index sentinel
        PREPARED_GTF_FILE


# Rule: prepare_reference - Downloads and prepares reference genome, GTF, and STAR index.
rule prepare_reference:
    message: "-~- Preparing reference genome, GTF, and STAR index... -~-"
    output:
        star_index_sentinel=os.path.join(REFERENCE_GENOME_DIR, "Genome"), # STAR Index sentinel
        gtf=PREPARED_GTF_FILE, # Path from config, e.g., "reference/TAIR10_prepared.gtf"
        fixed_fasta=os.path.join(config["reference_preparation"]["base_dir"], config["reference_preparation"]["fixed_fasta_name"])
    log:
        os.path.join(OUTPUT_DIR, "logs", "reference_preparation", "prepare_reference.log")
    conda: "env/rnaseq.yml"
    params:
        # Directories
        base_prep_dir=config["reference_preparation"]["base_dir"], # e.g., "reference_sources"
        star_index_dir=REFERENCE_GENOME_DIR, # e.g., "reference"
        # File names
        raw_fasta_leafname=config["reference_preparation"]["raw_fasta_name"],
        # fixed_fasta_leafname is implicitly defined by output.fixed_fasta
        raw_gff_leafname=config["reference_preparation"]["raw_gff_name"],
        # URLs
        fasta_url=config["reference_preparation"]["fasta_url"],
        gff_url=config["reference_preparation"]["gff_url"],
        # STAR params
        threads=config["general"]["threads"],
        sjdb_overhang=config["reference_preparation"]["star_sjdb_overhang"],
        sa_index_nbases=config["reference_preparation"]["star_genome_sa_index_nbases"],
        genome_load=config["reference_preparation"]["star_genome_load"]
    shell:
        # Create directories
        "set -e; "
        "mkdir -p {params.base_prep_dir} && mkdir -p {params.star_index_dir} > {log} 2>&1 && "
        # Define full paths for clarity in shell
        "RAW_FASTA_PATH='{params.base_prep_dir}/{params.raw_fasta_leafname}' && "
        "FIXED_FASTA_PATH='{output.fixed_fasta}' && "
        "RAW_GFF_PATH='{params.base_prep_dir}/{params.raw_gff_leafname}' && "
        "PREPARED_GTF_PATH='{output.gtf}' && "
        # Download files
        "echo 'Downloading source files...' >> {log} 2>&1 && "
        "curl -L {params.gff_url} -o $RAW_GFF_PATH >> {log} 2>&1 && "
        "curl -L {params.fasta_url} -o $RAW_FASTA_PATH >> {log} 2>&1 && "
        # Fix FASTA headers
        "echo 'Fixing FASTA headers...' >> {log} 2>&1 && "
        "python -c \""
        "import sys; "
        "fixed_fasta_out_path = '$FIXED_FASTA_PATH'; "
        "with open('$RAW_FASTA_PATH', 'r') as fi, open(fixed_fasta_out_path, 'w') as fo: "
        "    for line in fi: "
        "        if line.startswith('>'): "
        "            header_parts = line[1:].split(None, 1); "
        "            seq_id = header_parts[0]; "
        # Try to make a robust Chr prefixing
        "            if seq_id.lower().startswith('chr'): "
        "                new_id = 'Chr' + seq_id[3:] if len(seq_id) > 3 else 'Chr' + seq_id[0].upper() + seq_id[1:]; "
        "            elif len(seq_id) == 1 and seq_id.isalpha(): "
        "                new_id = 'Chr' + seq_id.upper(); "
        "            elif seq_id.isdigit(): "
        "                new_id = 'Chr' + seq_id; "
        "            else: " # Default to just prefixing with Chr if uncertain
        "                new_id = 'Chr' + seq_id; "
        "            fo.write('>' + new_id + '\\n'); " # Escaped newline
        "        else: "
        "            fo.write(line)"
        "\" >> {log} 2>&1 && "
        # Convert GFF to GTF
        "echo 'Converting GFF to GTF...' >> {log} 2>&1 && "
        "gffread $RAW_GFF_PATH -T -o $PREPARED_GTF_PATH >> {log} 2>&1 && "
        # Generate STAR genome index
        "echo 'Generating STAR genome index...' >> {log} 2>&1 && "
        "STAR --runThreadN {params.threads} --genomeDir {params.star_index_dir} --genomeFastaFiles $FIXED_FASTA_PATH --sjdbGTFfile $PREPARED_GTF_PATH --runMode genomeGenerate --sjdbOverhang {params.sjdb_overhang} --genomeSAindexNbases {params.sa_index_nbases} --genomeLoad {params.genome_load} >> {log} 2>&1"

# Rule: fastqdump - Downloads FASTQ files from SRA.
rule fastqdump:
    message: "-~- Downloading fastq files for sample {wildcards.sample}... -~-"
    output:
        fwd=os.path.join(OUTPUT_DIR, "{sample}", "RAW", "{sample}_pass_1.fastq.gz"),
        rev=os.path.join(OUTPUT_DIR, "{sample}", "RAW", "{sample}_pass_2.fastq.gz")
    log:
        os.path.join(OUTPUT_DIR, "logs", "fastqdump", "{sample}.log")
    conda: "env/rnaseq.yml"
    params:
        outdir=lambda wildcards: os.path.join(OUTPUT_DIR, wildcards.sample, "RAW")
    shell:
        "set -e; fastq-dump --outdir {params.outdir} --gzip --skip-technical --readids --read-filter pass --dumpbase --split-3 --clip {wildcards.sample} > {log} 2>&1"

# Rule: trimming - Trims adapter sequences and low-quality bases from raw sequencing reads using Trim Galore.
rule trimming:
    message: "-~- Trimming fastq files for sample {wildcards.sample}... -~-"
    input:
        fwd=os.path.join(OUTPUT_DIR, "{sample}", "RAW", "{sample}_pass_1.fastq.gz"),
        rev=os.path.join(OUTPUT_DIR, "{sample}", "RAW", "{sample}_pass_2.fastq.gz")
    output:
        fwd=os.path.join(OUTPUT_DIR, "{sample}", "TRIM", "{sample}_pass_1_val_1.fq.gz"),
        rev=os.path.join(OUTPUT_DIR, "{sample}", "TRIM", "{sample}_pass_2_val_2.fq.gz")
    log:
        os.path.join(OUTPUT_DIR, "logs", "trimming", "{sample}.log")
    conda: "env/rnaseq.yml"
    threads: config["general"]["threads"]
    params:
        r1_clip=config["params"]["trimming"]["three_prime_clip_R1"],
        r2_clip=config["params"]["trimming"]["three_prime_clip_R2"],
        max_n=config["params"]["trimming"]["max_n"],
        outdir=lambda wildcards: os.path.join(OUTPUT_DIR, wildcards.sample, "TRIM")
    shell:
        'set -e; trim_galore \
            --paired \
            --three_prime_clip_R1 {params.r1_clip} \
            --three_prime_clip_R2 {params.r2_clip} \
            --cores {threads} \
            --max_n {params.max_n} \
            --gzip \
            -o {params.outdir} \
            {input.fwd} \
            {input.rev} > {log} 2>&1'

# Rule: mapping - Aligns trimmed reads to the reference genome using STAR.
rule mapping:
    message: "-~- Mapping reads for sample {wildcards.sample} with STAR... -~-"
    input:
        fwd=os.path.join(OUTPUT_DIR, "{sample}", "TRIM", "{sample}_pass_1_val_1.fq.gz"),
        rev=os.path.join(OUTPUT_DIR, "{sample}", "TRIM", "{sample}_pass_2_val_2.fq.gz"),
        # Ensure reference genome is prepared before mapping
        ref_genome_sentinel=os.path.join(REFERENCE_GENOME_DIR, "Genome")
    output:
        bam=os.path.join(OUTPUT_DIR, "{sample}", "{sample}Aligned.sortedByCoord.out.bam")
    log:
        os.path.join(OUTPUT_DIR, "logs", "mapping", "{sample}.log")
    conda: "env/rnaseq.yml"
    threads: config["general"]["threads"]
    params:
        prefix=lambda wildcards: os.path.join(OUTPUT_DIR, wildcards.sample, wildcards.sample),
        genome_dir=REFERENCE_GENOME_DIR, # from config["paths"]["reference_genome_dir"]
        multimap_max=config["params"]["mapping_star"]["outFilterMultimapNmax"],
        intron_min=config["params"]["mapping_star"]["alignIntronMin"],
        intron_max=config["params"]["mapping_star"]["alignIntronMax"]
    shell:
        'set -e; STAR --runMode alignReads \
            --runThreadN {threads} \
            --readFilesCommand zcat \
            --outFilterMultimapNmax {params.multimap_max} \
            --alignIntronMin {params.intron_min} \
            --alignIntronMax {params.intron_max} \
            --genomeDir {params.genome_dir} \
            --readFilesIn {input.fwd} {input.rev} \
            --outSAMtype BAM SortedByCoordinate \
            --outFileNamePrefix {params.prefix} > {log} 2>&1'

# Rule: readscount - Counts reads mapped to genomic features (genes) using featureCounts.
rule readscount:
    message: "-~- Counting reads with featureCounts... -~-"
    input:
        bams=expand(os.path.join(OUTPUT_DIR, "{sample}", "{sample}Aligned.sortedByCoord.out.bam"), sample=SAMPLE_LIST),
        # Ensure GTF is prepared before counting
        gtf_annotation=PREPARED_GTF_FILE
    output:
        raw_counts=os.path.join(OUTPUT_DIR, "featureCount.cnt"),
        fixed_counts=os.path.join(OUTPUT_DIR, "featureCount.cnt.fixed")
    log:
        os.path.join(OUTPUT_DIR, "logs", "readscount", "featureCounts.log")
    conda: "env/rnaseq.yml"
    threads: config["general"]["threads"]
    params:
        # Use the GTF input from prepare_reference rule
        gtf_to_use=lambda wildcards, input: input.gtf_annotation,
        quality=config["params"]["readscount_featureCounts"]["quality_threshold"],
        strandedness=config["params"]["readscount_featureCounts"]["strandedness"],
        script=script_path("fix_featCnt_header")
    shell:
        "set -e; featureCounts -p -Q {params.quality} -M --fraction -s {params.strandedness} -T {threads} -o {output.raw_counts} -a {params.gtf_to_use} {input.bams} > {log} 2>&1 && \
         python {params.script} {output.raw_counts} >> {log} 2>&1"

# Rule: rundeg - Performs differential gene expression analysis.
rule rundeg:
    message: "-~- Running differential expression analysis... -~-"
    input:
        ftcount=os.path.join(OUTPUT_DIR, "featureCount.cnt.fixed")
    output:
        rd_out=directory(os.path.join(OUTPUT_DIR, config["rundeg-output"]))
    log:
        os.path.join(OUTPUT_DIR, "logs", "rundeg", "rundeg.log")
    conda: "env/rnaseq.yml"
    params:
        script=script_path("run_DE_analysis"),
        method=config["params"]["rundeg"]["method"],
        samples_file=config["params"]["rundeg"]["samples_file"],
        contrasts_file=config["params"]["rundeg"]["contrasts_file"]
    shell:
        "set -e; {params.script} \
            -m {input.ftcount} \
            --method {params.method} \
            --samples_file {params.samples_file} \
            --output {output.rd_out} \
            --contrasts {params.contrasts_file} > {log} 2>&1"

# Rule: TPM - Calculates Transcripts Per Million (TPM) values.
rule TPM:
    message: "-~- Calculating TPM values... -~-"
    input:
        ftcount=os.path.join(OUTPUT_DIR, "featureCount.cnt.fixed")
    output:
        tpm_intermediate=os.path.join(OUTPUT_DIR, "featureCount.cnt_for_tpm"),
        tpm_final=os.path.join(OUTPUT_DIR, "featureCount.cnt_for_tpm.tpm.tab")
    log:
        os.path.join(OUTPUT_DIR, "logs", "TPM", "tpm_calculation.log")
    conda: "env/rnaseq.yml"
    params:
        script=script_path("tpm_calculator")
    shell:
        "set -e; cut -f 1,6- {input.ftcount} | egrep -v '#' > {output.tpm_intermediate} 2>> {log} && \
         python {params.script} -count {output.tpm_intermediate} >> {log} 2>&1"

# Rule: diffexpr - Analyzes differential expression results and generates final matrices.
rule diffexpr:
    message: "-~- Analyzing differential expression results and generating matrices... -~-"
    input:
        tpm_table=os.path.join(OUTPUT_DIR, "featureCount.cnt_for_tpm.tpm.tab"),
        rundeg_analysis_dir=directory(os.path.join(OUTPUT_DIR, config["rundeg-output"]))
    output:
        final_matrix=os.path.join(OUTPUT_DIR, "Final_analysis", "diffExpr.P1e-3_C1.matrix")
    log:
        os.path.join(OUTPUT_DIR, "logs", "diffexpr", "analyze_diff_expr.log")
    conda: "env/rnaseq.yml"
    params:
        script=script_path("analyze_diff_expr"),
        script_cwd=os.path.join(OUTPUT_DIR, config["rundeg-output"]),
        matrix_path=os.path.relpath(
            os.path.join(OUTPUT_DIR, "featureCount.cnt_for_tpm.tpm.tab"),
            os.path.join(OUTPUT_DIR, config["rundeg-output"])
        ),
        final_dest_dir=os.path.relpath(
            os.path.join(OUTPUT_DIR, "Final_analysis"),
            os.path.join(OUTPUT_DIR, config["rundeg-output"])
        ),
        p_value=config["params"]["diffexpr"]["p_value_threshold"],
        fold_change=config["params"]["diffexpr"]["fold_change_threshold"],
        max_genes=config["params"]["diffexpr"]["max_genes_clust"]
    shell:
        "set -e; \
         cd {params.script_cwd} && \
         {params.script} --matrix {params.matrix_path} -P {params.p_value} -C {params.fold_change} --max_genes_clust {params.max_genes} >> {log} 2>&1 && \
         mkdir -p {params.final_dest_dir} >> {log} 2>&1 && \
         mv diffExpr.P1e-3_C1* {params.final_dest_dir}/ >> {log} 2>&1 && \
         mv *.subset {params.final_dest_dir}/ >> {log} 2>&1 && \
         mv *.samples {params.final_dest_dir}/ >> {log} 2>&1 && \
         mv *.matrix {params.final_dest_dir}/ >> {log} 2>&1"
