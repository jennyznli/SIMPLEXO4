import os

os.makedirs('logs/lsf',exist_ok=True)
os.makedirs('data/logs',exist_ok=True)
os.makedirs('data/bcftools',exist_ok=True)
os.makedirs('data/preprocess',exist_ok=True)

CHROMOSOMES_AUTOSOMAL=list(range(1,23))

BCF_INPUT=config['input']['chromosomes']

wildcard_constraints:
    CHR='[0-9]+'

rule select_variants_annotate:
    input:
        expand("data/bcftools/chr{CHR}.site-qc.bcf",CHR=CHROMOSOMES_AUTOSOMAL),
        expand("data/preprocess/chr{CHR}.annotation.pgen",CHR=CHROMOSOMES_AUTOSOMAL),
        expand("data/preprocess/chr{CHR}.annotation.pvar",CHR=CHROMOSOMES_AUTOSOMAL),
        expand("data/preprocess/chr{CHR}.annotation.psam",CHR=CHROMOSOMES_AUTOSOMAL),
        expand("data/preprocess/chr{CHR}.annotation.no_sample.vep.report.csv",CHR=CHROMOSOMES_AUTOSOMAL),
        "data/preprocess/all_chr.vep.report.csv"

rule bcftools_filter_per_chr:
    input:
        bcf=lambda wildcards: BCF_INPUT.format(CHR=wildcards.CHR)
    output:
        bcf="data/bcftools/chr{CHR}.site-qc.bcf",
        csi="data/bcftools/chr{CHR}.site-qc.bcf.csi"
    params:
        samples=config['input']['samples'],
        region_flag=("-R " + config['input']['targets']) if config.get('input',{}).get('targets') else ""
    log:
        "data/logs/chr{CHR}.bcftools_filter_per_chr.log"
    shell:
        """
        (bcftools view -S {params.samples} -a {params.region_flag} {input.bcf} |
        bcftools annotate --set-id '%CHROM\\_%POS\\_%REF\\_%ALT' |
        bcftools +fill-tags - -- -t VAF,AC,AC_Het,AC_Hom,AC_Hemi,AF,AN,NS,MAF,ExcHet,F_MISSING,HWE |
        bcftools filter -s LowCallRate -e 'INFO/F_MISSING > 0.05' -m + |
        bcftools view -f PASS -W=csi -Ob -o {output.bcf}) 2>&1 | tee {log}
        """

rule bcftools_no_sample_vcf:
    input:
        bcf="data/bcftools/chr{CHR}.site-qc.bcf"
    output:
        vcf="data/bcftools/chr{CHR}.site-qc.no_sample.vcf"
    log:
        "data/logs/chr{CHR}.bcftools_no_sample_vcf.log"
    shell:
        "(bcftools view -G -Ov -o {output.vcf} {input.bcf}) 2>&1 | tee {log}"

rule first_variant_annotation:
    input:
        "data/bcftools/chr{CHR}.site-qc.no_sample.vcf"
    output:
        "data/preprocess/chr{CHR}.annotation.no_sample.vep.vcf"
    threads: 16
    resources:
        mem_mb=48000
    log:
        "data/logs/chr{CHR}.first_variant_annotation.log"
    shell:
        """
        (export SINGULARITY_TMPDIR=/scratch/$USER/sing_tmp
        export SINGULARITY_CACHEDIR=/scratch/$USER/sing_cache

        mkdir -p /scratch/$USER/bin/chr{wildcards.CHR}
        cp /appl/samtools-1.21/bin/samtools /scratch/$USER/bin/chr{wildcards.CHR}/samtools
        chmod +x /scratch/$USER/bin/chr{wildcards.CHR}/samtools

        LOFTEE_REPO="/home/jennyzli/.vep/Plugins/loftee"
        LOFTEE_DIR="${{LOFTEE_REPO}}/GRCh38"

        singularity run --pwd "$PWD" -B "$PWD":"$PWD" -H "$PWD":"$PWD" \
        --bind /home/jennyzli/resources:/opt/vep/resources \
        --bind /home/jennyzli/.vep:/opt/vep/.vep \
        --bind /scratch/$USER/bin/chr{wildcards.CHR}/samtools:/usr/local/bin/samtools \
        --bind ${{LOFTEE_REPO}}:${{LOFTEE_REPO}} \
        --bind ${{LOFTEE_REPO}}/LoF.pm:/plugins/LoF.pm \
        --bind ${{LOFTEE_REPO}}/ancestral.pm:/plugins/ancestral.pm \
        --bind ${{LOFTEE_REPO}}/context.pm:/plugins/context.pm \
        --bind ${{LOFTEE_REPO}}/de_novo_donor.pl:/plugins/de_novo_donor.pl \
        --bind ${{LOFTEE_REPO}}/extended_splice.pl:/plugins/extended_splice.pl \
        --bind ${{LOFTEE_REPO}}/gerp_dist.pl:/plugins/gerp_dist.pl \
        --bind ${{LOFTEE_REPO}}/loftee_splice_utils.pl:/plugins/loftee_splice_utils.pl \
        --bind ${{LOFTEE_REPO}}/splice_site_scan.pl:/plugins/splice_site_scan.pl \
        --bind ${{LOFTEE_REPO}}/svm.pl:/plugins/svm.pl \
        --bind ${{LOFTEE_REPO}}/TissueExpression.pm:/plugins/TissueExpression.pm \
        --bind ${{LOFTEE_REPO}}/utr_splice.pl:/plugins/utr_splice.pl \
        --bind ${{LOFTEE_REPO}}/maxEntScan:/plugins/maxEntScan \
        /appl/containers/ensembl-vep_release_116.0.sif vep \
        --dir /opt/vep/.vep \
        -i $PWD/{input} \
        -o $PWD/{output} \
        --fork {threads} \
        --force_overwrite \
        --offline \
        --cache \
        --format vcf \
        --vcf --everything --canonical --mane \
        --assembly GRCh38 \
        --species homo_sapiens \
        --fasta /opt/vep/resources/Genomes/Human/hg38/fa/Homo_sapiens_assembly38.fasta \
        --vcf_info_field ANN \
        --plugin NMD \
        --plugin REVEL,/opt/vep/.vep/revel/revel_grch38.tsv.gz \
        --plugin SpliceAI,snv=/opt/vep/.vep/spliceai/spliceai_scores.raw.snv.hg38.vcf.gz,indel=/opt/vep/.vep/spliceai/spliceai_scores.raw.indel.hg38.vcf.gz \
        --plugin gnomADc,/opt/vep/.vep/gnomAD/gnomad.v3.1.1.hg38.genomes.gz \
        --plugin UTRAnnotator,/opt/vep/.vep/Plugins/UTRannotator/uORF_5UTR_GRCh38_PUBLIC.txt \
        --custom /opt/vep/.vep/clinvar/vcf_GRCh38/clinvar.autogvp.vcf.gz,ClinVar,vcf,exact,0,CLNSIG,CLNREVSTAT,CLNDN,AutoGVP \
        --plugin AlphaMissense,file=/opt/vep/.vep/alphamissense/AlphaMissense_GRCh38.tsv.gz \
        --plugin MaveDB,file=/opt/vep/.vep/mavedb/MaveDB_variants.tsv.gz \
        --plugin LoF,loftee_path:${{LOFTEE_REPO}},human_ancestor_fa:${{LOFTEE_DIR}}/human_ancestor.fa.gz,conservation_file:${{LOFTEE_DIR}}/loftee.sql,gerp_bigwig:${{LOFTEE_DIR}}/gerp_conservation_scores.homo_sapiens.GRCh38.bw \
        --cache_version 116) 2>&1 | tee {log}
        """

rule plink2_annotation_pgen:
    input:
        bcf="data/bcftools/chr{CHR}.site-qc.bcf",
        sex_file=config['input']['sex_info']
    output:
        pgen="data/preprocess/chr{CHR}.annotation.pgen",
        pvar="data/preprocess/chr{CHR}.annotation.pvar",
        psam="data/preprocess/chr{CHR}.annotation.psam"
    params:
        output_prefix="data/preprocess/chr{CHR}.annotation"
    threads: 16
    log:
        "data/logs/chr{CHR}.plink2_annotation_pgen.log"
    shell:
        """
        (plink2 --threads {threads} --bcf {input.bcf} --update-sex {input.sex_file} --double-id --vcf-half-call reference --make-pgen --out {params.output_prefix}) 2>&1 | tee {log}
        """

rule parse_annotation_no_sample_vep:
    input:
        "data/preprocess/chr{CHR}.annotation.no_sample.vep.vcf"
    output:
        "data/preprocess/chr{CHR}.annotation.no_sample.vep.report.csv"
    resources:
        mem_mb=32000,
        lsf_err="logs/lsf/parse_annotation_no_sample_vep.chr{CHR}.e",
        lsf_out="logs/lsf/parse_annotation_no_sample_vep.chr{CHR}.o"
    log:
        "data/logs/chr{CHR}.parse_annotation_no_sample_vep.log"
    shell:
        """
        (source $(conda info --base)/etc/profile.d/conda.sh
        conda activate vep_parser
        python scripts/vep_vcf_parser.py -i {input} -o {output} -m no_sample) 2>&1 | tee {log}
        """

rule combine_vep_reports:
    input:
        expand("data/preprocess/chr{CHR}.annotation.no_sample.vep.report.csv",CHR=CHROMOSOMES_AUTOSOMAL)
    output:
        "data/preprocess/all_chr.vep.report.csv"
    log:
        "data/logs/combine_vep_reports.log"
    shell:
        """
        (head -1 {input[0]} > {output}
        for f in {input}; do tail -n +2 "$f" >> {output}; done) 2>&1 | tee {log}
        """
