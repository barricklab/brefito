import glob
import os

brefito_config = {key.upper(): value for key, value in config.items()}
BREFITO_PACKAGE_PATH = brefito_config.get("BREFITO_PACKAGE_PATH", "")
TREE_UTILS = BREFITO_PACKAGE_PATH + "/workflow/scripts/tree_utils.pl"

# Detect samples from 01_breseq_initial_gd directory
_initial_gd_dir = "01_breseq_initial_gd"
samples = sorted([os.path.splitext(os.path.basename(f))[0]
                  for f in glob.glob(_initial_gd_dir + "/*.gd")])

# Ancestor GD files in working directory (e.g. Anc-_0gen_REL606.gd)
ancestor_files = sorted(glob.glob("Anc*.gd"))

# Reference and mask paths
REF_DIR = "references"
REF_GBK = REF_DIR + "/REL606.gbk"
PROPHAGE_GD = REF_DIR + "/prophage-amplifications.gd"

_nomask = config.get("NOMASK", False)
MASK_GD = REF_DIR + ("/empty.mask.gd" if _nomask else "/REL606.L20.G15.P0.M35.mask.gd")

_LTEE_REF_URL = "https://raw.githubusercontent.com/barricklab/LTEE-Ecoli/master/reference"

# Build target list
_all_targets = []
if samples:
    _all_targets += expand("04_final_normalized_gd/{sample}.gd", sample=samples)
    _all_targets += expand("mutated_genomes/{sample}.gff", sample=samples)
    _all_targets += ["compare_normalized.html", "initial.count.csv", "final.count.csv"]
    _all_targets += expand("05_normalized_masked_gd/{sample}.gd", sample=samples)
    _all_targets += expand("06_normalized_masked_no_IS_adjacent_gd/{sample}.gd", sample=samples)
    _all_targets += ["compare_normalized_masked.html", "final_masked.count.csv"]
    if ancestor_files:
        _all_targets += ["07_phylogeny/tree.rerooted.rescaled.tre"]


rule all_curate_ltee_clones:
    input: _all_targets
    default_target: True


# ── Download reference files from LTEE-Ecoli repository ─────────────────────

rule download_LTEE_ref_gbk:
    output: REF_GBK
    log:
        "logs/download-LTEE-ref-gbk.log"
    conda: "../envs/download.yml"
    shell:
        """
        mkdir -p {REF_DIR} > {log} 2>&1
        wget -q -O {output} "{_LTEE_REF_URL}/REL606.gbk" >> {log} 2>&1
        """

rule download_LTEE_mask_gd:
    output: REF_DIR + "/REL606.L20.G15.P0.M35.mask.gd"
    log:
        "logs/download-LTEE-mask-gd.log"
    conda: "../envs/download.yml"
    shell:
        """
        mkdir -p {REF_DIR} > {log} 2>&1
        wget -q -O {output} "{_LTEE_REF_URL}/REL606.L20.G15.P0.M35.mask.gd" >> {log} 2>&1
        """

rule download_LTEE_prophage_gd:
    output: PROPHAGE_GD
    log:
        "logs/download-LTEE-prophage-gd.log"
    conda: "../envs/download.yml"
    shell:
        """
        mkdir -p {REF_DIR} > {log} 2>&1
        wget -q -O {output} "{_LTEE_REF_URL}/prophage-amplifications.gd" >> {log} 2>&1
        """

rule download_LTEE_empty_mask_gd:
    output: REF_DIR + "/empty.mask.gd"
    log:
        "logs/download-LTEE-empty-mask-gd.log"
    conda: "../envs/download.yml"
    shell:
        """
        mkdir -p {REF_DIR} > {log} 2>&1
        wget -q -O {output} "{_LTEE_REF_URL}/empty.mask.gd" >> {log} 2>&1
        """


# ── Stub creation (ancient() ensures existing files are never overwritten) ───

rule create_LTEE_header:
    input:
        gd = ancient("01_breseq_initial_gd/{sample}.gd")
    output:
        "00_header/{sample}.gd"
    log:
        "logs/create-LTEE-header-{sample}.log"
    conda: "../envs/breseq_LTEE.yml"
    shell:
        """
        gdtools SUBTRACT -o {output}.tmp {input.gd} {input.gd} > {log} 2>&1
        gdtools NOT-EVIDENCE -o {output} {output}.tmp >> {log} 2>&1
        rm {output}.tmp
        """

rule create_LTEE_curate_add:
    input:
        gd = ancient("01_breseq_initial_gd/{sample}.gd")
    output:
        "02_curate_add/{sample}.gd"
    log:
        "logs/create-LTEE-curate-add-{sample}.log"
    conda: "../envs/breseq_LTEE.yml"
    shell:
        """
        gdtools SUBTRACT -o {output}.tmp {input.gd} {input.gd} > {log} 2>&1
        gdtools NOT-EVIDENCE -o {output} {output}.tmp >> {log} 2>&1
        rm {output}.tmp
        """

rule create_LTEE_curate_remove:
    input:
        gd = ancient("01_breseq_initial_gd/{sample}.gd")
    output:
        "02_curate_remove/{sample}.gd"
    log:
        "logs/create-LTEE-curate-remove-{sample}.log"
    conda: "../envs/breseq_LTEE.yml"
    shell:
        """
        gdtools SUBTRACT -o {output}.tmp {input.gd} {input.gd} > {log} 2>&1
        gdtools NOT-EVIDENCE -o {output} {output}.tmp >> {log} 2>&1
        rm {output}.tmp
        """


# ── Per-sample curation pipeline ─────────────────────────────────────────────

rule curate_LTEE_gd:
    input:
        initial    = "01_breseq_initial_gd/{sample}.gd",
        subtracts  = "02_curate_remove/{sample}.gd",
        add        = "02_curate_add/{sample}.gd",
        header     = "00_header/{sample}.gd",
        reference  = REF_GBK
    output:
        "03_curated/{sample}.gd"
    log:
        "logs/curate-LTEE-clones-{sample}.log"
    conda: "../envs/breseq_LTEE.yml"
    shell:
        """
        gdtools VALIDATE -r {input.reference} {input.initial} {input.add} {input.subtracts} {input.header} > {log} 2>&1
        gdtools SUBTRACT -o {output}.tmp1 {input.initial} {input.subtracts} >> {log} 2>&1
        gdtools UNION    -o {output}.tmp2 {input.add} {output}.tmp1 >> {log} 2>&1
        gdtools REHEADER -o {output} {input.header} {output}.tmp2 >> {log} 2>&1
        rm {output}.tmp1 {output}.tmp2
        """

rule normalize_LTEE_gd:
    input:
        gd        = "03_curated/{sample}.gd",
        reference = REF_GBK
    output:
        "04_final_normalized_gd/{sample}.gd"
    log:
        "logs/normalize-LTEE-clones-{sample}.log"
    conda: "../envs/breseq_LTEE.yml"
    shell:
        """
        gdtools NORMALIZE -a -r {input.reference} -o {output} {input.gd} > {log} 2>&1
        """

rule mask_LTEE_gd:
    input:
        gd        = "04_final_normalized_gd/{sample}.gd",
        reference = REF_GBK,
        prophage  = PROPHAGE_GD,
        mask      = MASK_GD,
        ancestors = ancestor_files
    output:
        "05_normalized_masked_gd/{sample}.gd"
    log:
        "logs/mask-LTEE-clones-{sample}.log"
    conda: "../envs/breseq_LTEE.yml"
    shell:
        """
        gdtools REMOVE   -c type==CON -o {output}.tmp1 {input.gd} > {log} 2>&1
        gdtools SUBTRACT -o {output}.tmp2 {output}.tmp1 {input.prophage} {input.ancestors} >> {log} 2>&1
        gdtools MASK     -v -s -o {output} {output}.tmp2 {input.mask} >> {log} 2>&1
        rm {output}.tmp1 {output}.tmp2
        """

rule no_is_adjacent_LTEE_gd:
    input:
        "05_normalized_masked_gd/{sample}.gd"
    output:
        "06_normalized_masked_no_IS_adjacent_gd/{sample}.gd"
    log:
        "logs/no-is-adjacent-LTEE-clones-{sample}.log"
    conda: "../envs/breseq_LTEE.yml"
    shell:
        """
        gdtools REMOVE -e -c TYPE!=UN -c adjacent!=UNDEFINED -o {output} {input} > {log} 2>&1
        """

rule apply_LTEE_mutations:
    input:
        gd        = "04_final_normalized_gd/{sample}.gd",
        reference = REF_GBK
    output:
        "mutated_genomes/{sample}.gff"
    log:
        "logs/apply-LTEE-mutations-{sample}.log"
    conda: "../envs/breseq_LTEE.yml"
    shell:
        """
        gdtools APPLY -r {input.reference} -f GFF3 -o {output} {input.gd} > {log} 2>&1
        """


# ── Aggregate analysis rules ─────────────────────────────────────────────────

rule compare_LTEE_normalized:
    input:
        gd_files  = expand("04_final_normalized_gd/{sample}.gd", sample=samples),
        reference = REF_GBK
    output:
        "compare_normalized.html"
    log:
        "logs/compare-LTEE-normalized.log"
    conda: "../envs/breseq_LTEE.yml"
    shell:
        """
        gdtools COMPARE -p -r {input.reference} -o {output} {input.gd_files} > {log} 2>&1
        """

rule compare_LTEE_normalized_masked:
    input:
        gd_files  = expand("05_normalized_masked_gd/{sample}.gd", sample=samples),
        reference = REF_GBK
    output:
        "compare_normalized_masked.html"
    log:
        "logs/compare-LTEE-normalized-masked.log"
    conda: "../envs/breseq_LTEE.yml"
    shell:
        """
        gdtools COMPARE -p -r {input.reference} -o {output} {input.gd_files} > {log} 2>&1
        """

rule count_LTEE_initial:
    input:
        gd_files  = expand("01_breseq_initial_gd/{sample}.gd", sample=samples),
        reference = REF_GBK
    output:
        "initial.count.csv"
    log:
        "logs/count-LTEE-initial.log"
    conda: "../envs/breseq_LTEE.yml"
    shell:
        """
        gdtools COUNT -o {output} -r {input.reference} {input.gd_files} > {log} 2>&1
        """

rule count_LTEE_final:
    input:
        gd_files  = expand("04_final_normalized_gd/{sample}.gd", sample=samples),
        reference = REF_GBK
    output:
        "final.count.csv"
    log:
        "logs/count-LTEE-final.log"
    conda: "../envs/breseq_LTEE.yml"
    shell:
        """
        gdtools COUNT -o {output} -r {input.reference} {input.gd_files} > {log} 2>&1
        """

rule count_LTEE_final_masked:
    input:
        gd_files  = expand("05_normalized_masked_gd/{sample}.gd", sample=samples),
        reference = REF_GBK
    output:
        "final_masked.count.csv"
    log:
        "logs/count-LTEE-final-masked.log"
    conda: "../envs/breseq_LTEE.yml"
    shell:
        """
        gdtools COUNT -o {output} -r {input.reference} {input.gd_files} > {log} 2>&1
        """

rule LTEE_phylogeny:
    input:
        sample_gds   = expand("05_normalized_masked_gd/{sample}.gd", sample=samples),
        ancestor_gds = ancestor_files,
        reference    = REF_GBK
    output:
        tre          = "07_phylogeny/tree.tre",
        genotypes    = "07_phylogeny/tree.genotypes.txt",
        sample_key   = "07_phylogeny/tree.sample.key.txt",
        mutation_key = "07_phylogeny/tree.mutation.key.txt",
        rerooted     = "07_phylogeny/tree.rerooted.tre",
        rescaled     = "07_phylogeny/tree.rerooted.rescaled.tre"
    log:
        "logs/LTEE-phylogeny.log"
    conda: "../envs/breseq_LTEE.yml"
    shell:
        """
        mkdir -p 07_phylogeny
        gdtools PHYLOGENY -p -a -o 07_phylogeny/tree -r {input.reference} {input.ancestor_gds} {input.sample_gds} > {log} 2>&1
        perl {TREE_UTILS} ROOT-ANCESTOR -i {output.tre} -o {output.rerooted} >> {log} 2>&1
        perl {TREE_UTILS} SCALE-PHYLIP  -i {output.rerooted} -o {output.rescaled} -p {output.genotypes} >> {log} 2>&1
        mkdir -p 07_phylogeny/discrepancies
        perl {TREE_UTILS} DISCREPANCIES -i {output.rerooted} -p 07_phylogeny/tree -o 07_phylogeny/discrepancies/tree >> {log} 2>&1
        """
