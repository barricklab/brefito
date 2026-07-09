#sample info not used but this loads the BRESEQ_ENV variable as well
# This workflow derives its samples from the gd/ directory, so it does not
# require sample information from data.csv / the read/assembly directories.
SAMPLE_INFO_REQUIRED = False
try: sample_info
except NameError:
    include: "load-sample-info.smk"

import glob
import os
import re
import sys

brefito_config = {key.upper(): value for key, value in config.items()}
BREFITO_PACKAGE_PATH = brefito_config.get("BREFITO_PACKAGE_PATH", "")
TREE_UTILS = BREFITO_PACKAGE_PATH + "/workflow/scripts/tree_utils.pl"

# Detect samples from gd directory
_initial_gd_dir = "gd"
samples = sorted([os.path.splitext(os.path.basename(f))[0]
                  for f in glob.glob(_initial_gd_dir + "/*.gd")])

# If specific sample names were given on the brefito command line, restrict the
# run to just building the curated mutant genome (mutants/{sample}.gff) for those
# samples. In this mode we skip all aggregate steps (count, compare, phylogeny).
_requested_samples = None
if 'SAMPLES' in brefito_config:
    _requested_samples = brefito_config['SAMPLES'].split("_,_")
    samples = [s for s in samples if s in _requested_samples]

# The 00_header / 01_curate_add / 02_curate_remove files are created (and then
# hand-curated) by the separate 'curate-LTEE-clones-init' workflow. This workflow
# consumes them as pre-existing inputs and never regenerates them, so that a
# --forceall run cannot clobber the hand-edited curation files. Fail early with a
# clear message if any are missing.
_init_dirs = ["00_header", "01_curate_add", "02_curate_remove"]
_missing_init_files = [
    os.path.join(d, s + ".gd")
    for s in samples
    for d in _init_dirs
    if not os.path.exists(os.path.join(d, s + ".gd"))
]
if _missing_init_files:
    sys.stderr.write(
        "Error: The following curation input files are missing:\n"
        + "".join("  " + f + "\n" for f in _missing_init_files)
        + "These files are created by the 'curate-LTEE-clones-init' workflow.\n"
        + "Please run:\n"
        + "    brefito curate-LTEE-clones-init\n"
        + "first (then edit the 01_curate_add/ and 02_curate_remove/ files as needed) "
        + "before running curate-LTEE-clones.\n"
    )
    sys.exit(1)

# Ancestor GD files in working directory (e.g. Anc-_0gen_REL606.gd)
ancestor_files = sorted(glob.glob("Anc*.gd"))

# If no ancestor GD file is present, try to determine which one to download
# automatically from the barricklab/LTEE-Ecoli repository based on the name
# of the directory brefito is being run in.
_LTEE_ANCESTOR_BASE_URL = "https://raw.githubusercontent.com/barricklab/LTEE-Ecoli/master"
_ANCESTOR_DOWNLOAD_FILE = None
_ANCESTOR_DOWNLOAD_URL = None

if samples and not ancestor_files and _requested_samples is None:
    _cwd_name = os.path.basename(os.getcwd())
    if re.search(r"(Ara-|A-)\d+", _cwd_name):
        _ANCESTOR_DOWNLOAD_FILE = "Anc-_0gen_REL606.gd"
        _ANCESTOR_DOWNLOAD_URL = _LTEE_ANCESTOR_BASE_URL + "/LTEE-clone-curated/" + _ANCESTOR_DOWNLOAD_FILE
    elif re.search(r"(Ara\+|A\+)\d+", _cwd_name):
        _ANCESTOR_DOWNLOAD_FILE = "Anc+_0gen_REL607.gd"
        _ANCESTOR_DOWNLOAD_URL = _LTEE_ANCESTOR_BASE_URL + "/LTEE-clone-curated/" + _ANCESTOR_DOWNLOAD_FILE
    elif re.search(r"MAE|MA", _cwd_name):
        _ANCESTOR_DOWNLOAD_FILE = "Anc+_REL1207.gd"
        _ANCESTOR_DOWNLOAD_URL = _LTEE_ANCESTOR_BASE_URL + "/MAE-clone-curated/" + _ANCESTOR_DOWNLOAD_FILE
    else:
        sys.stderr.write(
            "Error: No Anc*.gd ancestor genome diff file was found in this directory, "
            "and the directory name '" + _cwd_name + "' does not match a recognized "
            "LTEE naming convention (Ara-#, A-#, Ara+#, A+#, MA, MAE), so an ancestor "
            "genome diff cannot be downloaded automatically.\n"
            "Please place the appropriate ancestor Anc*.gd file in this directory "
            "before running curate-LTEE-clone.\n"
        )
        sys.exit(1)

    ancestor_files = [_ANCESTOR_DOWNLOAD_FILE]

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
    if _requested_samples is not None:
        # Only build the mutant genome for the requested samples. Snakemake pulls in
        # the full header/curate/normalize/apply chain automatically; the counting,
        # comparison, and phylogeny steps are intentionally omitted.
        _all_targets += expand("mutants/{sample}.gff", sample=samples)
    else:
        _all_targets += expand("04_final_normalized_gd/{sample}.gd", sample=samples)
        _all_targets += expand("mutants/{sample}.gff", sample=samples)
        _all_targets += ["output/compare_normalized.html", "output/count.initial.csv", "output/count.final.csv"]
        _all_targets += expand("05_normalized_masked_gd/{sample}.gd", sample=samples)
        _all_targets += expand("06_normalized_masked_no_IS_adjacent_gd/{sample}.gd", sample=samples)
        _all_targets += ["output/compare_normalized_masked.html", "output/count.final_masked.csv"]
        if ancestor_files:
            _all_targets += ["output/final.tre"]
            _all_targets += ["output/compare_normalized_masked.discrepancies.html"]


rule all_curate_ltee_clones:
    input: _all_targets
    default_target: True


# ── Download reference files from LTEE-Ecoli repository ─────────────────────

rule download_LTEE_ref_gbk:
    output: REF_GBK
    log:
        "logs/download/LTEE-ref-gbk.log"
    conda: "../envs/download.yml"
    shell:
        """
        mkdir -p {REF_DIR} > {log} 2>&1
        wget -q -O {output} "{_LTEE_REF_URL}/REL606.gbk" >> {log} 2>&1
        """

rule download_LTEE_mask_gd:
    output: REF_DIR + "/REL606.L20.G15.P0.M35.mask.gd"
    log:
        "logs/download/LTEE-mask-gd.log"
    conda: "../envs/download.yml"
    shell:
        """
        mkdir -p {REF_DIR} > {log} 2>&1
        wget -q -O {output} "{_LTEE_REF_URL}/REL606.L20.G15.P0.M35.mask.gd" >> {log} 2>&1
        """

rule download_LTEE_prophage_gd:
    output: PROPHAGE_GD
    log:
        "logs/download/LTEE-prophage-gd.log"
    conda: "../envs/download.yml"
    shell:
        """
        mkdir -p {REF_DIR} > {log} 2>&1
        wget -q -O {output} "{_LTEE_REF_URL}/prophage-amplifications.gd" >> {log} 2>&1
        """

rule download_LTEE_empty_mask_gd:
    output: REF_DIR + "/empty.mask.gd"
    log:
        "logs/download/LTEE-empty-mask-gd.log"
    conda: "../envs/download.yml"
    shell:
        """
        mkdir -p {REF_DIR} > {log} 2>&1
        wget -q -O {output} "{_LTEE_REF_URL}/empty.mask.gd" >> {log} 2>&1
        """

if _ANCESTOR_DOWNLOAD_FILE:
    rule download_LTEE_ancestor_gd:
        output: _ANCESTOR_DOWNLOAD_FILE
        log:
            "logs/download/LTEE-ancestor-gd.log"
        conda: "../envs/download.yml"
        shell:
            """
            wget -q -O {output} "{_ANCESTOR_DOWNLOAD_URL}" > {log} 2>&1
            """


# ── Per-sample curation pipeline ─────────────────────────────────────────────
# Note: the 00_header / 01_curate_add / 02_curate_remove inputs of curate_LTEE_gd
# are produced by the separate 'curate-LTEE-clones-init' workflow (see the guard at
# the top of this file), not by any rule here.

rule curate_LTEE_gd:
    input:
        initial    = "gd/{sample}.gd",
        subtracts  = "02_curate_remove/{sample}.gd",
        add        = "01_curate_add/{sample}.gd",
        header     = "00_header/{sample}.gd",
        reference  = REF_GBK
    output:
        "03_curated/{sample}.gd"
    log:
        "logs/curate/curate-LTEE-clones-{sample}.log"
    conda: BRESEQ_ENV
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
        "logs/curate/normalize-LTEE-clones-{sample}.log"
    conda: BRESEQ_ENV
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
        "logs/curate/mask-LTEE-clones-{sample}.log"
    conda: BRESEQ_ENV
    shell:
        """
        gdtools REMOVE   -c type==CON -o {output}.tmp1 {input.gd} > {log} 2>&1
        gdtools SUBTRACT -o {output}.tmp2 {output}.tmp1 {input.prophage} {input.ancestors} >> {log} 2>&1
        gdtools MASK     -v --mask-mode SMALL -o {output} {output}.tmp2 {input.mask} >> {log} 2>&1
        rm {output}.tmp1 {output}.tmp2
        """

rule no_is_adjacent_LTEE_gd:
    input:
        "05_normalized_masked_gd/{sample}.gd"
    output:
        "06_normalized_masked_no_IS_adjacent_gd/{sample}.gd"
    log:
        "logs/curate/no-is-adjacent-LTEE-clones-{sample}.log"
    conda: BRESEQ_ENV
    shell:
        """
        gdtools REMOVE -e -c TYPE!=UN -c adjacent!=UNDEFINED -o {output} {input} > {log} 2>&1
        """

rule apply_LTEE_mutations:
    input:
        gd        = "04_final_normalized_gd/{sample}.gd",
        reference = REF_GBK
    output:
        "mutants/{sample}.gff"
    log:
        "logs/curate/apply-LTEE-mutations-{sample}.log"
    conda: BRESEQ_ENV
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
        "output/compare_normalized.html"
    log:
        "logs/curate/compare-LTEE-normalized.log"
    conda: BRESEQ_ENV
    shell:
        """
        gdtools COMPARE -p -r {input.reference} -o {output} {input.gd_files} > {log} 2>&1
        """

rule compare_LTEE_normalized_masked:
    input:
        gd_files  = expand("05_normalized_masked_gd/{sample}.gd", sample=samples),
        reference = REF_GBK
    output:
        "output/compare_normalized_masked.html"
    log:
        "logs/curate/compare-LTEE-normalized-masked.log"
    conda: BRESEQ_ENV
    shell:
        """
        gdtools COMPARE -p -r {input.reference} -o {output} {input.gd_files} > {log} 2>&1
        """

rule mark_discrepancies_LTEE_compare:
    input:
        html          = "output/compare_normalized_masked.html",
        discrepancies = "07_phylogeny/discrepancies",   # directory() output of postprocess rule
    output:
        "output/compare_normalized_masked.discrepancies.html"
    run:
        import os, glob, re, html as _html

        # 1. Collect the set of discrepant keys from the discrepancy filenames.
        #    Current tree_utils.pl names files: tree.<i>.<mut_name>.<i>.tre
        #    We treat every "."-delimited field of each basename as a candidate key,
        #    so a data-key equal to the index <i> OR the mutation name <mut_name>
        #    will match. --- ADJUST THIS ONE BLOCK once the finalized breseq
        #    discrepancy-filename format is known. ---
        discrepant_keys = set()
        for path in glob.glob(os.path.join(input.discrepancies, "*")):
            base = os.path.basename(path)
            base_no_ext = os.path.splitext(base)[0]
            discrepant_keys.add(base_no_ext)                 # whole basename sans extension
            discrepant_keys.update(base_no_ext.split("."))   # each dotted field

        # 2. Rewrite each <tr ... data-key="K" ...> whose K is discrepant, injecting a
        #    light-red background. Uses a CSS <style> block keyed on data-key so we never
        #    have to merge into existing style/bgcolor attributes on the row.
        with open(input.html) as fh:
            doc = fh.read()

        present_keys = set(re.findall(r'data-key="([^"]*)"', doc))
        hit_keys = sorted(k for k in present_keys if k in discrepant_keys)

        if hit_keys:
            selectors = ", ".join(
                'tr[data-key="%s"]' % _html.escape(k, quote=True) for k in hit_keys
            )
            style_block = (
                "\n<style>%s { background-color: #ffb3b3 !important; }</style>\n"
                % selectors
            )
            # Inject before </head> if present, else at the very top of the document.
            if "</head>" in doc:
                doc = doc.replace("</head>", style_block + "</head>", 1)
            else:
                doc = style_block + doc

        with open(output[0], "w") as fh:
            fh.write(doc)


rule count_LTEE_initial:
    input:
        gd_files  = expand("gd/{sample}.gd", sample=samples),
        reference = REF_GBK
    output:
        "output/count.initial.csv"
    log:
        "logs/curate/count-LTEE-initial.log"
    conda: BRESEQ_ENV
    shell:
        """
        gdtools COUNT -o {output} -r {input.reference} {input.gd_files} > {log} 2>&1
        """

rule count_LTEE_final:
    input:
        gd_files  = expand("04_final_normalized_gd/{sample}.gd", sample=samples),
        reference = REF_GBK
    output:
        "output/count.final.csv"
    log:
        "logs/curate/count-LTEE-final.log"
    conda: BRESEQ_ENV
    shell:
        """
        gdtools COUNT -o {output} -r {input.reference} {input.gd_files} > {log} 2>&1
        """

rule count_LTEE_final_masked:
    input:
        gd_files  = expand("05_normalized_masked_gd/{sample}.gd", sample=samples),
        reference = REF_GBK
    output:
        "output/count.final_masked.csv"
    log:
        "logs/curate/count-LTEE-final-masked.log"
    conda: BRESEQ_ENV
    shell:
        """
        gdtools COUNT -o {output} -r {input.reference} {input.gd_files} > {log} 2>&1
        """

rule LTEE_phylogeny_tree:
    input:
        sample_gds   = expand("05_normalized_masked_gd/{sample}.gd", sample=samples),
        ancestor_gds = ancestor_files,
        reference    = REF_GBK
    output:
        tre          = "07_phylogeny/tree.tre",
        genotypes    = "07_phylogeny/tree.genotypes.txt",
        sample_key   = "07_phylogeny/tree.sample.key.txt",
        mutation_key = "07_phylogeny/tree.mutation.key.txt"
    log:
        "logs/curate/LTEE-phylogeny.log"
    conda: BRESEQ_ENV
    shell:
        """
        mkdir -p 07_phylogeny
        gdtools PHYLOGENY -p -a -o 07_phylogeny/tree -r {input.reference} {input.ancestor_gds} {input.sample_gds} > {log} 2>&1
        """

rule LTEE_phylogeny_postprocess:
    input:
        tre          = "07_phylogeny/tree.tre",
        genotypes    = "07_phylogeny/tree.genotypes.txt",
        reference    = REF_GBK
    output:
        rerooted     = "07_phylogeny/tree.rerooted.tre",
        rescaled     = "07_phylogeny/tree.rerooted.rescaled.tre",
        final        = "output/final.tre"
    log:
        "logs/curate/LTEE-phylogeny-postprocess.log"
    conda: "../envs/bioperl.yml"
    shell:
        """
        perl {TREE_UTILS} ROOT-ANCESTOR -i {input.tre} -o {output.rerooted} > {log} 2>&1
        perl {TREE_UTILS} SCALE-PHYLIP  -i {output.rerooted} -o {output.rescaled} -p {input.genotypes} >> {log} 2>&1
        mkdir -p 07_phylogeny/discrepancies
        perl {TREE_UTILS} DISCREPANCIES -i {output.rerooted} -p 07_phylogeny/tree -o 07_phylogeny/discrepancies/tree >> {log} 2>&1
        cp {output.rescaled} {output.final} >> {log} 2>&1
        """
