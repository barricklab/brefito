#sample info not used but this loads the BRESEQ_ENV variable as well
# This workflow derives its samples from the gd/ directory, so it does not
# require sample information from data.csv / the read/assembly directories.
SAMPLE_INFO_REQUIRED = False
try: sample_info
except NameError:
    include: "load-sample-info.smk"

import glob
import os

brefito_config = {key.upper(): value for key, value in config.items()}

# Detect samples from gd directory
_initial_gd_dir = "gd"
samples = sorted([os.path.splitext(os.path.basename(f))[0]
                  for f in glob.glob(_initial_gd_dir + "/*.gd")])

# If specific sample names were given on the brefito command line, restrict the run
# to just those samples.
if 'SAMPLES' in brefito_config:
    _requested_samples = brefito_config['SAMPLES'].split("_,_")
    samples = [s for s in samples if s in _requested_samples]

# Build target list: the three curation stub files for every sample. These are then
# hand-edited by the user and consumed (but never regenerated) by curate-LTEE-clones.
_all_targets = []
if samples:
    _all_targets += expand("00_header/{sample}.gd", sample=samples)
    _all_targets += expand("01_curate_add/{sample}.gd", sample=samples)
    _all_targets += expand("02_curate_remove/{sample}.gd", sample=samples)


rule all_curate_ltee_clones_init:
    input: _all_targets
    default_target: True


# ── Stub creation (ancient() ensures existing files are never overwritten) ───

rule create_LTEE_header:
    input:
        gd = ancient("gd/{sample}.gd")
    output:
        "00_header/{sample}.gd"
    log:
        "logs/curate/create-LTEE-header-{sample}.log"
    conda: BRESEQ_ENV
    shell:
        """
        gdtools SUBTRACT -o {output}.tmp {input.gd} {input.gd} > {log} 2>&1
        gdtools NOT-EVIDENCE -o {output} {output}.tmp >> {log} 2>&1
        rm {output}.tmp
        """

rule create_LTEE_curate_add:
    input:
        gd = ancient("gd/{sample}.gd")
    output:
        "01_curate_add/{sample}.gd"
    log:
        "logs/curate/create-LTEE-curate-add-{sample}.log"
    conda: BRESEQ_ENV
    shell:
        """
        gdtools SUBTRACT -o {output}.tmp {input.gd} {input.gd} > {log} 2>&1
        gdtools NOT-EVIDENCE -o {output} {output}.tmp >> {log} 2>&1
        rm {output}.tmp
        """

rule create_LTEE_curate_remove:
    input:
        gd = ancient("gd/{sample}.gd")
    output:
        "02_curate_remove/{sample}.gd"
    log:
        "logs/curate/create-LTEE-curate-remove-{sample}.log"
    conda: BRESEQ_ENV
    shell:
        """
        gdtools SUBTRACT -o {output}.tmp {input.gd} {input.gd} > {log} 2>&1
        gdtools NOT-EVIDENCE -o {output} {output}.tmp >> {log} 2>&1
        rm {output}.tmp
        """
