========================================================================================
    nf-exomiser-pipeline (Cohort Merging & Local Exomiser Execution)
========================================================================================
    Pipeline Description:
    1. Reads a samplesheet, validating that affected samples have HPO terms.
    2. Merges all VCFs into a single cohort-wide VCF.
    3. Filters the merged VCF, annotates de novo mutations, and splits it back into
       single-sample VCFs, annotating parental genotypes.
    4. Filters for AFFECTED samples only.
    5. The final annotated VCFs for affected samples are processed with Exomiser and VEP.
----------------------------------------------------------------------------------------
