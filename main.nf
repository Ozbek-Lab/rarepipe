#!/usr/bin/env nextflow
import java.text.SimpleDateFormat
nextflow.enable.dsl=2

/*
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
*/
"""
touch "NONE.fq"
"""


params.input_merged_vcf = "NONE.fq"  // Path to pre-merged cohort VCF (if available)
params.input_merged_vcf_index = "NONE.fq.csi"  // Path to index of pre-merged cohort VCF (if available)
params.skip_merge = true  // If true, will not merge the samples into a single VCF
// --- Log pipeline parameters ---
log.info """
          E X O M I S E R - N F   P I P E L I N E (Cohort Run with Pedigree)
          ==================================================================
          samplesheet                    : ${params.samplesheet}
          output dir                     : ${params.outdir}
          template YAML                  : ${params.template_yaml}
          exomiser YAML python path      : ${params.exomiser_yaml_python_script_path}
          Exomiser JAR                   : ${params.exomiser_jar_path}
          Exomiser data dir              : ${params.exomiser_data_dir_host}
          ClinVar reference file         : ${params.clinvar_reference_file}
          merged VCF                     : ${params.input_merged_vcf}
          merged VCF index               : ${params.input_merged_vcf_index}

          ---
          """
          .stripIndent()

workflow {
    main:
        channel
            .fromPath(params.samplesheet)
            //.ifEmpty { error "Samplesheet file not found: ${params.samplesheet}" }
            .splitCsv(header:true, sep:'\t')
            .map { row ->
                // VALIDATION: If sample is affected (2), it MUST have HPO terms.
                if (row.affected == '2' && !row.hpo?.trim()) {
                    error "Validation Error: Sample '${row.sample_id}' is affected but has no HPO terms provided in the samplesheet."
                }
                // Create a metadata map for each sample
                [
                    famid: row.famid,
                    sample_id: row.sample_id,
                    father_id: row.father_id,
                    mother_id: row.mother_id,
                    sex: row.sex,
                    affected: row.affected,
                    vcf: file(row.vcf, checkIfExists:true),
                    hpo: row.hpo
                ]
            }
            .collect()
            .set { cohort_ch }

        if(params.input_merged_vcf != "NONE.fq"){ // If a provided merged VCF is provided, filter the samples that exist in the provided merged VCF

            log.info "Filtering samples that exist in the provided merged VCF: ${params.input_merged_vcf} from the samplesheet samples list: ${params.samplesheet}."
            FILTER_EXISTING_SAMPLES(
                cohort_ch,
                file(params.input_merged_vcf, checkIfExists:true),
                file(params.input_merged_vcf_index, checkIfExists:true)
            )

            samples_to_merge_list_ch = FILTER_EXISTING_SAMPLES.out.samples_to_merge_list
                .splitCsv(header:true, sep:'\t')
                .map { row ->
                    // Create a metadata map for each sample
                    [
                        sample_id: row.sample_id,
                        vcf: row.vcf
                    ]
                }
                .collect()


            if(params.skip_merge == false){
                MERGE_COHORT_VCF(
                    samples_to_merge_list_ch,
                    file(params.input_merged_vcf),
                    file(params.input_merged_vcf_index)
                )
                merged_vcf_cohort = MERGE_COHORT_VCF.out.merged_vcf_cohort_path
                merged_vcf_cohort_index = MERGE_COHORT_VCF.out.merged_vcf_cohort_index_path
            }
            else{
                SKIP_MERGED_VCF(
                    file(params.input_merged_vcf),
                    file(params.input_merged_vcf_index)
                )
                merged_vcf_cohort = SKIP_MERGED_VCF.out.merged_vcf_path
                merged_vcf_cohort_index = SKIP_MERGED_VCF.out.merged_vcf_index_path
            }
        }
        else{ // If no provided merged VCF, create a new merged VCF from individual sample VCFs
            log.info "No provided merged VCF, will create new merged VCF from individual sample VCFs."
            // 2.2 Merge VCFs for the entire cohort (affected + unaffected)
            CREATE_MERGED_VCF(cohort_ch)
            // Extract file paths from tuple (meta_list, vcf, vcf_index) - process emits once, so channel broadcasts to all subscribers
            merged_vcf_cohort = CREATE_MERGED_VCF.out.merged_vcf_cohort_path
            merged_vcf_cohort_index = CREATE_MERGED_VCF.out.merged_vcf_cohort_index_path
        }

        // 3. Split the merged VCF and annotate each sample with parental genotypes
        affected_samples_ch = cohort_ch.flatMap { all_meta ->
            all_meta.collect { sample_meta ->
                tuple(sample_meta, all_meta)
            }
        }
        // 4. FILTER for AFFECTED samples before Annotation
        .filter { sample_meta, _all_meta ->
            sample_meta.affected == '2'
        }
        //affected_samples_ch.view { all_meta -> "AFFECTED SAMPLE: ${all_meta.sample_id}" }

        EXTRACT_AND_ANNOTATE_SAMPLE(
            affected_samples_ch,
            merged_vcf_cohort,
            merged_vcf_cohort_index
        )

        // 5. Run Exomiser and subsequent steps ONLY on the filtered affected samples
        CREATE_YAML(
            EXTRACT_AND_ANNOTATE_SAMPLE.out,
            file(params.template_yaml, checkIfExists:true),
            file(params.exomiser_yaml_python_script_path, checkIfExists:true)
        )
        
        RUN_EXOMISER(CREATE_YAML.out.yaml_and_vcf)

        // The explicit call below is REMOVED.
        GET_CLINVAR_REF(
            file(params.clinvar_reference_file).getParent()
        )

        // This call implicitly runs GET_CLINVAR_REF once and combines its
        // single output with every item from RUN_EXOMISER.out
        RUN_VEP_ANNOTATION(
            GET_CLINVAR_REF.out,
            RUN_EXOMISER.out.exomiser_results,
            params.vep_data_dir
        )

        RUN_VEP_SPLIT(
            RUN_VEP_ANNOTATION.out.vep_vcf
        )

        RUN_EXOMISER_SPLIT(
            RUN_VEP_SPLIT.out.vep_tsv,
            file(params.exomiser_split_python_script_path, checkIfExists:true)
        )
}

// --- Process Definitions ---
// Processes MERGE_COHORT_VCF and EXTRACT_AND_ANNOTATE_SAMPLE remain unchanged
// as they need to run on the full cohort. All other processes are also identical
// but will now only be executed for the filtered, affected samples.

process FILTER_EXISTING_SAMPLES {
    // Exists upon encountering any sample already present in the provided merged VCF
    input:
    val(meta_list)
    path(input_merged_vcf)
    path(input_merged_vcf_index)

    output:
    path("samples_to_merge_list.txt"), emit: samples_to_merge_list

    script:
    def sample_ids = meta_list.collect { it.sample_id }.join('" "')
    def vcf_files = meta_list.collect { it.vcf }.join('" "')
    def meta_list_size = meta_list.size()
    """
    # Gather existing samples from the provided merged VCF (one per line, trimmed)
    bcftools query -l ${input_merged_vcf} > existing_samples.txt

    # Build list of sample IDs from the samplesheet (per-sample VCFs)
    # that are NOT already present in the provided merged VCF
    SAMPLE_IDS=( "${sample_ids}" )
    VCF_FILES=( "${vcf_files}" )
    echo "sample_id\tvcf\n" > samples_to_merge_list.txt
    for sample_index in {0..${meta_list_size-1}}; do
        if ! grep -qFx "\${SAMPLE_IDS[sample_index]}" existing_samples.txt 2>/dev/null; then
            echo "\${SAMPLE_IDS[sample_index]}\t\${VCF_FILES[sample_index]}\n" >> samples_to_merge_list.txt
        fi
    done
    """
}

process SKIP_MERGED_VCF {
    tag "Skipping merging of samples into merged VCF and using provided merged VCF"
    publishDir "${params.outdir}/pipeline_info/merged_vcfs", mode: 'copy', pattern: "*.vcf.gz*", overwrite: true
    input:
    path(input_merged_vcf)
    path(input_merged_vcf_index)

    output:
    path(input_merged_vcf), emit: merged_vcf_path
    path(input_merged_vcf_index), emit: merged_vcf_index_path
    script:
    """
    touch ${input_merged_vcf}
    touch ${input_merged_vcf_index}
    """
}
process MERGE_COHORT_VCF {
    tag "Merging ${samples_to_merge_list.size()} samples into merged VCF"
    publishDir "${params.outdir}/pipeline_info/merged_vcfs", mode: 'copy', pattern: "*.vcf.gz*", overwrite: true

    input:
    val(samples_to_merge_list)
    path(input_merged_vcf)
    path(input_merged_vcf_index)

    output:
    path("cohort_${timestamp}.merged.vcf.gz"), emit: merged_vcf_cohort_path
    path("cohort_${timestamp}.merged.vcf.gz.tbi"), emit: merged_vcf_cohort_index_path

    script:
    def vcf_files = samples_to_merge_list.collect { it.vcf }.join(' ')
    // if(samples_to_merge_list.size() > 0){ vcf_files = samples_to_merge_list.collect { it.vcf }.join(' '); }

    def merged_file = input_merged_vcf.getName() == "NONE.fq" ? "" : input_merged_vcf.name

    Date now = new Date();
    SimpleDateFormat timestamp_formatter = new SimpleDateFormat("yyyy-MM-dd_HH-mm-ss");
    timestamp = timestamp_formatter.format(now)

    log.info "New sample IDs found, will create new merged VCF from individual sample VCFs."
    """
        echo "Processing with provided merged VCF: ${input_merged_vcf}"
        
        echo "Creating new merged VCF from individual sample VCFs"
        for vcf in ${vcf_files}; do bcftools index -f \$vcf; done
        bcftools merge -m none --threads 16 -Oz -o cohort.merged.preprocessed.vcf.gz ${merged_file} ${vcf_files}


    bcftools view cohort.merged.preprocessed.vcf.gz | \\
        bcftools norm -m-any | \\
        bcftools +fill-tags -- -t VAF | \\
        bcftools +setGT -- -t q -n . -i 'FORMAT/DP <= 3' 2>/dev/null | \\
        bcftools +setGT -- -t q -n . -i 'FORMAT/DP < 6 & FORMAT/VAF < 1' 2>/dev/null | \\
        bcftools +setGT -- -t q -n . -i '6 <= FORMAT/DP & FORMAT/DP < 24 & FORMAT/VAF < 0.35' 2>/dev/null | \\
        bcftools +setGT -- -t q -n . -i '24 <= FORMAT/DP & FORMAT/VAF< 0.17' 2>/dev/null | \\
        bcftools +fill-tags -- -t AC,AC_Hom,AC_Het,AC_Hemi,AF,AN,ExcHet,HWE | \\
        bcftools view -Oz -o cohort.merged.AC.vcf.gz && \\
    bcftools view cohort.merged.AC.vcf.gz -e "INFO/AC>6 || INFO/AN=0" -Oz -o cohort_${timestamp}.merged.vcf.gz

    bcftools index -t cohort_${timestamp}.merged.vcf.gz
    """
}

process CREATE_MERGED_VCF {
    tag "Creating new merged VCF with ${meta_list.size()} samples."
    publishDir "${params.outdir}/pipeline_info/merged_vcfs", mode: 'copy', pattern: "*.vcf.gz*", overwrite: true

    input:
    val(meta_list)

    output:
    path("cohort_${timestamp}.merged.vcf.gz"), emit: merged_vcf_cohort_path
    path("cohort_${timestamp}.merged.vcf.gz.tbi"), emit: merged_vcf_cohort_index_path

    script:
    def vcf_files = meta_list.collect { it.vcf }.join(' ')
    Date now = new Date();
    SimpleDateFormat timestamp_formatter = new SimpleDateFormat("yyyy-MM-dd_HH-mm-ss");
    timestamp = timestamp_formatter.format(now)

    """    
    echo "Creating new merged VCF from individual sample VCFs"

    for vcf in ${vcf_files}; do bcftools index -f \$vcf; done
    bcftools merge -m none --threads 16 -Oz -o cohort.merged.preprocessed.vcf.gz ${vcf_files}
    bcftools view cohort.merged.preprocessed.vcf.gz | \\
        bcftools norm -m-any | \\
        bcftools +fill-tags -- -t VAF | \\
        bcftools +setGT -- -t q -n . -i 'FORMAT/DP <= 3' 2>/dev/null | \\
        bcftools +setGT -- -t q -n . -i 'FORMAT/DP < 6 & FORMAT/VAF < 1' 2>/dev/null | \\
        bcftools +setGT -- -t q -n . -i '6 <= FORMAT/DP & FORMAT/DP < 24 & FORMAT/VAF < 0.35' 2>/dev/null | \\
        bcftools +setGT -- -t q -n . -i '24 <= FORMAT/DP & FORMAT/VAF< 0.17' 2>/dev/null | \\
        bcftools +fill-tags -- -t AC,AC_Hom,AC_Het,AC_Hemi,AF,AN,ExcHet,HWE | \\
        bcftools view -Oz -o cohort.merged.AC.vcf.gz && \\
    bcftools view cohort.merged.AC.vcf.gz -e "INFO/AC>6 || INFO/AN=0" -Oz -o cohort_${timestamp}.merged.vcf.gz

    bcftools index -t cohort_${timestamp}.merged.vcf.gz
    """
}
process EXTRACT_AND_ANNOTATE_SAMPLE {
    tag "${meta.sample_id}"
    input:
    tuple val(meta), val(all_meta)
    path(merged_vcf)
    path(merged_vcf_index)

    output:
    tuple val(meta), path("${meta.sample_id}.final.vcf.gz"), path("${meta.sample_id}.final.vcf.gz.tbi")

    // MODIFIED: Script now uses original parent VCFs for annotation
    script:
    def sample_id = meta.sample_id
    def father_id = meta.father_id
    def mother_id = meta.mother_id

    // Find the original VCF files for the parents from the full cohort metadata list
    def father_meta = all_meta.find { it.sample_id == father_id }
    def mother_meta = all_meta.find { it.sample_id == mother_id }
    def father_vcf_path = father_meta ? father_meta.vcf : ''
    def mother_vcf_path = mother_meta ? mother_meta.vcf : ''

    """
    echo '##FORMAT=<ID=maternal_GT,Number=1,Type=String,Description="Maternal genotype">' > mat_hdr.txt
    echo '##FORMAT=<ID=paternal_GT,Number=1,Type=String,Description="Paternal genotype">' > pat_hdr.txt
    # Extract this sample's variants from the merged cohort VCF
    bcftools view -I -s $sample_id -Oz -o current.vcf.gz $merged_vcf
    bcftools view -h current.vcf.gz > oldheader
    echo "\$(sed \\\$d oldheader; cat mat_hdr.txt; cat pat_hdr.txt; sed -n \\\$p oldheader)"> newheader
    bcftools reheader -h newheader current.vcf.gz | bcftools view -Oz -o reheadered.vcf.gz
    mv reheadered.vcf.gz current.vcf.gz

    # Annotate with paternal genotype from the ORIGINAL father's VCF
    if [ "$father_id" != "0" ]; then
        bcftools index -f "${father_vcf_path}"
        bcftools query -f '%CHROM\\t%POS\\t[%GT]\\n' -s "$father_id" "${father_vcf_path}" | bgzip -c > paternal_gt.txt.gz
        tabix -s1 -b2 -e2 paternal_gt.txt.gz
        bcftools annotate -a paternal_gt.txt.gz -h pat_hdr.txt -c CHROM,POS,FORMAT/paternal_GT current.vcf.gz -Oz -o annotated.tmp.vcf.gz
        mv annotated.tmp.vcf.gz current.vcf.gz
    fi

    # Annotate with maternal genotype from the ORIGINAL mother's VCF
    if [ "$mother_id" != "0" ]; then
        bcftools index -f "${mother_vcf_path}"
        bcftools query -f '%CHROM\\t%POS\\t[%GT]\\n' -s "$mother_id" "${mother_vcf_path}" | bgzip -c > maternal_gt.txt.gz
        tabix -s1 -b2 -e2 maternal_gt.txt.gz
        bcftools annotate -a maternal_gt.txt.gz -h mat_hdr.txt -c CHROM,POS,FORMAT/maternal_GT current.vcf.gz -Oz -o annotated.tmp.vcf.gz
        mv annotated.tmp.vcf.gz current.vcf.gz
    fi

    mv current.vcf.gz ${sample_id}.final.vcf.gz
    bcftools index -t ${sample_id}.final.vcf.gz
    """
}

process CREATE_YAML {
    tag "${meta.sample_id}"
    publishDir "${params.outdir}/pipeline_info/generated_yamls", mode: 'copy', pattern: "*.exomiser.yaml", overwrite: true

    input:
    tuple val(meta), path(annotated_vcf), path(annotated_vcf_index)
    path template_yaml_file
    path python_script

    output:
    tuple val(meta), path("*.exomiser.yaml"), path(annotated_vcf), emit: yaml_and_vcf

    script:
    def sample_id = meta.sample_id
    def hpo_terms = meta.hpo
    """
    python3 "${python_script}" \\
        -t "${template_yaml_file}" \\
        -v "${annotated_vcf.name}" \\
        -H "${hpo_terms}" \\
        -o "${sample_id}.exomiser.yaml" \\
        -s "${sample_id}"
    """
}

process RUN_EXOMISER {
    tag "${meta.sample_id}"
    publishDir path: "${params.outdir}/${meta.sample_id}", source: "results/${meta.sample_id}", mode: 'copy', overwrite: true, failIfExists: false

    input:
    tuple val(meta), path(generated_yaml), path(original_vcf_staged)

    output:
    tuple val(meta), path("results/${meta.sample_id}"), path("results/${meta.sample_id}/${meta.sample_id}.vcf.gz"), emit: exomiser_results

    script:
    def sample_id = meta.sample_id
    """
    java -jar "${params.exomiser_jar_path}" \\
        --analysis "${generated_yaml.name}" \\
        --exomiser.data-directory="${params.exomiser_data_dir_host}/" \\
        --exomiser.hg38.data-version="${params.exomiser_data_version}" \\
        --exomiser.phenotype.data-version="${params.exomiser_phenotype_version}" \\
        --exomiser.hg38.remm-path="${params.exomiser_remm_path_host}" \\
        --exomiser.hg38.cadd-snv-path="${params.exomiser_cadd_snv_path_host}" \\
        --exomiser.hg38.cadd-in-del-path="${params.exomiser_cadd_indel_path_host}"
    """
}

process GET_CLINVAR_REF{
    publishDir "${clinvar_ref_dir}", mode: 'copy', overwrite: true
    input:
    path clinvar_ref_dir

    output:
    tuple path("clinvar.vcf.gz"), path("clinvar.vcf.gz.tbi"), emit: clinvar_vcf
    
    script:
    """
    REMOTE_URL="https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz"
    LOCAL_FILE="${params.clinvar_reference_file}"
    
    if [ ! -f "\$LOCAL_FILE" ]; then
        curl -L -o clinvar.vcf.gz "\$REMOTE_URL"
        curl -L -o clinvar.vcf.gz.tbi "\${REMOTE_URL}.tbi"
    else
        ln -s \$LOCAL_FILE clinvar.vcf.gz
        ln -s \${LOCAL_FILE}.tbi clinvar.vcf.gz.tbi
    fi
    """
}

process RUN_VEP_ANNOTATION{
    tag "${meta.sample_id}"
    container 'ensemblorg/ensembl-vep:release_113.4'

    input:
    tuple path(clinvar_vcf), path(clinvar_tbi)
    tuple val(meta), path(exomiser_dir), path(exomiser_vcf)
    path vep_data_dir

    output:
    tuple val(meta), path("${meta.sample_id}.vep.vcf.gz"), emit: vep_vcf

    script:
    def sample_id = meta.sample_id
    """
    vep \\
        --fork 16 --buffer_size 50000 -i $exomiser_vcf \\
        --dir_cache $vep_data_dir --offline --cache --merged --force_overwrite \\
        --dir_plugins $vep_data_dir/plugins --format vcf \\
        --af_1kg --af_gnomade --af_gnomadg --check_existing \\
        --per_gene --pick_order mane_select \\
        --compress_output gzip --vcf --output_file ${sample_id}.vep.vcf.gz \\
        --custom file=$clinvar_vcf,short_name=ClinVar,format=vcf,type=exact,coords=0,fields=CLNSIG%CLNREVSTAT%CLNDN \\
        --custom file=$vep_data_dir/ucsc/repeatmasker.GRCh38.bed.gz,short_name=repeatmasker,format=bed,type=overlap,coords=1 \\
        --custom file=$vep_data_dir/promoterai/promoterAI_tss500.GRCh38.vcf.gz,short_name=PromoterAI,format=vcf,type=exact,coords=0,fields=PromoterAIScore \\
        --custom file=$vep_data_dir/primateai3d/PrimateAI-3D.hg38.vcf.gz,short_name=PrimateAI3D,format=vcf,type=exact,coords=0,fields=PrimateAI3DScore%PrimateAI3DPrediction \\
        --plugin SpliceAI,snv=$vep_data_dir/spliceai/spliceai_scores.raw.snv.hg38.vcf.gz,indel=$vep_data_dir/spliceai/spliceai_scores.raw.indel.hg38.vcf.gz \\
        --plugin dbNSFP,$vep_data_dir/dbNSFP/dbNSFP5.2a_grch38.gz,AlphaMissense_pred,AlphaMissense_rankscore,AlphaMissense_score,CADD_raw,CADD_phred,MetaRNN_pred,MetaRNN_score,RegeneronMe_ALL_AC,MutationTaster_pred,MutationTaster_score,ClinPred_pred
    """
}

process RUN_VEP_SPLIT{
    tag "${meta.sample_id}"
    publishDir "${params.outdir}/${meta.sample_id}", mode: 'copy', pattern: "*.vep.tsv", overwrite: true

    input:
    tuple val(meta), path(vep_vcf)

    output:
    tuple val(meta), path("${meta.sample_id}.vep.tsv"), emit: vep_tsv

    script:
    def sample_id = meta.sample_id
    """
    columns="[%SAMPLE]\\t%CHROM\\t%POS\\t%REF\\t%ALT\\t%ID\\t%FILTER\\t%Exomiser[\\t%GT\\t%maternal_GT\\t%paternal_GT\\t%VAF\\t%AD\\t%DP]\\t%AC_Het\\t%AC_Hom\\t%AC_Hemi\\t%HWE\\t%ExcHet\\t%AN\\t%AC\$(bcftools +split-vep -l ${vep_vcf} | cut -f2 | sed 's/^/\\\\\\\\t%/' | tr -d '\n' | xargs)"
    header=\$(echo "\$columns" | sed "s/%//g;s/\\[//g;s/\\]//g")

    (echo -e \$header; bcftools view ${vep_vcf} | \\
        bcftools +split-vep -f "\$columns\\n" -d -A tab) > ${sample_id}.vep.tsv
    """
}

process RUN_EXOMISER_SPLIT{
    tag "${meta.sample_id}"
    publishDir "${params.outdir}/${meta.sample_id}", mode: 'copy', pattern: "*.vep.populated.tsv", overwrite: true

    input:
    tuple val(meta), path(vep_tsv)
    path python_script

    output:
    path "${meta.sample_id}.vep.populated.tsv", emit: populated_tsv

    script:
    def sample_id = meta.sample_id
    """
    python3 $python_script $vep_tsv "${sample_id}.vep.populated.tsv"
    """
}