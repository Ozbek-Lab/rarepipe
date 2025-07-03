#!/usr/bin/env python3

import argparse
import os
import yaml

def create_exomiser_yaml(template_file, vcf_file_path, hpo_ids, output_yaml_file, sample_id_for_output):
    """
    Creates an Exomiser analysis YAML file from a template.
    The VCF path written to the YAML is just the filename (relative path).
    The outputDirectory and outputFileName are based on sample_id_for_output.
    """
    try:
        with open(template_file, 'r') as f:
            template_data = yaml.safe_load(f)
    except FileNotFoundError:
        print(f"Error: Template file not found at {template_file}")
        raise # Propagate error to Nextflow
    except yaml.YAMLError as e:
        print(f"Error parsing template YAML file {template_file}: {e}")
        raise # Propagate error to Nextflow

    # Get just the VCF filename
    vcf_filename = os.path.basename(vcf_file_path)
    
    # Update template data
    if 'analysis' not in template_data:
        template_data['analysis'] = {}
    template_data['analysis']['vcf'] = vcf_filename # Relative path (just the filename)
    template_data['analysis']['hpoIds'] = hpo_ids.split() # Split space-separated string

    if 'outputOptions' not in template_data:
        template_data['outputOptions'] = {}
    # Use the provided sample_id_for_output for naming
    template_data['outputOptions']['outputFileName'] = sample_id_for_output
    # This directory will be created by Exomiser relative to its working directory inside Docker.
    template_data['outputOptions']['outputDirectory'] = os.path.join('results', sample_id_for_output)

    try:
        with open(output_yaml_file, 'w') as f:
            yaml.dump(template_data, f, sort_keys=False, indent=2)
    except IOError as e:
        print(f"Error writing output YAML file {output_yaml_file}: {e}")
        raise # Propagate error to Nextflow
    
    # No return value needed for this version, and no special print to stdout

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Create Exomiser analysis YAML file with relative VCF path.')
    parser.add_argument('-t', '--template', help='Path to template YAML file', required=True)
    parser.add_argument('-v', '--vcf', help='Path to VCF file (Nextflow will stage this)', required=True)
    parser.add_argument('-o', '--output', help='Path to output YAML file (e.g., sample_id.exomiser.yaml)', required=True)
    parser.add_argument('-H', '--hpo_ids', help='HPO IDs separated by spaces (e.g., "HP:0001234 HP:0005678").', required=True)
    parser.add_argument('-s', '--sample_id', help='Sample ID to use for output naming in the YAML.', required=True)
    
    args = parser.parse_args()

    try:
        create_exomiser_yaml(args.template, args.vcf, args.hpo_ids, args.output, args.sample_id)
        # Script will exit with 0 on success, non-zero on error if exceptions are raised
    except Exception as e:
        print(f"An error occurred in create_exomiser_yaml: {e}")
        exit(1) # Ensure Nextflow sees a failure

