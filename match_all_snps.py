# Second script: match_all_snps.py
import os
import gzip
import json
import glob
import logging
import argparse

# Column mappings (same as first script)
COLUMN_MAPPINGS = {
    'mrmega': {
        'id': 0, 'chr': 1, 'pos': 2, 'ref': 3, 'alt': 4,
        'beta': 5, 'se': 6, 'pval': 7, 'aaf': 12,
        'n': 16, 'n_study': 18
    },
    'gwama': {
        'id': 0, 'chr': 1, 'pos': 2, 'ref': 3, 'alt': 4,
        'beta': 5, 'se': 6, 'pval': 7, 'aaf': 12,
        'n': 16, 'n_study': 18
    }
}

def setup_logging(output_dir, software):
    log_file = os.path.join(output_dir, f'snp_matching_{software}.log')
    logging.basicConfig(
        filename=log_file,
        level=logging.INFO,
        format='%(asctime)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S',
        force=True  # Reset logging for each software type
    )

def process_gz_file(file_path, original_snp_data, software_type):
    snp_data = {}
    matched_snps = set()
    
    filename = os.path.basename(file_path)
    phenotype = filename.split('.')[0]
    cohort = filename.split('.')[1]
    
    mapping = COLUMN_MAPPINGS[software_type]
    
    with gzip.open(file_path, 'rt') as f:
        next(f)
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) < max(mapping.values()) + 1:
                continue
            
            try:
                snp_id = fields[mapping['id']]
                if snp_id not in original_snp_data:
                    continue
                
                matched_snps.add(snp_id)
                
                snp_info = {
                    'chromosome': fields[mapping['chr']],
                    'position': fields[mapping['pos']],
                    'ref_allele': fields[mapping['ref']],
                    'alt_allele': fields[mapping['alt']],
                    'beta': fields[mapping['beta']],
                    'se': fields[mapping['se']],
                    'p_value': fields[mapping['pval']],
                    'aaf': fields[mapping['aaf']],
                    'n': fields[mapping['n']],
                    'n_study': fields[mapping['n_study']]
                }
                
                if snp_id not in snp_data:
                    snp_data[snp_id] = {}
                
                if phenotype not in snp_data[snp_id]:
                    snp_data[snp_id][phenotype] = {}
                
                snp_data[snp_id][phenotype][cohort] = snp_info
                
            except (IndexError, ValueError) as e:
                print(f"Error processing line in {file_path}: {str(e)}")
                continue
    
    return snp_data, matched_snps

def process_software_type(input_dir, consolidated_dir, output_dir, software_type):
    # Setup logging for this software type
    setup_logging(output_dir, software_type)
    logging.info(f"Starting processing for {software_type}")

    # Load the original consolidated SNP data
    consolidated_file = os.path.join(consolidated_dir, f'consolidated_snp_data_{software_type}.json.gz')
    with gzip.open(consolidated_file, 'rt') as f:
        original_snp_data = json.load(f)
    
    new_consolidated_data = {}
    all_matched_snps = set()

    # Find relevant tabix files
    pattern = os.path.join(input_dir, f"*.{software_type}_pval_*.gz")
    gz_files = glob.glob(pattern)
    gz_files = [f for f in gz_files if not f.endswith('.tbi')]

    if not gz_files:
        print(f"No tabix files found for {software_type}")
        return

    print(f"\nProcessing {software_type} files...")
    total_files = len(gz_files)
    for index, file_path in enumerate(gz_files, 1):
        print(f"Processing file {index}/{total_files}: {os.path.basename(file_path)}")
        try:
            file_data, matched_snps = process_gz_file(file_path, original_snp_data, software_type)
            all_matched_snps.update(matched_snps)
            
            for snp_id, snp_info in file_data.items():
                if snp_id not in new_consolidated_data:
                    new_consolidated_data[snp_id] = {}
                new_consolidated_data[snp_id].update(snp_info)
            
            logging.info(f"File: {os.path.basename(file_path)}")
            logging.info(f"Matched SNPs in this file: {len(matched_snps)}")
            logging.info("--------------------")
        except Exception as e:
            error_msg = f"Error processing {file_path}: {str(e)}"
            print(error_msg)
            logging.error(error_msg)

    # Save matched data
    output_file = os.path.join(output_dir, f'matched_snp_data_{software_type}.json.gz')
    with gzip.open(output_file, 'wt') as f:
        json.dump(new_consolidated_data, f)

    print(f"{software_type} data saved to {output_file}")
    print(f"Total matched SNPs for {software_type}: {len(all_matched_snps)}")
    
    logging.info(f"Total matched SNPs for {software_type}: {len(all_matched_snps)}")
    logging.info("Processing completed")

def main():
    parser = argparse.ArgumentParser(description='Match SNPs from tabix files for both MRMEGA and GWAMA.')
    parser.add_argument('-i', '--input_dir', required=True,
                        help='Input directory containing the tabix files')
    parser.add_argument('-c', '--consolidated_dir', required=True,
                        help='Directory containing the consolidated SNP data files')
    parser.add_argument('-o', '--output_dir', required=True,
                        help='Output directory for matched data')
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    
    # Process both software types
    for software_type in ['mrmega', 'gwama']:
        process_software_type(args.input_dir, args.consolidated_dir, args.output_dir, software_type)

if __name__ == "__main__":
    main()