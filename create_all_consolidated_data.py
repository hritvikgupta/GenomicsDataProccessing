# First script: create_all_consolidated_data.py
import os
import gzip
import json
import glob
import argparse

# Column mapping definitions
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

def process_gz_file(file_path, software_type):
    """Process a gzipped file based on its software type."""
    snp_data = {}
    mapping = COLUMN_MAPPINGS[software_type]
    
    filename = os.path.basename(file_path)
    phenotype = filename.split('.')[0]
    cohort = filename.split('.')[1]
    
    print(f"Processing {phenotype} - {cohort} ({software_type})")
    
    with gzip.open(file_path, 'rt') as f:
        header = next(f)
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) < max(mapping.values()) + 1:
                continue
            
            try:
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
                
                snp_id = fields[mapping['id']]
                
                if snp_id not in snp_data:
                    snp_data[snp_id] = {}
                
                if phenotype not in snp_data[snp_id]:
                    snp_data[snp_id][phenotype] = {}
                
                snp_data[snp_id][phenotype][cohort] = snp_info
                
            except (IndexError, ValueError) as e:
                print(f"Error processing line in {file_path}: {str(e)}")
                continue
    
    return snp_data

def process_software_type(input_dir, output_dir, software_type):
    consolidated_data = {}
    pattern = os.path.join(input_dir, f"*.{software_type}.sumstats.txt.gz")
    gz_files = glob.glob(pattern)
    
    if not gz_files:
        print(f"No files found for {software_type}")
        return
        
    print(f"\nProcessing {software_type} files...")
    
    for file_path in gz_files:
        try:
            file_data = process_gz_file(file_path, software_type)
            consolidated_data.update(file_data)
        except Exception as e:
            print(f"Error processing {file_path}: {str(e)}")

    output_file = os.path.join(output_dir, f'consolidated_snp_data_{software_type}.json.gz')
    with gzip.open(output_file, 'wt') as f:
        json.dump(consolidated_data, f)

    print(f"{software_type} data saved to {output_file}")
    print(f"Total SNPs processed for {software_type}: {len(consolidated_data)}")

def main():
    parser = argparse.ArgumentParser(description='Process GWAS summary statistics files for both MRMEGA and GWAMA.')
    parser.add_argument('-i', '--input_dir', required=True, 
                        help='Input directory containing the summary statistics files')
    parser.add_argument('-o', '--output_dir', required=True, 
                        help='Output directory for processed files')
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    
    # Process both software types
    for software_type in ['mrmega', 'gwama']:
        process_software_type(args.input_dir, args.output_dir, software_type)

if __name__ == "__main__":
    main()