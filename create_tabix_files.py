import subprocess
import os
import base64
import glob
import argparse

# Column mapping definitions
COLUMN_MAPPINGS = {
    'mrmega': {
        'id': 0,            # #ID
        'chr': 1,          # CHR
        'pos': 2,          # POS
        'ref': 3,          # REF
        'alt': 4,          # ALT
        'beta': 5,         # BETA
        'se': 6,           # SE
        'pval': 7,         # P
        'aaf': 12,         # AAF
        'n': 16,           # N
        'n_study': 18      # N_STUDY
    },
    'gwama': {
        'id': 0,           # #ID
        'chr': 1,          # CHR
        'pos': 2,          # POS
        'ref': 3,          # REF
        'alt': 4,          # ALT
        'beta': 5,         # BETA
        'se': 6,           # SE
        'pval': 7,         # P
        'aaf': 12,         # AAF
        'n': 16,           # N
        'n_study': 18      # N_STUDY
    }
}

def create_tabix_file(input_file, output_file, software_type, pval_min=None, pval_max=None):
    mapping = COLUMN_MAPPINGS[software_type]
    sorted_file = f"{output_file}.sorted"
    filtered_file = f"{output_file}.filtered"
    print(f"Sorting data for {output_file}...")

    # Modify sorting command to use correct p-value column based on software type
    pval_col = mapping['pval'] + 1  # Add 1 because sort command uses 1-based indexing
    if input_file.endswith('.gz'):
        sort_cmd = f"(zcat {input_file} | head -n 1 && zcat {input_file} | tail -n +2 | sort -k2,2V -k3,3n -k{pval_col},{pval_col}n) > {sorted_file}"
    else:
        sort_cmd = f"(head -n 1 {input_file} && tail -n +2 {input_file} | sort -k2,2V -k3,3n -k{pval_col},{pval_col}n) > {sorted_file}"
    subprocess.run(sort_cmd, shell=True, check=True)

    # Filter based on p-value range
    print(f"Filtering data with p-value in range {pval_min} to {pval_max}...")

    with open(sorted_file, 'r') as sorted_fh, open(filtered_file, 'w') as filtered_fh:
        header = sorted_fh.readline()
        filtered_fh.write(header)  # Keep the header
        for line in sorted_fh:
            fields = line.strip().split('\t')
            try:
                pval = float(fields[mapping['pval']])  # Use correct p-value column from mapping
                if (pval_min is None or pval >= pval_min) and (pval_max is None or pval <= pval_max):
                    filtered_fh.write(line)
            except (ValueError, IndexError) as e:
                print(f"Skipping line due to error: {e}")
    
    # Compress with bgzip
    output_gz_file = f"{output_file}.gz"
    print(f"Compressing filtered file to {output_gz_file} using bgzip...")
    subprocess.run(['bgzip', '-c', filtered_file], stdout=open(output_gz_file, 'wb'), check=True)

    # Create Tabix index
    print(f"Creating Tabix index for {output_gz_file}...")
    subprocess.run(['tabix', '-s', '2', '-b', '3', '-e', '3', '-S', '1', output_gz_file], check=True)

    # Convert to base64
    print("Encoding files to base64...")
    with open(output_gz_file, "rb") as f:
        encoded_gz = base64.b64encode(f.read()).decode('utf-8')
    with open(f"{output_gz_file}.tbi", "rb") as f:
        encoded_tbi = base64.b64encode(f.read()).decode('utf-8')

    print(f"Compressed file and index created: {output_gz_file} and {output_gz_file}.tbi")

    # Clean up
    os.remove(sorted_file)
    os.remove(filtered_file)

    return encoded_gz, encoded_tbi

def main():
    parser = argparse.ArgumentParser(description='Process GWAS summary statistics files and create tabix indexes.')
    parser.add_argument('-i', '--input_folder', required=True,
                        help='Input folder containing the summary statistics files')
    parser.add_argument('-o', '--output_folder', required=True,
                        help='Output folder for processed files')
    parser.add_argument('-s', '--software', choices=['mrmega', 'gwama', 'both'], default='both',
                        help='Software type to process (mrmega, gwama, or both) (default: both)')
    parser.add_argument('-p', '--phenotype', 
                        help='Specific phenotype to process (if not specified, processes all phenotypes)')
    parser.add_argument('-c', '--cohort',
                        help='Specific cohort to process (if not specified, processes all cohorts)')
    
    args = parser.parse_args()

    os.makedirs(args.output_folder, exist_ok=True)

    # Define p-value ranges
    pval_ranges = [
        (None, 1e-5),  # p-values <= 1e-4
        (None, 1e-1)   # p-values <= 1e-1
    ]

    # Define software types to process
    software_types = ["mrmega", "gwama"] if args.software == 'both' else [args.software]

    for software in software_types:
        # Construct file pattern
        if args.phenotype and args.cohort:
            pattern = f"{args.phenotype}.{args.cohort}.{software}.sumstats.txt.gz"
        elif args.phenotype:
            pattern = f"{args.phenotype}.*.{software}.sumstats.txt.gz"
        elif args.cohort:
            pattern = f"*.{args.cohort}.{software}.sumstats.txt.gz"
        else:
            pattern = f"*.{software}.sumstats.txt.gz"
        
        full_pattern = os.path.join(args.input_folder, pattern)
        sumstats_files = glob.glob(full_pattern)
        
        if not sumstats_files:
            print(f"No files found matching pattern: {pattern}")
            continue
        
        for input_file in sumstats_files:
            filename = os.path.basename(input_file)
            components = filename.split('.')
            phenotype = components[0]
            cohort = components[1]
            
            for pval_min, pval_max in pval_ranges:
                if pval_max is None:
                    pval_range_str = f"{pval_min}_and_above"
                elif pval_min is None:
                    pval_range_str = f"up_to_{pval_max}"
                else:
                    pval_range_str = f"{pval_min}_to_{pval_max}"
                
                output_file = os.path.join(
                    args.output_folder,
                    f"{phenotype}.{cohort}.{software}_pval_{pval_range_str}"
                )
                
                print(f"Processing {phenotype}-{cohort}-{software} for p-value range {pval_min} to {pval_max if pval_max else 'above'}...")
                encoded_gz, encoded_tbi = create_tabix_file(input_file, output_file, software, pval_min, pval_max)

if __name__ == "__main__":
    main()
