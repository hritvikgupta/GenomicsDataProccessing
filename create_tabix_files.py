import subprocess
import os
import base64
import glob
import argparse

def create_tabix_file(input_file, output_file, pval_min=None, pval_max=None):
    # Add counters for debugging
    total_lines = 0
    filtered_lines = 0
    error_lines = 0
    
    # Step 1: First count lines in input file
    print(f"Counting lines in input file {input_file}...")
    if input_file.endswith('.gz'):
        wc_cmd = f"gunzip -c {input_file} | wc -l"
    else:
        wc_cmd = f"wc -l < {input_file}"
    
    try:
        total_input_lines = int(subprocess.check_output(wc_cmd, shell=True).decode().strip())
        print(f"Total input lines: {total_input_lines}")
    except subprocess.CalledProcessError as e:
        print(f"Error counting input lines: {e}")
        raise

    # Step 2: Sort the input file
    sorted_file = f"{output_file}.sorted"
    filtered_file = f"{output_file}.filtered"
    print(f"Sorting data for {output_file}...")

    # Modified sort command to sort by CHR, POS, and P columns (2, 3, and 8)
    if input_file.endswith('.gz'):
        sort_cmd = f"(gunzip -c {input_file} | head -n 1 && gunzip -c {input_file} | tail -n +2 | sort -k2,2V -k3,3n -k8,8n) > {sorted_file}"
    else:
        sort_cmd = f"(head -n 1 {input_file} && tail -n +2 {input_file} | sort -k2,2V -k3,3n -k8,8n) > {sorted_file}"
    
    try:
        subprocess.run(sort_cmd, shell=True, check=True)
        sorted_size = os.path.getsize(sorted_file)
        print(f"Sorted file size: {sorted_size} bytes")
    except subprocess.CalledProcessError as e:
        print(f"Error during sorting: {e}")
        raise

    # Step 3: Filter the sorted file
    print(f"Filtering data with p-value in range {pval_min} to {pval_max}...")
    
    with open(sorted_file, 'r') as sorted_fh, open(filtered_file, 'w') as filtered_fh:
        header = sorted_fh.readline()
        filtered_fh.write(header)
        filtered_lines += 1  # Count header
        
        for line_num, line in enumerate(sorted_fh, 1):
            total_lines += 1
            fields = line.strip().split('\t')
            try:
                if len(fields) < 8:  # Verify we have enough fields
                    print(f"Warning: Line {line_num} has only {len(fields)} fields")
                    error_lines += 1
                    continue
                    
                pval = float(fields[7])  # P-value is in column 8 (index 7)
                
                # Print some sample p-values for debugging
                if line_num <= 5:
                    print(f"Sample p-value at line {line_num}: {pval}")
                
                # Validate p-value is between 0 and 1
                if not (0 <= pval <= 1):
                    error_lines += 1
                    if error_lines <= 5:
                        print(f"Invalid p-value at line {line_num}: {pval}")
                    continue
                
                if (pval_min is None or pval >= pval_min) and (pval_max is None or pval <= pval_max):
                    filtered_fh.write(line)
                    filtered_lines += 1
                    
            except (ValueError, IndexError) as e:
                error_lines += 1
                if error_lines <= 5:  # Only print first 5 errors
                    print(f"Error at line {line_num}: {e}")
                    print(f"Problematic line: {line[:100]}...")
    
    print(f"\nFiltering Statistics:")
    print(f"Total lines processed: {total_lines}")
    print(f"Lines after filtering: {filtered_lines}")
    print(f"Error lines: {error_lines}")
    
    # Rest of the function remains the same...
    # Step 4: Compress with bgzip
    output_gz_file = f"{output_file}.gz"
    print(f"Compressing filtered file to {output_gz_file} using bgzip...")
    subprocess.run(['bgzip', '-c', filtered_file], stdout=open(output_gz_file, 'wb'), check=True)
    
    # Create Tabix index
    print(f"Creating Tabix index for {output_gz_file}...")
    try:
        subprocess.run(['tabix', '-s', '2', '-b', '3', '-e', '3', '-S', '1', output_gz_file], check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error creating tabix index: {e}")
        raise

    # Base64 encoding
    print("Encoding files to base64...")
    with open(output_gz_file, "rb") as f:
        encoded_gz = base64.b64encode(f.read()).decode('utf-8')
    with open(f"{output_gz_file}.tbi", "rb") as f:
        encoded_tbi = base64.b64encode(f.read()).decode('utf-8')

    print(f"\nFinal file sizes:")
    print(f"Original input: {os.path.getsize(input_file)} bytes")
    print(f"Sorted file: {os.path.getsize(sorted_file)} bytes")
    print(f"Filtered file: {os.path.getsize(filtered_file)} bytes")
    print(f"Compressed output: {os.path.getsize(output_gz_file)} bytes")
    print(f"Index file: {os.path.getsize(f'{output_gz_file}.tbi')} bytes")

    # Clean up temporary files
    os.remove(sorted_file)
    os.remove(filtered_file)

    return encoded_gz, encoded_tbi

def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description='Process GWAS summary statistics files and create tabix indexes.')
    
    # Required arguments
    parser.add_argument('-i', '--input_folder', required=True,
                        help='Input folder containing the summary statistics files')
    parser.add_argument('-o', '--output_folder', required=True,
                        help='Output folder for processed files')
    
    # Optional arguments
    parser.add_argument('-s', '--software', choices=['mrmega', 'gwama', 'both'], default='both',
                        help='Software type to process (mrmega, gwama, or both) (default: both)')
    parser.add_argument('-p', '--phenotype', 
                        help='Specific phenotype to process (if not specified, processes all phenotypes)')
    parser.add_argument('-c', '--cohort',
                        help='Specific cohort to process (if not specified, processes all cohorts)')
    
    args = parser.parse_args()

    # Create output directory if it doesn't exist
    os.makedirs(args.output_folder, exist_ok=True)

    # Define the p-value ranges
    pval_ranges = [
        (None, 1e-4),  # p-values <= 1e-4
        (None, 1e-1)   # p-values <= 1e-1
    ]

    # Define which software types to process
    if args.software == 'both':
        software_types = ["mrmega", "gwama"]
    else:
        software_types = [args.software]

    # Process each software type
    for software in software_types:
        # Construct the file pattern based on arguments
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
        
        # Process each file
        for input_file in sumstats_files:
            # Extract phenotype and cohort from filename
            filename = os.path.basename(input_file)
            components = filename.split('.')
            phenotype = components[0]
            cohort = components[1]
            
            # Process each p-value range
            for pval_min, pval_max in pval_ranges:
                if pval_max is None:
                    pval_range_str = f"{pval_min}_and_above"
                elif pval_min is None:
                    pval_range_str = f"up_to_{pval_max}"
                else:
                    pval_range_str = f"{pval_min}_to_{pval_max}"
                
                # Construct the output filename
                output_file = os.path.join(
                    args.output_folder,
                    f"{phenotype}.{cohort}.{software}_pval_{pval_range_str}"
                )
                
                print(f"Processing {phenotype}-{cohort}-{software} for p-value range {pval_min} to {pval_max if pval_max else 'above'}...")
                encoded_gz, encoded_tbi = create_tabix_file(input_file, output_file, pval_min, pval_max)

if __name__ == "__main__":
    main()
