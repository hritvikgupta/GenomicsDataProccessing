import subprocess
import os
import base64
import glob
import argparse

def create_tabix_file(input_file, output_file, pval_min=None, pval_max=None):
    # Step 1: Sort the input file using Unix sort
    sorted_file = f"{output_file}.sorted"
    filtered_file = f"{output_file}.filtered"
    print(f"Sorting data for {output_file}...")

    # Modify the sorting command to sort by chromosome, position, and pval
    if input_file.endswith('.gz'):
        sort_cmd = f"(zcat {input_file} | head -n 1 && zcat {input_file} | tail -n +2 | sort -k2,2V -k3,3n -k15,15n) > {sorted_file}"
    else:
        sort_cmd = f"(head -n 1 {input_file} && tail -n +2 {input_file} | sort -k2,2V -k3,3n -k15,15n) > {sorted_file}"
    subprocess.run(sort_cmd, shell=True, check=True)

    # Step 2: Filter the sorted file based on p-value range
    print(f"Filtering data with p-value in range {pval_min} to {pval_max}...")

    with open(sorted_file, 'r') as sorted_fh, open(filtered_file, 'w') as filtered_fh:
        header = sorted_fh.readline()
        filtered_fh.write(header)  # Keep the header
        for line in sorted_fh:
            fields = line.strip().split('\t')
            try:
                pval = float(fields[14])  # Assuming p-value is in the 15th column (index 14)

                # Check if the pval falls within the specified range
                if (pval_min is None or pval >= pval_min) and (pval_max is None or pval <= pval_max):
                    filtered_fh.write(line)
            except (ValueError, IndexError) as e:
                print(f"Skipping line due to error: {e}")
    
    # Step 3: Compress the filtered file using bgzip
    output_gz_file = f"{output_file}.gz"
    print(f"Compressing filtered file to {output_gz_file} using bgzip...")
    subprocess.run(['bgzip', '-c', filtered_file], stdout=open(output_gz_file, 'wb'), check=True)

    # Step 4: Create a Tabix index
    print(f"Creating Tabix index for {output_gz_file}...")

    # Tabix indexing on chromosome (-s 2), position (-b 3, -e 3)
    subprocess.run(['tabix', '-s', '2', '-b', '3', '-e', '3', '-S', '1', output_gz_file], check=True)

    # Step 5: Convert files to base64 to store in Cosmos DB (or wherever needed)
    print("Encoding files to base64...")
    with open(output_gz_file, "rb") as f:
        encoded_gz = base64.b64encode(f.read()).decode('utf-8')
    with open(f"{output_gz_file}.tbi", "rb") as f:
        encoded_tbi = base64.b64encode(f.read()).decode('utf-8')

    print(f"Compressed file and index created: {output_gz_file} and {output_gz_file}.tbi")

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
        (None, 1e-5),  # p-values <= 1e-4
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
