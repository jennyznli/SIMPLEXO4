import os
import pandas as pd
import argparse

def get_args():
    p = argparse.ArgumentParser(description='Preprocess data for regenie analysis')
    p.add_argument('--vep-file', required=True, help='Single VEP CSV file containing all chromosome annotations')
    p.add_argument('--samples', required=True, help='Samples list file (one IID per line)')
    p.add_argument('--controls', help='Controls list file (one IID per line). Required unless --pheno-file is given.')
    p.add_argument('--cases', help='Cases list file (one IID per line). Required unless --pheno-file is given.')
    p.add_argument('--pheno-file', help='Pre-built multi-phenotype file with a header row (must include an IID column, '
                                         'FID optional, plus one or more phenotype columns coded 1/0/NA). If given, '
                                         'this is used directly instead of building a single STATUS column from --cases/--controls.')
    p.add_argument('--step1-covariates', help='Step 1 covariates file (IID and additional columns)')
    p.add_argument('--step2-covariates', help='Step 2 covariates file (IID and additional columns)')
    p.add_argument('--gene-consequence-exclude', default='',
                    help='TSV: GENE<TAB>CONSEQUENCE_SUBSTRING (e.g. CHEK2 missense_variant). '
                         'Drops matching sites from the annotation/set/mask files before pathogenic/VUS filtering.')
    p.add_argument('--negative-control-blacklist', default='',
                    help='Optional file with variant IDs (one per line) to exclude from the M4 synonymous negative-control mask.')
    p.add_argument('--mask-file', default='',
                    help='Optional custom regenie mask definition file (lines of "Mx category1,category2,..."). '
                         'Categories must be a subset of {pathogenic, vus, synonymous} (the labels written to '
                         'annotation.txt). If omitted, falls back to the built-in M1 pathogenic / M2 pathogenic,vus default.')
    p.add_argument('-O', '--output-prefix', default='', help='Output prefix (optional)')
    args = p.parse_args()
    if not args.pheno_file and not (args.cases and args.controls):
        p.error('must provide either --pheno-file, or both --cases and --controls')
    return args

def load_gene_consequence_exclude(path):
    """Load GENE<TAB>CONSEQUENCE_SUBSTRING exclusion rules. Match is case-insensitive substring on
    Variant.Consequence, so compound VEP terms like missense_variant&splice_region_variant still match."""
    rules = []
    if not path:
        return rules
    if not os.path.exists(path):
        print(f"Warning: gene-consequence exclude not found: {path} (continuing without)")
        return rules
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split('\t')
            if len(parts) < 2:
                parts = line.split()
            if len(parts) < 2:
                print(f"Warning: skip bad exclude rule: {line}")
                continue
            gene = parts[0].strip().upper()
            cons = parts[1].strip().lower()
            if gene and cons:
                rules.append((gene, cons))
    return rules

def matches_gene_consequence_exclude(gene, consequence, rules):
    if not rules:
        return False
    g = str(gene).strip().upper()
    c = str(consequence).strip().lower() if pd.notna(consequence) else ''
    if not g or not c:
        return False
    for gene_u, cons_sub in rules:
        if g == gene_u and cons_sub in c:
            return True
    return False

VALID_MASK_CATEGORIES = {'pathogenic', 'vus', 'synonymous'}

def load_mask_file(path, available_categories):
    """Load a custom regenie mask definition file (lines of 'Mx category1,category2,...').
    Validates each category against VALID_MASK_CATEGORIES (typo/format check) and warns if a
    referenced category has zero variants in this run (mask would come out empty)."""
    if not path:
        return None
    if not os.path.exists(path):
        print(f"Warning: mask file not found: {path} (falling back to built-in default mask)")
        return None
    lines = []
    with open(path) as f:
        for raw_line in f:
            line = raw_line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split(None, 1)
            if len(parts) != 2:
                print(f"Warning: skipping malformed mask line: {raw_line!r}")
                continue
            mask_name, cats_str = parts
            categories = [c.strip() for c in cats_str.split(',') if c.strip()]
            unknown = [c for c in categories if c not in VALID_MASK_CATEGORIES]
            if unknown:
                print(f"Warning: mask {mask_name} references unrecognized category(ies) {unknown} "
                      f"(valid: {sorted(VALID_MASK_CATEGORIES)}) -- possible typo")
            empty = [c for c in categories if c in VALID_MASK_CATEGORIES and c not in available_categories]
            if empty:
                print(f"Warning: mask {mask_name} references category(ies) {empty} with zero variants "
                      f"in this run -- mask may come out empty")
            lines.append(line)
    return lines

def is_synonymous_consequence(value):
    if pd.isna(value):
        return False
    text = str(value).strip().lower()
    return bool(text) and 'synonymous_variant' in text

def load_sample_list(filepath):
    """Load sample list from file, handling potential headers"""
    df = pd.read_csv(filepath, sep='\t', header=None)
    # Check if first row looks like a header
    if len(df) > 0 and df.iloc[0, 0] in ['IID', 'FID', '#IID', '#FID']:
        samples = df.iloc[1:, 0].tolist()
        print(f"  Detected header in {filepath}, skipping first row")
    else:
        samples = df[0].tolist()
    return samples

def load_pheno_file(filepath):
    """Load a pre-built multi-phenotype file. Must have a header row with an 'IID' column
    (FID added if missing) plus one or more phenotype columns coded 1/0/NA."""
    df = pd.read_csv(filepath, sep=r'\s+')
    if 'IID' not in df.columns:
        raise ValueError(f"Phenotype file {filepath} must have a header row including an 'IID' column "
                          f"(found columns: {list(df.columns)})")
    if 'FID' not in df.columns:
        df.insert(0, 'FID', df['IID'])
    return df

def main():
    args = get_args()
    
    # Set up output prefix
    prefix = f"{args.output_prefix}." if args.output_prefix else ""
    
    # Load samples
    print("Loading sample lists...")
    all_samples = set(load_sample_list(args.samples))
    print(f"  Total samples: {len(all_samples)}")

    pheno_df = None
    if args.pheno_file:
        print(f"\nLoading multi-phenotype file: {args.pheno_file}")
        pheno_df = load_pheno_file(args.pheno_file)
        pheno_cols = [c for c in pheno_df.columns if c not in ('FID', 'IID')]
        print(f"  Phenotype columns: {', '.join(pheno_cols)}")
        final_samples = all_samples & set(pheno_df['IID'])
        print(f"  Samples in both --samples and --pheno-file: {len(final_samples)}")
        valid_controls = set()
        valid_cases = set()
    else:
        all_controls = set(load_sample_list(args.controls))
        all_cases = set(load_sample_list(args.cases))
        print(f"  Total controls: {len(all_controls)}")
        print(f"  Total cases: {len(all_cases)}")

        # Get intersection of samples with controls and cases
        valid_controls = all_samples & all_controls
        valid_cases = all_samples & all_cases
        final_samples = valid_controls | valid_cases

        print(f"\nAfter intersection:")
        print(f"  Valid controls: {len(valid_controls)}")
        print(f"  Valid cases: {len(valid_cases)}")
        print(f"  Final samples (union): {len(final_samples)}")

        # Check for overlap between cases and controls
        overlap = valid_controls & valid_cases
        if overlap:
            print(f"  WARNING: {len(overlap)} samples appear in both cases and controls!")

    # Check that VEP file exists
    if not os.path.exists(args.vep_file):
        print(f"Error: VEP file does not exist: {args.vep_file}")
        return

    print(f"\nProcessing VEP file: {args.vep_file}")

    # Process VEP CSV file
    required_columns = ['ID', 'Gene', 'Variant.LoF_level', 'Variant.Consequence']

    print("  Reading VEP file...")
    df = pd.read_csv(args.vep_file, usecols=required_columns)
    print(f"  Loaded {len(df)} total variants")

    # Gene x consequence exclusion (e.g. CHEK2 missense_variant), applied before LoF-level filtering
    exclude_rules = load_gene_consequence_exclude(args.gene_consequence_exclude)
    if exclude_rules:
        print(f"  Gene x consequence exclude rules: {len(exclude_rules)}")
        for g, c in exclude_rules:
            print(f"    {g} + {c}")
        hit = df.apply(lambda r: matches_gene_consequence_exclude(r['Gene'], r['Variant.Consequence'], exclude_rules), axis=1)
        excluded_df = df[hit]
        df = df[~hit]
        excl_file = f'{prefix}gene_consequence_excluded.csv'
        excluded_df.to_csv(excl_file, index=False)
        print(f"  Excluded {len(excluded_df)} gene x consequence sites from all masks; wrote {excl_file}")
    else:
        print("  No gene x consequence exclude rules applied")

    # Filter for pathogenic and VUS variants (levels 1 and 2)
    filtered_df = df[df['Variant.LoF_level'].isin([1, 2])][['ID', 'Gene', 'Variant.LoF_level']]
    print(f"  Filtered to {len(filtered_df)} pathogenic/VUS variants")
    
    # Deduplicate: for each ID (variant), keep only ONE row with the most pathogenic annotation
    # Sort by LoF_level (1 is most pathogenic, 2 is less pathogenic)
    filtered_df = filtered_df.sort_values('Variant.LoF_level')
    
    # Check for duplicates before deduplication
    duplicate_count = filtered_df['ID'].duplicated().sum()
    if duplicate_count > 0:
        print(f"  Found {duplicate_count} duplicate variant IDs")
    
    # Keep only the first occurrence of each ID (most pathogenic due to sorting)
    combined_df = filtered_df.drop_duplicates(subset=['ID'], keep='first')
    print(f"  Total unique variants after deduplication: {len(combined_df)}")
    
    # Split back into pathogenic and VUS
    pathogenic_df = combined_df[combined_df['Variant.LoF_level'] == 1][['ID', 'Gene']]
    vus_df = combined_df[combined_df['Variant.LoF_level'] == 2][['ID', 'Gene']]

    print(f"  Pathogenic variants: {len(pathogenic_df)}")
    print(f"  VUS variants: {len(vus_df)}")

    # M4: synonymous negative-control mask. Built from the same (post gene-consequence-exclude) df,
    # excluding anything already claimed by pathogenic/VUS and any blacklisted (e.g. potential MNP-pair) IDs.
    blacklist_ids = set()
    if args.negative_control_blacklist and os.path.exists(args.negative_control_blacklist):
        with open(args.negative_control_blacklist) as f:
            blacklist_ids = {line.strip() for line in f if line.strip()}
        print(f"  Loaded {len(blacklist_ids)} blacklisted variant IDs for M4 negative control")
    elif args.negative_control_blacklist:
        print(f"  Warning: negative-control blacklist not found: {args.negative_control_blacklist} (continuing without)")

    claimed_ids = set(pathogenic_df['ID']) | set(vus_df['ID'])
    syn_base = df[df['Variant.Consequence'].apply(is_synonymous_consequence)]
    syn_base = syn_base[~syn_base['ID'].isin(claimed_ids)]
    if blacklist_ids:
        syn_base = syn_base[~syn_base['ID'].isin(blacklist_ids)]
    syn_df = syn_base[['ID', 'Gene']].drop_duplicates(subset=['ID'])
    print(f"  Synonymous negative-control variants: {len(syn_df)}")

    # Write annotation file
    annotation_file = f'{prefix}annotation.txt'
    with open(annotation_file, 'w') as ann_f:
        for _, row in pathogenic_df.iterrows():
            ann_f.write(f"{row['ID']} {row['Gene']} pathogenic\n")
        for _, row in vus_df.iterrows():
            ann_f.write(f"{row['ID']} {row['Gene']} vus\n")
        for _, row in syn_df.iterrows():
            ann_f.write(f"{row['ID']} {row['Gene']} synonymous\n")

    print(f"\nWritten {len(pathogenic_df) + len(vus_df) + len(syn_df)} unique variant annotations to {annotation_file}")

    # Write set file (group by gene) -- pathogenic + VUS + synonymous all need to appear here
    # so regenie can build M1-M4 masks from the categories declared in annotation.txt
    set_file = f'{prefix}set.txt'
    set_variants_df = pd.concat([pathogenic_df, vus_df, syn_df], ignore_index=True).drop_duplicates(subset=['ID'])
    gene_groups = set_variants_df.groupby('Gene')['ID'].agg(list)

    print(f"  Writing gene sets for {len(gene_groups)} genes...")
    with open(set_file, 'w') as set_f:
        for gene, variants in gene_groups.items():
            # Get chromosome and position from first variant
            first_variant = variants[0]
            chr_pos = first_variant.split('_')[:2]
            chr_name = chr_pos[0]
            pos = int(chr_pos[1])

            # Format variant list
            variant_str = ','.join(variants)
            set_f.write(f"{gene} {chr_name} {pos} {variant_str}\n")

    print(f"Written gene sets to {set_file}")

    # Create mask file: custom --mask-file if provided (validated against categories actually
    # present in this run), else the built-in default.
    available_categories = set()
    if len(pathogenic_df) > 0:
        available_categories.add('pathogenic')
    if len(vus_df) > 0:
        available_categories.add('vus')
    if len(syn_df) > 0:
        available_categories.add('synonymous')

    custom_mask_lines = load_mask_file(args.mask_file, available_categories)

    mask_file = f'{prefix}mask.txt'
    with open(mask_file, 'w') as f:
        if custom_mask_lines is not None:
            print(f"  Using custom mask file: {args.mask_file}")
            for line in custom_mask_lines:
                f.write(f"{line}\n")
        else:
            f.write('M1 pathogenic\n')
            # f.write('M2 vus\n')
            f.write('M2 pathogenic,vus\n')
            # f.write('M4 synonymous\n')

    print(f"Written mask definitions to {mask_file}")
    
    # Create covariate files for Step 1 and Step 2
    print("\nCreating covariate files...")
    final_samples_list = sorted(list(final_samples))
    
    # Process Step 1 covariates
    if args.step1_covariates:
        print(f"  Processing Step 1 covariates: {args.step1_covariates}")
        step1_covar_df = pd.read_csv(args.step1_covariates, sep='\s+')
        
        # Ensure IID column exists
        if 'IID' not in step1_covar_df.columns:
            print("Error: Step 1 covariates file must contain 'IID' column")
            return
        
        # Filter to only include final samples
        step1_covar_df = step1_covar_df[step1_covar_df['IID'].isin(final_samples)]
        
        # Add FID if not present
        if 'FID' not in step1_covar_df.columns:
            step1_covar_df.insert(0, 'FID', step1_covar_df['IID'])
        
        print(f"    Retained {len(step1_covar_df)} samples with {len(step1_covar_df.columns)-2} covariates")
    else:
        # Create minimal covar file with just FID and IID
        step1_covar_df = pd.DataFrame({'FID': final_samples_list, 'IID': final_samples_list})
        print(f"    No Step 1 covariates provided, created minimal file with FID/IID only")

    # Write Step 1 covar file
    step1_covar_file = f'{prefix}covar.step1.txt'
    step1_covar_df.to_csv(step1_covar_file, sep=' ', index=False, na_rep='NA')
    print(f"  Written Step 1 covariates to {step1_covar_file}")
    
    # Process Step 2 covariates
    if args.step2_covariates:
        print(f"  Processing Step 2 covariates: {args.step2_covariates}")
        step2_covar_df = pd.read_csv(args.step2_covariates, sep='\s+')
        
        # Ensure IID column exists
        if 'IID' not in step2_covar_df.columns:
            print("Error: Step 2 covariates file must contain 'IID' column")
            return
        
        # Filter to only include final samples
        step2_covar_df = step2_covar_df[step2_covar_df['IID'].isin(final_samples)]
        
        # Add FID if not present
        if 'FID' not in step2_covar_df.columns:
            step2_covar_df.insert(0, 'FID', step2_covar_df['IID'])
        
        print(f"    Retained {len(step2_covar_df)} samples with {len(step2_covar_df.columns)-2} covariates")
    else:
        # Create minimal covar file with just FID and IID
        step2_covar_df = pd.DataFrame({'FID': final_samples_list, 'IID': final_samples_list})
        print(f"    No Step 2 covariates provided, created minimal file with FID/IID only")
    
    # Write Step 2 covar file
    step2_covar_file = f'{prefix}covar.step2.txt'
    step2_covar_df.to_csv(step2_covar_file, sep=' ', index=False, na_rep='NA')
    print(f"  Written Step 2 covariates to {step2_covar_file}")
    
    # Create phenotype file
    print("\nCreating phenotype file...")
    pheno_file = f'{prefix}pheno.txt'

    if args.pheno_file:
        # Multi-phenotype mode: filter the pre-built pheno file down to final_samples,
        # keep FID/IID first, preserve whatever phenotype columns it already has.
        pheno_df = pheno_df[pheno_df['IID'].isin(final_samples)].sort_values('IID').reset_index(drop=True)
        other_cols = [c for c in pheno_df.columns if c not in ('FID', 'IID')]
        pheno_df = pheno_df[['FID', 'IID'] + other_cols]

        pheno_df.to_csv(pheno_file, sep=' ', index=False, na_rep='NA')

        print(f"  Written multi-phenotype file to {pheno_file}")
        print(f"    Samples: {len(pheno_df)}")
        for col in other_cols:
            n_cases = (pheno_df[col] == 1).sum()
            n_controls = (pheno_df[col] == 0).sum()
            n_missing = pheno_df[col].isna().sum()
            print(f"    {col}: {n_cases} cases, {n_controls} controls, {n_missing} missing")
    else:
        # Single-phenotype mode: STATUS column, 0 for controls, 1 for cases
        pheno_df = pd.DataFrame({'FID': final_samples_list, 'IID': final_samples_list})
        pheno_df['STATUS'] = pheno_df['IID'].apply(lambda x: 0 if x in valid_controls else 1)

        pheno_df.to_csv(pheno_file, sep=' ', index=False, na_rep='NA')

        cases_in_pheno = (pheno_df['STATUS'] == 1).sum()
        controls_in_pheno = (pheno_df['STATUS'] == 0).sum()
        print(f"  Written phenotype file to {pheno_file}")
        print(f"    Cases: {cases_in_pheno}")
        print(f"    Controls: {controls_in_pheno}")

    print("\nPreprocessing complete!")

if __name__ == '__main__':
    main()
