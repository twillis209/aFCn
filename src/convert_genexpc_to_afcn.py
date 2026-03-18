#!/usr/bin/env python3
"""
Convert genexpc formatted files to aFCn format.

genexpc format:
- Expression: haplotype-level (sample_id_0, sample_id_1) with diploid expression values
- eQTL: haplotype-level genotypes with variant IDs as columns (chr:pos:ref:alt)

aFCn format:
- Expression: diploid-level (sample_id) with gene expression values (gzipped CSV)
- eQTL: tab-separated file with gene_id and variant_id columns
- VCF: standard phased VCF file (bgzipped and tabix-indexed)
"""

import pandas as pd
import numpy as np
import argparse
import gzip
import sys
from pathlib import Path


def convert_expression_file(genexpc_exp_path, afcn_exp_path, gene_id, convert_ln_to_log2=True, samples=None):
    """
    Convert genexpc haplotype-level expression to aFCn diploid-level expression.
    
    genexpc format:
        haplotype_id,diploid_expression
        HG00096_0,-0.119
        HG00096_1,-0.119
    
    aFCn format (gzipped):
        Name,HG00096,HG00097,...
        ENSG00000227232,5.2,6.1,...
    
    Note: aFCn requires log2-transformed expression data.
    If genexpc data is in natural log (ln) scale, it will be converted to log2.
    """
    print(f"Reading expression file: {genexpc_exp_path}")
    exp_df = pd.read_csv(genexpc_exp_path, index_col='haplotype_id')
    
    # Extract sample IDs (remove _0 and _1 suffixes)
    sample_data = {}
    
    for hap_id, row in exp_df.iterrows():
        if hap_id.endswith('_0'):
            sample_id = hap_id[:-2]
            if samples is not None and sample_id not in samples:
                continue
            expr_val = row['diploid_expression']
            
            # Convert from natural log to log2 if needed
            # genexpc uses: log(exp(diploid_expression))
            # aFCn requires: log2 scale
            # Conversion: log2(x) = ln(x) / ln(2)
            if convert_ln_to_log2:
                expr_val = expr_val / np.log(2)
            
            sample_data[sample_id] = [expr_val]
    
    # Create output dataframe efficiently with pd.concat
    output_df = pd.DataFrame({'Name': [gene_id]})
    sample_df = pd.DataFrame(sample_data)
    output_df = pd.concat([output_df, sample_df], axis=1)
    
    # Write gzipped CSV
    print(f"Writing expression file: {afcn_exp_path}")
    if convert_ln_to_log2:
        print("  Note: Converting from natural log (ln) to log2 scale (required by aFCn)")
    with gzip.open(afcn_exp_path, 'wt') as f:
        output_df.to_csv(f, index=False)
    
    num_samples = len(sample_data)
    print(f"Converted {num_samples} samples")
    return list(sample_data.keys())


def convert_eqtl_file(genexpc_eqtl_path, afcn_eqtl_path, gene_id, vcf_output_path=None):
    """
    Convert genexpc eQTL file to aFCn eQTL format and optionally create VCF stub.
    
    genexpc format:
        haplotype_id,3:38079424:A:G,3:37995316:A:G,...
        HG00096_0,0,0,...
        HG00096_1,0,0,...
    
    aFCn eQTL format (tab-separated):
        gene_id	variant_id
        ENSG00000227232	chr3_38079424_A_G_b38
    
    Returns: list of variant IDs and dictionary of haplotypes for VCF creation
    """
    print(f"Reading eQTL file: {genexpc_eqtl_path}")
    eqtl_df = pd.read_csv(genexpc_eqtl_path, index_col='haplotype_id')
    
    # Get variant IDs and convert format: chr:pos:ref:alt -> chr#_pos_ref_alt_b38
    genexpc_variants = eqtl_df.columns.tolist()
    afcn_variants = []
    
    for var in genexpc_variants:
        parts = var.split(':')
        if len(parts) == 4:
            chrom, pos, ref, alt = parts
            # Add 'chr' prefix if not present
            if not chrom.startswith('chr'):
                chrom = f'chr{chrom}'
            afcn_var = f"{chrom}_{pos}_{ref}_{alt}_b38"
            afcn_variants.append(afcn_var)
        else:
            print(f"Warning: Skipping malformed variant ID: {var}")
    
    # Create eQTL output file
    eqtl_output_df = pd.DataFrame({
        'gene_id': [gene_id] * len(afcn_variants),
        'variant_id': afcn_variants
    })
    
    print(f"Writing eQTL file: {afcn_eqtl_path}")
    eqtl_output_df.to_csv(afcn_eqtl_path, sep='\t', index=False)
    print(f"Converted {len(afcn_variants)} variants")
    
    # If VCF output requested, create the genotype data structure
    if vcf_output_path:
        return afcn_variants, genexpc_variants, eqtl_df
    
    return afcn_variants, None, None


def create_vcf_from_eqtl(genexpc_eqtl_path, vcf_output_path, afcn_variants, genexpc_variants, eqtl_df, samples=None):
    """
    Create a phased VCF file from genexpc eQTL haplotype data.
    
    VCF format:
        #CHROM  POS     ID                    REF  ALT  QUAL  FILTER  INFO  FORMAT  HG00096  HG00097
        chr3    38079424 chr3_38079424_A_G_b38  A    G    .     .       .     GT      0|0      0|1
    
    Returns path to the compressed and indexed VCF file.
    """
    import pysam
    
    print(f"Creating VCF file: {vcf_output_path}")
    
    # Extract sample IDs (remove _0 and _1 suffixes)
    haplotype_ids = eqtl_df.index.tolist()
    sample_ids = sorted(list(set([hap_id[:-2] for hap_id in haplotype_ids if hap_id.endswith(('_0', '_1'))])))
    if samples is not None:
        sample_ids = [s for s in sample_ids if s in samples]    
    # Parse variant information
    vcf_lines = []
    vcf_lines.append("##fileformat=VCFv4.2")
    vcf_lines.append("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">")
    
    # Header line
    header = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT'] + sample_ids
    vcf_lines.append('\t'.join(header))
    
    # Data lines - need to sort by chromosome and position
    variant_data = []
    for afcn_var, genexpc_var in zip(afcn_variants, genexpc_variants):
        # Parse aFCn variant: chr3_38079424_A_G_b38
        parts = afcn_var.rsplit('_', 1)[0].split('_')  # Remove build suffix first
        chrom = parts[0]
        pos = int(parts[1])  # Convert to int for sorting
        ref = parts[2]
        alt = parts[3]
        
        # Build genotype fields for each sample
        genotypes = []
        for sample_id in sample_ids:
            hap0_id = f"{sample_id}_0"
            hap1_id = f"{sample_id}_1"
            
            # Get haplotype values (0 or 1)
            if hap0_id in eqtl_df.index and hap1_id in eqtl_df.index:
                h0 = int(eqtl_df.loc[hap0_id, genexpc_var])
                h1 = int(eqtl_df.loc[hap1_id, genexpc_var])
                genotypes.append(f"{h0}|{h1}")
            else:
                genotypes.append(".|.")
        
        # Store data for sorting
        variant_data.append((chrom, pos, afcn_var, ref, alt, genotypes))
    
    # Sort by chromosome and position
    variant_data.sort(key=lambda x: (x[0], x[1]))
    
    # Create VCF lines from sorted data
    for chrom, pos, afcn_var, ref, alt, genotypes in variant_data:
        vcf_line = [chrom, str(pos), afcn_var, ref, alt, '.', '.', '.', 'GT'] + genotypes
        vcf_lines.append('\t'.join(vcf_line))
    
    # Determine output paths
    if vcf_output_path.endswith('.gz'):
        vcf_temp = vcf_output_path[:-3]  # Remove .gz for temp file
        vcf_compressed = vcf_output_path
    else:
        vcf_temp = vcf_output_path
        vcf_compressed = vcf_output_path + '.gz'
    
    # Write temporary uncompressed VCF
    print(f"Writing VCF: {vcf_temp}")
    with open(vcf_temp, 'w') as f:
        f.write('\n'.join(vcf_lines))
    
    print(f"Converted {len(sample_ids)} samples and {len(afcn_variants)} variants to VCF")
    
    # Compress with pysam.tabix_compress
    print(f"Compressing with pysam.tabix_compress...")
    try:
        pysam.tabix_compress(vcf_temp, vcf_compressed, force=True)
        print(f"Created: {vcf_compressed}")
        
        # Remove uncompressed file
        Path(vcf_temp).unlink()
    except Exception as e:
        print(f"Error compressing VCF: {e}")
        print(f"Uncompressed VCF available at: {vcf_temp}")
        return vcf_temp
    
    # Index with pysam.tabix_index
    print(f"Indexing with pysam.tabix_index...")
    try:
        pysam.tabix_index(vcf_compressed, preset='vcf', force=True)
        print(f"Created index: {vcf_compressed}.tbi")
    except Exception as e:
        print(f"Error indexing VCF: {e}")
        print(f"VCF is compressed but not indexed: {vcf_compressed}")
    
    return vcf_compressed


def main():
    parser = argparse.ArgumentParser(
        description='Convert genexpc formatted files to aFCn format',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Convert expression and eQTL files
  %(prog)s --exp-in genexpc_exp.csv --exp-out afcn_exp.csv.gz \\
           --eqtl-in genexpc_eqtl.csv --eqtl-out afcn_eqtl.tsv \\
           --gene-id ENSG00000227232
  
  # Also create VCF file
  %(prog)s --exp-in genexpc_exp.csv --exp-out afcn_exp.csv.gz \\
           --eqtl-in genexpc_eqtl.csv --eqtl-out afcn_eqtl.tsv \\
           --vcf-out output.vcf --gene-id ENSG00000227232
        """
    )
    
    parser.add_argument('--exp-in', required=True,
                        help='Input genexpc expression CSV file')
    parser.add_argument('--exp-out', required=True,
                        help='Output aFCn expression CSV.GZ file')
    parser.add_argument('--eqtl-in', required=True,
                        help='Input genexpc eQTL CSV file')
    parser.add_argument('--eqtl-out', required=True,
                        help='Output aFCn eQTL TSV file')
    parser.add_argument('--vcf-out', required=False, default=None,
                        help='Output VCF file path (default: auto-generated from eqtl-out path)')
    parser.add_argument('--gene-id', required=True,
                        help='Gene ID to use (e.g., ENSG00000227232)')
    parser.add_argument('--no-log-conversion', action='store_true',
                        help='Do not convert from ln to log2 (if your data is already in log2)')
    parser.add_argument('--training-test-csv', required=False, default=None,
                        help='CSV with SampleID and is_training columns; restricts output to training samples')
    parser.add_argument('--skip-vcf', action='store_true',
                        help='Do not create VCF file (not recommended - aFCn requires VCF)')
    
    args = parser.parse_args()
    
    # Auto-generate VCF output path if not provided
    if args.vcf_out is None and not args.skip_vcf:
        # Generate from eqtl_out path: test_eqtl.tsv -> test_eqtl.vcf.gz
        eqtl_path = Path(args.eqtl_out)
        vcf_name = eqtl_path.stem + '.vcf.gz'
        args.vcf_out = str(eqtl_path.parent / vcf_name)
    
    # Validate inputs
    if not Path(args.exp_in).exists():
        print(f"Error: Expression file not found: {args.exp_in}", file=sys.stderr)
        sys.exit(1)
    
    if not Path(args.eqtl_in).exists():
        print(f"Error: eQTL file not found: {args.eqtl_in}", file=sys.stderr)
        sys.exit(1)
    
    # Load training samples if requested
    training_samples = None
    if args.training_test_csv is not None:
        tt_df = pd.read_csv(args.training_test_csv, index_col='SampleID')
        training_samples = set(tt_df.index[tt_df['is_training']].tolist())
        print(f"Filtering to {len(training_samples)} training samples from {args.training_test_csv}")

    # Convert expression file
    print("\n=== Converting Expression File ===")
    sample_ids = convert_expression_file(
        args.exp_in, args.exp_out, args.gene_id, 
        convert_ln_to_log2=not args.no_log_conversion,
        samples=training_samples,
    )
    
    # Convert eQTL file
    print("\n=== Converting eQTL File ===")
    afcn_variants, genexpc_variants, eqtl_df = convert_eqtl_file(
        args.eqtl_in, args.eqtl_out, args.gene_id, 
        vcf_output_path=args.vcf_out if not args.skip_vcf else None
    )
    
    # Create VCF (unless skipped)
    vcf_file = None
    if not args.skip_vcf and eqtl_df is not None:
        print("\n=== Creating VCF File ===")
        vcf_file = create_vcf_from_eqtl(
            args.eqtl_in, args.vcf_out, 
            afcn_variants, genexpc_variants, eqtl_df,
            samples=training_samples,
        )
    
    print("\n=== Conversion Complete ===")
    print(f"✓ Expression file: {args.exp_out}")
    print(f"✓ eQTL file: {args.eqtl_out}")
    if vcf_file:
        print(f"✓ VCF file: {vcf_file}")
        print(f"✓ VCF index: {vcf_file}.tbi")
        print(f"\nReady to run aFCn:")
        print(f"  python3 /path/to/afcn.py \\")
        print(f"    --vcf {vcf_file} \\")
        print(f"    --expr {args.exp_out} \\")
        print(f"    --eqtl {args.eqtl_out} \\")
        print(f"    --output results.csv \\")
        print(f"    --nthreads 4")
    else:
        print("\nNote: VCF file was not created. You will need a VCF file to run aFCn.")


if __name__ == "__main__":
    main()
