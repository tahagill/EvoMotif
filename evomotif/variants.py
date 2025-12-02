"""
Variant analysis module.

Analyzes pathogenic variants from ClinVar and tests for enrichment
in conserved motifs.
"""

import numpy as np
import pandas as pd
import logging
from typing import Dict, List, Optional, Tuple
from pathlib import Path
from scipy.stats import fisher_exact, chi2_contingency
from statsmodels.stats.proportion import proportion_confint
import json

logger = logging.getLogger(__name__)


class VariantAnalyzer:
    """Analyze variant enrichment in motifs."""
    
    def __init__(self):
        """Initialize variant analyzer."""
        self.logger = logging.getLogger(__name__)
    
    def load_clinvar_variants(
        self,
        clinvar_file: Path,
        protein_id: Optional[str] = None
    ) -> pd.DataFrame:
        """
        Load and filter ClinVar variants.
        
        Expects TSV/CSV with columns:
        - position: residue position
        - clinical_significance: pathogenic/benign/VUS
        - variant: amino acid change
        
        Args:
            clinvar_file: Path to ClinVar data file
            protein_id: Optional protein ID to filter
            
        Returns:
            DataFrame of variants
        """
        # Determine file format
        if clinvar_file.suffix == '.csv':
            df = pd.read_csv(clinvar_file)
        else:
            df = pd.read_csv(clinvar_file, sep='\t')
        
        # Filter by protein if specified
        if protein_id and 'protein_id' in df.columns:
            df = df[df['protein_id'] == protein_id]
        
        # Standardize clinical significance categories
        df['clinical_significance'] = df['clinical_significance'].str.lower()
        
        self.logger.info(f"Loaded {len(df)} variants from ClinVar")
        return df
    
    def classify_variants(
        self,
        variants: pd.DataFrame
    ) -> Dict[str, pd.DataFrame]:
        """
        Classify variants into pathogenic, benign, and uncertain.
        
        Args:
            variants: DataFrame of variants
            
        Returns:
            Dictionary of classified variant DataFrames
        """
        # Define classification keywords
        pathogenic_keywords = ['pathogenic', 'likely pathogenic']
        benign_keywords = ['benign', 'likely benign']
        
        classified = {}
        
        # Pathogenic variants
        classified['pathogenic'] = variants[
            variants['clinical_significance'].str.contains('|'.join(pathogenic_keywords))
        ]
        
        # Benign variants
        classified['benign'] = variants[
            variants['clinical_significance'].str.contains('|'.join(benign_keywords))
        ]
        
        # VUS (everything else)
        classified['vus'] = variants[
            ~variants['clinical_significance'].str.contains(
                '|'.join(pathogenic_keywords + benign_keywords)
            )
        ]
        
        self.logger.info(
            f"Classified variants: "
            f"{len(classified['pathogenic'])} pathogenic, "
            f"{len(classified['benign'])} benign, "
            f"{len(classified['vus'])} VUS"
        )
        
        return classified
    
    def map_variants_to_motifs(
        self,
        variants: pd.DataFrame,
        motif_positions: List[Tuple[int, int]],
        alignment_to_structure: Optional[Dict[int, int]] = None
    ) -> pd.DataFrame:
        """
        Determine which variants fall within motifs.
        
        Args:
            variants: DataFrame of variants
            motif_positions: List of (start, end) tuples for motifs
            alignment_to_structure: Optional mapping for position conversion
            
        Returns:
            DataFrame with added 'in_motif' column
        """
        variants = variants.copy()
        variants['in_motif'] = False
        variants['motif_id'] = None
        
        for idx, row in variants.iterrows():
            position = row['position']
            
            # Check if position is in any motif
            for motif_id, (start, end) in enumerate(motif_positions):
                if start <= position < end:
                    variants.at[idx, 'in_motif'] = True
                    variants.at[idx, 'motif_id'] = motif_id
                    break
        
        n_in_motif = variants['in_motif'].sum()
        self.logger.info(
            f"{n_in_motif}/{len(variants)} variants fall within motifs"
        )
        
        return variants
    
    def test_variant_enrichment(
        self,
        variants: pd.DataFrame,
        total_positions: int,
        motif_length: int
    ) -> Dict[str, float]:
        """
        Test if pathogenic variants are enriched in motifs.
        
        Uses Fisher's exact test.
        
        Args:
            variants: DataFrame with 'clinical_significance' and 'in_motif'
            total_positions: Total number of positions in protein
            motif_length: Total length of all motifs combined
            
        Returns:
            Dictionary of enrichment statistics
        """
        # Create contingency table
        # Rows: pathogenic vs non-pathogenic
        # Cols: in motif vs outside motif
        
        pathogenic = variants[
            variants['clinical_significance'].str.contains('pathogenic')
        ]
        non_pathogenic = variants[
            ~variants['clinical_significance'].str.contains('pathogenic')
        ]
        
        # Count variants
        path_in_motif = pathogenic['in_motif'].sum()
        path_outside = len(pathogenic) - path_in_motif
        nonpath_in_motif = non_pathogenic['in_motif'].sum()
        nonpath_outside = len(non_pathogenic) - nonpath_in_motif
        
        # Contingency table
        table = np.array([
            [path_in_motif, path_outside],
            [nonpath_in_motif, nonpath_outside]
        ])
        
        # Fisher's exact test
        odds_ratio, p_value = fisher_exact(table, alternative='greater')
        
        # Calculate variant density (variants per residue)
        motif_density = path_in_motif / motif_length if motif_length > 0 else 0
        outside_density = path_outside / (total_positions - motif_length)
        
        # Relative risk
        risk_motif = path_in_motif / (path_in_motif + nonpath_in_motif) if (path_in_motif + nonpath_in_motif) > 0 else 0
        risk_outside = path_outside / (path_outside + nonpath_outside) if (path_outside + nonpath_outside) > 0 else 0
        relative_risk = risk_motif / risk_outside if risk_outside > 0 else 0
        
        # 95% confidence interval for odds ratio
        # Using log transformation
        if odds_ratio > 0:
            log_or = np.log(odds_ratio)
            se_log_or = np.sqrt(
                1/path_in_motif + 1/path_outside + 
                1/nonpath_in_motif + 1/nonpath_outside
            ) if all([path_in_motif, path_outside, nonpath_in_motif, nonpath_outside]) else 0
            
            ci_lower = np.exp(log_or - 1.96 * se_log_or)
            ci_upper = np.exp(log_or + 1.96 * se_log_or)
        else:
            ci_lower = ci_upper = 0
        
        results = {
            'contingency_table': table.tolist(),
            'odds_ratio': float(odds_ratio),
            'odds_ratio_ci_lower': float(ci_lower),
            'odds_ratio_ci_upper': float(ci_upper),
            'p_value': float(p_value),
            'relative_risk': float(relative_risk),
            'pathogenic_in_motif': int(path_in_motif),
            'pathogenic_outside': int(path_outside),
            'variant_density_motif': float(motif_density),
            'variant_density_outside': float(outside_density),
            'enrichment_fold': float(motif_density / outside_density) if outside_density > 0 else 0
        }
        
        self.logger.info(
            f"Enrichment test: OR={odds_ratio:.2f}, "
            f"p={p_value:.4f}, "
            f"fold={results['enrichment_fold']:.2f}x"
        )
        
        return results
    
    def permutation_test_enrichment(
        self,
        variants: pd.DataFrame,
        motif_positions: List[Tuple[int, int]],
        total_positions: int,
        n_permutations: int = 10000
    ) -> Dict[str, float]:
        """
        Permutation test for variant enrichment.
        
        Randomly shuffle variant positions and recalculate enrichment.
        
        Args:
            variants: DataFrame of variants
            motif_positions: List of motif (start, end) tuples
            total_positions: Total protein length
            n_permutations: Number of permutations
            
        Returns:
            Dictionary with permutation test results
        """
        # Observed statistic: number of pathogenic variants in motifs
        pathogenic = variants[
            variants['clinical_significance'].str.contains('pathogenic')
        ]
        
        observed_in_motif = 0
        for _, var in pathogenic.iterrows():
            for start, end in motif_positions:
                if start <= var['position'] < end:
                    observed_in_motif += 1
                    break
        
        # Permutation test
        random_counts = []
        motif_length = sum(end - start for start, end in motif_positions)
        
        for _ in range(n_permutations):
            # Randomly place variants
            random_positions = np.random.choice(
                total_positions,
                size=len(pathogenic),
                replace=False
            )
            
            # Count how many fall in motifs
            count_in_motif = 0
            for pos in random_positions:
                for start, end in motif_positions:
                    if start <= pos < end:
                        count_in_motif += 1
                        break
            
            random_counts.append(count_in_motif)
        
        random_counts = np.array(random_counts)
        
        # Calculate p-value
        p_value = np.sum(random_counts >= observed_in_motif) / n_permutations
        
        results = {
            'observed': int(observed_in_motif),
            'expected': float(np.mean(random_counts)),
            'expected_std': float(np.std(random_counts)),
            'p_value': float(p_value),
            'z_score': float(
                (observed_in_motif - np.mean(random_counts)) / np.std(random_counts)
            ) if np.std(random_counts) > 0 else 0
        }
        
        self.logger.info(
            f"Permutation test: observed={observed_in_motif}, "
            f"expected={results['expected']:.1f}, "
            f"p={p_value:.4f}"
        )
        
        return results
    
    def analyze_motif_specific_variants(
        self,
        variants: pd.DataFrame,
        motif_positions: List[Tuple[int, int]]
    ) -> List[Dict]:
        """
        Analyze variant enrichment for each motif individually.
        
        Args:
            variants: DataFrame with variants
            motif_positions: List of motif positions
            
        Returns:
            List of per-motif statistics
        """
        results = []
        
        for motif_id, (start, end) in enumerate(motif_positions):
            # Find variants in this motif
            motif_variants = variants[
                (variants['position'] >= start) &
                (variants['position'] < end)
            ]
            
            pathogenic_count = motif_variants[
                motif_variants['clinical_significance'].str.contains('pathogenic')
            ].shape[0]
            
            benign_count = motif_variants[
                motif_variants['clinical_significance'].str.contains('benign')
            ].shape[0]
            
            vus_count = motif_variants[
                ~motif_variants['clinical_significance'].str.contains('pathogenic|benign')
            ].shape[0]
            
            motif_result = {
                'motif_id': motif_id,
                'start': start,
                'end': end,
                'length': end - start,
                'total_variants': len(motif_variants),
                'pathogenic': pathogenic_count,
                'benign': benign_count,
                'vus': vus_count,
                'pathogenic_density': pathogenic_count / (end - start)
            }
            
            results.append(motif_result)
        
        return results
    
    def save_variant_analysis(
        self,
        output_path: Path,
        enrichment_results: Dict,
        permutation_results: Dict,
        motif_specific_results: List[Dict],
        classified_variants: Dict[str, pd.DataFrame]
    ):
        """
        Save comprehensive variant analysis results.
        
        Args:
            output_path: Path to save JSON file
            enrichment_results: Fisher's test results
            permutation_results: Permutation test results
            motif_specific_results: Per-motif statistics
            classified_variants: Classified variant DataFrames
        """
        data = {
            'enrichment_test': enrichment_results,
            'permutation_test': permutation_results,
            'motif_specific': motif_specific_results,
            'variant_counts': {
                'pathogenic': len(classified_variants.get('pathogenic', [])),
                'benign': len(classified_variants.get('benign', [])),
                'vus': len(classified_variants.get('vus', []))
            }
        }
        
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with open(output_path, 'w') as f:
            json.dump(data, f, indent=2)
        
        self.logger.info(f"Variant analysis saved to {output_path}")
    
    def generate_variant_report(
        self,
        variants: pd.DataFrame,
        motif_positions: List[Tuple[int, int]],
        output_path: Path
    ):
        """
        Generate human-readable variant report.
        
        Args:
            variants: DataFrame of variants
            motif_positions: List of motif positions
            output_path: Path to save report
        """
        with open(output_path, 'w') as f:
            f.write("# Variant Analysis Report\n\n")
            
            # Summary statistics
            f.write("## Summary\n\n")
            f.write(f"Total variants: {len(variants)}\n")
            
            pathogenic = variants[
                variants['clinical_significance'].str.contains('pathogenic')
            ]
            f.write(f"Pathogenic variants: {len(pathogenic)}\n")
            
            in_motif = variants[variants['in_motif'] == True]
            f.write(f"Variants in motifs: {len(in_motif)}\n\n")
            
            # Per-motif breakdown
            f.write("## Motif-Specific Analysis\n\n")
            for motif_id, (start, end) in enumerate(motif_positions):
                motif_vars = variants[
                    (variants['position'] >= start) &
                    (variants['position'] < end)
                ]
                
                f.write(f"### Motif {motif_id + 1} (positions {start}-{end})\n\n")
                f.write(f"- Total variants: {len(motif_vars)}\n")
                
                path_vars = motif_vars[
                    motif_vars['clinical_significance'].str.contains('pathogenic')
                ]
                f.write(f"- Pathogenic: {len(path_vars)}\n")
                
                if len(path_vars) > 0:
                    f.write("\nPathogenic variants in this motif:\n")
                    for _, var in path_vars.iterrows():
                        f.write(f"  - Position {var['position']}: {var.get('variant', 'N/A')}\n")
                
                f.write("\n")
        
        self.logger.info(f"Variant report saved to {output_path}")
