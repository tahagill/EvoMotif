"""
Example: Correlate conservation scores with AlphaFold confidence.

This demonstrates how to validate that conserved motifs align with
high-confidence structural regions in AlphaFold models.
"""

from evomotif.structure import StructureMapper
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
from pathlib import Path


def analyze_conservation_confidence_correlation(
    structure_path: Path,
    conservation_scores: np.ndarray,
    alignment_to_structure: dict,
    chain_id: str = 'A'
):
    """
    Correlate conservation scores with AlphaFold pLDDT confidence.
    
    Args:
        structure_path: Path to AlphaFold PDB model
        conservation_scores: Array of conservation scores (0-1)
        alignment_to_structure: Mapping from alignment positions to structure residues
        chain_id: PDB chain identifier
    """
    # Load structure
    mapper = StructureMapper()
    structure = mapper.load_structure(structure_path)
    
    # Extract AlphaFold confidence scores
    confidence = mapper.get_alphafold_confidence(structure, chain_id=chain_id)
    
    # Map conservation to structure
    conservation_map = mapper.map_conservation_to_structure(
        conservation_scores,
        alignment_to_structure
    )
    
    # Find overlapping residues
    common_residues = set(confidence.keys()) & set(conservation_map.keys())
    
    if not common_residues:
        print("❌ No overlapping residues between conservation and structure")
        return
    
    # Extract paired values
    plddt_values = [confidence[res] for res in sorted(common_residues)]
    conservation_values = [conservation_map[res] for res in sorted(common_residues)]
    
    # Calculate correlation
    correlation, p_value = pearsonr(conservation_values, plddt_values)
    
    # Print results
    print("\n" + "="*60)
    print("Conservation-Confidence Correlation Analysis")
    print("="*60)
    print(f"Total residues analyzed: {len(common_residues)}")
    print(f"Pearson correlation: {correlation:.3f} (p={p_value:.2e})")
    print()
    
    if correlation > 0.5:
        print("✅ Strong positive correlation!")
        print("   → Conserved regions are structurally confident")
        print("   → AlphaFold model validates evolutionary conservation")
    elif correlation > 0.3:
        print("✓ Moderate positive correlation")
        print("   → Some alignment between conservation and structure")
    else:
        print("⚠️ Weak correlation")
        print("   → Conservation may be driven by non-structural factors")
        print("   → Check for functional sites in flexible regions")
    
    # Identify high-conservation, high-confidence residues (motif candidates)
    high_cons_conf = [
        res for res in common_residues
        if conservation_map[res] > 0.7 and confidence[res] > 70
    ]
    
    print(f"\nHigh conservation + high confidence residues: {len(high_cons_conf)}")
    if high_cons_conf:
        print(f"  Residues: {sorted(high_cons_conf)[:10]}" + 
              (" ..." if len(high_cons_conf) > 10 else ""))
        print("  → Strong candidates for functional motifs")
    
    # Identify high-conservation, low-confidence residues (interesting cases)
    high_cons_low_conf = [
        res for res in common_residues
        if conservation_map[res] > 0.7 and confidence[res] < 50
    ]
    
    print(f"\nHigh conservation + low confidence residues: {len(high_cons_low_conf)}")
    if high_cons_low_conf:
        print(f"  Residues: {sorted(high_cons_low_conf)}")
        print("  → May be functional sites in disordered regions")
        print("  → Consider experimental validation")
    
    # Create visualization
    plt.figure(figsize=(10, 6))
    
    plt.subplot(1, 2, 1)
    plt.scatter(conservation_values, plddt_values, alpha=0.6, s=30)
    plt.xlabel('Conservation Score', fontsize=12)
    plt.ylabel('AlphaFold pLDDT', fontsize=12)
    plt.title(f'Correlation: r={correlation:.3f}', fontsize=14)
    plt.axhline(70, color='red', linestyle='--', alpha=0.3, label='pLDDT threshold')
    plt.axvline(0.7, color='blue', linestyle='--', alpha=0.3, label='Conservation threshold')
    plt.legend()
    plt.grid(alpha=0.3)
    
    plt.subplot(1, 2, 2)
    residue_numbers = sorted(common_residues)
    plt.plot(residue_numbers, 
             [conservation_map[r] for r in residue_numbers],
             label='Conservation', alpha=0.7, linewidth=2)
    plt.plot(residue_numbers,
             [confidence[r]/100 for r in residue_numbers],  # Normalize pLDDT to 0-1
             label='pLDDT (normalized)', alpha=0.7, linewidth=2)
    plt.xlabel('Residue Number', fontsize=12)
    plt.ylabel('Score', fontsize=12)
    plt.title('Conservation vs Confidence Profile', fontsize=14)
    plt.legend()
    plt.grid(alpha=0.3)
    
    plt.tight_layout()
    output_file = 'alphafold_conservation_correlation.png'
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    print(f"\n📊 Plot saved to: {output_file}")
    
    return {
        'correlation': correlation,
        'p_value': p_value,
        'n_residues': len(common_residues),
        'high_conf_high_cons': high_cons_conf,
        'low_conf_high_cons': high_cons_low_conf
    }


# Example usage
if __name__ == '__main__':
    print("""
Example Usage:
--------------

from example_alphafold import analyze_conservation_confidence_correlation
from pathlib import Path
import numpy as np

# Your data
structure_path = Path('AF-P04637-F1-model_v4.pdb')  # AlphaFold model
conservation = np.array([0.85, 0.92, 0.78, ...])    # From your analysis
mapping = {0: 1, 1: 2, 2: 3, ...}                    # Alignment to structure

# Run analysis
results = analyze_conservation_confidence_correlation(
    structure_path,
    conservation,
    mapping,
    chain_id='A'
)

# Interpretation:
# - High correlation (>0.5): Conservation driven by structural constraint
# - Low correlation (<0.3): Conservation may reflect functional sites
# - High cons + high conf: Core functional motifs (targets for drugs)
# - High cons + low conf: Functional sites in flexible regions (allosteric sites)
    """)
