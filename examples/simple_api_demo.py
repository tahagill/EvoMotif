"""
EvoMotif Simple API Example

This demonstrates the new high-level API that makes EvoMotif
incredibly easy to use - just 2 lines of code!
"""

import evomotif

# ============================================================================
# Example 1: Basic Analysis (Simplest Possible)
# ============================================================================

print("Example 1: Basic Analysis")
print("-" * 60)

# Just two lines - that's it!
results = evomotif.analyze_protein(
    protein_name="ubiquitin",
    email="your.email@example.com"  # Replace with your email
)

# Print summary
print(results.summary())

# Access results
print(f"\nNumber of motifs found: {len(results.motifs)}")
print(f"Number of conserved positions: {len(results.conserved_positions)}")
print(f"Mean conservation score: {results.data['mean_conservation']:.3f}")


# ============================================================================
# Example 2: With 3D Structure
# ============================================================================

print("\n\n" + "="*70)
print("Example 2: Analysis with 3D Structure Mapping")
print("-" * 60)

results = evomotif.analyze_protein(
    protein_name="p53",
    email="your.email@example.com",
    pdb_id="1TUP",  # Add structure visualization
    max_sequences=30  # Use fewer sequences for faster demo
)

print(results.summary())


# ============================================================================
# Example 3: Custom Parameters
# ============================================================================

print("\n\n" + "="*70)
print("Example 3: Custom Analysis Parameters")
print("-" * 60)

results = evomotif.analyze_protein(
    protein_name="insulin",
    email="your.email@example.com",
    output_dir="./my_results/insulin",  # Custom output location
    max_sequences=100,                   # More sequences
    min_conservation=0.65,               # Lower threshold
    threads=8,                           # Use more CPU cores
    verbose=True                         # Show detailed logging
)

print(results.summary())


# ============================================================================
# Example 4: Using the Pipeline Object (More Control)
# ============================================================================

print("\n\n" + "="*70)
print("Example 4: Using Pipeline Object")
print("-" * 60)

# Create pipeline object for more control
pipeline = evomotif.EvoMotifPipeline()

# Check dependencies before running
deps = pipeline.check_dependencies()
print("Dependency status:")
for tool, available in deps.items():
    status = "✓" if available else "✗"
    print(f"  {status} {tool}")

# Run analysis
if all(deps.values()):
    results = pipeline.analyze(
        protein_name="AKT1",
        email="your.email@example.com",
        pdb_id="3CQU"
    )
    
    print(results.summary())
    
    # Export results to JSON
    results.export_json(results.output_dir / "custom_export.json")
    print(f"\n📄 Results exported to: {results.output_dir / 'custom_export.json'}")
else:
    print("\n⚠️  Some dependencies are missing. Please install them first.")


# ============================================================================
# Example 5: Accessing Individual Results
# ============================================================================

print("\n\n" + "="*70)
print("Example 5: Accessing Individual Results")
print("-" * 60)

# Run analysis
results = evomotif.analyze_protein(
    protein_name="ubiquitin",
    email="your.email@example.com",
    max_sequences=30
)

# Access specific data
print(f"Protein analyzed: {results.protein}")
print(f"Number of sequences: {results.n_sequences}")

# Get top conserved positions
if results.conserved_positions:
    print("\nTop 10 conserved positions:")
    sorted_positions = sorted(
        results.conserved_positions,
        key=lambda x: x['conservation'],
        reverse=True
    )
    
    for i, pos in enumerate(sorted_positions[:10], 1):
        print(f"  {i:2d}. Position {pos['position']:3d}: "
              f"{pos['amino_acid']} (conservation={pos['conservation']:.3f})")

# Get file paths
print("\nGenerated files:")
for file_type, filepath in results.files.items():
    print(f"  {file_type}: {filepath}")


# ============================================================================
# Example 6: Error Handling
# ============================================================================

print("\n\n" + "="*70)
print("Example 6: Error Handling")
print("-" * 60)

try:
    results = evomotif.analyze_protein(
        protein_name="nonexistent_protein_xyz123",
        email="your.email@example.com",
        max_sequences=10
    )
except RuntimeError as e:
    print(f"❌ Analysis failed: {e}")
    print("This is expected for proteins that don't exist in NCBI")


# ============================================================================
# Comparison: Old vs New API
# ============================================================================

print("\n\n" + "="*70)
print("API Comparison")
print("="*70)

print("""
OLD API (Still works for power users):
---------------------------------------
from evomotif.retrieval import SequenceRetriever
from evomotif.alignment import SequenceAligner
from evomotif.conservation import ConservationScorer
from evomotif.motif_discovery import MotifDiscoverer
from evomotif.stats import StatisticalAnalyzer
# ... 100+ more lines of boilerplate ...


NEW API (Recommended for most users):
--------------------------------------
import evomotif

results = evomotif.analyze_protein("p53", "user@email.com")
print(results.summary())

# Done! 🎉


Both APIs work! Use the simple API for quick analyses,
or use individual modules for maximum control.
""")


print("\n" + "="*70)
print("✅ All examples complete!")
print("="*70)
