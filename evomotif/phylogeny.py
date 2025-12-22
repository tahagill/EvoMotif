"""
Phylogenetic analysis module.

Builds phylogenetic trees and analyzes motif evolutionary origin
using maximum likelihood and ancestral state reconstruction.
"""

import numpy as np
import logging
from typing import Dict, List, Optional, Tuple
from pathlib import Path
from Bio import Phylo
from Bio.Align import MultipleSeqAlignment
from ete3 import Tree
import subprocess

# Optional imports for visualization (requires PyQt)
try:
    from ete3.treeview import TreeStyle, NodeStyle
    VISUALIZATION_AVAILABLE = True
except ImportError:
    VISUALIZATION_AVAILABLE = False

logger = logging.getLogger(__name__)


class PhylogeneticAnalyzer:
    """Perform phylogenetic analysis and motif evolution inference."""
    
    def __init__(self, fasttree_path: str = "fasttree"):
        """Initialize phylogenetic analyzer.
        
        Args:
            fasttree_path: Path to FastTree executable (default: "fasttree")
        """
        self.logger = logging.getLogger(__name__)
        self.fasttree_path = fasttree_path
    
    def build_tree_fasttree(
        self,
        alignment_file: Path,
        output_tree: Path,
        model: str = "JTT"
    ) -> Tree:
        """
        Build phylogenetic tree using FastTree.
        
        Args:
            alignment_file: Path to alignment file
            output_tree: Path to save Newick tree
            model: Evolutionary model (JTT, WAG, LG)
            
        Returns:
            ETE3 Tree object
        """
        cmd = [self.fasttree_path]
        
        # Model selection
        if model.upper() == "WAG":
            cmd.append("-wag")
        elif model.upper() == "LG":
            cmd.append("-lg")
        # JTT is default
        
        cmd.append(str(alignment_file))
        
        self.logger.info(f"Building tree with FastTree (model: {model})")
        
        try:
            with open(output_tree, 'w') as outfile:
                result = subprocess.run(
                    cmd,
                    stdout=outfile,
                    stderr=subprocess.PIPE,
                    text=True,
                    check=True
                )
            
            self.logger.info(f"Tree saved to {output_tree}")
            
            # Load tree with ETE3
            tree = Tree(str(output_tree))
            
            return tree
            
        except FileNotFoundError:
            self.logger.error(
                f"FastTree not found at '{self.fasttree_path}'. "
                "Install: http://www.microbesonline.org/fasttree/ "
                "or specify correct path in __init__(fasttree_path='...')"
            )
            raise
        except subprocess.CalledProcessError as e:
            self.logger.error(f"FastTree failed: {e.stderr}")
            raise
    
    def build_tree_iqtree(
        self,
        alignment_file: Path,
        output_prefix: Path,
        model: str = "JTT",
        bootstrap: int = 1000,
        threads: int = 4
    ) -> Tree:
        """
        Build phylogenetic tree using IQ-TREE with bootstrap support.
        
        Args:
            alignment_file: Path to alignment file
            output_prefix: Prefix for output files
            model: Evolutionary model
            bootstrap: Number of bootstrap replicates
            threads: Number of CPU threads
            
        Returns:
            ETE3 Tree object
        """
        cmd = [
            "iqtree",
            "-s", str(alignment_file),
            "-m", model,
            "-bb", str(bootstrap),
            "-nt", str(threads),
            "-pre", str(output_prefix)
        ]
        
        self.logger.info(
            f"Building tree with IQ-TREE "
            f"(model: {model}, bootstrap: {bootstrap})"
        )
        
        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                check=True
            )
            
            tree_file = Path(str(output_prefix) + ".treefile")
            self.logger.info(f"Tree saved to {tree_file}")
            
            # Load tree
            tree = Tree(str(tree_file))
            
            return tree
            
        except FileNotFoundError:
            self.logger.error("IQ-TREE not found. Install: http://www.iqtree.org/")
            raise
        except subprocess.CalledProcessError as e:
            self.logger.error(f"IQ-TREE failed: {e.stderr}")
            raise
    
    def annotate_tree_with_motifs(
        self,
        tree: Tree,
        motif_presence: Dict[str, Dict[int, bool]]
    ) -> Tree:
        """
        Annotate tree leaves with motif presence/absence.
        
        Args:
            tree: ETE3 Tree object
            motif_presence: Dict mapping {seq_id: {motif_id: present}}
            
        Returns:
            Annotated Tree object
        """
        for leaf in tree.iter_leaves():
            leaf_name = leaf.name
            
            if leaf_name in motif_presence:
                # Add motif features
                for motif_id, present in motif_presence[leaf_name].items():
                    leaf.add_feature(f"motif_{motif_id}", int(present))
        
        self.logger.info("Tree annotated with motif presence/absence")
        return tree
    
    def reconstruct_ancestral_states(
        self,
        tree: Tree,
        motif_id: int
    ) -> Tree:
        """
        Reconstruct ancestral states for motif presence using parsimony.
        
        Args:
            tree: Annotated ETE3 Tree object
            motif_id: Motif ID to analyze
            
        Returns:
            Tree with ancestral states
        """
        feature_name = f"motif_{motif_id}"
        
        # Fitch parsimony algorithm (simplified)
        def fitch_up(node):
            """Upward pass: assign possible states to internal nodes."""
            if node.is_leaf():
                # Leaf nodes have known states
                state = getattr(node, feature_name, 0)
                node.add_feature(f"{feature_name}_states", {state})
            else:
                # Internal nodes: intersection of children states
                child_states = [
                    getattr(child, f"{feature_name}_states", {0})
                    for child in node.children
                ]
                
                intersection = set.intersection(*child_states)
                if intersection:
                    node.add_feature(f"{feature_name}_states", intersection)
                else:
                    union = set.union(*child_states)
                    node.add_feature(f"{feature_name}_states", union)
        
        def fitch_down(node, parent_state=None):
            """Downward pass: resolve ambiguous states."""
            states = getattr(node, f"{feature_name}_states", {0})
            
            if parent_state is not None and parent_state in states:
                chosen_state = parent_state
            else:
                # Choose most parsimonious state (prefer 0 for absence)
                chosen_state = min(states)
            
            node.add_feature(feature_name, chosen_state)
            
            # Recurse to children
            for child in node.children:
                fitch_down(child, chosen_state)
        
        # Run Fitch algorithm - correct ETE3 API usage
        for node in tree.traverse("postorder"):
            fitch_up(node)
        fitch_down(tree)
        
        self.logger.info(f"Reconstructed ancestral states for motif {motif_id}")
        return tree
    
    def identify_motif_origin(
        self,
        tree: Tree,
        motif_id: int
    ) -> Optional[str]:
        """
        Identify the node where motif first appeared.
        
        Args:
            tree: Tree with reconstructed ancestral states
            motif_id: Motif ID to analyze
            
        Returns:
            Name/description of the node where motif originated
        """
        feature_name = f"motif_{motif_id}"
        origin_node = None
        
        # Find the most recent common ancestor of all motif-bearing leaves
        leaves_with_motif = [
            leaf for leaf in tree.iter_leaves()
            if getattr(leaf, feature_name, 0) == 1
        ]
        
        if not leaves_with_motif:
            self.logger.warning(f"No leaves have motif {motif_id}")
            return None
        
        if len(leaves_with_motif) == 1:
            origin_node = leaves_with_motif[0]
        else:
            origin_node = tree.get_common_ancestor(leaves_with_motif)
        
        # Get taxonomic info if available
        origin_info = self._get_node_taxonomy(origin_node, tree)
        
        self.logger.info(f"Motif {motif_id} originated at: {origin_info}")
        return origin_info
    
    def _get_node_taxonomy(self, node: Tree, tree: Tree) -> str:
        """
        Get taxonomic description of a node.
        
        Args:
            node: Tree node
            tree: Full tree
            
        Returns:
            Taxonomic description string
        """
        if node.is_leaf():
            return node.name
        
        # Get all leaf descendants
        descendants = node.get_leaf_names()
        
        # Simple heuristic: find common taxonomic groups
        # In real implementation, would query taxonomy database
        return f"MRCA of {len(descendants)} sequences"
    
    def calculate_motif_age(
        self,
        tree: Tree,
        motif_id: int
    ) -> float:
        """
        Estimate motif age based on branch lengths.
        
        Args:
            tree: Tree with reconstructed states
            motif_id: Motif ID
            
        Returns:
            Evolutionary distance (branch length sum) from origin to root
        """
        feature_name = f"motif_{motif_id}"
        
        # Find origin node
        leaves_with_motif = [
            leaf for leaf in tree.iter_leaves()
            if getattr(leaf, feature_name, 0) == 1
        ]
        
        if not leaves_with_motif:
            return 0.0
        
        if len(leaves_with_motif) == 1:
            origin_node = leaves_with_motif[0]
        else:
            origin_node = tree.get_common_ancestor(leaves_with_motif)
        
        # Calculate distance from origin node to root
        root = tree.get_tree_root()
        age = tree.get_distance(origin_node, root)
        
        return age
    
    def visualize_tree(
        self,
        tree: Tree,
        output_path: Path,
        motif_id: Optional[int] = None,
        show_bootstrap: bool = True
    ):
        """
        Visualize phylogenetic tree with motif annotations.
        
        Requires PyQt4 or PyQt5 for GUI rendering.
        
        Args:
            tree: ETE3 Tree object
            output_path: Path to save visualization
            motif_id: Optional motif ID to highlight
            show_bootstrap: Show bootstrap support values
        """
        if not VISUALIZATION_AVAILABLE:
            self.logger.error(
                "Tree visualization requires PyQt4 or PyQt5. "
                "Install with: pip install PyQt5"
            )
            raise ImportError("TreeStyle and NodeStyle not available")
        
        # Create tree style
        ts = TreeStyle()
        ts.show_leaf_name = True
        ts.show_branch_length = False
        ts.show_branch_support = show_bootstrap
        
        # Style nodes based on motif presence
        if motif_id is not None:
            feature_name = f"motif_{motif_id}"
            
            for node in tree.traverse():
                nstyle = NodeStyle()
                
                if hasattr(node, feature_name):
                    if getattr(node, feature_name) == 1:
                        nstyle["bgcolor"] = "#ff6b6b"  # Red for motif present
                        nstyle["size"] = 10
                    else:
                        nstyle["bgcolor"] = "#e0e0e0"  # Gray for absent
                        nstyle["size"] = 5
                
                node.set_style(nstyle)
        
        # Render tree
        tree.render(str(output_path), tree_style=ts)
        self.logger.info(f"Tree visualization saved to {output_path}")
    
    def analyze_motif_coevolution(
        self,
        tree: Tree,
        motif_ids: List[int]
    ) -> np.ndarray:
        """
        Analyze co-evolution patterns between motifs.
        
        Args:
            tree: Annotated tree
            motif_ids: List of motif IDs to analyze
            
        Returns:
            Correlation matrix of motif co-occurrence
        """
        # Extract motif presence matrix
        n_motifs = len(motif_ids)
        n_leaves = len(list(tree.iter_leaves()))
        
        presence_matrix = np.zeros((n_leaves, n_motifs))
        
        for i, leaf in enumerate(tree.iter_leaves()):
            for j, motif_id in enumerate(motif_ids):
                feature_name = f"motif_{motif_id}"
                presence_matrix[i, j] = getattr(leaf, feature_name, 0)
        
        # Calculate correlation
        correlation = np.corrcoef(presence_matrix.T)
        
        self.logger.info("Calculated motif co-evolution correlation matrix")
        return correlation
