"""
Sequence retrieval module using NCBI Entrez API.

Handles fetching protein sequences from multiple species,
filtering, and redundancy removal.
"""

from typing import List, Dict, Optional, Set
from Bio import Entrez, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import logging
import time
from pathlib import Path

logger = logging.getLogger(__name__)


class SequenceRetriever:
    """Retrieve and process protein sequences from NCBI."""
    
    def __init__(self, email: str, api_key: Optional[str] = None):
        """
        Initialize the retriever.
        
        Args:
            email: Email for NCBI Entrez (required)
            api_key: Optional API key for higher rate limits
        """
        Entrez.email = email
        if api_key:
            Entrez.api_key = api_key
        self.logger = logging.getLogger(__name__)
    
    def search_protein(
        self,
        protein_name: str,
        organisms: Optional[List[str]] = None,
        max_sequences: int = 1000,
        reviewed_only: bool = True
    ) -> List[str]:
        """
        Search for protein sequences in NCBI.
        
        Args:
            protein_name: Name of the protein to search
            organisms: List of organism names to include (None = all)
            max_sequences: Maximum number of sequences to retrieve
            reviewed_only: Only retrieve reviewed/SwissProt sequences
            
        Returns:
            List of NCBI protein IDs
        """
        # Build query
        query_parts = [f'"{protein_name}"[Protein Name]']
        
        if reviewed_only:
            query_parts.append('"reviewed"[Filter]')
        
        if organisms:
            org_query = " OR ".join([f'"{org}"[Organism]' for org in organisms])
            query_parts.append(f"({org_query})")
        
        query = " AND ".join(query_parts)
        
        self.logger.info(f"Searching with query: {query}")
        
        try:
            handle = Entrez.esearch(
                db="protein",
                term=query,
                retmax=max_sequences,
                sort="relevance"
            )
            record = Entrez.read(handle)
            handle.close()
            
            ids = record["IdList"]
            self.logger.info(f"Found {len(ids)} sequences")
            return ids
            
        except Exception as e:
            self.logger.error(f"Search failed: {e}")
            raise
    
    def fetch_sequences(
        self,
        id_list: List[str],
        batch_size: int = 100
    ) -> List[SeqRecord]:
        """
        Fetch sequences from NCBI in batches.
        
        Args:
            id_list: List of NCBI protein IDs
            batch_size: Number of sequences to fetch per request
            
        Returns:
            List of SeqRecord objects
        """
        sequences = []
        
        for i in range(0, len(id_list), batch_size):
            batch_ids = id_list[i:i + batch_size]
            self.logger.info(f"Fetching batch {i//batch_size + 1}: {len(batch_ids)} sequences")
            
            try:
                handle = Entrez.efetch(
                    db="protein",
                    id=batch_ids,
                    rettype="fasta",
                    retmode="text"
                )
                batch_records = list(SeqIO.parse(handle, "fasta"))
                handle.close()
                
                sequences.extend(batch_records)
                
                # Be nice to NCBI servers
                time.sleep(0.34)  # ~3 requests/second
                
            except Exception as e:
                self.logger.warning(f"Failed to fetch batch: {e}")
                continue
        
        self.logger.info(f"Successfully fetched {len(sequences)} sequences")
        return sequences
    
    def filter_sequences(
        self,
        sequences: List[SeqRecord],
        min_length: int = 50,
        max_length: Optional[int] = None,
        remove_fragments: bool = False
    ) -> List[SeqRecord]:
        """
        Filter sequences based on length and quality criteria.
        
        Args:
            sequences: List of SeqRecord objects
            min_length: Minimum sequence length
            max_length: Maximum sequence length (None = no limit)
            remove_fragments: Remove sequences marked as fragments
            
        Returns:
            Filtered list of SeqRecord objects
        """
        filtered = []
        
        for seq in sequences:
            # Length filters
            if len(seq.seq) < min_length:
                continue
            if max_length and len(seq.seq) > max_length:
                continue
            
            # Fragment filter
            if remove_fragments:
                desc_lower = seq.description.lower()
                if any(word in desc_lower for word in ['fragment', 'partial']):
                    continue
            
            filtered.append(seq)
        
        self.logger.info(f"Filtered {len(sequences)} -> {len(filtered)} sequences")
        return filtered
    
    def remove_duplicates(
        self,
        sequences: List[SeqRecord],
        identity_threshold: float = 0.95
    ) -> List[SeqRecord]:
        """
        Remove duplicate sequences based on sequence identity.
        
        Currently uses exact match. The identity_threshold parameter
        is reserved for future implementation of clustering (e.g., CD-HIT).
        For now, only 100% identical sequences are removed.
        
        Args:
            sequences: List of SeqRecord objects
            identity_threshold: Similarity threshold (currently unused; exact match only)
            
        Returns:
            Deduplicated list of SeqRecord objects
        """
        seen_seqs: Set[str] = set()
        unique = []
        
        for seq in sequences:
            seq_str = str(seq.seq).upper()
            if seq_str not in seen_seqs:
                seen_seqs.add(seq_str)
                unique.append(seq)
        
        self.logger.info(f"Removed {len(sequences) - len(unique)} duplicates")
        return unique
    
    def save_sequences(self, sequences: List[SeqRecord], output_path: Path):
        """
        Save sequences to FASTA file.
        
        Args:
            sequences: List of SeqRecord objects
            output_path: Path to output FASTA file
        """
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with open(output_path, 'w') as f:
            SeqIO.write(sequences, f, "fasta")
        self.logger.info(f"Saved {len(sequences)} sequences to {output_path}")
    
    def retrieve_all(
        self,
        protein_name: str,
        output_path: Path,
        organisms: Optional[List[str]] = None,
        max_sequences: int = 1000,
        reviewed_only: bool = False,
        **kwargs
    ) -> List[SeqRecord]:
        """
        Complete retrieval pipeline.
        
        Args:
            protein_name: Name of protein to search
            output_path: Path to save sequences
            organisms: List of organisms (None = all)
            max_sequences: Maximum sequences to retrieve
            reviewed_only: Only search reviewed/SwissProt entries
            **kwargs: Additional arguments for filtering
            
        Returns:
            List of processed SeqRecord objects
        """
        # Search
        ids = self.search_protein(protein_name, organisms, max_sequences, reviewed_only)
        
        if not ids:
            raise ValueError(f"No sequences found for {protein_name}")
        
        # Fetch
        sequences = self.fetch_sequences(ids)
        
        # Filter
        sequences = self.filter_sequences(sequences, **kwargs)
        
        # Remove duplicates
        sequences = self.remove_duplicates(sequences)
        
        # Save
        self.save_sequences(sequences, output_path)
        
        return sequences
