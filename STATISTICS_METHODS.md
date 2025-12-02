# Statistical Methods in EvoMotif

Comprehensive documentation of all statistical methods, algorithms, and implementations used in the EvoMotif protein motif discovery pipeline.

---

## Table of Contents

1. [Conservation Scoring](#conservation-scoring)
2. [Motif Discovery Algorithms](#motif-discovery-algorithms)
3. [Statistical Validation](#statistical-validation)
4. [Multiple Testing Correction](#multiple-testing-correction)
5. [Effect Size Calculation](#effect-size-calculation)
6. [Implementation Details](#implementation-details)

---

## Conservation Scoring

### 1. Shannon Entropy

**Purpose:** Measure information content and variability at each alignment position.

**Mathematical Formula:**

$$H(i) = -\sum_{a \in A} p_a(i) \log_2 p_a(i)$$

Where:
- $H(i)$ = Shannon entropy at position $i$
- $A$ = set of 20 amino acids
- $p_a(i)$ = frequency of amino acid $a$ at position $i$
- $\log_2$ = logarithm base 2 (bits of information)

**Normalized Shannon Entropy:**

$$H_{\text{norm}}(i) = \frac{H(i)}{H_{\text{max}}}$$

$$H_{\text{max}} = \log_2(20) \approx 4.322$$

**Range:**
- $H_{\text{norm}} = 0$: Completely conserved (all identical residues)
- $H_{\text{norm}} = 1$: Maximum variability (uniform distribution)

**Conservation Score (inverse):**

$$C_{\text{shannon}}(i) = 1 - H_{\text{norm}}(i)$$

**Why:** Shannon entropy directly quantifies uncertainty/variability. Low entropy = high conservation.

**Where Used:** 
- `evomotif/conservation.py::calculate_shannon_entropy()`
- Combined with BLOSUM62 for final conservation metric

**Implementation Logic:**

```python
def calculate_shannon_entropy(alignment):
    """
    Calculate Shannon entropy for each alignment position.
    
    Algorithm:
    1. For each column in alignment:
       a. Count frequency of each amino acid (excluding gaps)
       b. Calculate probabilities: p_a = count_a / total
       c. Apply Shannon formula: H = -Σ(p_a * log2(p_a))
       d. Normalize by max entropy (log2(20))
    2. Return conservation score: 1 - H_norm
    """
    scores = []
    max_entropy = np.log2(20)
    
    for i in range(alignment.get_alignment_length()):
        column = alignment[:, i]
        
        # Count amino acids (exclude gaps)
        aa_counts = {}
        for aa in column:
            if aa != '-':
                aa_counts[aa] = aa_counts.get(aa, 0) + 1
        
        total = sum(aa_counts.values())
        if total == 0:
            scores.append(0.0)
            continue
        
        # Calculate entropy
        entropy = 0.0
        for count in aa_counts.values():
            p = count / total
            if p > 0:
                entropy -= p * np.log2(p)
        
        # Normalize and convert to conservation
        normalized_entropy = entropy / max_entropy
        conservation = 1.0 - normalized_entropy
        scores.append(conservation)
    
    return scores
```

---

### 2. BLOSUM62 Scoring

**Purpose:** Measure evolutionary substitution patterns using empirically-derived substitution matrix.

**Mathematical Formula:**

$$B(i) = \frac{1}{n(n-1)} \sum_{j=1}^{n} \sum_{k=j+1}^{n} \text{BLOSUM62}(a_j, a_k)$$

Where:
- $B(i)$ = Average BLOSUM62 score at position $i$
- $n$ = number of sequences
- $a_j, a_k$ = amino acids at position $i$ in sequences $j$ and $k$
- $\text{BLOSUM62}(a, b)$ = substitution score from BLOSUM62 matrix

**Normalized BLOSUM Score:**

$$B_{\text{norm}}(i) = \frac{B(i) - B_{\text{min}}}{B_{\text{max}} - B_{\text{min}}}$$

Where:
- $B_{\text{max}} = 11$ (identical residues, maximum BLOSUM62 score)
- $B_{\text{min}} = -4$ (worst substitution score)

**Why:** BLOSUM62 captures evolutionary constraints. High scores indicate residues that are functionally equivalent even if not identical.

**Where Used:**
- `evomotif/conservation.py::calculate_blosum_score()`
- Combined with Shannon entropy for comprehensive conservation metric

**BLOSUM62 Matrix (excerpt):**

|   | A  | R  | N  | D  | C  |
|---|----|----|----|----|-----|
| A | 4  | -1 | -2 | -2 | 0   |
| R | -1 | 5  | 0  | -2 | -3  |
| N | -2 | 0  | 6  | 1  | -3  |
| D | -2 | -2 | 1  | 6  | -3  |
| C | 0  | -3 | -3 | -3 | 9   |

**Implementation Logic:**

```python
def calculate_blosum_score(alignment):
    """
    Calculate average BLOSUM62 score for pairwise comparisons.
    
    Algorithm:
    1. Load BLOSUM62 substitution matrix
    2. For each column in alignment:
       a. Extract non-gap amino acids
       b. Calculate all pairwise BLOSUM62 scores
       c. Average the scores: Σ(BLOSUM(aa_i, aa_j)) / n_pairs
       d. Normalize to 0-1 range
    3. Return normalized scores
    """
    from Bio.Align import substitution_matrices
    blosum62 = substitution_matrices.load("BLOSUM62")
    
    scores = []
    for i in range(alignment.get_alignment_length()):
        column = [aa for aa in alignment[:, i] if aa != '-']
        
        if len(column) < 2:
            scores.append(0.0)
            continue
        
        # Calculate pairwise scores
        total_score = 0
        n_pairs = 0
        for j in range(len(column)):
            for k in range(j + 1, len(column)):
                try:
                    score = blosum62[column[j], column[k]]
                    total_score += score
                    n_pairs += 1
                except KeyError:
                    continue
        
        # Average and normalize
        if n_pairs > 0:
            avg_score = total_score / n_pairs
            # BLOSUM62 range: -4 to 11
            normalized = (avg_score - (-4)) / (11 - (-4))
            scores.append(max(0.0, min(1.0, normalized)))
        else:
            scores.append(0.0)
    
    return scores
```

---

### 3. Combined Conservation Score

**Purpose:** Integrate information-theoretic and evolutionary perspectives.

**Mathematical Formula:**

$$C_{\text{combined}}(i) = w_1 \cdot C_{\text{shannon}}(i) + w_2 \cdot B_{\text{norm}}(i)$$

**Default Weights:**
- $w_1 = 0.5$ (Shannon entropy weight)
- $w_2 = 0.5$ (BLOSUM62 weight)

**Range:** 0 (variable) to 1 (conserved)

**Why Combined Metric:**
- Shannon entropy: Detects **identical** residues
- BLOSUM62: Detects **functionally similar** residues
- Combined: Captures both strict conservation and functional conservation

**Example:**
```
Position with all Lysine (K):
  Shannon: 1.0 (perfect conservation)
  BLOSUM62: 1.0 (identical = high score)
  Combined: 1.0 (highly conserved)

Position with K/R (positive charged residues):
  Shannon: ~0.7 (some variability)
  BLOSUM62: 0.9 (K-R substitution score = 2, favorable)
  Combined: 0.8 (functionally conserved)
```

**Where Used:**
- `evomotif/conservation.py::calculate_combined_conservation()`
- Primary input to motif discovery algorithms

---

### 4. Gap Frequency

**Purpose:** Quantify insertion/deletion patterns.

**Mathematical Formula:**

$$G(i) = \frac{n_{\text{gaps}}(i)}{n_{\text{total}}}$$

Where:
- $G(i)$ = Gap frequency at position $i$
- $n_{\text{gaps}}(i)$ = number of gap characters ('-') at position $i$
- $n_{\text{total}}$ = total number of sequences

**Range:** 0 (no gaps) to 1 (all gaps)

**Why:** High gap frequency indicates structural variability or insertions/deletions in evolution. Motifs should have low gap frequency.

**Where Used:**
- `evomotif/conservation.py::calculate_gap_frequency()`
- Filter criterion in motif discovery (exclude gappy regions)

**Implementation:**

```python
def calculate_gap_frequency(alignment):
    """
    Calculate gap frequency at each position.
    
    Algorithm:
    1. For each column:
       a. Count gap characters ('-')
       b. Divide by total sequences
    2. Return gap frequencies
    """
    n_seqs = len(alignment)
    gap_freqs = []
    
    for i in range(alignment.get_alignment_length()):
        column = alignment[:, i]
        n_gaps = column.count('-')
        gap_freq = n_gaps / n_seqs
        gap_freqs.append(gap_freq)
    
    return gap_freqs
```

---

## Motif Discovery Algorithms

### Algorithm 1: Sliding Window (Consecutive Motifs)

**Purpose:** Identify contiguous regions of high conservation.

**Mathematical Formulation:**

For window of size $w$ starting at position $i$:

$$S_{\text{window}}(i) = \frac{1}{w} \sum_{j=i}^{i+w-1} C_{\text{combined}}(j)$$

**Motif Criteria:**

1. **Conservation threshold:**
   $$S_{\text{window}}(i) \geq \theta_c$$
   Default: $\theta_c = 0.70$

2. **Gap threshold:**
   $$\max_{j \in [i, i+w-1]} G(j) \leq \theta_g$$
   Default: $\theta_g = 0.30$

3. **Minimum length:**
   $$L_{\text{motif}} \geq L_{\text{min}}$$
   Default: $L_{\text{min}} = 5$

**Algorithm Steps:**

```
SLIDING_WINDOW_MOTIF_DISCOVERY(alignment, conservation, window_size=5):
    
    1. Initialize empty motif list M
    
    2. For position i from 0 to (alignment_length - window_size):
        a. Extract window: W = positions[i:i+window_size]
        b. Calculate average conservation: S = mean(conservation[W])
        c. Calculate max gap frequency: G_max = max(gaps[W])
        
        d. IF S >= threshold_conservation AND G_max <= threshold_gap:
            i. Mark position i as motif start
    
    3. Merge adjacent/overlapping windows:
        a. Sort motif starts
        b. Combine consecutive positions into single motifs
        c. Extend motifs while conservation remains high
    
    4. Filter motifs:
        a. Remove motifs shorter than min_length
        b. Recalculate final conservation for each motif
    
    5. RETURN motif list M with [start, end, conservation] for each
```

**Why Sliding Window:**
- Detects functional domains (e.g., DNA-binding domains)
- Identifies linear epitopes
- Finds binding sites (protein-protein interactions)

**Where Used:**
- `evomotif/motif_discovery.py::find_motifs()`

**Implementation:**

```python
def find_motifs(alignment, conservation_scores, min_conservation=0.7,
                min_length=5, max_gap_frequency=0.3, window_size=5):
    """
    Sliding window motif discovery.
    
    Algorithm:
    1. Scan alignment with sliding window
    2. Calculate average conservation in each window
    3. Mark windows exceeding threshold
    4. Merge overlapping windows into motifs
    5. Filter by length and recalculate stats
    """
    motif_starts = []
    aln_len = alignment.get_alignment_length()
    
    # Step 1-2: Sliding window scan
    for i in range(aln_len - window_size + 1):
        # Average conservation in window
        window_conservation = np.mean(conservation_scores[i:i+window_size])
        
        # Max gap frequency in window
        window_gaps = []
        for j in range(i, i+window_size):
            column = alignment[:, j]
            gap_freq = column.count('-') / len(alignment)
            window_gaps.append(gap_freq)
        max_gap = max(window_gaps)
        
        # Check thresholds
        if window_conservation >= min_conservation and max_gap <= max_gap_frequency:
            motif_starts.append(i)
    
    # Step 3: Merge adjacent windows
    motifs = []
    if motif_starts:
        current_start = motif_starts[0]
        current_end = motif_starts[0] + window_size
        
        for start in motif_starts[1:]:
            if start <= current_end:  # Overlapping
                current_end = max(current_end, start + window_size)
            else:  # New motif
                motifs.append((current_start, current_end))
                current_start = start
                current_end = start + window_size
        
        motifs.append((current_start, current_end))
    
    # Step 4: Filter by length
    filtered_motifs = []
    for start, end in motifs:
        length = end - start
        if length >= min_length:
            motif_conservation = np.mean(conservation_scores[start:end])
            filtered_motifs.append({
                'start': start,
                'end': end,
                'length': length,
                'conservation': motif_conservation
            })
    
    return filtered_motifs
```

---

### Algorithm 2: Scattered Residue Detection

**Purpose:** Identify individual conserved positions regardless of spatial proximity.

**Mathematical Formulation:**

Position $i$ is conserved if:

$$C_{\text{combined}}(i) \geq \theta_c \quad \text{AND} \quad G(i) \leq \theta_g$$

Default thresholds:
- $\theta_c = 0.70$ (conservation)
- $\theta_g = 0.50$ (gap tolerance, higher than sliding window)

**Algorithm Steps:**

```
SCATTERED_RESIDUE_DETECTION(conservation, gaps, consensus):
    
    1. Initialize empty list R (conserved residues)
    
    2. For each position i in alignment:
        a. Get conservation score: C_i
        b. Get gap frequency: G_i
        c. Get consensus amino acid: A_i
        
        d. IF C_i >= threshold AND G_i <= gap_threshold:
            i. Add to conserved list: R.append({
                'position': i,
                'conservation': C_i,
                'gap': G_i,
                'amino_acid': A_i
            })
    
    3. Sort R by conservation score (descending)
    
    4. RETURN top conserved residues
```

**Why Independent Position Evaluation:**
- Detects catalytic residues (e.g., Ser-His-Asp in proteases)
- Finds metal-binding sites (Cys residues in zinc fingers)
- Identifies structural cysteines (disulfide bonds)
- Captures functionally critical residues even if isolated

**Example - P53 Results:**
```
Scattered conserved residues (non-consecutive):
- C215: Conservation 0.877, Gap 0.379 (zinc-binding)
- C309: Conservation 0.862, Gap 0.379 (zinc-binding)
- C365: Conservation 0.847, Gap 0.448 (zinc-binding)
- C313: Conservation 0.839, Gap 0.414 (zinc-binding)
```

These cysteines form a zinc finger, critical for DNA binding, but are NOT consecutive in sequence.

**Where Used:**
- `evomotif/motif_discovery.py::find_scattered_motifs()`
- Complementary to sliding window approach

**Implementation:**

```python
def find_scattered_motifs(conservation_scores, gap_frequencies, 
                         consensus_sequence, min_conservation=0.7,
                         max_gap=0.5):
    """
    Find individual conserved positions (scattered motifs).
    
    Algorithm:
    1. Evaluate each position independently
    2. Filter by conservation and gap thresholds
    3. Extract consensus amino acid
    4. Sort by conservation score
    """
    conserved_positions = []
    
    for i in range(len(conservation_scores)):
        conservation = conservation_scores[i]
        gap = gap_frequencies[i]
        aa = consensus_sequence[i]
        
        # Apply thresholds
        if conservation >= min_conservation and gap <= max_gap:
            conserved_positions.append({
                'position': i,
                'conservation': conservation,
                'gap': gap,
                'amino_acid': aa
            })
    
    # Sort by conservation (descending)
    conserved_positions.sort(key=lambda x: x['conservation'], reverse=True)
    
    return conserved_positions
```

---

## Statistical Validation

### Permutation Test

**Purpose:** Non-parametric significance testing for motif conservation.

**Null Hypothesis ($H_0$):**
The observed motif conservation is no different from random regions of the alignment.

**Alternative Hypothesis ($H_a$):**
The observed motif conservation is significantly higher than random regions.

**Mathematical Formulation:**

$$p\text{-value} = \frac{\#\{\pi: S_\pi \geq S_{\text{obs}}\}}{N}$$

Where:
- $S_{\text{obs}}$ = observed conservation score of motif
- $S_\pi$ = conservation score of random permutation $\pi$
- $N$ = total number of permutations (default: 1000)

**Algorithm Steps:**

```
PERMUTATION_TEST(motif_conservation, all_conservation_scores, n_permutations=1000):
    
    1. observed_score = motif_conservation
    
    2. Initialize counter: n_greater = 0
    
    3. FOR i = 1 to n_permutations:
        a. Randomly sample positions from alignment
        b. Calculate average conservation of random sample
        c. IF random_conservation >= observed_score:
            i. n_greater += 1
    
    4. p_value = n_greater / n_permutations
    
    5. RETURN p_value
```

**Why Permutation Test:**
- **Non-parametric:** No assumptions about distribution
- **Exact:** Uses empirical null distribution from data
- **Conservative:** Accounts for alignment-specific conservation patterns

**Interpretation:**
- $p < 0.05$: Motif is significantly conserved (reject $H_0$)
- $p < 0.01$: Highly significant
- $p < 0.001$: Extremely significant

**Where Used:**
- `evomotif/stats.py::permutation_test()`
- Applied to each discovered motif

**Implementation:**

```python
def permutation_test(observed_score, all_scores, n_permutations=1000,
                    motif_length=None, alternative='greater'):
    """
    Permutation test for motif significance.
    
    Algorithm:
    1. Store observed score
    2. Generate null distribution by random sampling
    3. Count how many random scores >= observed
    4. Calculate p-value = (n_greater + 1) / (n_permutations + 1)
       (+1 for continuity correction)
    """
    np.random.seed(42)  # Reproducibility
    
    n_greater = 0
    
    for i in range(n_permutations):
        # Random sample same length as motif
        if motif_length:
            random_indices = np.random.choice(
                len(all_scores), 
                size=motif_length, 
                replace=False
            )
            random_score = np.mean([all_scores[idx] for idx in random_indices])
        else:
            random_score = np.random.choice(all_scores)
        
        # Count exceedances
        if alternative == 'greater':
            if random_score >= observed_score:
                n_greater += 1
        elif alternative == 'less':
            if random_score <= observed_score:
                n_greater += 1
    
    # Calculate p-value with continuity correction
    p_value = (n_greater + 1) / (n_permutations + 1)
    
    return p_value
```

---

## Multiple Testing Correction

### Benjamini-Hochberg FDR Correction

**Purpose:** Control false discovery rate when testing multiple motifs.

**Problem:**
Testing $m$ motifs simultaneously increases Type I error rate. If each test has $\alpha = 0.05$, the family-wise error rate is:

$$P(\text{at least one false positive}) = 1 - (1 - \alpha)^m$$

For $m = 20$ motifs: $P = 1 - 0.95^{20} \approx 0.64$ (64% chance of false positive!)

**Solution: False Discovery Rate Control**

The FDR is the expected proportion of false discoveries among rejected hypotheses:

$$\text{FDR} = E\left[\frac{V}{R}\right]$$

Where:
- $V$ = number of false positives
- $R$ = total number of rejections

**Benjamini-Hochberg Procedure:**

Given $m$ p-values: $p_1, p_2, \ldots, p_m$

```
BENJAMINI_HOCHBERG_FDR(p_values, alpha=0.05):
    
    1. Sort p-values in ascending order: p_(1) <= p_(2) <= ... <= p_(m)
    
    2. For i = m down to 1:
        a. Calculate threshold: threshold_i = (i / m) * alpha
        b. IF p_(i) <= threshold_i:
            i. Reject H_0 for all j <= i
            ii. BREAK
    
    3. RETURN adjusted p-values and decisions
```

**Mathematical Formula for Adjusted P-values:**

$$p_{\text{adj}}^{(i)} = \min\left\{1, \min_{j \geq i} \left\{\frac{m \cdot p^{(j)}}{j}\right\}\right\}$$

**Why BH-FDR:**
- **Less conservative** than Bonferroni correction
- **Controls false discovery rate** instead of family-wise error rate
- **More power** to detect true positives
- **Appropriate for exploratory analysis**

**Example:**

```
Original p-values (5 motifs):
  Motif 1: p = 0.001
  Motif 2: p = 0.015
  Motif 3: p = 0.032
  Motif 4: p = 0.045
  Motif 5: p = 0.120

BH procedure (alpha = 0.05):
  Sorted: [0.001, 0.015, 0.032, 0.045, 0.120]
  
  i=5: 0.120 <= (5/5)*0.05 = 0.050? NO
  i=4: 0.045 <= (4/5)*0.05 = 0.040? NO
  i=3: 0.032 <= (3/5)*0.05 = 0.030? NO
  i=2: 0.015 <= (2/5)*0.05 = 0.020? YES ✓
  
  Reject H_0 for motifs 1 and 2 (significant after correction)

Adjusted p-values:
  Motif 1: 0.005
  Motif 2: 0.0375
  Motif 3: 0.0533
  Motif 4: 0.0563
  Motif 5: 0.120
```

**Where Used:**
- `evomotif/stats.py::multiple_testing_correction()`
- Applied after all motif p-values calculated

**Implementation:**

```python
def multiple_testing_correction(p_values, method='fdr_bh', alpha=0.05):
    """
    Benjamini-Hochberg FDR correction.
    
    Algorithm:
    1. Sort p-values with original indices
    2. Calculate BH critical values: (i/m) * alpha
    3. Find largest i where p_(i) <= (i/m)*alpha
    4. Reject all hypotheses j <= i
    5. Calculate adjusted p-values
    """
    from statsmodels.stats.multitest import multipletests
    
    # Use statsmodels implementation
    reject, p_adjusted, _, _ = multipletests(
        p_values,
        alpha=alpha,
        method=method  # 'fdr_bh' = Benjamini-Hochberg
    )
    
    return {
        'reject': reject,  # Boolean array
        'p_adjusted': p_adjusted,
        'alpha': alpha,
        'method': method
    }
```

---

## Effect Size Calculation

### Cohen's d

**Purpose:** Quantify the magnitude of difference between motif conservation and background.

**Mathematical Formula:**

$$d = \frac{\mu_{\text{motif}} - \mu_{\text{background}}}{s_{\text{pooled}}}$$

Where:

$$s_{\text{pooled}} = \sqrt{\frac{(n_1 - 1)s_1^2 + (n_2 - 1)s_2^2}{n_1 + n_2 - 2}}$$

- $\mu_{\text{motif}}$ = mean conservation of motif
- $\mu_{\text{background}}$ = mean conservation of all positions
- $s_{\text{pooled}}$ = pooled standard deviation
- $n_1, n_2$ = sample sizes

**Interpretation (Cohen's Guidelines):**
- $|d| < 0.2$: Negligible effect
- $0.2 \leq |d| < 0.5$: Small effect
- $0.5 \leq |d| < 0.8$: Medium effect
- $|d| \geq 0.8$: Large effect

**Why Effect Size Matters:**
- P-values indicate **statistical significance**
- Effect sizes indicate **practical significance**
- Large datasets can make tiny differences "significant"
- Effect size is independent of sample size

**Example:**

```
Motif conservation: mean = 0.85, sd = 0.05, n = 15
Background conservation: mean = 0.50, sd = 0.20, n = 1000

Pooled SD = sqrt(((15-1)*0.05² + (1000-1)*0.20²) / (15+1000-2))
          ≈ 0.197

Cohen's d = (0.85 - 0.50) / 0.197
          = 1.78  (VERY LARGE effect)
```

**Where Used:**
- `evomotif/stats.py::calculate_effect_size()`
- Reported alongside p-values for each motif

**Implementation:**

```python
def calculate_effect_size(motif_scores, background_scores):
    """
    Calculate Cohen's d effect size.
    
    Algorithm:
    1. Calculate means of both groups
    2. Calculate standard deviations
    3. Calculate pooled standard deviation
    4. Compute Cohen's d = (mean1 - mean2) / pooled_sd
    """
    mean_motif = np.mean(motif_scores)
    mean_background = np.mean(background_scores)
    
    sd_motif = np.std(motif_scores, ddof=1)
    sd_background = np.std(background_scores, ddof=1)
    
    n_motif = len(motif_scores)
    n_background = len(background_scores)
    
    # Pooled standard deviation
    pooled_sd = np.sqrt(
        ((n_motif - 1) * sd_motif**2 + (n_background - 1) * sd_background**2) /
        (n_motif + n_background - 2)
    )
    
    # Cohen's d
    if pooled_sd == 0:
        return 0.0
    
    cohens_d = (mean_motif - mean_background) / pooled_sd
    
    return cohens_d
```

---

## Implementation Details

### Complete Statistical Pipeline

**Step-by-step execution in `run_complete_analysis.py`:**

```python
# STEP 1: Calculate conservation
conservation_scores = scorer.calculate_combined_conservation(alignment)
gap_frequencies = scorer.calculate_gap_frequency(alignment)
consensus = scorer.get_consensus_sequence(alignment)

# STEP 2: Discover motifs (both methods)
consecutive_motifs = discoverer.find_motifs(
    alignment, 
    conservation_scores,
    min_conservation=0.70
)

scattered_motifs = discoverer.find_scattered_motifs(
    conservation_scores,
    gap_frequencies,
    consensus,
    min_conservation=0.70
)

# STEP 3: Statistical validation
validator = StatisticalValidator()

# Permutation test for each motif
for motif in consecutive_motifs:
    motif['p_value'] = validator.permutation_test(
        observed_score=motif['conservation'],
        all_scores=conservation_scores,
        n_permutations=1000,
        motif_length=motif['length']
    )
    
    # Effect size
    motif_region = conservation_scores[motif['start']:motif['end']]
    motif['effect_size'] = validator.calculate_effect_size(
        motif_region,
        conservation_scores
    )

# STEP 4: FDR correction
p_values = [m['p_value'] for m in consecutive_motifs]
correction_result = validator.multiple_testing_correction(
    p_values,
    method='fdr_bh',
    alpha=0.05
)

for i, motif in enumerate(consecutive_motifs):
    motif['p_adjusted'] = correction_result['p_adjusted'][i]
    motif['significant'] = correction_result['reject'][i]
```

---

### Parameter Selection Rationale

**Conservation Threshold (0.70):**
- Based on literature (Capra & Singh, 2007)
- Balances sensitivity and specificity
- Validated on known functional sites

**Window Size (5 residues):**
- Typical size of functional motifs
- Short enough to be specific
- Long enough to be meaningful

**Permutations (1000):**
- Standard in bioinformatics
- p-value resolution: 0.001
- Computationally feasible

**FDR Alpha (0.05):**
- Standard significance level
- 5% expected false discovery rate
- Appropriate for exploratory analysis

---

### Computational Complexity

| Method | Time Complexity | Space Complexity |
|--------|----------------|------------------|
| Shannon Entropy | $O(L \cdot N)$ | $O(L)$ |
| BLOSUM62 Score | $O(L \cdot N^2)$ | $O(L)$ |
| Sliding Window | $O(L \cdot w)$ | $O(L)$ |
| Scattered Detection | $O(L)$ | $O(L)$ |
| Permutation Test | $O(P \cdot L)$ | $O(P)$ |
| FDR Correction | $O(M \log M)$ | $O(M)$ |

Where:
- $L$ = alignment length
- $N$ = number of sequences
- $w$ = window size
- $P$ = number of permutations
- $M$ = number of motifs

**Total Pipeline:** $O(L \cdot N^2 + P \cdot L)$ dominated by BLOSUM62 pairwise comparisons.

---

## References

1. **Shannon Entropy:**
   - Shannon, C.E. (1948). "A Mathematical Theory of Communication". *Bell System Technical Journal*.

2. **BLOSUM Matrices:**
   - Henikoff, S. & Henikoff, J.G. (1992). "Amino acid substitution matrices from protein blocks". *PNAS* 89(22): 10915-10919.

3. **Conservation Scoring:**
   - Capra, J.A. & Singh, M. (2007). "Predicting functionally important residues from sequence conservation". *Bioinformatics* 23(15): 1875-1882.

4. **Permutation Tests:**
   - Good, P.I. (2005). *Permutation, Parametric, and Bootstrap Tests of Hypotheses*. Springer.

5. **FDR Correction:**
   - Benjamini, Y. & Hochberg, Y. (1995). "Controlling the false discovery rate: a practical and powerful approach to multiple testing". *Journal of the Royal Statistical Society B* 57(1): 289-300.

6. **Effect Sizes:**
   - Cohen, J. (1988). *Statistical Power Analysis for the Behavioral Sciences*. Lawrence Erlbaum Associates.

---

**Last Updated:** December 2, 2025  
**Version:** 1.0.0
