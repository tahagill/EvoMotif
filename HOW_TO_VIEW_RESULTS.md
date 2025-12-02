# 📊 How to View EvoMotif Results

## 🧬 3D Structure Visualizations (Interactive!)

### The EASIEST Way - Already Built-In! ✅

The `*_structure_3d.html` files are **interactive 3D molecular viewers**. Just open them in any web browser!

**Files Generated:**
- `complete_analysis_results/p53/p53_structure_3d.html` (527KB)
- `complete_analysis_results/BRCA1/BRCA1_structure_3d.html` (3.5MB)
- `complete_analysis_results/AKT1/AKT1_structure_3d.html` (283KB)

**How to Open:**

```bash
# Option 1: Firefox
firefox complete_analysis_results/AKT1/AKT1_structure_3d.html

# Option 2: Chrome
google-chrome complete_analysis_results/p53/p53_structure_3d.html

# Option 3: Any browser
xdg-open complete_analysis_results/BRCA1/BRCA1_structure_3d.html

# Option 4: Just double-click the HTML file in your file manager!
```

**What You'll See:**
- ✅ Rotating 3D protein structure
- ✅ Conserved residues highlighted in RED (shown as sticks)
- ✅ Protein backbone in gray (cartoon view)
- ✅ Interactive controls:
  - **Left-click + drag** = Rotate
  - **Right-click + drag** = Zoom
  - **Middle-click + drag** = Pan
  - **Scroll wheel** = Zoom in/out

---

## 🌳 Phylogenetic Trees

### Option 1: Online Tree Viewers (EASIEST!)

Upload your `.nwk` files to these free online tools:

**🌐 iTOL (Interactive Tree of Life) - RECOMMENDED**
- Website: https://itol.embl.de/upload.cgi
- Upload: `complete_analysis_results/AKT1/AKT1_tree.nwk`
- Features: Beautiful visualizations, export as PDF/PNG/SVG

**🌐 Phylo.io**
- Website: http://phylo.io
- Drag & drop your `.nwk` file
- Features: Simple, clean, fast

**🌐 ETE Toolkit TreeView**
- Website: http://etetoolkit.org/treeview/
- Paste Newick content or upload file
- Features: Scientific-quality figures

### Option 2: Desktop Software

**FigTree (RECOMMENDED for publications)**
```bash
# Install
sudo apt-get install figtree

# Open tree
figtree complete_analysis_results/p53/p53_tree.nwk
```

**Dendroscope**
```bash
# Download from: https://software-ab.informatik.uni-tuebingen.de/download/dendroscope/welcome.html
# Then open your .nwk file
```

### Option 3: Python/Jupyter Notebook

**Quick Viewer (using ete3):**
```python
from ete3 import Tree

# Load tree
tree = Tree("complete_analysis_results/AKT1/AKT1_tree.nwk")

# Display ASCII art
print(tree)

# Or render to file
tree.render("AKT1_tree.png")
```

**Using BioPython:**
```python
from Bio import Phylo
import matplotlib.pyplot as plt

# Load and display
tree = Phylo.read("complete_analysis_results/p53/p53_tree.nwk", "newick")
Phylo.draw(tree)
plt.show()
```

---

## 📁 PDB Structure Files

The raw PDB files can be opened with:

### PyMOL (Professional - FREE for Education)
```bash
# Install
sudo apt-get install pymol

# Open
pymol complete_analysis_results/AKT1/3CQU.pdb
```

**In PyMOL:**
- Color by conservation (manually)
- Create publication-quality figures
- Advanced structural analysis

### UCSF Chimera/ChimeraX (FREE)
```bash
# Download from: https://www.cgl.ucsf.edu/chimerax/
chimerax complete_analysis_results/BRCA1/1JM7.pdb
```

### VMD (Visual Molecular Dynamics - FREE)
```bash
# Install
sudo apt-get install vmd

# Open
vmd complete_analysis_results/p53/1TUP.pdb
```

### Online Viewers
- **RCSB PDB 3D View**: https://www.rcsb.org/3d-view/
- **Mol* Viewer**: https://molstar.org/viewer/
- **NGL Viewer**: http://nglviewer.org/ngl/

---

## 📊 Conservation Data

### JSON Files

**View conservation scores:**
```bash
# Pretty print JSON
jq '.' complete_analysis_results/AKT1/AKT1_conservation.json | less

# Get top 10 conserved positions
jq '. | sort_by(-.conservation) | .[:10]' complete_analysis_results/p53/conserved_positions.json
```

**Plot in Python:**
```python
import json
import matplotlib.pyplot as plt

# Load conservation data
with open('complete_analysis_results/AKT1/conserved_positions.json') as f:
    data = json.load(f)

# Extract positions and scores
positions = [d['position'] for d in data]
conservation = [d['conservation'] for d in data]

# Plot
plt.figure(figsize=(15, 5))
plt.plot(positions, conservation, 'b-', alpha=0.7)
plt.fill_between(positions, conservation, alpha=0.3)
plt.axhline(y=0.7, color='r', linestyle='--', label='Threshold')
plt.xlabel('Position')
plt.ylabel('Conservation Score')
plt.title('AKT1 Conservation Profile')
plt.legend()
plt.grid(True, alpha=0.3)
plt.savefig('AKT1_conservation_plot.png', dpi=300, bbox_inches='tight')
plt.show()
```

---

## 🎨 Creating Publication Figures

### 1. 3D Structure Figure

**From HTML (Screenshot):**
- Open the HTML file in browser
- Rotate to desired view
- Take screenshot (Print Screen or browser dev tools)

**From PyMOL (Professional Quality):**
```pymol
# Open PyMOL
pymol complete_analysis_results/AKT1/3CQU.pdb

# PyMOL commands:
hide everything
show cartoon
color gray, all
# Manually select conserved residues and color red
select conserved, resi 234+250+345+631+710+765+766+854
show sticks, conserved
color red, conserved
bg_color white
ray 1200, 800
png AKT1_structure_figure.png
```

### 2. Phylogenetic Tree Figure

**From iTOL:**
1. Upload `.nwk` file to https://itol.embl.de
2. Customize colors, labels, bootstrap values
3. Export as PDF/SVG/PNG (publication quality)

**From FigTree:**
1. Open in FigTree
2. Adjust layout (radial, rectangular, etc.)
3. File → Export PNG (set DPI to 300)

### 3. Conservation Plot

Use the Python code above to create a conservation profile plot.

---

## 🚀 Quick Start Examples

### View Everything for AKT1
```bash
# 3D structure (interactive)
firefox complete_analysis_results/AKT1/AKT1_structure_3d.html &

# Conservation scores
jq '.[:20]' complete_analysis_results/AKT1/conserved_positions.json

# Tree (upload to iTOL)
cat complete_analysis_results/AKT1/AKT1_tree.nwk
# Then paste at: https://itol.embl.de/upload.cgi

# Summary
jq '.' complete_analysis_results/AKT1/AKT1_summary.json | less
```

### View Everything for P53
```bash
# Open 3D structure
xdg-open complete_analysis_results/p53/p53_structure_3d.html

# Show tree in ASCII
python3 << 'EOF'
from Bio import Phylo
tree = Phylo.read("complete_analysis_results/p53/p53_tree.nwk", "newick")
Phylo.draw_ascii(tree)
EOF

# Top conserved positions
jq 'sort_by(-.conservation) | .[:10]' complete_analysis_results/p53/conserved_positions.json
```

---

## 📦 Install Recommended Tools

```bash
# Browser (usually already installed)
sudo apt-get install firefox

# Tree viewers
sudo apt-get install figtree

# Structure viewers
sudo apt-get install pymol
# or
wget https://www.cgl.ucsf.edu/chimerax/...  # ChimeraX

# Python visualization
pip install matplotlib seaborn biopython ete3
```

---

## 💡 Pro Tips

### 1. Interactive HTML is Self-Contained
The `*_structure_3d.html` files contain EVERYTHING:
- 3D viewer (py3Dmol)
- Structure coordinates
- Styling and colors
- No internet required after download!

### 2. Convert Newick to Image Quickly
```bash
# Python one-liner
python3 -c "from Bio import Phylo; import matplotlib.pyplot as plt; tree = Phylo.read('complete_analysis_results/AKT1/AKT1_tree.nwk', 'newick'); Phylo.draw(tree); plt.savefig('tree.png', dpi=300, bbox_inches='tight')"
```

### 3. Batch View All Structures
```bash
# Open all HTML files
for html in complete_analysis_results/*/*.html; do
    firefox "$html" &
    sleep 2
done
```

### 4. Create Animated Structure
Use PyMOL scripting to rotate and save frames, then combine with ffmpeg.

---

## 🎯 What Each File Type Is For

| File Type | Purpose | Best Viewer |
|-----------|---------|-------------|
| `*.html` | Interactive 3D viewer | **Web browser** |
| `*.nwk` | Phylogenetic tree | **iTOL** or FigTree |
| `*.pdb` | Raw structure coordinates | PyMOL, ChimeraX |
| `*.json` | Conservation data | Python, jq |
| `*.fasta` | Sequence alignments | Jalview, MEGA |

---

## 🆘 Troubleshooting

### HTML file won't open
- Try different browser (Firefox, Chrome, Edge)
- Check file size (should be >100KB)
- Open in VS Code and check for errors

### Tree looks messy
- Use online viewer (iTOL handles large trees better)
- Try different layout (radial vs rectangular)
- Increase figure size

### Can't see conserved residues in HTML
- Make sure you have conserved positions (check JSON file)
- File might have failed during generation
- Check console in browser (F12) for errors

---

**Need help?** Check the tool's GitHub issues or documentation!

🎨 **Happy visualizing!** 🧬
