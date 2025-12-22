# EvoMotif Quick Reference Card

## 🚀 The New Way (2 Lines!)

```python
import evomotif
results = evomotif.analyze_protein("p53", "your@email.com")
```

---

## 📊 What Changed?

| Aspect | Before | After |
|--------|--------|-------|
| **Import** | 7+ modules | 1 module |
| **Code** | 150+ lines | 2 lines |
| **Learning** | Hours | Minutes |
| **Pattern** | Custom | Standard library style |

---

## ✅ Yes, it works like a standard library!

### Like NumPy:
```python
import numpy as np
arr = np.array([1, 2, 3])
mean = np.mean(arr)
```

### Like Pandas:
```python
import pandas as pd
df = pd.read_csv('data.csv')
summary = df.describe()
```

### Like Scikit-learn:
```python
from sklearn.ensemble import RandomForestClassifier
model = RandomForestClassifier()
model.fit(X, y)
predictions = model.predict(X_test)
```

### Like EvoMotif (NEW!):
```python
import evomotif
results = evomotif.analyze_protein("p53", "user@email.com")
summary = results.summary()
```

---

## 📖 Common Usage Patterns

### Basic Analysis
```python
import evomotif
results = evomotif.analyze_protein("ubiquitin", "user@email.com")
print(results.summary())
```

### With 3D Structure
```python
results = evomotif.analyze_protein(
    "p53", 
    "user@email.com",
    pdb_id="1TUP"
)
```

### Custom Parameters
```python
results = evomotif.analyze_protein(
    "BRCA1",
    "user@email.com",
    max_sequences=100,
    min_conservation=0.65,
    threads=8,
    verbose=True
)
```

### Check Dependencies
```python
pipeline = evomotif.EvoMotifPipeline()
deps = pipeline.check_dependencies()
for tool, available in deps.items():
    print(f"{tool}: {'✓' if available else '✗'}")
```

### Access Results
```python
results.protein              # Protein name
results.n_sequences          # Number of sequences
results.motifs               # List of motifs
results.conserved_positions  # Conserved residues
results.summary()            # Human-readable summary
results.export_json(path)    # Export to JSON
results.get_file('alignment') # Get file path
```

---

## 🎯 Key Features

✅ **Single import** - `import evomotif`  
✅ **Clean API** - One function does everything  
✅ **Type hints** - Full IDE support  
✅ **Result objects** - Easy data access  
✅ **Error handling** - Clear, helpful messages  
✅ **Automatic checks** - Dependencies validated  
✅ **Progress feedback** - Know what's happening  
✅ **Documentation** - `help()` works everywhere  
✅ **Jupyter ready** - Tab completion + inline docs  
✅ **Standard pattern** - Like NumPy, Pandas, sklearn  

---

## 🔧 For Power Users

The old modular way still works:

```python
from evomotif.retrieval import SequenceRetriever
from evomotif.alignment import SequenceAligner
# ... full control over every step
```

---

## 📚 Documentation

- **Quick Start**: [GETTING_STARTED.md](GETTING_STARTED.md)
- **Examples**: [examples/simple_api_demo.py](examples/simple_api_demo.py)
- **Full Guide**: [USER_GUIDE.md](USER_GUIDE.md)
- **Modules**: [PIPELINE_GUIDE.md](PIPELINE_GUIDE.md)

---

## 💡 Compare to Other Tools

| Tool | Lines | Complexity |
|------|-------|------------|
| MEME Suite (CLI) | Command-line | High |
| Pfam (Database) | Web interface | Medium |
| EvoMotif (OLD) | 150+ | High |
| **EvoMotif (NEW)** | **2** | **Low** |

---

**Bottom Line:** EvoMotif now works exactly like NumPy, Pandas, scikit-learn, and other standard Python libraries! 🎉
