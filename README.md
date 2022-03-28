# CIARA <br /> (Cluster Independent Algorithm for the identification of RAre cell types)

Python implementation of the CIARA algorithm that integrates into scanpy's analysis with the AnnData format.

The package can be installed via pip:

`python -m pip install ciara_python`

*Note: The package only works on UNIX / MacOS operating systems, not on Windows systems, due to the copy-on-write multiprocessing setup used.*

## Tutorial

The functions are designed to work with scanpy's AnnData objects. For an interactive tutorial check out the Human Gastrula IPython Notebook which is also part of this repository.

First you load your dataset as a scanpy object and after normal preprocessing calculate the knn-network:

```python
import scanpy as sc

pbmc = sc.datasets.pbmc3k()
sc.pp.filter_cells(pbmc, min_genes=50)
sc.pp.filter_genes(pbmc, min_cells=10)
sc.pp.log1p(pbmc)

sc.pp.pca(pbmc)
sc.pp.neighbors(pbmc)
```

The CIARA package contains the two main functions `get_full_background()` and `ciara()` which should be imported via:

```python
from ciara_python import get_ background_full, ciara
```

Afterwards, the background genes get marked by running the `get_full_background()` function on your scanpy dataset. This adds the boolean column 'CIARA_background' to your `pbmc.var` AnnData slot, where relevant background genes are marked.

```python
get_background_full(pbmc, threshold=1, n_cells=2, n_cells_high=5)
```

Finally, the `ciara()` function is run on the dataset. This adds the column 'CIARA_p_value' to your `pbmc.var` object, where the calculated p_values for each of the previously marked background genes are stored.

```python
ciara(pbmc, n_cores=4, p_value=0.001, odds_ratio=2, approximation=True, local_region=1)
```

We can then extract the top markers or rare cells (lowest CIARA_p_value) and plot them onto the UMAP:

```python
from matplotlib.pyplot import rc_context

sc.tl.umap(pbmc)

top_markers = pbmc.var.nsmallest(4, ["CIARA_p_value"],)

with rc_context({'figure.figsize': (3, 3)}):
    sc.pl.umap(pbmc, color=top_markers.index.tolist())
```

## R package

Link to R package:

https://github.com/ScialdoneLab/CIARA/
