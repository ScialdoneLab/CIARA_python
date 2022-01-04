import numpy as np

def get_background_full(norm_adata, threshold, n_cells, n_cells_high):

    thr_per_gene = np.sum(norm_adata.X > threshold, axis=0)
    genes_filter = np.logical_and(thr_per_gene >= n_cells, thr_per_gene <= n_cells_high)
    print("Background genes: " + str(np.sum(genes_filter)))

    return genes_filter.tolist()
