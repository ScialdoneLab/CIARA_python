import numpy as np
from scipy.stats import fisher_exact
import scipy.sparse as sp
import multiprocessing
from functools import partial

def perform_fisher(nn_gene_expression, binary_expression, p_value, odds_ratio=2):

    p_value_nn = 1

    if np.sum(nn_gene_expression) != 0:
        input_fisher = np.array([[np.sum(nn_gene_expression), np.sum(binary_expression)-np.sum(nn_gene_expression)],
                                [np.sum(~nn_gene_expression), np.sum(~binary_expression)-np.sum(~nn_gene_expression)]])
        oddsr_test, p_test = fisher_exact(input_fisher, alternative = 'greater')
        if p_test < p_value and oddsr_test > odds_ratio:
            p_0 = np.sum(nn_gene_expression) / len(nn_gene_expression)
            p_value_nn = p_test

    return p_value_nn

def ciara_gene(gene_idx, p_value, odds_ratio, local_region, approximation):

    gene_expression = gene_expressions_g[gene_idx]
    binary_expression = gene_expression > np.median(gene_expression)

    if approximation:
        knn_subset = np.nditer(np.where(binary_expression))
    else:
        knn_subset = range(knn_matrix_g.shape[0])

    p_values_nn = np.array([])
    for cell in knn_subset:
        nn_cell = knn_matrix_g[cell, :]
        nn_gene_expression = binary_expression[nn_cell==1]
        p_value_sub = perform_fisher(nn_gene_expression, binary_expression , p_value, odds_ratio)
        p_values_nn = np.append(p_values_nn, p_value_sub)

    if np.sum(p_values_nn<p_value) >= local_region:
        p_value_gene = np.min(p_values_nn)
    else:
        p_value_gene = 1

    return p_value_gene

def ciara(norm_adata, n_cores, p_value, odds_ratio, local_region, approximation):

    multiprocessing.set_start_method("fork", force=True)

    background = norm_adata.X[:, norm_adata.var["CIARA_background"]]
    if sp.issparse(background):
        background = background.toarray()
    global gene_expressions_g
    gene_expressions_g = [background[:,i].flatten() for i in range(np.shape(background)[1])]
    global knn_matrix_g
    knn_matrix_g = norm_adata.obsp["connectivities"].toarray()

    pool = multiprocessing.Pool(n_cores)
    chunksize, extra = divmod(len(gene_expressions_g), 4 * n_cores)
    if extra:
        chunksize += 1
    print("\n## Running on " + str(n_cores) + " cores with a chunksize of " + str(chunksize))
    temp = partial(ciara_gene, p_value=p_value, odds_ratio=odds_ratio, local_region=local_region, approximation=approximation)
    results = pool.map(func=temp, iterable=range(len(gene_expressions_g)), chunksize=chunksize)
    pool.close()
    pool.join()


    p_values_output = [np.NAN for i in range(len(norm_adata.var_names))]
    for index, gene_pos in enumerate(np.where(norm_adata.var["CIARA_background"])[0]):
        p_values_output[gene_pos] = results[index]

    norm_adata.var.drop(columns="CIARA_p_value", inplace=True, errors='ignore')
    norm_adata.var.insert(0, "CIARA_p_value", p_values_output)

    print('\n---- Finished sucessfully! ----')

    return
