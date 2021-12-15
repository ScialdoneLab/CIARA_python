import numpy as np
from scipy.stats import fisher_exact, entropy
import multiprocessing
from functools import partial

def perform_sub_entropy(nn_gene_expression, binary_expression, p_value, odds_ratio=2):

    entropy_nn, p_value_nn = 1, 1 
    
    if np.sum(nn_gene_expression) > 0:
        input_fisher = np.array([[np.sum(nn_gene_expression), np.sum(binary_expression)-np.sum(nn_gene_expression)],
                                [np.sum(~nn_gene_expression), np.sum(~binary_expression)-np.sum(~nn_gene_expression)]])
        oddsr_test, p_test = fisher_exact(input_fisher, alternative = 'greater')
        if p_test < p_value and oddsr_test > odds_ratio:
            p_0 = np.sum(nn_gene_expression) / len(nn_gene_expression)
            entropy_nn = entropy([p_0, 1 - p_0], base=2)
            p_value_nn = p_test

    return entropy_nn, p_value_nn

def entropy_mixing_gene(gene_expression, knn_matrix, p_value, odds_ratio, local_region, approximation):

    binary_expression = gene_expression > np.median(gene_expression)
    
    if approximation:
        knn_matrix = knn_matrix[binary_expression,:]
    
    entropies_nn = np.array([])
    p_values_nn = np.array([])
    for cell in range(knn_matrix.shape[0]):
        nn_cell = knn_matrix[cell, :]
        nn_gene_expression = binary_expression[nn_cell==1]
        entropy_sub, p_value_sub = perform_sub_entropy(nn_gene_expression, binary_expression , p_value, odds_ratio)
        entropies_nn = np.append(entropies_nn, entropy_sub)
        p_values_nn = np.append(p_values_nn, p_value_sub)

    if np.sum(entropies_nn < 1) >= local_region:
        entropy_gene = np.min(entropies_nn)
        p_value_gene = np.min(p_values_nn)
    else:
        entropy_gene = 1
        p_value_gene = 1

    return entropy_gene, p_value_gene

def entropy_mixing(norm_adata, knn_matrix, n_cores, p_value, odds_ratio, local_region, approximation):
    
    assert((knn_matrix.index.values == knn_matrix.columns.values).all)
    assert((norm_adata.obs_names == knn_matrix.index.values).all)
    #knn_matrix is pandas dataframe in same orientation as Anndata.X (cellsxgenes)
    #knn_matrix and Anndata.X should contain the cells in same order
    
    background = norm_adata.X[:, norm_adata.var["EOM_background"]]
    gene_expressions = [background[:,i].flatten() for i in range(np.shape(background)[1])]
    knn_matrix = knn_matrix.to_numpy()
    
    pool = multiprocessing.Pool(n_cores)
    temp = partial(entropy_mixing_gene, knn_matrix=knn_matrix, p_value=p_value, odds_ratio=odds_ratio, local_region=local_region, approximation=approximation)
    results = pool.map(func=temp, iterable=gene_expressions, chunksize=50)
    pool.close()
    pool.join()
    
    entropies_output = [np.NAN for i in range(len(norm_adata.var_names))]
    p_values_output = [np.NAN for i in range(len(norm_adata.var_names))]
    for index, gene_pos in enumerate(np.where(norm_adata.var["EOM_background"])[0]):
        entropies_output[gene_pos] = results[index][0]
        p_values_output[gene_pos] = results[index][1]
    
    print('\n---- Finished sucessfully! ----')

    return entropies_output, p_values_output
