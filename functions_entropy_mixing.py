import numpy as np
from scipy.stats import fisher_exact, entropy

def get_background_full(norm_matrix, threshold, n_cells, n_cells_high):

    #norm_matrix is Anndata object
    thr_per_gene = np.sum(norm_matrix.X > threshold, axis=2)
    genes_filter = thr_per_gene >= n_cells & thr_per_gene <= n_cells_high
    genes_important = norm_matrix.var_names[genes_filter.values]

    return genes_important.tolist()

def perform_sub_entropy(obs_names, cell, binary_expression, knn_matrix_sub, p_value, odds_ratio=2):

    nn_cell = obs_names[knn_matrix_sub[cell, :]]
    nn_gene_expression = binary_expression[set(obs_names) & set(nn_cell)]
    entropy_nn, p_value_nn = 1, 1

    if np.sum(nn_gene_expression) > 0:
        input_fisher = np.array([[np.sum(nn_gene_expression), np.sum(binary_expression)-np.sum(nn_gene_expression)],
                                [np.sum(~nn_gene_expression), np.sum(~binary_expression)-np.sum(~nn_gene_expresssion)]])
        oddsr_test, p_test = fisher_exact(input_fisher, alternative = 'greater')
        if p_test < p_value & oddsr_test > odds_ratio:
            p_0 = np.sum(nn_gene_expression) / len(nn_gene_expression)
            entropy_nn = entropy([p_0, 1 - p_0])
            p_value_nn = p_test

    return entropy_nn, p_value_nn

def entropy_mixing_gene(obs_names, knn_matrix, gene_expression, p_value, odds_ratio, local_region=2, approximation):

    binary_expression = gene_expression > np.median(gene_expression)

    if approximation:
        sub_feature = obs_names[binary_expression]
    else:
        sub_feature = obs_names

    knn_matrix_sub = knn_matrix[sub_feature, :]

    entropies_nn = []
    p_values_nn = []
    for cell in sub_feature:
        entropy_sub, p_value_sub = perform_sub_entropy(obs_names, cell, binary_expression, knn_matrix_sub, p_value, odds_ratio)
        entropies_nn.append(entropy_sub)
        p_values_nn.append(p_value_sub)

    if np.sum(entropy_nn < 1) >= local_region:
        entropy_gene = min(entropies_nn)
        p_value_gene = min(p_values_nn)
    else:
        entropy_gene = 1
        p_value_gene = 1

    return entropy_gene, p_value_gene

def entropy_mixing(norm_matrix, knn_matrix, background, cores_number, p_value, odds_ratio, local_region, approximation):

    #norm_matrix is Anndata object
    #knn_matrix is numpy array in same orientation as Anndata (cellsxgenes)
    obs_names = norm_matrix.obs_names.to_numpy()
    results = {} #dictionary
    for number, gene in enumerate(background):
        gene_expression = norm_matrix.X[gene].to_numpy()
        print('Gene: ' + gene + ' number: ' + (number + 1)
        results[gene] = entropy_mixing_gene(obs_names, knn_matrix, gene_expression, p_value, odds_ratio, local_region, approximation)

    return results
