import numpy as np
# in each connectivity map each [1,2,3,4] is a 'yarn'


def get_evals(huckel):
    evals, evecs = np.linalg.eig(huckel)
    return(sorted(evals))


def gen_lin_yarn(length):
    return ([n for n in range(length)])


def gen_cyc_yarn(length):
    return([n % length for n in range(length + 1)])


def true2one(a):
    return (1 if a else 0)


def adj_list_to_adj_matrix(adjency_list):
    leng = len(adjency_list)
    adgency_array = [[true2one(n in neighbors) for n in range(leng)]
                     for neighbors in adjency_list]
    return(np.matrix(adgency_array))


def yarn_to_adj_list(connection_map):
    adjency1list = []
    for n in range(max(connection_map) + 1):
        adjency1list.append([])
    for index, node in enumerate(connection_map):
        if index + 1 < len(connection_map):
            adjency1list[node].append(connection_map[index + 1])
            adjency1list[connection_map[index + 1]].append(node)
    return(adjency1list)


def yarn_list_to_adj_list(connection_map):
    listofmax = [max(single_list) for single_list in connection_map]
    adjency_list = []
    for n in range(max(listofmax) + 1):
        adjency_list.append([])
    for single_list in connection_map:
        single_adjency_list = yarn_to_adj_list(single_list)
        for node, neighbors in enumerate(single_adjency_list):
            adjency_list[node] += neighbors
    return (adjency_list)


def merge_degeneracies(evals):
    current_eval_index = -1
    merged_evals = []
    current_eval = False
    for eigenvalue in evals:
        if (current_eval and current_eval < eigenvalue + 10 **
                (-12)) and (current_eval > eigenvalue - 10**(-12)):
            merged_evals[current_eval_index].append(eigenvalue)
        else:
            current_eval = eigenvalue
            current_eval_index += 1
            merged_evals.append([])
            merged_evals[current_eval_index].append(eigenvalue)
    return(merged_evals)


def sorted_degeneracies(mer_deg):
    sort_deg = [[eigval[0], len(eigval)] for eigval in mer_deg]
    return(sort_deg)


def arbitrary_to_numeral(arb_list):
    flatten_list = [node for yarn in arb_list for node in yarn]
    unique_flatten_list = list(set(flatten_list))
    numeral_list = [[unique_flatten_list.index(
        node) for node in yarn] for yarn in arb_list]
    return (numeral_list)


def kron_delta(a, b):
    if a == b:
        return(1)
    else:
        return(0)


def connectivity_to_eigenvalues(yarn_list):
    num_yarn_list = arbitrary_to_numeral(yarn_list)
    adj_list = yarn_list_to_adj_list(num_yarn_list)
    adj_matrix = adj_list_to_adj_matrix(adj_list)
    huckel_matrix = -1 * adj_matrix
    evals = get_evals(huckel_matrix)
    merged_eval = merge_degeneracies(evals)
    sorted_eig = sorted_degeneracies(merged_eval)
    return(sorted_eig)


def poly_ene_huckel(n):
    if n < 1:
        raise ValueError('polyene size must be 1 or greater')
    yarn = gen_lin_yarn(n)
    yarn_list = [yarn]
    return(connectivity_to_eigenvalues(yarn_list))


def cyclicpolyenehuckel(n):
    if n < 1:
        raise ValueError('polyene size must be 1 or greater')
    yarn = gen_cyc_yarn(n)
    yarn_list = [yarn]
    return(connectivity_to_eigenvalues(yarn_list))
