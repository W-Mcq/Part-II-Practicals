import numpy as np

i = "ana"
print(i)

M = np.ndarray((3,3))
print(M)

propene = np.matrix([[0,1,0],
	[1,0,1],
	[0,1,0]])

print(propene)

evals, evecs = np.linalg.eig(propene)

print(sorted(evals))

def get_evals(huckel):
    evals, evecs = np.linalg.eig(huckel)
    return(sorted (evals))

cyclopropene = np.matrix([
    [0,1,1],
    [1,0,1],
    [1,1,0]])

outy = get_evals(cyclopropene)

print(outy)

oops = [(i == 2) for i in range(4)]
print(oops)

def linear_connectivity(atom, neighbor):
    return (-1 if atom == neighbor+1 or atom == neighbor -1 else 0)

oops2 = [linear_connectivity(2, i) for i in range(5)]
print(oops2)

def poly_ene_huckel(n):
    initial_array = [[linear_connectivity(j, i) for i in range(n)] for j in range (n)]
    return(np.matrix(initial_array))

print(poly_ene_huckel(7))

def cyclic_connectivity(atom, neighbor,length):
    return (-1 if atom == (neighbor+1)%length or atom == (neighbor -1)%length else 0)
def cyclicpolyenehuckel(n):
    initial_array = [[cyclic_connectivity(j, i, n) for i in range(n)] for j in range (n)]
    return(np.matrix(initial_array))
print(cyclicpolyenehuckel(6))

print(np.linalg.eig(cyclicpolyenehuckel(6)))

def true2one(a):
    return (1 if a else 0)

def adjency_list_to_adgency_matrix(adjency_list):
    leng = len(adjency_list)
    adgency_array = [[true2one(n in neighbors) for n in range(leng)] for neighbors in adjency_list]
    return(np.matrix(adgency_array))

test_adj_list = [[1,4],[0,4,2],[1,3],[2,4],[3,1,0]]
print( adjency_list_to_adgency_matrix(test_adj_list))

def single_connect_list_to_adjency_list(connection_map):
    adjency1list = []
    for n in range(max(connection_map)+1):
        adjency1list.append([])
    for index, node in enumerate(connection_map):
        if index+1 < len(connection_map):
            adjency1list[node].append(connection_map[index+1])
            adjency1list[connection_map[index+1]].append(node)
    return(adjency1list)

print(2)
humma = single_connect_list_to_adjency_list([0,1,2,3,0,2,1,3])
habba = adjency_list_to_adgency_matrix(humma)
print (humma)
print (habba)

def multi_connect_list_to_adjency_list(connection_map):
    listofmax = [max(single_list) for single_list in connection_map]
    adjency_list = []
    for n in range(max(listofmax)+1):
        adjency_list.append([])
    for single_list in connection_map:
        single_adjency_list = single_connect_list_to_adjency_list(single_list)
        for node, neighbors in enumerate(single_adjency_list):
            adjency_list[node] += neighbors
    return (adjency_list)
        
huny = multi_connect_list_to_adjency_list([[0,1,2,0],[0,3,2],[1,3]])
hobo = adjency_list_to_adgency_matrix(huny)

print(huny)
print(hobo)
cube_adgency_list = multi_connect_list_to_adjency_list([[0,1,2,3,0],[4,5,6,7,4],[0,4],[1,5],[2,6],[3,7]])
cube_matrix_list = adjency_list_to_adgency_matrix(cube_adgency_list)
print(cube_matrix_list)

octahedron_adgency_list = multi_connect_list_to_adjency_list([[0,1,2,0],[3,4,5,3],[0,4,1,5,2,3,0]])
print(octahedron_adgency_list)
octahedron_adgency_matrix = adjency_list_to_adgency_matrix(octahedron_adgency_list)
print(octahedron_adgency_matrix)
print(get_evals(-1*octahedron_adgency_matrix))

dodecahedron_connectivity_list = [[1,2,3,4,5,1],[6,7,8,9,10,11,12,13,14,15,6],[0,16,17,18,19,0],
                                  [1,8],[2,10],[3,12],[4,14],[5,6],[0,13],[19,11],[18,9],[17,7],[16,15]] ##connectivity determined using a numbered 20-sided die
dodecahedron_adgency_list = multi_connect_list_to_adjency_list(dodecahedron_connectivity_list)
print(dodecahedron_adgency_list)
dodecahedron_adgency_matrix = adjency_list_to_adgency_matrix(dodecahedron_adgency_list)
print(dodecahedron_adgency_matrix)
print(get_evals(-1*dodecahedron_adgency_matrix))

icosahedron_connectivity_list = [[0,1,2,3,4,0],[5,6,7,8,9,5],
                                 [10,0],[10,1],[10,2],[10,3],[10,4],
                                 [11,5],[11,6],[11,7],[11,8],[11,9],
                                 [0,5,1,6,2,7,3,8,4,9,0]]
icosahedron_adgency_list = multi_connect_list_to_adjency_list(icosahedron_connectivity_list)
print(icosahedron_adgency_list)
icosahedron_adgency_matrix = adjency_list_to_adgency_matrix(icosahedron_adgency_list)
print(icosahedron_adgency_matrix)
print(get_evals(-1*icosahedron_adgency_matrix))

bucky_connectivity_list = [[0,1,2,3,4,0],
                           [5,6,7,8,9,5],
                           [10,11,12,13,14,10],
                           [15,16,17,18,19,15],
                           [20,21,22,23,24,20],
                           [25,26,27,28,29,25],
                           [30,31,32,33,34,30],
                           [35,36,37,38,39,35],
                           [40,41,42,43,44,40],
                           [45,46,47,48,49,45],
                           [50,51,52,53,54,50],
                           [55,56,57,58,59,55],
                           [0,46],[1,43],[2,56],[3,11],[4,5],
                           [6,47],[7,54],[8,16],[9,10],
                           [12,57],[13,21],[14,15],
                           [17,53],[18,26],[19,20],
                           [22,58],[23,31],[24,25],
                           [27,52],[28,36],[29,30],
                           [32,59],[33,41],[34,35],
                           [37,51],[38,49],[39,40],
                           [42,55],[44,45],[48,50]]
bucky_adgency_list = multi_connect_list_to_adjency_list(bucky_connectivity_list)
print(bucky_adgency_list)
bucky_adgency_matrix = adjency_list_to_adgency_matrix(bucky_adgency_list)
print(bucky_adgency_matrix)
print(get_evals(-1*bucky_adgency_matrix))

def merge_degeneracies(evals):
    current_eval_index = -1
    merged_evals = []
    current_eval = -10000.00
    for eigenvalue in evals:
        if (current_eval < eigenvalue + 10**(-14)) and (current_eval > eigenvalue - 10**(-14)):
            merged_evals[current_eval_index].append(eigenvalue)
        else:
            current_eval = eigenvalue
            current_eval_index += 1
            merged_evals.append([])
            merged_evals[current_eval_index].append(eigenvalue)
    return(merged_evals)

benzene = get_evals(cyclicpolyenehuckel(6))
print(merge_degeneracies(benzene))

def sorted_degeneracies(mer_deg):
    sort_deg = [[eigval[0],len(eigval)] for eigval in mer_deg]
    return(sort_deg)
print(sorted_degeneracies(merge_degeneracies(benzene)))

def arbitrary_to_numeral(arb_list):
    flatten_list = [node for yarn in arb_list for node in yarn]
    unique_flatten_list = list(set(flatten_list))
    numeral_list = [[unique_flatten_list.index(node) for node in yarn] for yarn in arb_list]
    return (numeral_list)

test_arb_list =[['a','ab',2,31,'a'],['ab',31]]
print(arbitrary_to_numeral(test_arb_list))
