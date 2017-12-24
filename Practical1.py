import huckellib as hl

print(hl.gen_lin_yarn(6))
print(hl.gen_cyc_yarn(6))

print(hl.cyclicpolyenehuckel(6))

print(hl.poly_ene_huckel(6))

yarn_list = [list(yarn) for yarn in ['abcda','bd']]
print(yarn_list)

eig = hl.connectivity_to_eigenvalues(yarn_list)

print(eig)
print('-----------------------------------------------------')
benzene_eig = hl.connectivity_to_eigenvalues([list('abcdefa')])
for eig in benzene_eig:
    print(eig)


print('-----------------------------------------------------')
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

for eigen in hl.connectivity_to_eigenvalues(bucky_connectivity_list):
    print(eigen)
