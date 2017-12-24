import huckellib as hl
import sys
import getopt


def main(argv):
    eigenvalues = []
    fout = ''
    output = ''
    try:
        opts, args = getopt.getopt(argv, "hl:c:f:")
    except getopt.GetoptError:
        sys.exit('huckel.py -l <linear length> -c <cyclic length> -f <filename>')
    if opts == []:
        print('huckel.py -l <linear length> -c <cyclic length> -f <filename>')
        sys.exit()
    opt, arg = opts[0]

    if opt == '-h':
        print('huckel.py -l <linear length> -c <cyclic length> -f <filename>')
        sys.exit()
    elif opt == '-l':  # linear polyene
        try:
            eigenvalues = hl.poly_ene_huckel(int(arg))
        except TypeError as err:
            sys.exit(err)
        except ValueError as err:
            sys.exit(err)
        fout = 'out_linear_polyene_c' + arg + '.txt'
    elif opt == '-c':  # cyclic polyene
        try:
            eigenvalues = hl.cyclicpolyenehuckel(int(arg))
        except TypeError as err:
            sys.exit(err)
        except ValueError as err:
            sys.exit(err)
        fout = 'out_cyclic_polyene_c' + arg + '.txt'
    elif opt == '-f':  # calculate based on a given file with a connection map
        try:
            with open(arg, 'r') as f:
                yarns = f.read()
                yarns = ''.join(yarns.split())
                yarn_list = [yarn.split(',') for yarn in yarns.split(';')]
                eigenvalues = hl.connectivity_to_eigenvalues(yarn_list)
        except FileNotFoundError:
            sys.exit('Please provide an existing file')
        fout = 'out_' + arg
    else:
        print('huckel.py -l <linear length> -c <cyclic length> -f <filename>')
        sys.exit()

    output += '{:^{width}} {} {:^{width}}'.format(
        'Eigenvalue', ',', 'Degeneracy', width=15) + '\n'
    for eig in eigenvalues:
        output += '{e[0].real:>{width}.{prec}f} {delim} {e[1]:>{width}d}'.format(
            e=eig, delim=',', width=15, prec=9) + '\n'
    with open(fout, 'w') as file:
        print(output)
        file.write(output)


if __name__ == "__main__":
    main(sys.argv[1:])
