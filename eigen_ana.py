#!/home/hliu/anaconda2/bin/python
import sys
import argparse
from researchcode.plotting.plot_set import *

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument(dest='file')
    parser.add_argument('--plot',
                        dest='plot',
                        action='store_true')
    parser.add_argument('--stat',
                        dest='stat',
                        action='store_true')
    parser.add_argument('--index',
                        dest='index')
    args = parser.parse_args()
    f = open(args.file, 'r')
    lines = [l for l in f.readlines() if l.startswith(' ')]
    lines_num = [float(l.strip().split()[1]) for l in lines]
    eigen_vals = np.array(lines_num)

    if args.plot:
        plt.plot(range(len(eigen_vals)), eigen_vals)
        plt.xlabel('Index of Eigenvalues')
        plt.ylabel('Variance')
        plt.tight_layout()
        plt.savefig('eignval.png')
    
    if args.stat:
        print(eigen_vals[:int(int(args.index))].sum()/eigen_vals.sum())
