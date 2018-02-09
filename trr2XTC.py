import glob
import os

trr = glob.glob('*.trr')

for t in trr:
    t_name = t.split('.')[0]
    os.system('gmx trjconv -f %s -o %s.xtc' % (t, t_name))
