import sys
import pandas as pd

prefix = sys.argv[1].split('.')[0]
lines = open(sys.argv[1], 'r').readlines()

def cleanLine(l):
    l = l.replace('-Side', '')
    l = l.replace('-Main', '')
    l = l.replace('%', '')
    l = [i for i in l.split() if i != sys.argv[2]]
    return l

new_lines = []
for l in lines[2:]:
    new_lines.append(cleanLine(l))

df = pd.DataFrame(new_lines)
df.columns = ['ResInfo', 'Percent']
df.Percent = df.Percent.astype('float')

result = df.groupby('ResInfo')['Percent'].apply(sum)
result.sort(ascending=False)

result_l = ['{:<8}{:5.2f}'.format(index, result.ix[index]) for index in result.index]
s = '\n'.join(result_l)
out = open(prefix+'_sum.dat','w')
out.write(s)
