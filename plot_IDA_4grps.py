# -*- coding: utf-8 -*-
"""
Created on Wed Jan 17 16:43:39 2018

@author: hliu
"""

IDA_plot = IDA.reindex(columns=IDA.columns.drop(['1','2','3'],2))
IDA_plot.columns = IDA_plot.columns.droplevel(2)
IDA_plot = IDA_plot.dropna()
IDA_plot.swaplevel(i=0,j=2, axis=1)['open'].max().max()
IDA_plot.swaplevel(i=0,j=2, axis=1)['open'].min().min()
IDA_plot.swaplevel(i=0,j=2, axis=1)['twist'].max().max()
IDA_plot.swaplevel(i=0,j=2, axis=1)['twist'].min().min()
xmin = 38.55
xmax = 54.39
ymin = 19.95
ymax = 70.32

ngrid = 50
x = IDA_plot['apo']['A']['open']
y = IDA_plot['apo']['A']['twist']
x_binNum, x_range = np.histogram(x, ngrid, range=(x.min(), x.max()))
y_binNum, y_range = np.histogram(y, ngrid, range=(y.min(), y.max()))
X,Y = np.meshgrid(x_range[:-1], y_range[:-1])
Z1, Z2 = np.meshgrid(x_binNum, y_binNum)
Z = Z1+Z2
plt.contour(X,Y,Z, colors='r')

def getBin(data, binNum, binRange=None):
    x = data[:, 0]
    y = data[:, 1]
    if binRange != None:
        x_range, y_range = binRange
        x_start = min(x_range); x_end = max(x_range)
        y_start = min(y_range); y_end = max(y_range)
        x_len = x_end - x_start
        y_len = y_end - y_start
    else:
        x_start = min(x); x_end = max(x)
        y_start = min(y); y_end = max(y)
        x_len = x_end - x_start
        y_len = y_end - y_start
    n_Xgrid, n_Ygrid = binNum
    Bins = np.zeros((n_Xgrid, n_Ygrid))


##########################################
# Find out structures near local minimum#
########################################

IDA_idx = IDA.reindex(columns=IDA.columns.drop(['1','2','3'],2))
IDA_idx.columns = IDA_idx.columns.droplevel(2)

p = sns.kdeplot(IDA_idx.apo.B.twist.dropna(), shade=True);x,y = p.get_lines()[0].get_data()
x[y.argmax()]
from scipy.signal import find_peaks_cwt
indexes = find_peaks_cwt(y, np.arange(1,10))


p = sns.kdeplot(IDA_idx.holo.A.twist.dropna(), shade=True);x,y = p.get_lines()[0].get_data()
x[y.argmax()]
from scipy.signal import find_peaks_cwt
indexes = find_peaks_cwt(y, np.arange(1,10))

p = sns.kdeplot(IDA_idx.holo.B.twist.dropna(), shade=True)
x,y = p.get_lines()[0].get_data()
x[y.argmax()]
from scipy.signal import find_peaks_cwt
indexes = find_peaks_cwt(y, np.arange(1,10))
# apo B first peak idx 48, value 36.5, traj idx 1004; second peak idx 90, value 55.7, traj idx 6782
# holo A first peak idx 52, value 41.77, traj idx 3899; second peak idx 90, value 54.05, traj idx 739
# holo B first peak idx 62, value 50.96; second peak idx 78, value 56.55