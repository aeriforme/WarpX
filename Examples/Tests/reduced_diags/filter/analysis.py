#!/usr/bin/env python3
import re 
import numpy as np
import openpmd_api as io
import yt 
from scipy.constants import m_e, c
import matplotlib.pyplot as plt 


# function that extracts the bin centers from the first line of the histogram diagnostic 
def get_bins(fname):
    my_bins = []
    with open(fname) as f:
        for line in f: 
            if line.startswith('#'):
                data = line.split()
                for d in data[2:]:
                    m = re.search('=(.*)\(\)', d)
                    my_bins.append(m.group(1))
    return np.asarray(my_bins, dtype=float)

filename = 'diags/diag1000000'
ds = yt.load(filename)    
ad = ds.all_data()
    
px = ad['electrons', 'particle_momentum_x'].to_ndarray()      
py = ad['electrons', 'particle_momentum_y'].to_ndarray() 
pz = ad['electrons', 'particle_momentum_z'].to_ndarray()
w = ad['electrons', 'particle_weight'].to_ndarray()
x = ad['electrons', 'particle_position_x'].to_ndarray()

ux = px / (m_e*c)


fig, ax = plt.subplots(ncols=1, nrows=1, figsize=(1500./300., 1000./300.), dpi=300)

fname = 'diags/reducedfiles/histo_ele.txt'
mybins = get_bins(fname) # get bin centers 
red = np.loadtxt(fname)
ax.plot(mybins, red[0,2:], label='red diags - unfil') 

# get bin edges 
db =  mybins[1]-mybins[0]
mybins = mybins - 0.5*db 
mybins = np.append(mybins, mybins[-1]+0.5*db)

H, b = np.histogram(ux, bins=mybins, weights=w)
ax.plot(b[:-1],H, label='full diags - unfil', lw=2)

fname_fil = 'diags/reducedfiles/histo_ele_fil.txt'
mybins = get_bins(fname_fil)
red_fil = np.loadtxt(fname_fil)
ax.plot(mybins, red_fil[0,2:], label='red diags - fil') 

# get bin edges
db =  mybins[1]-mybins[0]
mybins = mybins - 0.5*db 
mybins = np.append(mybins, mybins[-1]+0.5*db)
H_fil, b_fil = np.histogram(ux[x>0], bins=mybins, weights=w[x>0])
ax.plot(b_fil[:-1],H_fil, label='full diags - fil', lw=2)

ax.legend()

plt.tight_layout()
plt.savefig('histo.png', dpi=300, bbox_inches='tight')
plt.close("all")
