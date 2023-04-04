import os 
import sys 
import re
import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import m_e, c, hbar, e 
import openpmd_api as io
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import LogNorm, Normalize
from matplotlib import use, cm
import matplotlib.colors

E_crit = m_e**2*c**3/(e*hbar)
B_crit = m_e**2*c**2/(e*hbar)

# extract numbers from a string 
def find_num_in_line(line):
    items = re.findall('-?\ *[0-9]+\.?[0-9]*(?:[Ee]\ *-?\ *[0-9]+)?', line)
    fitems = [float(it) for it in items]
    if len(fitems)==1:
        return fitems[0]
    else:
        return fitems

# get input parameters from warpx_used_inputs 
with open('./warpx_used_inputs', 'rt') as f:
    lines = f.readlines()
    for line in lines:
        if 'my_constants.nx' in line:
            nx = find_num_in_line(line)
        if 'my_constants.ny' in line:
            ny = find_num_in_line(line)
        if 'my_constants.nz' in line:
            nz = find_num_in_line(line)
        if 'my_constants.Lx' in line:
            Lx = find_num_in_line(line)
        if 'my_constants.Ly' in line:
            Ly = find_num_in_line(line)
        if 'my_constants.Lz' in line:
            Lz = find_num_in_line(line)                        
        if 'beam_e.density' in line:
            n1 = find_num_in_line(line)   
        if 'beam_p.density' in line:
            n2 = find_num_in_line(line)   
        if 'beam_e.ux' in line:
            u1x = find_num_in_line(line) 
        if 'beam_e.uy' in line:
            u1y = find_num_in_line(line) 
        if 'beam_e.uz' in line:
            u1z = find_num_in_line(line) 
        if 'beam_p.ux' in line:
            u2x = find_num_in_line(line) 
        if 'beam_p.uy' in line:
            u2y = find_num_in_line(line) 
        if 'beam_p.uz' in line:
            u2z = find_num_in_line(line) 
        if 'beam_e.num_particles_per_cell_each_dim' in line:
            aux1, aux2, aux3 = find_num_in_line(line)      
            Ne = aux1*aux2*aux3
        if 'beam_p.num_particles_per_cell_each_dim' in line:
            aux1, aux2, aux3 = find_num_in_line(line)      
            Np = aux1*aux2*aux3
        if 'particles.E_external_particle' in line:
            Ex, Ey, Ez = find_num_in_line(line)      
        if 'particles.B_external_particle' in line:
            Bx, By, Bz = find_num_in_line(line)      
            
dV = Lx/nx * Ly/ny * Lz/nz

    
def chi(ux, uy, uz, Ex=Ex, Ey=Ey, Ez=Ez, Bx=Bx, By=By, Bz=Bz):
    #print(ux, uy, uz, Ex, Ey, Ez, Bx, By, Bz)
    gamma = np.sqrt(1.+ux**2+uy**2+uz**2)
    vx = ux / gamma * c
    vy = uy / gamma * c
    vz = uz / gamma * c
    tmp1x = Ex + vy*Bz - vz*By 
    tmp1y = Ey - vx*Bz + vz*Bx
    tmp1z = Ez + vx*By - vy*Bx
    tmp2 = (Ex*vx + Ey*vy + Ez*vz)/c
    
    chi = gamma/E_crit*np.sqrt(tmp1x**2+tmp1y**2+tmp1z**2 - tmp2**2)
    return chi 


def luminosity():
    series = io.Series("diags/diag1/openpmd_%T.bp",io.Access.read_only) 
    iterations = np.asarray(series.iterations)
    lumi = []
    for n,ts in enumerate(iterations):
        it = series.iterations[ts]
        rho1 = it.meshes["rho_beam_e"]
        dV = np.prod(rho1.grid_spacing)
        rho1 = it.meshes["rho_beam_e"][io.Mesh_Record_Component.SCALAR].load_chunk()
        rho2 = it.meshes["rho_beam_p"][io.Mesh_Record_Component.SCALAR].load_chunk()
        q1 = np.unique(it.particles["beam_e"]["charge"][io.Mesh_Record_Component.SCALAR].load_chunk())
        q2 = np.unique(it.particles["beam_p"]["charge"][io.Mesh_Record_Component.SCALAR].load_chunk())
        series.flush()
        n1 = rho1/q1
        n2 = rho2/q2
        #print(np.where((n1>0), n1, 0.))
        print(n1)
        l = 2*np.sum(rho1/q1*rho2/q2)*dV
        print('llllllllllllllllllllll ', l, np.sum(n1*n2), np.sum(n1*n1),np.sum(n2*n2))
        lumi.append(l)
    return lumi



print('chi of electrons ----------------------------------------')
print('theory:', chi(u1x, u1y, u1z))

fname='diags/reducedfiles/ParticleExtrema_beam_e.txt'
chimax = np.loadtxt(fname)[:,19]
print('ParticleExtrema diag = ', chimax)

fname='diags/reducedfiles/ParticleExtrema_beam_e.txt'
chimin = np.loadtxt(fname)[:,18]
print('ParticleExtrema diag = ', chimin)

fname='diags/reducedfiles/ColliderRelevant_beam_e_beam_p.txt'
chimin = np.loadtxt(fname)[:,6]
print('ColliderRelevant diag = ', chimin)

fname='diags/reducedfiles/ColliderRelevant_beam_e_beam_p.txt'
chiave = np.loadtxt(fname)[:,7]
print('ColliderRelevant diag = ', chiave)


print('luminosity')


print('from PIC data', luminosity())

fname='diags/reducedfiles/ColliderRelevant_beam_e_beam_p.txt'
lum = np.loadtxt(fname)[:,2]
print('ColliderRelevant diag = ', lum)


print('theory from input', 2.*n1*n2*Lx*Ly*Lz)













