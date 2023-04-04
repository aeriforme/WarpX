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
        if 'warpx.cfl' in line:
            cfl = find_num_in_line(line)
        if 'max_step' in line:
            num_steps = find_num_in_line(line)
        if 'beam_e.multiple_particles_ux' in line:
            uex1, uex2, uex3 = find_num_in_line(line)
        if 'beam_e.multiple_particles_uy' in line:
            uey1, uey2, uey3 = find_num_in_line(line)
        if 'beam_e.multiple_particles_uz' in line:
            uez1, uez2, uez3 = find_num_in_line(line)   
        if 'beam_e.multiple_particles_weight' in line:
            we1, we2, we3 = find_num_in_line(line)                                
        if 'beam_p.multiple_particles_ux' in line:
            upx1, upx2, upx3 = find_num_in_line(line)
        if 'beam_p.multiple_particles_uy' in line:
            upy1, upy2, upy3 = find_num_in_line(line)
        if 'beam_p.multiple_particles_uz' in line:
            upz1, upz2, upz3 = find_num_in_line(line)       
        if 'beam_p.multiple_particles_weight' in line:
            wp1, wp2, wp3 = find_num_in_line(line)                      
        if 'particles.E_external_particle ' in line:
            Ex, Ey, Ez = find_num_in_line(line)      
        if 'particles.B_external_particle ' in line:
            Bx, By, Bz = find_num_in_line(line)      

pex1, pex2, pex3 = uex1*m_e*c, uex2*m_e*c, uex3*m_e*c
pey1, pey2, pey3 = uey1*m_e*c, uey2*m_e*c, uey3*m_e*c
pez1, pez2, pez3 = uez1*m_e*c, uez2*m_e*c, uez3*m_e*c    

ppx1, ppx2, ppx3 = upx1*m_e*c, upx2*m_e*c, upx3*m_e*c
ppy1, ppy2, ppy3 = upy1*m_e*c, upy2*m_e*c, upy3*m_e*c
ppz1, ppz2, ppz3 = upz1*m_e*c, upz2*m_e*c, upz3*m_e*c    

ge1, ge2, ge3 = np.sqrt(1.+uex1**2+uey1**2+uez1**2), np.sqrt(1.+uex2**2+uey2**2+uez2**2), np.sqrt(1.+uex3**2+uey3**2+uez3**2) 
gp1, gp2, gp3 = np.sqrt(1.+upx1**2+upy1**2+upz1**2), np.sqrt(1.+upx2**2+upy2**2+upz2**2), np.sqrt(1.+upx3**2+upy3**2+upz3**2) 

Ee1, Ee2, Ee3 = m_e*c**2*(ge1-1.), m_e*c**2*(ge2-1.), m_e*c**2*(ge3-1.)
Ep1, Ep2, Ep3 = m_e*c**2*(gp1-1.), m_e*c**2*(gp2-1.), m_e*c**2*(gp3-1.)

    
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
        print(np.where((n1>0), n1, 0.))
        
        l = 2*np.sum(rho1/q1*rho2/q2)*dV
        print('llllllllllllllllllllll ', l, np.sum(n1*n2), np.sum(n1*n1),np.sum(n2*n2))
        lumi.append(l)
    return lumi



print('chi of electrons anal = ', chi(uex1, uey1, uez1), chi(uex2, uey2, uez2), chi(uex3, uey3, uez3))
print('chi max of electrons ----------------------------------------')
print('theory:', np.max( (chi(uex1, uey1, uez1), chi(uex2, uey2, uez2), chi(uex3, uey3, uez3))) )
fname='diags/reducedfiles/ParticleExtrema_beam_e.txt'
chimax = np.loadtxt(fname)[:,19]
print('ParticleExtrema diag = ', chimax)
fname='diags/reducedfiles/ColliderRelevant_beam_e_beam_p.txt'
chimax = np.loadtxt(fname)[:,8]
print('ColliderRelevant diag = ', chimax)


print('chi min of electrons ----------------------------------------')
print('theory:', np.min( (chi(uex1, uey1, uez1), chi(uex2, uey2, uez2), chi(uex3, uey3, uez3))))
fname='diags/reducedfiles/ParticleExtrema_beam_e.txt'
chimin = np.loadtxt(fname)[:,18]
print('ParticleExtrema diag = ', chimin)
fname='diags/reducedfiles/ColliderRelevant_beam_e_beam_p.txt'
chimin = np.loadtxt(fname)[:,6]
print('ColliderRelevant diag = ', chimin)


print('chi ave of electrons ----------------------------------------')
print('theory:', np.average( (chi(uex1, uey1, uez1), chi(uex2, uey2, uez2), chi(uex3, uey3, uez3)), weights = (we1, we2, we3) ))
fname='diags/reducedfiles/ColliderRelevant_beam_e_beam_p.txt'
chiave = np.loadtxt(fname)[:,7]
print('ColliderRelevant diag = ', chiave)



print('luminosityyyyy')
print('theory', luminosity())
fname='diags/reducedfiles/ColliderRelevant_beam_e_beam_p.txt'
lum = np.loadtxt(fname)[:,2]
print('ColliderRelevant diag = ', lum)
















