
#!/usr/bin/python

# Standard Python modules
from pyns.standard import *

import matplotlib.pylab as pylab
import numpy as np
import matplotlib
from scipy.optimize import fsolve

# PyNS modules
from pyns.constants          import *
from pyns.operators          import *
from pyns.discretization     import *
from pyns.display            import plot, write
from pyns.physical           import properties
from pyns.physical.constants import G

# membrane aux functions
from p_v_sat import *
# from latent_heat import *

plt.close('all')

#==========================================================================
#
# Define problem
#
#==========================================================================

AIR = 0
H2O = 1
FIL = 2

u_in   = 0.12 # m/s
t_in   = 80   # °C
a_salt = 90.0 # g/l
t_cold = 20   # °C

# Node coordinates for both domains
xn = (nodes(0,    0.1, 150), nodes(0,    0.1, 150), nodes(0,           0.1, 150))
yn = (nodes(-0.009, 0, 30),  nodes(0,  0.004, 12),  nodes(-0.01,    -0.009,  3))
zn = (nodes(0,    0.05, 75), nodes(0,    0.05, 75), nodes(0,           0.05, 75))

# Cell coordinates 
xc = (avg(xn[AIR]), avg(xn[H2O]), avg(xn[FIL]))
yc = (avg(yn[AIR]), avg(yn[H2O]), avg(yn[FIL]))
zc = (avg(zn[AIR]), avg(zn[H2O]), avg(zn[FIL]))

# Cell dimensions
cell = [cartesian_grid(xn[AIR],yn[AIR],zn[AIR]),  \
        cartesian_grid(xn[H2O],yn[H2O],zn[H2O]),  \
        cartesian_grid(xn[FIL],yn[FIL],zn[FIL])]

nx,ny,nz, dx,dy,dz = (cell[AIR][0], cell[H2O][0], cell[FIL][0]),  \
                     (cell[AIR][1], cell[H2O][1], cell[FIL][1]),  \
                     (cell[AIR][2], cell[H2O][2], cell[FIL][2]),  \
                     (cell[AIR][3], cell[H2O][3], cell[FIL][3]),  \
                     (cell[AIR][4], cell[H2O][4], cell[FIL][4]),  \
                     (cell[AIR][5], cell[H2O][5], cell[FIL][5])
rc,ru,rv,rw =        (cell[AIR][6], cell[H2O][6], cell[FIL][6]),  \
                     (cell[AIR][7], cell[H2O][7], cell[FIL][7]),  \
                     (cell[AIR][8], cell[H2O][8], cell[FIL][8]),  \
                     (cell[AIR][9], cell[H2O][9], cell[FIL][9])

# Set physical properties					 
prop = [properties.air(rc[AIR]),    \
        properties.water(rc[H2O]),  \
        properties.water(rc[FIL])]
					 
# Set physical properties temperature dependent:
#prop = [properties.air(round((t_in+t_cold)/2,-1),rc[AIR]), \
#        properties.water(t_in,rc[H2O]),                    \
#        properties.water(t_cold,rc[FIL])]

rho, mu, cap, kappa = (prop[AIR][0], prop[H2O][0], prop[FIL][0]),  \
                      (prop[AIR][1], prop[H2O][1], prop[FIL][1]),  \
                      (prop[AIR][2], prop[H2O][2], prop[FIL][2]),  \
                      (prop[AIR][3], prop[H2O][3], prop[FIL][3])
                      
diff = (ones(rc[AIR])*4.0E-4, ones(rc[H2O])*1.99E-09)

#h_d = [ 0.0, latent_heat(t_in), latent_heat(t_cold) ]
h_d = [ 0.0, 2359E3, 2359E3]

M_H2O  = 18E-3      # kg/mol
M_AIR  = 28E-3      # kg/mol
M_salt = 58.4428E-3 # kg/mol      

R = 8.314
pi = 3.1415

# Membrane properties
membrane = namedtuple('membrane', 'd kap eps tau r p t pv j t_old')
  # d is thickness in m, kap is thermal conductivity in W/mK
  # eps is porosity, tau is tortuosity
  # r is pore radius

mem = membrane(166E-6,   \
                 0.2,    \
                 0.687,  \
                 2,      \
                 0.1E-6, \
                 zeros((nx[AIR],1,nz[AIR])), \
                 zeros((nx[AIR],1,nz[AIR])), \
                 zeros((nx[AIR],1,nz[AIR])), \
                 zeros((nx[AIR],1,nz[AIR])), \
                 zeros((nx[AIR],1,nz[AIR])))
 
# Create unknowns; names, positions and sizes
uf    = [Unknown('face-u-vel',    X, ru[AIR], DIRICHLET),  \
         Unknown('face-u-vel',    X, ru[H2O], DIRICHLET),  \
         Unknown('face-u-vel',    X, ru[FIL], DIRICHLET)]
vf    = [Unknown('face-v-vel',    Y, rv[AIR], DIRICHLET),  \
         Unknown('face-v-vel',    Y, rv[H2O], DIRICHLET),  \
         Unknown('face-v-vel',    Y, rv[FIL], DIRICHLET)]
wf    = [Unknown('face-w-vel',    Z, rw[AIR], DIRICHLET),  \
         Unknown('face-w-vel',    Z, rw[H2O], DIRICHLET),  \
         Unknown('face-w-vel',    Z, rw[FIL], DIRICHLET)]
p     = [Unknown('pressure',      C, rc[AIR], NEUMANN),  \
         Unknown('pressure',      C, rc[H2O], NEUMANN),  \
         Unknown('pressure',      C, rc[FIL], NEUMANN)]
t     = [Unknown('temperature',   C, rc[AIR], NEUMANN),  \
         Unknown('temperature',   C, rc[H2O], NEUMANN),  \
         Unknown('temperature',   C, rc[FIL], NEUMANN)]
a     = [Unknown('concentration', C, rc[AIR], NEUMANN),  \
         Unknown('concentration', C, rc[H2O], NEUMANN),  \
         Unknown('concentration', C, rc[FIL], NEUMANN)]
p_tot = [Unknown('pressure-tot',  C, rc[AIR], NEUMANN),  \
         Unknown('pressure-tot',  C, rc[H2O], NEUMANN)]

# just for air
p_v =[Unknown('vapor_pressure',C, rc[AIR], NEUMANN)]
M  = [Unknown('molar mass',    C, rc[AIR], NEUMANN)]

# Specify boundary conditions

for k in range(0,nz[H2O]):
  uf[H2O].bnd[W].val[:1,:,k] = par(u_in, yn[H2O])
uf[H2O].bnd[E].typ[:1,:,:] = OUTLET 
uf[H2O].bnd[E].val[:1,:,:] = u_in

for c in range(AIR,H2O):
  for j in (B,T):
    uf[c].bnd[j].typ[:] = NEUMANN     
    vf[c].bnd[j].typ[:] = NEUMANN     
    wf[c].bnd[j].typ[:] = NEUMANN     
  
t[AIR].bnd[S].typ[:,:1,:] = DIRICHLET  
t[AIR].bnd[S].val[:,:1,:] = t_cold
t[AIR].bnd[N].typ[:,:1,:] = DIRICHLET  
t[AIR].bnd[N].val[:,:1,:] = 70

t[H2O].bnd[W].typ[:1,:,:] = DIRICHLET
t[H2O].bnd[W].val[:1,:,:] = t_in
t[H2O].bnd[S].typ[:,:1,:] = DIRICHLET  
t[H2O].bnd[S].val[:,:1,:] = 70
 
t[FIL].bnd[S].typ[:,:1,:] = DIRICHLET
t[FIL].bnd[S].val[:,:1,:] = t_cold
t[FIL].bnd[N].typ[:,:1,:] = DIRICHLET  
t[FIL].bnd[N].val[:,:1,:] = t_cold

t_int_mem = t[H2O].bnd[S].val[:,:1,:] # temporary

a[H2O].bnd[W].typ[:1,:,:] = DIRICHLET
a[H2O].bnd[W].val[:1,:,:] = a_salt/rho[H2O][:1,:,:]

M[AIR].bnd[S].typ[:,:1,:] = DIRICHLET
M[AIR].bnd[S].val[:,:1,:] = M[AIR].val[:,:1,:]

p_v[AIR].bnd[S].typ[:,:,:] = DIRICHLET
p_v[AIR].bnd[S].val[:,:,:] = p_v_sat(t[FIL].bnd[N].val[:,:,:])
p_v[AIR].bnd[N].typ[:,:,:] = DIRICHLET
p_v[AIR].bnd[N].val[:,:,:] = p_v_sat(t[H2O].bnd[S].val[:,:,:])
p_v[AIR].bnd[S].typ[:,:,:] = DIRICHLET
p_v[AIR].bnd[S].val[:,:,:] = p_v_sat(t[FIL].bnd[N].val[:,:,:])

t[AIR].val[:,:,:] = round((t_in+t_cold)/2,-1)
t[H2O].val[:,:,:] = 70
t[FIL].val[:,:,:] = t_cold

a[AIR].val[:,:,:] = p_v_sat(t[AIR].val[:,:,:])*1E-5*M_H2O/M_AIR
M[AIR].val[:,:,:] = 1/((1-a[AIR].val[:,:,:])/M_AIR + a[AIR].val[:,:,:]/M_H2O)
a[H2O].val[:,:,:] = a_salt/rho[H2O][:,:,:]
 
for c in range(AIR,FIL):
  adj_n_bnds(p[c])
  adj_n_bnds(t[c])
  adj_n_bnds(a[c])
  
# max values in domain:
t_max = 0.0
t_min = 100.0

for c in range(W,T):
  t_max = np.amax([t_max, np.amax(t[H2O].bnd[c].val)])
  t_min = np.amin([t_min, np.amin(t[FIL].bnd[c].val)]) 
  
  # Time-stepping parameters
dt  =    0.0005  # time step
ndt =    1       # number of time steps
dt_plot = ndt    # plot frequency

obst = [zeros(rc[AIR]), zeros(rc[H2O]),zeros(rc[FIL])]

#==========================================================================
#
# Solution algorithm
#
#==========================================================================

#-----------
#
# Time loop 
#
#-----------
for ts in range(1,ndt+1):

  write.time_step(ts)
 
  #------------------
  # Store old values
  #------------------
  for c in range(AIR,FIL):
    a[c].old[:]  = a[c].val
    t[c].old[:]  = t[c].val
    uf[c].old[:] = uf[c].val
    vf[c].old[:] = vf[c].val
    wf[c].old[:] = wf[c].val
  mem.t_old[:] = mem.t
  t_int_mem_old = t_int_mem
  
  #calculate rho
  t_min_rho = 20
  t_max_rho = 80
  rho_min = 1.205
  rho_max = 1.000
  rho[AIR][:,:,:] = (t[AIR].val[:,:,:] - t_min_rho)/(t_max_rho - t_min_rho) * \
                     (rho_max - rho_min) + rho_min
    
  #------------------------
  # Partial vapor pressure & molar mass
  #------------------------  
  
  # AIR domain values  
  M[AIR].val[:,:,:] = 1/((1-a[AIR].val[:,:,:])/M_AIR + a[AIR].val[:,:,:]/M_H2O)  
  p_v[AIR].val[:,:,:] = a[AIR].val[:,:,:] *M[AIR].val[:,:,:]/M_H2O * (p_tot[AIR].val[:,:,:] +1E5) 
  
  # Liquid film & AIR S bnd values 
  M[AIR].bnd[S].val[:,:,:] = np.power(((1-a[AIR].bnd[S].val[:,:,:])/M_AIR + a[AIR].bnd[S].val[:,:,:]/M_H2O),(-1))
  p_v[AIR].bnd[S].val[:,:,:] = a[AIR].bnd[S].val[:,:,:] *M[AIR].bnd[S].val[:,:,:]/M_H2O * (p_tot[AIR].val[:,:1,:] +1E5)
  
  t_int_film = t_sat(p_v[AIR].bnd[S].val[:,:,:])
  print('t_int_film = ' + '%3.4f' %np.mean(t_int_film))
  
  t[AIR].bnd[S].val[:,:,:] = t_int_film
  t[FIL].bnd[N].val[:,:,:] = t_int_film
  
  m_out = (2* kappa[AIR][:,:1,:] / dy[AIR][:,:1,:] *    \
                   (t_int_film - t[AIR].val[:,:1,:])    \
          + 2*kappa[FIL][:,-1:,:] / dy[FIL][:,-1:,:] *  \
                  (t_int_film - t[FIL].val[:,-1:,:]))   \
          * dx[AIR][:,:1,:] * dz[AIR][:,:1,:] / h_d[FIL]  # kg/s
  
  print('m_out = ' + '%3.4e' %np.mean(m_out))  
    
  #------------------------
  # Membrane
  #------------------------   
  
  # Compute new temperature of the membrane  
  
  kappa_mem = mem.kap*(1-mem.eps) + kappa[AIR][:,-1:,:]*mem.eps 
  const_mem_1 = kappa_mem *dy[AIR][:,-1:,:] /kappa[AIR][:,-1:,:] /mem.d
  mem.t[:,:1,:] = ((1.0+const_mem_1)*t_int_mem + t[AIR].val[:,-1:,:]) / (2.0+const_mem_1)
  
  mem.t[:,:1,:] = mem.t + 273.15;
  mem.p[:,:1,:] = (p_tot[H2O].val[:,:1,:] + p_tot[AIR].val[:,-1:,:]) /2.0 + 1E5
  mem.pv[:,:1,:] = (p_v[AIR].bnd[N].val[:,:1,:] + p_v[AIR].val[:,-1:,:]) /2.0
  
  # Diffusion Coefficients
  C_K = 2.0*mem.eps*mem.r/(3.0*mem.tau*mem.d)*np.power(8.0*M_H2O/(mem.t*R*pi),0.5)
  C_M = mem.eps*mem.p*diff[AIR][:,-1:,:]/(mem.tau*R*mem.t*(mem.p-mem.pv))
  C_T = 1.0/(1.0/C_K + 1.0/C_M)
  
  # Jump condition at membrane
  lhs_lin_mem = (2.0*kappa[H2O][:,:1,:]/dy[H2O][:,:1,:] \
     + 1.0/(dy[AIR][:,-1:,:]/(2.0*kappa[AIR][:,-1:,:]) + mem.d/kappa_mem)) \
     * mem.eps * dx[AIR][:,-1:,:] * dz[AIR][:,-1:,:] / h_d[H2O]
     
  lhs_fun_mem = C_T*dx[AIR][:,-1:,:]*dz[AIR][:,-1:,:]
  
  rhs_mem = C_T*dx[AIR][:,-1:,:]*dz[AIR][:,-1:,:]*p_v[AIR].val[:,-1:,:] \
    + mem.eps*dx[AIR][:,-1:,:]*dz[AIR][:,-1:,:]/h_d[H2O] \
    * ( 2.0*kappa[H2O][:,:1,:]/dy[H2O][:,:1,:]*t[H2O].val[:,:1,:] \
      + 1.0/(dy[AIR][:,-1:,:]/(2.0*kappa[AIR][:,-1:,:])+mem.d/kappa_mem) \
      * t[AIR].val[:,-1:,:])
  
  for ii in range(0,nx[AIR]):
    for kk in range(0,nz[AIR]):
      jump_cond_mem = lambda t: lhs_lin_mem[ii,:1,kk]*t + lhs_fun_mem[ii,:1,kk] *\
        p_v_sat_salt(t, a[H2O].val[ii,:1,kk], M_salt) - rhs_mem[ii,:1,kk]
      t_int_mem[ii,:1,kk] = fsolve(jump_cond_mem, t[H2O].val[ii,-1:,kk])
  
  print('t_int_mem = ' + '%3.4f' %np.mean(t_int_mem))
  
  t[H2O].bnd[S].val[:,:1,:] = t_int_mem
  const_mem_2 = 2*kappa[AIR][:,-1:,:]*mem.d/kappa_mem/dy[AIR][:,-1:,:];
  t[AIR].bnd[N].val[:,:1,:] = (t_int_mem + const_mem_2 *t[AIR].val[:,-1:,:])/(1+const_mem_2)
  p_v[AIR].bnd[N].val[:,:,:]= p_v_sat_salt(t_int_mem, a[H2O].val[:,:1,:], M_salt)
  mem.j[:,:,:] = C_T[:,:,:] *dx[AIR][:,-1:,:]*dz[AIR][:,-1:,:]*(p_v[AIR].bnd[N].val[:,:,:]-p_v[AIR].val[:,-1:,:])  
  
  #------------------------
  # Concentration
  #------------------------
  
  q_a = [zeros(rc[AIR]),
         zeros(rc[H2O])]
  dv  = [dx[AIR]*dy[AIR]*dz[AIR], 
         dx[H2O]*dy[H2O]*dz[H2O]]  
  q_a[AIR][:,-1:,:] = mem.j [:,:1,:] / dv[AIR][:,-1:,:]
  q_a[AIR][:,:1,:]  = - m_out[:,:1,:] / dv[AIR][:,:1,:] 
  
  # in case of v[H2O].bnd[S].val ~= 0 correct convection into membrane 
  for c in range(AIR,H2O):
    calc_t(a[c], (uf[c],vf[c],wf[c]), rho[c], diff[c],  \
           dt, (dx[c],dy[c],dz[c]), obst[c],
           source = q_a[c])

  #------------------------
  # Temperature (enthalpy)
  #------------------------

  q_t = [zeros(rc[AIR]),
         zeros(rc[H2O]), 
         zeros(rc[FIL])]
  dv  = [dx[AIR]*dy[AIR]*dz[AIR], 
         dx[H2O]*dy[H2O]*dz[H2O], 
         dx[FIL]*dy[FIL]*dz[FIL]]  
  q_t[H2O][:,:1,:]  = -h_d[H2O]*mem.j [:,:1,:] / dv[H2O][:,:1,:]
  q_t[FIL][:,-1:,:] =  h_d[FIL]*m_out[:,:1,:] / dv[FIL][:,-1:,:]
  
  for c in range(AIR,FIL):
    calc_t(t[c], (uf[c],vf[c],wf[c]), (rho[c]*cap[c]), kappa[c],  \
           dt, (dx[c],dy[c],dz[c]), obst[c],
           source = q_t[c])

  for c in range(AIR,FIL):
    t[c].val[t[c].val > t_max] = t_max
    t[c].val[t[c].val < t_min] = t_min

  #-----------------------
  # Momentum conservation
  #-----------------------
  for c in (AIR,H2O):
    g_v = -G * avg(Y, rho[c])
  
    ef = zeros(ru[c]), g_v, zeros(rw[c])
    
    calc_uvw((uf[c],vf[c],wf[c]), (uf[c],vf[c],wf[c]), rho[c], mu[c],  \
             dt, (dx[c],dy[c],dz[c]), obst[c],
             pressure = p_tot[c],
             force = ef)
  
  #----------
  # Pressure
  #----------
  for c in (AIR,H2O):
    calc_p(p[c], (uf[c],vf[c],wf[c]), rho[c],  \
           dt, (dx[c],dy[c],dz[c]), obst[c])
  
    p_tot[c].val = p_tot[c].val + p[c].val
  
  #---------------------
  # Velocity correction
  #---------------------
  for c in (AIR,H2O):
    corr_uvw((uf[c],vf[c],wf[c]), p[c], rho[c],  \
             dt, (dx[c],dy[c],dz[c]), obst[c])
 
  # Compute volume balance for checking 
  for c in (AIR,H2O):
    err = vol_balance((uf[c],vf[c],wf[c]),  \
                      (dx[c],dy[c],dz[c]), obst[c])
    print('Maximum volume error after correction: %12.5e' % abs(err).max())

  # Check the CFL number too 
  for c in (AIR,H2O):
    cfl = cfl_max((uf[c],vf[c],wf[c]), dt, (dx[c],dy[c],dz[c]))
    print('Maximum CFL number: %12.5e' % cfl)

#==========================================================================
#
# Visualisation
#
#==========================================================================
#%%
  if ts % dt_plot == 0:
    plt.close('all')
    
    z_pos = 10
    
    xc = avg(xn[AIR])
    yc = np.append(avg(yn[FIL]), avg(yn[AIR]),axis=0)
    yc = np.append(yc, avg(yn[H2O]),axis=0)
    
    t_plot=np.append(t[FIL].val[:,:,z_pos],t[AIR].val[:,:,z_pos],axis=1)
    t_plot=np.append(t_plot, t[H2O].val[:,:,z_pos],axis=1)
    t_plot=transpose(t_plot)
    p_plot=np.append(p[FIL].val[:,:,z_pos],p[AIR].val[:,:,z_pos],axis=1)
    p_plot=np.append(p_plot, p[H2O].val[:,:,z_pos],axis=1)
    p_plot=transpose(p_plot)
    a_plot=np.append(a[FIL].val[:,:,z_pos],a[AIR].val[:,:,z_pos],axis=1)
    a_plot=np.append(a_plot, a[H2O].val[:,:,z_pos],axis=1)
    a_plot=transpose(a_plot)
    
    plt.figure
    plt.subplot(2,2,1)
    levels_t=linspace( t_plot.min(), t_plot.max(), 11)
    norm_t=cm.colors.Normalize( vmax=t_plot.max(), vmin=t_plot.min() )
    cax_t=plt.contourf(xc,yc,t_plot,levels_t,cmap='rainbow',norm=norm_t)
    cbar_t=plt.colorbar(cax_t)
    plt.title('Temperature')
    plt.xlabel('x [m]')
    #plt.ylim([-1E1,1E1])
    plt.ylabel('y [m]' )
    
    plt.subplot(2,2,2)
    matplotlib.rcParams['contour.negative_linestyle'] = 'solid'
    cax_p=plt.contourf(xc,yc,p_plot,cmap='rainbow')
    cax_p2=plt.contour(xc,yc,p_plot,colors='k')
    plt.clabel(cax_p2, fontsize=12, inline=1)
    cbar_p = plt.colorbar(cax_p)
    plt.title('Pressure Correction')
    plt.xlabel('x [m]')
    #plt.ylim([-1E1,1E1])
    plt.ylabel('y [m]' )
    
    plt.subplot(2,2,3)
    cax_a=plt.contourf(xc,yc,a_plot,cmap='rainbow')
    cbar_a=plt.colorbar(cax_a)
    plt.title('Concentration')
    plt.xlabel('x [m]')
    #plt.ylim([-1E1,1E1])
    plt.ylabel('y [m]' )
    
    pylab.show()

    #for c in (AIR,H2O):
    #  plot_isolines(t[c].val, (uf[c],vf[c],wf[c]), (xn[c],yn[c],zn[c]), Z)
     # plot_isolines(p_tot[c], (uf[c],vf[c],wf[c]), (xn[c],yn[c],zn[c]), Z)
