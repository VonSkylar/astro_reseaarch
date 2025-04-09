import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, MaskedColumn
from astropy import units as u
import astropy.coordinates as coord
import scipy.stats as stats
from galpy.potential import MiyamotoNagaiPotential,\
     NFWPotential,HernquistPotential,MWPotential2014
from galpy.potential.mwpotentials import Irrgang13I
from galpy.potential import plotRotcurve
from galpy.potential import plotEscapecurve
from galpy.potential import vesc
from galpy.orbit import Orbit
from galpy.util.conversion import get_physical
from matplotlib import rcParams
rcParams.update({'font.size':13})

R_sun = 8.2  # kpc
Z_sun = 20.8 # pc
[U,V,W] = [11.1,12.24,7.25] # radial,tangential,vertival
rotv = 230 # clockwise,along with V
AU_km = u.AU.to(u.km,1)
yr_s = u.yr.to(u.s,1)

#定义银河坐标系，考虑太阳的位置和运动。
gc_frame = coord.Galactocentric(\
    galcen_distance=R_sun*u.kpc,\
    galcen_v_sun=[U,V+rotv,W]*(u.km / u.s),\
    z_sun=Z_sun*u.pc)
cir_v = 230 # km/s
sun_r = 8.2 # kpc

#目标天体的观测数据
ra   = 25.278694
dec  = 16.532547
D    = 37.28
pmra = 1.461# mas/yr（mas代表毫角秒）
pmde = -1.564# mas/yr
rv   = 384.3

#将观测数据转换为坐标
c    = coord.SkyCoord(ra=ra*u.degree, dec=dec*u.degree,\
                        distance=D*u.kpc,\
                        pm_ra_cosdec=pmra*u.mas/u.yr,\
                        pm_dec=pmde*u.mas/u.yr,\
                        radial_velocity=rv*u.km/u.s,\
                        frame='icrs')
o = Orbit(c)
ts = np.linspace(0,-1e9,1001)*u.yr # a total time of 1 Gyr,
                                   # with 1000 steps
o.integrate(ts,MWPotential2014)# integrate the orbit
# get the orbit at each time step
data = [[o.x(t_),o.y(t_),o.z(t_)]for t_ in ts]
coor = np.array(data)
coor = coor.T


# 按照测量误差进行MC sampling
sampleNo = 10
De   = 3.42
pmrae= 0.117
pmdee= 0.093
rve  = 4.0
mu  = [pmra,pmde,rv,D]
sx2 = pmrae**2
sy2 = pmdee**2
sz2 = rve**2
su2 = De**2
Sigma = np.array([[sx2,0  ,0  ,0  ],
                  [0  ,sy2,0  ,0  ],
                  [0  ,0  ,sz2,0  ],
                  [0  ,0  ,0  ,su2]])
s = np.random.multivariate_normal(mu,Sigma,sampleNo)
                
fig = plt.figure(figsize=(14,6.5))
plt.subplots_adjust(left=0.05,right=0.98,top=0.88,bottom=0.08,
                    hspace=0.2,wspace=0.2)
grid = plt.GridSpec(1,2)
# 画轨道的投影图
ax1 = plt.subplot(grid[0,0])
plt.xlabel('x/kpc')
plt.ylabel('z/kpc')

ax2 = plt.subplot(grid[0,1])
plt.xlabel('y/kpc')
ax2.yaxis.set_ticklabels([]) # 不显示刻度

# 画不考虑误差的轨道
ax1.plot(-coor[0],coor[2],zorder=10)
ax2.plot(coor[1],coor[2],zorder=10)
ax1.scatter(-coor[0][0],coor[2][0],zorder=10,c='r')
ax2.scatter(coor[1][0],coor[2][0],zorder=10,c='r')
ax1.scatter(0,0,c='black')
ax2.scatter(0,0,c='black')

#考虑误差
ras    = np.array([ra  for j in range(sampleNo)])
decs   = np.array([dec for j in range(sampleNo)])
pmras  = np.array([s[j][0] for j in range(sampleNo)])
pmdecs = np.array([s[j][1] for j in range(sampleNo)])
rvs    = np.array([s[j][2] for j in range(sampleNo)])
Ds     = np.array([s[j][3] for j in range(sampleNo)])
cs = coord.SkyCoord(ra=ras*u.degree, dec=decs*u.degree,\
                    distance=Ds*u.kpc,\
                    pm_ra_cosdec=pmras*u.mas/u.yr,\
                    pm_dec=pmdecs*u.mas/u.yr,\
                    radial_velocity=rvs*u.km/u.s,\
                    frame='icrs')
os = Orbit(cs)
os.integrate(ts,MWPotential2014)
index = 0

# 画MC sampling的轨道（有误差）
for o in os:
    index += 1
    print(index)
    data = [[o.x(t_),o.y(t_),o.z(t_)]for t_ in ts]

    coor = np.array(data)
    coor = coor.T
    ax1.plot(-coor[0],coor[2],c='grey',alpha=0.2)
    ax2.plot(coor[1],coor[2],c='grey',alpha=0.2)
plt.show()
