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

def radec_to_deg(ra, dec):
    # 将赤经从 "hh:mm:ss.ss" 转换为度数
    ra_h, ra_m, ra_s = map(float, ra.split(':'))
    ra_deg = 15 * (ra_h + ra_m / 60 + ra_s / 3600)  # 1小时 = 15度
    
    # 将赤纬从 "dd:mm:ss.ss" 转换为度数
    dec_d, dec_m, dec_s = map(float, dec.split(':'))
    dec_deg = dec_d + dec_m / 60 + dec_s / 3600
    
    return ra_deg, dec_deg

R_sun = 8.178  # kpc
Z_sun = 25 # pc
[U,V,W] = [7.01,10.13,4.95] # radial,tangential,vertival
rotv = 225 # clockwise,along with V
AU_km = u.AU.to(u.km,1)
yr_s = u.yr.to(u.s,1)

#定义银河坐标系，考虑太阳的位置和运动。
gc_frame = coord.Galactocentric(\
    galcen_distance=R_sun*u.kpc,\
    galcen_v_sun=[U,V+rotv,W]*(u.km / u.s),\
    z_sun=Z_sun*u.pc)
cir_v = 225 # km/s
sun_r = 8.178 # kpc

#目标天体的观测数据
ra = '14:43:25.76'  # 赤经（字符串格式）
dec = '+14:53:36.3'  # 赤纬（字符串格式）
ra, dec = radec_to_deg(ra, dec)
D    = 2.90
pmra = -46.914# mas/yr（mas代表毫角秒）
pmde = -1.465# mas/yr
rv   = 194.25

#将观测数据转换为坐标
c    = coord.SkyCoord(ra=ra*u.degree, dec=dec*u.degree,\
                        distance=D*u.kpc,\
                        pm_ra_cosdec=pmra*u.mas/u.yr,\
                        pm_dec=pmde*u.mas/u.yr,\
                        radial_velocity=rv*u.km/u.s,\
                        frame='icrs')
o = Orbit(c)
ts = np.linspace(0, -5e7, int(5e7 / 0.1e6) + 1) * u.yr  # 5 Gyr, 0.1 Myr 步长 ,
o.integrate(ts,MWPotential2014)# integrate the orbit
# get the orbit at each time step
data = [[o.x(t_),o.y(t_),o.z(t_)]for t_ in ts]
coor = np.array(data)
coor = coor.T


# 按照测量误差进行MC sampling
sampleNo = 10
De   = 0.5
pmrae= 0.121
pmdee= 0.096
rve  = 1.97
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
                
fig = plt.figure(figsize=(14, 6.5))
plt.subplots_adjust(left=0.05, right=0.98, top=0.88, bottom=0.08,
                    hspace=0.2, wspace=0.2)
grid = plt.GridSpec(1, 3)  # 增加一个子图

# 画 x-y 投影图
ax1 = plt.subplot(grid[0, 0])
plt.xlabel('x/kpc')
plt.ylabel('y/kpc')

# 画 y-z 投影图
ax2 = plt.subplot(grid[0, 1])
plt.xlabel('y/kpc')
plt.ylabel('z/kpc')

# 画 x-z 投影图
ax3 = plt.subplot(grid[0, 2])
plt.xlabel('x/kpc')
plt.ylabel('z/kpc')

# 画不考虑误差的轨道
ax1.plot(-coor[0], coor[1], zorder=10)
ax2.plot(coor[1], coor[2], zorder=10)
ax3.plot(-coor[0], coor[2], zorder=10)

# 标记太阳的位置
ax1.scatter(-R_sun, 0, zorder=10, c='red', label='Sun')  # 太阳的位置
ax2.scatter(0, 0, zorder=10, c='red', label='Sun')  # 太阳的位置
ax3.scatter(-R_sun, 0, zorder=10, c='red', label='Sun')  # 太阳的位置

# 标记银河系中心的位置
ax1.scatter(0, 0, c='black', label='Galactic Center')
ax2.scatter(0, 0, c='black', label='Galactic Center')
ax3.scatter(0, 0, c='black', label='Galactic Center')

# 添加图例
ax1.legend(loc='upper right')
ax2.legend(loc='upper right')
ax3.legend(loc='upper right')

# 考虑误差
ras = np.array([ra for j in range(sampleNo)])
decs = np.array([dec for j in range(sampleNo)])
pmras = np.array([s[j][0] for j in range(sampleNo)])
pmdecs = np.array([s[j][1] for j in range(sampleNo)])
rvs = np.array([s[j][2] for j in range(sampleNo)])
Ds = np.array([s[j][3] for j in range(sampleNo)])
cs = coord.SkyCoord(ra=ras * u.degree, dec=decs * u.degree,
                    distance=Ds * u.kpc,
                    pm_ra_cosdec=pmras * u.mas / u.yr,
                    pm_dec=pmdecs * u.mas / u.yr,
                    radial_velocity=rvs * u.km / u.s,
                    frame='icrs')
os = Orbit(cs)
os.integrate(ts, MWPotential2014)
index = 0

# 画 MC sampling 的轨道（有误差）
for o in os:
    index += 1
    print(index)
    data = [[o.x(t_), o.y(t_), o.z(t_)] for t_ in ts]
    coor = np.array(data)
    coor = coor.T
    ax1.plot(-coor[0], coor[1], c='grey', alpha=0.2)
    ax2.plot(coor[1], coor[2], c='grey', alpha=0.2)
    ax3.plot(-coor[0], coor[2], c='grey', alpha=0.2)

plt.show()
