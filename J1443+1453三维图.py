import numpy as np
import matplotlib.pyplot as plt
from astropy import units as u
import astropy.coordinates as coord
from galpy.orbit import Orbit
from galpy.potential import MWPotential2014
from mpl_toolkits.mplot3d import Axes3D

# 定义银河系参数
R_sun = 8.178  # 太阳到银河系中心的距离，单位 kpc
Z_sun = 25  # 太阳相对于银河系平面的偏移量，单位 pc
rotv = 225  # 银河系的旋转速度，单位 km/s

# 创建三维图形
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')

# 设置三维图的坐标轴标签
ax.set_xlabel('x (kpc)')
ax.set_ylabel('y (kpc)')
ax.set_zlabel('z (kpc)')

# 设置坐标轴范围
ax.set_xlim(-25, 25)
ax.set_ylim(-25, 25)
ax.set_zlim(-20, 5)

# 添加太阳圈和银河系盘边缘的标记
theta = np.linspace(0, 2 * np.pi, 100)
sun_circle_x = R_sun * np.cos(theta)
sun_circle_y = R_sun * np.sin(theta)
ax.plot(sun_circle_x, sun_circle_y, 0, color='black', linestyle='dotted', label='Solar Circle (R=8.178 kpc)')

disk_edge_x = 25 * np.cos(theta)
disk_edge_y = 25 * np.sin(theta)
ax.plot(disk_edge_x, disk_edge_y, 0, color='black', linestyle='dotted', label='MW Disk Edge (R=25 kpc)')

# 标记太阳的位置
ax.scatter(-R_sun, 0, 0, c='blue', label='Sun', s=50)

# 标记银河系中心的位置
ax.scatter(0, 0, 0, c='red', label='Galactic Center', s=50)

# 目标天体的观测数据
ra = '14:43:25.76'  # 赤经（字符串格式）
dec = '+14:53:36.3'  # 赤纬（字符串格式）
distance = 2.90  # 距离，单位 kpc
pmra = -46.914  # 赤经自行，单位 mas/yr
pmdec = -1.465  # 赤纬自行，单位 mas/yr
rv = 194.25  # 径向速度，单位 km/s

# 将赤经和赤纬转换为带有单位的值
ra = coord.Angle(ra, unit=u.hourangle)  # 赤经单位：小时
dec = coord.Angle(dec, unit=u.degree)  # 赤纬单位：度

# 将目标天体的观测数据转换为坐标
target_coord = coord.SkyCoord(ra=ra, dec=dec, distance=distance*u.kpc,
                              pm_ra_cosdec=pmra*u.mas/u.yr,
                              pm_dec=pmdec*u.mas/u.yr,
                              radial_velocity=rv*u.km/u.s,
                              frame='icrs')

# 创建目标天体的 Orbit 对象
target_orbit = Orbit(target_coord)

# 对目标天体的轨道进行积分
ts = np.linspace(0, -5e7, int(5e7 / 0.1e6) + 1) * u.yr  # 时间步长
target_orbit.integrate(ts, MWPotential2014)

# 获取目标天体的轨道数据
target_data = [[target_orbit.x(t_), target_orbit.y(t_), target_orbit.z(t_)] for t_ in ts]
target_coor = np.array(target_data)
target_coor = target_coor.T

# 绘制目标天体的轨道
ax.scatter(-target_coor[0][0], target_coor[1][0], target_coor[2][0], c='black', s=50)
ax.plot(-target_coor[0], target_coor[1], target_coor[2], label='Target Orbit', c='black', zorder=10, linewidth=2)

# 在目标天体轨道末尾添加实心箭头
arrow_index = -1  # 最后一个点
ax.quiver(-target_coor[0][arrow_index], target_coor[1][arrow_index], target_coor[2][arrow_index],
          -target_coor[0][arrow_index] + target_coor[0][arrow_index-1],
          target_coor[1][arrow_index] - target_coor[1][arrow_index-1],
          target_coor[2][arrow_index] - target_coor[2][arrow_index-1],
          color='black', length=5, normalize=True, arrow_length_ratio=0.3)

# 使用 Orbit.from_name 加载人马座矮星系的轨道
sgr_orbit = Orbit.from_name('Sagittarius Dwarf Spheroidal')

# 对 Sgr dSph 的轨道进行积分
sgr_orbit.integrate(ts, MWPotential2014)

# 获取 Sgr dSph 的轨道数据
sgr_data = [[sgr_orbit.x(t_), sgr_orbit.y(t_), sgr_orbit.z(t_)] for t_ in ts]
sgr_coor = np.array(sgr_data)
sgr_coor = sgr_coor.T

# 绘制人马座矮星系的轨道
ax.scatter(-sgr_coor[0][0], sgr_coor[1][0], sgr_coor[2][0], c='green', s=50)
ax.plot(-sgr_coor[0], sgr_coor[1], sgr_coor[2], label='Sgr dSph Orbit', c='green', zorder=10, linewidth=2)

# 在 Sgr dSph 轨道末尾添加实心箭头
arrow_index = -1  # 最后一个点
ax.quiver(-sgr_coor[0][arrow_index], sgr_coor[1][arrow_index], sgr_coor[2][arrow_index],
          -sgr_coor[0][arrow_index] + sgr_coor[0][arrow_index-1],
          sgr_coor[1][arrow_index] - sgr_coor[1][arrow_index-1],
          sgr_coor[2][arrow_index] - sgr_coor[2][arrow_index-1],
          color='green', length=5, normalize=True, arrow_length_ratio=0.3)

# 添加图例
ax.legend(loc='upper right')

# 调整视角
ax.view_init(elev=15, azim=-45)  # 调整视角，默认角度为仰角 30 度，方位角 45 度

# 显示图像
plt.tight_layout()
plt.show()