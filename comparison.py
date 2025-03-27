import numpy as np
import matplotlib.pyplot as plt
from astropy import units as u
import astropy.coordinates as coord
from galpy.orbit import Orbit
from galpy.potential import MWPotential2014, SolidBodyRotationWrapperPotential
from galpy.potential import DehnenBarPotential, LogarithmicHaloPotential
from mpl_toolkits.mplot3d import Axes3D

# 定义银河系参数
R_sun = 8.178  # 太阳到银河系中心的距离，单位 kpc
Z_sun = 25 / 1000  # 太阳相对于银河系平面的偏移量，单位 kpc
rotv = 225  # 银河系的旋转速度，单位 km/s

# 创建 galpy 棒旋结构势场
# 使用 Dehnen 棒状势场模型作为近似
bar_strength = 0.3  # 棒的相对强度
bar_angle = 25.0 * np.pi/180  # 棒与太阳-银心连线的夹角，单位弧度
bar_pattern_speed = 39.0  # 棒的模式速度，单位 km/s/kpc
bar_length = 3.5  # 棒的长度，单位 kpc

# 创建 Dehnen 棒势场
dpn = DehnenBarPotential(
    omegab=bar_pattern_speed, 
    rb=bar_length,
    Af=bar_strength,
    alpha=0.01,  # 时间开启强度参数
    tsteady=1.0  # 棒结构达到稳定的时间，单位 Gyr
)

# 创建包含棒结构的完整势场模型
bar_potential = SolidBodyRotationWrapperPotential(pot=MWPotential2014 + [dpn],
                                                 omega=bar_pattern_speed,
                                                 pa=bar_angle)

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

# 创建人马座矮星系的轨道
sgr_orbit = Orbit.from_name('Sagittarius Dwarf Spheroidal')

# 创建图形
fig = plt.figure(figsize=(10, 6))

# 在 MWPotential2014 中计算轨道
ax1 = fig.add_subplot(121, projection='3d')
ax1.set_title('Orbit in MWPotential2014')
ax1.set_xlabel('x (kpc)')
ax1.set_ylabel('y (kpc)')
ax1.set_zlabel('z (kpc)')
ax1.set_xlim(-25, 25)
ax1.set_ylim(-25, 25)
ax1.set_zlim(-20, 5)

# 绘制太阳圈和银河系盘
theta = np.linspace(0, 2 * np.pi, 100)
sun_circle_x = R_sun * np.cos(theta)
sun_circle_y = R_sun * np.sin(theta)
ax1.plot(sun_circle_x, sun_circle_y, 0, color='black', linestyle='dotted')

# 标记太阳和银河系中心
ax1.scatter(-R_sun, 0, 0, c='blue', label='Sun', s=50)
ax1.scatter(0, 0, 0, c='red', label='Galactic Center', s=50)

# 计算 MWPotential2014 中的轨道
#ts = np.linspace(0, -5, 1001) * u.Gyr  # 时间步长
ts = np.linspace(0, -5e7, int(5e7 / 0.1e6) + 1) * u.yr  # 时间步长
target_orbit.integrate(ts, MWPotential2014)
sgr_orbit.integrate(ts, MWPotential2014)

# 绘制目标天体的轨道
target_data = [[target_orbit.x(t_), target_orbit.y(t_), target_orbit.z(t_)] for t_ in ts]
target_coor = np.array(target_data)
target_coor = target_coor.T
ax1.scatter(-target_coor[0][0], target_coor[1][0], target_coor[2][0], c='black', s=50, label='Target initial position')
ax1.plot(-target_coor[0], target_coor[1], target_coor[2], label='Target Orbit (MWPotential2014)', c='black', linewidth=2)

# 绘制人马座矮星系的轨道
sgr_data = [[sgr_orbit.x(t_), sgr_orbit.y(t_), sgr_orbit.z(t_)] for t_ in ts]
sgr_coor = np.array(sgr_data)
sgr_coor = sgr_coor.T
ax1.scatter(-sgr_coor[0][0], sgr_coor[1][0], sgr_coor[2][0], c='green', s=50, label='Sgr dSph initial position')
ax1.plot(-sgr_coor[0], sgr_coor[1], sgr_coor[2], label='Sgr dSph Orbit (MWPotential2014)', c='green', linewidth=2)

# 在包含棒旋结构的势场中计算轨道
ax2 = fig.add_subplot(122, projection='3d')
ax2.set_title('Orbit in MW + Bar Potential')
ax2.set_xlabel('x (kpc)')
ax2.set_ylabel('y (kpc)')
ax2.set_zlabel('z (kpc)')
ax2.set_xlim(-25, 25)
ax2.set_ylim(-25, 25)
ax2.set_zlim(-20, 5)

# 复制太阳圈和标记
ax2.plot(sun_circle_x, sun_circle_y, 0, color='black', linestyle='dotted')
ax2.scatter(-R_sun, 0, 0, c='blue', label='Sun', s=50)
ax2.scatter(0, 0, 0, c='red', label='Galactic Center', s=50)

# 绘制棒结构的示意图（简单的椭圆）
bar_ellipse_a = bar_length
bar_ellipse_b = bar_length * 0.4  # 假设宽度是长度的40%
t = np.linspace(0, 2*np.pi, 100)
x_bar = bar_ellipse_a * np.cos(t)
y_bar = bar_ellipse_b * np.sin(t)
# 旋转坐标
x_rot = x_bar * np.cos(bar_angle) - y_bar * np.sin(bar_angle)
y_rot = x_bar * np.sin(bar_angle) + y_bar * np.cos(bar_angle)
ax2.plot(x_rot, y_rot, 0, 'r-', linewidth=2, alpha=0.5, label='Bar (schematic)')

# 重新设置初始条件
target_orbit_bar = Orbit(target_coord)
sgr_orbit_bar = Orbit.from_name('Sagittarius Dwarf Spheroidal')

# 计算在棒旋势场中的轨道
target_orbit_bar.integrate(ts, bar_potential)
sgr_orbit_bar.integrate(ts, bar_potential)

# 绘制目标天体在棒旋势场中的轨道
target_bar_data = [[target_orbit_bar.x(t_), target_orbit_bar.y(t_), target_orbit_bar.z(t_)] for t_ in ts]
target_bar_coor = np.array(target_bar_data)
target_bar_coor = target_bar_coor.T
ax2.scatter(-target_bar_coor[0][0], target_bar_coor[1][0], target_bar_coor[2][0], c='navy', s=50, label='Target initial position')
ax2.plot(-target_bar_coor[0], target_bar_coor[1], target_bar_coor[2], label='Target Orbit (with Bar)', c='navy', linewidth=2)

# 绘制人马座矮星系在棒旋势场中的轨道
sgr_bar_data = [[sgr_orbit_bar.x(t_), sgr_orbit_bar.y(t_), sgr_orbit_bar.z(t_)] for t_ in ts]
sgr_bar_coor = np.array(sgr_bar_data)
sgr_bar_coor = sgr_bar_coor.T
ax2.scatter(-sgr_bar_coor[0][0], sgr_bar_coor[1][0], sgr_bar_coor[2][0], c='darkgreen', s=50, label='Sgr dSph initial position')
ax2.plot(-sgr_bar_coor[0], sgr_bar_coor[1], sgr_bar_coor[2], label='Sgr dSph Orbit (with Bar)', c='darkgreen', linewidth=2)

# 添加图例和调整视角
for ax in [ax1, ax2]:
    ax.legend(loc='upper right')
    ax.view_init(elev=15, azim=-45)

# 添加解释性标题
plt.suptitle('Comparison of Orbit Integration: Axisymmetric vs. Barred Potential', fontsize=16)

# 显示图像
plt.tight_layout()
plt.legend(framealpha=0.5) 
plt.legend(loc='best')
plt.subplots_adjust(top=0.9)
plt.show()

# 计算和绘制轨道差异
fig2, axs = plt.subplots(3, 2, figsize=(10, 6))

# 目标天体轨道的差异
R_mw = np.sqrt(target_coor[0]**2 + target_coor[1]**2)
R_bar = np.sqrt(target_bar_coor[0]**2 + target_bar_coor[1]**2)
t_Gyr = np.abs(ts.to(u.Gyr).value)  # 时间，单位 Gyr

# 左列：目标天体
axs[0, 0].plot(t_Gyr, R_bar - R_mw)
axs[0, 0].set_ylabel('ΔR (kpc)')
axs[0, 0].set_title('Target Orbit Differences')
axs[0, 0].axhline(0, color='k', linestyle='--', linewidth=0.8)

axs[1, 0].plot(t_Gyr, target_bar_coor[2] - target_coor[2])
axs[1, 0].set_ylabel('Δz (kpc)')
axs[1, 0].axhline(0, color='k', linestyle='--', linewidth=0.8)

# 总的空间差异
total_diff_target = np.sqrt((target_bar_coor[0] - target_coor[0])**2 + 
                     (target_bar_coor[1] - target_coor[1])**2 + 
                     (target_bar_coor[2] - target_coor[2])**2)
axs[2, 0].plot(t_Gyr, total_diff_target)
axs[2, 0].set_ylabel('Total difference (kpc)')
axs[2, 0].set_xlabel('Time (Gyr)')

# 人马座矮星系轨道的差异
R_sgr_mw = np.sqrt(sgr_coor[0]**2 + sgr_coor[1]**2)
R_sgr_bar = np.sqrt(sgr_bar_coor[0]**2 + sgr_bar_coor[1]**2)

# 右列：人马座矮星系
axs[0, 1].plot(t_Gyr, R_sgr_bar - R_sgr_mw)
axs[0, 1].set_ylabel('ΔR (kpc)')
axs[0, 1].set_title('Sgr dSph Orbit Differences')
axs[0, 1].axhline(0, color='k', linestyle='--', linewidth=0.8)

axs[1, 1].plot(t_Gyr, sgr_bar_coor[2] - sgr_coor[2])
axs[1, 1].set_ylabel('Δz (kpc)')
axs[1, 1].axhline(0, color='k', linestyle='--', linewidth=0.8)

# 总的空间差异
total_diff_sgr = np.sqrt((sgr_bar_coor[0] - sgr_coor[0])**2 + 
                 (sgr_bar_coor[1] - sgr_coor[1])**2 + 
                 (sgr_bar_coor[2] - sgr_coor[2])**2)
axs[2, 1].plot(t_Gyr, total_diff_sgr)
axs[2, 1].set_ylabel('Total difference (kpc)')
axs[2, 1].set_xlabel('Time (Gyr)')

fig2.suptitle('Orbit Differences: Axisymmetric vs. Barred Potential', fontsize=16)
plt.tight_layout()
plt.legend(loc='best')
plt.subplots_adjust(top=0.9)
plt.show()