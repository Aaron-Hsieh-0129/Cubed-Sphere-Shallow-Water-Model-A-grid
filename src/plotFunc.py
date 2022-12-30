import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import cartopy.util as cutil
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER


######## Should Be Tuned ########### 
left, right, split = 8000, 12000, 21
wind = 20
skip_car = 12
scale_car = 4000
skip = 4
scale = 1000
skip_sph = 4
scale_sph = 3000

left_zeta, right_zeta, split_zeta = -16 * 10 ** (-5), 16 * 10 ** (-5), 17
######## Should Be Tuned ###########

DX = DY = 2
NX = NY = int(90 / DX)
DT = 180
LEAP = 5

cmap = cm.viridis
wind_color = "#5c5959"
DPI = 100
fs = 15

def plotOnCubeWindMul(t):
    x, y = np.loadtxt("../outputs/grids/x.txt").reshape(6, NX, NY), np.loadtxt("../outputs/grids/y.txt").reshape(6, NX, NY)
    val = np.loadtxt(f"../outputs/h/h_{t*LEAP}.txt").reshape(6, NX, NY)
    u = np.loadtxt(f"../outputs/u/u_{t*LEAP}.txt").reshape(6, NX, NY)
    v = np.loadtxt(f"../outputs/v/v_{t*LEAP}.txt").reshape(6, NX, NY)

    fig = plt.figure(figsize=(18,10))
    ax5 = fig.add_subplot(3,4,2)
    ax4 = fig.add_subplot(3,4,5)
    ax1 = fig.add_subplot(3,4,6)
    ax2 = fig.add_subplot(3,4,7)
    ax3 = fig.add_subplot(3,4,8)
    ax6 = fig.add_subplot(3,4,10)
    ax5.set_title(f"t = {t * LEAP * DT / 60} min", fontsize=fs)

    cs1 = ax1.contourf(x[0], y[0], val[0], levels=np.linspace(left, right, split), extend='both', cmap=cmap)
    ax2.contourf(x[1], y[1], val[1], levels=np.linspace(left, right, split), extend='both', cmap=cmap)
    ax3.contourf(x[2], y[2], val[2], levels=np.linspace(left, right, split), extend='both', cmap=cmap)
    ax4.contourf(x[3], y[3], val[3], levels=np.linspace(left, right, split), extend='both', cmap=cmap)
    ax5.contourf(x[4], y[4], val[4], levels=np.linspace(left, right, split), extend='both', cmap=cmap)
    ax6.contourf(x[5], y[5], val[5], levels=np.linspace(left, right, split), extend='both', cmap=cmap)
    cb_ax1 = fig.add_axes([0.9235, 0.1, 0.015, 0.78])
    fig.colorbar(cs1, cax=cb_ax1, ticks=np.linspace(left, right, split))
    
    Q = ax1.quiver(x[0, ::skip, ::skip], y[0, ::skip, ::skip], u[0, ::skip, ::skip], v[0, ::skip, ::skip], angles='xy', units="width", scale=scale, color=wind_color)
    ax2.quiver(x[1, ::skip, ::skip], y[1, ::skip, ::skip], u[1, ::skip, ::skip], v[1, ::skip, ::skip], angles='xy', units="width", scale=scale, color=wind_color)
    ax3.quiver(x[2, ::skip, ::skip], y[2, ::skip, ::skip], u[2, ::skip, ::skip], v[2, ::skip, ::skip], angles='xy', units="width", scale=scale, color=wind_color)
    ax4.quiver(x[3, ::skip, ::skip], y[3, ::skip, ::skip], u[3, ::skip, ::skip], v[3, ::skip, ::skip], angles='xy', units="width", scale=scale, color=wind_color)
    ax5.quiver(x[4, ::skip, ::skip], y[4, ::skip, ::skip], u[4, ::skip, ::skip], v[4, ::skip, ::skip], angles='xy', units="width", scale=scale, color=wind_color)
    ax6.quiver(x[5, ::skip, ::skip], y[5, ::skip, ::skip], u[5, ::skip, ::skip], v[5, ::skip, ::skip], angles='xy', units="width", scale=scale, color=wind_color)

    ax1.quiverkey(Q, 0.7, 0.9, wind, f"{wind}" + r'$ \frac{m}{s}$', labelpos='E', coordinates='figure')  
    ax2.quiverkey(Q, 0.7, 0.9, wind, f"{wind}" + r'$ \frac{m}{s}$', labelpos='E', coordinates='figure')  
    ax3.quiverkey(Q, 0.7, 0.9, wind, f"{wind}" + r'$ \frac{m}{s}$', labelpos='E', coordinates='figure')  
    ax4.quiverkey(Q, 0.7, 0.9, wind, f"{wind}" + r'$ \frac{m}{s}$', labelpos='E', coordinates='figure')  
    ax5.quiverkey(Q, 0.7, 0.9, wind, f"{wind}" + r'$ \frac{m}{s}$', labelpos='E', coordinates='figure')  
    ax6.quiverkey(Q, 0.7, 0.9, wind, f"{wind}" + r'$ \frac{m}{s}$', labelpos='E', coordinates='figure')   
    
    plt.savefig(f"../graphs/h/curvilinear/{int(t/LEAP/2)}.png", dpi=DPI)
    plt.close()
    return

def plotOnSphereWindMul(t):
    x, y = np.loadtxt("../outputs/grids/lon.txt").reshape(6, NX, NY), np.loadtxt("../outputs/grids/lat.txt").reshape(6, NX, NY)
    val = np.loadtxt(f"../outputs/h/h_{t*LEAP}.txt").reshape(6, NX, NY)
    x = x * 180 / np.pi
    y = y * 180 / np.pi

    x[0, :, :NX//2] = x[0, :, :NX//2] - 360
    x[4, :, :NX//2] = x[4, :, :NX//2] - 360
    x[5, :, :NX//2] = x[5, :, :NX//2] - 360

    u = np.loadtxt(f"../outputs/u_lon_lat/u_lon_lat_{t*LEAP}.txt").reshape(6, NX, NY)
    v = np.loadtxt(f"../outputs/v_lon_lat/v_lon_lat_{t*LEAP}.txt").reshape(6, NX, NY)

    plt.figure(figsize=(18,8))
    plt.xlabel("LON")
    plt.ylabel("LAT")
    plt.title(f"t = {t * LEAP * DT / 60} min", fontsize=fs)

    plt.contourf(x[0], y[0], val[0], levels=np.linspace(left, right, split), extend='both', cmap=cmap)
    plt.contourf(x[1], y[1], val[1], levels=np.linspace(left, right, split), extend='both', cmap=cmap)
    plt.contourf(x[2, :, 0:NX//2], y[2, :, 0:NX//2], val[2, :, 0:NX//2], levels=np.linspace(left, right, split), extend='both', cmap=cmap)
    plt.contourf(x[2, :, NX//2:]-360, y[2, :, NX//2:], val[2, :, NX//2:], levels=np.linspace(left, right, split), extend='both', cmap=cmap)
    plt.contourf(x[3]-360, y[3], val[3], levels=np.linspace(left, right, split), extend='both', cmap=cmap)
    plt.contourf(x[4], y[4], val[4], levels=np.linspace(left, right, split), extend='both', cmap=cmap)
    plt.contourf(x[5], y[5], val[5], levels=np.linspace(left, right, split), extend='both', cmap=cmap)
    cbar = plt.colorbar(pad=0.05)
    cbar.set_ticks(np.linspace(left, right, split))
    
    plt.quiver(x[0][::skip_sph, ::skip_sph], y[0][::skip_sph, ::skip_sph], u[0][::skip_sph, ::skip_sph], v[0][::skip_sph, ::skip_sph], angles='xy', units="width", scale=scale_sph, color=wind_color)
    plt.quiver(x[1][::skip_sph, ::skip_sph], y[1][::skip_sph, ::skip_sph], u[1][::skip_sph, ::skip_sph], v[1][::skip_sph, ::skip_sph], angles='xy', units="width", scale=scale_sph, color=wind_color)
    plt.quiver(x[2, :, 0:NX//2][::skip_sph, ::skip_sph], y[2, :, 0:NX//2][::skip_sph, ::skip_sph], u[2, :, 0:NX//2][::skip_sph, ::skip_sph], v[2, :, 0:NX//2][::skip_sph, ::skip_sph], angles='xy', units="width", scale=scale_sph, color=wind_color)
    plt.quiver(x[2, :, NX//2:][::skip_sph, ::skip_sph]-360, y[2, :, NX//2:][::skip_sph, ::skip_sph], u[2, :, NX//2:][::skip_sph, ::skip_sph], v[2, :, NX//2:][::skip_sph, ::skip_sph], angles='xy', units="width", scale=scale_sph, color=wind_color)
    plt.quiver(x[3][::skip_sph, ::skip_sph]-360, y[3][::skip_sph, ::skip_sph], u[3][::skip_sph, ::skip_sph], v[3][::skip_sph, ::skip_sph], angles='xy', units="width", scale=scale_sph, color=wind_color)
    plt.quiver(x[4][::skip_sph, ::skip_sph], y[4][::skip_sph, ::skip_sph], u[4][::skip_sph, ::skip_sph], v[4][::skip_sph, ::skip_sph], angles='xy', units="width", scale=scale_sph, color=wind_color)
    Q = plt.quiver(x[5][::skip_sph, ::skip_sph], y[5][::skip_sph, ::skip_sph], u[5][::skip_sph, ::skip_sph], v[5][::skip_sph, ::skip_sph], angles='xy', units="width", scale=scale_sph, color=wind_color)
    qk = plt.quiverkey(Q, 0.7, 0.9, wind, f"{wind}" + r'$ \frac{m}{s}$', labelpos='E', coordinates='figure')
        
    plt.savefig(f"../graphs/h/sphere/{int(t/LEAP/2)}.png", dpi=DPI)
    plt.close()
    return

def plotOnSphereMul(t):
    x, y = np.loadtxt("../outputs/grids/lon.txt").reshape(6, NX, NY), np.loadtxt("../outputs/grids/lat.txt").reshape(6, NX, NY)
    val = np.loadtxt(f"../outputs/h/h_{t*LEAP}.txt").reshape(6, NX, NY)
    x = x * 180 / np.pi
    y = y * 180 / np.pi

    x[0, :, :NX//2] = x[0, :, :NX//2] - 360
    x[4, :, :NX//2] = x[4, :, :NX//2] - 360
    x[5, :, :NX//2] = x[5, :, :NX//2] - 360

    u = np.loadtxt(f"../outputs/u_lon_lat/u_lon_lat_{t*LEAP}.txt").reshape(6, NX, NY)
    v = np.loadtxt(f"../outputs/v_lon_lat/v_lon_lat_{t*LEAP}.txt").reshape(6, NX, NY)

    plt.figure(figsize=(18,8))
    plt.xlabel("LON")
    plt.ylabel("LAT")
    plt.title(f"t = {t * LEAP * DT / 60} min", fontsize=fs)

    plt.contourf(x[0], y[0], val[0], levels=np.linspace(left, right, split), extend='both', cmap=cmap)
    plt.contourf(x[1], y[1], val[1], levels=np.linspace(left, right, split), extend='both', cmap=cmap)
    plt.contourf(x[2, :, 0:NX//2], y[2, :, 0:NX//2], val[2, :, 0:NX//2], levels=np.linspace(left, right, split), extend='both', cmap=cmap)
    plt.contourf(x[2, :, NX//2:]-360, y[2, :, NX//2:], val[2, :, NX//2:], levels=np.linspace(left, right, split), extend='both', cmap=cmap)
    plt.contourf(x[3]-360, y[3], val[3], levels=np.linspace(left, right, split), extend='both', cmap=cmap)
    plt.contourf(x[4], y[4], val[4], levels=np.linspace(left, right, split), extend='both', cmap=cmap)
    plt.contourf(x[5], y[5], val[5], levels=np.linspace(left, right, split), extend='both', cmap=cmap)
    cbar = plt.colorbar(pad=0.05)
    cbar.set_ticks(np.linspace(left, right, split))
      
    plt.savefig(f"../graphs/h/sphere/{int(t/LEAP/2)}.png", dpi=DPI)
    plt.close()
    return

def plotOnCubeMul(t):
    x, y = np.loadtxt("../outputs/grids/x.txt").reshape(6, NX, NY), np.loadtxt("../outputs/grids/y.txt").reshape(6, NX, NY)
    val = np.loadtxt(f"../outputs/h/h_{t*LEAP}.txt").reshape(6, NX, NY)
    u = np.loadtxt(f"../outputs/u/u_{t*LEAP}.txt").reshape(6, NX, NY)
    v = np.loadtxt(f"../outputs/v/v_{t*LEAP}.txt").reshape(6, NX, NY)

    fig = plt.figure(figsize=(18,10))
    ax5 = fig.add_subplot(3,4,2)
    ax4 = fig.add_subplot(3,4,5)
    ax1 = fig.add_subplot(3,4,6)
    ax2 = fig.add_subplot(3,4,7)
    ax3 = fig.add_subplot(3,4,8)
    ax6 = fig.add_subplot(3,4,10)
    ax5.set_title(f"t = {t * LEAP * DT / 60} min", fontsize=fs)

    cs1 = ax1.contourf(x[0], y[0], val[0], levels=np.linspace(left, right, split), extend='both', cmap=cmap)
    ax2.contourf(x[1], y[1], val[1], levels=np.linspace(left, right, split), extend='both', cmap=cmap)
    ax3.contourf(x[2], y[2], val[2], levels=np.linspace(left, right, split), extend='both', cmap=cmap)
    ax4.contourf(x[3], y[3], val[3], levels=np.linspace(left, right, split), extend='both', cmap=cmap)
    ax5.contourf(x[4], y[4], val[4], levels=np.linspace(left, right, split), extend='both', cmap=cmap)
    ax6.contourf(x[5], y[5], val[5], levels=np.linspace(left, right, split), extend='both', cmap=cmap)
    cb_ax1 = fig.add_axes([0.9235, 0.1, 0.015, 0.78])
    fig.colorbar(cs1, cax=cb_ax1, ticks=np.linspace(left, right, split))

    plt.savefig(f"../graphs/h/curvilinear/{int(t/LEAP/2)}.png", dpi=DPI)
    plt.close()
    return

def plotSphereCartopy(t):
    lon = np.loadtxt("../outputs/grids/lon.txt").reshape(6, NX, NY).flatten() * 180 / np.pi
    lat = np.loadtxt("../outputs/grids/lat.txt").reshape(6, NX, NY).flatten() * 180 / np.pi
    h = np.loadtxt(f"../outputs/h/h_{t*LEAP}.txt").reshape(6, NX, NY).flatten()

    map_proj = ccrs.PlateCarree(central_longitude=0.)
    data_crs = ccrs.PlateCarree()
    fig, ax = plt.subplots(nrows=1,ncols=1,figsize=(12,8),subplot_kw={'projection': map_proj},dpi=DPI)

    pos = ax.get_position()
    pset = [0.07, 0.05, 0.8, 0.9]
    ax.set_position(pset)

    x, y, _ = map_proj.transform_points(data_crs, lon, lat).T
    mask = np.invert(np.logical_or(np.isinf(x), np.isinf(y)))
    x = np.compress(mask, x)
    y = np.compress(mask, y)

    cs = ax.tricontourf(x,y,h[mask],extend="both",cmap=cmap,levels=np.linspace(left, right, split))
    ax.set_global()
    cb_ax1 = fig.add_axes([0.92, 0.16, 0.012, 0.67])
    fig.colorbar(cs, cax=cb_ax1, ticks=np.linspace(left, right, split))

    ax.coastlines(resolution='110m',color='k', lw=0.2, zorder=13)
    gl = ax.gridlines(draw_labels=True)

    ax.set_title(f"t = {t * LEAP * DT / 60} min", fontsize=fs)
    
    plt.savefig(f"../graphs/h/sphere_cartopy/{int(t/LEAP/2)}.png", dpi=DPI)
    plt.close()

def plotSphereWindCartopy(t):
    lon = np.loadtxt("../outputs/grids/lon.txt").reshape(6, NX, NY).flatten() * 180 / np.pi
    lat = np.loadtxt("../outputs/grids/lat.txt").reshape(6, NX, NY).flatten() * 180 / np.pi
    h = np.loadtxt(f"../outputs/h/h_{t*LEAP}.txt").reshape(6, NX, NY).flatten()
    u = np.loadtxt(f"../outputs/u_lon_lat/u_lon_lat_{t*LEAP}.txt").reshape(6, NX, NY).flatten()
    v = np.loadtxt(f"../outputs/v_lon_lat/v_lon_lat_{t*LEAP}.txt").reshape(6, NX, NY).flatten()

    map_proj = ccrs.PlateCarree(central_longitude=0.)
    data_crs = ccrs.PlateCarree()
    fig, ax = plt.subplots(nrows=1,ncols=1,figsize=(12,8),subplot_kw={'projection': map_proj},dpi=DPI)

    pos = ax.get_position()
    pset = [0.07, 0.05, 0.8, 0.9]
    ax.set_position(pset)

    x, y, _ = map_proj.transform_points(data_crs, lon, lat).T
    mask = np.invert(np.logical_or(np.isinf(x), np.isinf(y)))
    x = np.compress(mask, x)
    y = np.compress(mask, y)

    cs = ax.tricontourf(x,y,h[mask],extend="both",cmap=cmap,levels=np.linspace(left, right, split))
    ax.set_global()
    cb_ax1 = fig.add_axes([0.92, 0.16, 0.012, 0.67])
    fig.colorbar(cs, cax=cb_ax1, ticks=np.linspace(left, right, split))

    Q = ax.quiver(x[::skip_car], y[::skip_car], u[mask][::skip_car], v[mask][::skip_car], angles='xy', units="width", scale=scale_sph, color=wind_color)
    ax.quiverkey(Q, 0.85, 0.85, wind, f"{wind}" + r'$ \frac{m}{s}$', labelpos='E', coordinates='figure')

    ax.coastlines(resolution='110m',color='k', lw=0.2, zorder=13)
    gl = ax.gridlines(draw_labels=True)

    ax.set_title(f"t = {t * LEAP * DT / 60} min", fontsize=fs)
    
    plt.savefig(f"../graphs/h/sphere_cartopy/{int(t/LEAP/2)}.png", dpi=DPI)
    plt.close()

def plotWind():
    #### plot curvilinear wind
    u = np.loadtxt("../outputs/u/u_0.txt").reshape(6, NX, NY)
    v = np.loadtxt("../outputs/v/v_0.txt").reshape(6, NX, NY)
    x, y = np.loadtxt("../outputs/grids/x.txt").reshape(6, NX, NY), np.loadtxt("../outputs/grids/y.txt").reshape(6, NX, NY)

    fig = plt.figure(figsize=(18, 10))
    ax5 = fig.add_subplot(342)
    ax4 = fig.add_subplot(345)
    ax1 = fig.add_subplot(346)
    ax2 = fig.add_subplot(347)
    ax3 = fig.add_subplot(348)
    ax6 = fig.add_subplot(3,4,10)

    ax1.quiver(x[0][::skip, ::skip], y[0][::skip, ::skip], u[0][::skip, ::skip], v[0][::skip, ::skip], angles='xy', units="width", scale=scale)
    ax2.quiver(x[1][::skip, ::skip], y[1][::skip, ::skip], u[1][::skip, ::skip], v[1][::skip, ::skip], angles='xy', units="width", scale=scale)
    ax3.quiver(x[2][::skip, ::skip], y[2][::skip, ::skip], u[2][::skip, ::skip], v[2][::skip, ::skip], angles='xy', units="width", scale=scale)
    ax4.quiver(x[3][::skip, ::skip], y[3][::skip, ::skip], u[3][::skip, ::skip], v[3][::skip, ::skip], angles='xy', units="width", scale=scale)
    Q = ax5.quiver(x[4][::skip, ::skip], y[4][::skip, ::skip], u[4][::skip, ::skip], v[4][::skip, ::skip], angles='xy', units="width", scale=scale)
    ax6.quiver(x[5][::skip, ::skip], y[5][::skip, ::skip], u[5][::skip, ::skip], v[5][::skip, ::skip], angles='xy', units="width", scale=scale)
    qk = plt.quiverkey(Q, 0.48, 0.975, 30, r'$30 \frac{m}{s}$', labelpos='E', coordinates='figure')
    ax5.set_title("Curvilinear Coordinate")
    plt.tight_layout()
    plt.savefig("../graphs/wind/curvilinear.png", dpi=150)
    plt.close()

    #### plot sphere wind
    u = np.loadtxt("../outputs/u_lon_lat/u_lon_lat_0.txt").reshape(6, NX, NY)
    v = np.loadtxt("../outputs/v_lon_lat/v_lon_lat_0.txt").reshape(6, NX, NY)
    lon = np.loadtxt("../outputs/grids/lon.txt").reshape(6, NX, NY) * 180. / np.pi
    lat = np.loadtxt("../outputs/grids/lat.txt").reshape(6, NX, NY) * 180. / np.pi

    lon[0, :, :NX//2] = lon[0, :, :NX//2] - 360
    lon[4, :, :NX//2] = lon[4, :, :NX//2] - 360
    lon[5, :, :NX//2] = lon[5, :, :NX//2] - 360

    plt.figure(figsize=(18, 12))
    plt.title("Spherical Coordinate")

    plt.quiver(lon[0][::skip_sph, ::skip_sph], lat[0][::skip_sph, ::skip_sph], u[0][::skip_sph, ::skip_sph], v[0][::skip_sph, ::skip_sph], angles='xy', units="width", scale=scale_sph)
    plt.quiver(lon[1][::skip_sph, ::skip_sph], lat[1][::skip_sph, ::skip_sph], u[1][::skip_sph, ::skip_sph], v[1][::skip_sph, ::skip_sph], angles='xy', units="width", scale=scale_sph)
    plt.quiver(lon[2, :, 0:NX//2][::skip_sph, ::skip_sph], lat[2, :, 0:NX//2][::skip_sph, ::skip_sph], u[2, :, 0:NX//2][::skip_sph, ::skip_sph], v[2, :, 0:NX//2][::skip_sph, ::skip_sph], angles='xy', units="width", scale=scale_sph)
    plt.quiver(lon[2, :, NX//2:][::skip_sph, ::skip_sph]-360, lat[2, :, NX//2:][::skip_sph, ::skip_sph], u[2, :, NX//2:][::skip_sph, ::skip_sph], v[2, :, NX//2:][::skip_sph, ::skip_sph], angles='xy', units="width", scale=scale_sph)
    plt.quiver(lon[3][::skip_sph, ::skip_sph]-360, lat[3][::skip_sph, ::skip_sph], u[3][::skip_sph, ::skip_sph], v[3][::skip_sph, ::skip_sph], angles='xy', units="width", scale=scale_sph)
    plt.quiver(lon[4][::skip_sph, ::skip_sph], lat[4][::skip_sph, ::skip_sph], u[4][::skip_sph, ::skip_sph], v[4][::skip_sph, ::skip_sph], angles='xy', units="width", scale=scale_sph)
    Q = plt.quiver(lon[5][::skip_sph, ::skip_sph], lat[5][::skip_sph, ::skip_sph], u[5][::skip_sph, ::skip_sph], v[5][::skip_sph, ::skip_sph], angles='xy', units="width", scale=scale_sph)
    qk = plt.quiverkey(Q, 0.7, 0.9, wind, f"{wind}" + r'$ \frac{m}{s}$', labelpos='E', coordinates='figure')
            
    plt.savefig("../graphs/wind/spherical.png", dpi=DPI)
    plt.close()

    ### plot cartopy wind
    u = np.loadtxt("../outputs/u_lon_lat/u_lon_lat_0.txt").reshape(6, NX, NY).flatten()
    v = np.loadtxt("../outputs/v_lon_lat/v_lon_lat_0.txt").reshape(6, NX, NY).flatten()
    lon = np.loadtxt("../outputs/grids/lon.txt").reshape(6, NX, NY).flatten() * 180. / np.pi
    lat = np.loadtxt("../outputs/grids/lat.txt").reshape(6, NX, NY).flatten() * 180. / np.pi

    map_proj = ccrs.PlateCarree(central_longitude=0.)
    data_crs = ccrs.PlateCarree()
    fig, ax = plt.subplots(nrows=1,ncols=1,figsize=(12,8),subplot_kw={'projection': map_proj},dpi=DPI)

    x, y, _ = map_proj.transform_points(data_crs, lon, lat).T
    mask = np.invert(np.logical_or(np.isinf(x), np.isinf(y)))
    x = np.compress(mask, x)
    y = np.compress(mask, y)

    Q = ax.quiver(x[::skip_car], y[::skip_car], u[mask][::skip_car], v[mask][::skip_car], angles='xy', units="width", scale=scale_car)
    ax.quiverkey(Q, 0.85, 0.85, wind, f"{wind}" + r'$ \frac{m}{s}$', labelpos='E', coordinates='figure')

    plt.savefig("../graphs/wind/spherical_cartopy.png", dpi=DPI)
    plt.close()



def plotSphereCartopyZeta(t):
    u = np.loadtxt(f"../outputs/u/u_{t*LEAP}.txt").reshape(6, NX, NY)
    v = np.loadtxt(f"../outputs/v/v_{t*LEAP}.txt").reshape(6, NX, NY)
    x, y = np.loadtxt("../outputs/grids/x.txt").reshape(6, NX, NY), np.loadtxt("../outputs/grids/y.txt").reshape(6, NX, NY)

    zeta = (((v[:, :, 2:] - v[:, :, :-2]) / ((x[:, :, 2:] - x[:, :, :-2]) / 2))[:, 1:-1, :] - ((u[:, 2:, :] - u[:, :-2, :]) / ((y[:, 2:, :] - y[:, :-2, :]) / 2))[:, :, 1:-1]).flatten()

    lon = np.loadtxt("../outputs/grids/lon.txt").reshape(6, NX, NY)[:, 1:-1, 1:-1].flatten() * 180 / np.pi
    lat = np.loadtxt("../outputs/grids/lat.txt").reshape(6, NX, NY)[:, 1:-1, 1:-1].flatten() * 180 / np.pi

    map_proj = ccrs.PlateCarree(central_longitude=0.)
    data_crs = ccrs.PlateCarree()
    fig, ax = plt.subplots(nrows=1,ncols=1,figsize=(12,8),subplot_kw={'projection': map_proj},dpi=DPI)

    pos = ax.get_position()
    pset = [0.07, 0.05, 0.8, 0.9]
    ax.set_position(pset)

    x, y, _ = map_proj.transform_points(data_crs, lon, lat).T
    mask = np.invert(np.logical_or(np.isinf(x), np.isinf(y)))
    x = np.compress(mask, x)
    y = np.compress(mask, y)

    cs = ax.tricontourf(x, y, zeta[mask], extend="both", cmap=cmap, levels=np.linspace(left_zeta, right_zeta, split_zeta))
    ax.set_global()
    cb_ax1 = fig.add_axes([0.92, 0.16, 0.012, 0.67])
    fig.colorbar(cs, cax=cb_ax1, ticks=np.linspace(left_zeta, right_zeta, split_zeta))

    ax.coastlines(resolution='110m',color='k', lw=0.2, zorder=13)
    gl = ax.gridlines(draw_labels=True)

    ax.set_title(f"t = {t * LEAP * DT / 60} min", fontsize=fs)
    
    plt.savefig(f"../graphs/zeta/{int(t/LEAP/2)}.png", dpi=DPI)
    plt.close()

def plotOnCubeZeta(t):
    u = np.loadtxt(f"../outputs/u/u_{t*LEAP}.txt").reshape(6, NX, NY)
    v = np.loadtxt(f"../outputs/v/v_{t*LEAP}.txt").reshape(6, NX, NY)
    x, y = np.loadtxt("../outputs/grids/x.txt").reshape(6, NX, NY), np.loadtxt("../outputs/grids/y.txt").reshape(6, NX, NY)

    val = (((v[:, :, 2:] - v[:, :, :-2]) / ((x[:, :, 2:] - x[:, :, :-2]) / 2))[:, 1:-1, :] - ((u[:, 2:, :] - u[:, :-2, :]) / ((y[:, 2:, :] - y[:, :-2, :]) / 2))[:, :, 1:-1])


    fig = plt.figure(figsize=(18,10))
    ax5 = fig.add_subplot(3,4,2)
    ax4 = fig.add_subplot(3,4,5)
    ax1 = fig.add_subplot(3,4,6)
    ax2 = fig.add_subplot(3,4,7)
    ax3 = fig.add_subplot(3,4,8)
    ax6 = fig.add_subplot(3,4,10)
    ax5.set_title(f"t = {t * LEAP * DT / 60} min", fontsize=fs)

    cs1 = ax1.contourf(x[0, 1:-1, 1:-1], y[0, 1:-1, 1:-1], val[0], levels=np.linspace(left_zeta, right_zeta, split_zeta), extend='both', cmap=cmap)
    ax2.contourf(x[1, 1:-1, 1:-1], y[1, 1:-1, 1:-1], val[1], levels=np.linspace(left_zeta, right_zeta, split_zeta), extend='both', cmap=cmap)
    ax3.contourf(x[2, 1:-1, 1:-1], y[2, 1:-1, 1:-1], val[2], levels=np.linspace(left_zeta, right_zeta, split_zeta), extend='both', cmap=cmap)
    ax4.contourf(x[3, 1:-1, 1:-1], y[3, 1:-1, 1:-1], val[3], levels=np.linspace(left_zeta, right_zeta, split_zeta), extend='both', cmap=cmap)
    ax5.contourf(x[4, 1:-1, 1:-1], y[4, 1:-1, 1:-1], val[4], levels=np.linspace(left_zeta, right_zeta, split_zeta), extend='both', cmap=cmap)
    ax6.contourf(x[5, 1:-1, 1:-1], y[5, 1:-1, 1:-1], val[5], levels=np.linspace(left_zeta, right_zeta, split_zeta), extend='both', cmap=cmap)
    cb_ax1 = fig.add_axes([0.9235, 0.1, 0.015, 0.78])
    fig.colorbar(cs1, cax=cb_ax1, ticks=np.linspace(left_zeta, right_zeta, split_zeta))

    plt.savefig(f"../graphs/zeta/{int(t/LEAP/2)}.png", dpi=DPI)
    plt.close()
    return