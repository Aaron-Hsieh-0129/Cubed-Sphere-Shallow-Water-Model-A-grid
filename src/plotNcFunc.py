import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import cartopy.util as cutil
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import xarray as xr


######## Should Be Tuned ########### 
left, right, split = 8000, 10500, 21
wind = 20
skip_car = 55
scale_car = 4000
skip = 10
scale = 1000
skip_sph = 12
scale_sph = 3000

left_zeta, right_zeta, split_zeta = -16 * 10 ** (-5), 16 * 10 ** (-5), 17

DT = 45
DX = DY = 0.5

LEAP = 400
######## Should Be Tuned ###########

NX = NY = int(90 / DX)

cmap = cm.viridis
wind_color = "#5c5959"
DPI = 100
fs = 15

def plotOnCubeWindMul(t):
    path = f"../outputs/nc/{t}.nc"
    x = xr.open_dataset("../outputs/nc/grid.nc")['x_local'][:, 1:-1, 1:-1].to_numpy(), 
    y = xr.open_dataset("../outputs/nc/grid.nc")['y_local'][:, 1:-1, 1:-1].to_numpy()
    val = xr.open_dataset(path)['h'][:, 1:-1, 1:-1].to_numpy()
    u = xr.open_dataset(path)['u'][:, 1:-1, 1:-1].to_numpy()
    v = xr.open_dataset(path)['v'][:, 1:-1, 1:-1].to_numpy()

    fig = plt.figure(figsize=(18,10), dpi=DPI)
    ax5 = fig.add_subplot(3,4,2)
    ax4 = fig.add_subplot(3,4,5)
    ax1 = fig.add_subplot(3,4,6)
    ax2 = fig.add_subplot(3,4,7)
    ax3 = fig.add_subplot(3,4,8)
    ax6 = fig.add_subplot(3,4,10)
    ax5.set_title(f"t = {t * DT / 60} min", fontsize=fs)

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
    
    plt.savefig(f"../graphs/h/curvilinear/{int(t/LEAP)}.png", dpi=DPI)
    plt.close()
    return

def plotSphereWindCartopy(t):
    path = f"../outputs/nc/{t}.nc"
    lon = xr.open_dataset(f"../outputs/nc/grid.nc")['lon_sphere'][:, 1:-1, 1:-1].to_numpy() * 180 / np.pi
    lon = lon.flatten()

    lat = xr.open_dataset(f"../outputs/nc/grid.nc")['lat_sphere'][:, 1:-1, 1:-1].to_numpy().flatten() * 180 / np.pi
    h = xr.open_dataset(path)['h'][:, 1:-1, 1:-1].to_numpy().flatten()
    u = xr.open_dataset(path)['u_lonlat'][:, 1:-1, 1:-1].to_numpy().flatten()
    v = xr.open_dataset(path)['v_lonlat'][:, 1:-1, 1:-1].to_numpy().flatten()

    map_proj = ccrs.PlateCarree(central_longitude=0.)
    data_crs = ccrs.PlateCarree()
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(12,8), subplot_kw={'projection': map_proj}, dpi=DPI)

    pos = ax.get_position()
    pset = [0.07, 0.05, 0.8, 0.9]
    ax.set_position(pset)

    x, y, _ = map_proj.transform_points(data_crs, lon, lat).T
    mask = np.invert(np.logical_or(np.isinf(x), np.isinf(y)))
    x = np.compress(mask, x)
    y = np.compress(mask, y)

    cs = ax.tricontourf(x, y, h[mask], extend="both", cmap=cmap, levels=np.linspace(left, right, split))
    ax.set_global()
    cb_ax1 = fig.add_axes([0.92, 0.16, 0.012, 0.67])
    fig.colorbar(cs, cax=cb_ax1, ticks=np.linspace(left, right, split))

    Q = ax.quiver(x[::skip_car], y[::skip_car], u[mask][::skip_car], v[mask][::skip_car], angles='xy', units="width", scale=scale_sph, color=wind_color)
    ax.quiverkey(Q, 0.85, 0.85, wind, f"{wind}" + r'$ \frac{m}{s}$', labelpos='E', coordinates='figure')

    # ax.coastlines(resolution='110m',color='k', lw=0.2, zorder=13)
    gl = ax.gridlines(draw_labels=True)

    ax.set_title(f"t = {t * DT / 60} min", fontsize=fs)
    
    plt.savefig(f"../graphs/h/sphere_cartopy/{int(t/LEAP)}.png", dpi=DPI)
    plt.close()


def plotSphereCartopyZeta(t):
    u = xr.open_dataset(f"../outputs/nc/{t}.nc")['u'][:, 1:-1, 1:-1].to_numpy()
    v = xr.open_dataset(f"../outputs/nc/{t}.nc")['v'][:, 1:-1, 1:-1].to_numpy()
    x = xr.open_dataset("../outputs/nc/grid.nc")['x_local'][:, 1:-1, 1:-1].to_numpy()
    y = xr.open_dataset("../outputs/nc/grid.nc")['y_local'][:, 1:-1, 1:-1].to_numpy()

    zeta = (((v[:, 2:, :] - v[:, :-2, :]) / ((x[:, 2:, :] - x[:, :-2, :]) / 2))[:, :, 1:-1] - ((u[:, :, 2:] - u[:, :, :-2]) / ((y[:, :, 2:] - y[:, :, :-2]) / 2))[:, 1:-1, :]).flatten()

    lon = xr.open_dataset(f"../outputs/nc/grid.nc")['lon_sphere'][:, 2:-2, 2:-2].to_numpy().flatten() * 180 / np.pi
    lat = xr.open_dataset(f"../outputs/nc/grid.nc")['lat_sphere'][:, 2:-2, 2:-2].to_numpy().flatten() * 180 / np.pi

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

    # ax.coastlines(resolution='110m',color='k', lw=0.2, zorder=13)
    gl = ax.gridlines(draw_labels=True)

    ax.set_title(f"t = {t * DT / 60} min", fontsize=fs)
    
    plt.savefig(f"../graphs/zeta/sphere_cartopy/{int(t/LEAP)}.png", dpi=DPI)
    plt.close()

def plotOnCubeZeta(t):
    u = xr.open_dataset(f"../outputs/nc/{t}.nc")['u'][:, 1:-1, 1:-1].to_numpy()
    v = xr.open_dataset(f"../outputs/nc/{t}.nc")['v'][:, 1:-1, 1:-1].to_numpy()
    x = xr.open_dataset("../outputs/nc/grid.nc")['x_local'][:, 1:-1, 1:-1].to_numpy(), 
    y = xr.open_dataset("../outputs/nc/grid.nc")['y_local'][:, 1:-1, 1:-1].to_numpy()

    val = (((v[:, 2:, :] - v[:, :-2, :]) / ((x[:, 2:, :] - x[:, :-2, :]) / 2))[:, :, 1:-1] - ((u[:, :, 2:] - u[:, :, :-2]) / ((y[:, :, 2:] - y[:, :, :-2]) / 2))[:, 1:-1, :])


    fig = plt.figure(figsize=(18,10))
    ax5 = fig.add_subplot(3,4,2)
    ax4 = fig.add_subplot(3,4,5)
    ax1 = fig.add_subplot(3,4,6)
    ax2 = fig.add_subplot(3,4,7)
    ax3 = fig.add_subplot(3,4,8)
    ax6 = fig.add_subplot(3,4,10)
    ax5.set_title(f"t = {t * DT / 60} min", fontsize=fs)

    x = x[:, 1:-1, 1:-1]
    y = y[:, 1:-1, 1:-1]
    cs1 = ax1.contourf(x[0], y[0], val[0], levels=np.linspace(left_zeta, right_zeta, split_zeta), extend='both', cmap=cmap)
    ax2.contourf(x[1], y[1], val[1], levels=np.linspace(left_zeta, right_zeta, split_zeta), extend='both', cmap=cmap)
    ax3.contourf(x[2], y[2], val[2], levels=np.linspace(left_zeta, right_zeta, split_zeta), extend='both', cmap=cmap)
    ax4.contourf(x[3], y[3], val[3], levels=np.linspace(left_zeta, right_zeta, split_zeta), extend='both', cmap=cmap)
    ax5.contourf(x[4], y[4], val[4], levels=np.linspace(left_zeta, right_zeta, split_zeta), extend='both', cmap=cmap)
    ax6.contourf(x[5], y[5], val[5], levels=np.linspace(left_zeta, right_zeta, split_zeta), extend='both', cmap=cmap)
    cb_ax1 = fig.add_axes([0.9235, 0.1, 0.015, 0.78])
    fig.colorbar(cs1, cax=cb_ax1, ticks=np.linspace(left_zeta, right_zeta, split_zeta))

    plt.savefig(f"../graphs/zeta/curvilinear/{int(t/LEAP)}.png", dpi=DPI)
    plt.close()
    return