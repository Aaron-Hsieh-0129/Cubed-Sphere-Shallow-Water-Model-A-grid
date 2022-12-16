import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

DX = DY = 2
NX = NY = int(90 / DX)
DT = 360
LEAP = 10
cmap = cm.viridis

def plotOnCubeWindMul(t):
    left, right, split = 0, 2000, 21

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
    ax5.set_title(f"t = {t * 10 * DT / 60} min", fontsize=14)

    cs1 = ax1.contourf(x[0], y[0], val[0], levels=np.linspace(left, right, split), extend='both', cmap=cmap)
    ax2.contourf(x[1], y[1], val[1], levels=np.linspace(left, right, split), extend='both', cmap=cmap)
    ax3.contourf(x[2], y[2], val[2], levels=np.linspace(left, right, split), extend='both', cmap=cmap)
    ax4.contourf(x[3], y[3], val[3], levels=np.linspace(left, right, split), extend='both', cmap=cmap)
    ax5.contourf(x[4], y[4], val[4], levels=np.linspace(left, right, split), extend='both', cmap=cmap)
    ax6.contourf(x[5], y[5], val[5], levels=np.linspace(left, right, split), extend='both', cmap=cmap)
    cb_ax1 = fig.add_axes([0.9235, 0.1, 0.015, 0.78])
    fig.colorbar(cs1, cax=cb_ax1, ticks=np.linspace(left, right, split))

    Q = ax1.quiver(x[0, ::2, ::2], y[0, ::2, ::2], u[0, ::2, ::2], v[0, ::2, ::2], angles='xy', units="width", scale=100)
    ax2.quiver(x[1, ::2, ::2], y[1, ::2, ::2], u[1, ::2, ::2], v[1, ::2, ::2], angles='xy', units="width", scale=100)
    ax3.quiver(x[2, ::2, ::2], y[2, ::2, ::2], u[2, ::2, ::2], v[2, ::2, ::2], angles='xy', units="width", scale=100)
    ax4.quiver(x[3, ::2, ::2], y[3, ::2, ::2], u[3, ::2, ::2], v[3, ::2, ::2], angles='xy', units="width", scale=100)
    ax5.quiver(x[4, ::2, ::2], y[4, ::2, ::2], u[4, ::2, ::2], v[4, ::2, ::2], angles='xy', units="width", scale=100)
    ax6.quiver(x[5, ::2, ::2], y[5, ::2, ::2], u[5, ::2, ::2], v[5, ::2, ::2], angles='xy', units="width", scale=100)

    ax1.quiverkey(Q, 0.7, 0.9, 10, r'$10 \frac{m}{s}$', labelpos='E', coordinates='figure')  
    ax2.quiverkey(Q, 0.7, 0.9, 10, r'$10 \frac{m}{s}$', labelpos='E', coordinates='figure')  
    ax3.quiverkey(Q, 0.7, 0.9, 10, r'$10 \frac{m}{s}$', labelpos='E', coordinates='figure')  
    ax4.quiverkey(Q, 0.7, 0.9, 10, r'$10 \frac{m}{s}$', labelpos='E', coordinates='figure')  
    ax5.quiverkey(Q, 0.7, 0.9, 10, r'$10 \frac{m}{s}$', labelpos='E', coordinates='figure')  
    ax6.quiverkey(Q, 0.7, 0.9, 10, r'$10 \frac{m}{s}$', labelpos='E', coordinates='figure')   
    
    plt.savefig(f"../graphs/h/curvilinear/{int(t/10)}.png", dpi=150)
    plt.close()
    return

def plotOnSphereWindMul(t):
    x, y = np.loadtxt("../outputs/grids/lon.txt").reshape(6, NX, NY), np.loadtxt("../outputs/grids/lat.txt").reshape(6, NX, NY)
    val = np.loadtxt(f"../outputs/h/h_{t*LEAP}.txt").reshape(6, NX, NY)
    x = x * 180 / np.pi
    y = y * 180 / np.pi

    u = np.loadtxt(f"../outputs/u/u_{t*LEAP}.txt").reshape(6, NX, NY)
    v = np.loadtxt(f"../outputs/v/v_{t*LEAP}.txt").reshape(6, NX, NY)

    left, right, split = 0, 2000, 21
    plt.figure(figsize=(18,8))
    plt.xlabel("LON")
    plt.ylabel("LAT")
    plt.title(f"t = {t * 10 * DT / 60} min", fontsize=14)

    # cmap = cm.twilight_shifted
    plt.contourf(x[0], y[0], val[0], levels=np.linspace(left, right, split), extend='both', cmap=cmap)
    plt.contourf(x[1], y[1], val[1], levels=np.linspace(left, right, split), extend='both', cmap=cmap)
    plt.contourf(x[2, :, 0:NX//2], y[2, :, 0:NX//2], val[2, :, 0:NX//2], levels=np.linspace(left, right, split), extend='both', cmap=cmap)
    plt.contourf(x[2, :, NX//2:]-360, y[2, :, NX//2:], val[2, :, NX//2:], levels=np.linspace(left, right, split), extend='both', cmap=cmap)
    plt.contourf(x[3]-360, y[3], val[3], levels=np.linspace(left, right, split), extend='both', cmap=cmap)
    plt.contourf(x[4], y[4], val[4], levels=np.linspace(left, right, split), extend='both', cmap=cmap)
    plt.contourf(x[5], y[5], val[5], levels=np.linspace(left, right, split), extend='both', cmap=cmap)
    cbar = plt.colorbar(pad=0.05)
    cbar.set_ticks(np.linspace(left, right, split))
    
    # plt.quiver(x[0][::2, ::2], y[0][::2, ::2], u[0][::2, ::2], v[0][::2, ::2], angles='xy', units="width", scale=100)
    # plt.quiver(x[1][::2, ::2], y[1][::2, ::2], u[1][::2, ::2], v[1][::2, ::2], angles='xy', units="width", scale=100)
    # plt.quiver(x[2, :, 0:NX//2][::2, ::2], y[2, :, 0:NX//2][::2, ::2], u[2, :, 0:NX//2][::2, ::2], v[2, :, 0:NX//2][::2, ::2], angles='xy', units="width", scale=100)
    # plt.quiver(x[2, :, NX//2:][::2, ::2]-360, y[2, :, NX//2:][::2, ::2], u[2, :, NX//2:][::2, ::2], v[2, :, NX//2:][::2, ::2], angles='xy', units="width", scale=100)
    # plt.quiver(x[3][::2, ::2]-360, y[3][::2, ::2], u[3][::2, ::2], v[3][::2, ::2], angles='xy', units="width", scale=100)
    # plt.quiver(x[4][::2, ::2], y[4][::2, ::2], u[4][::2, ::2], v[4][::2, ::2], angles='xy', units="width", scale=100)
    # Q = plt.quiver(x[5][::2, ::2], y[5][::2, ::2], u[5][::2, ::2], v[5][::2, ::2], angles='xy', units="width", scale=100)
    # qk = plt.quiverkey(Q, 0.7, 0.9, 0.5, r'$5 \frac{m}{s}$', labelpos='E', coordinates='figure')
        
    plt.savefig(f"../graphs/h/sphere/{int(t/10)}.png", dpi=150)
    plt.close()
    return

def plotOnSphereMul(t):
    x, y = np.loadtxt("../outputs/grids/x.txt").reshape(6, NX, NY), np.loadtxt("../outputs/grids/y.txt").reshape(6, NX, NY)
    val = np.loadtxt(f"../outputs/h/h_{t*LEAP}.txt").reshape(6, NX, NY)

    left, right, split = -1200, 1200, 13
    plt.figure(figsize=(18,8))
    plt.xlabel("LON")
    plt.ylabel("LAT")
    plt.title(f"t = {t * 10 * DT / 60} min", fontsize=14)

    plt.contourf(x[0], y[0], val[0], levels=np.linspace(left, right, split), extend='both', cmap=cmap)
    plt.contourf(x[1], y[1], val[1], levels=np.linspace(left, right, split), extend='both', cmap=cmap)
    plt.contourf(x[2, :, 0:NX//2], y[2, :, 0:NX//2], val[2, :, 0:NX//2], levels=np.linspace(left, right, split), extend='both', cmap=cmap)
    plt.contourf(x[2, :, NX//2:]-360, y[2, :, NX//2:], val[2, :, NX//2:], levels=np.linspace(left, right, split), extend='both', cmap=cmap)
    plt.contourf(x[3]-360, y[3], val[3], levels=np.linspace(left, right, split), extend='both', cmap=cmap)
    plt.contourf(x[4], y[4], val[4], levels=np.linspace(left, right, split), extend='both', cmap=cmap)
    plt.contourf(x[5], y[5], val[5], levels=np.linspace(left, right, split), extend='both', cmap=cmap)
    cbar = plt.colorbar(pad=0.05)
    cbar.set_ticks(np.linspace(left, right, split))

    plt.savefig(f"../graphs/h/sphere/{int(t/10)}.png", dpi=150)
    plt.close()
    return

def plotOnCubeMul(t):
    x, y = np.loadtxt("../outputs/grids/x.txt").reshape(6, NX, NY), np.loadtxt("../outputs/grids/y.txt").reshape(6, NX, NY)
    val = np.loadtxt(f"../outputs/h/h_{t*LEAP}.txt").reshape(6, NX, NY)

    fig = plt.figure(figsize=(18,10))
    ax5 = fig.add_subplot(3,4,2)
    ax4 = fig.add_subplot(3,4,5)
    ax1 = fig.add_subplot(3,4,6)
    ax2 = fig.add_subplot(3,4,7)
    ax3 = fig.add_subplot(3,4,8)
    ax6 = fig.add_subplot(3,4,10)
    ax5.set_title(f"t = {t * 10 * DT / 60} min", fontsize=14)

    left, right, split = 0, 2000, 13
    cs1 = ax1.contourf(x[0], y[0], val[0], levels=np.linspace(left, right, split), extend='both', cmap=cmap)
    ax2.contourf(x[1], y[1], val[1], levels=np.linspace(left, right, split), extend='both', cmap=cmap)
    ax3.contourf(x[2], y[2], val[2], levels=np.linspace(left, right, split), extend='both', cmap=cmap)
    ax4.contourf(x[3], y[3], val[3], levels=np.linspace(left, right, split), extend='both', cmap=cmap)
    ax5.contourf(x[4], y[4], val[4], levels=np.linspace(left, right, split), extend='both', cmap=cmap)
    ax6.contourf(x[5], y[5], val[5], levels=np.linspace(left, right, split), extend='both', cmap=cmap)
    cb_ax1 = fig.add_axes([0.9235, 0.1, 0.015, 0.78])
    fig.colorbar(cs1, cax=cb_ax1, ticks=np.linspace(left, right, split))

    plt.savefig(f"../graphs/h/curvilinear/{int(t/10)}.png", dpi=150)
    plt.close()
    return