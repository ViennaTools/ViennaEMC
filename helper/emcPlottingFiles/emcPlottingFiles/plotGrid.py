import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm
from matplotlib.colors import Normalize

# plot a grid structure (either 3d or 2d)
def plotGrid(data, zLabel="", addColorBar=True):
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    if data.ndim == 3: # 3d case ( TODO find better way to plot)
        Z, Y, X = np.mgrid[0:data.shape[0], 0:data.shape[1], 0:data.shape[2]]
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        scat = ax.scatter(X, Y, Z, c=data.flatten('C'), alpha=0.5)
        if addColorBar:
            fig.colorbar(scat, shrink=0.5, aspect=5)
    else:
        x = np.linspace(0, data.shape[1], data.shape[1])
        y = np.linspace(0, data.shape[0], data.shape[0])
        X, Y = np.meshgrid(x, y)
        surf = ax.plot_surface(X, Y, data, cmap=cm.viridis,
                            linewidth=0, antialiased=False)
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel(zLabel)
        if addColorBar:
            fig.colorbar(surf)
    plt.tight_layout()
    

    
