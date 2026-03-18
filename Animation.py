import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

try:
    df = pd.read_csv('orbit.csv')
except FileNotFoundError:
    print("Error: orbit.csv not found. Please run your C simulation first!")
    exit()

EV_CONV = 1.60217663e-19

ex = df["ex"].to_numpy()
ey = df["ey"].to_numpy()
px = df["px"].to_numpy()
py = df["py"].to_numpy()
ke = df["ke"].to_numpy()
pe = df["pe"].to_numpy()
te = df["te"].to_numpy()
steps = df["step"].to_numpy()
te = te/EV_CONV
fig , (ax1,ax2) = plt.subplots(1,2)

ax1.set_xlim(np.min(ex) - 1e-11,np.max(ex)+1e-11)
ax1.set_ylim(np.min(ey)-1e-11,np.max(ey)+1e-11)
ax1.set_aspect("equal")
ax1.set_title("Bohr's model electron and proton path")
ax1.axis("off")

proton_plot, = ax1.plot([],[],'ro',markersize = 10,label = 'Proton')
electron_plot, = ax1.plot([],[],'bo',markersize=5,label='Electron')
electron_trail, = ax1.plot([], [], 'b-', alpha=0.2, linewidth=1)

ax2.set_title("Total Energy")
ax2.set_xlabel("Steps")
ax2.set_xlim(steps[0], steps[-1])
te_mean = np.mean(te)
ax2.set_ylim(te_mean - 0.05, te_mean + 0.05)
ax2.set_ylabel("Energy(ev)")
ax2.grid(True)
te_plot, = ax2.plot([],[],'g-',linewidth = 2)

def update(frame):
    electron_plot.set_data([ex[frame]],[ey[frame]])
    proton_plot.set_data([px[frame]],[py[frame]])
    start = max(0, frame - 200)
    electron_trail.set_data(ex[start:frame], ey[start:frame])

    te_plot.set_data(steps[:frame],te[:frame])

    return electron_plot, proton_plot, te_plot ,electron_trail

ani = FuncAnimation(fig,update,frames=range(0,len(df),20),interval=30,blit = True)

plt.tight_layout()
plt.show()
