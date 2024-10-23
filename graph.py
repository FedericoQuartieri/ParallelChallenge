import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Caricare i dati dai file CSV
parallel_df = pd.read_csv('tempi_3d_interpolated_parallel.csv')
serial_df = pd.read_csv('tempi_3d_interpolated_serial.csv')

# Estrarre i valori da entrambe le tabelle
cutoff_values = np.unique(parallel_df['Cutoff'].values)
size_values = np.unique(parallel_df['Array Size'].values)

# Creare griglie per i valori di cutoff e size
X, Y = np.meshgrid(cutoff_values, size_values)

# Creare griglie per i valori di tempo parallelo e seriale
Z_parallel = parallel_df.pivot(index='Array Size', columns='Cutoff', values='Parallel Time').values
Z_serial = serial_df.pivot(index='Array Size', columns='Cutoff', values='Serial Time').values

# Calcolare lo speedup (Serial Time / Parallel Time) in ogni punto
Z_speedup = np.divide(Z_serial, Z_parallel, out=np.zeros_like(Z_serial), where=Z_parallel!=0)

# Creare la figura e l'asse 3D
fig = plt.figure(figsize=(10, 7))
ax = fig.add_subplot(111, projection='3d')

# Grafico della superficie dello speedup
surf_speedup = ax.plot_surface(X, Y, Z_speedup, cmap='coolwarm', alpha=0.8)

# Aggiungere etichette e titolo
ax.set_xlabel('Cutoff')
ax.set_ylabel('Array Size')
ax.set_zlabel('Speedup (Serial/Parallel)')
ax.set_title('3D Plot of Speedup (Serial Time / Parallel Time)')

# Aggiungere la barra dei colori per lo speedup
fig.colorbar(surf_speedup, ax=ax, shrink=0.5, aspect=5, label="Speedup")

# Mostrare il grafico
plt.savefig('plot_3d_cutoff.png', format='png', dpi=300)
