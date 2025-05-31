import matplotlib.pyplot as plt
import numpy as np
import matplotlib.animation as animation
from utils import update_pos
from DNA import DNA

def animate_dna(dna: DNA, steps=100, dt=0.03, save_path="outputs/dna_folding.mp4"):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    line1, = ax.plot([], [], [], 'ro-')
    line2, = ax.plot([], [], [], 'bo-')
    lines = [ax.plot([], [], [], 'k--')[0] for _ in dna.strand1.bases]

    def init():
        ax.set_xlim(-10, 10)
        ax.set_ylim(-10, 10)
        ax.set_zlim(-15, 15)
        return [line1, line2] + lines

    def update(frame):
        for base in dna.strand1.bases + dna.strand2.bases:
            base.force = np.zeros(3)

        dna.strand1.backbone_forces()
        dna.strand2.backbone_forces()
        dna.base_pair_forces()
        dna.compute_electrostatics()
        dna.stacking_forces()
        # dna.apply_dihedral_forces()

        update_pos(dna.strand1.bases + dna.strand2.bases, dt=dt)

      
        offset = np.array([20, 0, 0])  

        xs1 = [b.position[0] - offset[0] for b in dna.strand1.bases]
        ys1 = [b.position[1] - offset[1] for b in dna.strand1.bases]
        zs1 = [b.position[2] - offset[2] for b in dna.strand1.bases]

        xs2 = [b.position[0] - offset[0] for b in dna.strand2.bases]
        ys2 = [b.position[1] - offset[1] for b in dna.strand2.bases]
        zs2 = [b.position[2] - offset[2] for b in dna.strand2.bases]

        line1.set_data(ys1, zs1)          
        line1.set_3d_properties(xs1)    

        line2.set_data(ys2, zs2)
        line2.set_3d_properties(xs2)


        for i, (b1, b2) in enumerate(zip(dna.strand1.bases, dna.strand2.bases)):
            x_pair = [b1.position[1] - offset[1], b2.position[1] - offset[1]]
            y_pair = [b1.position[2] - offset[2], b2.position[2] - offset[2]]
            z_pair = [b1.position[0] - offset[0], b2.position[0] - offset[0]]
            lines[i].set_data(x_pair, y_pair)
            lines[i].set_3d_properties(z_pair)


        return [line1, line2] + lines

    ani = animation.FuncAnimation(fig, update, init_func=init, frames=steps, interval=50, blit=False)
    # ani.save(save_path, writer='ffmpeg', fps=30)
    plt.show()
