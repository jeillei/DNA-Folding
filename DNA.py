import numpy as np
from utils import noise, rotational_steps
from physics import spring, electro, apply_dihedral_force

class Base:
    def __init__(self, type, position):
        self.type = type
        self.charge = -1
        self.position = np.array(position, dtype=float)
        self.velocity = np.zeros(3)
        self.force = np.zeros(3)

class Strand:
    pair_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

    def __init__(self, sequence, angle_variation=0.05):
        self.bases = []
        current_pos = np.zeros(3)
        stepper = rotational_steps(angle_variation=angle_variation)

        for b_char in sequence:
            base = Base(b_char, current_pos)
            self.bases.append(base)
            step = stepper.step() + noise() / 10
            current_pos += step

    def complement(self):
        return ''.join(self.pair_dict[base.type] for base in self.bases)

    def backbone_forces(self, k=10, dflt=0.34):
        for i in range(len(self.bases) - 1):
            b1, b2 = self.bases[i], self.bases[i + 1]
            f = spring(b1.position, b2.position, k, dflt)
            b1.force += f
            b2.force -= f

class DNA:
    def __init__(self, strand: Strand, strand_dist=2.0, angle_variation=0.05):
        self.strand1 = strand
        comp_seq = strand.complement()

        base_dir = strand.bases[1].position - strand.bases[0].position
        base_dir /= np.linalg.norm(base_dir)
        arbitrary = np.array([0, 0, 1]) if not np.allclose(np.cross(base_dir, [0, 0, 1]), 0) else [0, 1, 0]
        offset_axis = np.cross(base_dir, arbitrary)
        offset_axis /= np.linalg.norm(offset_axis)
        fixed_offset = offset_axis * strand_dist

        self.strand2 = Strand.__new__(Strand)
        self.strand2.bases = []
        stepper = rotational_steps(angle_variation=angle_variation)
        current_pos = strand.bases[0].position + fixed_offset

        for b_char in comp_seq:
            base = Base(b_char, current_pos)
            self.strand2.bases.append(base)
            step = stepper.step() + noise() / 10
            current_pos += step

    def base_pair_forces(self, k=5, dflt=2.0):
        pair_strength = {
                ('A', 'T'): 2.0,
                ('T', 'A'): 2.0,
                ('C', 'G'): 3.0,
                ('G', 'C'): 3.0
            }

        for b1, b2 in zip(self.strand1.bases, self.strand2.bases):
            k = pair_strength.get((b1.type, b2.type), 0.0)
            f = spring(b1.position, b2.position, k=k, dflt=2.0)
            b1.force += f
            b2.force -= f


    def compute_electrostatics(self, k_e=-6):
        all_bases = self.strand1.bases + self.strand2.bases
        for i in range(len(all_bases)):
            for j in range(i + 1, len(all_bases)):
                if abs(i - j) <= 1:
                    continue
                f = electro(all_bases[i].position, all_bases[j].position, k_e=k_e)
                all_bases[i].force += f
                all_bases[j].force -= f

    def stacking_forces(self, k=20, dflt=0.34):
        for strand in [self.strand1, self.strand2]:
            for i in range(len(strand.bases) - 1):
                b1, b2 = strand.bases[i], strand.bases[i + 1]
                f = spring(b1.position, b2.position, k, dflt)
                b1.force += f
                b2.force -= f

    def apply_dihedral_forces(self):
        for strand in [self.strand1, self.strand2]:
            for i in range(len(strand.bases) - 3):
                apply_dihedral_force(*strand.bases[i:i+4])
