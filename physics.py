import numpy as np

def spring(pos1, pos2, k, dflt):
    vec = pos2-pos1
    dist = np.linalg.norm(vec)
    direction = vec / dist if dist != 0 else np.zeros(3)
    return k * (dist-dflt) * direction

def electro(pos1, pos2, q1=-1, q2=-1, k_e=6, lambda_D=1.0, min_dist=0.2, max_force=1000):
    vec = pos2 - pos1
    dist = np.linalg.norm(vec)

    if dist < min_dist:
        dist = min_dist

    direction = vec / dist
    screened = np.exp(-dist / lambda_D)
    force = k_e * q1 * q2 * screened / (dist**2)
    
    if abs(force) > max_force:
        force = max_force * np.sign(force)

    return force * direction

def dihedral_angle(a, b, c, d):
    b1 = b.position - a.position
    b2 = c.position - b.position
    b3 = d.position - c.position

    n1 = np.cross(b1, b2)
    n2 = np.cross(b2, b3)

    n1 /= np.linalg.norm(n1) + 1e-8
    n2 /= np.linalg.norm(n2) + 1e-8
    b2_unit = b2 / (np.linalg.norm(b2) + 1e-8)

    x = np.dot(n1, n2)
    y = np.dot(np.cross(n1, n2), b2_unit)

    return np.arctan2(y, x)

def apply_dihedral_force(a, b, c, d, k_dihedral=0.4, preferred_phi=np.pi / 5):
    phi = dihedral_angle(a, b, c, d)
    delta_phi = phi - preferred_phi
    torque = -k_dihedral * delta_phi

    axis = c.position - b.position
    axis /= np.linalg.norm(axis) + 1e-8

    perp_a = np.cross(axis, b.position - a.position)
    perp_d = np.cross(axis, d.position - c.position)

    perp_a /= np.linalg.norm(perp_a) + 1e-8
    perp_d /= np.linalg.norm(perp_d) + 1e-8

    force_mag = torque / 1.0

    f_a = force_mag * perp_a
    f_d = -force_mag * perp_d

    a.force += f_a
    b.force -= f_a * 0.5
    c.force -= f_d * 0.5
    d.force += f_d

def brownian(base, strength= 0):
    motion = np.random.normal(0, strength, size = 3)
    base.force += motion