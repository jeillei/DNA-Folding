import numpy as np

def noise():
    vec = np.random.randn(3)
    vec /= np.linalg.norm(vec)
    return vec

def rand_unit_vec():
    vec = np.random.randn(3)
    return vec / np.linalg.norm(vec)

class rotational_steps:
    def __init__(self, initial_dir=np.array([1.0, 0, 0]), angle_variation=0.1):
        self.direction = initial_dir / np.linalg.norm(initial_dir)
        self.angle_var = angle_variation
        self.rng = np.random.default_rng()

    def step(self, step_size=0.34):
        angle = self.rng.normal(0, self.angle_var)
        axis = self.rng.standard_normal(3)
        axis -= axis.dot(self.direction) * self.direction / np.linalg.norm(self.direction)**2
        axis /= np.linalg.norm(axis)

        cos_a = np.cos(angle)
        sin_a = np.sin(angle)
        cross = np.cross(axis, self.direction)
        dot = np.dot(axis, self.direction)
        new_dir = (
            self.direction * cos_a +
            cross * sin_a +
            axis * dot * (1 - cos_a)
        )
        new_dir /= np.linalg.norm(new_dir)
        self.direction = new_dir
        return new_dir * step_size

def update_pos(bases, dt=0.03, damping=0.98):
    for base in bases:
        base.velocity += base.force * dt
        base.position += base.velocity * dt
        base.velocity *= damping
        base.force = np.zeros(3)
