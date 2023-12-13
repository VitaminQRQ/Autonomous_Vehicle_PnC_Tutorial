import math

from . import data_struct as struct
from . import utils as utils

class Vehicle:
    def __init__(self, L, initial_state: struct.State):
        self.L_ = L
        self.x_ = initial_state.pose_.x_
        self.y_ = initial_state.pose_.y_
        self.theta_ = utils.normalize_angle(initial_state.pose_.theta_)
        self.v_ = initial_state.v_

    def update_state(self, delta, a, dt):
        delta = min(delta,  35 * math.pi / 180)
        delta = max(delta, -35 * math.pi / 180)
        
        a = min(a,  5)
        a = max(a, -5)
        
        self.x_ += self.v_ * math.cos(self.theta_) * dt
        self.y_ += self.v_ * math.sin(self.theta_) * dt
        self.theta_ += self.v_ / self.L_ * math.tan(delta) * dt
        self.theta_ = utils.normalize_angle(self.theta_)
        self.v_ += a * dt
        self.v_ = max(0, self.v_)
        
        kappa = math.tan(delta) / self.L_
        
        new_pose  = struct.Transform(self.x_, self.y_, self.theta_)
        new_state = struct.State(new_pose, self.v_, a, kappa)
         
        return new_state