'''
用于存储坐标信息的结构体
'''
class Transform:
    def __init__(self, x: float, y: float, theta: float):
        self.x_ = x
        self.y_ = y
        self.theta_ = theta
        
'''
用于存储静态障碍物信息
'''
class StaticObstacle:
    def __init__(self, pose: Transform, radius: float):
        self.pose_   = pose
        self.radius_ = radius
        
'''
用于存储 Frenet 坐标系下的车辆状态信息
'''
class FrenetState:
    def __init__(self, s, s_dot, s_ddot, d, d_dot, d_ddot):
        self.s_      = s
        self.s_dot_  = s_dot
        self.s_ddot_ = s_ddot
        self.d_      = d
        self.d_dot_  = d_dot
        self.d_ddot_ = d_ddot

'''
用于存储车辆状态信息的结构体
'''
class State:
    def __init__(self, pose: Transform, v: float, a: float, kappa: float):
        self.pose_ = pose
        self.v_    = v
        self.a_    = a
        self.kappa_ = kappa