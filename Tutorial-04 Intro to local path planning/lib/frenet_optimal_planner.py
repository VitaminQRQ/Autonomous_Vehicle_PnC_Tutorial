import copy
import math
import numpy as np
import matplotlib.pyplot as plt

from scipy.spatial import cKDTree
from typing import List, Tuple
from scipy.interpolate import BSpline

from . import utils
from . import data_struct as struct
from . import param_parser as param_parser
'''
用于补充参考线的航向、曲率信息，并转化至 frenet 坐标系
'''
class ReferencePath:
    def __init__(
        self,
        waypoint_x: np.mat, # 参考路径点的 x 坐标
        waypoint_y: np.mat, # 参考路径点的 y 坐标
        resolution: float
        ):       

        # 计算参考路径累计弧长
        waypoint_s = np.zeros_like(waypoint_x)
        for i in range(1, len(waypoint_x)):
            prev_x, prev_y = waypoint_x[i - 1], waypoint_y[i - 1]
            curr_x, curr_y = waypoint_x[i], waypoint_y[i]
            waypoint_s[i] = waypoint_s[i - 1] + np.sqrt((curr_x - prev_x) ** 2 + (curr_y - prev_y) ** 2)
        
        if len(waypoint_s) < 2:
            print('waypoints_s: ', waypoint_s)
            raise ValueError("Waypoint size less than 2")
        
        # 记录路径基本信息
        self.total_length_ = waypoint_s[-1, 0]
        
        # 使用三次样条曲线拟合
        spline_x = BSpline(np.ravel(waypoint_s), np.ravel(waypoint_x), 3)
        spline_y = BSpline(np.ravel(waypoint_s), np.ravel(waypoint_y), 3)
        
        # 初始化参考路径的分辨率
        self.sample_resolution_ = resolution
        self.interp_size_ = (np.ceil(self.total_length_ / self.sample_resolution_) + 1).astype(int)
        
        self.interp_s_ = np.linspace(0, self.total_length_, self.interp_size_)
        
        # 计算插值后的 x y theta kappa dkappa
        self.interp_x_      = np.mat(spline_x(self.interp_s_)).T
        self.interp_y_      = np.mat(spline_y(self.interp_s_)).T
        
        self.interp_theta_  = self.calculate_theta() 
        self.interp_kappa_  = self.calculate_kappa()
        self.interp_dkappa_ = self.calculate_dkappa()

        interp_xy_stack = np.hstack((self.interp_x_, self.interp_y_)) 
        self.kd_tree_ = cKDTree(interp_xy_stack)
        
    def calculate_theta(self):
        # 初始化航向角为空列表
        interp_theta = np.zeros((self.interp_size_, 1))
        
        # 如果路径点数量为 2，则无法使用中心差分
        if self.interp_size_ == 2:
            dy = self.interp_y_[1, 0] - self.interp_y_[1, 0]
            dx = self.interp_x_[1, 0] - self.interp_x_[0, 0]
            theta = np.arctan2(dy, dx)
            
            interp_theta[0, 0] = theta
            interp_theta[1, 0] = theta
            
            return np.mat(interp_theta)
        
        # 使用中心差分方法计算航向角        
        for i in range(1, self.interp_size_-1):
            dy = self.interp_y_[i + 1, 0] - self.interp_y_[i - 1, 0]
            dx = self.interp_x_[i + 1, 0] - self.interp_x_[i - 1, 0]    
            
            theta = np.arctan2(dy, dx)
            
            if abs(abs(theta) - math.pi) < 1e-4:
                theta = math.pi
            
            interp_theta[i, 0] = theta
        
        return np.mat(interp_theta)
    
    def calculate_kappa(self):
        # 初始化所有曲率为 0
        interp_kappa = np.zeros((self.interp_size_, 1))
        
        # 使用中心差分法计算曲率
        for i in range(1, self.interp_size_ - 1):
            dtheta = utils.normalize_angle(self.interp_theta_[i + 1, 0] - self.interp_theta_[i - 1, 0])
            ds = self.interp_s_[i + 1] - self.interp_s_[i - 1]
            interp_kappa[i, 0] = dtheta / ds
        
        # 补全两头的航向信息    
        interp_kappa[ 0, 0] = copy.deepcopy(interp_kappa[ 1, 0])
        interp_kappa[-1, 0] = copy.deepcopy(interp_kappa[-2, 0])

        return np.mat(interp_kappa)

    def calculate_dkappa(self):
        # 初始化曲率变化率为 0
        interp_dkappa = np.zeros((self.interp_size_, 1))
        
        # 使用中心差分法计算曲率
        for i in range(1, self.interp_size_ - 1):
            dkappa = self.interp_kappa_[i + 1, 0] - self.interp_kappa_[i - 1, 0]
            ds = self.interp_s_[i + 1] - self.interp_s_[i - 1]
            interp_dkappa[i, 0] = dkappa / ds

        # 补全两头的曲率变化率信息
        interp_dkappa[0, 0]  = copy.deepcopy(interp_dkappa[1, 0])
        interp_dkappa[-1, 0] = copy.deepcopy(interp_dkappa[-2, 0])

        return np.mat(interp_dkappa)


'''
计算五次多项式轨迹
'''
class QuinticPolynomial:
    def __init__(
        self,
        p_start     : float,   # 起点的函数值 p(t)
        p_dot_start : float,   # 起点函数值的导数 p'(t)
        p_ddot_start: float,   # 起点函数值的二阶导 p''(t)
        p_end       : float,   # 终点的函数值 p(t)
        p_dot_end   : float,   # 终点函数值的导数 p'(t) = 0
        p_ddot_end  : float,   # 终点函数值的二阶导 p''(t) = 0 
        sample_time : np.array # 离散的时间
        ):

        t_end = sample_time[-1]

        matrix_a = np.mat([[1, 0    , 0         , 0             , 0              , 0              ],
                           [0, 1    , 0         , 0             , 0              , 0              ],
                           [0, 0    , 2         , 0             , 0              , 0              ],
                           [1, t_end, t_end ** 2, t_end ** 3    , t_end ** 4     , t_end ** 5     ],
                           [0, 1    , 2 * t_end , 3 * t_end ** 2, 4 * t_end ** 3 , 5 * t_end ** 4 ],
                           [0, 0    , 2         , 6 * t_end     , 12 * t_end ** 2, 20 * t_end ** 3]])
        
        vector_b = np.mat([p_start, p_dot_start, p_ddot_start, p_end, p_dot_end, p_ddot_end]).T

        # calculate coefficients of quintic polynomial
        coeff_c = np.linalg.solve(matrix_a, vector_b)

        # d and temporal derivatives
        self.p_ = np.mat(coeff_c[0]
                         + coeff_c[1] * sample_time
                         + coeff_c[2] * sample_time ** 2
                         + coeff_c[3] * sample_time ** 3
                         + coeff_c[4] * sample_time ** 4
                         + coeff_c[5] * sample_time ** 5).T
        
        self.p_dot_ = np.mat(coeff_c[1]
                            + 2 * coeff_c[2] * sample_time
                            + 3 * coeff_c[3] * sample_time ** 2
                            + 4 * coeff_c[4] * sample_time ** 3
                            + 5 * coeff_c[5] * sample_time ** 4).T
        
        self.p_ddot_ = np.mat(2 * coeff_c[2]
                             + 6  * coeff_c[3] * sample_time
                             + 12 * coeff_c[4] * sample_time ** 2
                             + 20 * coeff_c[5] * sample_time ** 3).T
        
        self.p_dddot_ = np.mat(6 * coeff_c[3] 
                              + 24 * coeff_c[4] * sample_time 
                              + 60 * coeff_c[5] * sample_time ** 2).T

'''
计算四次多项式轨迹
'''        
class QuarticPolinomial:
    def __init__(
        self,
        p_start     : float, # 起点的函数值 p(t)
        p_dot_start : float, # 起点的函数导数 p'(t)
        p_ddot_start: float, # 起点得到函数二阶导 p''(t)
        p_dot_end   : float, # 终点的函数导数 p'(t) 
        p_ddot_end  : float, # 终点函数的二阶导 p''(t) = 0
        sample_time : np.array
        ):

        t_end = sample_time[-1]

        matrix_a = np.mat(
            [
                [1, 0, 0, 0, 0],
                [0, 1, 0, 0, 0],
                [0, 0, 2, 0, 0],
                [0, 1, 2 * t_end, 3 * t_end ** 2, 4 * t_end ** 3],
                [0, 0, 2, 6 * t_end, 12 * t_end ** 2],
            ]
        )
        vector_b = np.mat([p_start, p_dot_start, p_ddot_start, p_dot_end, p_ddot_end]).T

        # 计算四次多项式的参数
        coeff_c = np.linalg.solve(matrix_a, vector_b)

        # s and temporal derivatives
        self.p_ = np.mat(
            coeff_c[0]
            + coeff_c[1] * sample_time
            + coeff_c[2] * sample_time ** 2
            + coeff_c[3] * sample_time ** 3
            + coeff_c[4] * sample_time ** 4
        ).T
        self.p_dot_ = np.mat(
            coeff_c[1]
            + 2 * coeff_c[2] * sample_time
            + 3 * coeff_c[3] * sample_time ** 2
            + 4 * coeff_c[4] * sample_time ** 3
        ).T
        self.p_ddot_  = np.mat(2 * coeff_c[2] + 6 * coeff_c[3] * sample_time + 12 * coeff_c[4] * sample_time ** 2).T
        self.p_dddot_ = np.mat(6 * coeff_c[3] + 24 * coeff_c[4] * sample_time).T
        
'''
用于描述 Cartesian 坐标系下的局部轨迹
'''
class Trajectory:
    def __init__(
        self,
        longitudinal_curve: QuarticPolinomial,
        lateral_curve: QuinticPolynomial,
        t_array: np.array,
        reference_path: ReferencePath
        ):
        
        self.sample_time_ = t_array

        self.s_      = longitudinal_curve.p_
        self.s_dot_  = longitudinal_curve.p_dot_
        self.s_ddot_ = longitudinal_curve.p_ddot_
        self.s_dddot_ = longitudinal_curve.p_dddot_

        self.d_       = lateral_curve.p_
        self.d_dot_   = lateral_curve.p_dot_
        self.d_ddot_  = lateral_curve.p_ddot_
        self.d_dddot_ = lateral_curve.p_dddot_

        # 计算横 / 纵向轨迹的 Jerk 代价
        self.lateral_jerk_      = np.trapz(np.ravel(self.d_dddot_) ** 2, self.sample_time_)
        self.longitudinal_jerk_ = np.trapz(np.ravel(self.s_dddot_) ** 2, self.sample_time_)
        
        self.frenet_to_cartesian(reference_path)
    
    def frenet_to_cartesian(self, reference_path: ReferencePath):
        x_ref      = np.interp(self.s_, np.ravel(reference_path.interp_s_), np.ravel(reference_path.interp_x_))
        y_ref      = np.interp(self.s_, np.ravel(reference_path.interp_s_), np.ravel(reference_path.interp_y_))
        theta_ref  = np.interp(self.s_, np.ravel(reference_path.interp_s_), np.ravel(reference_path.interp_theta_))
        kappa_ref  = np.interp(self.s_, np.ravel(reference_path.interp_s_), np.ravel(reference_path.interp_kappa_))
        dkappa_ref = np.interp(self.s_, np.ravel(reference_path.interp_s_), np.ravel(reference_path.interp_dkappa_))
        
        self.x_ref_      = np.mat(x_ref)
        self.y_ref_      = np.mat(y_ref)
        self.theta_ref_  = np.mat(theta_ref)
        self.kappa_ref_  = np.mat(kappa_ref)
        self.dkappa_ref_ = np.mat(dkappa_ref)
        
        cartesian_x = np.ravel(self.x_ref_) - np.ravel(self.d_) * np.sin(np.ravel(self.theta_ref_))
        cartesian_y = np.ravel(self.y_ref_) + np.ravel(self.d_) * np.cos(np.ravel(self.theta_ref_))
        
        self.x_ = np.mat(cartesian_x).T
        self.y_ = np.mat(cartesian_y).T
    
        # 计算 d_prime
        d_prime = np.ravel(self.d_dot_) / np.ravel(self.s_dot_)
                
        # 计算 d(s) 的二阶导       
        d_pprime = (np.ravel(self.d_ddot_) - d_prime * np.ravel(self.s_ddot_)) / (np.ravel(self.s_dot_) ** 2)

        # 计算航向角增量
        delta_theta = np.arctan2(d_prime, (1 - np.ravel(kappa_ref) * np.ravel(self.d_)))
        
        # 计算航向角
        theta = np.ravel(theta_ref) + delta_theta
        self.theta_ = np.mat(theta).T
        
        # 计算车速
        v = np.ravel(self.s_dot_) * (1 - np.ravel(kappa_ref) * np.ravel(self.d_)) / np.cos(delta_theta)
        self.v_ = np.mat(v).T
        
        # 计算曲率
        kappa = (((d_pprime + (np.ravel(dkappa_ref) * np.ravel(self.d_) + np.ravel(kappa_ref) * d_prime) * np.tan(delta_theta))
            * (np.cos(delta_theta) ** 2) / (1 - np.ravel(kappa_ref) * np.ravel(self.d_))
            + np.ravel(kappa_ref)
        ) * np.cos(delta_theta) / (1 - np.ravel(kappa_ref) * np.ravel(self.d_)))    
        
        self.kappa_ = np.mat(kappa).T
        
        # 计算加速度
        a = np.ravel(self.d_ddot_) * (1 - np.ravel(kappa_ref) * np.ravel(self.d_)) / np.cos(delta_theta) + (
            np.ravel(self.s_dot_) ** 2
        ) / np.cos(delta_theta) * (
            (1 - np.ravel(kappa_ref) * np.ravel(self.d_))
            * np.tan(delta_theta)
            * (
                np.ravel(self.kappa_) * (1 - np.ravel(kappa_ref) * np.ravel(self.d_)) / np.cos(delta_theta)
                - np.ravel(kappa_ref)
            )
            - (np.ravel(dkappa_ref) * np.ravel(self.d_) + np.ravel(dkappa_ref) * d_prime)
        )
        
        self.a_ = np.mat(a).T
        return

'''
碰撞检测器
'''
class VehicleGeometry:
    # rl ------------- fl
    #    |           |
    #    |           |
    # rr ------------ fr
    def __init__(self, l: float, w: float):
        self.l_ = l
        self.w_ = w
        
        self.corner_fl_ = [ l/2,  w/2]
        self.corner_fr_ = [ l/2, -w/2]
        self.corner_rl_ = [-l/2,  w/2]
        self.corner_rr_ = [-l/2, -w/2]
        
        self.radius_ = math.sqrt((l/2)**2 + (w/2)**2)

    def check_collision_free(self, 
                             ego_x: float, 
                             ego_y: float, 
                             obstacle_state: struct.StaticObstacle) -> bool:
        
        obs_x = obstacle_state.pose_.x_
        obs_y = obstacle_state.pose_.y_
        
        dist  = utils.distance(ego_x, ego_y, obs_x, obs_y)
        
        # 如果车和障碍物之间的距离大于碰撞半径，则说明无碰撞
        if dist > self.radius_ + obstacle_state.radius_:
            return True
        else:
            return False
    
    def get_vehicle_shape(self, vehicle_pose: struct.Transform):
        ego_x     = vehicle_pose.x_
        ego_y     = vehicle_pose.y_
        ego_theta = vehicle_pose.theta_

        rotation_matrix = np.mat([[ math.cos(ego_theta),   math.sin(ego_theta)], 
                                  [-math.sin(ego_theta),  math.cos(ego_theta)]])

        original_bound = np.mat([self.corner_fl_, 
                                 self.corner_fr_, 
                                 self.corner_rr_, 
                                 self.corner_rl_, 
                                 self.corner_fl_])
        
        # Then, apply the rotation
        rotated_bound = original_bound * rotation_matrix

        rotated_bound[:, 0] += ego_x
        rotated_bound[:, 1] += ego_y

        return rotated_bound

        
'''
Frenet Optimal Planner
'''
class FrenetOptimalPlanner:
    def __init__(
        self, 
        reference_path: ReferencePath, 
        parameters: param_parser.Parameters
        ):
        self.min_terminal_t_ = parameters.min_terminal_t_ # 最小采样终止时间
        self.max_terminal_t_ = parameters.max_terminal_t_ # 最大采样终止时间
        self.num_terminal_t_ = parameters.num_terminal_t_ # 采样终止时间的个数
        self.num_terminal_t_ = int(self.num_terminal_t_)
        self.sample_step_    = parameters.sample_step_    # 采样步长
        
        self.min_terminal_v_ = parameters.min_terminal_v_ # 最小减速到 0
        self.max_terminal_v_ = parameters.max_terminal_v_ # 最快加速到 60 km/h
        self.num_terminal_v_ = parameters.num_terminal_t_ # 速度采样的个数
        self.num_terminal_v_ = int(self.num_terminal_v_)
        
        self.min_terminal_d_ = parameters.min_terminal_d_ # 左侧最大偏移量
        self.max_terminal_d_ = parameters.max_terminal_d_ # 右侧最大偏移量
        self.num_terminal_d_ = parameters.num_terminal_d_ # 横向采样个数 
        self.num_terminal_d_ = int(self.num_terminal_d_)
        
        self.lateral_jerk_coeff_    = parameters.lateral_jerk_coeff_    # 侧向加速度变化代价
        self.lateral_time_coeff_    = parameters.lateral_time_coeff_    # 换道时间代价
        self.lateral_deviate_coeff_ = parameters.lateral_deviate_coeff_ # 与参考线偏差代价
        
        self.longitudinal_jerk_coeff_    = parameters.longitudinal_jerk_coeff_    # 纵向加速度变化代价
        self.longitudinal_time_coeff_    = parameters.longitudinal_time_coeff_    # 加/减速时间代价
        self.longitudinal_deviate_coeff_ = parameters.longitudinal_deviate_coeff_ # 与目标车速偏差代价
        
        self.lateral_cost_coeff_      = parameters.lateral_cost_coeff_      # 横向代价
        self.longitudinal_cost_coeff_ = parameters.longitudinal_cost_coeff_ # 纵向代价
        
        self.target_v_ = self.max_terminal_v_ # 目标车速
        
        self.kappa_limit_ = parameters.kappa_limit_ # 曲率限制
        self.a_limit_     = parameters.a_limit_     # 加速度限制
        self.v_limit_     = parameters.v_limit_     # 速度限制
        
        # 碰撞检测距离阈值
        self.obs_check_thresh_ = self.v_limit_ * self.max_terminal_t_
        
        self.reference_path_ = reference_path
        
    def generate_trajectory_sample(self, frenet_state: struct.FrenetState) -> List[Trajectory]:
        longitudinal_sample_list  = []
        lateral_sample_list       = []
        cartesian_trajectory_list = []
        
        terminal_s_ddot = 0
        terminal_d_dot  = 0
        terminal_d_ddot = 0
        
        # 生成采样轨迹
        for terminal_t in np.linspace(self.min_terminal_t_, 
                                      self.max_terminal_t_, 
                                      self.num_terminal_t_):
            sample_time = np.arange(0, terminal_t + self.sample_step_, self.sample_step_)
            for terminal_s_dot in np.linspace(self.min_terminal_v_, 
                                              self.max_terminal_v_, 
                                              self.num_terminal_v_):
                # 如果始末都没速度，则跳过
                if frenet_state.s_dot_ == 0 and terminal_s_dot == 0:
                    continue
                
                longitudinal_sample = QuarticPolinomial(
                    frenet_state.s_, frenet_state.s_dot_, frenet_state.s_ddot_,
                    terminal_s_dot , terminal_s_ddot    , sample_time)
                
                longitudinal_sample_list.append(longitudinal_sample)
                
                for terminal_d in np.linspace(self.min_terminal_d_, 
                                              self.max_terminal_d_, 
                                              self.num_terminal_d_):    
                    lateral_sample = QuinticPolynomial(
                        frenet_state.d_, frenet_state.d_dot_, frenet_state.d_ddot_, 
                        terminal_d     , terminal_d_dot     , terminal_d_ddot, 
                        sample_time
                    )
                    
                    lateral_sample_list.append(lateral_sample)
                    
                    # 将横、纵向采样结果耦合，得到 Cartesian 坐标系下的轨迹
                    cartesian_trajectory = Trajectory(longitudinal_sample, 
                                                      lateral_sample, 
                                                      sample_time, 
                                                      self.reference_path_)
                    
                    cartesian_trajectory_list.append(cartesian_trajectory)
                            
        return cartesian_trajectory_list
                        
    def choose_optimal_trajectory(
        self, 
        cartesian_trajectory_list: list, 
        vehicle_geometry: VehicleGeometry, 
        obs_list: list
        ) -> Tuple[Trajectory, List[Trajectory]]:
        
        valid_trajectory = []
        trajectory_cost_list = []
        
        for cartesian_trajectory in cartesian_trajectory_list:
            # max_kappa = np.amax(cartesian_trajectory.kappa_)
            # max_a = np.amax(cartesian_trajectory.a_)
            # max_v = np.amax(cartesian_trajectory.v_) 
            min_v = np.amin(cartesian_trajectory.v_)
            # max_d_dot = np.amax(cartesian_trajectory.d_dot_)
            
            # 检查硬约束
            # if abs(max_kappa) > self.kappa_limit_ or max_a > self.a_limit_ or max_v > self.v_limit_ or min_v < 0:
            #     continue   
            if min_v < 0:
                continue
            
            # 检查碰撞约束
            collision_flag = 0
            for i in range(len(obs_list)):
                obs = obs_list[i]
                
                # 判断障碍物与自车的距离，如果距离比较远，则不进行碰撞检测
                obs_dist = utils.distance(obs.pose_.x_, obs.pose_.y_,
                                      cartesian_trajectory.x_[0, 0], 
                                      cartesian_trajectory.y_[0, 0])
                if obs_dist > self.obs_check_thresh_:
                    continue
                
                # 判断障碍物与自车在参考线上的距离差，如果距离比较短，则不进行碰撞检测
                _, obs_idx = self.reference_path_.kd_tree_.query([obs.pose_.x_, obs.pose_.y_])
                obs_s = self.reference_path_.interp_s_[obs_idx]
                if obs_s - cartesian_trajectory.s_[0] > self.obs_check_thresh_ or obs_s < cartesian_trajectory.s_[0]:
                    continue
                
                for j in range(len(cartesian_trajectory.s_)):
                    if not vehicle_geometry.check_collision_free(cartesian_trajectory.x_[j], 
                                                                 cartesian_trajectory.y_[j],
                                                                 obs):
                        collision_flag = 1
            
            if collision_flag:
                continue
                
            lateral_cost = (
                self.lateral_jerk_coeff_    * cartesian_trajectory.lateral_jerk_ + 
                self.lateral_time_coeff_    * cartesian_trajectory.sample_time_[-1] +
                self.lateral_deviate_coeff_ * cartesian_trajectory.d_[-1, 0] ** 2
                )
            
            longitudinal_cost = (
                self.longitudinal_jerk_coeff_ * cartesian_trajectory.longitudinal_jerk_ + 
                self.longitudinal_time_coeff_ * cartesian_trajectory.sample_time_[-1] + 
                self.longitudinal_deviate_coeff_  * (cartesian_trajectory.v_[-1, 0] - self.target_v_)**2
            )
            
            total_cost = (self.lateral_cost_coeff_ * lateral_cost 
                            + self.longitudinal_cost_coeff_ * longitudinal_cost)
            
            if math.isnan(total_cost):
                continue
            
            valid_trajectory.append(cartesian_trajectory)
            trajectory_cost_list.append(total_cost)
    
        if valid_trajectory:
            optimal_trajectory_cost = min(trajectory_cost_list)
            optimal_trajectory_index = trajectory_cost_list.index(optimal_trajectory_cost)
            optimal_trajectory = valid_trajectory[optimal_trajectory_index]
        else:
            optimal_trajectory = math.nan
        
        return optimal_trajectory, valid_trajectory
    
    def update_planning(
        self, 
        frenet_state: struct.FrenetState, 
        vehicle_geometry: VehicleGeometry, 
        obs_list: list
        ) -> Tuple[Trajectory, List[Trajectory], List[Trajectory]]:
        
        trajectory_list = self.generate_trajectory_sample(frenet_state)
        optimal_trajectory, valid_trajectory = self.choose_optimal_trajectory(trajectory_list, vehicle_geometry, obs_list)
        
        return optimal_trajectory, valid_trajectory, trajectory_list
        
'''
将笛卡尔坐标转换为 frenet 坐标
'''
def transform_cartesian_to_frenet(
    cartesian_state: struct.State, 
    reference_path: ReferencePath
    ) -> struct.FrenetState:

    dist, idx = reference_path.kd_tree_.query([cartesian_state.pose_.x_,
                                              cartesian_state.pose_.y_])
    
    # get x- and y-position along reference path
    x_reference  = reference_path.interp_x_[idx, 0]
    y_reference  = reference_path.interp_y_[idx, 0]
    xy_reference = np.array([x_reference, y_reference])

    # get s-coordinate along reference path
    s_reference = reference_path.interp_s_[idx]

    # get heading along reference path
    theta_reference = reference_path.interp_theta_[idx, 0]

    # get curvature along reference path
    kappa_reference = reference_path.interp_kappa_[idx, 0]

    # get first derivative of curvature along reference path
    dkappa_reference = reference_path.interp_dkappa_[idx, 0]

    # calculate difference in heading between trajectory and reference path
    delta_theta = cartesian_state.pose_.theta_ - theta_reference

    # current position of the car
    xy_state = np.array([cartesian_state.pose_.x_, cartesian_state.pose_.y_])

    # calculate distance vector
    d_vector = xy_state - xy_reference

    # calculate tangential vector
    tangential_vector = np.array([np.cos(theta_reference), np.sin(theta_reference)])

    # calculate sign of lateral position
    d_sign = np.sign(np.cross(a=tangential_vector, b=d_vector))

    # calculate quantities in frenèt coordinate system
    frenet_s = s_reference

    frenet_d = d_sign * dist
    
    frenet_s_dot = (cartesian_state.v_ * np.cos(delta_theta)) / (
        1 - kappa_reference * frenet_d
    )

    d_prime = (1 - kappa_reference * frenet_d) * np.tan(delta_theta)
    frenet_d_dot = frenet_s_dot * d_prime

    frenet_s_ddot = (
        np.cos(delta_theta)
        / (1 - kappa_reference * frenet_d)
        * (
            cartesian_state.a_
            - (frenet_s_dot ** 2)
            / np.cos(delta_theta)
            * (
                (1 - kappa_reference * frenet_d)
                * np.tan(delta_theta)
                * (
                    kappa_reference #cartesian_state.kappa_
                    * (1 - kappa_reference * frenet_d)
                    / np.cos(delta_theta)
                    - kappa_reference
                )
                - (dkappa_reference * frenet_d + kappa_reference * d_prime)
            )
        )
    )

    d_pprime = -(dkappa_reference * frenet_d + kappa_reference * d_prime) * np.tan(
        delta_theta
    ) + (1 - kappa_reference * frenet_d) / (np.cos(delta_theta) ** 2) * (
        kappa_reference #cartesian_state.kappa_
        * (1 - kappa_reference * frenet_d)
        / np.cos(delta_theta)
        - kappa_reference
    )
    frenet_d_ddot = (
        d_pprime * frenet_s_dot ** 2 + d_prime * frenet_s_ddot
    )

    frenet_s_dot = max(2, frenet_s_dot)
    
    frenet_state = struct.FrenetState(frenet_s, frenet_s_dot, frenet_s_ddot, 
                                      frenet_d, frenet_d_dot, frenet_d_ddot)
    
    return frenet_state