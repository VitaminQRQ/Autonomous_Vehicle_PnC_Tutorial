import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

import lib.frenet_optimal_planner as fop
import lib.utils as utils
import lib.data_struct as struct
import lib.vehicle_simulator as sim
import lib.controller as controller
import lib.param_parser as param_parser

def plot_trajectory(i, 
                    optimal_trajectory_list, 
                    valid_trajectory_list, 
                    actual_x_list, 
                    actual_y_list, 
                    actual_theta_list,
                    obs_index_list, 
                    obs_x, 
                    obs_y, 
                    obs_r, 
                    reference_x, 
                    reference_y, 
                    vehicle_geometry):
    
    valid_trajectory   = valid_trajectory_list[i]
    optimal_trajectory = optimal_trajectory_list[i]

    ego_x = actual_x_list[i]
    ego_y = actual_y_list[i]
    ego_theta = actual_theta_list[i]
    
    plt.clf()

    if valid_trajectory == []:
        return

    # 绘制车辆当前位置
    ego_bound_x, ego_bound_y = utils.plot_circle_helper(ego_x, ego_y, vehicle_geometry.radius_)
    plt.plot(ego_bound_x, ego_bound_y, color='tomato', linewidth=0.7)

    vehicle_pose = struct.Transform(ego_x, ego_y, ego_theta)
    ego_shape = vehicle_geometry.get_vehicle_shape(vehicle_pose)
    plt.plot(ego_shape[:, 0], ego_shape[:, 1], color='orange', linewidth=1)
    
    # 绘制障碍物
    for j in range(len(obs_index_list)):   # 使用不同的迭代变量
        obs_bound_x, obs_bound_y = utils.plot_circle_helper(obs_x[j], obs_y[j], obs_r[j])
        plt.plot(obs_bound_x, obs_bound_y, color='black', label='obstacle')    

    # 绘制参考路径
    plt.plot(reference_x, reference_y, label='reference_trajectory')

    # 绘制采样路径
    for trajectory in valid_trajectory:
        plt.plot(trajectory.x_, trajectory.y_, color='gray', linewidth=0.5)

    # 绘制最优路径
    plt.plot(optimal_trajectory.x_, optimal_trajectory.y_, color='red', linewidth=0.5)

    # 绘制已行驶路径
    plt.plot(actual_x_list[0 : i+1], actual_y_list[0 : i+1], color='orange')

    # 绘制
    plt.axis('equal')
    plt.grid(True)

    area = 30
    plt.xlim(ego_x - area, ego_x + area)
    plt.ylim(ego_y - area, ego_y + area)

if __name__ == '__main__':
    # 初始化算法参数
    parameters = param_parser.Parameters('./config/params.ini')
    
    # 通过 numpy 的 genfromtxt 读取 csv 文件，并跳过第一行
    trajectory_data = np.genfromtxt('./data/global_trajectory.csv', delimiter=';', skip_header=1)
    reference_trajectory = np.mat(trajectory_data)

    # 将数组分割，方便后续操作
    reference_s = reference_trajectory[:, 0] # [m]
    reference_x = reference_trajectory[:, 1] # [m]
    reference_y = reference_trajectory[:, 2] # [m]
    
    # 建立 ReferencePath
    ref_path = fop.ReferencePath(reference_x, reference_y, 0.1)    
                
    # 初始化车辆位置和碰撞检测器
    car_length = 3
    car_width  = 2
    
    initial_pose = struct.Transform(ref_path.interp_x_[0, 0], 
                                    ref_path.interp_y_[0, 0], 
                                    ref_path.interp_theta_[0, 0])

    initial_state = struct.State(initial_pose, v=5, a=0, kappa=0)
    vehicle_geometry = fop.VehicleGeometry(car_length, car_width)

    # 生成测试场景
    obs_index_gap  = 30
    obs_index_list = np.floor(np.arange(obs_index_gap, len(reference_s), obs_index_gap))
    obs_x = np.zeros_like(obs_index_list)
    obs_y = np.zeros_like(obs_index_list)
    obs_r = np.zeros_like(obs_index_list)

    obs_list = []

    for i in range(len(obs_index_list)):
        index = int(obs_index_list[i])
        obs_x[i] = reference_x[index] + np.random.uniform(-5, 5, 1)
        obs_y[i] = reference_y[index] + np.random.uniform(-5, 5, 1)
        obs_r[i] = np.random.uniform(0, 2, 1)
        
        obs_pose  = struct.Transform(obs_x[i], obs_y[i], 0)
        obs_state = struct.StaticObstacle(obs_pose, obs_r[i])

        obs_list.append(obs_state)
            
    # 设置车辆初始状态
    initial_x     = ref_path.interp_x_[0, 0]
    initial_y     = ref_path.interp_y_[0, 0]
    initial_theta = ref_path.interp_theta_[0, 0]

    initial_kappa = 0
    initial_v     = 0
    initial_a     = 0

    vehicle_pose  = struct.Transform(initial_x, initial_y, initial_theta)
    vehicle_state = struct.State(vehicle_pose, initial_v, initial_a, initial_kappa)

    # 创建控制器
    stanley = controller.Stanley(parameters.K_)
    pid = controller.PID(parameters.Kp_, 
                         parameters.Ki_, 
                         parameters.Kd_)
    
    # 创建车辆仿真器
    wheel_base  = 2
    sim_vehilce = sim.Vehicle(wheel_base, vehicle_state)

    # 设置仿真条件
    total_sim_time = 300
    sample_time    = 0.1
    sim_time       = np.arange(0, total_sim_time, sample_time)

    # 初始化记录轨迹的空数组
    actual_x_list = []
    actual_y_list = []
    actual_v_list = []
    actual_theta_list = []
    
    sample_trajectory_list  = []
    valid_trajectory_list   = []
    optimal_trajectory_list = []
    
    # 规划一个轨迹
    frenet_state = fop.transform_cartesian_to_frenet(vehicle_state, ref_path)
    planner      = fop.FrenetOptimalPlanner(ref_path, parameters)

    optimal_trajectory, valid_trajectory, trajectory_list = planner.update_planning(frenet_state, vehicle_geometry, obs_list)

    print("Processing ...")
    for i in range(len(sim_time)):
        if valid_trajectory != [] and len(optimal_trajectory.s_) > 2:
            index = 1
        else:
            print("\nNo valid path")
            break
        
        if optimal_trajectory.s_[0] >= reference_s[-1, 0]:
            actual_x_list.append(vehicle_state.pose_.x_)
            actual_y_list.append(vehicle_state.pose_.y_)
            actual_v_list.append(vehicle_state.v_)
            print("\nFinish")
            break
        
        actual_x_list.append(vehicle_state.pose_.x_)
        actual_y_list.append(vehicle_state.pose_.y_)
        actual_theta_list.append(vehicle_state.pose_.theta_)
        actual_v_list.append(vehicle_state.v_)
        
        valid_trajectory_list.append(valid_trajectory)
        optimal_trajectory_list.append(optimal_trajectory)
        sample_trajectory_list.append(trajectory_list)
        
        # 寻找控制器参考点
        front_x = vehicle_state.pose_.x_ + wheel_base * math.cos(vehicle_state.pose_.theta_)
        front_y = vehicle_state.pose_.y_ + wheel_base * math.sin(vehicle_state.pose_.theta_)
        
        dx = np.ravel(optimal_trajectory.x_) - front_x
        dy = np.ravel(optimal_trajectory.y_) - front_y
        dist = np.abs(dx + dy)
        min_dist_idx = np.argmin(dist)
        min_dist = np.min(dist)
        
        target_x     = optimal_trajectory.x_[min_dist_idx, 0]
        target_y     = optimal_trajectory.y_[min_dist_idx, 0]
        target_theta = optimal_trajectory.theta_[min_dist_idx, 0]
        target_v     = optimal_trajectory.v_[min_dist_idx, 0]
        target_pose  = struct.Transform(target_x, target_y, target_theta)
        
        # 计算控制量
        delta         = stanley.update_control(target_pose, vehicle_state)
        accel         = pid.update_control(target_v, vehicle_state.v_, sample_time)
        vehicle_state = sim_vehilce.update_state(delta, accel, sample_time)
        
        # 将笛卡尔坐标系转换到 Cartesian 坐标系
        dx = np.ravel(optimal_trajectory.x_) - vehicle_state.pose_.x_
        dy = np.ravel(optimal_trajectory.y_) - vehicle_state.pose_.y_
        dist = np.abs(dx**2 + dy**2)
        min_dist_idx = np.argmin(dist)
        min_dist = np.min(dist)
           
        # 从上一帧轨迹中寻找当前位置的投影
        # 并以此为起点规划下一帧轨迹 
        frenet_state = struct.FrenetState(optimal_trajectory.s_[min_dist_idx, 0], 
                                            optimal_trajectory.s_dot_[min_dist_idx, 0], 
                                            optimal_trajectory.s_ddot_[min_dist_idx, 0], 
                                            optimal_trajectory.d_[min_dist_idx, 0], 
                                            optimal_trajectory.d_dot_[min_dist_idx, 0], 
                                            optimal_trajectory.d_ddot_[min_dist_idx, 0])
        
        optimal_trajectory, valid_trajectory, trajectory_list = planner.update_planning(frenet_state, 
                                                                                        vehicle_geometry, 
                                                                                        obs_list)
        
        progress = (i + 1) / len(sim_time)
        bar = "#" * int(progress * 20)
        print("\rApproximate Progress: [{:<20}] {:.1f}%".format(bar, progress * 100), end="")
    print()
    
    print("Generating video ...")        
    fig = plt.figure()
    ani = animation.FuncAnimation(fig, 
                                  plot_trajectory, 
                                  frames=len(optimal_trajectory_list), 
                                  fargs=(
                                      optimal_trajectory_list, 
                                      valid_trajectory_list, 
                                      actual_x_list, 
                                      actual_y_list, 
                                      actual_theta_list,
                                      obs_index_list, 
                                      obs_x, 
                                      obs_y, 
                                      obs_r, 
                                      reference_x, 
                                      reference_y, 
                                      vehicle_geometry
                                      )
                                  )
    
    ani.save('./media/trajectory.gif', writer='pillow', fps=30)     
    print("Done!")