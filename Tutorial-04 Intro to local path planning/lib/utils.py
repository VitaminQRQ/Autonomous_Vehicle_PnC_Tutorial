# Utils 中存放的是一些常用的计算工具

import math
import numpy as np
import matplotlib.pyplot as plt

'''
将角度归一化到 [-pi, pi]
'''
def normalize_angle(angle: float):
    normalized_angle = math.fmod(angle + math.pi, 2.0 * math.pi)
    
    if (normalized_angle < 0.0):
        normalized_angle = normalized_angle + 2.0 * math.pi
    
    return normalized_angle - math.pi

'''
将一整个列表的角度都归一化（一维）
'''
def normalize_angle_list(angle_list):
    normalized_angles = np.fmod(angle_list + np.pi, 2.0 * np.pi)
    normalized_angles[normalized_angles < 0.0] += 2.0 * np.pi
    
    return normalized_angles - np.pi

'''
计算两个坐标点 (x1, y1), (x2, y2) 之间的距离
'''
def distance(x1: float, y1: float, x2: float, y2: float):
    dx = x2 - x1
    dy = y2 - y1
    
    return math.sqrt(dy**2 + dx**2)

'''
用于辅助绘制圆形
'''
def plot_circle_helper(center_x: float, center_y: float, radius: float):
    
    theta = np.linspace(0, 2*math.pi, 100)
    x = radius * np.cos(theta) + center_x
    y = radius * np.sin(theta) + center_y
    
    return x, y