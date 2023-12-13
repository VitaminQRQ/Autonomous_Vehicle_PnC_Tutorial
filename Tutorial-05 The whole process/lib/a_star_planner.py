import math
import cv2
import numpy as np
import matplotlib.pyplot as plt

class Node:
    def __init__(
        self, 
        position: tuple, 
        g: float, 
        f: float, 
        parent: 'Node'
        ) -> None:
        
        self.position = position
        self.g = g
        self.f = f
        self.parent_node = parent

class AStar:
    def __init__(
        self, 
        start_x: int, 
        start_y: int, 
        end_x  : int, 
        end_y  : int,
        map_image: np.ndarray,
        step_size    = 5,
        search_range = 40, 
        scale = 5
        ) -> None:
        
        self.start_x_   = start_x // scale
        self.start_y_   = start_y // scale
        self.end_x_     = end_x   // scale
        self.end_y_     = end_y   // scale
        self.scale_     = scale
        
        resize_height   = map_image.shape[0] // scale
        resize_width    = map_image.shape[1] // scale
        self.map_image_ = cv2.resize(map_image, (resize_width, resize_height))
        
        self.step_size_    = step_size
        self.search_range_ = search_range
        
    def search_path(self):
        start_node = Node((self.start_x_, self.start_y_), 0       , 0       , None)
        end_node   = Node((self.end_x_  , self.end_y_ ) , math.inf, math.inf, None)

        open_list = []
        close_list = [start_node]
        
        step_size = self.step_size_
        map_image = self.map_image_
        
        map_height = map_image.shape[0]
        map_width  = map_image.shape[1]
        
        while 1:
            current_node = close_list[-1]
            
            possible_moves = (
                [ 0        ,  step_size], 
                [ 0        , -step_size], 
                [ step_size,  0        ], 
                [-step_size,  0        ]
                )
            
            for move in possible_moves:
                g_cost = current_node.g + 1
                x = current_node.position[0] + move[0]
                if x < 0 or x >= map_width:
                    continue
                y = current_node.position[1] + move[1]
                if y < 0 or y >= map_height:
                    continue
                if map_image[y][x] == 0:
                    continue
                h_cost = abs(x - end_node.position[0]) + abs(y - end_node.position[1])
                f_cost = g_cost + h_cost

                if h_cost < self.search_range_:
                    step_size = 1
                new_node = Node((x, y), g_cost, f_cost, current_node)
                count = 0
                for node in open_list:
                    if node.position == new_node.position:
                        count += 1
                        if node.f > new_node.f:
                            node.g = new_node.g
                            node.f = new_node.f
                            node.parent_node = current_node
                for node in close_list:
                    if node.position == new_node.position:
                        count += 1
                if count == 0:
                    open_list.append(new_node)

            temp_node = None
                
            for i in range(len(open_list)):
                if temp_node is None or open_list[i].f < temp_node.f:
                    temp_node = open_list[i]
            
            if temp_node is None:
                print('No path found to the destination.')
                return []  
            
            for i in range(len(open_list)):
                if temp_node == open_list[i]:
                    open_list.pop(i)
                    break
            close_list.append(temp_node)
            
            if temp_node.position == end_node.position:
                print('End point found')
                break
        return close_list

    def retrieve_path(self, close_list):
        start_node_position = (self.start_x_, self.start_y_)
        path = []
        path.append(close_list[-1])
        point = path[-1]

        while 1:
            for node in close_list:
                if node.position == (point.parent_node.position if point.parent_node else None):
                    point = node
                    path.append(point)
            if point.position == start_node_position:
                print('Path search completed')
                break       
        return path 
    
    def update_planning(self):
        close_list = self.search_path()
        
        if close_list == []:
            return []
        
        path = self.retrieve_path(close_list)

        resize_path = []
        for node in path:
            resize_path.append((node.position[0] * self.scale_,
                                node.position[1] * self.scale_))
        
        return resize_path