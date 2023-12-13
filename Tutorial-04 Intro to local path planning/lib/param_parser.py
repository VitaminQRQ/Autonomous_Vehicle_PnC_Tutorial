import configparser

class Parameters:
    def __init__(self, file_path):
        config = configparser.ConfigParser()
        config.read(file_path)

        # FRENET_OPTIMAL_PLANNER parameters
        self.min_terminal_t_ = config.getfloat('FRENET_OPTIMAL_PLANNER', 'min_terminal_t')
        self.max_terminal_t_ = config.getfloat('FRENET_OPTIMAL_PLANNER', 'max_terminal_t')
        self.num_terminal_t_ = config.getfloat('FRENET_OPTIMAL_PLANNER', 'num_terminal_t')
        self.sample_step_    = config.getfloat('FRENET_OPTIMAL_PLANNER', 'sample_step_')
        self.min_terminal_v_ = config.getfloat('FRENET_OPTIMAL_PLANNER', 'min_terminal_v')
        self.max_terminal_v_ = config.getfloat('FRENET_OPTIMAL_PLANNER', 'max_terminal_v')
        self.num_terminal_v_ = config.getfloat('FRENET_OPTIMAL_PLANNER', 'num_terminal_v')
        self.min_terminal_d_ = config.getfloat('FRENET_OPTIMAL_PLANNER', 'min_terminal_d')
        self.max_terminal_d_ = config.getfloat('FRENET_OPTIMAL_PLANNER', 'max_terminal_d')
        self.num_terminal_d_ = config.getfloat('FRENET_OPTIMAL_PLANNER', 'num_terminal_d')
        
        self.v_limit_        = config.getfloat('FRENET_OPTIMAL_PLANNER', 'v_limit')
        self.a_limit_        = config.getfloat('FRENET_OPTIMAL_PLANNER', 'a_limit')
        self.kappa_limit_    = config.getfloat('FRENET_OPTIMAL_PLANNER', 'kappa_limit')
        
        self.lateral_jerk_coeff_         = config.getfloat('FRENET_OPTIMAL_PLANNER', 'lateral_jerk_coeff')
        self.lateral_time_coeff_         = config.getfloat('FRENET_OPTIMAL_PLANNER', 'lateral_time_coeff')
        self.lateral_deviate_coeff_      = config.getfloat('FRENET_OPTIMAL_PLANNER', 'lateral_deviate_coeff')
        self.longitudinal_jerk_coeff_    = config.getfloat('FRENET_OPTIMAL_PLANNER', 'longitudinal_jerk_coeff')
        self.longitudinal_time_coeff_    = config.getfloat('FRENET_OPTIMAL_PLANNER', 'longitudinal_time_coeff')
        self.longitudinal_deviate_coeff_ = config.getfloat('FRENET_OPTIMAL_PLANNER', 'longitudinal_deviate_coeff')
        self.lateral_cost_coeff_         = config.getfloat('FRENET_OPTIMAL_PLANNER', 'lateral_cost_coeff')
        self.longitudinal_cost_coeff_    = config.getfloat('FRENET_OPTIMAL_PLANNER', 'longitudinal_cost_coeff')

        # LONGITUDINAL_PID parameters
        self.Kp_ = config.getfloat('LONGITUDINAL_PID', 'kp')
        self.Ki_ = config.getfloat('LONGITUDINAL_PID', 'ki')
        self.Kd_ = config.getfloat('LONGITUDINAL_PID', 'kd')

        # LATERAL_STANLEY parameters
        self.K_ = config.getfloat('LATERAL_STANLEY', 'k')
