from numpy import array, sin, cos, pi, maximum, minimum, random, tan, sqrt
from numpy import arctan2, sign, abs, mean, argmax, argmin, transpose
from mathtools import fast_poisson_disk
from time import time

"""
This module provides the RailPlace class and a generating function rand_rail 
"""

def _rand_angle(cfg, x, y, alpha_seed):
    # Given a set of (0,1) random alpha_seed,
    # the function output a uniform constrained by x/y/anglelimit set of angles. 
    
    alpha_min = cfg.environment.rail_angle_mean - cfg.environment.rail_angle_limit
    alpha_max = cfg.environment.rail_angle_mean + cfg.environment.rail_angle_limit
    
    # Get min and max angles for that position
    # (so to fit inside the avaible area)
    local_limit = [sign(y)*arctan2(abs(cfg.environment.x_max - x), abs(y)),
                   -sign(y)*arctan2(abs(cfg.environment.x_min - x), abs(y))]

    
    alpha_max = minimum( maximum(local_limit[0],local_limit[1]), alpha_max)
    alpha_min = maximum( minimum(local_limit[0],local_limit[1]), alpha_min)

    alpha = (alpha_max - alpha_min)*alpha_seed + alpha_min # - alphalimit <= alpha < alphalimit

    return alpha


def rand_rail(cfg, N = 1, equidistant = True):
    """
    This function generates N random objects of the class RailPlace,
    conditioned to the turbine's constrains.
    """
    
    #P - primary rail,
    #S - secondary rail,
    #alpha - angle from the perpendicular to the primary rail

    alpha_min = cfg.environment.rail_angle_mean - cfg.environment.rail_angle_limit
    alpha_max = cfg.environment.rail_angle_mean + cfg.environment.rail_angle_limit
    
    # Randomness
    random.seed(int(time()+int((time()-int(time()))*10000)))
    if equidistant:
        delta_x = cfg.environment.x_max - cfg.environment.x_min
        delta_y = cfg.environment.y_max - cfg.environment.y_min
        limits = array([[cfg.environment.x_min, cfg.environment.x_max],
                        [cfg.environment.y_min, cfg.environment.y_max]])
        x,y = transpose(fast_poisson_disk(sqrt(delta_x*delta_y*1.0/(N*sqrt(2))),limits))
    else:
        x,y,alpha = random.rand(3,N)

        x = (cfg.environment.x_max - cfg.environment.x_min)*x + cfg.environment.x_min


    if alpha_min > 0: #Limits of the avaible area depends on alpha limits
        y_max = (cfg.environment.x_max - x)/tan(alpha_min)
        y_min = (cfg.environment.x_min - x)/tan(alpha_min)
        y_max = minimum(y_max,cfg.environment.y_max)
        y_min = maximum(y_min,cfg.environment.y_min)
    elif alpha_max < 0:
        y_max = (cfg.environment.x_min - x)/tan(alpha_max)
        y_min = (cfg.environment.x_max - x)/tan(alpha_max)
        y_max = minimum(y_max,cfg.environment.y_max)
        y_min = maximum(y_min,cfg.environment.y_min)
    else:
        y_max = cfg.environment.y_max
        y_min = cfg.environment.y_min

    if equidistant:
        x = x[(y<y_max) & (y>y_min)]
        y = y[(y<y_max) & (y>y_min)]
        alpha = random.rand(len(x))
    else:
        y = (y_max - y_min)*y + y_min

    alpha = _rand_angle(cfg, x, y, alpha)

    P,S = xya2ps(x,y,alpha)
    
    return [RailPlace((p,s,a)) for p,s,a in zip(P,S,alpha)]

def xya2ps(x,y,alpha):
    """
    Compute primary (p) and secoundary (s) from (x,y) and angle (alpha).
    """
    S = y/cos(alpha)
    P = x + S*sin(alpha)
    return P,S

class RailPlace:
    """
    This class provides a simple way to store and manipulate
    information of the configuration space of the rails.
    That is:
    position of the (end of) primary rail counting from origin (p),
    angle between the secondary rail and a perpendicular from the primary rail (alpha)
    lenght of the secondary rail (s)
    """
    
    def __init__(self,(p,s,alpha) = (0,0,0)):
        self.p = p
        self.s = s
        self.alpha = alpha

    def __repr__(self):
        return ("p = " + str(self.p*100) + "cm" +
                "\ns = " + str(self.s*100) + "cm" +
                "\nalpha = " + str(self.alpha) + "rad ("+ str(round(self.alpha*180/pi,1))+"deg)")
                
    def setPSAlpha(self,(p,s,alpha)):
        self.p, self.s, self.alpha = p, s, alpha

    def getPSAlpha(self):
        # Get (p,s,alpha) as array
        return array((self.p, self.s, self.alpha))

    def getXYZ(self, cfg = None):
        # Get the XYZ value of the end of the secondary rail
        if cfg is None:
            return array([ self.p - self.s*sin(self.alpha), self.s*cos(self.alpha), 0 ])
        return array([ self.p - self.s*sin(self.alpha), self.s*cos(self.alpha), cfg.environment.z_floor_level ])

