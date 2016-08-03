from turbine import Turbine
from numpy import array, sin, cos, pi, maximum, minimum, random, tan, arctan2




def rand_rail(turbine, N = 1):
    #P - primary rail,
    #S - secondary rail,
    #alpha - angle from the perpendicular to the primary rail

    alpha_min = turbine.environment.rail_angle_mean - turbine.environment.rail_angle_limit
    alpha_max = turbine.environment.rail_angle_mean + turbine.environment.rail_angle_limit
    
    # Randomness
    x,y,alpha = random.rand(3,N)
    
    x = (turbine.environment.x_max - turbine.environment.x_min)*x + turbine.environment.x_min

    if alpha_min > 0: #Limits depend on alpha limits
        ymax = (turbine.environment.x_max - x)/tan(alpha_min)
        ymin = (turbine.environment.x_min - x)/tan(alpha_min)
        ymax = minimum(ymax,turbine.environment.y_max)
        ymin = maximum(ymin,turbine.environment.y_min)
    elif alpha_max < 0:
        ymax = (turbine.environment.x_min - x)/tan(alpha_max)
        ymin = (turbine.environment.x_max - x)/tan(alpha_max)
        ymax = minimum(ymax,turbine.environment.y_max)
        ymin = maximum(ymin,turbine.environment.y_min)
    else:
        ymax = turbine.environment.y_max
        ymin = turbine.environment.y_min
        
    y = (ymax - ymin)*y + ymin
    
    
    # Get min and max angles for that position (so to fit inside the avaible area)
    local_limit = [arctan2(turbine.environment.x_max - x, y),
                   arctan2(turbine.environment.x_min - x, y)]
    
    alpha_max = minimum( maximum(local_limit[0],local_limit[1]), alpha_max)
    alpha_min = maximum( minimum(local_limit[0],local_limit[1]), alpha_min)
    
    alpha = (alpha_max - alpha_min)*alpha + alpha_min # - alphalimit <= alpha < alphalimit
    
    # S/P conversion
    S = y/cos(alpha)
    P = x + S*sin(alpha)
    
    return [RailPlace((p,s,a)) for p,s,a in zip(P,S,alpha)]
    

class RailPlace:
    
    def __init__(self,(p,s,alpha) = (0,0,0)):
        self.p = p
        self.s = s
        self.alpha = alpha

    def __repr__(self):
        return ("p = " + str(self.p*100) + "cm" +
                "\ns = " + str(self.s*100) + "cm" +
                "\nalpha = " + str(self.alpha) + "rad ("+ str(round(self.alpha*180/pi,1))+"°)")
                
    def setPSAlpha(self,(p,s,alpha)):
        self.p, self.s, self.alpha = p, s, alpha

    def getPSAlpha(self):
        return array((self.p, self.s, self.alpha))

    def getXYZ(self, turbine):
        return array([ self.p - self.s*sin(self.alpha), self.s*cos(self.alpha), turbine.environment.z_floor_level ])

