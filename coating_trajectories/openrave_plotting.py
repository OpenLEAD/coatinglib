from numpy import array

def plot_points(turbine, points, key, color):
    try:
        turbine.handles
    except AttributeError:
        turbine.handles={}
    
    if key in turbine.handles:
        temp = turbine.handles[key]
        temp.append(turbine.env.plot3(points=array(points)[:,0:3],pointsize=5,colors=color))
        turbine.handles[key] = temp
    else:    
        turbine.handles[key] = [turbine.env.plot3(points=array(points)[:,0:3],pointsize=5,colors=color)]
    
def plot_point(turbine, point, key, color):
    try:
        turbine.handles
    except AttributeError:
        turbine.handles={}
        
    if key in turbine.handles:
        temp = turbine.handles[key]
        temp.append(turbine.env.plot3(points=array(point)[0:3],pointsize=5,colors=color))
        turbine.handles[key] = temp
    else:
        turbine.handles[key] = [turbine.env.plot3(points=array(point)[0:3],pointsize=5,colors=color)]

def plot_points_array(turbine, pointsArray, key, color):
    try:
        turbine.handles
    except AttributeError:
        turbine.handles={}
        
    points = [item for sublist in pointsArray for item in sublist]
    plotPoints(turbine, points, key, color)

def remove_points(turbine, key):
    try:
        turbine.handles.pop(key,None)
        turbine.env.UpdatePublishedBodies()
    except AttributeError:
        turbine.handles={}
        
        
        
