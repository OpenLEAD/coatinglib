from numpy import array, random

class Visualizer:
    """
    QTcoin initializer for Openrave simulator visualization.
    """

    def __init__(self, environment):
        self.env = environment
        self.env.SetViewer('qtcoin')
        self.handles = {}

    def plot(self, points, key='point_'+str(random.uniform(1,10)), color=(0,0,0), pointsize = 5):
        """
        Method to plot points, array(points), array(array(points))

        Keyword arguments:
        points -- points to be plotted.
        key -- string, key dictionary, name of the points.
        color -- tuple (R,G,B) 
        """
        points = array(points)
        if len(points.shape)==1:
            points = points[0:3]
        if len(points.shape)==2:
            points = points[:,0:3]
        if len(points.shape)==3:
            points = [item for sublist in points for item in sublist]
            points = array(points)[:,0:3]
        
        if key in self.handles:
            temp = self.handles[key]
            temp.append(self.env.plot3(points=points, pointsize=pointsize, colors=color))
            self.handles[key] = temp
        else:
            self.handles[key] = [self.env.plot3(points=points, pointsize=pointsize, colors=color)]
        return            

    def remove_points(self, key):
        """
        Method to remove points

        Keyword arguments:
        key -- string, key dictionary, name of the points.
        """
        try:
            self.handles.pop(key,None)
            self.update()
        except AttributeError:
            raise AttributeError('There is not that key')

    def update(self):
        """
        Update the environment.
        """
            
        self.env.UpdatePublishedBodies()
        
