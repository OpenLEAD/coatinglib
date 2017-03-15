from numpy import array, random, reshape, c_

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

        if len(points)==0:
            return
        
        if key in self.handles:
            temp = self.handles[key]
            temp.append(self.env.plot3(points=points, pointsize=pointsize, colors=color))
            self.handles[key] = temp
        else:
            self.handles[key] = [self.env.plot3(points=points, pointsize=pointsize, colors=color)]
        return key         

    def plot_lists(self, lists_of_points, key='point_'+str(random.uniform(1,10)), color=(0,0,0), pointsize = 5):
        for list_of_points in lists_of_points:
            if len(list_of_points)>0:
                for list_of_point in list_of_points:
                    self.plot(list_of_points, key, color, pointsize)
        return key

    def plot_normal(self, points, key='normal_'+str(random.uniform(1,10)), color=(1,0,0), linewidth = 4):
        """
        Method to plot points, array(points), array(array(points))

        Keyword arguments:
        points -- points to be plotted.
        key -- string, key dictionary, name of the points.
        color -- tuple (R,G,B) 
        """
        points = array(points)
        if len(points.shape)==1:
            points = points.reshape(1,len(points))
        if len(points.shape)==3:
            points = [item for sublist in points for item in sublist]
            points = array(points)[:,0:6]
        
        if key in self.handles:
            temp = self.handles[key]
            temp.append(self.env.drawlinelist(points=reshape(c_[points[:,0:3], points[:,0:3]+0.05*points[:,3:6]],
                                                             (2*len(points),3)),
                                              linewidth=linewidth,
                                              colors=color))
            self.handles[key] = temp
        else:
            self.handles[key] = [self.env.drawlinelist(points=reshape(c_[points[:,0:3], points[:,0:3]+0.05*points[:,3:6]],
                                                             (2*len(points),3)),
                                              linewidth=linewidth,
                                              colors=color)]
        return key  

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

    def drawline(self, points, key='normal_'+str(random.uniform(1,10)), color=(1,0,0), linewidth = 4):
        if key in self.handles:
            temp = self.handles[key]
            temp.append(self.env.drawlinestrip(points=points, linewidth=3.0, colors=color))
            self.handles[key] = temp
        else:
            self.handles[key] = [self.env.drawlinestrip(points=points, linewidth=linewidth, colors=color)]
        return key
        
