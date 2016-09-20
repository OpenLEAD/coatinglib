from numpy import load, savez_compressed, zeros, ones, arange, r_, c_, outer, tile
from numpy import meshgrid, array, shape, sum, eye, dot, argsort, concatenate, sqrt
from numpy import argmax, argmin
from os import makedirs
import errno
from openravepy import RaveCreateCollisionChecker, matrixFromAxisAngle
from scipy.spatial import KDTree
import mathtools
from math import pi, ceil
from copy import copy, deepcopy
from lxml import etree as ET

class IterSurfaceError(Exception):    
    def __init__(self):
        Exception.__init__(self, "object is not a valid surface.")

class BladeModeling:
    """ BladeModeling class for blade modelling.

    Keyword arguments:
    name -- the name of the blade.
    model_type -- model object (e.g. RBF).
    env -- environment object
    blade -- body blade object.
    """

    turbine = None
    name = None
    blade = None
    _points = []
    
    def __init__(self, name, turbine, blade):
        self.turbine = turbine
        self._name = name
        self._blade = blade    

    def save_samples(self, directory_to_save):
        
        try:
            makedirs(directory_to_save)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise   
        savez_compressed(directory_to_save + self._name + '_points.npz', array=self._points)
        return

    def save_model(self, models, directory_to_save, iter_surface=None, model_index=[]):

        """ A method to save model and model info

        Keyword arguments:
        models -- the computed models.
        model_index -- when to switch between models .
        iter_surface -- the iter surface used to compute the trajectories.
        directory_to_save -- where to save the files.
        """

        try:
            makedirs(directory_to_save)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise

        model = ET.Element("model")
        ET.SubElement(model, "name").text = self._name

        if iter_surface is not None:
            surface = ET.SubElement(model, "iter_surface")
            ET.SubElement(surface, "type").text = iter_surface.name()
            ET.SubElement(surface, "Rn").text = str(iter_surface._Rn)

        ET.SubElement(model, "index").text = str(model_index)

        for i in range(0,len(models)):
            doc = ET.SubElement(model, "interpolation")
            ET.SubElement(doc, "type").text = models[i].model_type
            ET.SubElement(doc, "w").text = directory_to_save + 'w_' + str(i)
            ET.SubElement(doc, "points").text = directory_to_save + 'points_' + str(i)
            if model_index:
                ET.SubElement(doc, "index").text = str(model_index[i])
            savez_compressed(directory_to_save + 'w_' + str(i) + '.npz', array=models[i]._w)
            savez_compressed(directory_to_save + 'points_' + str(i) + '.npz', array=models[i]._points)
            if models[i].model_type == 'RBF':
                ET.SubElement(doc, "kernel").text = models[i]._kernel
                if models[i]._kernel=='gaussr':
                    parameters = ET.SubElement(model_type, "paremeters")
                    ET.SubElement(paremeters, "gauss_parameter").text = str(models[i].gausse)
                ET.SubElement(doc, "eps").text = str(models[i]._eps)
  
        tree = ET.ElementTree(model)
        tree.write(directory_to_save + "model.xml", pretty_print=True)
        return

    def save_trajectory(self, trajectories, xml_model, iter_surface,
                        tangent_step, directory_to_save):
        """ A method to save trajectory and trajectory info

        Keyword arguments:
        trajectories -- the computed trajectories.
        xml_model -- path to xml model file.
        iter_surface -- the iter surface used to compute the trajectories.
        tangent_step -- the tangent step used.
        directory_to_save -- where to save the files.
        """

        try:
            makedirs(directory_to_save)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise

        trajectory = ET.Element("trajectory")
        ET.SubElement(trajectory, "name").text = self._name
        ET.SubElement(trajectory, "tangent_step").text = str(tangent_step)
        
        ET.SubElement(trajectory, "xml_model").text = xml_model
        
        surface = ET.SubElement(trajectory, "iter_surface")
        ET.SubElement(surface, "type").text = iter_surface.name()
        ET.SubElement(surface, "Rn").text = str(iter_surface._Rn)
        ET.SubElement(surface, "stopR").text = str(iter_surface.stopR)
        ET.SubElement(surface, "coatingstep").text = str(iter_surface.coatingstep)

        savez_compressed(directory_to_save + 'trajectory.npz', array=trajectories)
        
        tree = ET.ElementTree(trajectory)
        tree.write(directory_to_save + "trajectory.xml", pretty_print=True)
        return

    def load_samples(self, directory_to_load = 'Blade/'):
        try:
            self._points = load(directory_to_load+self._name+'_points.npz')
            self._points = self._points['array']
        except IOError:
            raise IOError('Samples could not be loaded. Call method BladeModeling::sampling')
        return

    def load_model(self, xml_model):

        root = ET.parse(open(xml_model)).getroot()
        try:
            model._w = load(directory_to_load+model.model_type+'/'+model._name+'_w.npz')
            model._w = model._w['array']
            model._points = load(directory_to_load+model.model_type+'/'+model._name+'_points.npz')
            model._points = model._points['array']
        except IOError:
            raise IOError("Model could not be loaded.")
        return
    
    def load_trajectories(self, model, directory_to_load = 'Blade/Trajectory/'):
        try:
            self._trajectories = load(directory_to_load+model._name+'_trajectories.npz')
            self._trajectories = self._trajectories['array']
            self._trajectories = self._trajectories.tolist()
        except IOError:
            raise IOError("Trajectories could not be loaded.")
        return

    def sampling(self, delta = 0.005, min_distance_between_points=0.05):
        """ The sampling method is an algorithm for uniformly sampling objects.
        It computes the bounding box of the object (object inscribed by cube), uniformly
        samples the box's faces, and checks rays collision from cube's samples to the object.
        Rays with collision and its normal vectors are stored as object's samples.
        Finally, a kdtree is computed from object's samples and neighbooring points are removed  
        to generate an uniformly sampling on the object (filtering data by distance).

        This is a data-intensive computing and might freeze your computer.

        Keyword arguments:
        delta -- uniform distance between the cube's samples.
        min_distance_between_points -- uniform distance between the object's samples.
        min_distance_between_points >= delta
        directory_to_save -- directory to save samples (default is './Blade')
        """
        
        cc = RaveCreateCollisionChecker(self.turbine.env,'ode')
        if cc is not None:
                ccold = self.turbine.env.GetCollisionChecker()
                self.turbine.env.SetCollisionChecker(cc)
                cc = ccold
 
        self._blade.SetTransform(eye(4))
        self._blade.SetTransform(dot(self._blade.GetTransform(),
                                     matrixFromAxisAngle([0, -self.turbine.config.environment.blade_angle, 0]))
                                 )
        
        ab = self._blade.ComputeAABB()
        p = ab.pos()
        e = ab.extents()+0.01 # increase since origin of ray should be outside of object
        sides = array((
                     (e[0],0,0,-1,0,0,0,e[1],0,0,0,e[2]), #x
                     (-e[0],0,0,1,0,0,0,e[1],0,0,0,e[2]), #-x
                     (0,0,e[2],0,0,-1,e[0],0,0,0,e[1],0), #z
                     (0,0,-e[2],0,0,1,e[0],0,0,0,e[1],0), #-z
                     (0,e[1],0,0,-1,0,e[0],0,0,0,0,e[2]), #y
                     (0,-e[1],0,0,1,0,e[0],0,0,0,0,e[2])  #-y
                     ))
        maxlen = 2*sqrt(sum(e**2))+0.03
        self._points = zeros((0,6))
        for side in sides:
            ex = sqrt(sum(side[6:9]**2))
            ey = sqrt(sum(side[9:12]**2))
            XX,YY = meshgrid(r_[arange(-ex,-0.25*delta,delta),0,arange(delta,ex,delta)],
                                  r_[arange(-ey,-0.25*delta,delta),0,arange(delta,ey,delta)])
            localpos = outer(XX.flatten(),side[6:9]/ex)+outer(YY.flatten(),side[9:12]/ey)
            N = localpos.shape[0]
            rays = c_[tile(p+side[0:3],(N,1))+localpos,maxlen*tile(side[3:6],(N,1))]
            collision, info = self.turbine.env.CheckCollisionRays(rays, self._blade)
            # make sure all normals are the correct sign: pointing outward from the object)
            newinfo = info[collision,:]
            if len(newinfo) > 0:
                  newinfo[sum(rays[collision,3:6]*newinfo[:,3:6],1)>0,3:6] *= -1
                  self._points = r_[self._points,newinfo]

        self._points = self.filter_by_distance(self._points, min_distance_between_points)
        return             
   
    def filter_by_distance(self, points, r):
        """
        The filter_by_distance method is an algorithm to delete the nearest neighbors
        points, inside a distance threshold. 

        Keyword arguments:
        points -- points to be filtered array(array).
        r -- distance threshold. If r is None, points are sorted w.r.t. x axis.
        """

        points = points[argsort(points[:,0])]
        if r is None:
            return points
        
        Tree = KDTree(points[:,0:3])
        rays = []
        N = len(points)
        i=0
        I = ones((len(points), 1), dtype=bool)
        while True:
            if I[i]:
                rays.append(points[i])
                idx = Tree.query_ball_point(points[i,0:3],r)
                idx = array(idx)
                idx = idx[idx>i]
                for j in idx:I[j]=False
            i+=1
            if i==len(points):break
        return array(rays)
        
    def make_model(self, model, iter_surface=None, number_of_points_per_part = 5000,
                   intersection_between_divisions = 0.15):
        """
        The make_model method is an algorithm to generate a mathematical representation
        of the object (mesh). This method can be called after the sampling method,
        and never before. The model is generated by the model_type (constructor), e.g. RBF.

        This is a data-intensive computing and might freeze your computer.

        Keyword arguments:
        model -- model type for the implicit computation, e.g. the RBF model (rbf.py).
        iter_surface -- This is the surface to be iterated, as mathtools.sphere.
        If the model has more than 9000 points, this argument is needed.
        directory_to_save -- directory to save the model (default is 'Blade/')
        """

        model_points = [self._points]
        model_index = []
        
        if len(self._points)>9000:
            if iter_surface is None:
                raise ValueError("The number of points is "+str(len(self._points))+""", which is
                                bigger than 9000. It is not safe to make an unique model this big.
                                Create an iterate surface to partitionate the model
                                (check mathtools.IterSurface).""")
            elif not issubclass(iter_surface.__class__, mathtools.IterSurface):
                    raise IterSurfaceError()
            else:
                model_points, model_index = self.divide_model_points(self._points, iter_surface,
                                                                     number_of_points_per_part,
                                                                     intersection_between_divisions)                

        models = []
        for i in range(0,len(model_points)):
            model._points = model_points[i]
            model.make()
            models.append(model)
        return models, model_index

    def divide_model_points(self, points, iter_surface,
                            number_of_points_per_part = 5000,
                            intersection_between_divisions = 0.15):
        """
        This method divides the points for multiple model generation, e.g. multiple RBFs.
        It is required for objects with heavy sampling density, as the computation
        of large RBFs is not possible due to the computer memory capacity.

        Keyword arguments:
        points -- samples to be divided. Each division will produce a model.
        iter_surface -- This is the surface to be iterated, as mathtools.sphere.
        number_of_points_per_part -- number of samples per division (5000 default)
        intersection_between_divisions -- coeficient that multiplies the number_of_points_per_part.
        The result is the samples which belong to two RBFs simultaneously. 
        """

        if number_of_points_per_part >= len(points):
            return [points], []

        model_points = []
        points_distance = iter_surface.f_array(points)
        points = [x for (y,x) in sorted(zip(points_distance, points))]
        points = array(points)
        number_of_points = len(points)
        model_points.append(points[0:number_of_points_per_part])
        model_index = [iter_surface.f(points[number_of_points_per_part-1])]
        counter = 1
        while True:
            k = (1-intersection_between_divisions)*number_of_points_per_part*counter
            if int(k+number_of_points_per_part) >= number_of_points:
                model_points.append(points[int(k):])
                model_index.append(iter_surface.f(points[-1]))
                return model_points, model_index    
            else:
                model_points.append(points[int(k):int(k+number_of_points_per_part)])
                model_index.append(iter_surface.f(points[int(k+number_of_points_per_part)-1]))
            counter += 1
                            
        

    def compute_initial_point(self, model, iter_surface, trajectories):
        """
        The compute_initial_point computes the initial point to start the generating trajectories
        algorithm. If trajectories were loaded, the initial point is the last computed point projected
        Keyword arguments:
        point_on_surfaces -- initial point on both surfaces.
        iter_surface -- surface to be iterated, as mathtools.sphere.
        step -- it must be small, e.g. 1e-3. Otherwise the method will fail.
        """
                
        if len(trajectories)>0:
            last_computed_point = trajectories[-1][-1]
            iter_surface.find_iter(last_computed_point)
            iter_surface.update()
            return mathtools.curvepoint(model, iter_surface, last_computed_point[0:3])
        else:    
            initial_point = self._points[argmax(iter_surface.f_array(self._points))]
            iter_surface.findnextparallel(initial_point)
            return mathtools.curvepoint(model, iter_surface, initial_point[0:3])

    def draw_parallel(self, point_on_surfaces, model, iter_surface, step):
        """
        The draw_parallel generates one coating trajectory. The trajectory is
        the intersection between two surfaces: the blade model, and the surface
        to be iterated, e.g. spheres. The algorithm
        follows the marching method, documentation available in:
        http://www.mathematik.tu-darmstadt.de/~ehartmann/cdgen0104.pdf
        page 94 intersection surface - surface.
        The function stops when a closed curve is found, i.e., when it is found a
        point near the initial point (<step distance).

        This is a data-intensive computing and might freeze your computer.

        Keyword arguments:
        point_on_surfaces -- initial point on both surfaces.
        iter_surface -- surface to be iterated, as mathtools.sphere.
        step -- it must be small, e.g. 1e-3. Otherwise the method will fail.
        """
        
        initial_point = copy(point_on_surfaces)
        counter = 0
        trajectory = [point_on_surfaces]
        while True:
            tan = mathtools.surfaces_tangent(trajectory[-1], iter_surface)
            next_point_on_surfaces = mathtools.curvepoint(model,
                                                          iter_surface,
                                                          trajectory[-1][0:3]-tan*step)
            dP = abs(next_point_on_surfaces[0:3]-initial_point[0:3])
            if counter==0:
                if max(dP)<=step:counter+=1
            else:
                if max(dP)<=step:
                    try:
                        p0 = trajectory[0][0:3]; p1 = trajectory[1][0:3]; p2 = trajectory[2][0:3]
                        if (max(abs(p0-p1))<=step) and (max(abs(p0-p2))>=step):
                            trajectory.append(next_point_on_surfaces)
                            return trajectory
                    except IndexError:
                        raise IndexError('Step is too big and the function terminated soon.')
            trajectory.append(next_point_on_surfaces)  

    def generate_trajectories(self, model, iter_surface, step = 1e-3):
        """
        Method generate the coating trajectories. The trajectories are
        the intersection between two surfaces: the blade model, and the surface
        to be iterated, e.g. spheres. The algorithm
        follows the marching method, documentation available in:
        http://www.mathematik.tu-darmstadt.de/~ehartmann/cdgen0104.pdf
        page 94 intersection surface - surface.

        This is a data-intensive computing and might freeze your computer.

        Keyword arguments:
        iter_surface -- surface to be iterated, as mathtools.sphere.
        step -- it must be small, e.g. 1e-3. Otherwise the method will fail.
        directory_to_save -- directory to save the trajectories (default is './Blade/Trajectory')
        """

        trajectories = []
        if not issubclass(iter_surface.__class__, mathtools.IterSurface):
            raise IterSurfaceError()

        self._blade.SetTransform(eye(4))
        self._blade.SetTransform(dot(self._blade.GetTransform(),
                                     matrixFromAxisAngle([0, -self.turbine.config.environment.blade_angle, 0]))
                                 )
        point_on_surfaces = self.compute_initial_point(model, iter_surface, trajectories)
 
        while iter_surface.criteria():
            trajectories.append(self.draw_parallel(point_on_surfaces, model, iter_surface, step))
            p0=trajectories[-1][-1]
            iter_surface.update()
            point_on_surfaces = mathtools.curvepoint(model, iter_surface, p0[0:3])
        return trajectories
