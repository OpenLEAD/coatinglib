from numpy import load, savez_compressed, zeros, ones, arange
from numpy import meshgrid, array, shape, sum, eye, dot, argsort
from numpy import argmax, argmin, savetxt, mean, loadtxt, cross, linalg
from numpy import r_, c_, outer, tile, concatenate, sqrt, linspace
from os import makedirs
from os.path import join
import errno
from openravepy import RaveCreateCollisionChecker, matrixFromAxisAngle
import mathtools
from math import pi, ceil, isnan
from copy import copy, deepcopy
from lxml import etree as ET
import ast
from rbf import RBF

class BladeModeling:
    """ BladeModeling class for blade modelling.

    Keyword arguments:
    name -- the name of the blade.
    turbine -- turbine object.
    blade -- body blade object.
    """

    turbine = None
    blade = None
    
    def __init__(self, turbine, blade):
        self.turbine = turbine
        self._blade = blade
        self.trajectories = []
        self.points = []
        self.models = []
        self.models_index = []
        self.model_iter_surface = None
        self.samples_delta = turbine.config.model.samples_delta
        self.min_distance_between_points = turbine.config.model.min_distance_between_points
        self.intersection_between_divisions = turbine.config.model.intersection_between_divisions
        self.number_of_points_per_model = turbine.config.model.number_of_points_per_model
        self.trajectory_step = turbine.config.model.trajectory_step
        self.trajectory_iter_surface = None
        self.gap = turbine.config.coating.parallel_gap

    def save_samples(self, directory_to_save, name):
        """
        A method to save object samples and samples's info.

        Keyword arguments:
        delta -- the used uniform distance between the cube's samples parameter.
        min_distance_between_points -- the used uniform distance between the
        object's samples parameter.
        """

        delta = self.samples_delta
        min_distance_between_points = self.min_distance_between_points
        
        if delta is None:
            raise TypeError("Samples were not generated. samples_delta is None.")
            
        try:
            makedirs(directory_to_save)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise

        samples = ET.Element("samples")
        ET.SubElement(samples, "name").text = name
        ET.SubElement(samples, "number_of_samples").text = str(len(self.points))
        ET.SubElement(samples, "delta").text = str(delta)
        ET.SubElement(samples, "min_distance_between_points").text = str(min_distance_between_points)
        ET.SubElement(samples, "path").text = join(directory_to_save, 'samples.csv')
            
        savetxt(join(directory_to_save, 'samples.csv'), self.points, delimiter = ',')

        tree = ET.ElementTree(samples)
        tree.write(join(directory_to_save, "samples.xml"), pretty_print=True)
        
        return

    def save_model(self, directory_to_save, name):
        """
        A method to save model and model info

        Keyword arguments:
        models -- the computed models.
        model_index -- when to switch between models.
        iter_surface -- the iter surface used to compute the trajectories.
        directory_to_save -- where to save the files.
        """

        models = self.models
        iter_surface = self.model_iter_surface
        model_index = self.models_index
        intersection_between_divisions = self.intersection_between_divisions
        number_of_points_per_model = self.number_of_points_per_model

        try:
            makedirs(directory_to_save)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise

        model = ET.Element("model")
        ET.SubElement(model, "name").text = name

        surface = ET.SubElement(model, "iter_surface")
        if iter_surface is not None:   
            ET.SubElement(surface, "type").text = iter_surface.name()
            ET.SubElement(surface, "Rn").text = str(iter_surface._Rn)
            ET.SubElement(surface, "switching_parameters").text = str(model_index)
        else: ET.SubElement(surface, "type").text = 'None'

        ET.SubElement(model, "intersection_between_divisions").text = str(intersection_between_divisions)
        ET.SubElement(model, "number_of_points_per_model").text = str(number_of_points_per_model)

        for i in range(0,len(models)):
            doc = ET.SubElement(model, "interpolation")
            ET.SubElement(doc, "type").text = models[i].model_type
            ET.SubElement(doc, "points").text = join(directory_to_save,'points_' + str(i) + ".csv")
            savetxt(join(directory_to_save, 'points_' + str(i) + '.csv'), models[i]._points, delimiter = ',')
            
            if models[i].model_type == 'RBF':
                ET.SubElement(doc, "w").text = join(directory_to_save, 'w_' + str(i) + ".csv")
                savetxt(join(directory_to_save, 'w_' + str(i) + '.csv'), models[i]._w, delimiter = ',')
                ET.SubElement(doc, "kernel").text = models[i]._kernel
                ET.SubElement(doc, "eps").text = str(models[i]._eps)
                if models[i]._kernel=='gaussr': ET.SubElement(doc, "gauss_parameter").text = str(models[i].gausse)
                
  
        tree = ET.ElementTree(model)
        tree.write(join(directory_to_save, "model.xml"), pretty_print=True)
        return

    def save_trajectory(self, xml_model, directory_to_save, name):
        """
        A method to save trajectory and trajectory info. Two type of files are saved:
        .npz and n .csv files

        Keyword arguments:
        trajectories -- the computed trajectories.
        xml_model -- path to xml model file.
        iter_surface -- the iter surface used to compute the trajectories.
        trajectory_step -- the tangent step used.
        directory_to_save -- where to save the files.
        """

        trajectories = self.trajectories
        trajectory_step = self.trajectory_step
        iter_surface = self.trajectory_iter_surface

        try:
            makedirs(directory_to_save)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise

        trajectory = ET.Element("trajectory")
        ET.SubElement(trajectory, "name").text = name
        ET.SubElement(trajectory, "trajectory_step").text = str(trajectory_step)
        
        ET.SubElement(trajectory, "xml_model").text = xml_model
        
        surface = ET.SubElement(trajectory, "iter_surface")
        ET.SubElement(surface, "type").text = iter_surface.name()
        ET.SubElement(surface, "Rn").text = str(iter_surface._Rn0)
        ET.SubElement(surface, "stopR").text = str(iter_surface.stopR)
        ET.SubElement(surface, "coatingstep").text = str(iter_surface.coatingstep)

        csv_files = ET.SubElement(trajectory, "csv_files")
        iter_surface._Rn = iter_surface._Rn0
        for i in range(0,len(trajectories)):
            ETfile = ET.SubElement(csv_files, "file")
            ET.SubElement(ETfile, "path").text = join(directory_to_save,'trajectory_' + str(i) + '.csv')
            iter_surface.update()
            ET.SubElement(ETfile, "iter_R").text = str(iter_surface._Rn)
            savetxt(join(directory_to_save, 'trajectory_' + str(i) + '.csv'), trajectories[i], delimiter = ',')

        ET.SubElement(trajectory, "npz_file").text = join(directory_to_save,'trajectory.npz')
        savez_compressed(join(directory_to_save,'trajectory.npz'), array=trajectories)
        
        tree = ET.ElementTree(trajectory)
        tree.write(join(directory_to_save, "trajectory.xml"), pretty_print=True)
        return

    def load_samples(self, xml_file):
        """
        A method to load the samples and the used sampling configurations.

        Keyword arguments:
        xml_file -- the xml file path of the samples.
        """
        
        try:
            xml = ET.parse(open(xml_file))
        except IOError:
            raise IOError('Samples could not be loaded. Call method BladeModeling::sampling')

        self.points = loadtxt(xml.find('path').text, delimiter=',')
        if len(self.points.shape)!=2 or self.points.shape[1]!=6:
            raise TypeError("samples's shape must be (n,6). [x,y,z, nx,ny,nz]")

        self.samples_delta = float(xml.find('delta').text)
        self.min_distance_between_points = float(xml.find('min_distance_between_points').text)
        return

    def load_model(self, xml_file):
        """
        A method to load the models and the used modeling configurations.

        Keyword arguments:
        xml_file -- the xml file path of the model.
        """
        
        try:
            xml = ET.parse(open(xml_file))
        except IOError:
            raise IOError('Model could not be loaded. Call method BladeModeling::make_model')

        if xml.find('iter_surface').find('type').text != 'None':
            self.model_iter_surface = getattr(mathtools, xml.find('iter_surface').find('type').text)(float(xml.find('iter_surface').find('Rn').text),0,0)
            self.models_index = ast.literal_eval(xml.find('iter_surface').find('switching_parameters').text)
        else: self.model_iter_surface = None

        self.intersection_between_divisions = float(xml.find('intersection_between_divisions').text)
        self.number_of_points_per_model = int(xml.find('number_of_points_per_model').text)

        self.models = []
        for xml_model in xml.findall('interpolation'):
            if xml_model.find('type').text == 'RBF':
                points = loadtxt(xml_model.find('points').text, delimiter=',')
                eps = float(xml_model.find('eps').text)
                kernel = xml_model.find('kernel').text
                gausse = None
                if xml_model.find("gauss_parameter") is not None: gausse = float(xml_model.find("gauss_parameter").text)
                model = RBF(kernel, points, eps, gausse)
                model._w = loadtxt(xml_model.find('w').text, delimiter=',')
                self.models.append(model)
            else: raise TypeError('This type of interpolation was not yet implemented')
        return
    
    def load_trajectory(self, xml_file):
        """
        A method to load the trajectories and the used trajectory configurations.

        Keyword arguments:
        xml_file -- the xml file path of the trajectories.
        Note that the load_trajectories calls the load_model method.
        """
        
        try:
            xml = ET.parse(open(xml_file))
        except IOError:
            raise IOError('Trajectory could not be loaded. Call method BladeModeling::generate_trajectories')

        self.trajectory_step = float(xml.find('trajectory_step').text)

        self.trajectory_iter_surface = getattr(mathtools, xml.find('iter_surface').find('type').text)(float(xml.find('iter_surface').find('Rn').text),
                                                                                                      float(xml.find('iter_surface').find('stopR').text),
                                                                                                      float(xml.find('iter_surface').find('coatingstep').text)
                                                                                                      )
        self.load_model(xml.find('xml_model').text)
        self.trajectories = (load(xml.find('npz_file').text)['array']).tolist()
        return

    def validate_models(self, models, samples):
        """
        After the load_model call, the user may not know if the models are valid.
        This method verifies if the model are valid w.r.t. the samples. 

        Keyword arguments:
        models -- the model objects (e.g. RBF).
        samples -- samples of the object.
        """

        return
            
            

    def sampling(self):
        """
        The sampling method is an algorithm for uniformly sampling objects.
        It computes the bounding box of the object (object inscribed by cube), uniformly
        samples the box's faces, and checks rays collision from cube's samples to the object.
        Rays with collision and its normal vectors are stored as object's samples.
        Finally, a kdtree is computed from object's samples and neighbooring points are removed  
        to generate an uniformly sampling on the object (filtering data by distance).

        This is a data-intensive computing and might freeze your computer.
        """

        delta = self.samples_delta
        min_distance_between_points = self.min_distance_between_points
        
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
                     (0,0,-e[2],0,0,1,e[0],0,0,0,e[1],0), #-z
                     (0,0,e[2],0,0,-1,e[0],0,0,0,e[1],0), #z
                     (0,e[1],0,0,-1,0,e[0],0,0,0,0,e[2]), #y
                     (0,-e[1],0,0,1,0,e[0],0,0,0,0,e[2])  #-y
                     ))
        maxlen = 2*sqrt(sum(e**2))+0.03
        self.points = zeros((0,6))
        for side in sides:
            ex = sqrt(sum(side[6:9]**2))
            ey = sqrt(sum(side[9:12]**2))
            XX,YY = meshgrid(r_[arange(-ex,-0.25*delta,delta),0,arange(delta,ex,delta)],
                             r_[arange(-ey,-0.25*delta,delta),0,arange(delta,ey,delta)])
            localpos = outer(XX.flatten(),side[6:9]/ex)+outer(YY.flatten(),side[9:12]/ey)
            N = localpos.shape[0]
            rays = c_[tile(p+side[0:3],(N,1))+localpos,
                      maxlen*tile(side[3:6],(N,1))]
            collision, info = self.turbine.env.CheckCollisionRays(rays, self._blade)
            # make sure all normals are the correct sign: pointing outward from the object)
            newinfo = info[collision,:]
            if len(newinfo) > 0:
                  newinfo[sum(rays[collision,3:6]*newinfo[:,3:6],1)>0,3:6] *= -1
                  self.points = r_[self.points,newinfo]

        self.points = mathtools.filter_by_distance(self.points, min_distance_between_points)
        self.samples_delta = delta
        self.min_distance_between_points = min_distance_between_points
        return 
        
    def make_model(self, model, model_iter_surface = None):
        """
        The make_model method is an algorithm to generate a mathematical representation
        of the object (mesh). This method can be called after the sampling method,
        and never before. The model is generated by the model_type (constructor), e.g. RBF.

        This is a data-intensive computing and might freeze your computer.

        Keyword arguments:
        model -- model type for the implicit computation, e.g. the RBF model (rbf.py).
        iter_surface -- This is the surface to be iterated, as mathtools.sphere.
        If the model has more than 9000 points, this argument is needed.
        """

        number_of_points_per_model = self.number_of_points_per_model
        intersection_between_divisions = self.intersection_between_divisions

        model_points = [self.points]
        models_index = []
        models = []
        
        if len(self.points)>7300:
            if model_iter_surface is None:
                raise ValueError("The number of points is "+str(len(self.points))+""", which is
                                bigger than 9000. It is not safe to make an unique model this big.
                                Create an iterate surface to partitionate the model
                                (check mathtools.IterSurface).""")
            elif not issubclass(model_iter_surface.__class__, mathtools.IterSurface):
                    raise TypeError("Object is not a valid surface.")
            else:
                model_points, models_index = self.divide_model_points(self.points, model_iter_surface,
                                                                     number_of_points_per_model,
                                                                     intersection_between_divisions)                

        for i in range(0,len(model_points)):
            temp_model = deepcopy(model)
            temp_model._points = model_points[i]
            temp_model.make()
            models.append(temp_model)

        self.models = models
        self.models_index = models_index
        self.model_iter_surface = model_iter_surface
        self.intersection_between_divisions = intersection_between_divisions
        return models, models_index

    def divide_model_points(self, points, iter_surface, number_of_points_per_part,
                            intersection_between_divisions):
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

        if intersection_between_divisions<=0:
            raise ValueError('intersection_between_divisions must be bigger than 0')

        if number_of_points_per_part<1000:
            raise ValueError('number_of_points_per_division must be bigger than 1000')
        
        model_points = []
        points_distance = iter_surface.f_array(points)
        points = points[argsort(points_distance)]
        number_of_points = len(points)
        model_points.append(points[0:number_of_points_per_part])
        model_index = [(iter_surface.f(points[0]),
                        iter_surface.f(points[number_of_points_per_part-1]))]
        counter = 1
        while True:
            k = (1-intersection_between_divisions)*number_of_points_per_part*counter
            if int(k+number_of_points_per_part) >= number_of_points:
                model_points.append(points[int(k):])
                model_index.append((iter_surface.f(points[int(k)]),
                                    iter_surface.f(points[-1])))
                return model_points, model_index    
            else:
                model_points.append(points[int(k):int(k+number_of_points_per_part)])
                model_index.append((iter_surface.f(points[int(k)]),
                                    iter_surface.f(points[int(k+number_of_points_per_part)-1])))
            counter += 1

    def compute_initial_point(self, iter_surface, trajectories):
        """
        The compute_initial_point computes the initial point to start the generating trajectories
        algorithm. If trajectories were loaded, the initial point is the last computed point projected
        
        Keyword arguments:
        point_on_surfaces -- initial point on both surfaces.
        iter_surface -- surface to be iterated, as mathtools.sphere.
        step -- it must be small, e.g. 1e-3. Otherwise the method will fail.
        """

        point = []
        if len(trajectories)>0:
            last_computed_point = trajectories[-1][-1]
            iter_surface.find_iter(last_computed_point)
            iter_surface.update()
            point = last_computed_point
        else:    
            initial_point = self.points[argmax(iter_surface.f_array(self.points))]
            iter_surface.findnextparallel(initial_point)
            point = initial_point
            
        model = self.select_model(point)
        return mathtools.curvepoint(model, iter_surface, point[0:3])

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

    def draw_meridian_by_point(self, origin, meridian_step = 0.001):
        
        iter_surface = self.trajectory_iter_surface
        stopR = iter_surface.stopR
        meridian = []

        # Going up
        point = array(origin)
        iter_surface.find_iter(point)
        iter_surface.coatingstep = meridian_step
        while iter_surface.criteria():
            model = self.select_model(point)
            df = model.df(point)
            df = df/linalg.norm(df)
            point[3:6] = df
            meridian += [point]
            if abs(dot(df, point[0:3]/linalg.norm(point[0:3]))) >= 1-1e-1:
                break
            tan = cross(point[0:3]/linalg.norm(point[0:3]), df)
            tan = cross(df, tan)
            iter_surface.update()
            point[0:3] = point[0:3] + tan*meridian_step
            model = self.select_model(point)
            point = mathtools.curvepoint(model,
                                         iter_surface,
                                         point[0:3])

        
        #Going down
        point = array(origin)
        iter_surface.find_iter(point)
        iter_surface.coatingstep = -meridian_step
        iter_surface.stopR = iter_surface._Rn0
        while not iter_surface.criteria():
            model = self.select_model(point)

            df = model.df(point)
            df = df/linalg.norm(df)
            point[3:6] = df
            meridian += [point]
            if abs(dot(df, point[0:3]/linalg.norm(point[0:3]))) >= 1-1e-1:
                break
            tan = cross(point[0:3]/linalg.norm(point[0:3]), df)
            tan = cross(df, tan)
            iter_surface.update()
            point[0:3] = point[0:3] + tan*meridian_step
            model = self.select_model(point)
            point = mathtools.curvepoint(model,
                                         iter_surface,
                                         point[0:3])

            
        iter_surface.stopR = stopR
                
        return meridian
        

    def draw_meridians(self, parallel, meridian_step, parallel_step):
        
        meridians = []

        for i in linspace(0, len(parallel), parallel_step).astype(int)[:-1]:
            point = parallel[i]
            meridian = self.draw_meridian_by_point(point, meridian_step)
            meridians += [meridian]
            
        return meridians

    def draw_meridians_old(self, parallel, miridian_step, parallel_step, parallel_init = 0):
        
        meridians = []
        iter_surface = self.trajectory_iter_surface
        stopR = iter_surface.stopR
        Rn0 = iter_surface._Rn0

        for i in linspace(parallel_init, len(parallel), parallel_step).astype(int)[:-1]:
            meridian = []
            point = parallel[i]
            iter_surface.find_iter(point)
            iter_surface.coatingstep = miridian_step
            iter_surface.stopR = stopR
            while iter_surface.criteria():
                model = self.select_model(point)
                df = model.df(point)
                df = df/linalg.norm(df)
                point[3:6] = df
                meridian.append(point)
                if abs(dot(df, point[0:3]/linalg.norm(point[0:3]))) >= 1-1e-1:
                    break
                tan = cross(point[0:3]/linalg.norm(point[0:3]), df)
                tan = cross(df, tan)
                iter_surface.update()
                point[0:3] = point[0:3] + tan*miridian_step
                model = self.select_model(point)
                point = mathtools.curvepoint(model,
                                             iter_surface,
                                             point[0:3])

            
            point = parallel[i]
            iter_surface.find_iter(point)
            iter_surface.coatingstep = -miridian_step
            iter_surface.stopR = Rn0
            while not iter_surface.criteria():
                model = self.select_model(point)

                df = model.df(point)
                df = df/linalg.norm(df)
                point[3:6] = df
                meridian.append(point)
                if abs(dot(df, point[0:3]/linalg.norm(point[0:3]))) >= 1-1e-1:
                    break
                tan = cross(point[0:3]/linalg.norm(point[0:3]), df)
                tan = cross(df, tan)
                iter_surface.update()
                point[0:3] = point[0:3] + tan*miridian_step
                model = self.select_model(point)
                point = mathtools.curvepoint(model,
                                             iter_surface,
                                             point[0:3])
            meridians.append(meridian)
        return meridians

    def generate_trajectories(self, iter_surface):
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
        """

        step = self.trajectory_step

        trajectories = deepcopy(self.trajectories)
        
        if not issubclass(iter_surface.__class__, mathtools.IterSurface):
            raise TypeError("Object is not a valid surface.")

        if self.model_iter_surface is not None:
            if iter_surface.name() != self.model_iter_surface.name():
                raise TypeError("Object iter_surface must have same type of model_iter_surface.")

        if not self.models:
            raise IndexError("Object models is empty. Load or create a model before generate trajectories.")

        self._blade.SetTransform(eye(4))
        self._blade.SetTransform(dot(self._blade.GetTransform(),
                                     matrixFromAxisAngle([0, -self.turbine.config.environment.blade_angle, 0]))
                                 )

        point_on_surfaces = self.compute_initial_point(iter_surface, trajectories)
        self.trajectory_iter_surface = iter_surface
 
        while iter_surface.criteria():
            model = self.select_model(point_on_surfaces)
            trajectories.append(self.draw_parallel(point_on_surfaces, model, iter_surface, step))
            p0=trajectories[-1][-1]
            iter_surface.update()
            point_on_surfaces = mathtools.curvepoint(model, iter_surface, p0[0:3])
            self.trajectories = trajectories
        return trajectories

    def select_model(self, point):
        """
        From a list of models generated for a highly sampled object,
        this method selects the model that interpolates the parallel to be computed.

        Keyword arguments:
        point -- a point of the parallel.
        """
        
        models = array(self.models)
        models_index = array(self.models_index)
        intersection_between_divisions = self.intersection_between_divisions
        model_iter_surface = self.model_iter_surface
        
        if len(models) == 1:
            return models[0]

        f_point = model_iter_surface.f(point)
        models = models[(models_index[:,0] <= f_point+1e-5) & (models_index[:,1] >= f_point-1e-5)]

        if len(models) == 1:
            return models[0]
        
        if len(models)==0:
            raise ValueError("Point:"+str(point)+", with value:"+str(f_point)+", is not in range of the RBF models.\n "+str(self.models_index))

        
        f_points = float('Inf')
        for i in range(0,len(models)):
            dif = abs(sum(model_iter_surface.f_array(models[i]._points)-f_point)/len(models[i]._points))
            if dif < f_points:
                f_points = dif
            if dif > f_points:
                return models[i-1]
        return models[-1]

    def filter_trajectory(self, interpolation='linear', new_step = None, max_error = None):

        if interpolation == 'linear':
            f = mathtools.distance_point_line_3d

        if new_step is None:
            new_step = 20*self.trajectory_step

        if new_step < self.trajectory_step:
            return self.trajectory_step

        if max_error is None:
            max_error = self.gap/2
        
        D = []
        new_trajectories = []
        input_step = new_step
        for trajectory in self.trajectories:
            d = 100
            while d > max_error:
                d_list = []
                index_factor = int(new_step/self.trajectory_step)
                index_selection = arange(0, len(trajectory), index_factor)
                for i in range(0,len(index_selection)):
                    if i==len(index_selection)-1:
                        for j in range(index_selection[i], len(trajectory)):
                            d_list.append(f(array(trajectory[index_selection[i]])[0:3],
                                            array(trajectory[index_selection[0]])[0:3],
                                            array(trajectory[j])[0:3]))
                    else:
                        for j in range(index_selection[i], index_selection[i+1]):
                            d_list.append(f(array(trajectory[index_selection[i]])[0:3],
                                       array(trajectory[index_selection[i+1]])[0:3],
                                       array(trajectory[j])[0:3]))
                d = max(d_list)
                new_step = new_step - self.trajectory_step
            D.append(d_list)
            new_trajectories.append((array(trajectory)[arange(0, len(trajectory), index_factor)]).tolist())
            new_step = input_step
        return new_trajectories, D

    def filter_trajectory_opt(self, interpolation='linear', max_error = None):

        if interpolation == 'linear':
            f = mathtools.distance_point_line_3d

        if max_error is None:
            max_error = self.gap/2
        
        new_trajectories = []
        N = len(self.trajectories)
        counter = 0
        for trajectory in self.trajectories:
            i = 0
            new_trajectory = []
            while i < len(trajectory):
                new_trajectory.append(trajectory[i])
                temp_i=i+1
                for j in range(i+2, len(trajectory)):
                    for k in range(i,j):
                        d = f(array(trajectory[i])[0:3],
                              array(trajectory[j])[0:3],
                              array(trajectory[k])[0:3])
                        if d >= max_error:
                            break
                    else:
                        temp_i = j
                        continue
                    break
                i = temp_i
            new_trajectories.append(new_trajectory)
            counter+=1
            print('iter %3i / %3i' % (counter,N))
        return new_trajectories

    def rotate_models(self, T):
        for model in self.models:
            model._points = array(mathtools.rotate_trajectories([model._points], T)[0])
        return

    def find_borders(self, init_parallel, end_parallel):
        neg_border = []
        pos_border = []
        for i in range(init_parallel, end_parallel):
            rays = mathtools.filter_by_distance(array(self.trajectories[i]), None, 0.9995, True)
            neg_border.append(rays[rays[:,1]<0])
            pos_border.append(rays[rays[:,1]>0])
        return neg_border, pos_border
