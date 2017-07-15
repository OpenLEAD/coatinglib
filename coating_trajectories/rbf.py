from numpy import dot, array, sqrt, log, exp, vstack
from numpy import sum, zeros, ndarray, linalg
from copy import copy

## @file
# @brief This contains functions to make the RBF model
# @author Renan S. Freitas
# @bug No known bugs


##@page RBF
class RBF:
    """
    Radial basis function class for blade modeling.
    The RBF interpolation is an implicit reconstruction of an object.
    The algorithm is described in pages 35-40 of course available in link:
    http://graphics.stanford.edu/courses/cs468-10-fall/LectureSlides/04_Surface_Reconstruction.pdf
    Or in paper: Reconstruction and Representation of 3D Objects with Radial Basis
    Functions, J.C. Carr, R.K. Beatson.

    Keyword arguments: 
    kernel -- radial basis functions: r3, logr, gaussr. (r = abs(ci-cj))
    points -- of the object to be modeled.
    """

    def __init__(self, kernel, points=[], eps=1e-3, gausse = 1e-1):

        if not isinstance(kernel, basestring):
            raise ValueError('kernel is not a string')

        if not isinstance(points, ndarray):
            if isinstance(points, list): points = array(points)
            else: raise ValueError('points is not a list or numpy.ndarray')
        
        self._kernel = kernel
        self._phi_dict = {'r3':self._r3,'logr':self._logr,'gaussr':self._gaussr}
        self._dphi_dict = {'r3':self._dr3,'logr':self._dlogr,'gaussr':self._dgaussr}
        self._eps = eps
        self._w = []
        self._points = points
        self.model_type = 'RBF'
        self.gausse = 1e-1
       
        try:
            self._phi_dict[self._kernel]
        except KeyError:
            raise
        
    def _r3(self, ci, cj):
        c = cj-ci
        c = sqrt(sum(c*c,1))
        return c*c*c

    def _dr3(self, ci, cj):
        c = -(cj-ci)
        nc = 3*sqrt(sum(c*c,1))
        b = zeros((len(c),3))
        b[:,0] = nc*c[:,0]*self._w
        b[:,1] = nc*c[:,1]*self._w
        b[:,2] = nc*c[:,2]*self._w
        return sum(b,0)

    def _logr(self, ci, cj):
        c = -(cj-ci)
        c2 = sum(c*c,1)
        logc = log(sqrt(c2))
        logc[logc==-float('Inf')]=0
        return c2*logc

    def _dlogr(self, ci, cj):
        c = -(cj-ci)
        logc2 = log(sum(c*c,1))
        logc2[logc2==-float('Inf')]=0

        b = zeros((len(c),3))
        b[:,0] = (logc2*c[:,0]+c[:,0])*self._w
        b[:,1] = (logc2*c[:,1]+c[:,1])*self._w
        b[:,2] = (logc2*c[:,2]+c[:,2])*self._w
        return sum(b,0)

    def _gaussr(self, ci, cj):
        c = -(cj-ci)
        e = self.gausse
        return exp(-e*e*sum(c*c,1))

    def _dgaussr(self, ci, cj):
        e = self.gausse; e2 = e*e 
        k = -2*e2
        c = -(cj-ci)
        expc2 = exp(-e2*sum(c*c,1))

        b = zeros((len(c),3))
        b[:,0] = (expc2*c[:,0]*k)*self._w
        b[:,1] = (expc2*c[:,1]*k)*self._w
        b[:,2] = (expc2*c[:,2]*k)*self._w
        return sum(b,0)

    def _phi(self, ci, cj):
        return self._phi_dict[self._kernel](ci, cj)

    def _dphi(self, ci,cj):
        return self._dphi_dict[self._kernel](ci, cj)

    def f(self, c):
        """
        f(c) = sum(phi(ci-cj)), i samples of the object. It computes the RBF.
        If f(c) == 0, the point belongs to the object.
        If f(c) < 0, the point is inside the object.
        If f(c) > 0, the point is outside the object.

        Keyword arguments: 
        c -- point to evaluate the function.
        """

        c = array(c[0:3])
        return dot(self._w, self._phi(c, self._points[:,0:3]))
        
    def df(self, c):
        """
        It computes the gradient of the RBF.

        Keyword arguments:
        c -- point to evaluate the gradient of the function.
        """
        c = array(c[0:3])
        return self._dphi(c, self._points[:,0:3])

    def _pointsaugment(self, points):
        b = copy(points)
        b[:,0:3]=b[:,0:3]+b[:,3:6]*self._eps   
        return b   

    def make(self):
        """
        Make the implicit reconstruction with RBF.

        This is a data-intensive computing and might freeze your computer.
        """
        if len(self._points)>0:
            self._points = vstack((self._points,
                                   self._pointsaugment(self._points)))
            N = len(self._points)
            K = zeros((N, N))
            d = []
            for i in range(0,N):
                K[i,:] = self._phi(self._points[i,0:3], self._points[:,0:3])
                if i>=N/2:
                    d.append(self._eps)
                else: d.append(0)   
            self._w = linalg.solve(K,d)
        else:
            raise ValueError("There are no points to make a RBF.")
        
