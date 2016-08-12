from numpy import dot, array, sqrt, log, exp, vstack, load, savez_compressed, sum, zeros
from scipy.linalg import solve
from os import makedirs
import errno
from copy import copy


class RBF:
    """Radial basis function class for blade modeling.
       The kernel must be specified (string): r3, logr, gaussr or iqr.
       The name of the object must be specified (string): e.g. jiraublade.
       The points of the object to be modeled, must be specified if first time.
       RBF(name,kernel,points,eps)."""

    def __init__(self, name, kernel, points=[], eps=0.001):
        self._kernel = kernel
        self._phi_dict = {'r3':self._r3,'logr':self._logr,'gaussr':self._gaussr,'iqr':self._iqr}
        self._dphi_dict = {'r3':self._dr3,'logr':self._dlogr,'gaussr':self._dgaussr,'iqr':self._diqr}
        self._name = name+'_'+kernel
        self._eps = eps
        self._w = []
        self._points = points
        self.model_type = 'RBF'
        
    def _r3(self, ci, cj):
        c = cj-ci
        c = sqrt(sum(c*c,1))
        return c*c*c

    def _dr3(self, ci, cj):
        c = cj-ci
        nc = 3*sqrt(sum(c*c,1))
        b = zeros((len(c),3))
        b[:,0] = -nc*c[:,0]*self._w
        b[:,1] = -nc*c[:,1]*self._w
        b[:,2] = -nc*c[:,2]*self._w
        return sum(b,0)

    def _logr(self, ci, cj):
        c = ci-cj
        c2 = dot(c, c)
        if c2==0:
            return 0
        c = sqrt(c2)
        return c2*log(c)

    def _dlogr(self, ci, cj):
        c = ci-cj
        c2 = dot(c, c)
        if c2==0:
            return array([0,0,0])
        return array([c[0]*log(c2) + c[0],
                      c[1]*log(c2) + c[1],
                      c[2]*log(c2) + c[2]])

    def _gaussr(self, ci, cj):
        e = 0.1
        c = ci-cj
        c2 = e**2*dot(c, c)
        return exp(-c2)

    def _dgaussr(self, ci, cj):
        e = 0.1; e2 = e*e 
        k = -2*e2
        c = ci-cj
        c2 = e2*dot(c, c)
        oexp = exp(-c2)
        return array([k*c[0]*oexp, k*c[1]*oexp, k*c[2]*oexp])

    def _iqr(self, ci, cj):
        e = 0.001
        c = ci-cj
        c2 = e**2*dot(c, c)
        return 1/(1+c2)

    def _diqr(self, ci,cj):
        e2 = 0.001**2; k = -2*e2
        c = ci-cj
        c2 = 1/((e2*dot(c,c)+1)**2)
        return array([k*c[0]*c2,k*c[1]*c2, k*c[2]*c2])        

    def _phi(self, ci, cj):
        try:
            return self._phi_dict[self._kernel](ci, cj)
        except:
            raise #ValueError('kernel is not in list.')
        return

    def _dphi(self, ci,cj):
        try:
            return self._dphi_dict[self._kernel](ci, cj)
        except:
            raise #ValueError('kernel is not in list.')
        return

    def f(self, c):
        c = array(c[0:3])
        return dot(self._w, self._phi(c, self._points[:,0:3]))
        
    def df(self, c):
        c = array(c[0:3])
        return self._dphi(c, self._points[:,0:3])

    def _pointsaugment(self):
        b = copy(self._points)
        b[:,0:3]=b[:,0:3]+b[:,3:6]*self._eps
        self._points = vstack((self._points,b))   
        return    

    def make(self):
        print 'RBF::make -  Warning: this is a data-intensive computing and might freeze your computer.'
        if len(self._points)>0:
            self._pointsaugment()
            N = len(self._points)
            K = zeros((N, N))
            d = []
            for i in range(0,N):
                K[i,:] = self._phi(self._points[i,0:3],self._points[:,0:3])
                if i>=N/2:
                    d.append(self._eps)
                else: d.append(0)
            try:    
                self._w = solve(K,d)
            except:
                raise
        else:
            print "RBF::make - There are no points to make a RBF. Please create points with BladeModeling::sampling."
        
