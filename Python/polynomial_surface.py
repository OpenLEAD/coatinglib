from numpy import *
from mayavi import mlab
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import implicit
from scipy.linalg import solve
import tree_time
#from sympy import *

Rays1 = load('blade_sampling/blade_crop_fast.npz')
Rays1  = Rays1['array']
Rays2 = load('blade_sampling/blade_crop_fast2.npz')
Rays2  = Rays2['array']

rR = concatenate((Rays1,Rays2))


def vector2(x,y,z):
    a = [1,
         x, y, z,
         x**2, y**2, z**2, x*y, x*z, y*z
         ]
    return a

def vector3(x,y,z):
    a = [1,
         x, y, z,
         x**2, y**2, z**2, x*y, x*z, y*z, 
         x**3, y**3, z**3, x**2*y, x**2*z, y**2*x, y**2*z, z**2*x, z**2*y,
         x*y*z]
    return a

def vector4(x,y,z):
    return [1, z, z**2, z**3, z**4, y, y*z, y*z**2, y*z**3, y**2, y**2*z, y**2*z**2, y**3, y**3*z, y**4, x, x*z, x*z**2, x*z**3, x*y, x*y*z, x*y*z**2, x*y**2, x*y**2*z, x*y**3, x**2, x**2*z, x**2*z**2, x**2*y, x**2*y*z, x**2*y**2, x**3, x**3*z, x**3*y, x**4]

def vector5(x,y,z):
    return [1, z, z**2, z**3, z**4, z**5, y, y*z, y*z**2, y*z**3, y*z**4, y**2, y**2*z, y**2*z**2, y**2*z**3, y**3, y**3*z, y**3*z**2, y**4, y**4*z, y**5, x, x*z, x*z**2, x*z**3, x*z**4, x*y, x*y*z, x*y*z**2, x*y*z**3, x*y**2, x*y**2*z, x*y**2*z**2, x*y**3, x*y**3*z, x*y**4, x**2, x**2*z, x**2*z**2, x**2*z**3, x**2*y, x**2*y*z, x**2*y*z**2, x**2*y**2, x**2*y**2*z, x**2*y**3, x**3, x**3*z, x**3*z**2, x**3*y, x**3*y*z, x**3*y**2, x**4, x**4*z, x**4*y, x**5]

def dvector5(x,y,z):
    dx = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, z, z**2, z**3, z**4, y, y*z, y*z**2, y*z**3, y**2, y**2*z, y**2*z**2, y**3, y**3*z, y**4, 2*x, 2*x*z, 2*x*z**2, 2*x*z**3, 2*x*y, 2*x*y*z, 2*x*y*z**2, 2*x*y**2, 2*x*y**2*z, 2*x*y**3, 3*x**2, 3*x**2*z, 3*x**2*z**2, 3*x**2*y, 3*x**2*y*z, 3*x**2*y**2, 4*x**3, 4*x**3*z, 4*x**3*y, 5*x**4]
    dy = [0, 0, 0, 0, 0, 0, 1, z, z**2, z**3, z**4, 2*y, 2*y*z, 2*y*z**2, 2*y*z**3, 3*y**2, 3*y**2*z, 3*y**2*z**2, 4*y**3, 4*y**3*z, 5*y**4, 0, 0, 0, 0, 0, x, x*z, x*z**2, x*z**3, 2*x*y, 2*x*y*z, 2*x*y*z**2, 3*x*y**2, 3*x*y**2*z, 4*x*y**3, 0, 0, 0, 0, x**2, x**2*z, x**2*z**2, 2*x**2*y, 2*x**2*y*z, 3*x**2*y**2, 0, 0, 0, x**3, x**3*z, 2*x**3*y, 0, 0, x**4, 0]
    dz = [0, 1, 2*z, 3*z**2, 4*z**3, 5*z**4, 0, y, 2*y*z, 3*y*z**2, 4*y*z**3, 0, y**2, 2*y**2*z, 3*y**2*z**2, 0, y**3, 2*y**3*z, 0, y**4, 0, 0, x, 2*x*z, 3*x*z**2, 4*x*z**3, 0, x*y, 2*x*y*z, 3*x*y*z**2, 0, x*y**2, 2*x*y**2*z, 0, x*y**3, 0, 0, x**2, 2*x**2*z, 3*x**2*z**2, 0, x**2*y, 2*x**2*y*z, 0, x**2*y**2, 0, 0, x**3, 2*x**3*z, 0, x**3*y, 0, 0, x**4, 0, 0]
    a=[dx,dy,dz]
    return a

def vectorn(n):
#    x=Symbol('x')
#    y=Symbol('y')
#    z=Symbol('z')
    
    n=n+1
    a=[]
    for i in range(0,n):
        for j in range(0,n-i):
            for k in range(0,n-i-j):
                a.append(x**i*y**j*z**k)

    return a            
                        
            

def dvector3(x,y,z):
    dx = [0,
         1, 0, 0,
         2*x, 0, 0, y, z, 0,
         3*x**2, 0, 0, 2*x*y, 2*x*z, y**2, 0, z**2, 0,
         y*z]

    dy = [0,
          0, 1, 0,
          0, 2*y, 0, x, 0, z,
          0, 3*y**2, 0, x**2, 0, 2*y*x, 2*y*z, 0, z**2,
          x*z]

    dz = [0,
          0, 0, 1,
          0, 0, 2*z, 0, x, y,
          0, 0, 3*z**2, 0, x**2, 0, y**2, 2*z*x, 2*z*y,
          x*y]
    a = [dx,dy,dz]
    return a

def dvector4(x,y,z):
    dx = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, z, z**2, z**3, y, y*z, y*z**2, y**2, y**2*z, y**3, 2*x, 2*x*z, 2*x*z**2, 2*x*y, 2*x*y*z, 2*x*y**2, 3*x**2, 3*x**2*z, 3*x**2*y, 4*x**3]
    dy = [0, 0, 0, 0, 0, 1, z, z**2, z**3, 2*y, 2*y*z, 2*y*z**2, 3*y**2, 3*y**2*z, 4*y**3, 0, 0, 0, 0, x, x*z, x*z**2, 2*x*y, 2*x*y*z, 3*x*y**2, 0, 0, 0, x**2, x**2*z, 2*x**2*y, 0, 0, x**3, 0]
    dz = [0, 1, 2*z, 3*z**2, 4*z**3, 0, y, 2*y*z, 3*y*z**2, 0, y**2, 2*y**2*z, 0, y**3, 0, 0, x, 2*x*z, 3*x*z**2, 0, x*y, 2*x*y*z, 0, x*y**2, 0, 0, x**2, 2*x**2*z, 0, x**2*y, 0, 0, x**3, 0, 0]
    a=[dx,dy,dz]
    return a

def matrix(rays):
    m = []
    for ray in rays:
        m.append(vector3(ray[0],ray[1],ray[2]))
    return array(m)

def matrix2(rays):
    m = []
    S = []
    for ray in rays:
        a=vector5(ray[0],ray[1],ray[2])
        da=dvector5(ray[0],ray[1],ray[2])
        b = [a,da[0],da[1],da[2]]
        s = [0, ray[3], ray[4], ray[5]]
        S.extend(s)
        m.extend(b)
    return m, S

def polynomial_surface(rays):
    m = matrix(rays)
    M = dot(transpose(m),m)
    w, v = linalg.eig(M)
    v = transpose(v)
    return v[w==min(w)],min(w)

def polynomial_surface2(rays):
    m, S = matrix2(rays)
    A = dot(transpose(m),m)
    b = dot(transpose(m),S)
    x = solve(A, b)
    return x


##v, w = polynomial_surface(rR)
##v=v[0]

v=[]

##mlab.clf()
##x, y, z = mgrid[-2:0:50j, -3:-1:50j, -2:0:50j]
#values = x+y+z
#values = v[0]+v[1]*x+v[2]*y+v[3]*z+v[4]*x**2+v[5]*y**2+v[6]*z**2+v[7]*x*y+v[8]*x*z+v[9]*y*z+v[10]*x**3+v[11]*y**3+v[12]*z**3+v[13]*x**2*y+v[14]*x**2*z+v[15]*y**2*x+v[16]*y**2*z+v[17]*z**2*x+v[18]*z**2*y+v[19]*x*y*z
##mlab.contour3d(x, y, z, values, contours=[2])
##mlab.axes()

##x=-1.1819988005749862..-0.48598792067278029
##y=-3.5753923468214932..-1.0554489912321572
##z=-1.4602732920175647..-0.1333696140894175

def fn3(x,y,z):
    value0 = v[0]
    value1 = v[1]*x+v[2]*y+v[3]*z
    value2 = v[4]*x**2+v[5]*y**2+v[6]*z**2+v[7]*x*y+v[8]*x*z+v[9]*y*z
    value3 = v[10]*x**3+v[11]*y**3+v[12]*z**3+v[13]*x**2*y+v[14]*x**2*z+v[15]*y**2*x+v[16]*y**2*z+v[17]*z**2*x+v[18]*z**2*y+v[19]*x*y*z
    value = value0+value1+value2+value3
    return value

def fn2(x,y,z):
    value0 = v[0]
    value1 = v[1]*x+v[2]*y+v[3]*z
    value = value0+value1
    return value

def fn4(x,y,z):
    
    return dot(vector4(x,y,z),transpose(v))


def f4(x,y,z):
    return 0.0158409735155322*x**4 - 0.0455175545677363*x**3*y + 0.229708256297771*x**3*z - 0.0123881037382856*x**3 - 0.0404958971530159*x**2*y**2 + 0.0370550760761227*x**2*y*z - 0.577768075781302*x**2*y - 0.155580530030559*x**2*z**2 + 0.0362405492287055*x**2*z - 0.131409797907689*x**2 + 0.0202057385466043*x*y**3 + 0.0426930375374075*x*y**2*z + 0.0696863140265359*x*y**2 - 0.0505449727608702*x*y*z**2 + 0.225686420684704*x*y*z - 0.0384970823430439*x*y + 0.105768716128252*x*z**3 - 0.118067494638723*x*z**2 - 1.21492908245443*x*z - 0.18276378366575*x + 0.0166396304986572*y**4 + 0.00338711654878079*y**3*z + 0.216501859969158*y**3 + 0.0811414830256883*y**2*z**2 + 0.02769438794406*y**2*z + 0.942772232605902*y**2 - 0.0290415751360407*y*z**3 + 0.627763785746382*y*z**2 + 0.0884479216928731*y*z + 1.61784537809966*y + 0.0109771138288731*z**4 - 0.178048926480733*z**3 + 1.50165747078441*z**2 + 0.196218414740315*z + 0.76794299165963

def f5(x,y,z):
    return 0.0297659116455479*x**5 + 0.181062045560249*x**4*y - 0.128661031864172*x**4*z + 0.307369024428665*x**4 + 0.0187730062188317*x**3*y**2 - 0.352382984812728*x**3*y*z - 0.0165274665944978*x**3*y + 0.109314604473024*x**3*z**2 - 0.285345405421894*x**3*z - 0.0953176916491287*x**3 + 0.0408304246877191*x**2*y**3 + 0.00447133758293315*x**2*y**2*z + 0.195842840630828*x**2*y**2 + 0.27146542716375*x**2*y*z**2 + 0.233595210793662*x**2*y*z - 0.260349111071952*x**2*y - 0.0589021386012211*x**2*z**3 + 0.187272006426444*x**2*z**2 + 0.662311023454343*x**2*z - 0.100094166057317*x**2 - 0.0432138132832048*x*y**4 + 0.0400826659272909*x*y**3*z - 0.42553574962281*x*y**3 - 0.102140941691896*x*y**2*z**2 + 0.409711840025748*x*y**2*z - 1.565282909628*x*y**2 - 0.0334481817041418*x*y*z**3 - 0.733013704495247*x*y*z**2 + 1.27796791340411*x*y*z - 2.51747409763771*x*y - 0.020191776904239*x*z**4 + 0.203163668332096*x*z**3 - 1.18680358227468*x*z**2 - 0.369578257596795*x*z - 1.47497490230636*x - 0.00520797176161978*y**5 + 0.040846567888235*y**4*z - 0.0651785345972887*y**4 - 0.00582434668496111*y**3*z**2 + 0.429749743471307*y**3*z - 0.294262656848997*y**3 + 0.0302244925660516*y**2*z**3 + 0.00973115339906459*y**2*z**2 + 1.6524767421243*y**2*z - 0.600269372065211*y**2 - 0.0149898305108897*y*z**4 + 0.183269458830206*y*z**3 + 0.374946543237288*y*z**2 + 2.72744600165394*y*z - 0.577089308079262*y + 0.00729646367142089*z**5 - 0.065684684064366*z**4 + 0.185426421490161*z**3 + 1.25471056201527*z**2 + 1.66800538058776*z - 0.370560555814104

def main():
    global rR
    global v
    v = polynomial_surface2(rR)

    ax = implicit.plot_implicit(f5,bbox=(-5,5,-5,5,-5,5))
    ##
    #rR = tree_time.tree_time(rR,0.05)
    #rR=array(rR)
    Axes3D.scatter(ax,rR[:,0],rR[:,1],rR[:,2])
    ##
    plt.show()