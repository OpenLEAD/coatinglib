from numpy import *
from mayavi import mlab
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import implicit
from scipy.linalg import solve
#import tree_time
#from sympy import *

Rays1 = load('coated_points/pos4_black.npz')
Rays1  = Rays1['array']
Rays2 = load('coated_points/pos4_blue.npz')
Rays2  = Rays2['array']

#rR=Rays1
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
        a=vector4(ray[0],ray[1],ray[2])
        da=dvector4(ray[0],ray[1],ray[2])
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
v = array([ 6.01287452e-01, 4.38723449e-01, 1.43084162e+00,
                -1.62499201e-01, 9.11754926e-04, 1.35783742e+00,
                4.18625568e-01, 5.48893964e-01, -2.88737812e-02,
                8.21810883e-01, 1.60021716e-01, 6.75498267e-02,
                2.02815593e-01, 2.10435937e-02, 1.76268132e-02,
                -4.42642101e-01, -1.11486732e+00, -1.44688518e-01,
                7.06505573e-02, -4.29640954e-01, 2.66671401e-01,
                -6.87355941e-02, -1.30855344e-01, 4.67678861e-02,
                -1.16512312e-02, -1.85217968e-01, 5.48742656e-02,
                -1.15640905e-01, -6.01464874e-01, 3.69218793e-02,
                -3.81872876e-02, 1.38044769e-03, 1.51381408e-01,
                -4.33149498e-02, 3.88079293e-02])

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

def f4behind_filtered(x,y,z):
    return 0.23256524000168*x**4 - 0.0803593464907042*x**3*y - 0.383734254134677*x**3*z + 0.00436393131989541*x**3 + 0.0507931779830319*x**2*y**2 + 0.0634094586409759*x**2*y*z - 0.230454925983146*x**2*y + 0.492399768743622*x**2*z**2 - 0.0841779752966015*x**2*z - 0.0719708729835943*x**2 - 0.0173621871966343*x*y**3 - 0.0075263728926647*x*y**2*z - 0.0648891554000706*x*y**2 - 0.0824741231493418*x*y*z**2 + 0.101379769705647*x*y*z - 0.0305238793240705*x*y - 0.201578287070596*x*z**3 - 0.0138180109505892*x*z**2 - 0.909489575895646*x*z - 0.159044612063701*x + 0.0232495546616985*y**4 - 0.00512568852492615*y**3*z + 0.255063023472764*y**3 + 0.0832234479565902*y**2*z**2 - 0.0839605876912045*y**2*z + 1.01291311017395*y**2 - 0.0275741095916358*y*z**3 + 0.590020095137534*y*z**2 - 0.250204195750542*y*z + 1.67203834312611*y + 0.0716075331987429*z**4 - 0.237276989300336*z**3 + 1.30487684236424*z**2 + 0.0472995776370099*z + 0.750897097381243

def f4behind_3mm(x,y,z):
    return -0.0194242633552399*x**4 + 0.341146158309173*x**3*y - 0.0933102548831795*x**3*z - 0.0943235598570498*x**3 + 0.263214686468164*x**2*y**2 - 0.358423951545163*x**2*y*z + 0.920029193808266*x**2*y + 0.344663121655037*x**2*z**2 + 0.532945098232743*x**2*z + 0.58340150810563*x**2 - 0.0331170392490096*x*y**3 - 0.14886096292297*x*y**2*z - 0.0767252298546826*x*y**2 - 0.0311011259699964*x*y*z**2 - 0.5043143718347*x*y*z + 0.156603097998775*x*y - 0.189876806925384*x*z**3 - 0.698038218804208*x*z**2 - 0.603034682066102*x*z - 0.113579708087724*x + 0.00766169965664368*y**4 - 0.000852081771972584*y**3*z + 0.0797448256363611*y**3 + 0.0722629763165681*y**2*z**2 - 0.0703089122856204*y**2*z + 0.321256788562008*y**2 - 0.010854988146435*y*z**3 + 0.43152774654703*y*z**2 - 0.178331727877306*y*z + 0.532913320217161*y + 0.0511681485098616*z**4 - 0.0320708568210133*z**3 + 0.794314248916288*z**2 + 0.289653534545425*z + 0.06112422436045

def f4all_3mm(x,y,z):
    return 0.0388079293*x**4 - 0.0433149498*x**3*y + 0.151381408*x**3*z + 0.00138044769*x**3 - 0.0381872876*x**2*y**2 + 0.0369218793*x**2*y*z - 0.601464874*x**2*y - 0.115640905*x**2*z**2 + 0.0548742656*x**2*z - 0.185217968*x**2 - 0.0116512312*x*y**3 + 0.0467678861*x*y**2*z - 0.130855344*x*y**2 - 0.0687355941*x*y*z**2 + 0.266671401*x*y*z - 0.429640954*x*y + 0.0706505573*x*z**3 - 0.144688518*x*z**2 - 1.11486732*x*z - 0.442642101*x + 0.0176268132*y**4 + 0.0210435937*y**3*z + 0.202815593*y**3 + 0.0675498267*y**2*z**2 + 0.160021716*y**2*z + 0.821810883*y**2 - 0.0288737812*y*z**3 + 0.548893964*y*z**2 + 0.418625568*y*z + 1.35783742*y + 0.000911754926*z**4 - 0.162499201*z**3 + 1.43084162*z**2 + 0.438723449*z + 0.601287452

def f4coated4_50mm(x,y,z):
    return 0.210366913650167*x**4 + 0.607743786899794*x**3*y - 0.707234609953402*x**3*z + 0.573075139580213*x**3 + 0.203164874363018*x**2*y**2 - 0.421908589577347*x**2*y*z + 1.31121211241278*x**2*y + 1.24086914637841*x**2*z**2 + 0.065885415118794*x**2*z + 1.09974168847929*x**2 - 0.014150506173351*x*y**3 - 0.125701424698597*x*y**2*z + 0.037300424487786*x*y**2 - 0.0419025582651391*x*y*z**2 - 0.724750622468922*x*y*z + 0.630146164462021*x*y - 0.690968507486295*x*z**3 - 0.12442556368759*x*z**2 - 0.496555882470579*x*z + 0.0809323387905233*x + 0.0113323258784397*y**4 - 0.0247689443465589*y**3*z + 0.118782418819956*y**3 + 0.0687662568798315*y**2*z**2 - 0.247141943680235*y**2*z + 0.461628174707939*y**2 - 0.0442022398229372*y*z**3 + 0.415072461077496*y*z**2 - 0.633447995669211*y*z + 0.752767098639409*y + 0.186666432940899*z**4 - 0.296366521474835*z**3 + 0.565807556196118*z**2 + 0.168658949987404*z + 0.109139797748772


def main():
    global rR
    global v
    v = polynomial_surface2(rR)

    ax = implicit.plot_implicit(f4behind_filtered,bbox=(-5,5,-5,5,-5,5))

    #rR = tree_time.tree_time(rR,0.05)
    rR=array(rR)
    Axes3D.scatter(ax,rR[:,0],rR[:,1],rR[:,2])

    plt.show()
