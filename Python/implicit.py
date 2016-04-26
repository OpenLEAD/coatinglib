from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
import numpy as np
import time


def plot_implicit(fn, bbox=(-2.5,2.5,-2.5,2.5,-2.5,2.5),rc=100,ns=15,colors='b',ax=[]):
    ''' create a plot of an implicit function
    fn  ...implicit function (plot where fn==0)
    bbox ..the x,y,and z limits of plotted interval'''
    xmin, xmax, ymin, ymax, zmin, zmax = bbox
    if not ax:
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
    A = np.linspace(xmin, xmax, rc) # resolution of the contour
    B = np.linspace(xmin, xmax, ns) # number of slices
    A1,A2 = np.meshgrid(A,A) # grid on which the contour is plotted

    i=0
    t=time.time()
    initialTime = t
    N=len(B)
    print len(A1)
    for z in B: # plot contours in the XY plane
        X,Y = A1,A2
        Z = fn(X,Y,z)
        cset = ax.contour(X, Y, Z+z, [z], zdir='z',colors=colors)
        # [z] defines the only level to plot for this contour for this value of z
        i+=1
        if i%1==0:
            print str(i)+'/'+str(N)
            print 'elapsed time= '+str(time.time()-t)
            t=time.time()        
            print 'total Time = '+str(time.time()- initialTime)
    i=0
    t=time.time()
    initialTime = t
    for y in B: # plot contours in the XZ plane
        X,Z = A1,A2
        Y = fn(X,y,Z)
        cset = ax.contour(X, Y+y, Z, [y], zdir='y',colors=colors)

        i+=1
        if i%1==0:
            print str(i)+'/'+str(N)
            print 'elapsed time= '+str(time.time()-t)
            t=time.time()        
            print 'total Time = '+str(time.time()- initialTime)

    i=0
    t=time.time()
    initialTime = t
    for x in B: # plot contours in the YZ plane
        Y,Z = A1,A2
        X = fn(x,Y,Z)
        cset = ax.contour(X+x, Y, Z, [x], zdir='x',colors=colors)

        i+=1
        if i%1==0:
            print str(i)+'/'+str(N)
            print 'elapsed time= '+str(time.time()-t)
            t=time.time()        
            print 'total Time = '+str(time.time()- initialTime)

    # must set plot limits because the contour will likely extend
    # way beyond the displayed level.  Otherwise matplotlib extends the plot limits
    # to encompass all values in the contour.
    ax.set_zlim3d(zmin,zmax)
    ax.set_xlim3d(xmin,xmax)
    ax.set_ylim3d(ymin,ymax)

    #plt.show()
    return ax

def plot_implicit_save(fn, bbox=(-2.5,2.5,-2.5,2.5,-2.5,2.5),rc=100,ns=15,colors='b',ax=[]):
    xmin, xmax, ymin, ymax, zmin, zmax = bbox
    
    if not ax:
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

    A = np.linspace(xmin, xmax, rc) # resolution of the contour
    B = np.linspace(xmin, xmax, ns) # number of slices
    A1,A2 = np.meshgrid(A,A) # grid on which the contour is plotted

    i=0
    t=time.time()
    initialTime = t
    N=len(B)
    for z in B: # plot contours in the XY plane
        X,Y = A1,A2
        Z = fn(X,Y,z)
        np.savez_compressed('implicit/Z/'+str(i)+'.npz', array=Z)
        i+=1

    i=0
    for y in B: # plot contours in the XZ plane
        X,Z = A1,A2
        Y = fn(X,y,Z)
        np.savez_compressed('implicit/Y/'+str(i)+'.npz', array=Y)
        i+=1

    i=0
    for x in B: # plot contours in the YZ plane
        Y,Z = A1,A2
        X = fn(x,Y,Z)
        np.savez_compressed('implicit/X/'+str(i)+'.npz', array=X)
        i+=1
    return ax

def plot_implicit_load(folder,bbox=(-2.5,2.5,-2.5,2.5,-2.5,2.5),rc=100,ns=15,colors='b',ax=[]):
    xmin, xmax, ymin, ymax, zmin, zmax = bbox
    if not ax:
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
    A = np.linspace(xmin, xmax, rc) # resolution of the contour
    B = np.linspace(xmin, xmax, ns) # number of slices
    A1,A2 = np.meshgrid(A,A) # grid on which the contour is plotted

    i=0
    for z in B: # plot contours in the XY plane
        X,Y = A1,A2
        Z = np.load(folder+'Z/'+str(i)+'.npz')
        Z = Z['array']
        cset = ax.contour(X, Y, Z+z, [z], zdir='z',colors=colors)
        i+=1
    i=0

    for y in B: # plot contours in the XZ plane
        X,Z = A1,A2
        Y = np.load(folder+'Y/'+str(i)+'.npz')
        Y = Y['array']
        cset = ax.contour(X, Y+y, Z, [y], zdir='y',colors=colors)

        i+=1

    i=0
    for x in B: # plot contours in the YZ plane
        Y,Z = A1,A2
        X = np.load(folder+'X/'+str(i)+'.npz')
        X = X['array']
        cset = ax.contour(X+x, Y, Z, [x], zdir='x',colors=colors)
        i+=1

    ax.set_zlim3d(zmin,zmax)
    ax.set_xlim3d(xmin,xmax)
    ax.set_ylim3d(ymin,ymax)

    return ax
