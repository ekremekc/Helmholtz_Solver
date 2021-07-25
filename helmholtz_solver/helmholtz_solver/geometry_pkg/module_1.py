
"""
This script defines the geometry for Rijke Tube to pass it to the GMSH python API

There are two methods to define the geometry. This script defines the geometry by using 
dictionary for the bspline edge. The other method uses list for BSpline definition.
Numbering of edges are quite important, please see the instructions.

Actually, using dictionaries is much more convenient, because it tags the edge key which eases to 
accessing edge parameters while building mesh parameters in main_geometry.py

LAST UPDATE: USE DICTIONARIES!
"""


def g(p0, p1, n):
    # n is the number of intervals
    pts = []
    for i in range(n + 1):
        x = p0[0] + i*(p1[0] - p0[0])/n
        y = p0[1] + i*(p1[1] - p0[1])/n
        pts.append([x, y])
    return pts


def h(p0, p1, n):
    """
    the same as g, but with additional points in the first and last intervals
    """
    # n is the number of intervals
    pts = g(p0, p1, n)

    def lininterp():

        # (linearly) interpolate between pts[0] and pts[1]
        x = (pts[0][0] + pts[1][0])/2
        y = (pts[0][1] + pts[1][1])/2
        pts.insert(1, [x, y])

        # (linearly) interpolate between pts[-2] and pts[-1]
        x = (pts[-2][0] + pts[-1][0])/2
        y = (pts[-2][1] + pts[-1][1])/2
        pts.insert(-1, [x, y])

    lininterp()

    return pts


def f():
    """
    This function builds the geometry structure of the desired geometry
    for mesh with main_geometry.py file

    Returns
    -------
    pts : list
        Main points, used for generating control points.
    ll : list
        list of lines for GMSH API.
    physical_ll : list
        list of physical lines for GMSH API.
    ctrl_pts : dict
        control points for BSpline edge.
    lcar : float
        characteristic length for GMSH API

    """

    # elementary points

    p0 = [0., - .0235]
    p1 = [1., - .0235]
    p2 = [1., .0235]
    p3 = [0., .0235]

    pts  = [p0, p1, p2, p3]

    """
    
            POINTS    
    
         3__________ 2
         |          |   
         |          |   
         |__________|   
         0           1  
         
         ELEMENTARY LINES
    
           ____3______ 
          |          |  
         0|          | 1
          |__________|   
               2       
        
        PHYSICAL CURVES   
    
           ____4______ 
          |          |  
         1|          | 2
          |__________|   
               3     
         
         * Elementary lines should be numbered as counter-clockwise
    """
    # elementary lines
    
    l0 = [3, 0]  # inlet 
    l1 = [1, 2]  # outlet

    """ Order of the l0, l1, l2.. is not important
    """
    ll = [l0, l1]
    
    # physical_lines

    physical_ll = [[0], [1]]

    # control points
    """ For unified BSplines (from point 0 to point 2), the way of doing this is;
        1) Build control points from 0 to 1
        2) Build control points from 1 to 2 by EXCLUDING intersection point on point 1
        3) sum the lists of points
    make sure the intersection points are not doubled in the final list for this BSpline.
    see module_1_b
    """
    pts3 = h(p0, p1, 10)
    pts4 = h(p2, p3, 10)

    ctrl_pts = {3: pts3, 4: pts4}
    # ctrl_pts = {3: pts4, 4: pts3} # Also works

    lcar = 9.4e-3

    return pts, ll, physical_ll, ctrl_pts, lcar


if __name__ == '__main__':

    import numpy as np
    import matplotlib.pyplot as plt

    s = (pts, ll, physical_ll, ctrl_pts, lcar) = f()

    for t in s:
        print(t)
        print()
                                     
    ctrl_pts_3 = np.array([np.array(foo) for foo in ctrl_pts[3]])
    # ctrl_pts_3 = np.concatenate([np.array(foo) for foo in ctrl_pts[3]])

    ctrl_pts_4 = np.array([np.array(foo) for foo in ctrl_pts[4]])
    # ctrl_pts_4 = np.concatenate([np.array(foo) for foo in ctrl_pts[4]])

    plt.plot(ctrl_pts_3[:, 0], ctrl_pts_3[:, 1], 's', markersize=3)
    plt.plot(ctrl_pts_4[:, 0], ctrl_pts_4[:, 1], 's', markersize=3)
    plt.axis('equal')
    plt.show()
