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

def my_func(p0, p1, n, start=False, end=False):

    pts1 = g(p0, p1, n)
    

    if not start:
        pts1.pop(0)

    if not end:
        pts1.pop(-1)

    return pts1

def f():

    """
          3 _____ 2
           |    |   if you want to obtain bspline from 0 to 1,
     5_____|    |
     |     4    |   start looping from point 1 like;
     |__________|   pts  = [p1, p2, p3, p0]
     0           1  

     EDGES  

            __3__ 
            |   |   
      __5___|4  |2
    6|          |   
     |__________|   
          1        
    
    PHYSICAL CURVES   

             __2__ 
            |   |   
      __3___|4  |4
    4|          |   
     |__________|   
          1       

    """

    # elementary points

    p0 = [0., 0.]
    p1 = [2., 0.]
    p2 = [2., 2.]
    p3 = [1., 2.]
    p4 = [1., 1.]
    p5 = [0., 1.]

    pts  = [p0, p1, p2, p3, p4, p5] 

    # elementary lines

    l0 = [0, 1]  
    l1 = [2, 3]  
    l2 = [4, 5]
    
    

    ll = [ l0, l1, l2]
    
    # BSpline Lines
    bsp1 = [1,2]
    bsp2 = [3,4]
    bsp3 = [5,0]
    bspline_loops = {2: bsp1, 4: bsp2, 6:bsp3}
    
    # physical_lines

    physical_ll = [[0],[1],[2]]

    # control points
    
    pts2 = g(p1, p2, 10)
    pts4 = g(p3, p4, 10)
    pts6 = g(p5, p0, 10)
    
    ctrl_pts = {2:pts2 , 4:pts4, 6:pts6}
    
    # perturb the control point
    # ctrl_pts[4][1]+=0.01
    # ctrl_pts[12][0]-=0.01
    # ctrl_pts[20][1]-=0.01
    lcar = 1e-1

    return pts, ll, physical_ll, ctrl_pts, lcar

if __name__ == '__main__':

    import numpy as np
    import matplotlib.pyplot as plt

    s = (pts, ll, physical_ll, ctrl_pts, lcar, bspline_loops) = f()

    for t in s:
        print(t)
        print()

    ctrl_pts_3 = np.array([np.array(ctrl_pts[key]) for key in ctrl_pts])
    
    # ctrl_pts_3 = np.array([np.array(foo) for foo in ctrl_pts[3]])
    # ctrl_pts_3 = np.concatenate([np.array(foo) for foo in ctrl_pts[3]])

    # ctrl_pts_4 = np.array([np.array(foo) for foo in ctrl_pts[4]])
    # ctrl_pts_4 = np.concatenate([np.array(foo) for foo in ctrl_pts[4]])

    [plt.plot(np.array(ctrl_pts[key])[:, 0], np.array(ctrl_pts[key])[:, 1], 's', markersize=3) for key in ctrl_pts]
    # plt.plot(ctrl_pts_4[:, 0], ctrl_pts_4[:, 1], 's', markersize=3)
    # plt.axis('equal')
    # plt.show()
    # 
    """
    
    
    
    
    """
