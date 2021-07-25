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

    # elementary points

    p0 = [0., 0.]
    p1 = [.4, 0.]
    p2 = [.4, .1]
    p3 = [0., .1]

    pts  = [p0, p1, p2, p3]

    """
    
    BOTTOM AND RIGHT EDGES OF THE TUBE PERTURBED
    
    
            POINTS    
    
         3__________ 2
         |          |   
         |          |   
         |__________|   
         0           1  
         
         ELEMENTARY LINES
    
           ____1______ 
          |          |  
         0|          | 2
          |__________|   
               2       
        
        PHYSICAL CURVES   
    
           ____2______ 
          |          |  
         1|          | 3
          |__________|   
               3     
         
         * Elementary lines should be numbered as counter-clockwise
    """
    # elementary lines
    
    l0 = [3, 0]  # inlet 
    l1 = [2, 3]  # outlet


    ll = [l0, l1]
    
    # physical_lines

    physical_ll = [[0], [1]]

    # control points

    pts3_bottom = g(p0, p1, 10)
    pts3_right  = my_func(p1, p2, 5, end = True)
    pts3 = pts3_bottom + pts3_right
    ctrl_pts = {3: pts3}
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

    # ctrl_pts_4 = np.array([np.array(foo) for foo in ctrl_pts[4]])
    # ctrl_pts_4 = np.concatenate([np.array(foo) for foo in ctrl_pts[4]])

    plt.plot(ctrl_pts_3[:, 0], ctrl_pts_3[:, 1], 's', markersize=3)
    # plt.plot(ctrl_pts_4[:, 0], ctrl_pts_4[:, 1], 's', markersize=3)
    plt.axis('equal')
    plt.show()
