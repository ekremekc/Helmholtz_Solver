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
     3__________ 2
     |          |   if you want to obtain bspline from 0 to 1,
     |          |   start looping from point 1 like;
     |__________|   pts  = [p1, p2, p3, p0]
     0           1  
    """
    
    # elementary points

    p0 = [0., 0.]
    p1 = [.4, 0.]
    p2 = [.4, .1]
    p3 = [0., .1]

    pts  = [p1, p2, p3, p0]

    # elementary lines

    l0 = [0, 1]  # inlet
    l1 = [1, 2]  # outlet
    l2 = [2, 3]  # top

    ll = [l0, l1, l2]

    # physical_lines

    physical_ll = [[0], [1], [2]]

    # control points

    pts3 = g(p0, p1, 10)
    # pts4 = h(p2, p3, 10)

    # ctrl_pts = {3: pts3}
    ctrl_pts = pts3
    
    # perturb the control point
    ctrl_pts[4][1]+=0.01
    
    lcar = 1e-2
    return pts, ll, physical_ll, ctrl_pts, lcar

# def f():

#     # elementary points

#     p0 = [0., - .0235]
#     p1 = [1., - .0235]
#     p2 = [1., .0235]
#     p3 = [0., .0235]

#     pts  = [p0, p1, p2, p3]

#     # elementary lines

#     l0 = [3, 0]  # inlet
#     l1 = [2, 3]  # outlet

#     ll = [l0, l1]

#     # physical_lines

#     physical_ll = [[0], [1]]

#     # control points

#     pts3 = h(p0, p1, 10)
#     pts4 = h(p1, p2, 10)

#     ctrl_pts = {3: pts3, 4: pts4}

#     lcar = 9.4e-3

#     return pts, ll, physical_ll, ctrl_pts, lcar

if __name__ == '__main__':

    import numpy as np
    import matplotlib.pyplot as plt

    s = (pts, ll, physical_ll, ctrl_pts, lcar) = f()

    for t in s:
        print(t)
        print()

    ctrl_pts_3 = np.array([np.array(foo) for foo in ctrl_pts])
    
    # ctrl_pts_3 = np.array([np.array(foo) for foo in ctrl_pts[3]])
    # ctrl_pts_3 = np.concatenate([np.array(foo) for foo in ctrl_pts[3]])

    # ctrl_pts_4 = np.array([np.array(foo) for foo in ctrl_pts[4]])
    # ctrl_pts_4 = np.concatenate([np.array(foo) for foo in ctrl_pts[4]])

    plt.plot(ctrl_pts_3[:, 0], ctrl_pts_3[:, 1], 's', markersize=3)
    # plt.plot(ctrl_pts_4[:, 0], ctrl_pts_4[:, 1], 's', markersize=3)
    plt.axis('equal')
    plt.show()
    
    """
    
    
    
    
    """
