# def g(p0, p1, n):
#     # n is the number of intervals
#     pts = []
#     for i in range(n + 1):
#         x = p0[0] + i*(p1[0] - p0[0])/n
#         y = p0[1] + i*(p1[1] - p0[1])/n
#         pts.append([x, y])
#     return pts
#
#
# def h(p0, p1, n):
#     """
#     the same as g, but with additional points in the first and last intervals
#     """
#     # n is the number of intervals
#     pts = g(p0, p1, n)
#
#     def lininterp():
#
#         # (linearly) interpolate between pts[0] and pts[1]
#         x = (pts[0][0] + pts[1][0])/2
#         y = (pts[0][1] + pts[1][1])/2
#         pts.insert(1, [x, y])
#
#         # (linearly) interpolate between pts[-2] and pts[-1]
#         x = (pts[-2][0] + pts[-1][0])/2
#         y = (pts[-2][1] + pts[-1][1])/2
#         pts.insert(-1, [x, y])
#
#     lininterp()
#
#     return pts


def my_func(p0, p1, n, start_interp=1, end_interp=1, start=False, end=False):

    pts1 = []
    for i in range(n + 1):
        x = p0[0] + i * (p1[0] - p0[0]) / n
        y = p0[1] + i * (p1[1] - p0[1]) / n
        pts1.append([x, y])

    def start_lininterp():

        # (linearly) interpolate between pts[0] and pts[1]
        x1 = (pts1[0][0] + pts1[1][0])/2
        y1 = (pts1[0][1] + pts1[1][1])/2
        pts1.insert(1, [x1, y1])

    def end_lininterp():

        # (linearly) interpolate between pts[-2] and pts[-1]
        x1 = (pts1[-2][0] + pts1[-1][0])/2
        y1 = (pts1[-2][1] + pts1[-1][1])/2
        pts1.insert(-1, [x1, y1])

    for _ in range(start_interp):
        start_lininterp()

    for _ in range(end_interp):
        end_lininterp()

    if not start:
        pts1.pop(0)

    if not end:
        pts1.pop(-1)

    return pts1


def f():

    """
        _____________               ___________________
        |            \              |                  |
    r1  |             \_____________|                  | r3
        |                   ^                          | 
    _  _|_  _  _  _  _  _  _r2_  _  _  _  _  _  _  _  _|_  _  _ 
        |                                              |
        |              _____________                   |
        |             /             |                  |
        |____________/              |__________________|
             l1       l2    l3               l4
    """

    l_1 = 96.
    l_2 = 67.5
    l_3 = 78.
    l_4 = 400.

    r_1 = 32.5
    r_2 = 11.
    r_3 = 35.

    # nondimensionalization

    l_1 /= l_4
    l_2 /= l_4
    l_3 /= l_4

    r_1 /= l_4
    r_2 /= l_4
    r_3 /= l_4

    l_4 /= l_4

    # elementary points

    p0 = [0., - r_1]
    p1 = [l_1 + l_2 + l_3 + l_4, - r_3]

    p2 = [l_1 + l_2 + l_3 + l_4,   r_3]
    p3 = [0.,   r_1]

    p4 = [l_1, - r_1]
    p5 = [l_1 + l_2, - r_2]
    p6 = [l_1 + l_2 + l_3, - r_2]
    p7 = [l_1 + l_2 + l_3, - r_3]

    p8 = [l_1 + l_2 + l_3, r_3]
    p9 = [l_1 + l_2 + l_3, r_2]
    p10 = [l_1 + l_2, r_2]
    p11 = [l_1, r_1]

    pts = [p5, p6, p7, p1, p2, p8, p9, p10]
    # pts = [p0, p4, p5, p6, p7, p1, p2, p8, p9, p10, p11, p3, p0]

    # elementary lines

    l4 = [0, 1]
    l5 = [1, 2]
    l6 = [2, 3]
    l1 = [3, 4]
    l7 = [4, 5]
    l8 = [5, 6]
    l9 = [6, 7]

    ll = [l4, l5, l6, l1, l7, l8, l9]

    # physical_lines

    physical_ll = [[0, 1, 2], [3], [4, 5, 6]]

    # control points

    # to be improved

    pts10 = my_func(p10, p11, 4, start=True, end=True)
    pts11 = my_func(p11, p3, 6, end_interp=2)
    pts0 = my_func(p3, p0, 4, start_interp=2, end_interp=2)
    pts2 = my_func(p0, p4, 6, start_interp=2)
    pts3 = my_func(p4, p5, 4, start=True, end=True)

    ctrl_pts = pts10 + pts11 + pts0 + pts2 + pts3

    lcar = 5e-3

    return pts, ll, physical_ll, ctrl_pts, lcar


if __name__ == '__main__':

    import numpy as np
    import matplotlib.pyplot as plt

    from geometry_pkg import bspline_module

    BSpline = bspline_module.BSpline

    pts, ll, physical_ll, ctrl_pts, lcar = f()

    pts = np.array(pts)
    ctrl_pts = np.array(ctrl_pts)

    bspl = BSpline(ctrl_pts, lcar=lcar)
    bspl.create_curve()
    bspl_pts = bspl.pts

    marker_style = dict(color='k', linestyle='',
                        marker='o', fillstyle='none', markeredgewidth=0.5, markersize=3)

    plt.plot(pts[:, 0], pts[:, 1], 'dimgray', Linewidth=0.5)
    plt.plot(bspl_pts[:, 0], bspl_pts[:, 1], '--k', Linewidth=1)
    plt.plot(ctrl_pts[:, 0], ctrl_pts[:, 1], **marker_style)
    plt.axis('equal')
    plt.show()
