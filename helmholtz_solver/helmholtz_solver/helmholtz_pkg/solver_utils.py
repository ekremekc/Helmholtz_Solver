import numpy as np
import dolfin as dolf
import matplotlib.pyplot as plt
import matplotlib.tri as tri
from mpl_toolkits.axes_grid1 import make_axes_locatable

def absolute_eigenfunction(p):
    """
    Calculates the absolute value of the given eigenfunction

    Parameters
    ----------
    p : dolfin.dolfin.function.function.Function
        Calculated unsplitted eigenfunction

    Returns
    -------
    p_abs : dolfin.dolfin.function.function.Function
        Calculated absolute eigenfunction

    """
    
    p_r,p_i = p.split(True)
    p_abs = p_r.copy()
    p_abs.vector()[:] = np.abs(p_r.vector()[:]+1j*p_i.vector()[:])
    return p_abs

def phase_eigenfunction(p):
    """
    Calculates the phase value of the given eigenfunction

    Parameters
    ----------
    p : dolfin.dolfin.function.function.Function
        Calculated unsplitted eigenfunction

    Returns
    -------
    p_abs : dolfin.dolfin.function.function.Function
        Calculated phase eigenfunction

    """
    
    p_r,p_i = p.split(True)
    p_abs = p_r.copy()
    p_abs.vector()[:] = np.angle(p_r.vector()[:]+1j*p_i.vector()[:])
    return p_abs

def barplot(p, save = False, plotname = "p",  plotformat = ".pdf"):
    """
    Plot the eigenfunction with scalarbar

    Parameters
    ----------
    p : dolfin.dolfin.function.function.Function
        Eigenfunction solution
    save : Boolean, optional
        Saving option. The default is False.
    plotname : String, optional
        name of the plot. The default is False.
    plotformat : String, optional
        plotting format option. The default is ".pdf".

    """
    
    fig = dolf.plot(p)
    plt.colorbar(fig)
    
    if save:
        
        filename = plotname + plotformat
        plt.savefig(filename)
        
def plot_2d(p, mesh, save = False, plotname = "2D",  plotformat = ".pdf", colorbar = True,
            orientation = 'V'):
    """
    Generate 2D nice plots

    Parameters
    ----------
    p : Variable we want to plot
        dolfin.function.function.Function.
    mesh : dolfin.cpp.mesh.Mesh
        mesh of the domain.
    save : Boolean, optional
        Saving option. The default is False.
    plotname : String, optional
        name of the plot. The default is False.
    plotformat : String, optional
        plotting format option. The default is ".pdf".

    Returns
    -------
    TYPE
        DESCRIPTION.

    """
    fig, ax = plt.subplots()
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    def mesh2triang(mesh):
        xy = mesh.coordinates()
        return tri.Triangulation(xy[:, 0], xy[:, 1], mesh.cells())
    
    # plt.gca().set_aspect('equal')
    if isinstance(p, dolf.Function):
        mesh = p.function_space().mesh()
        if (mesh.geometry().dim() != 2):
            raise(AttributeError)
        if p.vector().size() == mesh.num_cells():
            C = p.vector().array()
            tpc = ax.tripcolor(mesh2triang(mesh), C)
        else:
            C = p.compute_vertex_values(mesh)
            tpc = ax.tripcolor(mesh2triang(mesh), C, shading='gouraud')
            
    if colorbar == True:
        if orientation == "H":
            divider = make_axes_locatable(ax)
            cax = divider.new_vertical(size="5%", pad=0.7, pack_start=True)
            fig.add_axes(cax)
            fig.colorbar(tpc, cax=cax, orientation="horizontal")
        elif orientation == "V":
            fig.colorbar(tpc)
    elif isinstance(p, dolf.Mesh):
        if (p.geometry().dim() != 2):
            raise(AttributeError)
        plt.triplot(mesh2triang(p), color='k')
    
    if save:
        
        filename = plotname + plotformat
        plt.savefig(filename)
    