import numpy as np
import dolfin as dolf
import matplotlib.pyplot as plt


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
        
# def save_result(p, saveformat = ".pvd"):
    