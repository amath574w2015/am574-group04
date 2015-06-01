

def mapc2p(xc,yc):
    """
    Specifies the mapping to curvilinear coordinates -- should be consistent
    with mapc2p.f
    """
    from numpy import where, cos, pi
    import math

    # Assumes 0 <= xc <= 30 and -1 <= yc <= 1

    # Nozzle geometry
    rad1 = 0.75-0.25*cos(pi/2 - pi*(xc-10.)/10.)
    radius = where(abs(xc-15)<10, rad1, 1.)
    xp = xc
    yp = yc * radius
    return xp,yp

def plot_grid(mx=30,my=10):
    from pylab import linspace, meshgrid, plot, figure
    xc1 = linspace(3, 28, mx+1)
    yc1 = linspace(-1, 1, my+1)
    xc,yc = meshgrid(xc1,yc1)
    xp,yp = mapc2p(xc,yc)
    figure(figsize=(15,7))
    plot(xp,yp,'k')
    plot(xp.T,yp.T,'k')
    
