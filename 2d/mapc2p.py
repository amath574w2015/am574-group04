

def mapc2p(xc,yc):
    """
    Specifies the mapping to curvilinear coordinates -- should be consistent
    with mapc2p.f
    """
    from numpy import where, cos, pi

    # Assumes 0 <= xc <= 32 and -1 <= yc <= 1

    # Nozzle geometry
    rad1 = 1-0.3*(1 + cos(pi*(xc-16)/8.))
    radius = where(abs(xc-16)<8, rad1, 1.)
    xp = xc
    yp = yc * radius
    return xp,yp

def plot_grid(mx=32,my=10):
    from pylab import linspace, meshgrid, plot, figure
    xc1 = linspace(0, 32, mx+1)
    yc1 = linspace(-1, 1, my+1)
    xc,yc = meshgrid(xc1,yc1)
    xp,yp = mapc2p(xc,yc)
    figure(figsize=(15,5))
    plot(xp,yp,'k')
    plot(xp.T,yp.T,'k')
    
