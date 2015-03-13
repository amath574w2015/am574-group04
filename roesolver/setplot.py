""" 
Set up the plot figures, axes, and items to be done for each frame.

This module is imported by the plotting routines and then the
function setplot is called to set the plot parameters.
    
""" 

#--------------------------
def setplot(plotdata):
#--------------------------
    import numpy as np
    import math
    import matplotlib.pyplot as plt
    """ 
    Specify what is to be plotted at each frame.
    Input:  plotdata, an instance of clawpack.visclaw.data.ClawPlotData.
    Output: a modified version of plotdata.
    
    """ 
    def density(current_data):
        import math

        q = current_data.q
        mx = len(q[1,:])
        xlower = 6.0
        dx = (26 - xlower)/mx
	radius = []
	area = []
        for i in range(mx):
		xcell = xlower + (i +1.0 - 0.5)*dx 
		if (abs(xcell-16.0) < 8.0):
			radius.append(1.0-0.30 * (1.0+math.cos(math.pi*(xcell-16.0)/8.0)))
		else:
			radius.append(1.0)
        
		area.append(math.pi*radius[i]**2.0)
	density = q[0,:]/area[:] 
        return density

    plotdata.clearfigures()  # clear any old figures,axes,items data

    # Figure for q[0]
    plotfigure = plotdata.new_plotfigure(name='Density, Momentum, and Energy', figno=1)
    plotfigure.kwargs = {'figsize':(10,12)}

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(5,1,1)'   # top figure
    plotaxes.xlimits = 'auto'
    plotaxes.ylimits = 'auto' #[0.0,5.0]
    plotaxes.title = 'Density'

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = 0 #density 
    plotitem.plotstyle = '-o'
    plotitem.color = 'b'


    # Figure for q[1]
    def vel(current_data):
        q = current_data.q
        u = q[1,:] / q[0,:]
        return u


    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(5,1,2)'   # bottom figure
    plotaxes.xlimits = 'auto'
    plotaxes.ylimits = 'auto'  #[-2.0,15.0]
    plotaxes.title = 'Velocity'

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = vel
    plotitem.plotstyle = '-o'
    plotitem.color = 'b'

    # Figure for q[2]

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(5,1,3)'   # 3rd figure
    plotaxes.xlimits = 'auto'
    plotaxes.ylimits = 'auto' #[0.5,35.0]
    plotaxes.title = 'Energy'

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = 2
    plotitem.plotstyle = '-o'
    plotitem.color = 'b'


    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(5,1,4)'   # bottom figure
    plotaxes.xlimits = 'auto'
    plotaxes.ylimits = 'auto' #[0.5,35.0]
    plotaxes.title = 'Pressure'

    def pressure(current_data):
        q = current_data.q
        gamma = 1.4
        p = (gamma-1.) * (q[2,:] - 0.5*q[1,:]**2/ q[0,:])
        return p

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = pressure
    plotitem.plotstyle = '-o'
    plotitem.color = 'b'

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(5,1,5)'   # bottom figure
    plotaxes.xlimits = 'auto'
    plotaxes.ylimits = 'auto' #[0.5,35.0]
    plotaxes.title = 'Mach'

    def Mach(current_data):
        q = current_data.q
        gamma = 1.4
        p = (gamma-1.) * (q[2,:] - 0.5*q[1,:]**2 / q[0,:])
        M = (q[1,:]/q[0,:]) / (gamma*p/q[0,:])**0.5 
        return M
    
    def add_shock_location(current_data):
        from pylab import plot
        xshock = 17.01
        plot([xshock, xshock], [0, 2], 'r')  
	

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = Mach
    plotitem.plotstyle = '-o'
    plotitem.color = 'b'
    plotaxes.afteraxes = add_shock_location
 
    # Parameters used only when creating html and/or latex hardcopy
    # e.g., via clawpack.visclaw.frametools.printframes:

    plotdata.printfigs = True                # print figures
    plotdata.print_format = 'png'            # file format
    plotdata.print_framenos = 'all'          # list of frames to print
    plotdata.print_fignos = 'all'            # list of figures to print
    plotdata.html = True                     # create html files of plots?
    plotdata.html_homelink = '../README.html'
    plotdata.latex = True                    # create latex file of plots?
    plotdata.latex_figsperline = 5          # layout of plots
    plotdata.latex_framesperline = 1         # layout of plots
    plotdata.latex_makepdf = False           # also run pdflatex?

    return plotdata

    
