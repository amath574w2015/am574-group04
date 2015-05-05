u""" 
Set up the plot figures, axes, and items to be done for each frame.

This module is imported by the plotting routines and then the
function setplot is called to set the plot parameters.
    
"""
import numpy as np
import math
import matplotlib.pyplot as plt

def ArMach(M,x):
	g = 1.40
	g_m = g - 1.0
	g_a = g + 1.0
	radius = 0.750-0.250*math.sin(math.pi*(x-10)/10)
	area = radius**2.0*math.pi
	at = (0.750-0.250)**2.0*math.pi
	return area/at - (g_a/2.0)**(-g_a/(2.0*g_m))* (1.0 + g_m/2.0*M**2.0)**(g_a/(2.0*g_m))/M


def analtMach():
	from pylab import plot
	import scipy.optimize
	import numpy as np
	x = np.linspace(0.0,30.0,num=50)
	M = []
	for i in range(len(x)):
	   if(x[i] < 15.0): 
	       Mlow = 0.01
	       Mhigh = 1.0
	   else:
	       Mlow = 1.0
	       Mhigh = 4.0
	   if x[i] < 5.0:
	       xx = 5.0
	   elif x[i] > 25.0:
	       xx = 25.0
	   else:
	       xx = x[i]
	   M.append(scipy.optimize.brentq(ArMach, Mlow, Mhigh, args=(xx)))
	plt.hold(True)
	plot(x,M, 'r',linewidth='2')

def analtrho(current_data):
	from pylab import plot
	import scipy.optimize
	import numpy as np
	q = current_data.q
	x = np.linspace(0.0,30.0,num=50)
	M = []
	rho = []
	g = 1.4
	g_m = g - 1.0
	g_a = g + 1.0
	for i in range(len(x)):
	   if(x[i] < 15.0):
	       Mlow = 0.01
	       Mhigh = 1.0
	   else:
	       Mlow = 1.0
	       Mhigh = 4.0
	   if x[i] < 5.0:
	       xx = 5.0
	   elif x[i] > 25.0:
	       xx = 25.0
	   else:
	       xx = x[i]

	   M.append(scipy.optimize.brentq(ArMach, Mlow, Mhigh, args=(xx)))
	   rho.append(q[0,0]*(1.0+g_m/2.0*M[i]**2.0)**(-1/g_m))
	plot(x,rho, 'r',linewidth='2')


def analtp(current_data):
	from pylab import plot
	import scipy.optimize
	import numpy as np
	q = current_data.q
	x = np.linspace(0.0,30.0,num=50)
	M = []
	p = []
	g = 1.4
	g_m = g - 1.0
	g_a = g + 1.0 
	pt = (g_m) * (q[2,0] - 0.5*q[1,0]**2/ q[0,0])
	for i in range(len(x)):
	   if(x[i] < 15.0):
	       Mlow = 0.01
	       Mhigh = 1.0
	   else:
	       Mlow = 1.0
	       Mhigh = 4.0
	   if x[i] < 5.0:
	       xx = 5.0
	   elif x[i] > 25.0:
	       xx = 25.0
	   else: 
	       xx = x[i]
	   M.append(scipy.optimize.brentq(ArMach, Mlow, Mhigh, args=(xx)))
	   p.append(pt*(1.0+g_m/2.0*M[i]**2.0)**(-g/g_m))
	plot(x,p,'r',linewidth='2')


def analtvel(current_data):
	from pylab import plot
	import scipy.optimize
	import numpy as np
	q = current_data.q
	x = np.linspace(0.0,30.0,num=50)
	M = []
	p = []
	rho = []
	u = []
	g = 1.4
	g_m = g - 1.0
	g_a = g + 1.0
	pt = (g_m) * (q[2,0] - 0.5*q[1,0]**2/ q[0,0])
	for i in range(len(x)):
	   if(x[i] < 15.0):
	       Mlow = 0.01
	       Mhigh = 1.0
	   else:
	       Mlow = 1.0
	       Mhigh = 4.0
	   if x[i] < 5.0:
	       xx = 5.0
	   elif x[i] > 25.0:
	       xx = 25.0
	   else:
	       xx = x[i]
	   M.append(scipy.optimize.brentq(ArMach, Mlow, Mhigh, args=(xx)))
	   p.append(pt*(1.0+g_m/2.0*M[i]**2.0)**(-g/g_m))
	   rho.append(q[0,0]*(1.0+g_m/2.0*M[i]**2.0)**(-1/g_m))
	   u.append(M[i]*math.sqrt(g*p[i]/rho[i]))
	plot(x,u,'r',linewidth='2')

 

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
    plotitem.plot_var = 0  
    plotitem.plotstyle = '-o'
    plotitem.color = 'b'
    plotaxes.afteraxes = analtrho


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
    plotaxes.afteraxes = analtvel

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
    plotaxes.afteraxes = analtp

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
        xshock = 15.00
        plt.hold(True)
        plot([xshock, xshock], [0, 3], 'r')  
        plot([0, 30],[1, 1],'--r')	
 

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = Mach
    plotitem.plotstyle = '-o'
    plotitem.color = 'b'
    plotaxes.afteraxes = add_shock_location
#    plotaxes.afteraxes = analtMach

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

    
