u""" 
Set up the plot figures, axes, and items to be done for each frame.

This module is imported by the plotting routines and then the
function setplot is called to set the plot parameters.
    
"""
import numpy as np
import math
import matplotlib.pyplot as plt

def ArMach(M,at,x):
	g = 1.40
	g_m = g - 1.0
	g_a = g + 1.0
	radius = 0.750-0.250*math.sin(math.pi*(x-10)/10)
	area = radius**2.0*math.pi
	#at = (0.750-0.250)**2.0*math.pi
	return area/at - (g_a/2.0)**(-g_a/(2.0*g_m))* (1.0 + g_m/2.0*M**2.0)**(g_a/(2.0*g_m))/M


def analtMach(current_data):
	from pylab import plot
	import scipy.optimize
	import numpy as np

#	Adjust this tolerance to get an accurate shock location
	tol = 5*10**-5.0

	q = current_data.q
	x = np.linspace(0.0,30.0,num=100)
	g = 1.4
	g_m = g - 1.0
	g_a = g + 1.0 
	rhot = 1.0
	pt =  g_m * (q[2,0] - 0.5*q[1,0]**2/ q[0,0])
	patm = 4.0*0.31 #g_m * (q[2,len(q[2])-1] - 0.5*q[1,len(q[2])-1]**2/ q[0,len(q[2])-1]) 
	at = (0.750-0.250)**2.0*math.pi
	M = []
	p = []
	rho = []
	u = []

#	Open files to write
	Mwrite = open('M_analt','wt')
	pwrite = open('p_analt','wt')
	rhowrite = open('rho_analt','wt')
	uwrite = open('u_analt','wt')
	xwrite = open('x_analt','wt')

	Msup = scipy.optimize.brentq(ArMach, 1.0, 4.0, args=(at, 25))
	psup = pt*(1.0 + g_m/2.0*Msup**2.0)**(-g/g_m)
	Msub = scipy.optimize.brentq(ArMach, 0.01, 1.0, args=(at,5))
	psub = pt*(1.0 + g_m/2.0*Msub**2.0)**(-g/g_m)

#	Flow with a shock in nozzle
	if(patm > psup and patm < psub):
	   xs_guess = 20.0

 	   for j in range(50000):
		Ms = scipy.optimize.brentq(ArMach, 1.0, 4.0, args = (at,xs_guess))
		As = (0.750-0.250*math.sin(math.pi*(xs_guess-10)/10))**2.0*math.pi
		pexitt_pt = ( (g_a *Ms**2.0)/( g_m*Ms**2.0 + 2.0))**(g/g_m)*(g_a/(2.0*g*Ms**2.0-g_m))**(1.0/g_m)
		M2s = math.sqrt( (2.0 + g_m*Ms**2.0) / (2.0*g*Ms**2.0 - g_m) )
		Astar_2 = As * math.sqrt(((2.0/g_a)*(1.0 + g_m/2.0*M2s**2.0))**(-g_a/g_m)) *M2s

		Me = scipy.optimize.brentq(ArMach, 0.01, 1.0, args = (Astar_2, 25))
		patm_pexitt = (1.0 + g_m/2.0*Me**2.0)**(-g/g_m)

#	        Solve for guessed patm by patm = (patm/pt2) * (pt2/pt) * pt
		patm_guess = patm_pexitt * pexitt_pt * pt

		err = abs(patm_guess - patm)
		if err <= tol:
		    xs = xs_guess
		    break
		elif patm_guess > patm:
		    xs_guess = xs_guess + 0.0001
		elif patm_guess < patm:
		    xs_guess = xs_guess - 0.0001

#	   Shock can't be found so assume pure isentropic
	   if j == 49999:
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
			M.append(scipy.optimize.brentq(ArMach, Mlow, Mhigh, args=(at, xx)))
			p.append(pt*(1.0+g_m/2.0*M[i]**2.0)**(-g/g_m))
			rho.append(rhot*(1.0+g_m/2.0*M[i]**2.0)**(-1.0/g_m))
			u.append(M[i]*math.sqrt(g*p[i]/rho[i]))

	#		Write out data to file
			Mwrite.write(str(M[i]))
			Mwrite.write('\n')
			pwrite.write(str(p[i]))
			pwrite.write('\n')
			rhowrite.write(str(rho[i]))
			rhowrite.write('\n')
			uwrite.write(str(u[i]))
			uwrite.write('\n')
			xwrite.write(str(x[i]))
			xwrite.write('\n')

		    plt.hold(True)
		    plot(x,M, 'r',linewidth='2')


	   else:
		   x1 = np.linspace(0.0,xs,num=50)
		   x2 = np.linspace(xs,30,num=50)
		   M1 = []
		   M2 = []
		   p1 = []
		   p2 = []
		   u1 = []
	  	   u2 = []
		   rho1 = []
		   rho2 = []
		   pexitt = pexitt_pt * pt

#		   Determine the total density after shock, rhot2
		   rhos = rhot * (1.0 + g_m/2.0*Ms**2.0)**(-1.0/g_m)
		   rho2s = rhos * (g_a * Ms**2.0) / (g_m*Ms**2.0 + 2.0)
		   rhot2 = rho2s *  (1.0 + g_m/2.0*M2s**2.0)**(1.0/g_m)

#  	           Create Mach structure before the shock
		   for i in range(len(x1)):
			if(x1[i] < 15.0): 
			   Mlow = 0.01
			   Mhigh = 1.0
			else:
			   Mlow = 1.0
			   Mhigh = 4.0

			if x1[i] < 5.0:
			   xx = 5.0
			else:
			   xx = x1[i]
			M1.append(scipy.optimize.brentq(ArMach, Mlow, Mhigh, args=(at, xx)))
			p1.append(pt*(1.0+g_m/2.0*M1[i]**2.0)**(-g/g_m))
			rho1.append(rhot*(1.0+g_m/2.0*M1[i]**2.0)**(-1.0/g_m))
			u1.append(M1[i]*math.sqrt(g*p1[i]/rho1[i]))

#			Write out data to file
			Mwrite.write(str(M1[i]))
			Mwrite.write('\n')
			pwrite.write(str(p1[i]))
			pwrite.write('\n')
			rhowrite.write(str(rho1[i]))
			rhowrite.write('\n')
			uwrite.write(str(u1[i]))
			uwrite.write('\n')
			xwrite.write(str(x1[i]))
			xwrite.write('\n')

		   M1.append(scipy.optimize.brentq(ArMach, 0.01, 1.0, args=(Astar_2, x2[0])))
		   x1=list(x1)
		   x1.append(x1[i])

#  	           Create Mach structure after the shock
		   for i in range(len(x2)):
			Mlow = 0.01
			Mhigh = 1.0
	
			if x2[i] > 25.0:
			   xx = 25.0
			else:
			   xx = x2[i]
			rho
			M2.append(scipy.optimize.brentq(ArMach, Mlow, Mhigh, args=(Astar_2, xx)))
			p2.append(pexitt*(1.0+g_m/2.0*M2[i]**2.0)**(-g/g_m))
			rho2.append(rhot2*(1.0+g_m/2.0*M2[i]**2.0)**(-1.0/g_m))
			u2.append(M2[i]*math.sqrt(g*p2[i]/rho2[i]))

#			Write out data to file
			Mwrite.write(str(M2[i]))
			Mwrite.write('\n')
			pwrite.write(str(p2[i]))
			pwrite.write('\n')
			rhowrite.write(str(rho2[i]))
			rhowrite.write('\n')
			uwrite.write(str(u2[i]))
			uwrite.write('\n')
			xwrite.write(str(x2[i]))
			xwrite.write('\n')

		 
		   plt.hold(True)
		   plot(x1,M1, 'r',linewidth='2')
		   plot(x2,M2,'r',linewidth='2')
#	Pure isentropic flow with subsonic outlet
	elif( patm > psub):
	   for i in range(len(x)):
		Mlow = 0.01
		Mhigh = 1.0
		if x[i] < 5.0:
		   xx = 5.0
		elif x[i] > 25.0:
		   xx = 25.0
		else:
		   xx = x[i]
		M.append(scipy.optimize.brentq(ArMach, Mlow, Mhigh, args=(at, xx)))
		p.append(pt*(1.0+g_m/2.0*M[i]**2.0)**(-g/g_m))
		rho.append(rhot*(1.0+g_m/2.0*M[i]**2.0)**(-1.0/g_m))
		u.append(M[i]*math.sqrt(g*p[i]/rho[i]))

#		Write out data to file
		Mwrite.write(str(M[i]))
		Mwrite.write('\n')
		pwrite.write(str(p[i]))
		pwrite.write('\n')
		rhowrite.write(str(rho[i]))
		rhowrite.write('\n')
		uwrite.write(str(u[i]))
		uwrite.write('\n')
		xwrite.write(str(x[i]))
		xwrite.write('\n')

	   plt.hold(True)
           plot(x,M, 'r',linewidth='2.5')	
#	Pure isentropic flow with supersonic outlet
	else:
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
		M.append(scipy.optimize.brentq(ArMach, Mlow, Mhigh, args=(at, xx)))
		p.append(pt*(1.0+g_m/2.0*M[i]**2.0)**(-g/g_m))
		rho.append(rhot*(1.0+g_m/2.0*M[i]**2.0)**(-1.0/g_m))
		u.append(M[i]*math.sqrt(g*p[i]/rho[i]))

#		Write out data to file
		Mwrite.write(str(M[i]))
		Mwrite.write('\n')
		pwrite.write(str(p[i]))
		pwrite.write('\n')
		rhowrite.write(str(rho[i]))
		rhowrite.write('\n')
		uwrite.write(str(u[i]))
		uwrite.write('\n')
		xwrite.write(str(x[i]))
		xwrite.write('\n')

	   plt.hold(True)
           plot(x,M, 'r',linewidth='2.5')

def analtrho(current_data):
	from pylab import plot
	import scipy.optimize
	import numpy as np
#	q = current_data.q
#	x = np.linspace(0.0,30.0,num=50)
#	M = []
#	rho = []
#	g = 1.4
#	g_m = g - 1.0
#	g_a = g + 1.0
#	at = (0.750-0.250)**2.0*math.pi
#	for i in range(len(x)):
#	   if(x[i] < 15.0):
#	       Mlow = 0.01
#	       Mhigh = 1.0
#	   else:
#	       Mlow = 1.0
#	       Mhigh = 4.0
#	   if x[i] < 5.0:
#	       xx = 5.0
#	   elif x[i] > 25.0:
#	       xx = 25.0
#	   else:
#	       xx = x[i]

#	   M.append(scipy.optimize.brentq(ArMach, Mlow, Mhigh, args=(at, xx)))
#	   rho.append(q[0,0]*(1.0+g_m/2.0*M[i]**2.0)**(-1/g_m))
	rhoopen = open('rho_analt','rt')
	xopen = open('x_analt','rt')
	rho = []
	x = []
	for i in range(100):
	   x.append(float(xopen.readline()))
	   rho.append(float(rhoopen.readline()))

	plot(x,rho, 'r',linewidth='2.5')


def analtp(current_data):
	from pylab import plot
	import scipy.optimize
	import numpy as np
#	q = current_data.q
#	x = np.linspace(0.0,30.0,num=50)
#	at = (0.750-0.250)**2.0*math.pi
#	M = []
#	p = []
#	g = 1.4
#	g_m = g - 1.0
#	g_a = g + 1.0 
#	pt = (g_m) * (q[2,0] - 0.5*q[1,0]**2/ q[0,0])
#	for i in range(len(x)):
#	   if(x[i] < 15.0):
#	       Mlow = 0.01
#	       Mhigh = 1.0
#	   else:
#	       Mlow = 1.0
#	       Mhigh = 4.0
#	   if x[i] < 5.0:
#	       xx = 5.0
#	   elif x[i] > 25.0:
#	       xx = 25.0
#	   else: 
#	       xx = x[i]
#	   M.append(scipy.optimize.brentq(ArMach, Mlow, Mhigh, args=(at,xx)))
#	   p.append(pt*(1.0+g_m/2.0*M[i]**2.0)**(-g/g_m))
	popen = open('p_analt','rt')
	xopen = open('x_analt','rt')
	p = []
	x = []
	for i in range(100):
	   x.append(float(xopen.readline()))
	   p.append(float(popen.readline()))

	plot(x,p,'r',linewidth='2.5')


def analtvel(current_data):
	from pylab import plot
	import scipy.optimize
	import numpy as np
#	q = current_data.q
#	at = (0.750-0.250)**2.0*math.pi
#	x = np.linspace(0.0,30.0,num=50)
#	M = []
#	p = []
#	rho = []
#	u = []
#	g = 1.4
#	g_m = g - 1.0
#	g_a = g + 1.0
#	pt = (g_m) * (q[2,0] - 0.5*q[1,0]**2/ q[0,0])
#	for i in range(len(x)):
#	   if(x[i] < 15.0):
#	       Mlow = 0.01
#	       Mhigh = 1.0
#	   else:
#	       Mlow = 1.0
#	       Mhigh = 4.0
#	   if x[i] < 5.0:
#	       xx = 5.0
#	   elif x[i] > 25.0:
#	       xx = 25.0
#	   else:
#	       xx = x[i]
#	   M.append(scipy.optimize.brentq(ArMach, Mlow, Mhigh, args=(at,xx)))
#	   p.append(pt*(1.0+g_m/2.0*M[i]**2.0)**(-g/g_m))
#	   rho.append(q[0,0]*(1.0+g_m/2.0*M[i]**2.0)**(-1/g_m))
#	   u.append(M[i]*math.sqrt(g*p[i]/rho[i]))
	uopen = open('u_analt','rt')
	xopen = open('x_analt','rt')
	u = []
	x = []
	for i in range(100):
	   x.append(float(xopen.readline()))
	   u.append(float(uopen.readline()))

	plot(x,u,'r',linewidth='2.5')

 

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
    plotfigure = plotdata.new_plotfigure(name='Mach, Density, Velocity, and Energy', figno=1)
    plotfigure.kwargs = {'figsize':(14,16)}

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(5,1,1)'   # bottom figure
    plotaxes.xlimits = 'auto'
    plotaxes.ylimits = 'auto' #[0.5,35.0]
    plotaxes.title = 'Mach'

    def Mach(current_data):
        q = current_data.q
        gamma = 1.4
        p = (gamma-1.) * (q[2,:] - 0.5*q[1,:]**2 / q[0,:])
        M = (q[1,:]/q[0,:]) / (gamma*p/q[0,:])**0.5 
        return M
     

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = Mach
    plotitem.plotstyle = '-o'
    plotitem.color = 'b'
    plotaxes.afteraxes = analtMach
    
    # Figure for q[0]

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(5,1,2)'   # top figure
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
    plotaxes.axescmd = 'subplot(5,1,3)'   # bottom figure
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
    plotaxes.axescmd = 'subplot(5,1,4)'   # 3rd figure
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
    plotaxes.axescmd = 'subplot(5,1,5)'   # bottom figure
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

    
