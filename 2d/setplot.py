""" 
Set up the plot figures, axes, and items to be done for each frame.

This module is imported by the plotting routines and then the
function setplot is called to set the plot parameters.
    
""" 

from mapc2p import mapc2p
import numpy as np

#--------------------------
def setplot(plotdata):
#--------------------------
    
    """ 
    Specify what is to be plotted at each frame.
    Input:  plotdata, an instance of pyclaw.plotters.data.ClawPlotData.
    Output: a modified version of plotdata.
    
    """ 


    from clawpack.visclaw import colormaps

    plotdata.clearfigures()  # clear any old figures,axes,items data
    

    # Figure for pcolor plot
    plotfigure = plotdata.new_plotfigure(name='Density and Velocities', figno=0)
    plotfigure.kwargs = {'figsize':(14,20)}


    # Figure for q[0]
    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(3,1,1)'
    plotaxes.xlimits = 'auto'
    plotaxes.ylimits = 'auto'
    plotaxes.title = 'q[0]'
    plotaxes.scaled = True

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = 0
    plotitem.pcolor_cmap = colormaps.red_yellow_blue
    plotitem.pcolor_cmin = 0.
    plotitem.pcolor_cmax = 1.
    plotitem.add_colorbar = True
    plotitem.celledges_show = 0
    plotitem.patchedges_show = 1
    plotitem.MappedGrid = True
    plotitem.mapc2p = mapc2p
    plotitem.show = True       # show on plot?

    def uvel(current_data):
        q = current_data.q
        u = q[1,:,:] / q[0,:,:]
        return u

    def vvel(current_data):
        q = current_data.q
        v = q[2,:,:] / q[0,:,:]
        return v


    # FIgure for q[1]/q[0]
    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(3,1,2)'
    plotaxes.xlimits = 'auto'
    plotaxes.ylimits = 'auto'
    plotaxes.title = 'u-vel'
    plotaxes.scaled = True

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = uvel
    plotitem.pcolor_cmap = colormaps.red_yellow_blue
    plotitem.pcolor_cmin = 0.
    plotitem.pcolor_cmax = 4.
    plotitem.add_colorbar = True
    plotitem.celledges_show = 0
    plotitem.patchedges_show = 1
    plotitem.MappedGrid = True
    plotitem.mapc2p = mapc2p
    plotitem.show = True       # show on plot?

    # FIgure for q[2]/q[0]
    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(3,1,3)'
    plotaxes.xlimits = 'auto'
    plotaxes.ylimits = 'auto'
    plotaxes.title = 'v-vel'
    plotaxes.scaled = True

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = vvel
    plotitem.pcolor_cmap = colormaps.red_yellow_blue
    plotitem.pcolor_cmin = 0.
    plotitem.pcolor_cmax = 0.3
    plotitem.add_colorbar = True
    plotitem.celledges_show = 0
    plotitem.patchedges_show = 1
    plotitem.MappedGrid = True
    plotitem.mapc2p = mapc2p
    plotitem.show = True       # show on plot?

#   Figure for pcolor plot 2
    plotfigure = plotdata.new_plotfigure(name='q[3] and pressure', figno=1)
    plotfigure.kwargs = {'figsize':(14,20)}


    # Figure for q[3]
    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(2,1,1)'
    plotaxes.xlimits = 'auto'
    plotaxes.ylimits = 'auto'
    plotaxes.title = 'q[3]'
    plotaxes.scaled = True

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = 3
    plotitem.pcolor_cmap = colormaps.red_yellow_blue
    plotitem.pcolor_cmin = 0.
    plotitem.pcolor_cmax = 10.
    plotitem.add_colorbar = True
    plotitem.celledges_show = 0
    plotitem.patchedges_show = 1
    plotitem.MappedGrid = True
    plotitem.mapc2p = mapc2p
    plotitem.show = True       # show on plot?

    def pressure(current_data):
        q = current_data.q
        gamma = 1.4
        p = (gamma-1.) * (q[3,:,:] - 0.5*(q[1,:,:]**2+q[2,:,:]**2)/ q[0,:,:])
        return p

    # Figure for pressure
    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(2,1,2)'
    plotaxes.xlimits = 'auto'
    plotaxes.ylimits = 'auto'
    plotaxes.title = 'Pressure'
    plotaxes.scaled = True

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = pressure
    plotitem.pcolor_cmap = colormaps.red_yellow_blue
    plotitem.pcolor_cmin = 0.
    plotitem.pcolor_cmax = 4.
    plotitem.add_colorbar = True
    plotitem.celledges_show = 0
    plotitem.patchedges_show = 1
    plotitem.MappedGrid = True
    plotitem.mapc2p = mapc2p
    plotitem.show = True       # show on plot?

    def Mach(current_data):
        q = current_data.q
        gamma = 1.4
        p = (gamma-1.) * (q[3,:,:] - 0.5*(q[1,:,:]**2+q[2,:,:]**2)/q[0,:,:])
        M = (q[1,:,:]/q[0,:,:]) / (gamma*p/q[0,:,:])**0.5
        return M


    
    # Figure for contour plot
#    plotfigure = plotdata.new_plotfigure(name='contour', figno=2)

    # Set up for axes in this figure:
#    plotaxes = plotfigure.new_plotaxes()
#    plotaxes.axescmd = 'subplot(3,1,3)'
#    plotaxes.xlimits = 'auto'
#    plotaxes.ylimits = 'auto'
#    plotaxes.title = 'q[0]'
#    plotaxes.scaled = True

    # Set up for item on these axes:
#    plotitem = plotaxes.new_plotitem(plot_type='2d_contour')
#    plotitem.plot_var = 0
#    plotitem.contour_levels = np.linspace(-0.9, 0.9, 10)
#    plotitem.contour_colors = 'k'
#    plotitem.patchedges_show = 1
#    plotitem.MappedGrid = True
#    plotitem.mapc2p = mapc2p
#    plotitem.show = True       # show on plot?
    

    # Parameters used only when creating html and/or latex hardcopy
    # e.g., via pyclaw.plotters.frametools.printframes:

    plotdata.printfigs = True                # print figures
    plotdata.print_format = 'png'            # file format
    plotdata.print_framenos = 'all'          # list of frames to print
    plotdata.print_fignos = 'all'            # list of figures to print
    plotdata.html = True                     # create html files of plots?
    plotdata.html_homelink = '../README.html'   # pointer for top of index
    plotdata.latex = True                    # create latex file of plots?
    plotdata.latex_figsperline = 4           # layout of plots
    plotdata.latex_framesperline = 4         # layout of plots
    plotdata.latex_makepdf = False           # also run pdflatex?

    return plotdata

    
