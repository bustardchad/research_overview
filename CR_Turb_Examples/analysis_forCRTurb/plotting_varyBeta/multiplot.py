import yt
from yt.units import dimensions
from yt.units import erg
from yt.units.yt_array import YTQuantity
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid

yt.enable_parallelism()

edenstocgs = 6.54e-11
denstocgs = 6.85e-27
gammatocgs = edenstocgs*denstocgs/1.67e-24 # ecr * rho/1.67e-24
prestocgs = edenstocgs
temptocgs = prestocgs/denstocgs


Myr = 1.
kpc = 1.

def crscale(field,data):
          return np.abs(data['user_out_var8'])/(0.667*3.0856e21)

def _pcr(field, data):
 return ds.arr(data['Ec'].v*edenstocgs*0.3333, 'g/(cm*s**2)')

def _densratio(field, data):
  return data['density']/(YTQuantity(1,"g/cm**3")*0.03066)

times = range(7,8,1)

for i in times:
 fn = "../"
 base = 'cr.out1.'  # beta = 1
 ds = yt.load(fn+base+str(i).zfill(5)+'.athdf')

 fn2 = "../"
 base = 'cr.out2.'
 ds2 = yt.load(fn2+base+str(i).zfill(5)+'.athdf')

 fn3 = "../../../../beta10/res256/streaming/"   # beta = 10
 base = 'cr.out1.'
 ds3 = yt.load(fn3+base+str(i).zfill(5)+'.athdf')

 fn4 = "../../../../beta10/res256/streaming/"
 base = 'cr.out2.'
 ds4 = yt.load(fn4+base+str(i).zfill(5)+'.athdf')


 ds.add_field(('gas', u'pcr'), function = _pcr, units="g/(cm*s**2)", dimensions=dimensions.pressure)
 ds.add_field(('gas', u'densratio'), function = _densratio, units="",display_name=r"$\rho / \rho_{i}$", dimensions=dimensions.density)
 ds2.add_field(("gas","crscale"), function=crscale,display_name=r"$L_{CR}/L_{0}$",units="")

 ds3.add_field(('gas', u'pcr'), function = _pcr, units="g/(cm*s**2)", dimensions=dimensions.pressure)
 ds3.add_field(('gas', u'densratio'), function = _densratio, units="",display_name=r"$\rho / \rho_{i}$", dimensions=dimensions.density)
 ds4.add_field(("gas","crscale"), function=crscale,display_name=r"$L_{CR}/L_{0}$",units="")
 
 grad_fields = ds.add_gradient_fields(("gas","pcr"))
 grad_fields3 = ds3.add_gradient_fields(("gas","pcr"))

 # We don't need to do the same for the pressure field because yt already
 # has pressure gradient fields. Now, define the "degree of hydrostatic
 # equilibrium" field.


 def _misalign(field, data):
   # Remember that g is the negative of the potential gradient
   gx = data[("gas", "pcr_gradient_x")]
   bx = data[("gas", "magnetic_field_x")]
   gy = data[("gas", "pcr_gradient_y")]
   by = data[("gas", "magnetic_field_y")]
   gz = data[("gas", "pcr_gradient_z")]
   bz = data[("gas", "magnetic_field_z")]
   noalign = np.sqrt(gx * gx + gy * gy + gz * gz)*data['magnetic_field_magnitude']
   align = np.sqrt((gx*bx)**2.0 + (gy*by)**2.0 + (gz*bz)**2.0)
   h = align/noalign
   return h


 ds.add_field(
   ("gas", "misalign"),
   function=_misalign,
   units="",
   take_log=True,
   display_name=r"$\frac{|v_{A} \cdot \nabla P_{CR}|}{|v_{A} \nabla P_{CR}|}$",
   sampling_type="cell",
 )
  
 ds3.add_field(
   ("gas", "misalign"),
   function=_misalign,
   units="",
   take_log=True,
   display_name=r"$\frac{|v_{A} \cdot \nabla P_{CR}|}{|v_{A} \nabla P_{CR}|}$",
   sampling_type="cell",
 )



 fig = plt.figure()
 
 # See http://matplotlib.org/mpl_toolkits/axes_grid/api/axes_grid_api.html
 # These choices of keyword arguments produce a four panel plot that includes
 # four narrow colorbars, one for each plot.  Axes labels are only drawn on the
 # bottom left hand plot to avoid repeating information and make the plot less
 # cluttered.
 grid = AxesGrid(fig, (0.075,0.075,0.85,0.85),
                 nrows_ncols = (3, 2),
                 axes_pad = 0.35,
                 label_mode = "L",
                 share_all = True,
                 cbar_location="right",
                 cbar_mode="edge",
                 cbar_size="3%",
                 cbar_pad="1%",
 		direction='row')
 
 
 fsize = 26
 
 fields = ['densratio', 'misalign', 'crscale']
 
 for j in range(3):
   plotvar = fields[j]
   if (j==0):
     varmin = 0.3
     varmax = 3.0
     slc = yt.SlicePlot(ds, 'z', plotvar,fontsize=fsize)
     slc.set_zlim(plotvar, varmin, varmax)
     slc.set_xlabel('x (kpc)')
     slc.set_ylabel('y (kpc)')
     slc.hide_axes()
     slc.set_cmap(field=plotvar, cmap='kelp')
     slc.set_zlim(plotvar, varmin, varmax)
     slc.set_background_color(plotvar)
     slc.set_log(plotvar,False)
     slc.annotate_line_integral_convolution(("gas", "magnetic_field_x"), ("gas", "magnetic_field_y"), kernellen = 100, lim=(0.5, 0.65),alpha = 0.6)
     slc.annotate_title(r"$\beta \sim 1$")
     plot = slc.plots[plotvar]
     plot.figure = fig
     plot.axes = grid[0].axes
     plot.cax = grid.cbar_axes[0]
     slc._setup_plots()
 
     slc = yt.SlicePlot(ds3, 'z', plotvar,fontsize=fsize)
     slc.set_zlim(plotvar, varmin, varmax)
     slc.set_xlabel('x (kpc)')
     slc.set_ylabel('y (kpc)')
     slc.hide_axes()
     slc.set_cmap(field=plotvar, cmap='kelp')
     slc.set_zlim(plotvar, varmin, varmax)
     slc.set_background_color(plotvar)
     slc.set_log(plotvar,False)
     slc.annotate_line_integral_convolution(("gas", "magnetic_field_x"), ("gas", "magnetic_field_y"), kernellen = 100, lim=(0.5, 0.65),alpha = 0.6)
     slc.annotate_title(r"$\beta \sim 10$")
     plot = slc.plots[plotvar]
     plot.figure = fig
     plot.axes = grid[1].axes
     plot.cax = grid.cbar_axes[j]
     slc._setup_plots()
   if (j==1):
     varmin = 1.e-2
     varmax = 1.e0
    # varmin = 1.e3
    # varmax = 1.e6
     slc = yt.SlicePlot(ds, 'z', plotvar,fontsize=fsize)
     slc.set_zlim(plotvar, varmin, varmax)
     slc.set_cmap(field=plotvar, cmap='kamae')
     slc.set_xlabel('x (kpc)')
     slc.set_ylabel('y (kpc)')
     slc.hide_axes()
     slc.set_background_color(plotvar)
     slc.set_log(plotvar, True)
     plot = slc.plots[plotvar]
     plot.figure = fig
     plot.axes = grid[2].axes
     plot.cax = grid.cbar_axes[j]
     slc._setup_plots()
 
     slc = yt.SlicePlot(ds3, 'z', plotvar,fontsize=fsize)
     slc.set_zlim(plotvar, varmin, varmax)
     slc.set_cmap(field=plotvar, cmap='kamae')
     slc.set_xlabel('x (kpc)')
     slc.set_ylabel('y (kpc)')
     slc.hide_axes()
     slc.set_background_color(plotvar)
     slc.set_log(plotvar, True)
     plot = slc.plots[plotvar]
     plot.figure = fig
     plot.axes = grid[3].axes
     plot.cax = grid.cbar_axes[j]
     slc._setup_plots()
   if (j==2):
     varmin = 1.e-1
     varmax = 1.e3
     slc = yt.SlicePlot(ds2, 'z', plotvar,fontsize=fsize)
     slc.set_zlim(plotvar, varmin, varmax)
     slc.set_cmap(field=plotvar, cmap='hot_r')
     slc.set_xlabel('x (kpc)')
     slc.set_ylabel('y (kpc)')
     slc.hide_axes()
     slc.set_background_color(plotvar)
     slc.set_log(plotvar, True)
     plot = slc.plots[plotvar]
     plot.figure = fig
     plot.axes = grid[4].axes
     plot.cax = grid.cbar_axes[j]
     slc._setup_plots()
 
     slc = yt.SlicePlot(ds4, 'z', plotvar,fontsize=fsize)
     slc.set_zlim(plotvar, varmin, varmax)
     slc.set_cmap(field=plotvar, cmap='hot_r')
     slc.set_xlabel('x (kpc)')
     slc.set_ylabel('y (kpc)')
     slc.hide_axes()
     slc.set_background_color(plotvar)
     slc.set_log(plotvar, True)
     plot = slc.plots[plotvar]
     plot.figure = fig
     plot.axes = grid[5].axes
     plot.cax = grid.cbar_axes[j]
 #    slc._setup_plots()

     slc._setup_plots()
     #fig.set_size_inches(14,14)
     fig.set_size_inches(12,16)
 
 
 # For each plotted field, force the SlicePlot to redraw itself onto the AxesGrid
 # axes.
 #for i, field in enumerate(fields):
 #    plot = slc.plots[field]
 #    plot.figure = fig
 #    plot.axes = grid[i].axes
 #    plot.cax = grid.cbar_axes[i]
 
 # Finally, redraw the plot on the AxesGrid axes.
 plt.savefig('multiplot_beta1_beta10_'+str(i).zfill(5)+'.pdf')
