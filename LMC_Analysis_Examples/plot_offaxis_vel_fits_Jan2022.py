import yt
import trident
from trident import LightRay
import aplpy
import numpy as np
#import matplotlib as plt
from yt.units.yt_array import YTQuantity
from yt import YTArray
import h5py
import matplotlib.pyplot as plt
from astropy.io import fits
#ts = yt.load("windCRs_hdf5_plt_cnt_0025")
#ts = yt.DatasetSeries("windCRs_hdf5_plt_cnt_0005")

from scipy.optimize import curve_fit

def func(x,a,b,c):
    return a * np.exp(-b * x) + c

xdata = [0.025,0.091,0.223,0.364,0.496,0.709,0.84,1.043,1.316,1.569,1.893,2.227,2.642,3.107,3.573,4.038,4.464,4.879,5.294,5.698,6.144,6.579,6.953,7.308,7.763,8.381,9.008,9.524,9.975]

ydata = [21.671,25.196,30.809,35.64,40.47,45.822,49.086,53.264,58.355,62.533,66.449,69.582,72.715,75.718,77.546,79.373,80.287,81.07,81.593,81.984,82.115,82.245,82.245,81.984,82.115,81.854,81.723,81.07,80.94]


x_arr = np.arange(0.01,10,0.01)
popt, pcov = curve_fit(func, xdata, ydata)
#print(popt)
plt.plot(xdata,ydata,label="From Salem2015")
plt.plot(x_arr, func(x_arr, *popt), 'r-', label='fit: a=%5.3f, b=%5.3f, c=%5.3f' % tuple(popt))
plt.xlabel("x (kpc)")
plt.ylabel("Velocity (km/s)")
plt.legend()
plt.savefig("Rotation_Velocity_Fit_Salem2015.pdf")
plt.close()




yt.enable_parallelism()
ts = yt.DatasetSeries("windCRs_hdf5_plt_cnt_0050*",parallel=False)
#ts = yt.DatasetSeries("RP_Streaming",parallel=False)
rho_ex = []
times = []
def _ramPres(field, data):
        velrot_y = -7.5e6*data['x']/np.sqrt((data['x']**2.0 + data['y']**2.0))*YTQuantity(1,'cm/s')
        velrot_x = 7.5e6*data['y']/np.sqrt((data['x']**2.0 + data['y']**2.0))*YTQuantity(1,'cm/s')
        val = data['density']*((data['velocity_x'] - velrot_x)**2.0 + (data['velocity_y'] - velrot_y)**2.0 + (data['velocity_z'])**2.0)
        return val

def _ramPres_2(field, data):
        rad = np.sqrt((data['x']**2.0 + data['y']**2.0))
        velrot_y = 1.e5*func(rad,*popt)*data['x']/rad*YTQuantity(1,'cm/s')
        velrot_x = 1.e5*func(rad,*popt)*data['y']/rad*YTQuantity(1,'cm/s')
        val = 1.65e-28*YTQuantity(1,'g/cm**3')*((220.0*1e5*data['ones']*YTQuantity(1,'cm/s') - velrot_x)**2.0 + (159.0*1e5*data['ones']*YTQuantity(1,'cm/s') - velrot_y)**2.0 + (192.0*1e5*data['ones']*YTQuantity(1,'cm/s'))**2.0)
        return val

yt.add_field(("gas","ramPres"), function=_ramPres_2,units="dyne/cm**2")



for ds in ts.piter():
        dd = ds.all_data()


        ds.periodicity=(True,True,True)
        ds.coordinates.x_axis[1] = 0
        ds.coordinates.y_axis[1] = 2

        print(ds.derived_field_list)
        L = [-0.561,0.099,0.822] # vector normal to cutting plane
        north_vector = [0,0,1]

       # cut_data = dd.cut_region(['(obj["igm "] > 0.95) & (obj["z"].in_units("kpc") > -5) & (obj["z"].in_units("kpc") < 0.0)'])
        cut_data = dd.cut_region('(obj["ism "] > 0.6)')

       # proj = yt.ProjectionPlot(ds, "z", "H_p0_number_density")
       # proj = yt.OffAxisProjectionPlot(ds, L, "H_p1_number_density", width=(40, 'kpc'),fontsize=20)
       # proj = yt.OffAxisProjectionPlot(ds, L, "H_nuclei_density", width=(60, 'kpc'),fontsize=20)
       # proj = yt.OffAxisProjectionPlot(ds, L, "ramPres",center = (0,0,0), width=(40, 'kpc'))
        proj = yt.SlicePlot(ds, "z", "ramPres",center = (0,0,0), width=(40, 'kpc'),data_source=cut_data)
       # proj.set_zlim("H_nuclei_density",1.E18, 1.E22)
       # proj.set_zlim("H_p1_number_density",1.E18, 1.E21)
       # proj.set_zlim("H_p0_number_density",1.E-5, 1.E0)
        proj.set_cmap(field="ramPres", cmap='kamae')
        proj.set_colorbar_label("ramPres", r"P$_{ram}$ (dyne/cm$^{2}$)")
        proj.set_zlim("ramPres",1.E-13, 4.E-13)
        proj.annotate_timestamp()
        proj.save()

        """
      #  d.set_cmap(field="H_nuclei_density", cmap='viridis')
      #  d.set_colorbar_label("H_nuclei_density", "Column Density (cm$^{-2}$)")
      #  d.annotate_particles((10, 'Mpc'),minimum_mass=1.9E37,p_size=0.8)
       # d.set_width((30,'kpc'))
       # d.set_zlim('H_nuclei_density',1e17,1e22)
       # d.annotate_timestamp()
      #  d.annotate_contour('H_nuclei_density',ncont = 1,clim=[5e17,5e17],take_log = True, label = True,
       #         plot_args={"colors": "white", "linewidths": 0.5})
       # d.annotate_contour('H_nuclei_density',ncont = 1,clim=[1e18,1e18],take_log = True, label = True,
       #         plot_args={"colors": "cyan", "linewidths": 0.6})
       # d.annotate_contour('H_nuclei_density',ncont = 1,clim=[1e19,1e19],take_log = True, label = True,
       #         plot_args={"colors": "yellow", "linewidths": 0.7})
       # d.annotate_contour('H_nuclei_density',ncont = 1,clim=[1e20,1e20],take_log = True, label = True,
       #         plot_args={"colors": "red", "linewidths": 0.8})
       # d.annotate_contour('density',ncont = 3,clim=[3e-6,5e-6],take_log = True, label = True,
       #         plot_args={"linewidths": 1})
       # d.annotate_contour('density',ncont = 2,clim=[5e-6,5e-5],take_log = True, label = True,
       #         plot_args={"linewidths": 1})
        #d.annotate_contour('density',ncont = 1,clim=[1e-6],take_log = True, label = True,
        #        plot_args={"colors": "yellow",
        #                   "linewidths": 2})
        #d.annotate_contour('density',ncont = 1,clim=[1e-5],take_log = True, label = True,
        #        plot_args={"colors": "red",
        #                   "linewidths": 2})
       # d.annotate_text((0, 15), 'Many winds launched after 400 Myrs', coord_system='plot')
       # d.save()
        proj.save()

        p = yt.SlicePlot(ds, "z", "pressure")
        p.set_zlim('pressure',1.E-14,1.E-11)
        p.annotate_timestamp()
        p.save()
        t = yt.SlicePlot(ds, "z", "temperature")
  #      t.set_zlim('temperature',1.E4,1.E8)
        t.annotate_timestamp()
        t.save()

        v = yt.SlicePlot(ds, "z", "velocity_magnitude")
        v.set_zlim('velocity_magnitude',1.E5,1.E8)
        v.annotate_timestamp()
        v.annotate_quiver("velocity_x","velocity_y",16)
 #       v.set_logscale('False')
        v.save()
        rho_ex.append(dd.quantities.extrema("pressure"))
        times.append(ds.current_time.in_units("s"))
#        print(dd)

       """
       # prj_fits = yt.FITSOffAxisProjection(ds,L,'H_p0_number_density',width=(40, 'kpc'))
       # prj_fits = yt.FITSOffAxisProjection(ds,L,'H_p0_number_density',width=(20, 'kpc'))
        prj_fits = yt.FITSOffAxisProjection(ds,L,'ramPres',center = (0,0,0),width=(20, 'kpc'),weight_field="ones",data_source=cut_data)
       # sky_center = [82.245,-69.5] # in degrees
       # sky_center = [78.76,-69.19] # in degrees, what I previously sent to Yong
        sky_center = [79.0,-68.68] # in degrees, what Jack told me to use
       # sky_center = [82.245-16.19,-69.5+16.19] # in degrees
       # sky_center = [30.,45.] # in degrees
        sky_scale = (4123.71, "arcsec/kpc") # could also use a YTQuantity
        prj_fits.create_sky_wcs(sky_center, sky_scale, ctype=["RA---TAN","DEC--TAN"],replace_old_wcs=True)
       # print(prj_fits["density"].header)
        prj_fits.writeto("LMC_ram_RP_rotate_Jan2022_newSkyCenter.fits", clobber=True)
       # prj_fits.write_fits("LMC_ram.fits", clobber=True,sky_center=sky_center,sky_scale=sky_scale)
        image_data = fits.getdata('LMC_ram_RP_rotate_Jan2022_newSkyCenter.fits')
        print(type(image_data))
        print(image_data.shape)

        fig = aplpy.FITSFigure("LMC_ram_RP_rotate_Jan2022_newSkyCenter.fits")
        fig.add_grid()
 #       fig.set_title(r'H I Column Density ($cm^{-2}$)',size=24)
       # fig.show_colorscale(cmap="kamae_r",stretch = 'linear') #,vmin=7.5E-14, vmax = 1.5E-13)
        fig.show_colorscale(cmap="kamae_r",stretch = 'linear', vmin=1E-13, vmax = 3.5E-13)
        #fig.set_axis_labels_font(family="serif", size=20)
        fig.tick_labels.set_xposition("top")
        #fig.set_tick_labels_font(family="serif", size=18)
        fig.ticks.set_xspacing(8.0)
        fig.ticks.set_yspacing(8.0)
        fig.add_colorbar()

       # fig.tick_labels.set_xformat('hh:mm')
        fig.tick_labels.set_xformat('dd')
        fig.tick_labels.set_yformat('dd')
       # fig.set_theme('pretty')
        fig.colorbar.show(log_format=False)
       # fig.colorbar.set_axis_label_text(r'H I Column Density ($cm^{-2}$)')
       # fig.colorbar.set_axis_label_text(r'H II Column Density ($cm^{-2}$)')
        fig.colorbar.set_font(size=20)
       # fig.set_title("Low Gas Mass LMC, Ram Pressure Only")
     #   fig.colorbar.set_label_properties(size=20)
       # fig.colorbar.set_ticks([1.e18,1.e19,1.e20,1.e21])
       # fig.set_theme('publication')
       # fig.savefig('LMC_mock_RP_rotate_nearSide_5kpc_95.pdf')
        fig.savefig('LMC_mock_RP_rotate_Jan2022_newSkyCenter.pdf')

        """

       # prj_fits = yt.FITSOffAxisProjection(ds,L,'H_p0_number_density',width=(40, 'kpc'))
       # prj_fits = yt.FITSOffAxisProjection(ds,L,'H_p0_number_density',width=(20, 'kpc'))
        prj_fits = yt.FITSOffAxisProjection(ds,L,'O_p5_number_density',width=(20, 'kpc'))
        sky_center = [82.245,-69.5] # in degrees
       # sky_center = [30.,45.] # in degrees
        sky_scale = (4123.71, "arcsec/kpc") # could also use a YTQuantity
        prj_fits.create_sky_wcs(sky_center, sky_scale, ctype=["RA---TAN","DEC--TAN"],replace_old_wcs=True)
       # print(prj_fits["density"].header)
        prj_fits.writeto("LMC_ram.fits", clobber=True)
       # prj_fits.write_fits("LMC_ram.fits", clobber=True,sky_center=sky_center,sky_scale=sky_scale)
        image_data = fits.getdata('LMC_ram.fits')
        print(type(image_data))
        print(image_data.shape)

        fig = aplpy.FITSFigure("LMC_ram.fits")
        fig.add_grid()
        fig.show_colorscale(cmap="viridis",stretch = 'log',vmin=1.E13, vmax = 1.E15)
        fig.set_axis_labels_font(family="serif", size=12)
        fig.set_tick_labels_font(family="serif", size=6)
        fig.add_colorbar()
       # fig.colorbar.set_axis_label_text(r'HI Column Density ($cm^{-2}$)')
        fig.colorbar.set_axis_label_text(r'O VI Column Density ($cm^{-2}$)')
       # fig.set_theme('publication')
        fig.savefig('LMC_mock.pdf')



        prj_fits = yt.FITSOffAxisProjection(ds,L,'H_nuclei_density',width=(15, 'kpc'))
        sky_center = [82.245,-69.5] # in degrees
       # sky_center = [30.,45.] # in degrees
        sky_scale = (4123.71, "arcsec/kpc") # could also use a YTQuantity
        prj_fits.create_sky_wcs(sky_center, sky_scale, ctype=["GLON","GLAT"],replace_old_wcs=True)
       # print(prj_fits["density"].header)
        prj_fits.writeto("LMC_ram_galactic.fits", clobber=True)
       # prj_fits.write_fits("LMC_ram.fits", clobber=True,sky_center=sky_center,sky_scale=sky_scale)
        image_data = fits.getdata('LMC_ram_galactic.fits')
        print(type(image_data))
        print(image_data.shape)

        fig = aplpy.FITSFigure("LMC_ram_galactic.fits")
        fig.add_grid()
        fig.show_colorscale(cmap="viridis",stretch = 'linear',vmin=1.E13, vmax = 1.E15)
        fig.set_axis_labels_font(family="serif", size=12)
        fig.set_tick_labels_font(family="serif", size=6)
        fig.add_colorbar()
        fig.colorbar.set_axis_label_text(r'Column Density ($cm^{-2}$)')
       # fig.set_theme('publication')
        fig.savefig('LMC_mock_galactic.pdf')
        """

 #       plt.show(image_data,cmap='gray')
 #       plt.colorbar()
 #       plt.savefig('LMC_ram.pdf')

rho_ex = np.array(rho_ex)
print(rho_ex)

 #       plt.show(image_data,cmap='gray')
 #       plt.colorbar()
 #       plt.savefig('LMC_ram.pdf')

rho_ex = np.array(rho_ex)
print(rho_ex)
