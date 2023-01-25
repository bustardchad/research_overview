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
import matplotlib as plt
from numpy import *
import matplotlib.cm as cm
import matplotlib.mlab as mlab
from yt.fields import interpolated_fields
from yt.fields.field_detector import \
    FieldDetector
from yt.utilities.linear_interpolators import \
    BilinearFieldInterpolator
ts = yt.DatasetSeries("windCRs_hdf5_plt_cnt_0050")
rho_ex = []
times = []

nenh = loadtxt('nenh.txt')
tempTxt = loadtxt('Temps.txt')
nHTxt = loadtxt('NH.txt')




for ds in ts:
 #       time = ds.current_time.in_units("Kyrs")
        dd = ds.all_data()

       # ds.add_interpolated_field(name="nenh_interp",units="",table_data=nenh,axes_data=[tempTxt,nHTxt],axes_fields = ['temperature','H_nuclei_density'])
        ds.periodicity=(True,True,True)
 #       print(ds.derived_field_list)
        ds.coordinates.x_axis[1] = 0
        ds.coordinates.y_axis[1] = 2

        def _nenh_interp(field, data):
          #  nH = 0.76*(data['number_density'])
          #  temp = data['temperature']


            # use yt's bilinear interpolator:
            interp = BilinearFieldInterpolator(nenh,boundaries=[tempTxt,nHTxt],field_names=['temperature','H_nuclei_density'], truncate=True)

            field_data = interp(data)
            field_data[field_data < 0] = 0.01
            return(field_data)

        ds.add_field(('gas','nenh_interp'), function=_nenh_interp, units="", take_log=False, display_name='$n_{e} / n_{H}$', sampling_type="cell")


        def _rotationMeasure(field, data):
          #  nH = 0.76*(data['number_density'])
          #  temp = data['temperature']


            # use yt's bilinear interpolator:
            interp = BilinearFieldInterpolator(nenh,boundaries=[tempTxt,nHTxt],field_names=['temperature','H_nuclei_density'], truncate=True)
            if (any(interp(data) < 0)):
                print(interp(data))
           # field_data = interp(data)
            field_data = 1.2
            return 0.812*field_data*(data['H_nuclei_density'])*((-0.561*data['magnetic_field_x'] + 0.099*data['magnetic_field_y'] + 0.822*data['magnetic_field_z'])/1.e-6)*(1./3.0856e18)



        ds.add_field(('gas','rotationMeasure'), function=_rotationMeasure, units="gauss/cm**3", take_log=False, display_name='Rotation Measure', sampling_type="cell")

        def _losBField(field,data):
            return (-0.561*data['magnetic_field_x'] + 0.099*data['magnetic_field_y'] + 0.822*data['magnetic_field_z'])/1.e-6

        ds.add_field(('gas','losBField'), function=_losBField, units="gauss", take_log=False, display_name='Line of Sight B Field', sampling_type="cell")



        L = [-0.561,0.099,0.822] # vector normal to cutting plane
        region = ds.box([0.0,-6.17E22,-3.0856E22],[1.2E23,6.17E22,3.0856E22])
       # d = yt.ProjectionPlot(ds, "z", "density",origin='native',center=[0.0,0.0,0.0])

       # d = yt.ProjectionPlot(ds, "z", "rotationMeasure",fontsize = 20)
        L = [-0.561,0.099,0.822] # vector normal to cutting plane
        d = yt.OffAxisProjectionPlot(ds, L, "nenh_interp",weight_field="density",center=[0.0,0.0,0.0], width=(10, 'kpc'),fontsize=20)
       # d = yt.OffAxisProjectionPlot(ds, L, "losBField", width=(30, 'kpc'),fontsize=20,weight_field='H_nuclei_density')
       # d = yt.ProjectionPlot(ds, "z", "density")
        d.set_cmap(field="nenh_interp", cmap='PuOr')
        d.set_colorbar_label("nenh_interp", "$n_{e}/n_{H}$")
        #d.set_width((60,'kpc'))
       # d.set_zlim('rotationMeasure',-200,200)
        d.annotate_timestamp()
        d.save()
        """

        prj_fits = yt.FITSOffAxisProjection(ds,L,'rotationMeasure',center = (0,0,0),width=(10, 'kpc'))
       # sky_center = [82.245,-69.5] # in degrees
        sky_center = [79.0,-68.68] # in degrees
        # sky_center = [82.245-16.19,-69.5+16.19] # in degrees
        # sky_center = [30.,45.] # in degrees
        sky_scale = (4123.71, "arcsec/kpc") # could also use a YTQuantity
        prj_fits.create_sky_wcs(sky_center, sky_scale, ctype=["RA---TAN","DEC--TAN"],replace_old_wcs=True)
        # print(prj_fits["density"].header)
        prj_fits.writeto("LMC_Ram_evolved_Faraday.fits", clobber=True)
        # prj_fits.write_fits("LMC_ram.fits", clobber=True,sky_center=sky_center,sky_scale=sky_scale)
        image_data = fits.getdata('LMC_Ram_evolved_Faraday.fits')
        print(type(image_data))
        print(image_data.shape)

        fig = aplpy.FITSFigure("LMC_Ram_evolved_Faraday.fits")
        fig.add_grid()
        fig.show_colorscale(cmap="RdBu",stretch = 'linear',vmin=-100, vmax = 100)
       # fig.set_axis_labels_font(family="serif", size=18)
        fig.tick_labels.set_xposition("top")
       # fig.set_tick_labels_font(family="serif", size=12)
        fig.ticks.set_xspacing(8.0)
        fig.ticks.set_yspacing(8.0)
        fig.add_colorbar()
        # fig.tick_labels.set_xformat('hh:mm')
       # fig.tick_labels.set_xformat('dd')
        fig.tick_labels.set_xformat('hh:mm')
        fig.tick_labels.set_yformat('dd')
        # fig.set_theme('pretty')
        fig.colorbar.set_axis_label_text(r'Rotation Measure $\phi$ (rad $ m^{-2}$)')
        # fig.colorbar.set_axis_label_text(r'H II Column Density ($cm^{-2}$)')
        # fig.colorbar.set_font(size=14)
        # fig.set_theme('publication')
        fig.savefig('LMC_mock_Faraday_Ram_evolved.pdf')
        """
