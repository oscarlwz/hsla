
import numpy as np
from astropy.table import Table, Column, vstack 
from bokeh.plotting import figure, show, output_file, vplot
from bokeh.models import BoxAnnotation


def plot_interactive(target): 
 
    data = Table.read(target+'_coadd_FUVM_final_all.fits.gz') 

    output_file(target+'_interactive.html', title=target+' Interactive Quick Look') 

    TOOLS = "pan,wheel_zoom,box_zoom,resize,reset,save"
    p1 = figure(plot_width=1000, tools=TOOLS, y_range=[0,3e-14])

    p1.line(data['WAVE'], data['FLUX'], color='#AA0000', legend='Flux', alpha=0.5)
    p1.line(data['WAVE'], data['ERROR'], color='#009900', legend='Error', alpha=0.5)
    mid_box = BoxAnnotation(plot=p1, left=1000, right=1500, bottom=0, top=1e-14, fill_alpha=0.1, fill_color='green')

    p1.title = target 
    p1.grid.grid_line_alpha=0.5
    p1.xaxis.axis_label = 'Wavelength'
    p1.yaxis.axis_label = 'Flux'
    
    show(p1)  # open a browser
