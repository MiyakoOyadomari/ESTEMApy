import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.patches as pat
from astropy.io import fits
import copy, optparse, os, sys, io
import re, string, pprint, math
import time
import traceback
import logging
import matplotlib.cm as cm
from ctypes import *
from numpy.ctypeslib import ndpointer
import scipy
from scipy.stats import norm
from numpy.random import *
import gc

global freq_if, ch_sep, c, convert_f
c = 2.99792458e8 # Speed of the light [m/s]

def parse_inp(filename):
    # form a hash of parameter names and values parsed from an input file.
    # Don't worry about parameter types (lists, strings, etc.) as we sort that
    # out after
    INPUTFILE = open(filename, "r")
    control = dict()
    
    # a few useful regular expressions
    newline = re.compile(r'\n')
    space = re.compile(r'\s')
    char = re.compile(r'\w') #a-xA-Z0-9
    comment = re.compile(r'#.*')
    
    # parse the input file assuming '=' is used to separate names from values
    for line in INPUTFILE:
        if char.match(line):
            line = comment.sub(r'', line) #remove a comment line
            line = line.replace("'", '')  #remove single quatation marks
            (param, value) = line.split('=')
            
            param = newline.sub(r'', param)
            param = param.strip()
            param = space.sub(r'', param)
            
            value = newline.sub(r'', value)
            value = value.strip()
            value = space.sub(r'', value)
            valuelist = value.split(',')
            #print param,'=',valuelist
            control[param] = valuelist
    INPUTFILE.close()
    return control

def chekin(control):
    global out_dir, output_file_name, pipelog, output_file ,CP_file_name,CP_Hist_name
    global raw_data_binary,savefile_Hist1,savefile_Hist2, initial_image_file
    global plain_file,flux_cal_file_name, CP_result_file, CF_info, RingWidth
    global tmask ,image_file_name, x_axis, y_axis
    global initial_center, initial_radius, area_size_x, area_size_y
    global cum_flux_up, cum_flux_lo, search_range, noise, itn ,cbmax
    out_dir = control.get('out_dir', [False])[0]
    output_file_name = control.get('output_file_name', [False])[0]
    output_file = out_dir + '/' + output_file_name
    pipelog = output_file + '.pipelog'
    raw_data_binary = output_file + '_raw_data.array'
    initial_image_file = output_file + '_initial_image.png'
    plain_file = output_file + '_plain_CP.txt'
    flux_cal_file_name = output_file + '_Radius_Report.txt'
    CP_file_name = out_dir + '/tmep_CenterPosition.txt'
    CP_Hist_name = out_dir + '/CP_Histgram.txt'
    savefile_Hist1 = out_dir + '/CP_Histgram.png'
    savefile_Hist2 = out_dir + '/CP_Histgram_CrossSum.png'
    CP_result_file = output_file + '_CenterPosition_Result.txt'
    CF_info = out_dir + '/converted_function.txt'
    RingWidth = output_file + '_RignWidth.txt'
    
    tmask = control.get('tmask', [False])

    image_file_name = control.get('image_file_name', [False])[0]
    x_axis = control.get('x_axis', [0])
    y_axis = control.get('y_axis', [0])
    initial_center = control.get('initial_center', [])
    for i in range (len(initial_center)):
        initial_center[i] = int(initial_center[i])

    initial_radius = int(control.get('initial_radius', [200])[0])
    area_size_xy = control.get('area_size_xy', [1000,1000])
    area_size_x = int(area_size_xy[0])
    area_size_y = int(area_size_xy[1])

    cum_flux_lo = float(control.get('cum_flux_lo', [0.02])[0])
    cum_flux_up = float(control.get('cum_flux_up', [0.98])[0])
    search_range = int(control.get('search_range', [20])[0])
    noise = float(control.get('noise', [])[0])
    itn = int(control.get('itn', [])[0])
    cbmax = float(control.get('cbmax', [])[0])



def plot_image(im_data_array,x_axis,y_axis, savefile):
    fig = plt.figure()  # fig,ax = plt.subplots()
    ax1 = fig.add_subplot(111)
    imp = ax1.imshow(im_data_array,cmap='jet', norm=colors.LogNorm())
    cbar = plt.colorbar(imp)
    plt.ylim(y_axis[0],y_axis[1])
    plt.xlim(x_axis[0],x_axis[1])
    plt.xlabel('channel')
    plt.ylabel('channel')
    cbar.set_label('Jy/beam km/s')
    # Save the map
    logger.info('# Save the map')
    logger.info(savefile)
    plt.savefig(savefile,dpi=300)
    plt.show()

def plot_circle_on_map(im_data_array,x_axis,y_axis,radius):
    fig = plt.figure()  # fig,ax = plt.subplots()
    ax1 = fig.add_subplot(111)
    imp = ax1.imshow(im_data_array,cmap='jet', norm=colors.LogNorm())
    cbar = plt.colorbar(imp)
    plt.ylim(y_axis[0],y_axis[1])
    plt.xlim(x_axis[0],x_axis[1])
    plt.xlabel('channel')
    plt.ylabel('channel')
    cbar.set_label('Jy/beam km/s')
    plt.scatter(x=[x_axis[1]/2], y=[y_axis[1]/2], c='black', s=70, marker='+')
    c = pat.Circle(xy = (x_axis[1]/2, y_axis[1]/2), radius = radius, color = "black",fill=False)
    ax1.add_patch(c)
    plt.show()

def plot_xy(x, y, xlabel, ylabel, xy_title):
    fig = plt.figure()  # fig,ax = plt.subplots()
    ax1 = fig.add_subplot(111)
    ax1.scatter(x,y)
    ax1.set_title(xy_title)
    ax1.set_xlabel(xlabel)
    ax1.set_ylabel(ylabel)
    plt.show()

def plot_map_graph(im_data_array,x_axis,y_axis,radius,x,y,radius_lo,radius_up,im_title,xy_title,savefile, pop,cmap, cblb, yl,vmin=None, vmax=None, norm=None):
    fig = plt.figure(figsize=(15,5))  # fig,ax = plt.subplots()
    ax1 = fig.add_subplot(121)

    if task_no == 7:
        imp = ax1.imshow(im_data_array,cmap='jet', norm=colors.LogNorm())
    else:
        imp = ax1.imshow(im_data_array,cmap=cmap, vmin=vmin, vmax=vmax, norm=norm)


    cbar = plt.colorbar(imp)
    plt.ylim(y_axis[0],y_axis[1])
    plt.xlim(x_axis[0],x_axis[1])
    plt.xlabel('channel')
    plt.ylabel('channel')
    ax1.set_title(im_title)
    cbar.set_label(cblb)
    #    cbar.set_label('Jy/beam km/s')
    plt.scatter(x=[x_axis[1]/2], y=[y_axis[1]/2], c='black', s=70, marker='+')
    c = pat.Circle(xy = (x_axis[1]/2, y_axis[1]/2), radius = radius, color = "black",fill=False)
    ax1.add_patch(c)
    
    ax2 = fig.add_subplot(122)
    ax2.scatter(x,y)
    ax2.set_title(xy_title)
    ax2.set_xlabel('radius (ch)')
    ax2.set_ylabel(yl)
    ax2.axvline(radius)
    ax2.axvline(radius_lo)
    ax2.axvline(radius_up)

    plt.savefig(io.BytesIO())
    
    ax1_pos = ax1.get_position()
    cax1_pos0 = cbar.ax.get_position()
    cax1_pos1 = [cax1_pos0.x0, ax1_pos.y0, cax1_pos0.x1 - cax1_pos0.x0, ax1_pos.y1 - ax1_pos.y0]
    cbar.ax.set_position(cax1_pos1)
    plt.savefig(savefile,dpi=300)
    if pop == 'Y':
        plt.show()
    plt.close()


def cal_flux(im_data,center_x,center_y,radius):
    cal_flux_list = []
    data_slice = im_data[center_y - radius:center_y + radius + 1, center_x - radius:center_x + radius + 1].copy()
    # making a array_inout;
    # 1: The element is inner the circle,  0: The element is out of the circle
    array_inout = np.zeros((radius*2+1,radius*2+1))  # making all "0" array
    
    i = 0; j = 0
    for i in range(radius*2+1):
        for j in range(radius*2+1):
            RR = (i - radius)**2 + (j - radius)**2
            if RR <= radius**2:
                array_inout[i][j] = 1

    # The sum flux inner the circle
    flux_in_circle = array_inout * data_slice
    if np.sum(flux_in_circle) == 0:
        inclusion_rate = 0
    else:
        inclusion_rate = np.sum(flux_in_circle) / np.sum(im_data)
    D_rate = abs(inclusion_rate - 0.5)
    cal_flux_list = [center_x, center_y, radius, np.sum(im_data),np.sum(flux_in_circle), inclusion_rate,D_rate]
    return cal_flux_list

def cal_radius(im_data,center_x,center_y,fig):
    cal_flux_results = []; cal_flux_list = []
    
    if task_no == 7:
        if os.path.exists(flux_cal_file_name):
            os.remove(flux_cal_file_name)
        flux_cal_file = open(flux_cal_file_name, 'w')
        flux_cal_file.write('# center_x(ch), center_y(ch), radius(ch),total_flux, flux_in_circle, rate, abs(rate-0.5) \n')


    i = 0
    while (initial_radius -1 > i):
        radius = initial_radius - i
        # Calculate the flux in the circle
        cal_flux_list = cal_flux(im_data,center_x,center_y,radius)
        cal_flux_results.append(cal_flux_list)
        if task_no == 7:
            str_cal_flux_list = str(cal_flux_list[0])+','+str(cal_flux_list[1])+','+str(cal_flux_list[2])+','+str(cal_flux_list[3])+','+str(cal_flux_list[4])+','+str(cal_flux_list[5])+','+str(cal_flux_list[6])
            flux_cal_file.write(str_cal_flux_list+'\n')
        i += 1

    if task_no == 7:
        flux_cal_file.close()
    
    The_Radius = list_minimum(cal_flux_results, 6, 2)
    lo_radius = list_lower_limit(cal_flux_results, 5, cum_flux_lo, 2)
    up_radius = list_lower_limit(cal_flux_results, 5, cum_flux_up, 2)
    if not up_radius:
        logger.error('center_xy = ('+str(center_x)+', '+str(center_y)+')')
        logger.error('The normalized cumulative flux does not reach the upper rate.')
        return 0
    else:
        D1_R = up_radius[1] - lo_radius[1]
        D2_R = The_Radius[1] - lo_radius[1]
        D2_D1 = float(D2_R) / float(D1_R)
        logger.info('The_Radius = ' + str(The_Radius[1]) + '(index = ' + str(The_Radius[0]) + ')')
        logger.info('lo_radius = ' + str(lo_radius) + ', up_radius = ' + str(up_radius))
        logger.info('D1_R = ' + str(D1_R) + ', D2_R = ' + str(D2_R) + ', R2/R1 = ' + str(D2_D1))
        if fig == 'Y':
            return The_Radius,lo_radius, up_radius, D1_R, D2_R,D2_D1,cal_flux_results
        else:
            return The_Radius,lo_radius, up_radius, D1_R, D2_R,D2_D1

def list_minimum(list_name, require_j, out_j):
    i = 0; temp_list = []; m = 0;
    for i in range(len(list_name)):
        temp_list.append(list_name[i][require_j])
        i += 1
    min_i = [k for k, v in enumerate(temp_list) if v == min(temp_list)]
    for m in range(len(min_i)):
        min_val = list_name[min_i[m]][out_j]
    return min_i, min_val

def list_lower_limit(list_name, require_j, require_val, out_j):
    flux_abs = 0.5
    i = 0
    for i in range(len(list_name)):
        if flux_abs > abs(list_name[i][require_j] - require_val):
            flux_abs = abs(list_name[i][require_j] - require_val)
            lower_limit_i = i
            out_val = list_name[i][out_j]
        i += 1
    return lower_limit_i, out_val


#    i = 0; temp_list = []
#    for i in range(len(list_name)):
#        if list_name[i][require_j] >= require_val:
#            temp_list.append(list_name[i])
#        i += 1
#    if not temp_list:
#        return 0
#    else:
#        lower_limit_i = temp_list.index(temp_list[-1])
#        out_val = list_name[lower_limit_i][out_j]
#        return lower_limit_i, out_val

def print_start(text):
    logger.info('\n')
    logger.info('*' * 20)
    logger.info('Starting tmask ' + text + ' at ' + today)

def print_end(text):
    logger.info('Ending tmask ' + text + ' at ' + today)

##################################################################
today = time.asctime()


control = parse_inp(sys.argv[1])
if not sys.argv[1]:
    logger.error('Insufficient arguments')
#control = parse_inp('test2_input.txt')
chekin(control)

for i in range(len(tmask)):
    tmask[i] = int(tmask[i])
if (len(tmask) == 1):
    tmask.append(999)

# delete the log file if we are starting again.
if tmask[0] <= 1 <= tmask[1]:
    if os.path.exists(pipelog):
        os.remove(pipelog)

#  Output Messages to not only the log-file but also the consol.
logger = logging.getLogger('PIPELOG')
logger.setLevel(10)
fh = logging.FileHandler(pipelog)
logger.addHandler(fh)
sh = logging.StreamHandler()
logger.addHandler(sh)
########################################

logger.info("\n")
logger.info('*' * 20)
logger.info('Pipeline started on ' + today)
logger.info('Using:', sys.argv[1])
#logger.info('Using:'+'r18316c_v2.txt' )
logger.info("inputs:\n")
logger.info(pprint.pformat(control))
logger.info("\n")
logger.info('doing tmask ' + str(tmask[0]) + ' to '+ str(tmask[1]))


######## tmask 1 ########
### Load a fits file created after "squash2" (4D matrix),
# and creat an array of image data.
if tmask[0] <= 1 <= tmask[1]:
    print_start('1, Load a fits file and creat an array of image data.')
    logger.info('Loading ' +  image_file_name)

    fits_in = fits.open(image_file_name)
    logger.info('fits_in size =' + str(len(fits_in)))
    logger.info('fits_in[0] =' + str(fits_in[0]))
    logger.info('fits_in[1] =' + str(fits_in[1]))

    hdu = fits_in[0]
    fits_data = hdu.data  # Array
    fits_header = hdu.header  # Dictionary
    logger.info('fits_data.shape = ' + str(fits_data.shape))
    raw_data=fits_data[0,0] # raw_data array
    # Save the raw data array in a binary file
    logger.info('# Save the raw data array in a binary file')
    logger.info('file: '+raw_data_binary +'.npy')
    np.save(raw_data_binary, raw_data)


    freq_if = fits_header.get('CRVAL3', [False])
    ch_sep = fits_header.get('CDELT3', [False])
    logger.info('freq_if = ' + str(freq_if) + ', ch_sep = ' + str(ch_sep))
    convert_f = ch_sep * c * 1e-3 / freq_if
    logger.info('convert function = '+str(convert_f))

    if os.path.exists(CF_info):
        os.remove(CF_info)
    
    CF = open(CF_info, 'w')
    CF.write(str(convert_f))
    CF.close()
    
    # Convert the unit of raw data to "Jy/beam km/s".
    logger.info('# Convert the unit of raw data to "Jy/beam km/s')
#    freq_if = fits_header.get('CRVAL3', [False])
#    ch_sep = fits_header.get('CDELT3', [False])
#    logger.info('freq_if = ' + str(freq_if) + ', ch_sep = ' + str(ch_sep))
#
#    convert_f = ch_sep * c * 1e-3 / freq_if
#    logger.info('convert function = '+str(convert_f))
    conv_data = raw_data * convert_f # Converted data array, unit: Jy/beam km/s
 
    # Make a natural logarithmic flux data array from "comv_data".
    min_flux = np.min(conv_data[conv_data > 0])
    max_flux = np.max(conv_data)
    logger.info('min_flux for conv_data: '+str(min_flux))
    logger.info('max_flux for conv_data: '+str(max_flux))
    temp_data = conv_data / min_flux
    log_data = np.where(temp_data > 0, np.log(temp_data), 0)
    min_flux = np.min(log_data[conv_data > 0])
    max_flux = np.max(log_data)
    logger.info('min_flux for log_data: '+str(min_flux))
    logger.info('max_flux for log_data: '+str(max_flux))

    del hdu; del fits_data ; del fits_in
#    del(raw_data);del(conv_data)
#    del(temp_data);del(log_data)
    gc.collect()
    print_end('1')
#######################################
# Laod the binary data when the tmask 1 is skipped.
if tmask[0] != 1:
    # Laod the raw data array
    logger.info('# Laod the raw data array when the tmask 1 is skipped.')
    raw_data_binary_in = raw_data_binary + '.npy'
    logger.info(raw_data_binary_in)
    raw_data = np.load(raw_data_binary_in)
    logger.info('raw_data.shape = '+str(raw_data.shape))

#    fits_in = fits.open(image_file_name)
#    hdu = fits_in[0]
#    fits_header = hdu.header  # Dictionary
#
#    logger.info('# Convert the unit of raw data to "Jy/beam km/s')
#    freq_if = fits_header.get('CRVAL3', [False])
#    ch_sep = fits_header.get('CDELT3', [False])
#    convert_f = ch_sep * c * 1e-3 / freq_if

    # Make a converted data
    CF = open(CF_info, 'r')
    convert_f = CF.readline()
    convert_f = float(convert_f)
    CF.close()
    
    logger.info('convert function = '+str(convert_f))
#    logger.info('freq_if = ' + str(freq_if) + ', ch_sep = ' + str(ch_sep))
    conv_data = raw_data * convert_f # Converted data array, unit: Jy/beam km/s
    
    
    # Make a natural logarithmic flux data array from "comv_data".
    min_flux = np.min(conv_data[conv_data > 0])
    max_flux = np.max(conv_data)
    logger.info('min_flux for conv_data: '+str(min_flux))
    logger.info('max_flux for conv_data: '+str(max_flux))
    temp_data = conv_data / min_flux
    log_data = np.where(temp_data > 0, np.log(temp_data), 0)
    min_flux = np.min(log_data[conv_data > 0])
    max_flux = np.max(log_data)
    logger.info('min_flux for log_data: '+str(min_flux))
    logger.info('max_flux for log_data: '+str(max_flux))


#    del fits_in; del hdu
#    gc.collect()

    cdll.LoadLibrary("./RINGsize2.7_C.so")
    libc = CDLL("./RINGsize2.7_C.so")
    logger.info('# include the C-library:RINGsize2.7_C.so')
    
    # Give the directory name to shared C-library
    logger.info('# Give the directory name to shared C-library')
    enc_str1 = plain_file.encode('utf-8')
    enc_str2 = CP_file_name.encode('utf-8')
    enc_str3 = CP_Hist_name.encode('utf-8')
    #    enc_str = out_dir.encode('utf-8')
    file_name1 = create_string_buffer(enc_str1)
    file_name2 = create_string_buffer(enc_str2)
    file_name3 = create_string_buffer(enc_str3)
    libc.share_string(file_name1,file_name2,file_name3)
    
    # Give variables of parameters to shared C-library
    logger.info('# Give variables of parameters to shared C-library')
    libc.share_int(c_int(initial_center[0]), c_int(initial_center[1]), c_int(initial_radius), c_int(search_range))
    libc.share_float(c_float(cum_flux_lo), c_float(cum_flux_up))



######## tmask 2 ########
### Plot the image map (conv_data)
if tmask[0] <= 2 <= tmask[1]:
    print_start('2, Plot the image map')
    
    
    if x_axis == [0]:
        x_axis = [0,conv_data[1].size - 1]
    else:
        for i in range(len(x_axis)):
            x_axis[i] = int(x_axis[i])

    if y_axis == [0]:
        y_axis = [0,conv_data[1].size - 1]
    else:
        for i in range(len(y_axis)):
            y_axis[i] = int(y_axis[i])

    plot_image(conv_data,x_axis,y_axis, initial_image_file)
    logger.info('x_axis(channel) = ' + str(x_axis))
    logger.info('y_axis(channel) = ' + str(y_axis))
    del conv_data
    gc.collect()
    print_end('2')


######## tmask 3 ########
# Natural logarithmic flux data
#  - Calculate the radius for the initial center position,
#    and create the circle on the maser map.
#
if tmask[0] <= 3 <= tmask[1]:
    print_start('3, Natural logarithmic flux data. Create the circle on the map')
    
    task_no = 3
#    conv_data = raw_data * convert_f # Converted data array, unit: Jy/beam km/s

    # Set paratemers of initial center position and radius.
    logger.info('# Set paratemers of initial center position and radius.')
    logger.info('initial_center = ' + str(initial_center))
    logger.info('initial _radius = ' +str(initial_radius))
    center_x = initial_center[0];center_y = initial_center[1]
    min_flux = np.min(log_data[log_data > 0])
    max_flux = np.max(log_data)
    logger.info('min_flux for log_data: '+str(min_flux))
    logger.info('max_flux for log_data: '+str(max_flux))

    
    cal_radius_result = cal_radius(log_data,center_x,center_y,fig='Y')
    
    if not cal_radius_result:
        pass
    else:
        The_Radius = cal_radius_result[0]
        lo_radius = cal_radius_result[1]; up_radius = cal_radius_result[2]
        D1_R = cal_radius_result[3]; D2_R = cal_radius_result[4]
        D2_D1 = cal_radius_result[5]
        cal_flux_results = cal_radius_result[6]
        
        # Parameters for maser map
        ## Slice the image array for plot area
        maser_map = log_data[center_y - int(area_size_y/2):center_y + int(area_size_y/2) + 1, center_x - int(area_size_x/2):center_x + int(area_size_x/2) + 1].copy()
        # The map axes
        x_axis = [0,area_size_x]; y_axis = [0,area_size_x]
        
        # Parameters for the Normalized number of spots versus radius.
        x_radius = []; y_rate = []
        j = 0
        for j in range(len(cal_flux_results)):
            x_radius.append(cal_flux_results[j][2]); y_rate.append(cal_flux_results[j][5])
            j += 1
        
        # Plot the map (popup and save a file)
        IM = 'Center position = ('+str(center_x)+', '+str(center_y)+'), Radius = '+str(The_Radius[1])
        RR = 'Radius = '+str(The_Radius[1])+', D1 = '+str(D1_R)+', D2 = '+str(D2_R)+', D2/D1 = '+str(D2_D1)
        savefile = output_file+'_log_CP'+str(center_x)+'.'+str(center_y)+'_R'+str(The_Radius[1])+'_D'+str(round(D2_D1,2))+'_'+str(cum_flux_up)+'.png'
        plot_map_graph(maser_map,x_axis,y_axis,The_Radius[1], x_radius, y_rate, lo_radius[1], up_radius[1], IM, RR,savefile, pop='Y',cmap='Reds', cblb='', yl='Cumulative logarithmic flux')
    del maser_map; del log_data
    gc.collect()
    print_end('3')


######## tmask 4 ########
# Search the center position using the logarithmic flux data (without random noise)
#
if tmask[0] <= 4 <= tmask[1]:
    print_start('4, Search the center position using the logarithmic flux data (without random noise)')
    
    # Set paratemers of initial center position and radius.
    logger.info('# Set paratemers of initial center position and radius.')
    logger.info('initial_center = ' + str(initial_center))
    logger.info('initial _radius = ' +str(initial_radius))
    logger.info('nise in Image cube =' + str(noise))
    logger.info('iteration number =' + str(itn))
    center_x = initial_center[0]; center_y = initial_center[1]




    if os.path.exists(plain_file):
        os.remove(plain_file)



    # type of pointer of matrix
    _DOUBLE_PP = ndpointer(dtype=np.uintp, ndim=1, flags='C')
    # types of ctypes variables in function "add_matrix"
    libc.CalculateRadius.argtypes = [_DOUBLE_PP, c_int32, c_float]
    # type of function return
    libc.CalculateRadius.restype = None

    tp = np.uintp
    mpp = (log_data.__array_interface__['data'][0] + np.arange(log_data.shape[0])*log_data.strides[0]).astype(tp)

    flux_total = c_float(log_data.sum())
    task_no = c_int(4)
    col = c_int(len(log_data))

    libc.CalculateRadius(mpp, task_no, flux_total)

    

    del log_data
    gc.collect()
    print_end('4')


######## tmask 5 ########
# Search for the optimum center position and radius by giving random noise.
#
if tmask[0] <= 5 <= tmask[1]:
    print_start('5, Search for the optimum center position and radius by giving random noise.')
    
    # Set paratemers of initial center position and radius.
    logger.info('# Set paratemers of initial center position and radius.')
    logger.info('initial_center = ' + str(initial_center))
    logger.info('initial _radius = ' +str(initial_radius))
    logger.info('nise in Image cube =' + str(noise))
    logger.info('iteration number =' + str(itn))
    center_x = initial_center[0]; center_y = initial_center[1]
    

    # Make a valid_pixel data
    raw_X = raw_data.shape[1]; raw_Y = raw_data.shape[0]
    valid_pixel = np.where(raw_data > 0, 1, 0)
    
    
    if os.path.exists(CP_file_name):
        os.remove(CP_file_name)

    # Iteration: Randam Noise
    for count in range(itn):
        logger.info('--- count:'+str(count))
        normals = np.reshape(np.random.normal(0 ,noise, raw_X*raw_Y),(raw_Y,raw_X))
#        print 'normals:', normals.strides[0]
#        print 'valid_pixel:', valid_pixel.strides[0]

        temp_data = (normals * valid_pixel) + raw_data
#        print 'temp_data:', temp_data.strides[0]
        temp_data = np.where(temp_data > 0, temp_data, 0)
        temp_data = np.array(temp_data, dtype=np.float32)
#        print 'temp_data:', temp_data.strides[0]

        logger.info('raw data: mean='+str(raw_data.mean())+', std='+str(raw_data.std())+', max='+str(np.max(raw_data))+', min='+str(np.min(raw_data)))
        
        logger.info('normals: mean='+str(normals.mean())+', std='+str(normals.std())+', max='+str(np.max(normals))+', min='+str(np.min(normals)))
        logger.info('temp_data mean='+str(temp_data.mean())+', std='+str(temp_data.std())+', max='+str(np.max(temp_data))+', min='+str(np.min(temp_data)))
        
        
        # Make a Converted data array
        conv_data = temp_data * convert_f # Converted data array, unit: Jy/beam km/s
        #    conv_data = raw_data * convert_f # Converted data array, unit: Jy/beam km/s
        
        del temp_data
        gc.collect
        # Make a natural logarithmic flux data array from "comv_data".
        min_flux = np.min(conv_data[conv_data > 0])
        max_flux = np.max(conv_data)
        logger.info('min_flux for conv_data: '+str(min_flux))
        logger.info('max_flux for conv_data: '+str(max_flux))
        temp_data = conv_data / min_flux
        log_data = np.where(temp_data > 0, np.log(temp_data), 0)
        min_flux = np.min(log_data[conv_data > 0])
        max_flux = np.max(log_data)
        logger.info('min_flux for log_data: '+str(min_flux))
        logger.info('max_flux for log_data: '+str(max_flux))
        
        
        
        del conv_data; del temp_data
        gc.collect()
        
        
        # type of pointer of matrix
        _DOUBLE_PP = ndpointer(dtype=np.uintp, ndim=1, flags='C')
        # types of ctypes variables in function "add_matrix"
        libc.CalculateRadius.argtypes = [_DOUBLE_PP, c_int32, c_float]
        # type of function return
        libc.CalculateRadius.restype = None
        
        tp = np.uintp
        mpp = (log_data.__array_interface__['data'][0] + np.arange(log_data.shape[0])*log_data.strides[0]).astype(tp)
        
        flux_total = c_float(log_data.sum())
        task_no = c_int(5)
        
        libc.CalculateRadius(mpp, task_no, flux_total)



    del log_data
    gc.collect()
    print_end('5')


######## tmask 6 ########
# Plot the histogram of center position obtained in tmask 5.
#
if tmask[0] <= 6 <= tmask[1]:
    print_start('6, Plot the histogram of center position obtained in tmask 5.')
    
    if os.path.exists(CP_Hist_name):
        os.remove(CP_Hist_name)


    libc.HistCal()


    # Make a scatter plot
    CP_X = []; CP_Y =[]; CP_Count = []; CP_CrossSum =[];
    CP_Hist = open(CP_Hist_name, 'r')
    for line in CP_Hist:
        temp_list = line.split(',')
        CP_X.append(int(temp_list[0]))
        CP_Y.append(int(temp_list[1]))
        CP_Count.append(int(temp_list[2]))
        CP_CrossSum.append(int(temp_list[3]))
                             
    CP_Hist.close()

    # Calcurate the standard deviation
    CP_CrossSum_max = CP_CrossSum[0]
    CP_X_mu = CP_X[0]; CP_Y_mu = CP_Y[0]
    CP_X_min = CP_X[0]; CP_X_max = CP_X[0]
    CP_Y_min = CP_Y[0]; CP_Y_max = CP_Y[0]
    CP_Count_sum = 0
    for i in range(len(CP_X)):
        CP_Count_sum = CP_Count_sum + CP_Count[i]
        if (CP_CrossSum_max < CP_CrossSum[i]):
            CP_CrossSum_max = CP_CrossSum[i]
            CP_X_mu = CP_X[i]
            CP_Y_mu = CP_Y[i]
        if (CP_X_min > CP_X[i]):
            CP_X_min = CP_X[i]
        if (CP_X_max < CP_X[i]):
            CP_X_max = CP_X[i]
        if (CP_Y_min > CP_Y[i]):
            CP_Y_min = CP_Y[i]
        if (CP_Y_max < CP_Y[i]):
            CP_Y_max = CP_Y[i]


    temp_Hist = np.reshape(CP_Count,(CP_X_max - CP_X_min +1, CP_Y_max - CP_Y_min +1))
#    print temp_Hist

    Xj = CP_X_mu - CP_X_min
    Yk = CP_Y_mu - CP_Y_min
    logger.info('CP_Count_sum='+ str(CP_Count_sum))
    logger.info('CP_CrossSum_max=' + str(CP_CrossSum_max))
    logger.info('CP_X_mu=' + str(CP_X_mu) + ', CP_Y_mu=' + str(CP_Y_mu))
    logger.info('CP_X_min=' + str(CP_X_min) + ', CP_X_max=' + str(CP_X_max))
    logger.info('CP_Y_min=' + str(CP_Y_min) + ', CP_Y_max=' + str(CP_Y_max))
    logger.info('Xj=' + str(Xj) + ', Yk=' + str(Yk) + ', CP_Count=' + str(temp_Hist[Xj][Yk]))

    if os.path.exists(CP_result_file):
        os.remove(CP_result_file)
        
    CP_result = open(CP_result_file, 'w')
    CP_result.write('Center Position X=' + str(CP_X_mu) + '\n')
    CP_result.write('Center Position Y=' + str(CP_Y_mu) + '\n')



    error_pix = 1
    while (error_pix < 5):
        j = 0; k = 0; Count_IN = 0
        for j in range(CP_X_max - CP_X_min +1):
            for k in range(CP_Y_max - CP_Y_min +1):
                RR = (k - Yk)**2 + (j - Xj)**2
                if (RR <= (error_pix**2)):
                    Count_IN = Count_IN + temp_Hist[j][k]
        Hist_rate = float(Count_IN) / float(CP_Count_sum)
#        q_value = 1 - (1.0 - Hist_rate)/2 # Percent Point Function
#        ppf = norm.ppf(q=q_value, loc=0, scale=1)
        a,b = norm.interval(alpha=Hist_rate, loc=0, scale=1)
        if (error_pix == 1):
            b1 = str(round(b,2)) + ' Sigma'
        elif (error_pix == 2):
            b2 = str(round(b,2)) + ' Sigma'
        elif (error_pix == 3):
            b3 = str(round(b,2)) + ' Sigma'
        logger.info('error_pixel:' + str(error_pix) + ', Count_IN=' + str(Count_IN) + ', Count_IN_rate=' + str(Hist_rate) + ', CI=' + str(a) + ', '+ str(b))

        CP_list = 'error_pixel:' + str(error_pix) + ', Count_IN=' + str(Count_IN) + ', Count_IN_rate=' + str(Hist_rate) + ', CI=' + str(a) + ', '+ str(b) + '\n'
        CP_result.write(CP_list)
        
        error_pix += 1

    CP_result.close()


    # Plot the Histgram of Cenpter Position
    fig = plt.figure()  # fig,ax = plt.subplots()
    ax = plt.axes()
    #ax = fig.add_subplot(111)

    plt.scatter(CP_X, CP_Y, s=100, c=CP_Count, cmap='Blues')
    plt.xlabel('channel')
    plt.ylabel('channel')
    plt.title('Histgram of Center Position')

    plt.colorbar()

    plt.scatter(x=[CP_X_mu], y=[CP_Y_mu], c='black', s=100, marker='+', linewidth= 1.5)
    e1 = pat.Circle(xy = (CP_X_mu, CP_Y_mu), radius = 1, linewidth= 1.5,color = "black",label=b1,fill=False)
    e2 = pat.Circle(xy = (CP_X_mu, CP_Y_mu), radius = 2, linewidth= 1.5,color = "orange",label=b2,fill=False)
    e3 = pat.Circle(xy = (CP_X_mu, CP_Y_mu), radius = 3, linewidth= 1.5,color = "green",label=b3,fill=False)
    ax.add_patch(e1)
    ax.add_patch(e2)
    ax.add_patch(e3)
    plt.legend(fontsize=11)
    plt.axis('scaled')
    ax.set_aspect('equal')
    plt.savefig(savefile_Hist1,dpi=300)
    plt.show()
    plt.close()

    # Plot the Histgram of Cenpter Position: Sum of Cross-area
    plt.scatter(CP_X, CP_Y, s=100, c=CP_CrossSum, cmap='Reds')
    plt.xlabel('channel')
    plt.ylabel('channel')
    plt.title('Histgram of Cenpter Position: Sum of Cross-area')
    plt.colorbar()
    
    plt.savefig(savefile_Hist2,dpi=300)
    plt.close()
    
    print_end('6')


######## tmask 7 ########
# Plot the Final Map
if tmask[0] <= 7 <= tmask[1]:
    print_start('7, Plot the Final Map')
    
    task_no = 7

    # CP_X_mu, CP_Y_mu
    CP_X = []; CP_Y =[]; CP_Count = []; CP_CrossSum =[];
    CP_Hist = open(CP_Hist_name, 'r')
    for line in CP_Hist:
        temp_list = line.split(',')
        CP_X.append(int(temp_list[0]))
        CP_Y.append(int(temp_list[1]))
        CP_Count.append(int(temp_list[2]))
        CP_CrossSum.append(int(temp_list[3]))

    CP_Hist.close()

    CP_CrossSum_max = CP_CrossSum[0]
    CP_X_mu = CP_X[0]; CP_Y_mu = CP_Y[0]
    for i in range(len(CP_X)):
        if (CP_CrossSum_max < CP_CrossSum[i]):
            CP_CrossSum_max = CP_CrossSum[i]
            CP_X_mu = CP_X[i]
            CP_Y_mu = CP_Y[i]


    conv_data = raw_data * convert_f # Converted data array, unit: Jy/beam km/s
    
    final_result = cal_radius(conv_data,CP_X_mu,CP_Y_mu,fig='Y')


    The_Radius = final_result[0]
    lo_radius = final_result[1]; up_radius = final_result[2]
    D1_R = final_result[3]; D2_R = final_result[4]
    D2_D1 = final_result[5]
    cal_flux_results = final_result[6]
    
    # Parameters for maser map
    ## Slice the image array for plot area
    maser_map = conv_data[CP_Y_mu - int(area_size_y/2):CP_Y_mu + int(area_size_y/2) + 1, CP_X_mu - int(area_size_x/2):CP_X_mu + int(area_size_x/2) + 1].copy()
    # The map axes
    x_axis = [0,area_size_x]; y_axis = [0,area_size_y]
    
    # Parameters for the Normalized number of spots versus radius.
    x_radius = []; y_rate = []
    j = 0
    for j in range(len(cal_flux_results)):
        x_radius.append(cal_flux_results[j][2]); y_rate.append(cal_flux_results[j][5])
        j += 1
    
    # Plot the map (popup and save a file)
    IM = 'Center position = ('+str(CP_X_mu)+', '+str(CP_Y_mu)+'), Radius = '+str(The_Radius[1])
    RR = 'Radius = '+str(The_Radius[1])+', D1 = '+str(D1_R)+', D2 = '+str(D2_R)+', D2/D1 = '+str(D2_D1)
    savefile = output_file+'_CP'+str(CP_X_mu)+'.'+str(CP_Y_mu)+'_R'+str(The_Radius[1])+'.png'
    plot_map_graph(maser_map,x_axis,y_axis,The_Radius[1], x_radius, y_rate, lo_radius[1], up_radius[1], IM, RR,savefile, pop='Y',cmap='jet', cblb='Jy/bwam km/s', yl='Cumulative logarithmic flux')


    # Make a Final map
    x_axis_num = x_axis[1]/100
    y_axis_num = y_axis[1]/100
    x_axis_label = []; y_axis_label = []
    xticklabels = []; yticklabels = []
    if (x_axis_num % 2 == 1):
        for i in range(x_axis_num):
            x_axis_label.append(x_axis[1]/2 - 100*((x_axis_num/2)-i))
            xticklabels.append(0 + 10*((x_axis_num/2)-i))
        for i in range(y_axis_num):
            y_axis_label.append(y_axis[1]/2 - 100*((y_axis_num/2)-i))
            yticklabels.append(0 - 10*((y_axis_num/2)-i))
    if (x_axis_num % 2 == 0):
        for i in range(x_axis_num +1):
            x_axis_label.append((x_axis[1]/2) - 100*((x_axis_num/2)-i))
            xticklabels.append(0 + 10*((x_axis_num/2)-i))
        for i in range(y_axis_num +1):
            y_axis_label.append((y_axis[1]/2) - 100*((y_axis_num/2)-i))
            yticklabels.append(0 - 10*((y_axis_num/2)-i))


    fig = plt.figure(figsize=(10,10))  # fig,ax = plt.subplots()
    ax1 = fig.add_subplot(111)
    imp = ax1.imshow(maser_map,cmap='jet', norm=colors.LogNorm(vmin=0.1,vmax=cbmax))

    cbar = plt.colorbar(imp)
    plt.ylim(y_axis[0],y_axis[1])
    plt.xlim(x_axis[0],x_axis[1])
    ax1.set_xticks(x_axis_label)
    ax1.set_yticks(y_axis_label)
    ax1.set_xticklabels(xticklabels)
    ax1.set_yticklabels(yticklabels)
    ax1.tick_params(labelsize=16)
    plt.xlabel('R.A. Offset [mas]')
    plt.ylabel('Del. Offset [mas]')
    ax1.set_title(IM)
    cbar.set_label('Jy/beam km/s')
    plt.scatter(x=[x_axis[1]/2], y=[y_axis[1]/2], c='black', s=100, marker='+', linewidth= 1.5)
    c = pat.Circle(xy = (x_axis[1]/2, y_axis[1]/2), radius = The_Radius[1], linewidth= 1.5,color = "black",fill=False)
    ax1.add_patch(c)

    plt.savefig(io.BytesIO())

    ax1_pos = ax1.get_position()
    cax1_pos0 = cbar.ax.get_position()
    cax1_pos1 = [cax1_pos0.x0, ax1_pos.y0, cax1_pos0.x1 - cax1_pos0.x0, ax1_pos.y1 - ax1_pos.y0]
    cbar.ax.set_position(cax1_pos1)
    plt.rcParams["font.size"] = 16
    #    savefile = output_file+'_FinalMap_CP'+str(CP_X_mu)+'.'+str(CP_Y_mu)+'_R'+str(The_Radius[1])+'.png'
    savefile = output_file+'_FinalMap_CP'+str(CP_X_mu)+'.'+str(CP_Y_mu)+'_R'+str(The_Radius[1])+'.pdf'
    plt.savefig(savefile,dpi=500)
    plt.show()
    plt.close()
    
    del maser_map; del conv_data; del raw_data
    gc.collect()
    print_end('7')



######## tmask 8 ########
# Calcurate ring width
if tmask[0] <= 8 <= tmask[1]:
    print_start('8, Calcurate ring width')

    # Pick up Flux rate in "xxxx_Radius_Report.txt" file
    flux_rate = [];R_temp = []
    flux_cal_file = open(flux_cal_file_name, 'r')
    flux_cal_file.readline()
    for line in flux_cal_file:
        temp_list = line.split(',')
        flux_rate.append(float(temp_list[5]))
        R_temp.append(int(temp_list[2]))

    # Flux rate: 0.9
    flux_abs = abs(flux_rate[0] - 0.9)
    for i in range(len(flux_rate)):
        if flux_abs > abs(flux_rate[i] - 0.9):
            flux_abs = abs(flux_rate[i] - 0.9)
            R09 = R_temp[i]
    logger.info("R="+str(R09)+",: at flux rate 0.9")

    # Flux rate: 0.9
    flux_abs = abs(flux_rate[0] - 0.8)
    for i in range(len(flux_rate)):
        if flux_abs > abs(flux_rate[i] - 0.8):
            flux_abs = abs(flux_rate[i] - 0.8)
            R08 = R_temp[i]
    logger.info("R="+str(R08)+",: at flux rate 0.8")


    # Flux rate: 0.75
    flux_abs = abs(flux_rate[0] - 0.75)
    for i in range(len(flux_rate)):
        if flux_abs > abs(flux_rate[i] - 0.75):
            flux_abs = abs(flux_rate[i] - 0.75)
            R075 = R_temp[i]
    logger.info("R="+str(R075)+",: at flux rate 0.75")

    # Flux rate: 0.25
    flux_abs = abs(flux_rate[0] - 0.25)
    for i in range(len(flux_rate)):
        if flux_abs > abs(flux_rate[i] - 0.25):
            flux_abs = abs(flux_rate[i] - 0.25)
            R025 = R_temp[i]
    logger.info("R="+str(R025)+",: at flux rate 0.25")

    # Flux rate: 0.2
    flux_abs = abs(flux_rate[0] - 0.2)
    for i in range(len(flux_rate)):
        if flux_abs > abs(flux_rate[i] - 0.2):
            flux_abs = abs(flux_rate[i] - 0.2)
            R02 = R_temp[i]
    logger.info("R="+str(R02)+",: at flux rate 0.2")

    # Flux rate: 0.1
    flux_abs = abs(flux_rate[0] - 0.1)
    for i in range(len(flux_rate)):
        if flux_abs > abs(flux_rate[i] - 0.1):
            flux_abs = abs(flux_rate[i] - 0.1)
            R01 = R_temp[i]
    logger.info("R="+str(R01)+",: at flux rate 0.1")

    flux_cal_file.close()

    RingWidth_file = open(RingWidth, 'w')
    RingWidth_file.write("R(F0.9) = " + str(R09) +"\n")
    RingWidth_file.write("R(F0.8) = " + str(R08) +"\n")
    RingWidth_file.write("R(F0.75) = " + str(R075) +"\n")
    RingWidth_file.write("R(F0.25) = " + str(R025) +"\n")
    RingWidth_file.write("R(F0.2) = " + str(R02) +"\n")
    RingWidth_file.write("R(F0.1) = " + str(R01) +"\n")
    RingWidth_file.write("W(F0.1-0.9) = " + str(R09-R01) +"\n")
    RingWidth_file.write("W(F0.2-0.8) = " + str(R08-R02) +"\n")
    RingWidth_file.write("W(F0.25-0.75) = " + str(R075-R025) +"\n")
    RingWidth_file.close()
    print_end('8')



