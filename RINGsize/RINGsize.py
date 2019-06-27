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
    return control

def chekin(control):
    global out_dir, output_file_name, pipelog, output_file
    global data_array_binary, initial_image_file,results_report_name,flux_cal_file_name
    global tmask ,image_file_name, x_axis, y_axis
    global initial_center, initial_radius, area_size_x, area_size_y
    global cum_flux_up, cum_flux_lo, search_range
    out_dir = control.get('out_dir', [False])[0]
    output_file_name = control.get('output_file_name', [False])[0]
    output_file = out_dir + '/' + output_file_name
    pipelog = output_file + '.pipelog'
    data_array_binary = output_file + '_conv_data.array'
    initial_image_file = output_file + '_initial_image.png'
    results_report_name = output_file + '_results_report.txt'
    flux_cal_file_name = output_file + '_flux_cal_file.txt'
    
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

def plot_map_graph(im_data_array,x_axis,y_axis,radius,x,y,radius_lo,radius_up,im_title,xy_title,savefile, pop):
    fig = plt.figure(figsize=(15,5))  # fig,ax = plt.subplots()
    ax1 = fig.add_subplot(121)
    imp = ax1.imshow(im_data_array,cmap='jet', norm=colors.LogNorm())
    cbar = plt.colorbar(imp)
    plt.ylim(y_axis[0],y_axis[1])
    plt.xlim(x_axis[0],x_axis[1])
    plt.xlabel('channel')
    plt.ylabel('channel')
    ax1.set_title(im_title)
    cbar.set_label('Jy/beam km/s')
    plt.scatter(x=[x_axis[1]/2], y=[y_axis[1]/2], c='black', s=70, marker='+')
    c = pat.Circle(xy = (x_axis[1]/2, y_axis[1]/2), radius = radius, color = "black",fill=False)
    ax1.add_patch(c)
    
    ax2 = fig.add_subplot(122)
    ax2.scatter(x,y)
    ax2.set_title(xy_title)
    ax2.set_xlabel('radius (ch)')
    ax2.set_ylabel('Normalised cumulative flux')
    #ax2.axhline(0.5)
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
        #if os.path.exists(flux_cal_file_name):
        #os.remove(flux_cal_file_name)
    
    #flux_cal_file = open(flux_cal_file_name, 'a+')
    #flux_cal_file.write('# center_x(ch), center_y(ch), radius(ch),total_flux, flux_in_circle, rate, abs(rate-0.5) \n')


    i = 0
    while (initial_radius -1 > i):
        radius = initial_radius - i
        # Calculate the flux in the circle
        cal_flux_list = cal_flux(conv_data,center_x,center_y,radius)
        cal_flux_results.append(cal_flux_list)
        #str_cal_flux_list = str(cal_flux_list[0])+','+str(cal_flux_list[1])+','+str(cal_flux_list[2])+','+str(cal_flux_list[3])+','+str(cal_flux_list[4])+','+str(cal_flux_list[5])+','+str(cal_flux_list[6])
        #flux_cal_file.write(str_cal_flux_list+'\n')
        i += 1

    #flux_cal_file.close()
    
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
    i = 0; temp_list = []
    for i in range(len(list_name)):
        if list_name[i][require_j] >= require_val:
            temp_list.append(list_name[i])
        i += 1
    if not temp_list:
        return 0
    else:
        lower_limit_i = temp_list.index(temp_list[-1])
        out_val = list_name[lower_limit_i][out_j]
        return lower_limit_i, out_val

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
    #logger.error('Insufficient arguments')
#control = parse_inp('r18316c_v2.txt')
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

    # Convert the unit of raw data to "Jy/beam km/s".
    logger.info('# Convert the unit of raw data to "Jy/beam km/s')
    freq_if = fits_header.get('CRVAL3', [False])
    ch_sep = fits_header.get('CDELT3', [False])
    logger.info('freq_if = ' + str(freq_if) + ', ch_sep = ' + str(ch_sep))

    raw_data=fits_data[0,0] # raw_data array
    #np.savetxt('./raw_data.txt', raw_data)
    convert_f = ch_sep * c * 1e-3 / freq_if
    logger.info('convert function = '+str(convert_f))
    conv_data = raw_data * convert_f # Converted data array, unit: Jy/beam km/s
    # Save the converted data array in a binary file
    logger.info('# Save the converted data array in a binary file')
    logger.info('file: '+data_array_binary +'.npy')
    np.save(data_array_binary, conv_data)
    print_end('1')
#######################################
# Laod the converted data array when the tmask 1 is skipped.
if tmask[0] != 1:
    logger.info('# Laod the converted data array when the tmask 1 is skipped.')
    data_array_binary_in = data_array_binary + '.npy'
    logger.info(data_array_binary_in)
    conv_data = np.load(data_array_binary_in)
    logger.info('conv_data.shape = '+str(conv_data.shape))


######## tmask 2 ########
### Plot the image map
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
    print_end('2')


######## tmask 3 ########
# Calculate the radius for the initial center position,
# and create the circle on the maser map.
if tmask[0] <= 3 <= tmask[1]:
    print_start('3, Calculate the radius for the initial center position')

    # Set paratemers of initial center position and radius.
    logger.info('# Set paratemers of initial center position and radius.')
    logger.info('initial_center = ' + str(initial_center))
    logger.info('initial _radius = ' +str(initial_radius))
    center_x = initial_center[0];center_y = initial_center[1]
    
    cal_radius_result = cal_radius(conv_data,center_x,center_y,fig='Y')
        
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
        maser_map = conv_data[center_y - int(area_size_y/2):center_y + int(area_size_y/2) + 1, center_x - int(area_size_x/2):center_x + int(area_size_x/2) + 1].copy()
        # The map axes
        x_axis = [0,area_size_x]; y_axis = [0,area_size_x]
            
        # Parameters for the normalised cumulative flux versus radius.
        x_radius = []; y_rate = []
        j = 0
        for j in range(len(cal_flux_results)):
            x_radius.append(cal_flux_results[j][2]); y_rate.append(cal_flux_results[j][5])
            j += 1

        # Plot the map (popup and save a file)
        IM = 'Center position = ('+str(center_x)+', '+str(center_y)+'), Radius = '+str(The_Radius[1])
        RR = 'Radius = '+str(The_Radius[1])+', D1 = '+str(D1_R)+', D2 = '+str(D2_R)+', D2/D1 = '+str(D2_D1)
        savefile = output_file+'_CP'+str(center_x)+'.'+str(center_y)+'_R'+str(The_Radius[1])+'_D'+str(round(D2_D1,2))+'_'+str(cum_flux_up)+'.png'
        plot_map_graph(maser_map,x_axis,y_axis,The_Radius[1], x_radius, y_rate, lo_radius[1], up_radius[1], IM, RR,savefile, pop='Y')
    print_end('3')



######## tmask 4 ########
# Search the optimum center position and radius
# 
if tmask[0] <= 4 <= tmask[1]:
    print_start('4, Search the optimum center position and radius')

    # Set paratemers of initial center position and radius.
    logger.info('# Set paratemers of initial center position and radius.')
    logger.info('initial_center = ' + str(initial_center))
    logger.info('initial _radius = ' +str(initial_radius))
    center_x = initial_center[0];center_y = initial_center[1]
    radius_file_name = output_file + '_center_radius.txt'
    
    if os.path.exists(radius_file_name):
        os.remove(radius_file_name)

    radius_file = open(radius_file_name, 'a+')
    radius_file.write('# center_x(ch), center_y(ch), radius(ch), D1_R, D2_R, D2_R/D1_R \n')

    if os.path.exists(results_report_name):
        os.remove(results_report_name)
    
    report_file = open(results_report_name, 'a+')
    report_file.write('# center_x(ch), center_y(ch), radius(ch), D1_R, D2_R, D2_R/D1_R \n')

    center_x_bg = initial_center[0] - int(search_range/2)
    center_x_en = initial_center[0] + int(search_range/2)
    center_y_bg = initial_center[1] - int(search_range/2)
    center_y_en = initial_center[1] + int(search_range/2)
    search_x = list(range(center_x_bg, center_x_en + 1))
    search_y = list(range(center_y_bg, center_y_en + 1))

    center_radius = []
    x = 0; y = 0
    for y in range(len(search_y)):
        for x in range(len(search_x)):
            center_x = search_x[x]; center_y = search_y[y]
            logger.info('searching matrix = ('+str(x)+', '+str(y)+'), center_xy = (' +str(center_x)+', '+str(center_y)+')' )
            
            cal_radius_result = cal_radius(conv_data,center_x,center_y,fig='N')
            if not cal_radius_result:
                pass
            else:
                The_Radius = cal_radius_result[0]
                lo_radius = cal_radius_result[1]; up_radius = cal_radius_result[2]
                D1_R = cal_radius_result[3]; D2_R = cal_radius_result[4]
                D2_D1 = cal_radius_result[5]
                temp_list = []

                temp_list = [center_x, center_y,The_Radius[1],D1_R,D2_R,D2_D1]
                center_radius.append(temp_list)
                str_center_radius = str(center_x)+','+str(center_y)+','+str(The_Radius[1])+','+str(D1_R)+','+str(D2_R)+','+str(D2_D1)
                radius_file.write(str_center_radius + '\n')

            x += 1
        y += 1

    D1_min = list_minimum(center_radius, 3, 3)
    logger.info('*************************')
    logger.info('Results report')
    for i in D1_min[0]:
        logger.info('center_xy = (' + str(center_radius[i][0]) + ', '+str(center_radius[i][1])+'), The radius = '+ str(center_radius[i][2])+ ', D1 = '+str(center_radius[i][3])+', D2 = '+str(center_radius[i][4])+', D2_D1 = '+str(center_radius[i][4]))
    
        str_report_file = str(center_radius[i][0]) + ','+str(center_radius[i][1])+','+ str(center_radius[i][2])+','+str(center_radius[i][3])+','+str(center_radius[i][4])+','+str(center_radius[i][5])
        report_file.write(str_report_file + '\n')

    radius_file.close()
    report_file.close()
    print_end('4')


######## tmask 5 ########
# Plot the final map
#
if tmask[0] <= 5 <= tmask[1]:
    print_start('5, Plot the final maps')

    # read the center position
    center_xy = []
    report_file = open(results_report_name, 'r')
    for line in report_file:
        if line[0] == '#':
            continue
        temp_list = line.split(',')
        center_xy.append([int(temp_list[0]),int(temp_list[1])])
    report_file.close()

    k = 0
    for k in range(len(center_xy)):
        # Set paratemers of center position and radius.
        center_x = center_xy[k][0]; center_y = center_xy[k][1]
        logger.info('# Set paratemers of center position and radius.')
        logger.info('center = (' + str(center_x)+', '+str(center_y)+')')
 
        cal_radius_result = cal_radius(conv_data,center_x,center_y,fig='Y')

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
            maser_map = conv_data[center_y - int(area_size_y/2):center_y + int(area_size_y/2) + 1, center_x - int(area_size_x/2):center_x + int(area_size_x/2) + 1].copy()
            # The map axes
            x_axis = [0,area_size_x]; y_axis = [0,area_size_x]
                
            # Parameters for the normalised cumulative flux versus radius.
            x_radius = []; y_rate = []
            j = 0
            for j in range(len(cal_flux_results)):
                x_radius.append(cal_flux_results[j][2]); y_rate.append(cal_flux_results[j][5])
                j += 1

            # Plot the map (popup and save a file)
            IM = 'Center position = ('+str(center_x)+', '+str(center_y)+'), Radius = '+str(The_Radius[1])
            RR = 'Radius = '+str(The_Radius[1])+', D1 = '+str(D1_R)+', D2 = '+str(D2_R)+', D2/D1 = '+str(D2_D1)
            savefile = output_file+'_CP'+str(center_x)+'.'+str(center_y)+'_R'+str(The_Radius[1])+'_D'+str(round(D2_D1,2))+'_'+str(cum_flux_up)+'.png'
            plot_map_graph(maser_map,x_axis,y_axis,The_Radius[1], x_radius, y_rate, lo_radius[1], up_radius[1], IM, RR,savefile, pop='N')

    print_end('5')

