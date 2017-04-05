'''
Contains functions that allows the visualization of a
system_design.BackscatterSystem instance. Can view a 3D model or a 
2D model in the XZ plane or YZ plane.
'''

import system_design as design
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib as mpl
import numpy as np


def system3D(system, title='Backscatter System', show_surface = True , 
                                                 show_fanbeam = True ,
                                                 show_Y_collimation = True,
                                                 show_Z_collimation = True ):

    '''
    Plot a 3D view of the entire system. Have the option to turn off certain
    parts of the system (object surface, Y and Z collimators).

    The image can be very slow to move around and observe if the Y collimators
    are plotted, since there are so many of them. 
    '''

    if isinstance(system, design.BackscatterSystem):

        pass

    else:

        print "ERROR: must pass an instance of" 

        print "system_design.BackscatterSystem"

    # =======================================================================
    # Create boundary limits for plot based on system dimensions
    # =======================================================================

    x,y,z = system.Detector.v_vector

    xx,yy,zz = system.Detector.v_vector + system.Detector.c_vector

    xmax  = x*1.1

    yspan = system.Detector.width*1.1

    zmin  = -system.Collimator.max_depth*1.1

    zmax  = zz*1.1

    limits = [[-1, xmax], [-yspan/2, yspan/2], [zmin, zmax]]

    # =======================================================================
    # Create figure instance and start adding the system components to it
    # =======================================================================

    fig = create_figure_3D(title, limits)

    # Add to flat object surface being imaged

    if show_surface:

        create_object_surface(fig, limits)

    # Add the fins of the collimation grid

    if show_Z_collimation:

        for fin in system.Collimator.z_fins:

            create_fin(fig, fin)

    if show_Y_collimation:

        for fin in system.Collimator.y_fins:

            create_fin(fig, fin)

    # Add X-ray fan beam

    if show_fanbeam:

        create_fan(fig, system)

    # Add Detector

    create_detector(fig, system.Detector)

    plt.show()

    return fig

def xray_spectrum(xray_object, norm = True, title = 'X-ray Spectrum'):
  '''
  Plots the X ray spectrum from an XRayFanBeam instance

  Will only work if the XRayFanBeam instance actually has a spectrum to plot
  and is not a mono-energetic source.

  Can choose to plot either the original spectrum or the normalized spectrum,
  default is the normalized spectrum.
  '''

  if isinstance(xray_object, design.BackscatterSystem):

    if isinstance(xray_object.XRayFanBeam, design.XRayFanBeam):

      if isinstance(xray_object.XRayFanBeam.energy, np.ndarray):

        if norm:

          spectrum = xray_object.XRayFanBeam.normalized_energy

        else:

          spectrum = xray_object.XRayFanBeam.energy

      else:

        print "\nX-ray energy source is not a spectrum; nothing to plot"

        print "XRayFanBeam.energy needs to be a numpy.ndarray"

        return

    else:

        print "\nBackscatterSystem instance does not contain an"

        print "XRayFanBeam instance"

        return
 
  elif isinstance(xray_object, design.XRayFanBeam):

    if isinstance(xray_object.energy, np.ndarray):

        if norm:

          spectrum = xray_object.normalized_energy

        else:

          spectrum = xray_object.energy

    else:

      print "\nX-ray energy source is not a spectrum; nothing to plot"

      print "XRayFanBeam.energy needs to be a numpy.ndarray"

      return

  else:

    print "\nError: Did not pass an instance of XRayFanBeam or an"

    print "instance of BackscatterSystem containing an XRayFanBeam."

    return

  figure = plt.figure()

  ax = figure.add_subplot(111)

  ax.set_title(title)

  ax.set_xlabel('Energy (MeV)')

  if norm:

    ax.set_ylabel('Probability')

  else:

    ax.set_ylabel('Frequency')

  plt.plot(spectrum[:,0], spectrum[:,1])

  plt.show()

  return figure, ax

def get_vertices(fin):

    corner1 = fin.v_vector

    corner2 = fin.v_vector + fin.a_vector + fin.b_vector + fin.c_vector

    side1 = np.vstack((corner1, 
                       corner1 + fin.b_vector,
                       corner1 + fin.b_vector + fin.c_vector,
                       corner1 + fin.c_vector ))

    side2 = np.vstack((corner1, 
                       corner1 + fin.a_vector,
                       corner1 + fin.a_vector + fin.b_vector,
                       corner1 + fin.b_vector ))

    side3 = np.vstack((corner1, 
                       corner1 + fin.a_vector,
                       corner1 + fin.a_vector + fin.c_vector,
                       corner1 + fin.c_vector ))

    side4 = np.vstack((corner2, 
                       corner2 - fin.b_vector,
                       corner2 - fin.b_vector - fin.c_vector,
                       corner2 - fin.c_vector ))

    side5 = np.vstack((corner2, 
                       corner2 - fin.a_vector,
                       corner2 - fin.a_vector - fin.b_vector,
                       corner2 - fin.b_vector ))

    side6 = np.vstack((corner2, 
                       corner2 - fin.a_vector,
                       corner2 - fin.a_vector - fin.c_vector,
                       corner2 - fin.c_vector ))
    
    vertices_list = [side1, side2, side3, side4, side5, side6] 

    return vertices_list

def create_figure_3D(title, limits):
    '''
    Creates a 3D figure instance 

    plotting limits are determiend by the maximum dimensions of the system 
    '''

    figure = plt.figure()

    ax = figure.add_subplot(111, projection='3d')

    xlims, ylims, zlims = limits

    ax.set_xlim( xlims )

    ax.set_ylim( ylims )

    ax.set_zlim( zlims )

    ax.set_title(title)

    return figure

def create_figure_2D(title, limits):
    '''
    Creates a 2D figure instance 

    plotting limits are determiend by the maximum dimensions of the system 
    '''

    figure = plt.figure()

    ax = figure.add_subplot(111)

    xlims, ylims = limits

    ax.set_xlim( xlims )

    ax.set_ylim( ylims )

    ax.set_title(title)

    return figure, ax

def create_object_surface(figure, limits):
    
    xlims, ylims, zlims = limits

    xmin, xmax = xlims

    ymin, ymax = ylims

    X = [xmin, xmin, xmax, xmax]

    Y = [ymin, ymax, ymax, ymin]

    Z = [0]*4

    vertices = zip(X, Y, Z)

    add_surface(figure, vertices, color='b', alpha=1.0)

def create_fin(figure, fin):
    
    vertices_list = get_vertices(fin)

    for vertices in vertices_list:

        add_surface(figure, vertices, color='k', alpha=.8)

def create_fan(figure, system):

    fan_angle_radians = system.XRayFanBeam.fan_angle * np.pi/180.

    xray_path = system.XRayFanBeam.surface_offset + system.Collimator.max_depth

    source = [0, 0, system.XRayFanBeam.surface_offset]

    bottom_neg_y = [ 0, 
                    -(xray_path)*np.tan(fan_angle_radians),
                    -system.Collimator.max_depth ]

    bottom_pos_y = [ 0, 
                     (xray_path)*np.tan(fan_angle_radians),
                    -system.Collimator.max_depth ]

    vertices = [source, bottom_neg_y, bottom_pos_y]

    add_surface(figure, vertices, color='r', alpha=0.7)

def create_detector(figure, detector):
       
    vertices_list = get_vertices(detector)

    for vertices in vertices_list:

        add_surface(figure, vertices, color='#22F0E6', alpha=1.0)

def add_surface(figure, vertices, color, alpha):

    surface = Poly3DCollection([vertices], facecolors = color, 
                                           edgecolors=color  ,
                                           alpha=alpha         )

    figure.gca().add_collection(surface)

def xz_system(system, title = 'Depth Collimation', show_ray_limits = True):
    '''
    Plot a 2D view of just the Z collimators. 

    Have the option to view rays being traced from the detector pixels
    to their respective maximum depths of view for each collimator

    Due to the limitations of plotting rectangle patches in matplotlib, 
    the view is from the angle of the detector, so everything is rotated 
    by the degree that the detector was set up at. 
    '''

    rectangle_patches = []

    lines = []

    angle = system.Detector.angle_radians

    rotate_matrix = np.array([[ np.sin(angle), np.cos(angle)],
                              [-np.cos(angle), np.sin(angle)] ])

    x,y,z = system.Detector.v_vector + system.Detector.a_vector

    xx,yy,zz = system.Detector.v_vector + system.Detector.c_vector + \
                                          system.Detector.a_vector

    x_r = np.dot(rotate_matrix, np.array([[-1. , x*1.1, x*1.1],
                                          [0., 0., z*1.1 ]
                                         ])
                )

    z_r = np.dot(rotate_matrix, np.array([[0. , 0., xx*1.1],
                                          [-system.Collimator.max_depth*1.1,
                                          system.XRayFanBeam.surface_offset*1.1,
                                          zz*1.1]
                                         ])
                 )

    xlims = [np.min([x_r[0], z_r[0]]), np.max([x_r[0], z_r[0]])]

    zlims = [np.min([x_r[1], z_r[1]]), np.max([x_r[1], z_r[1]])]

    xlims = [round(x) for x in xlims]

    zlims = [round(z) for z in zlims]

    if zlims[1] - zlims[0] > xlims[1] - xlims[0]:

        limits = [[xlims[0], xlims[0] + (zlims[1] - zlims[0])], zlims] 

    else:

        limits = [xlims, [zlims[0], zlims[0] + (xlims[1] - xlims[0])]]

    fig, ax = create_figure_2D(title, limits)

    # Detector

    det_x_start, det_y_start  = np.dot(rotate_matrix, 
                              [system.Detector.v_vector[0], 
                               system.Detector.v_vector[2]] 
                              )

    det_patch = patches.Rectangle( ( det_x_start , det_y_start )        , 
                                       system.Detector.thickness , 
                                       system.Detector.height    , 
                                       fc='#22F0E6'       ,   
                                       alpha=1.0          , 
                                  )

    rectangle_patches.append(det_patch)

    # Fins 

    for fin in system.Collimator.z_fins:

        f_x, f_y  = np.dot(rotate_matrix, [fin.v_vector[0], fin.v_vector[2]])

        f = patches.Rectangle((f_x, f_y), -fin.length, fin.thickness, 
                              fc='k', alpha=1.0 )

        rectangle_patches.append(f)

    # Fan Beam

    fan_x_start, fan_y_start = np.dot(rotate_matrix, [0 , 
                                      system.XRayFanBeam.surface_offset])

    fan_x_finish, fan_y_finish = np.dot(rotate_matrix, [0 , 
                                       round(-system.Collimator.max_depth*1.1)])
   
    lines.append( plt.Line2D( (fan_x_start, fan_x_finish), 
                              (fan_y_start, fan_y_finish), 
                              color='r', lw=1.5) )

    # Object Surface

    surf_x_start, surf_y_start = np.dot(rotate_matrix, [-1 , 0]) 


    surf_x_finish, surf_y_finish = np.dot(rotate_matrix,
                                 [round(system.Detector.v_vector[0]), 0] )
   
    lines.append( plt.Line2D( (surf_x_start, surf_x_finish), 
                              (surf_y_start, surf_y_finish), 
                              color='b', lw=1.5) )

    # If desired, create Line patches and add to axes

    if show_ray_limits:
 
        for i in range(system.Detector.num_pixel_H):

            depth = system.Collimator.depths[i]

            pixel = system.Collimator.starting_point + \
                    (i + 1)*np.array([ -system.Detector.pixel_height * \
                                       np.cos(system.Detector.angle_radians) ,
                                       0 ,
                                       system.Detector.pixel_height * \
                                       np.sin(system.Detector.angle_radians)
                                       ])

            ray_x_start, ray_y_start = np.dot(rotate_matrix, [pixel[0], 
                                                              pixel[2]])

            ray_x_finish, ray_y_finish = np.dot(rotate_matrix, 
                                                [0, -depth]    )

            lines.append( plt.Line2D( (ray_x_start, ray_x_finish), 
                                      (ray_y_start, ray_y_finish), 
                                       color='r') )

    # Add to axes

    for patch in rectangle_patches:

        ax.add_patch(patch)

    for line in lines:

        ax.add_line(line)


    ax.set_aspect('equal')

    plt.show()

    return fig, ax

def yz_system(system, row = 0, title = 'Surface Collimation', 
                                       show_ray_limits = True):
    '''
    Plot a 2D view of just the Y collimators, at a specific row in the
    detector.

    Have the option to view rays being traced from the detector pixels
    to their respective maximum depths of view for each collimator
    '''

    rectangle_patches = []

    lines = []

    yspan = (system.Detector.width + system.Collimator.collimator_thickness +\
             system.Collimator.y_resolution)


    for fin in system.Collimator.y_fins:

        if fin.index_number <= (row + 1)*(system.Detector.num_pixel_W + 1) \
              and fin.index_number > row*(system.Detector.num_pixel_W + 1):

            fin_length = fin.length

            break

    zspan = (system.Collimator.paths[row] + fin_length)*1.5

    limits = [ [-round(yspan/2.), round(yspan/2.)], 
               [-round(zspan/2.), round(zspan/2.)]  ]

    fig, ax = create_figure_2D(title, limits)

    # Detector

    det_y_start, det_z_start  = (-system.Detector.width/2., fin_length)

    det_patch = patches.Rectangle( ( det_y_start , det_z_start ) , 
                                       system.Detector.width , 
                                       system.Detector.thickness , 
                                       fc='#22F0E6'       ,   
                                       alpha=1.0          , 
                                  )

    rectangle_patches.append(det_patch)

    # Fins 

    for fin in system.Collimator.y_fins:

        if fin.index_number <= (row + 1)*(system.Detector.num_pixel_W + 1) \
              and fin.index_number > row*(system.Detector.num_pixel_W + 1):

            f_y = -fin.v_vector[1]

            f_z = fin.length

            f = patches.Rectangle((f_y, f_z), -fin.thickness, -fin.length, 
                              fc='k', alpha=1.0 )

            rectangle_patches.append(f)

            # If desired, create Line patches and add to axes

            if show_ray_limits and \
               fin.index_number < (row + 1)*(system.Detector.num_pixel_W + 1):

                y1_start, z1_start = ( f_y - fin.thickness, f_z )

                y1_finish  = f_y - (system.Collimator.y_resolution + \
                                    system.Detector.pixel_width + \
                                    system.Collimator.collimator_thickness)/2.

                z1_finish = -system.Collimator.paths[row]
    
                y2_start, z2_start = ( f_y - system.Detector.pixel_width, f_z )

                y2_finish, z2_finish = ( y1_finish + \
                                         system.Collimator.y_resolution,
                                         -system.Collimator.paths[row] )

                lines.append( plt.Line2D( (y1_start, y1_finish), 
                                          (z1_start, z1_finish), 
                                           color='r') )

                lines.append( plt.Line2D( (y2_start, y2_finish), 
                                          (z2_start, z2_finish), 
                                           color='r') )

    # Object Surface

    surf_y_start, surf_y_finish = limits[0]


    surf_z_start, surf_z_finish = 2*[-system.Collimator.paths[row]]
   
    lines.append( plt.Line2D( (surf_y_start, surf_y_finish), 
                              (surf_z_start, surf_z_finish), 
                              color='b', lw=1.5) )

    # Add to axes

    for patch in rectangle_patches:

        ax.add_patch(patch)

    for line in lines:

        ax.add_line(line)


    ax.set_aspect('equal')

    plt.show()

    return fig, ax

def plot_system_parameters(system, variable, values, fs = 20):
    '''
    Pass it a fully designed system and tell it which parameter to vary, along
        with a range of values to vary it by

    variable : Type string, attribute to vary inside system.

    sugested variable options: "num_pixel_H", "num_pixel_W", "pixel_height",        
                          "pixel_width", "angle"      , "collimator_thickness",  
                          "max_depth"  , "surface_offset" , "beam_offset",          
                          "y_resolution" , "length"

    values : (Type: array or list) range of values for variable.
    '''

    depths = []

    if variable in system.Detector.__dict__:

        parameter = system.Detector

    elif variable in system.Collimator.__dict__:

        parameter = system.Collimator

    elif variable in system.XRayFanBeam.__dict__:

        parameter = system.XRayFanBeam

    else:

        print "\nRequested variable does not exist in System\n"

        return

    for value in values:

        if variable == "angle":

            parameter.__setattr__(variable, value)

            parameter.__setattr__("angle_radians", value*np.pi/180.)

        elif variable == "angle_radians":

            parameter.__setattr__(variable, value)

            parameter.__setattr__("angle", value*180./np.pi)

        else:

            parameter.__setattr__(variable, value)          

        system.optimize_fan_beam()

        system.optimize_collimator()

        system.system_check()

        depths.append(system.depths)

        #print system.Collimator.max_z_length


    var_title = ' '.join(variable.split('_')).title()

    fig1, ax1 = plt.subplots()

    fig2, ax2 = plt.subplots()

    ax1.set_title('System View vs. ' + var_title, fontsize=1.5*fs)

    ax1.set_xlabel(var_title, fontsize=fs)

    ax1.set_ylabel('Depth (cm)', fontsize=fs)

    pix_centers = [depth["pix_center"] for depth in depths]

    full_pix_range = [depth["full_pix_range"] for depth in depths]

    pix_range = [depth["range"] for depth in depths]

    for i,(pixels,fpranges,ranges) in enumerate(zip(pix_centers,full_pix_range,
                                                                   pix_range)):

        #ax2.scatter([values[i]]*len(pixels), pixels ,
         #          facecolor='b', label='Pixel Center')

        lower_depth_fpr = [fpr[0] - pix for pix,fpr in zip(pixels,fpranges)]
        upper_depth_fpr = [pix - fpr[1] for pix,fpr in zip(pixels,fpranges)]

        ax1.errorbar([values[i]]*len(pixels), pixels, 
                         yerr=[lower_depth_fpr, upper_depth_fpr],
                         fmt='ob',ecolor='g', capthick=2)

        lower_depth_r = [r[0] - pix for pix,r in zip(pixels,ranges)]
        upper_depth_r = [pix - r[1] for pix,r in zip(pixels,ranges)]

        ax2.errorbar([values[i]]*len(pixels), pixels, 
                         yerr=[lower_depth_r,upper_depth_r],
                         fmt='ob',ecolor='g', capthick=2)

    



    #preferred_view_depths = [i*system.Detector.pixel_height + 
     #                        system.Detector.pixel_height/2 for i in \
      #                       range(system.Detector.num_pixel_H)]




    plt.show()


    return depths

    '''
    print "\n\nPixel centers"
    print system.depths["pix_center"]
    print "\nfull pixel ranges"
    print system.depths["full_pix_range"]
    print "\nPixel ranges"
    print system.depths["range"]
    print"\ny resolutions"
    print system.depths["y_res"]
    '''


