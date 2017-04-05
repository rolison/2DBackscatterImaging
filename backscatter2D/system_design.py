'''
Contains all the classes and functions necessary to develop an
X-ray Backscatter System that utilizes a fan beam X-ray source,
a 2D pixelated detector array, and a collimation grid that rests on the face
of the detector array. This allows for 2D image acquisition for a single 
position of the backscatter system. By scanning the system in one dimension,
it's possible to develop a 3D image of an object
'''

import numpy as np

class BackscatterSystem(object):
    '''
    Object that is a composition of all cmoponents in the backscatter system:
        1) X-Ray Fan Beam source
        2) 2D Pixelated detector array
        3) Collimation grid for detector array

    Parameters Needed from User:
    ============================

    step_size          : distance taken by each step in the X dimension as
                         the system scans over an object (cm)

    Default Values:
    ===============
    
    step_size          : None

    Parameters:
    ===========

    depths             : A dictionary containing useful information on the
                         system. Gets created after calling the system_check
                         method. Each Key conatinsa list that's as long as the
                         number of pixels in the height of the detector.

                         Keys:

                         pix_center     : depth position (cm) determined by
                                          tracing a line from the center of a
                                          pixel to where it intersects with the
                                          fan beam. This line is orthogonal to
                                          the surface of the pixel. A negative 
                                          value indicates it is above the 
                                          surface of the object being imaged.

                         full_pix_range : tuple (cm), containing the 
                                          (deepest depth, shallowest depth) 
                                          that the pixel has full sight of. 
                                          Full sight meaning any photons that
                                          scatter in this range can pass
                                          through any part of the pixel's face.
                                          A negative value indicates it is
                                          above the surface of the object
                                          being imaged.

                         range          : tuples (cm), containing the 
                                          (deepest depth, shallowest depth)
                                          that the pixel can possibly see. So
                                          this range is always greater than the
                                          full_pix_range. It includes areas of
                                          the fan beam where photon scatters 
                                          are only partially visible to the
                                          pixel face. A negative value 
                                          indicates it is above the surface of
                                          the object being imaged.

                         y_res          : resolution in the y dimension (cm).
                                          Determined by using the path length
                                          of a photon scattering from the fan
                                          beam and going straight to the
                                          center of the pixel. This path length
                                          is the length of the line traced when 
                                          calculating pix_center


    You have the option to fully build a backscatter system from this object
    or you can build as many components as you want separately and then add
    them to an instance of this object

    Examples:
    =========

        1) Build from scratch with this object:
                
            import system_design as sysd
          
            your_system = sysd.BackscatterSystem()

        2) Build a detector separately and add it to your system, and build
           the remaining components during instantiation of your system:

            import system_design as sysd

            your_detector = sysd.Detector(your_settings)

            your_system = sysd.BackscatterSystem(your_detector)

        3) Build all components beforehand and then put in a system

            import system_design as sysd

            your_xray = sysd.XRayFanBeam(your_xray_settings)

            your_detector = sysd.Detector(your_detector_settings)

            your_collimator = sysd.UniformCollimator(your_collimator_settings)

            your_system = sysd.BackscatterSystem(your_xray, your_detector,
                                                your_collimator)


    '''
    def __init__(self, beam_input=None, det_input=None, collimator_input=None,
                                                              step_size = None):

        self.step_size = step_size

        # ====================================================================
        # Load or set up X-ray Fan Beam
        # ====================================================================


        if isinstance(beam_input, XRayFanBeam):

            self.XRayFanBeam = beam_input

        else:

            responded = False

            print "\nX-Ray Fan Beam source doesn't exist, would you like to" +\
                  " build one?"

            print "Options: yes, no, use defaults \n"

            response = raw_input('>>  ')

            while not responded:

                if response == 'y' or response == 'yes':

                    self.XRayFanBeam = XRayFanBeam(build_manually = True)

                    responded = True

                elif response == 'n' or response == 'no':

                    print "\nYou can add an X-Ray Fan Beam source later by" +\
                          " creating an XRayFanBeam instance and" +\
                          " assigning it to your BackscatterSystem instance."

                    self.XRayFanBeam = None 

                    responded = True

                elif response == 'use defaults':

                    print "\nUsing defaults established in XRayFanBeam" +\
                          " Class."

                    print "See XRayFanBeam documentation to learn more."

                    self.XRayFanBeam = XRayFanBeam()

                    responded = True

                else:

                    print "\nYour response was not a valid option."

                    print  "Options: yes, no, use defaults \n"

                    response = raw_input('>>  ')                    

 
        # ====================================================================
        # Load or set up Detector
        # ====================================================================


        if isinstance(det_input, Detector):

            self.Detector = det_input

        else:

            responded = False

            print "\n2D Detector array doesn't exist, would you like to" +\
                  " build one?"

            print "Options: yes, no, use defaults \n"

            response = raw_input('>>  ')

            while not responded:

                if response == 'y' or response == 'yes':

                    self.Detector = Detector(build_manually = True)

                    responded = True

                elif response == 'n' or response == 'no':

                    print "\nYou can add a 2D detector array later by" +\
                          " creating a Detector instance and" +\
                          " assigning it to your BackscatterSystem instance."

                    self.Detector = None 

                    responded = True

                elif response == 'use defaults':

                    print "\nUsing defaults established in Detector" +\
                          " Class."

                    print "See Detector documentation to learn more."

                    self.Detector = Detector()

                    responded = True

                else:

                    print "\nYour response was not a valid option."

                    print  "Options: yes, no, use defaults \n"

                    response = raw_input('>>  ')                    


        # ====================================================================
        # Load or set up Uniform Collimation Grid
        # ====================================================================

        if isinstance(collimator_input, UniformCollimator):

            self.Collimator = collimator_input

        else:

            responded = False

            print "\nDetector Collimation grid doesn't exist, " +\
                  "would you like to  build one? "

            print "Options: yes, no, use defaults \n"

            response = raw_input('>>  ')

            while not responded:

                if response == 'y' or response == 'yes':

                    self.Collimator = UniformCollimator(build_manually = True)

                    responded = True

                    if self.Detector is not None:

                        responded = False

                        print "\nSince there is a Detector present, do you want"

                        print "to construct the collimation grid?" 

                        print  "Options: yes, no\n" 

                        response = raw_input('>>  ')                                

                        while not responded:

                            if response == 'y' or response == 'yes':

                                self.Collimator.construct_grid(self.Detector)

                                responded = True

                            elif response == 'n' or response == 'no':

                                responded = True

                            else:

                                print "\nYour response was not a valid option."

                                print  "Options: yes, no\n"

                                response = raw_input('>>  ')

                elif response == 'n' or response == 'no':

                    print "\nYou can add a detector collimation grid later " +\
                          "by creating a Collimator instance and " +\
                          "assigning it to your BackscatterSystem instance."

                    self.Collimator = None 

                    responded = True

                elif response == 'use defaults':

                    print "\nUsing defaults established in UniformCollimator"+ 
                          "Class"

                    print "See UniformCollimator documentation to learn more."

                    self.Collimator = UniformCollimator()

                    responded = True

                    if self.Detector is not None:

                        responded = False

                        print "\nSince there is a Detector present, do you want"

                        print "to construct the collimation grid?" 

                        print  "Options: yes, no\n" 

                        response = raw_input('>>  ')                                

                        while not responded:

                            if response == 'y' or response == 'yes':

                                self.Collimator.construct_grid(self.Detector)

                                responded = True

                            elif response == 'n' or response == 'no':

                                responded = True

                            else:

                                print "\nYour response was not a valid option."

                                print  "Options: yes, no\n"

                                response = raw_input('>>  ')

                else:

                    print "\nYour response was not a valid option."
                          
                    print "Options: yes, no, use defaults \n"

                    response = raw_input('>>  ')  

    def system_check(self):
        '''
        Analyzes the system and calculates the following features:

        depths         : A dictionary containg useful information on the
                         system. Gets created after calling the system_check
                         method. Each Key conatinsa list that's as long as the
                         number of pixels in the height of the detector.

                         Keys:

                         pix_center     : depth position (cm) determined by
                                          tracing a line from the center of a
                                          pixel to where it intersects with the
                                          fan beam. This line is orthogonal to
                                          the surface of the pixel. A negative 
                                          value indicates it is above the 
                                          surface of the object being imaged.

                         full_pix_range : tuples (cm), containing the 
                                          (deepest depth, shallowest depth) 
                                          that the pixel has full sight of. 
                                          Full sight meaning any photons that
                                          scatter in this range can pass
                                          through any part of the pixel's face.
                                          A negative value indicates it is
                                          above the surface of the object
                                          being imaged.

                         range          : tuples (cm), containing the 
                                          (deepest depth, shallowest depth)
                                          that the pixel can possibly see. So
                                          this range is always greater than the
                                          full_pix_range. It includes areas of
                                          the fan beam where photon scatters 
                                          are only partially visible to the
                                          pixel face. A negative value 
                                          indicates it is above the surface of
                                          the object being imaged.

                         y_res          : resolution in the y dimension (cm).
                                          Determined by using the path length
                                          of a photon scattering from the fan
                                          beam and going straight to the
                                          center of the pixel. This path length
                                          is the length of the line traced when 
                                          calculating pix_center

        '''
           
        self.depths = {'pix_center' : [], 'full_pix_range' : [], 
                       'range' : [], 'y_res' : []}

        pix_height = np.array([ -self.Detector.pixel_height*\
                                 np.cos(self.Detector.angle_radians),
                                 0 ,
                                 self.Detector.pixel_height*\
                                 np.sin(self.Detector.angle_radians)
                                ])

        pix_width = self.Detector.pixel_width - \
                    self.Collimator.collimator_thickness

        y_index = self.Detector.num_pixel_W + 1

        for i in range(self.Detector.num_pixel_H):

            # Top of Pixel
            Xp,Yp,Zp = self.Collimator.z_fins[i+1].v_vector

            # Bottom of Pixel
            xp,yp,zp = self.Collimator.z_fins[i].v_vector + \
                       self.Collimator.z_fins[i].a_vector

            # Center of Pixel
            xcp,ycp,zcp = self.Collimator.z_fins[i].v_vector + \
                          self.Collimator.z_fins[i].a_vector/2. + \
                          pix_height/2.

            # Collimator above pixel
            Xc,Yc,Zc = self.Collimator.z_fins[i+1].v_vector +\
                       self.Collimator.z_fins[i].c_vector

            # Collimator below pixel
            xc,yc,zc = self.Collimator.z_fins[i].v_vector + \
                       self.Collimator.z_fins[i].a_vector + \
                       self.Collimator.z_fins[i].c_vector

            # Calculate pix_center
            d = xcp/np.tan(self.Detector.angle_radians) - zcp

            path_length = xcp/np.sin(self.Detector.angle_radians)

            self.depths["pix_center"].append(d)

            # Calculate full_pix_range
            lower_limit = -find_intercept((xp,zp),(xc,zc)) 

            upper_limit = -find_intercept((Xp,Zp),(Xc,Zc)) 

            self.depths["full_pix_range"].append((lower_limit, upper_limit))

            # Calculate range
            lower_limit = -find_intercept((Xp,Zp),(xc,zc))

            upper_limit = -find_intercept((xp,zp),(Xc,Zc))

            self.depths["range"].append((lower_limit, upper_limit))

            # Calculate y_res
            length = self.Collimator.y_fins[y_index*i].length

            y_res = pix_width*(path_length - (length/2.))/(length/2.)

            self.depths["y_res"].append(y_res)

    def optimize_fan_beam(self, fan_angle=True):

        self.XRayFanBeam.optimize_with_detector(self.Detector, fan_angle)

    def optimize_collimator(self):

        self.Collimator.construct_grid(self.Detector)


class XRayFanBeam(object):
    '''
    Parameters Needed from User:
    ============================

    energy         : either pass it a single value or a file path containing
                     an X-ray energy distribution. It will then store it
                     as an array. See example below
                     (units: MeV)

    surface_offset : how high up (+Z) the source is positioned, where
                     x-rays are emitted

    fan_angle      : angle between vector pointing directly down (-Z) and 
                     vector pointing along edge of fan 
                     so fan_angle = total fan spread / 2
                     See figure below
                     (units: degrees)

                     fan_angle = 0 is a pencil beam


    Example energy spectrum:

    energy = array([[  1.00000000e-02,   1.24156743e+01],
                    [  1.10000000e-02,   1.04968200e+02],
                    [  1.20000000e-02,   3.89738325e+02]])

        The first column contains the energy values in MeV
        The second column contains the photon intensity or frequency

        The second column is what gets normalized and placed inside the
            normalized_energy attribute

    Fan Beam Figure: 

                   Point Source
                        X
                       /|\
                      / | \
                     /  |__<---- Fan Angle 
                    /   |   \
                   /    |    \ <- Edge of Fan Beam
                  /     |     \
                 /      |      \

                       /|\
                        |______ -Z Direction

    Default Values:
    ===============

    Energy         : 0.075 MeV

    surface_offset : 16. cm

    fan_angle      : 43.0 degrees

    Attributes:
    ===========

    normalized_energy : If the energy attribute is an array containing an
                        energy spectrum, then there will also be a copy of
                        that spectrum that has been normalized. This is 
                        necessary for importing the spectrum into MCNP

    '''

    def __init__(self, energy    = 0.075, surface_offset =16.0,
                       fan_angle = 43.0, build_manually = False):

        if build_manually == True:

            responded = False

            print '\nBuilding fan beam manually'

            print 'Are you importing the X-Ray energy spectra from a file?'

            print  "Options: yes, no\n"

            response = raw_input('>>  ') 

            while not responded:

                if response == 'y' or response == 'yes':

                    print "\nEnter file path\n"

                    file_path = raw_input('>>  ')

                    self.energy = convert_file_to_array(file_path)

                    responded = True

                elif response == 'n' or response == 'no':

                    print "\nEnter an X-ray energy: "

                    print "(in units of MeV)"

                    self.energy = get_value_from_user('float', float)

                    responded = True

                else:

                    print "\nYour response was not a valid option."

                    print  "Options: yes, no\n"

                    response = raw_input('>>  ')

            print "\nSurface offset of X-ray source point"

            print "(in centimeters)"

            self.surface_offset = get_value_from_user('float', float)

            print "\nAngle of the fan beam (the angle between the edge of"

            print "the beam and the -Z direction)"

            print "(in degrees)"

            self.fan_angle = get_value_from_user('float', float)

        else:

            if isinstance(energy, str): 

                self.energy = convert_file_to_array(energy)

            else:

                self.energy = energy

            self.surface_offset = surface_offset

            self.fan_angle = fan_angle

        if isinstance(self.energy , np.ndarray):

            self.normalized_energy = normalize_spectrum(self.energy)

        else: 

            self.normalized_energy = None

    def optimize_with_detector(self, detector, opt_fan_angle=True):

        #if not isinstance(detector, Detector):

        #    print "Error: Did not pass an instance of system_design.Detector"

        #    return

        if opt_fan_angle:

            fan_angle = np.arctan((detector.width/2.)/self.surface_offset)

            self.fan_angle = fan_angle*180/np.pi

        else:

            self.surface_offset = detector.width/(2.* \
                                  np.tan(self.fan_angle*np.pi/180))

class Detector(object):
    '''
    Object that contains all the necessary information to build a
    2D pixelated detector array.

    The "height" of the detector determines the resolution and number
    of pixels seeing into the object (Z dimension) being imaged.

    The "width" of the detector determines the resolution and number
    of pixels seeing accross the object (Y dimension) being imaged.

    The system therefore mechanically moves in the X dimension, scanning
    over the object to form a 3D image. 


                        2D Detector Array
        ___         ___ ___ ___ ___ ___ ___ ___
         |         |   |   |   |   |   |   |   |
         |         |___|___|___|___|___|___|___|
      height       |   |   |   |   |   |   |   |
   (num_pixel_H)   |___|___|___|___|___|___|___|
         |         |   |   |   |   |   |   |   |
        _|_        |___|___|___|___|___|___|___|
               
                   <---------- width ---------->
                           (num_pixel_W)
               
              


                          Individual Pixel 
                                ___
                      pixel    |   |
                      height   |___|

                               pixel 
                               width


    Parameters Needed from User:
    ============================

    num_pixel_H          : number of pixels in the height dimension of detector 
                           array

    num_pixel_W          : number of pixels in the width dimension of detector 
                           array

    pixel_height         : height of a single detector pixel (units: cm)

    pixel_width          : width of a single detector pixel (units: cm)
                           
    angle                : angle of detector array from object surface
                           EX: 0 degrees means detector is parallel to surface
                               90 degress means detector is parallel to fan beam
                           (units: degrees) 

    Default Values:
    ===============
    
    num_pixel_H     : 8 pixels

    num_pixel_W     : 30 pixels

    pixel_height    : 1.0 cm

    pixel_width     : 1.0 cm 

    angle           : 40.0 degrees

    Attributes:
    ===========

    height          : total height of the detector (units: cm)

    width           : total width of the detector (units: cm)

    angle_radians   : angle of the detector array in radians (units: radians)

    thickness       : how thick the detector is. This is held constant and
                      should be changed with caution, since it will affect
                      Tallying (image results) in the MCNP simulations.
                      (units: cm)

    The following attributes are modeled after the BOX Macrobody used in
    MCNP6. See MCNP manual for more information on Macrobodies. 
    See Collimator.construct_grid or Fin documentation for more information
    on these vectors and MCNP format.

    v_vector        : vector pointing from the origin to the corner of the
                      detector. This is helpful for plotting and MCNP input
                      generation (units: cm)

    a_vector        : vector pointing from v_vector to a corner of the
                      detector, defining its thickness. This is helpful for
                      plotting and MCNP input generation (units: cm)

    b_vector        : vector pointing from v_vector to a corner of the
                      detector, defining its width. This is helpful for
                      plotting and MCNP input generation (units: cm)

    c_vector        : vector pointing from v_vector to a corner of the
                      detector, defining its height (length). This is helpful
                      for plotting and MCNP input generation (units: cm)
    '''

    def __init__(self, num_pixel_H = 8    , num_pixel_W = 30     , 
                       pixel_height = 1.0 , pixel_width = 1.0    , 
                       angle = 40.0       , thickness = 0.35     ,
                                            build_manually = False):

        if build_manually == True:

            print "\nBuilding 2D Detector array manually"

            print "Number of pixels in the height dimension"

            self.num_pixel_H = get_value_from_user('integer', int)

            print "\nNumber of pixels in the width dimension"

            self.num_pixel_W = get_value_from_user('integer', int) 

            print "\nHeight of a single pixel in the detector array"

            print "(in centimeters)"

            self.pixel_height = get_value_from_user('float', float)

            print "\nWidth of a single pixel in the detector array"

            print "(in centimeters)"

            self.pixel_width = get_value_from_user('float', float)

            print "\nAngle of detector array from object surface"

            print "0 degrees means detector is parallel to object surface"

            print "90 degress means detector is parallel to fan beam"

            print "Recommended values are between 50 and 70 degrees"

            print "(in degrees)"

            self.angle = get_value_from_user('float', float)

            print "\nDetector thickness"

            print "(in centimeters)"

            self.thickness = get_value_from_user('float', float)


        else:

            self.num_pixel_H    = num_pixel_H

            self.num_pixel_W    = num_pixel_W

            self.pixel_height   = pixel_height

            self.pixel_width    = pixel_width

            self.angle          = angle

            self.thickness      = thickness

        self.height = self.num_pixel_H * self.pixel_height

        self.width  = self.num_pixel_W * self.pixel_width

        self.angle_radians = self.angle *np.pi/180

        self.v_vector = 0

        self.a_vector = 0

        self.b_vector = 0

        self.c_vector = 0

    def create_detector_vectors(self, starting_point, collimator_thickness):

        collimator_thickness_vec = np.array([
                        -collimator_thickness/2. * np.cos(self.angle_radians),
                         collimator_thickness/2. ,
                         collimator_thickness/2. * np.sin(self.angle_radians)
                                                ])

        self.v_vector = starting_point + collimator_thickness_vec

        self.a_vector = np.array([ self.thickness * np.sin(self.angle_radians),
                                   0 ,
                                   self.thickness * np.cos(self.angle_radians)
                                   ])

        self.b_vector = np.array([ 0, self.width, 0])

        self.c_vector = np.array([ -self.height * np.cos(self.angle_radians),
                                    0 , 
                                    self.height * np.sin(self.angle_radians)
                                   ])

class UniformCollimator(object):
    '''
    Object that contains information and modules necessry to build a
    uniform collimation grid to accompany the 2D pixelated detector array.

    This collimator will have a grid that is all a single length.
    You will need to run a system check using the method insdie the
    BackscatterSystem class to find out what the resolutions and 
    viewing range is for this collimator.

    In order to construct the grid, a fully built Detector instance 
    must be passed to the appropriate modules. Design data located 
    inside the Detector class are necessary to correctly design the 
    collimators.

    The grid is broken into two categories: Z collimators and Y collimators

    Z collimators are what provide depth discrimination (Z dimension in this
    framework) and are placed along the height of the detector, looking like
    fins (thin rectangular parrallelapipeds) extending in the Y dimension.
    They blind a detector pixel from X-ray backscatter signals generated at
    a specific depth in the object and deeper. Therefore the pixels only 
    receive signal from positions inside the object above that specified
    depth. There are as many Z collimators as there are pixels in the 
    "height" dimension of the detector, one for each row. 

    Y collimators are what provide the resolution and discrimination in
    signal across the object being imaged (Y dimension in this
    framework) and are placed along the width of the detector, looking like
    fins (thin recatngular parrallelapipeds) extending in the Z dimension.
    They blind a detector pixel from X-ray backscatter signals generated
    from other areas of the object that are supposed to be received by a pixel 
    closer to that position of scatter. These X-rays that scatter at wrong 
    angles, causing them to contribute signal to the wrong pixel, contribute 
    to the reduction in image resolution in the XY plane. There are as many
    Y collimators as there are pixels in the "width" dimension of the detector 
    plus one. 

    See Detector class documentation for diagrams and explanations of
    detector "height" and "width".


    Parameters Needed from User:
    ============================

    collimator_thickness : thickness of the detector collimators
                           (units: cm)

    length               : How long all teh fins will be inside the grid 
                           (units: cm)
                           
    surface_offset       : closest point/distance the system can be to
                           the object  (units: cm)

    beam_offset          : closest point/distance the system can be to 
                           the fanbeam (units: cm)

    material             : the type of material that the collimator is made
                           out of.

    Default Values:
    ===============

    collimator_thickness : 0.2 cm (2 mm) 

    length               : 5 cm 

    surface_offset       : 2.0 cm
    
    beam_offset          : 2.0 cm 

    material             : 'lead'

    Attributes:
    ===========

    depths                  : maximum depth a row of detector pixels can
                              see inside the object. Will be a list with the
                              max depth for each row in grid array (units: cm)

    y_resolution            : the resolution of the pixels in the Y dimension.
                              Will be a list with the y resolution for each row
                              in grid array (units: cm)

    z_fins                  : List containing each fin in the height direction.
                              Length of list will be equal to the number of 
                              pixels in the height dimension of the detector 
                              plus one.  

    y_fins                  : List containing each fin in the width direction.
                              Length of list will be equal to the number of 
                              pixels in the width dimension of the detector
                              plus one, multiplied by the number of pixels in
                              the height dimension of the detector.

    starting_point          : Vector pointing to the start point when creating
                              the list of fins for the collimation grid 
                              (units: cm)
    '''

    def __init__(self, collimator_thickness = 0.2  , length = 5.0      ,
                       surface_offset = 2.0        , beam_offset = 2.0 ,
                       material = 'lead'           , build_manually = False):

        if build_manually == True:

            print "\nBuilding collimator manually"

            print "Thickness of the collimation grid" 
            
            print "(be mindful that the thickness of the grid will be blocking"

            print "that much of each pixel face in the detector array)"

            print "(in centimeters)"

            self.collimator_thickness = get_value_from_user('float', float)

            print "\nLength of the the colllimators in the grid"

            print "(in centimeters)"

            self.length = get_value_from_user('float', float)

            print "\nClosest point/distance the system can be to the object"

            print " (offset distance of detector from surface of object)"

            print "(in centimeters)"

            self.surface_offset = get_value_from_user('float', float)

            print "\nClosest point/distance the system can be to the"

            print " X-ray fan beam source (offset distance of detector"

            print " from X-ray fan beam plane)"

            print "(in centimeters)"

            self.beam_offset = get_value_from_user('float', float)

            print "\nMaterial the collimator is made out of"

            self.material = raw_input('\n>>  ')

            # These attributes are filled using methods in this class

            self.z_fins = []

            self.y_fins = []

            self.depths = []

            self.starting_point = None

            #self.y_resolution = []

            self.max_depth = None

        else:

            self.collimator_thickness = collimator_thickness

            self.length         = length

            self.surface_offset = surface_offset

            self.beam_offset    = beam_offset

            self.material       = material 

            self.z_fins = []

            self.y_fins = []

            self.depths = []

            self.starting_point = None

            #self.y_resolution = []

            self.max_depth = None
        
    def construct_grid(self, detector):
        '''
        fins contain data in a format that follows that of 
        MCNP's BOX Macrobody (See MCNP Manual for more details)

        BOX v_x v_y v_z a_x a_y a_z b_x b_y b_z c_x c_y c_z

        where the v vector points to a corner of the box (from the origin).

        In this situation, a box is one of the fins that form the collimation
        grid.

        vectors a, b, and c each are orthogonal to eachother and point to
        the other corners of the box FROM THE STARTING POINT OF THE CORNER
        ESTABLISHED BY VECTOR v. To say it another way, for these the vectors
        the origin is located at the corner of the fin pointed at by vector v.

        z_fins will be a list containing each fin in the height direction.
        Length of list will be equal to the number of pixels in the height
        dimension of the detector plus one. 

        y_fins will be a list containing each fin in the width direction.
        Length of list will be equal to the number of pixels in the width
        dimension of the detector plus one, multiplied by the number of 
        pixels in the height dimension of the detector.

        a single fin is an object of the Fin class containing
        those 4 arrays which correspond to the vectors mentioned above, as well
        as other data, such as fin length,  width, and thickness.

        See Fin documentation for more information.

        Example:

            z_fins = [fin, fin, fin, fin , fin, ....]
        '''

        # ==================================================================
        # Clear the lists and variables
        # ==================================================================

        self.z_fins = []

        self.y_fins = []

        #self.depths = []

        # ==================================================================
        # Fill the self.starting_point attribute
        # ==================================================================

        self.find_starting_point(detector)

        # ==================================================================
        # this vector will be used to increment the starting point (v vector)
        #   to point to each fin needed in the height dimension. 
        #   It is determined by pixel height and detector angle.
        #   Called increment_xz because it only affects the vector in the 
        #   x and z dimension.
        # ==================================================================

        increment_xz = np.array([ 
                        -detector.pixel_height*np.cos(detector.angle_radians),
                        0 ,
                        detector.pixel_height*np.sin(detector.angle_radians) ])

        # ==================================================================
        # this vector will be used to increment the starting point (v vector)
        #   to point to each fin needed in the width dimension. 
        #   It is determined by pixel width. Called increment_y because it
        #   only affects the vector in the y dimension.
        # ==================================================================

        increment_y = np.array([ 0 , detector.pixel_width , 0 ])

        # ==================================================================
        # Start iterating over all the fins 
        # ==================================================================

        for i in range(detector.num_pixel_H + 1):

            # Update the position for the starting point of a Z fin

            current_point_xz = self.starting_point + i*increment_xz

            # Create a Z fin and append it to the list

            self.z_fins.append( self.create_z_fin( detector, i+1, 
                                                   current_point_xz) )

            # ===============================================================
            # After the first iteration, start creating the y fins too 
            # ===============================================================
            
            if i>0:

                # Start iterating over y fins

                for j in range(detector.num_pixel_W + 1):

                    current_point_xyz = current_point_xz + j*increment_y

                    self.y_fins.append( self.create_y_fin(  detector, 
                                    (i-1)*(detector.num_pixel_W+1) + (j + 1),
                                                           current_point_xyz ) 
                                        )

            # Need to calculate depth now that the fin was made 

            if i < detector.num_pixel_H:

                Xp,y,Zp = self.starting_point + (i+1)*increment_xz

                xl,y,zl = current_point_xz + self.z_fins[i].a_vector + \
                                                self.z_fins[i].c_vector

                self.depths.append( (Zp - zl*Xp/xl)/(Xp/xl - 1) )

        self.max_depth = self.depths[0]

    def find_starting_point(self, detector):
        '''
        finds vector starting point for the construction of the collimation
        grid

        Heavy use of trigonometry

        Also creates the vectors that define the dimensions of the Detector
        instance being used. 
        '''

        start_x = self.beam_offset + \
                  detector.height  * np.cos(detector.angle_radians) + \
                  self.length * np.sin(detector.angle_radians) + \
                  self.collimator_thickness * np.cos(detector.angle_radians)

        start_y = -(detector.width + self.collimator_thickness)/2.

        start_z = self.surface_offset + \
                  self.length * np.cos(detector.angle_radians) - \
                  self.collimator_thickness * np.sin(detector.angle_radians)

        self.starting_point = np.array([ start_x, start_y, start_z ])

        detector.create_detector_vectors(self.starting_point, 
                                         self.collimator_thickness)

    def create_z_fin(self, detector, index, position):

        # ===============================================================
        # Instantiate a Fin object and assign values to its attributes
        # ===============================================================

        fin = Fin('Z', index, self.collimator_thickness)

        fin.width = detector.width + self.collimator_thickness

        # Fin starting position (v_vector)

        fin.v_vector = position

        # Fin thickness vector (a_vector)

        fin.a_vector = np.array([
                    -self.collimator_thickness*np.cos(detector.angle_radians),
                     0 ,
                     self.collimator_thickness*np.sin(detector.angle_radians)
                                 ])

        # Fin width vector (b_vector)

        fin.b_vector = np.array([ 0 , 
                                  detector.width + self.collimator_thickness ,
                                  0  ])



        # Fin length vector (c_vector)

        fin.length = self.length

        fin.c_vector = np.array([ -self.length*np.sin(detector.angle_radians), 
                                   0 ,
                                  -self.length*np.cos(detector.angle_radians)])

        return fin

    def create_y_fin(self, detector, index, position):

        # ===============================================================
        # Create useful variables from input arguments
        # ===============================================================        

        actual_pixel_height = detector.pixel_height - self.collimator_thickness

        actual_pixel_width = detector.pixel_width - self.collimator_thickness

        # ===============================================================
        # Instantiate a Fin object and assign values to its attributes
        # ===============================================================

        fin = Fin('Y', index, self.collimator_thickness)

        fin.width = actual_pixel_height

        # Fin starting position (v_vector)

        fin.v_vector = position

        # Fin thickness vector (a_vector)

        fin.a_vector = np.array([ 0, self.collimator_thickness, 0 ])

        # Fin width vector (b_vector)

        fin.b_vector = np.array([
                            actual_pixel_height*np.cos(detector.angle_radians), 
                            0 ,
                            -actual_pixel_height*np.sin(detector.angle_radians),  
                                  ])

        fin.length = self.length

        fin.c_vector = np.array([ -self.length*np.sin(detector.angle_radians), 
                                   0 ,
                                  -self.length*np.cos(detector.angle_radians)])

        return fin 

class Fin(object):
    '''
    This is a class specifically to be used by the Collimator class. 

    The Collimator is a collimation grid that is composed of many 
    individual fins.  

    a Fin is characterized by its dimensions: thickness, width, and length

    thickness is determined by the Collimator class and will be the same for
    every Fin in a Backscatter System. 

    width is determined by detector dimensions:
        detector width if it's a horizontal fin collimating for
        depth (Z) discrimination

        detector pixel height if it's a vertical fin collimating for
        surface (Y) discrimination

    length is the factor that must be calculated based on the desired
    system properties such as maximum depth, number of pixels, etc.

                        Detector Face
                             / \
                              |

                         __________  ___ 
              thickness /_________/|  |
                        |         ||  | length 
                        |         ||  |
                        |_________|/ _|_

                        |----------|
                            width

                              |
                             \ /
                        Object Surface

    Attributes:
    ===========

    collimation_dimension    : The dimension in which the fin is trying to
                               collimate its detector pixel. Z is collimating
                               the depth dimension and Y is collimating accross
                               the surface dimension 

    index_number             : The index of this specific fin amongst the other
                               fins that share the same collimation_dimension.
                               1 (Not 0 like most python indexes) is the first
                               in the collection. 

    thickness                : Determined by the Collimator class and will be 
                               the same for every Fin in a Backscatter System.
                               (units: cm)

    width                    : Determined by detector width if it's a Fin
                               with collimation_dimension = Z

                               Determined by detector pixel height if it's
                               a Fin with collimation_dimension = Y
                               (units: cm)

    length                   : Must be calculated based on the desired
                               system properties such as maximum depth,
                               number of pixels, etc. (units: cm)

    The following attributes are modeled after the BOX Macrobody used in
    MCNP6. See MCNP manual for more information on Macrobodies. 
    See Collimator.construct_grid documentation for more information
    on fin attributes and MCNP format.

    v_vector                 : Base vector. Points from origin to a corner
                               of the box that defines this fin
                               (units: cm)

    a_vector                 : Dimension vector describing fin thickness. 
                               Provides direction and therefore magnitude (size)
                               of a side of the box. This vector's components 
                               are relative to the corner pointed at by v_vector 
                               (units: cm)

    b_vector                 : Dimension vector describing fin width. 
                               Provides direction and therefore magnitude (size)
                               of a side of the box. This vector's components 
                               are relative to the corner pointed at by v_vector
                               (units: cm)

    c_vector                 : Dimension vector describing fin length. 
                               Provides direction and therefore magnitude (size)
                               of a side of the box. This vector's components 
                               are relative to the corner pointed at by v_vector
                               (units: cm)
    '''
    def __init__(self, collimation_dimension, index_number, thickness):
        
        self.collimation_dimension = collimation_dimension

        self.index_number = index_number

        self.thickness = thickness

        self.width = 0

        self.length = 0

        self.v_vector = 0

        self.a_vector = 0

        self.b_vector = 0

        self.c_vector = 0

def get_value_from_user(dtype, convert_func):
    '''
    helper function when system components are being built manually
    and require user input. 

    Tries to convert the user input into the right data type based on
    the conversion function that was passed. If successful, it loads it
    into variable 

    Example: passing float, will try converting user input into a float

    If an error is raised, it will tell the user and ask them to try again. 
    '''
    correctly_input = False

    while not correctly_input:

        try:

            response = convert_func(raw_input('\n>>  '))

            correctly_input = True

        except ValueError:

            print '\nIncorrect, the input must be a ' + dtype

            print 'Try again: '

    return response

def convert_file_to_array(file_path):

    opened = False

    spectrum = []

    while not opened:
    
        try:

            f = open(file_path, 'r')

            opened = True

        except IOError:

            print '\nFile and/or path does not exist'

            print 'Enter a file path or quit\n'

            response = raw_input('>>  ')

            if response == 'q' or response == 'quit':

                return 75.0

            else:

                file_path = response

    contents = f.readlines()

    t =[]

    for line in contents:

        # Ignores lines beginning with '#' and any blank lines

        if '#' not in line and line != '\n':

            # Remove 'newline' character and split up values

            spectrum.append([float(value) for value in \
                             line.rstrip('\n').split()])

    return np.array(spectrum)

def normalize_spectrum(spectrum):
    
    norm_spectrum = np.copy(spectrum)

    norm_spectrum[:,1] = norm_spectrum[:,1]/norm_spectrum[:,1].sum()

    return norm_spectrum

def find_intercept(point1, point2):
    
    x1, y1 = point1

    x2, y2 = point2

    m = (y2 - y1)/(x2 - x1)

    b = y1 - m*x1

    return b

