'''
Contains classes and modules for the creation of a 2D backscatter system
for MCNP input files and the creation of images using said system.
'''

import system_design as design
import numpy as np
import os

# This is a dictionary that contains materials that can be used for the
#   Material class. Keys are common names for the material and Values are
#   a tuple containing the material density and then the ZAID identifier 
#   and fraction of that isotope in the material
#   (See MCNP manual for ore details about material card specifications)

LIST_OF_KNOWN_MATERIALS = {'lead' : ('-11.34', '82000 -1.00'),
                           'tungsten' : ('-19.30', '74000 -1.00'),
                           'air' : ('-0.0012041', '6000 -0.000124 ' + \
                                                  '7000 -0.755268 ' + \
                                                  '8000 -0.231781 ' + \
                                                  '18000 -0.012827')} 
# ==========================================================================
#           Classes 
# ==========================================================================

class Cell(object):
    '''
    Contains the properties describing a cell card in MCNP
    '''

    def __init__(self, name=None, material_num=None, material_density=None,
                                        surfaces=None, importance='P=1' ):

        self.name = name

        self.material_num = material_num

        self.material_density = material_density

        self.surfaces = surfaces

        self.importance = importance

class Surface(object):
    '''
    Contains the properties describing a surface card in MCNP
    '''

    def __init__(self, name=None, TR=None, type=None, type_variables=''):

        self.name = name

        self.TR = TR

        self.type = type

        self.type_variables = type_variables

class SourceDefinition(object):
    '''
    Contains the properties describing a source definiion (SDEF) card in MCNP
    '''

    def __init__(self, position = None, energy = None, particle = '2', 
                         vector = None, direction  = None, TR= None,
                         weight = None):

        self.name = 'SDEF'

        self.position = position

        self.energy = energy

        self.particle = particle

        self.vector = vector

        self.direction = direction

        self.TR = TR

        self.weight = weight

        self.distributions = {}

    def make_energy_distribution(self, energies):

        name = self.energy.split('d')[-1]

        SI = 'SI'+name + ' A ' + ' '.join(map(str, energies[:,0]))

        SP = 'SP'+name + ' ' + ' '.join(map(str, energies[:,1]))

        self.distributions.update({self.energy : [SI, SP]})

    def make_direction_distribution(self, angle):

        name = self.direction.split('d')[-1]

        angle_bins = [-1, round(np.cos(angle*np.pi/180), 5), 1]

        SI = 'SI'+name + ' ' + ' '.join(map(str, angle_bins))

        probabilities = [ 0,
                round((1 - angle_bins[0] - 1 + angle_bins[1])/2., 5), 
                round((1 - angle_bins[1] - 1 + angle_bins[2])/2., 5) ]

        SP = 'SP'+name + ' ' + ' '.join(map(str, probabilities))

        SB = 'SB'+name + ' 0 0 1'

        self.weight = str(round(1./probabilities[-1], 5))

        self.distributions.update({self.direction : [SI, SP, SB]})

class CoordinateTransformation(object):
    '''
    Contains the properties describing a coordinate transformation (TR) card 
    in MCNP
    '''

    def __init__(self, name = None, degrees = True, displacement_vector = None):

        self.name = name

        self.degrees = degrees # If True, then name has *TR instead of TR

        self.displacement_vector = displacement_vector

        if self.degrees:

            self.rotation_matrix = np.array([[0, 90, 90], 
                                             [90, 0, 90], 
                                             [90, 90, 0]])

        else:

            self.rotation_matrix = np.array([[1, 0, 0], 
                                             [0, 1, 0], 
                                             [0, 0, 1]])

        self.M = '1'

class Material(object):
    '''
    Contains the properties describing a material card in MCNP
    '''

    def __init__(self, description = None, name = None ):

        self.description = description

        self.name = name

        if LIST_OF_KNOWN_MATERIALS.has_key(self.description):
        
            self.composition = LIST_OF_KNOWN_MATERIALS[description][1]

        else:

            self.composition = None
        
class Tally(object):
    '''
    Contains the properties describing a material card in MCNP
    '''
    
    def __init__(self, type, name, radiation = 'P', type_variables = ''):

        self.type = type

        self.name = name

        self.radiation = radiation

        self.type_variables = type_variables
        
# ==========================================================================
#           Modules
# ==========================================================================

# ==========================================================================
#           Modules Specific to Creating this 2D Backscatter System
# ==========================================================================

def create_input_deck(system, file_path = 'system_template.i', 
                              TR = '1', mode = 'P', nps = 5000000):
#
    '''
    Takes the backscatter system and exports it into a readable format
    in an MCNP input file

    system is an instance of system_design.BackscatterSystem

    file_path is either a file name or a path to (and including) the file

    TR is the coordinate transformation name that will be used when creating
        surface cards and source definition cards.
    '''

    if not isinstance(system, design.BackscatterSystem):

        print "ERROR: must pass an instance of" 

        print "system_design.BackscatterSystem"

        return

    # ====================================================================
    # Either create a new file or overwrite one. Check if new directories
    #    were requested too.
    # ====================================================================

    if os.path.isfile(file_path):

        responded = False

        print "\nWARNING: This file already exists, do you wish to proceed?"

        print "Doing so will overwrite the file."

        print "(yes/no)\n"

        response = raw_input('>>  ')

        while not responded:

            if response == 'y' or response == 'yes':

                input_deck = open(file_path, 'w')

                responded = True

            elif response == 'n' or response == 'no':

                return

            else:

                    print "\nYour response was not a valid option."

                    print  "Options: yes, no \n"

                    response = raw_input('>>  ')

    elif len(file_path.split('/')) > 1:

        directory = file_path.rpartition('/')[0]

        if not os.path.isdir(directory):

            os.makedirs(directory)

            input_deck = open(file_path, 'w')

        else:

            input_deck = open(file_path, 'w')

    else:

        input_deck = open(file_path, 'w')

    # ====================================================================
    # Create Cell, Surface, and Material cards for system
    # ====================================================================

    # Collimator Data

    collimator_cells, collimator_surfaces, collimator_material = \
                                            create_collimation_cards(system, TR)

    # If it's a fan beam, need to create collimators that will shape the
    #   x-ray source into a fan beam

    if system.XRayFanBeam.fan_angle > 0.0:

        xray_collimator_cells, xray_collimator_surfaces = \
                                       create_xray_collimation_cards(system, TR)

    # Create Detector Cell where mesh tally will be placed inside of 

    detector_cell, detector_surface = create_detector_cards(system,
                                                            TR=str(int(TR) + 1))

    # Figure out maximum dimension for determining the size of the problem

    xmax, ymax, z = system.Detector.v_vector + system.Detector.a_vector

    x, y, zmax = system.Detector.v_vector + system.Detector.a_vector + \
                                            system.Detector.c_vector

    max_dimension = np.max([system.XRayFanBeam.surface_offset,
                            xmax, ymax, zmax ])

    problem_boundary_surface = Surface(name='99998', type='so', TR='' ,
                                type_variables=str(round(max_dimension*2., 3)))

    inside_boundary_material = Material(description='air' , name='900')

    inside_boundary_cell = create_inside_boundary_cell(collimator_cells +\
                                                       xray_collimator_cells +\
                                                       [detector_cell], 
                                                       problem_boundary_surface,
                                                       inside_boundary_material)

    outside_boundary_cell = Cell(name = '99999',
                                 material_num = '',
                                 material_density = '0',
                                 surfaces = problem_boundary_surface.name,
                                 importance = 'P=0' )

    # Source Definition (X Ray Fan Beam Data)

    source_def = create_source_definition(system, TR)

    # Coordinate Transformation (will be used when scanning over an object)

    transform = CoordinateTransformation(name = TR, 
                                         displacement_vector = '0 0 0')

    # Create Mesh Tally for Detector

    mesh_tally = create_detector_tally(system, TR=str(int(TR) + 1), 
                                                    energy_mesh=True)

    # Create Transform for detector 

    tally_transform = create_detector_transform(system, TR=str(int(TR) + 1))

    # ====================================================================
    # Start writing the contents for the input file
    # ====================================================================

    # Write Cell Cards
    # ================

    write_header('Cells', input_deck)
    
    write_subheader('Collimation Grid', input_deck)

    for cell in collimator_cells:

        write_cell_card(cell, input_deck)

    if system.XRayFanBeam.fan_angle > 0.0:

        write_subheader('X-Ray Fan Beam Collimators', input_deck)

        for cell in xray_collimator_cells:

            write_cell_card(cell, input_deck)

    write_subheader('Detector', input_deck)

    write_cell_card(detector_cell, input_deck)

    write_subheader('Phantom', input_deck)

    write_subheader('Outside of Problem', input_deck)

    write_cell_card(outside_boundary_cell, input_deck)

    write_subheader('Inside of Problem', input_deck)

    write_cell_card(inside_boundary_cell, input_deck)

    # Write Surface Cards
    # ===================

    input_deck.write('\n')

    write_header('Surfaces', input_deck)

    write_subheader('Collimation Grid', input_deck)

    for surface in collimator_surfaces:

        write_surface_card(surface, input_deck)

    if system.XRayFanBeam.fan_angle > 0.0:

        write_subheader('X-Ray Fan Beam Collimators', input_deck)

        for surface in xray_collimator_surfaces:

            write_surface_card(surface, input_deck)

    write_subheader('Detector', input_deck)

    write_surface_card(detector_surface, input_deck)
    
    write_subheader('Phantom', input_deck)

    write_subheader('Problem Boundary', input_deck)

    write_surface_card(problem_boundary_surface, input_deck)

    # Write Data Cards
    # ================

    input_deck.write('\n')

    write_header('Data', input_deck)

    # Write Coordinate Transformation Card

    write_subheader('Coordinate Transformation', input_deck)

    write_coordinate_transform_card(transform, input_deck)

    write_coordinate_transform_card(tally_transform, input_deck)

    # Write Material Cards

    write_subheader('Materials', input_deck)

    write_material_card(collimator_material, input_deck)

    write_material_card(inside_boundary_material, input_deck)

    # Write Source Definition

    write_subheader('Source Definition', input_deck)

    write_source_definition_card(source_def, input_deck)

    # Write Tallies

    write_subheader('Tallies', input_deck)

    write_tally_card(mesh_tally, input_deck)

    # Write Ending 

    input_deck.write('c\n')

    input_deck.write('MODE ' + mode + '\n')

    input_deck.write('NPS ' + str(nps) + '\n')

    input_deck.write('PRINT 110 160 161 162')

    input_deck.close()

def create_collimation_cards(system, TR):
    '''
    Creates surface, cell, and material cards for the collimation grid

    '''
    if not isinstance(system.Collimator, design.Collimator) and \
       not isinstance(system.Collimator, design.UniformCollimator):

        print "\nError: did not receive an instance of system_design.Collimator"

        print "or system_design.UniformCollimator"

        return

    coll = system.Collimator

    cells = []

    surfaces = []

    material = Material(description = coll.material.lower(), name = '200')

    num_z = coll.z_fins[-1].index_number
    
    for fin in coll.z_fins:

        surface = Surface(name = str(fin.index_number) ,
                          TR = TR , 
                          type = 'BOX')

        # Round values to 3 decimal places

        values =  map(lambda x: round(x, 3), np.hstack([fin.v_vector, 
                                                        fin.a_vector*0.99, 
                                                        fin.b_vector, 
                                                        fin.c_vector]))

        # Combine values into a string

        surface.type_variables = ' '.join(map(str, values))

        cell = Cell(name=str(fin.index_number),
                    material_num=material.name, 
             material_density=LIST_OF_KNOWN_MATERIALS[material.description][0],
                    surfaces='-'+surface.name)

        surfaces.append(surface)

        cells.append(cell)

    for fin in coll.y_fins:

        surface = Surface(name = str(num_z + fin.index_number) ,
                          TR = TR , 
                          type = 'BOX')

        # Round values to 3 decimal places

        values =  map(lambda x: round(x, 3), np.hstack([fin.v_vector, 
                                                        fin.a_vector, 
                                                        fin.b_vector*.99, 
                                                        fin.c_vector]))

        # Combine values into a string

        surface.type_variables = ' '.join(map(str, values))

        cell = Cell(name=str(num_z + fin.index_number),
                    material_num=material.name, 
             material_density=LIST_OF_KNOWN_MATERIALS[material.description][0],
                    surfaces='-'+surface.name)

        surfaces.append(surface)

        cells.append(cell)

    return cells, surfaces, material

def create_inside_boundary_cell(other_cells, surface, material):
    
    problem_cell = Cell(name=surface.name,
                        material_num=material.name, 
              material_density=LIST_OF_KNOWN_MATERIALS[material.description][0])
    
    problem_cell.surfaces = '-' + surface.name

    for cell in other_cells:

        problem_cell.surfaces += ' #' + cell.name

    return problem_cell

def create_xray_collimation_cards(system, TR):

    angle = system.XRayFanBeam.fan_angle *np.pi/180.

    zmin = system.XRayFanBeam.surface_offset - 2

    zmax = zmin + 1

    ymin = -np.tan(angle)*2 - 1

    ymax = np.abs(ymin)

    xmin = system.step_size/(system.XRayFanBeam.surface_offset + \
                                     system.Collimator.max_depth)

    xmax = ymax

    values = map(lambda x: round(x, 3), [xmin, xmax, ymin, ymax, zmin, zmax])

    surface1 = Surface(name = '99991', TR = TR, type = 'RPP')  

    surface1.type_variables = ' '.join(map(str, values))

    surface2 = Surface(name = '99992', TR = TR, type = 'RPP')

    values = map(lambda x: round(x, 3), [-xmax, -xmin, ymin, ymax, zmin, zmax])

    surface2.type_variables = ' '.join(map(str, values))

    cell1 = Cell(name = '99991', material_num = '', material_density = '0',
                           surfaces = '-'+surface1.name, importance = 'P=0')

    cell2 = Cell(name = '99992', material_num = '', material_density = '0',
                           surfaces = '-'+surface2.name, importance = 'P=0')

    return [cell1, cell2], [surface1, surface2]

def create_source_definition(system, TR):
    
    if not isinstance(system.XRayFanBeam, design.XRayFanBeam):

        print "\nError: did not receive an instance of " +\
              "system_design.XRayFanBeam"

        return

    fan = system.XRayFanBeam

    sdef = SourceDefinition(position = '0 0 '+str(fan.surface_offset),
                            vector = '0 0 -1', TR = TR)

    if isinstance(fan.energy, np.ndarray):

        sdef.energy = 'd1'

        sdef.make_energy_distribution(fan.normalized_energy)

    else:

        sdef.energy =  str(fan.energy)

    if fan.fan_angle > 0:

        sdef.direction = 'd2'

        sdef.make_direction_distribution(fan.fan_angle)

    else:

        sdef.direction = '1'

    return sdef

def create_detector_tally(system, TR, energy_mesh = False):

    if not isinstance(system.Detector, design.Detector):

        print "\nError: did not receive an instance of " +\
              "system_design.Detector"

        return

    det = system.Detector

    tally = Tally(type = 'FMESH', name = '14')

    tally.type_variables += 'GEOM=xyz ORIGIN=0 0 0'

    tally.type_variables += ' TR=' + TR + '\n'

    i = ' '.join([' '*8, 'IMESH=' + str(det.thickness) , 'IINTS=1\n'])

    j = ' '.join([' '*8, 'JMESH=' + str(det.width), 
                  'JINTS=' + str(det.num_pixel_W) + '\n'])

    k = ' '.join([' '*8, 'KMESH=' + str(det.height), 
                  'KINTS=' + str(det.num_pixel_H) + '\n'])

    tally.type_variables += ' '.join(['', i, j, k])

    if energy_mesh:

        if isinstance(system.XRayFanBeam.energy, np.ndarray):

            e = ' '*10 + 'EMESH=' + str(system.XRayFanBeam.energy[-2,0])

            e += ' EINTS=' + str(int(system.XRayFanBeam.energy[-2,0]/0.005)) +\
                 '\n'

        else: 

            e = ' '*10 + 'EMESH=' + str(system.XRayFanBeam.energy)

            e += ' EINTS=' + str(int(system.XRayFanBeam.energy/0.005)) + '\n'

    tally.type_variables += e

    return tally

def create_detector_transform(system, TR):
    
    transform = CoordinateTransformation(name = TR)

    values =  map(lambda x: round(x, 3), 
                  system.Detector.v_vector + np.array([0.01, 0 , 0.01]))

    transform.displacement_vector = ' '.join(map(str, values))

    angle = system.Detector.angle

    transform.rotation_matrix = np.array([[90-angle, 90, angle],
                                          [90, 0, 90],
                                          [180-angle, 90, 90-angle]])

    return transform

def create_detector_cards(system, TR):
    
    if not isinstance(system.Detector, design.Detector):

        print "\nError: did not receive an instance of " +\
              "system_design.Detector"

        return

    det = system.Detector

    xmin, ymin, zmin = det.v_vector

    surface = Surface(name = '99990', TR = TR, type = 'RPP')

    values = map(lambda x: round(x, 3), [0, det.thickness,
                                         0, det.width,
                                         0, det.height])

    surface.type_variables = ' '.join(map(str, values))

    cell = Cell(name = '99990', material_num = '', material_density = '0',
                           surfaces = '-'+surface.name)

    return cell, surface

# ==========================================================================
#           General Modules for Writing to an MCNP Input File
# ==========================================================================

def write_subheader(text, deck):

    text = ' ' + text + ' '

    deck.write('c\n')

    deck.write('c ' + text.center(78, '-') + '\n')

    deck.write('c\n')

def write_header(text, deck):

    text = ' ' + text + ' '
    
    deck.write('c ' + '-'*78 + '\n')

    deck.write('c ' + text.center(78, '-') + '\n')

    deck.write('c ' + '-'*78 + '\n')

    deck.write('c\n')

def write_cell_card(cell, deck):

    if not isinstance(cell, Cell):

        print "\nError: did not receive an instance of mcnp.Cell"

        return

    head = ' '.join([cell.name, cell.material_num, cell.material_density])

    tail = ' '.join([cell.surfaces, 'IMP:' + cell.importance])

    write_line(head, tail, deck)

def write_surface_card(surface, deck):

    if not isinstance(surface, Surface):

        print "\nError: did not receive an instance of mcnp.Surface"

        return

    head = ' '.join([surface.name, surface.TR, surface.type])

    tail = surface.type_variables

    write_line(head, tail, deck)

    #line = ' '.join([surface.name, surface.TR, 
    #                 surface.type, surface.type_variables])

    #deck.write(line + '\n')

def write_coordinate_transform_card(transform, deck):
    
    if not isinstance(transform, CoordinateTransformation):

        print "\nError: did not receive an instance of " +\
              "mcnp.CoordinateTransformation"

        return

    rotation = ' '.join(map(str, np.hstack(transform.rotation_matrix)))

    line = ' '.join(['TR' + transform.name, transform.displacement_vector, 
                      rotation, transform.M])

    if transform.degrees:

        deck.write('*' + line + '\n')

    else:

        deck.write(line + '\n')

def write_material_card(material, deck):

    if not isinstance(material, Material):

        print "\nError: did not receive an instance of " +\
              "mcnp.Material"

        return
    
    deck.write('c\nc ' + material.description.capitalize() + '\nc\n')

    name = 'm' + material.name + ' '

    zaids = material.composition.split()

    zaids.reverse()

    deck.write(name + ' '.join([ zaids.pop(), zaids.pop() ]) + '\n')

    for i in range(len(zaids)/2):

        zaid_and_fraction = ' '.join([ zaids.pop(), zaids.pop() ])

        deck.write( zaid_and_fraction.rjust(len(name)+len(zaid_and_fraction)) \
                    + '\n')

def write_source_definition_card(sdef, deck):
    
    if not isinstance(sdef, SourceDefinition):

        print "\nError: did not receive an instance of mcnp.SourceDefinition"

        return



    head = ' '.join(['SDEF',
                     'POS='+sdef.position,
                     'ERG='+sdef.energy])

    if sdef.weight:

        head += ' WGT='+sdef.weight



    tail = ' '.join(['PAR='+sdef.particle,
                     'VEC='+sdef.vector,
                     'DIR='+sdef.direction,
                     'TR='+sdef.TR])

    write_line(head, tail, deck)

    for distribution in sdef.distributions.keys():

        for line in sdef.distributions[distribution]:

            head, temp, tail = line.partition(' ')

            # In case the distribution is specified with a character,
            #   like A,D,L etc.

            if tail.split()[0].isalnum():

                head += ' ' + tail[0]

                tail = tail[1:]

            else:

                head += '  '

            write_line(head, tail, deck)

def write_tally_card(tally, deck):

    if not isinstance(tally, Tally):

        print "\nError: did not receive an instance of mcnp.Tally"

        return

    if 'FMESH' in tally.type.upper():

        deck.write(tally.type + tally.name + ':' + tally.radiation)

        deck.write(' ' + tally.type_variables)

    else:

        head = tally.type + tally.name + ':' + tally.radiation

        tail = tally.type_variables

        write_line(head, tail, deck)

def plot_basis(system):
    '''
    Helps you visualize collimation grid in the MCNP plotter
    by telling you what to input for the BASIS command in the 
    plotter.
    '''

    angle = system.Detector.angle_radians

    x = str(-np.cos(angle))

    y = '0'

    z = str(np.sin(angle))

    print 'For Y collimators:'

    print 'basis = 0 1 0 ' + ' '.join([x,y,z])

    print 'For Z collimators:'

    print 'basis = 1 0 0 0 0 1'

def write_line(head, tail, deck):

    # Check if line needs to be written over multiple lines

    if len(head + tail) > 80:

        pad = len(head) + 1

        # Find string length of element with longest name

        max_length = np.max(map(lambda x: len(x), tail.split())) + 1

        # determine how many surfaces can fit on a single line 

        num_per_line = (80 - pad)/max_length

        elements = tail.split()

        deck.write(head + ' ' + ' '.join( elements[0:num_per_line] ) + '\n')

        for i in range((len(elements)-num_per_line)/num_per_line + 1):

            current_row = ' '.join( elements[ (i+1)*num_per_line : \
                                              (i+1)*num_per_line + num_per_line] )

            deck.write( current_row.rjust( pad + len(current_row) ) + '\n')

    else:

        deck.write(' '.join([head, tail + ' \n']))
