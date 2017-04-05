import vtk
import system_design
import numpy as np

def system(system):

    
    # Create the graphics structure. The renderer renders into the render
    # window. The render window interactor captures mouse events and will
    # perform appropriate camera or actor manipulation depending on the
    # nature of the events.
    ren = vtk.vtkRenderer()
    renWin = vtk.vtkRenderWindow()
    renWin.AddRenderer(ren)
    iren = vtk.vtkRenderWindowInteractor()
    iren.SetRenderWindow(renWin)
     
    # Add the Z fins to the renderer
    for fin in system.Collimator.z_fins:

        fin_coords = [fin.v_vector, fin.b_vector, fin.a_vector, 
                                                          fin.c_vector] 

        finAct = create_hexahedron_actor(fin_coords)

        ren.AddActor(finAct)

    # Add the Y fins to the renderer
    for fin in system.Collimator.y_fins:

        fin_coords = [fin.v_vector, fin.a_vector, fin.b_vector, fin.c_vector] 

        finAct = create_hexahedron_actor(fin_coords)

        ren.AddActor(finAct)

    # Add the Detector to the renderer
    det_coords = [system.Detector.v_vector, system.Detector.b_vector,
                  system.Detector.c_vector, system.Detector.a_vector]

    ren.AddActor(create_hexahedron_actor(det_coords, color=(92, 144, 249)))

    # Add the X-ray Beam to the renderer
    beam_coords = get_beam_vectors(system)

    ren.AddActor(create_pyramid_actor(beam_coords))

    # Add plane showing Object surface
    center = (system.Detector.width, 0, 0 )
    point1 = (0, 4*system.Detector.width, 0)
    point2 = (4*system.Detector.width, 0, 0)

    ren.AddActor(create_plane_actor(center, point1, point2, 
                                        color=(100,100,100)))

    # Set the background and size
    ren.SetBackground(0.9, 0.9, 0.9)
    renWin.SetSize(500, 500)
     
    # This allows the interactor to initalize itself. It has to be
    # called before an event loop.
    iren.Initialize()
     
    # We'll zoom in a little by accessing the camera and invoking a "Zoom"
    # method on it.
    ren.ResetCamera()
    ren.GetActiveCamera().Zoom(1.5)
    renWin.Render()
     
    # Start the event loop.
    iren.Start()

    
def create_hexahedron_actor(hex_coords, color = (186, 200, 214)):

    if np.max(color) > 1.0:

        color = (color[0]/255., color[1]/255., color[2]/255.)

    src_vector, width_vector, thickness_vector, length_vector = hex_coords

    P0 = src_vector
    P1 = src_vector + width_vector
    P2 = src_vector + width_vector + length_vector
    P3 = src_vector + length_vector
    P4 = src_vector + thickness_vector
    P5 = src_vector + thickness_vector + width_vector
    P6 = src_vector + thickness_vector + width_vector + length_vector
    P7 = src_vector + thickness_vector + length_vector

    # Create the points
    points = vtk.vtkPoints()
    points.InsertNextPoint(P0)
    points.InsertNextPoint(P1)
    points.InsertNextPoint(P2)
    points.InsertNextPoint(P3)
    points.InsertNextPoint(P4)
    points.InsertNextPoint(P5)
    points.InsertNextPoint(P6)
    points.InsertNextPoint(P7)
     
    # Create a hexahedron from the points
    hexahedron = vtk.vtkHexahedron() 

    for i in range(8):
        hexahedron.GetPointIds().SetId(i,i)
     
    # Add the hexahedron to a cell array
    hexs = vtk.vtkCellArray()
    hexs.InsertNextCell(hexahedron)
     
    # Add the points and hexahedron to an unstructured grid
    uGrid = vtk.vtkUnstructuredGrid()
    uGrid.SetPoints(points)
    uGrid.InsertNextCell(hexahedron.GetCellType(), hexahedron.GetPointIds())
     
    # Create mapper and Actor
    mapper = vtk.vtkDataSetMapper()
    mapper.SetInputData(uGrid)
     
    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    actor.GetProperty().SetColor(color)

    return actor

def create_pyramid_actor(pyramid_coords, color = (239, 57, 57)):

    if np.max(color) > 1.0:

        color = (color[0]/255., color[1]/255., color[2]/255.)

    # Using vectors in the format of get_beam_vectors() function
    # source_vector, height_vector, step_vector, width_vector
    src_vector, h_vector, step_vector, w_vector = pyramid_coords

    P0 = src_vector
    P1 = src_vector + h_vector + step_vector + w_vector
    P2 = src_vector + h_vector + step_vector - w_vector
    P3 = src_vector + h_vector - step_vector - w_vector
    P4 = src_vector + h_vector - step_vector + w_vector

    # Create the points
    points = vtk.vtkPoints()
    points.InsertNextPoint(P0)
    points.InsertNextPoint(P1)
    points.InsertNextPoint(P2)
    points.InsertNextPoint(P3)
    points.InsertNextPoint(P4)
     
    # Create a pyramid from the points
    pyramid = vtk.vtkPyramid() 

    for i in range(8):
        pyramid.GetPointIds().SetId(i,i)
     
    # Add the pyramid to a cell array
    cells = vtk.vtkCellArray()
    cells.InsertNextCell(pyramid)
     
    # Add the points and pyramid to an unstructured grid
    uGrid = vtk.vtkUnstructuredGrid()
    uGrid.SetPoints(points)
    uGrid.InsertNextCell(pyramid.GetCellType(), pyramid.GetPointIds())
     
    # Create mapper and actor
    mapper = vtk.vtkDataSetMapper()
    mapper.SetInputData(uGrid)
     
    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    actor.GetProperty().SetColor(color)

    return actor

def get_beam_vectors(system):
    
    # Point where X-rays are generated
    source_vector = np.array([0.0, 0.0, system.XRayFanBeam.surface_offset])

    # Vector Pointing from source to max depth of imaging
    height_vector = np.array([0.0, 0.0, -system.XRayFanBeam.surface_offset -\
                                          system.Collimator.max_depth])

    # Vector pointing from source to half of system step size (x direction) 
    step_vector = np.array([system.step_size/2., 0.0, 0.0])

    # Vector pointing from source to half of detector width (y direction)
    angle = system.XRayFanBeam.fan_angle *np.pi/180.

    half_width = (system.XRayFanBeam.surface_offset + \
                  system.Collimator.max_depth) * np.tan(angle)

    width_vector = np.array([0.0, half_width, 0.0])

    return [source_vector, height_vector, step_vector, width_vector]

def create_plane_actor(center, point1, point2, color=(0,0,0)):

    if np.max(color) > 1.0:

        color = (color[0]/255., color[1]/255., color[2]/255.)

    source = vtk.vtkPlaneSource()
    source.SetPoint1(point1)
    source.SetPoint2(point2)
    source.SetCenter(center)
     
    # mapper
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputConnection(source.GetOutputPort())
     
    # actor
    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    actor.GetProperty().SetColor(color)

    return actor