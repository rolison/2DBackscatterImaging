This is a Python package intended to be used for the design and analysis of 
a 2D Backscatter Imaging system. This system is comprised of a fan beam X-ray 
source, a 2D pixelated radiation detector array, and a collimation grid. It is 
hypothesized that this design will allow for 3D (or at least depth sensitive) 
imaging of an object by measuring the X-rays that are scattered inside the 
object and back to the detector array. This technique requires only single sided 
access to said object, in contrast to CT imaging which requires full 360 degree 
access to the object. 

There will be 3 main capabilities for this package:

1) Easily and quickly design a 2D Backscatter imaging system by 
   just inputting simple system parameters (or use the default values). The code 
   will then calculate the required dimensions for the collimation grid in order
   to achieve the desired image resolution. 

2) visualize the system that was designed inside python 
   (using matplotlib and vtk libraries)

3) Export the designed system into a format that can be used 
   inside an MCNP6 input file. This will require the user to have access to 
   MCNP6, the export controlled Monte Carlo Radiation Transport code.

__What's New__: A simple tally_reader.py has been added that can read Cartesian
cooridnate Meshtallies that are output by MCNP simulations and store it in an
object. Can use methods like imshow() from matplotlib.pyplot to view 
backscatter images. image_utils.py is still under construction but contains handy functions for image analysis. 

__Current (Posible) Bug__: after closing the window when calling the system.
visualize() function, I recevie a segmentation fault 11 code. This may be due to 
my own python installation and other system properties so I do not know if other 
users will face this problem.

# Examples on How to Use:

## Create a backscatter system from scratch:

```
from backscatter2D import system_design
          
system = system_design.BackscatterSystem()
```

## Or build components separately and put them together later

You can either add your own settings for each component or just use the
defaults. Use the help() command for documentation on the settings you 
can give each component.

```
detector = system_design.Detector(angle=50.)

collimator = system_design.UniformCollimator(length=6.5)

collimator.construct_grid(detector)

system = system_design.BackscatterSystem(beam_input=None, det_input=detector, 
                                                  collimator_input=collimator)

```
This sets the step size the system will move when scanning across an object

`system.step_size = .5 `

## Visualize your newly created system

```
from backscatter2D import visualize

visualize.system(system)
```

You could also view a 2D plot of the detector and collimator.
Currently it creates this plot at an angle, since I'm unsure how to best
rotate all the shapes correctly inside matplotlib library. 

```
from backscatter2D import matplotting

figure, axes = matplotting.xz_system(system)
```

## Using an X-ray Spectrum and visualization

When creating an X-ray Fan Beam source, you can provide a file that contains
spectrum information and load it into the object. Open the 150kVp.spectrum file 
inside example_spectra/ in order to see what spectrum files should look like. 

```
xray = system_design.XRayFanBeam(energy='example_spectra/150kVp.spectrum')

matplotting.xray_spectrum(xray)
```

Or

```
system.XRayFanBeam = xray

matplotting.xray_spectrum(system)
```

If a file isn't provided, the XRayFanBeam class will default to a mono-energetic
source.

## Exporting your system to an MCNP input file

Just provide the system you've built as well as a file name/path you want it
to export the system to. 

```
from backscatter2D import mcnp

mcnp.create_input_deck(system, file_path='./sample_dir/sample_input.i')
```

If you have MCNP, you can now go plot the input file using MCNP to see that it worked.

