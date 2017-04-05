import numpy as np

#file_path = '../mesh_6B_particles'

class CartesianMeshTally(object):
    """
    Mesh Tally that uses x y z coordinates
    (as opposd to cylindrical coordinates)

    Initialize it by passing a file_path to a mesh tally output file
    """

    def __init__(self, file_path):
        '''
        num_histories

        tally_number

        x_bin_boundaries

        y_bin_boundaries

        z_bin_boundaries

        energy_bin_boundaries

        matrix  ( shape = (energy, x, z, y) )

        matrix_total  ( shape = (x,z,y) )

        energy

        X 

        Y

        Z

        flux

        relative_error
        '''

        f = open(file_path, 'r')

        lines = f.readlines()

        f.close()

        lines = [line.strip('\n') for line in lines]

        # grab meta data for the mesh tally

        for line in lines:

            if 'Number of histories' in line:

                self.num_histories = float(line.split()[-1])

            if 'Mesh Tally Number' in line:

                self.tally_number = line.split()[-1]

            if 'X direction:' in line:

                self.x_bin_boundaries = map(float, line.split()[2:])

            if 'Y direction:' in line:

                self.y_bin_boundaries = map(float, line.split()[2:])

            if 'Z direction:' in line:

                self.z_bin_boundaries = map(float, line.split()[2:])

            if 'Energy bin boundaries:' in line:

                self.energy_bin_boundaries = map(float, line.split()[3:])

            if 'Energy' in line and 'X' in line:

                start_pos = lines.index(line) + 1

        # Start grabing all the data for the mesh tally

        energy, X, Y, Z, flux, rel_err = [], [], [], [], [], []

        for line in lines[start_pos:]:

            # If there were energy bins, then there will be a "Total" over all
            #   energies at the end

            if 'Total' in line:

                total = True

                x,y,z,r,re = map(float, line.split()[1:])

                e = 1.00e36

            else: 

                total = False

                e,x,y,z,r,re = map(float, line.split())

            energy.append(e)

            X.append(x)

            Y.append(y)

            Z.append(z)

            flux.append(r)

            rel_err.append(re)

        # Get dimensions for the shape of teh arrays

        if total:

            size_E = len(self.energy_bin_boundaries)

        else:

            size_E = len(self.energy_bin_boundaries) - 1

        size_x = len(self.x_bin_boundaries) - 1

        size_y = len(self.y_bin_boundaries) - 1

        size_z = len(self.z_bin_boundaries) - 1

        shape = [size_E, size_x, size_y, size_z]
        # Create  Arrays

        self.energy = np.flipud(np.array(energy).reshape(shape).transpose())

        self.X = np.flipud(np.array(X).reshape(shape).transpose())

        self.Y = np.flipud(np.array(Y).reshape(shape).transpose())

        self.Z = np.flipud(np.array(Z).reshape(shape).transpose())

        self.flux = np.flipud(np.array(flux).reshape(shape).transpose())

        self.relative_error = \
                        np.flipud(np.array(rel_err).reshape(shape).transpose())
