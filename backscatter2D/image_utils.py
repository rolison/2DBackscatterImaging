import matplotlib.pyplot as plt
import numpy as np

def getVminVmax(image, num_std=3):
    '''
    Used for plotting images using the imshow() method from matplotlib.pyplot

    Returns a tuple containing values to be used for the vmin and vmax
    arguments inside imshow()

    Pass this method a 2D matrix ('image') and it will find the mean value
    for elements in the matrix that are either num_std above the average
    or below it.

    step 1) find mean value and standard deviation of elements in the entire
            2D array.

    step 2) find all elements that are the requested number of standard
            deviations (num_std) away from the mean. One array is below the
            mean (vmin values), the other array is above the 
            mean (vmax values).

    step 3) take the mean of those two new arrays, resulting in a single
            value for each vmin and vmax.

    step 4) check to see if either of those values are null (due to there being
            no elements that were above or below the requested number of 
            standard deviations). If it is null, just set the value to being 
            the maximum or minimum value in the original original array. 
    '''

    vmin = np.mean(image[image < image.mean() - num_std*image.std()])
    vmax = np.mean(image[image > image.mean() + num_std*image.std()])

    # In case the above search fails (i.e. there are no pixels with values
    #   greater than the requested number of standard deviations away from the
    #   mean)
    vmin = image.min() if np.isnan(vmin) else vmin 
    vmax = image.max() if np.isnan(vmax) else vmax

    return (vmin, vmax)

def filterFromImage(image_feature, image_background, num_std):
    '''
    Take the mean and std of image_background. Filter
    image_feature into 1's and 0's, based on whether the pixel intensity is
    greater than the mean + num_std*std's (=1) or is less than (=0)
    '''
    pass

    filtered_image = np.zeros(image_feature.shape)

    filtered_image[image_feature > image_background.mean() + 
                                   num_std*image_background.std()] = 1

    return filtered_image

def filterByRowFromImage(image_feature, image_background, num_std):
    '''
    For each row, take the mean and std of image_background. Filter
    image_feature into 1's and 0's, based on whether the pixel intensity is
    greater than the mean + num_std*std's (=1) or is less than (=0)
    '''

    filtered_image = np.zeros(image_feature.shape)

    for filtered_row, row_b, row_f in zip(filtered_image, image_background,
                                                              image_feature):
 
        if row_f.mean() - row_b.mean() > 0:
            
            filtered_row[row_f > row_b.mean() + num_std*row_b.std()] = 1

        else:

            filtered_row[row_f < row_b.mean() - num_std*row_b.std()] = 1

    return filtered_image

def filterByColumnFromImage(image_feature, image_background, num_std):
    '''
    For each column, take the mean and std of image_background. Filter
    image_feature into 1's and 0's, based on whether the pixel intensity is
    greater than the mean + num_std*std's (=1) or is less than (=0)
    '''

    im_f = image_feature.transpose()

    im_b = image_background.transpose()

    filtered_image = np.zeros(im_f.shape)

    for filtered_column, col_b, col_f in zip(filtered_image, im_b, im_f):

        if col_f.mean() - col_b.mean() > 0:
            
            filtered_column[col_f > col_b.mean() + num_std*col_b.std()] = 1

        else:

            filtered_column[col_f < col_b.mean() - num_std*col_b.std()] = 1


    return filtered_image.transpose()

def filterByRowFromSelf(image, num_std):
    '''
    For each row, take the mean and std of image. Filter
    itself into 1's and 0's, based on whether the pixel intensity is
    greater than the mean + num_std*std's (=1) or is less than (=0)
    '''

    pass
    
    filtered_image = np.zeros(image.shape)

    for filtered_row, im_row in zip(filtered_image, image):

        filtered_row[im_row > im_row.mean() + num_std*im_row.std()] = 1

    return filtered_image

def filterByColumnFromSelf(image, num_std):
    '''
    For each column, take the mean and std of image. Filter
    itself into 1's and 0's, based on whether the pixel intensity is
    greater than the mean + num_std*std's (=1) or is less than (=0)
    '''

    pass
    
    im_transp = image.transpose()

    filtered_image = np.zeros(im_transp.shape)

    for filtered_column, im_col in zip(filtered_image, im_transp):

        filtered_column[im_col > im_col.mean() + num_std*im_col.std()] = 1

    return filtered_image.transpose()

def filterFromSelf(image, num_std):
    '''
    Take the mean and std of image. Filter itself into 1's and 0's, 
    based on whether the pixel intensity is greater than the 
    mean + num_std*std's (=1) or is less than (=0)
    '''
    pass

def rescaleImage(matrix, range_values=255):
    '''
    pass it a 2D+ array and it will convert it's values into a linear range
    between 0 and range_values, for each array in the first dimension of the
    passed array (matrix). Each value is rounded up.

    Default is 255, after the 255 grayscale values in an image.

    Equation: scaled_element = scaling_factor*(element - min)/(max - min)
    '''
    new_matrix = np.copy(matrix)

    for i, depth_plane in enumerate(new_matrix):
        depth_plane=(depth_plane - depth_plane.min())/(depth_plane.max()-depth_plane.min())*range_values
        new_matrix[i] = depth_plane.round()

    return new_matrix

def plotFilterByRowFromImage(image_feature, image_background, num_std, fs=20, 
                                                              row=None):
    '''
    Makes a plot for each row, visually showing the filterByRowFromImage
    method.

    Either plot every row or specify a row to view.
    '''

    num_row, len_row = image_feature.shape

    if row is None:

        rows = range(num_row)

    else:

        rows = [row]

    for row in rows:

        fig = plt.figure('Filtering Image for Row ' + str(row))

        fig.clear()

        plt.plot(xrange(len_row), image_background[row], 'g',
                 lw=2, label='Background Image')

        plt.plot(xrange(len_row), image_feature[row], 'b',
                 lw=2, label='Image with Feature')

        mean = np.mean(image_background[row])

        std = np.std(image_background[row])

        plt.plot(xrange(len_row), [mean]*len_row, 'k',
                 lw=3, label='Mean Intensity of Background')

        plt.plot(xrange(len_row), [mean + num_std*std]*len_row, 'r', alpha=1,
                 lw=2, label=str(num_std)+' Standard Deviation of Background')

        plt.plot(xrange(len_row), [mean - num_std*std]*len_row, 'r', alpha=1,
                 lw=2)

        plt.xlabel('Pixel Location', fontsize=fs)
        plt.ylabel('Pixel Intensity (AU)', fontsize=fs)
        plt.yticks(fontsize=fs)
        plt.xticks(xrange(len_row), fontsize=fs)
        plt.xlim([0,len_row-1])

        ax = fig.gca()
        ax.legend(loc='best', fontsize=fs)

        ax.fill_between(xrange(len_row), [mean - num_std*std]*len_row, 
                        [mean + num_std*std]*len_row, 
                        facecolor='red', alpha=0.1)


    plt.show()

