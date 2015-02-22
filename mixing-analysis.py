"""
analysis.py
===========

Written by Scott Thomas 2012-02-01
Made available online 2015-02
See `LICENSE` for terms of re-use

Simple Python script used to test mixing along a microfluidic channel.

Inputs
------

Uses line profiles input from ImageJ. An option is also included to use data
from PMCapture, which is software from the fluorescence microscope. However,
this was never tested, as the fluorescence-based approach was a Plan B intended
for if the contrast-based approach did not work.

The line profiles should be provided in subfolders, one for each different type
of channel to be tested. File names should be as follows: `(x)mm_n.hst` where x
is the distance from the Y-junction in mm and n is the number of the profile at
that point. For example, to use 4 profiles taken at the 5mm mark, provide 4
ImageJ files named `5mm_1.hst`, `5mm_2.hst`, `5mm_3.hst`, `5mm_4.hst`. These
will be median-averaged to give a final profile.

Example folder structure
------------------------

The folder structure below would be used to compare two channels which had 4
measurements taken at regular 5mm intervals.

```
- root_folder/
  - analysis.py
  - channel1/
    - 0mm_1.hst
    - 0mm_2.hst
    - 0mm_3.hst
    - 0mm_4.hst
    - 5mm_1.hst
    - 5mm_2.hst
    - 5mm_3.hst
    - 5mm_4.hst
    - 10mm_1.hst
    - ...
  - channel2/
    - 0mm_1.hst
    - 0mm_2.hst
    - ...
```

Outputs
-------

After running the script, the median-averaged profiles will be presented to you
on screen. Click twice on each plot: once to define the bottom-left point of
the steep middle section, and the second time to define the top-right point.
Then close the plot. The gradient will be written to the file `results.txt` in
the channel folder.

`results.txt` contains one line for each measurement point with the measurement
point (e.g. 5 for 5 mm) followed by the calculated gradient. Simple! If you
want to modify this script to add more info, look at the end of the
`process_channel` function.
"""

from pylab import *
import os


def list_files(path='./'):
    """
    Makes a list of files in the currect (or specified) directory
    Arguments:      path = path to directory (default is current)
    Returns:        list of files in directory\
    """

    file_list = []
    directory_list = os.listdir(path)
    for filename in directory_list:
        file_list.append(filename)

    return file_list


def find_profile_files(path='./'):
    """
    Finds files with the extension .hst in the specified directory. Returns
    them as a list.
    """
    file_list = list_files(path)
    for filename in file_list[:]:
        ext = os.path.splitext(filename)[1]
        if ext != '.hst':
            file_list.remove(filename)

    return file_list


def determine_data_type(filename):
    """
    Figures out whether a given data file was generated using PMCapture
    (fluorescence microscope software) or ImageJ.
    Arguments:      filename
    Returns:        'FM' if from the fluorescence microscope, or 'IJ' if
    otherwise
    """

    file = open(filename, 'r')
    if file.readline() == '# PMCapture Pro Profile Data\n':
        data_type = 'FM'
    else:
        data_type = 'IJ'
    file.close()
    return data_type


def process_images(datalist, num_points=None, smoothing=False):
    """
    Import and process a list of cross-section data files.

    Arguments:
        datalist = list containing filenames of cross-sections from the
                   microscope
        num_points = number of points used (default: 1 pt = 1 pix on largest
                     cross-section)
        smoothing = whether to smooth the data or not (default not)

    Returns:
        positional coordinate array (range 0 - 1)
        averaged greyscale value array
        number of datapoints used
    """

    # Import the data files and extract relevant columns
    data = []
    len_slices = len(datalist)
    for i in range(len_slices):
        data.append([[], []])
        source = determine_data_type(datalist[i])
        datafile = open(datalist[i])
        if source == 'FM':
            col1 = 3
            col2 = 4
        if source == 'IJ':
            col1 = 0
            col2 = 1
        for line in datafile:
            splitline = line.split()
            if splitline != [] and splitline[0] != "#":
                data[i][0].append(float(splitline[col1]))
                data[i][1].append(float(splitline[col2]))
        datafile.close()

    positions = []
    values = []
    normcoords = []
    len_section = []
    for i in range(len_slices):
        # Find channel edges and remove values outside them
        position = array(data[i][0])
        value = array(data[i][1])
        start = float(position[value == min(value[:len(position)/2])])
        stop = float(position[value == min(value[len(position)/2:])])
        value[(position < start) | (position > stop)] = 0
        position[(position < start) | (position > stop)] = 0
        value = trim_zeros(value)
        position = trim_zeros(position)
        if num_points is None:
            len_section.append(size(position))

        # Normalise the distance across the channels on a 0-1 scale
        normcoord = (position - min(position))/(max(position)-min(position))
        positions.append(position)
        values.append(value)
        normcoords.append(normcoord)

    # figure out the number of points in the largest channel
    if num_points is None:
        num_points = max(len_section)

    # write arrays of equal size
    x = [n/float(num_points) for n in range(num_points)]
    ys = []
    for i in range(len_slices):
        y = []
        for n in range(num_points):
            y.append(linear_interpolate(normcoords[i], values[i], x[n]))
        ys.append(y)
    # Median-average all the cross-sections
    ys = array(ys)
    average = median(ys, axis=0)

    if smoothing:
        average = smooth(average[:], int(num_points/100))

    return x, average, num_points


def linear_interpolate(x, y, x_desired):
    """
    Given two lists or arrays x and y, linearly interpolates between adjacent
    x-values to give an estimated y-value for any particular x

    Arguments:
        x = 1-D array or list of floats
        y = 1-D array or list of floats containing values matching the
            positions in x
        x_desired = the x value to interpolate at
    Returns:
        interpolated y-value
    """

    index = 0
    while x[index] <= x_desired:
        index = index + 1

    x_left = x[index-1]
    x_right = x[index]

    y_interp = ((y[index-1] * (x_right - x_desired))
                + (y[index] * (x_desired - x_left))) / (x_right - x_left)

    return y_interp


def smooth(data, degree):
    """
    Smooths data using a moving triangular window. Leaves copies of data at the
    end.

    Arguments:
        data = list containing the data, assumed to be equally spaced
        degree = amount of smoothing to perform (size of triangle)

    Returns:
        smoothed array
    """

    triangle = array(range(degree) + [degree] + range(degree)[::-1]) + 1
    smoothed = []
    for i in range(degree, len(data) - degree*2):
        point = data[i:i+len(triangle)]*triangle
        smoothed.append(sum(point) / sum(triangle))
    smoothed = [smoothed[0]]*(degree + degree/2) + smoothed
    while len(smoothed) < len(data):
        smoothed.append(smoothed[-1])
    return smoothed


# This was written to trial a backgrouund subtraction procedure, which would
# allow systematic trends in the images (from slides, microscopes, or whatever)
# to be reduced by taking a picture of the channel before it was filled with
# liquid. However, this approach did not end up being used, as it would be
# necessary to either modify the images before using ImageJ on them, or to take
# cross-sections in exactly the same place each time and subtract the profiles
# from each other. Ultimately, comparing contrast gradients to the Y-junction
# gave the same effect.
def background_subtract(datalist, backgroundlist, source='IJ'):
    """ Given two lists of data files, subtract the second from the first
    Arguments:  datalist = list of files containing sample cross sections
                backgroundlist = list of files containing blank cross sections
                Source = 'FM' fluorescence microscope or 'IJ' ImageJ (default)
    Returns:    normalised position coordinate array
                background-subtracted grayscale value array
                number of datapoints used"""

    # Process the sample and background images
    x, y, n = process_images(sample, source)
    xb, yb, nb = process_images(background, source)

    # Truncate to the smaller one
    if n > nb:
        n = nb
        x = x[:nb]
        y = y[:nb]
    if n < nb:
        nb = n
        xb = xb[:n]
        yb = yb[:n]
    ys = (y - yb)/yb  # remove background effects

    return x, ys, n


def process_channel(channel_name):
    """
    Processes a single folder of files set up as described in the intro,
    writing the result to `results.txt`.

    Example:
        process_channel("channel1") will process the channel in the subfolder
        `channel1/`.
    """

    file_list = find_profile_files('./' + channel_name)
    sizes = set([f.split('mm_')[0] for f in file_list])
    groups = []
    for s in sizes:
        # scans through the different lengths (first part of filename) and
        # makes a list linking these to the repeats (second part of filename).
        # the end result is a list of lists containing grouped information
        # about the files to process
        groups.append([s, [f.split(s + 'mm_')[1][:-4]
                      for f in file_list
                      if s == f.split('mm_')[0]]])

    xs = {}
    ys = {}

    for g in groups:
        num_samples = len(g[1])
        prefix = './' + channel_name + '/' + g[0] + 'mm_'
        image_list = [prefix + g[1][i] + '.hst' for i in range(num_samples)]
        xout, yout, nout = process_images(image_list)
        xs[g[0]] = xout
        ys[g[0]] = yout

    xlist = []
    ylist = []
    for i in groups:
        # Plot the cross-section
        fig = figure()
        ax = fig.add_subplot(111)
        ax.plot(xs[i[0]], ys[i[0]], 'k-', linewidth=2)
        xlabel('Position across channel')
        ylabel('Pixel value')

        # Define event to handle clicks on the plot
        def onclick(event):
            xlist.append(event.xdata)
            ylist.append(event.ydata)
        # Connect event manager to plot and show it
        cid = fig.canvas.mpl_connect('button_press_event', onclick)
        show()

    # Calculate gradients using the click positions
    gradient = []
    for i in range(0, len(xlist), 2):
        gradient.append([float(groups[i/2][0]),
                        (ylist[i+1] - ylist[i])/(xlist[i+1] - xlist[i])])
    savetxt('./' + channel_name + '/' + 'results.txt', gradient)

    return
