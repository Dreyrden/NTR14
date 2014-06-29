#----------------------------------------------------------#
#---------------------- Plot Module -----------------------#
#----------------------------------------------------------#
"""
    module: plotmodule.py
    date:   2014_05_13
    author: ryne beeson
"""


#  function to decimate_flow flow
def decimate_flow(flow, pts):
    #  import statements
    from math import ceil
    import numpy as np
    
    #  check if flow argument is of type dict, which
    #+ indicates that there is a flow['x'] and flow['y']
    #+ list.
    #+ otherwise the argument should be a list of [[x, y, z],[...],..]
    if type(flow) == dict:
        lx = len(flow['x'])
        if lx > pts:
            xend = flow['x'][-1:][0]
            yend = flow['y'][-1:][0]
            step = int(ceil(lx/float(pts)))
            flow['x'] = flow['x'][0:lx:step]
            flow['y'] = flow['y'][0:lx:step]
            flow['x'].append(xend)
            flow['y'].append(yend)
    else:
        lx = len(flow)
        if lx > pts:
            yend = flow[-1]
            step = int(ceil(lx/float(pts)))
            flow = flow[0:lx:step]
            flow.append(yend)
    #  return the decimated flow argument
    return flow


#  function to generate a list of markers, colors for plotting
def generate_markers(numpts):
    markers = ['o', 'd', 'x', '+']
    colors  = ['b', 'r', 'g', 'y']
    mylist  = []
    count   = 1
    #  run a while-loop until all points have a marker
    #+ and a color associated with them
    while count <= numpts:
        #  select a marker
        for i in range(len(markers)):
            #  select a color
            for j in range(len(colors)):
                if count > numpts: return mylist
                else: count += 1; mylist.append([markers[i], colors[j]])
    
    #  return statement in if-statement above


#  extract the x, y, z cartesian positions from a flow matrix
#+ or from a list containing [x, y, z] lists
def extract_xyz(flow):
    #  check if flow argument is of type dict, which
    #+ indicates that there is a flow['x'] and flow['y']
    #+ list.
    #+ otherwise the argument should be a list of [[x, y, z],[...],..]
    x = []; y = []; z = []
    if type(flow) == dict:
        for row in flow['y']:
            x.append(row[0])
            y.append(row[1])
            z.append(row[2])
    else:
        for row in flow:
            x.append(row[0])
            y.append(row[1])
            z.append(row[2])
    #  return x, y, z
    return x, y, z


def histogram_plot(data, num_bins, \
                   fignum = 1, \
                   title = '', \
                   xlabel = '', \
                   ylabel = '', \
                   interactive = True, \
                   savefile = False, \
                   savepath = "/Users/rynebeeson/Desktop/", \
                   image_name = "histogram_plot.png"):
    
    #  import statments
    import numpy as np
    import matplotlib.pyplot  as plt
    import matplotlib.patches as patches
    import matplotlib.path    as path
    from matplotlib.pylab import savefig
    from time import sleep
    
    fig, ax = plt.subplots(1, 1, num = fignum)
    
    n, bins = np.histogram(data, num_bins)
    
    # get the corners of the rectangles for the histogram
    left = np.array(bins[:-1])
    right = np.array(bins[1:])
    bottom = np.zeros(len(left))
    top = bottom + n
    
    
    # we need a (numrects x numsides x 2) numpy array for the path helper
    # function to build a compound path
    XY = np.array([[left,left,right,right], [bottom,top,top,bottom]]).T
    
    # get the Path object
    barpath = path.Path.make_compound_path_from_polys(XY)
    
    # make a patch out of it
    patch = patches.PathPatch(barpath, facecolor='blue', \
                              edgecolor='gray', alpha=0.8)
    ax.add_patch(patch)
    
    # update the view limits
    ax.set_xlim(left[0], right[-1])
    ax.set_ylim(bottom.min(), top.max())

    # set title and axes labels
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    
    plt.grid()
    if interactive: plt.show()
    else: savefig(savepath + image_name)


#  LinePlot Class is a default class for
#+ creating 2D Plots
class LinePlot:
    'LinePlot class documentation string'

    #  initialization method
    def __init__(self, number):
        #  import statements
        import matplotlib.pyplot as plt
        self.fig_number = number
        self.axisbg = '1.0'
        self.fig = plt.figure(number)
        self.ax  = self.fig.add_subplot(111, axisbg = self.axisbg)
        
        #  additional plotting attributes
        self.xlabel = ''
        self.ylabel = ''
        self.title  = ''
        #  x, y axis containers for plot scaling
        self.xaxis = []
        self.yaxis = []
        #  x, y containers for plotting
        self.x = []
        self.y = []
        #  initialize min and max for setting axis
        self.lmin = 0.0
        self.lmax = 0.0

    #  2d plotting method
    def plot(self, color = 'b', marker = ''):
        # plot the flow
        self.ax.plot(self.x, self.y, c = color, marker = marker)
        #  find the min and max of the plot data
        lmins = [min(self.x), min(self.y)]
        lmaxs = [max(self.x), max(self.y)]
        lmin  = min(lmins)
        lmax  = max(lmaxs)
        #  find the extreme min and max of the plot data
        #+ and save to self attributes if
        #+ max is larger than current self.max or
        #+ min is smaller than current self.min
        if lmin < self.lmin: self.lmin = lmin
        if lmax > self.lmax: self.lmax = lmax

    #  set the axes
    def set_axis(self, aspect = 'equal'):
        #  if aspect is set to equal by default or user
        if aspect == 'equal':
            self.ax.set_xlim3d(self.lmin, self.lmax)
            self.ax.set_ylim3d(self.lmin, self.lmax)

    #  set the aspect ratio
    def set_aspect(self, aspect = 'equal'):
        if aspect == 'equal':
            self.ax.set_aspect(aspect)
        else:
            self.ax.set_aspect(aspect)

    #  set the legend
    def set_legend(self, list = ['']):
        self.ax.legend(list)

    #  draw the plot
    def draw(self):
        #  import statements
        import matplotlib.pyplot as plt
        plt.draw()

    #  draw a grid on plot
    def grid(self):
        #  import statements
        import matplotlib.pyplot as plt
        plt.grid()

    #  show the plot
    def show(self):
        #  import statements
        import matplotlib.pyplot as plt
        plt.show(self.fig_number)







#  ArcPlot Class is a default class for
#+ creating 3D Plots
class ArcPlot:
    'ArcPlot class documentation string'

    #  initialization method
    def __init__(self, number):
        #  import statements
        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d import Axes3D
        self.fig_number = number
        self.axisbg = '1.0'
        self.fig = plt.figure(number)
        self.ax  = self.fig.add_subplot(111, projection = '3d', \
                                        aspect = 'equal', axisbg = self.axisbg)
        #  additional plotting attributes
        self.xlabel = ''
        self.ylabel = ''
        self.zlable = ''
        self.title  = ''
        #  x, y, z axis containers for plot scaling
        self.xaxis = []
        self.yaxis = []
        self.zaxis = []
        #  x, y, z containers for plotting
        self.x = []
        self.y = []
        self.z = []
        #  initialize min and max for setting axis
        self.lmin = 0.0
        self.lmax = 0.0
    
    #  turn the flow into x, y, z lists and
    #+ reduce the number of data points to be plotted
    def flow2xyz(self, flow, number_of_points = 100):
        #  import statements
        from plotmodule import decimate_flow, extract_xyz
        #  decimate the flow
        flow = decimate_flow(flow, number_of_points)
        #  extract x, y, z
        self.x, self.y, self.z = extract_xyz(flow)
    
    #  3d arc plotting method
    def plot(self, color = 'b'):
        # plot the flow
        self.ax.plot(self.x, self.y, self.z, c = color)
        #  find the min and max of the plot data
        lmins = [min(self.x), min(self.y), min(self.z)]
        lmaxs = [max(self.x), max(self.y), max(self.z)]
        lmin  = min(lmins)
        lmax  = max(lmaxs)
        #  find the extreme min and max of the plot data
        #+ and save to self attributes if
        #+ max is larger than current self.max or
        #+ min is smaller than current self.min
        if lmin < self.lmin: self.lmin = lmin
        if lmax > self.lmax: self.lmax = lmax
    
    #  set the axes
    def set_axis(self, aspect = 'equal'):
        #  if aspect is set to equal by default or user
        if aspect == 'equal':
            self.ax.set_xlim3d(self.lmin, self.lmax)
            self.ax.set_ylim3d(self.lmin, self.lmax)
            self.ax.set_zlim3d(self.lmin, self.lmax)
    
    #  set the legend
    def set_legend(self, list = ['']):
        self.ax.legend(list)

    #  draw the arc plot
    def draw(self):
        #  import statements
        import matplotlib.pyplot as plt
        plt.draw()

    #  show the arc plot
    def show(self):
        #  import statements
        import matplotlib.pyplot as plt
        plt.show(self.fig_number)




    
    











#----------------------------------------------------------#
#---------------------- Plot Module -----------------------#
#----------------------------------------------------------#