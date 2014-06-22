#  write out a list to the filepath
def export_list(list, filepath):
    fileout = open(filepath, 'w')
    for item in list:
        fileout.write(str(item) + '\n')
    fileout.close()

#  extract an element from a list of lists/arrays and apply
#+ a multiplying factor
def extract_element(list, index, multiplier = 1.0):
    new_list = []
    for item in list:
        new_list.append(item[index] * multiplier)
    return new_list

#  extract consecutive elements from a list of lists/arrays and apply
#+ a multiplying factor
#  a list for appending can also be passed to the function
def extract_elements(thelist, start_index, end_index, new_list = None, \
                     multiplier = 1.0):
    #  import numpy.array incase list entries are lists
    from numpy import array
    #  if a list to append to is not given,
    #+ then create a new list
    if new_list == None: new_list = []
    for item in thelist:
        #  check to see if item is a list
        #+ and convert to array if it is
        if type(item) == type([]): islist = True; item = array(item)
        else: islist = False
        #  extract the elements
        extracted = item[start_index:end_index + 1] * multiplier
        #  convert back to list if necessary
        if islist: extracted = list(extracted)
        #  append the sublist/subarray
        new_list.append(extracted)
    return new_list

#  create one list of elements from a list of lists
def create_one_list(thelist, start_index, end_index, new_list = None, \
                    multiplier = 1.0):
    #  import numpy.array incase list entries are lists
    import numpy
    #  if a list to append to is not given,
    #+ then create a new list
    if new_list == None: new_list = []
    #  check if given a single list, else we have a list of lists
    #+ or a list of arrays
    if type(thelist[0]) != list and type(thelist[0]) != numpy.ndarray:
        #  create an appropriate slice
        item_slice = thelist[start_index:end_index + 1]
        for x in item_slice:
            new_list.append(x * multiplier)
    else:
        for item in thelist:
            #  create a slice
            item_slice = item[start_index:end_index + 1]
            for x in item_slice:
                new_list.append(x * multiplier)
    return new_list

#  convert r,v to an initial conditions list for use in the
#+ flow.py module
def rv2ic(r, v, mu, STM = False):
    #  import statements
    from numpy import zeros, eye
    #  if converting a planar problem
    if len(r) == 2:
        if STM: ic = zeros(21)
        else: ic = zeros(5)
        ic[0] = r[0]
        ic[1] = r[1]
        ic[2] = v[0]
        ic[3] = v[1]
        ic[4] = mu
        #  append the identity matrix
        if STM: ic[5:] = eye(4).reshape(16,)
    #  if converting a spatial problem
    elif len(r) == 3:
        if STM: ic = zeros(43)
        else: ic = zeros(7)
        ic[0] = r[0]
        ic[1] = r[1]
        ic[2] = r[2]
        ic[3] = v[0]
        ic[4] = v[1]
        ic[5] = v[2]
        ic[6] = mu
        #  append the identity matrix
        if STM: ic[7:] = eye(6).reshape(36,)
    #  return
    return ic

