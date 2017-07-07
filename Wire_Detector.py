# -*- coding: utf-8 -*-
"""
Created on Sun Jun 04 20:04:49 2017

@author: Paige
"""

from scipy import ndimage, misc
import matplotlib.pyplot as plt
import numpy as np

### Functions -----------------------------------------------------------------

def filter_pic(picture, sigma, t2, scale, tilt):
    pic = ndimage.imread(picture, mode  = 'L', flatten = True)                     # import picture
    pic = ndimage.gaussian_filter(pic, sigma)                                       # apply Gaussian filter
                                 
    pic_grad = np.gradient(pic)                                                     # take gradient          
    Gx = pic_grad[0]
    Gy = pic_grad[1]
    mag = np.sqrt(Gx**2 + Gy**2)                                                    # magnitude of gradient
    #theta = np.arctan(Gy, Gx)                                                       # orientation of gradient
    
    thresh2 = np.empty(np.shape(mag))
    
    for i in range(np.shape(mag)[0]):                                               # double threshold
        for j in range(np.shape(mag)[1]):
            if mag[i, j] < t2:
                thresh2[i, j] = 0
            else:
                thresh2[i, j] = mag[i, j]
    return thresh2*intensity

def maximum(picture, row, column, L_or_R):
    if L_or_R == "Right":
        while picture[row, column] <= picture[row, column + 1]:
            column += 1
    else:
        while picture[row, column] <= picture[row, column - 1]:
            column -= 1
    return column

def outline(picture):
    outline_pic = plt.figure(3)
    outline_pic_ax = outline_pic.add_subplot(111)
    outline_pic_ax.imshow(picture)
    outline_pic_ax.axis("off")
    outline_pic_ax.scatter(r, h, color = 'r')
    outline_pic_ax.scatter(l, h, color = 'g')
    plt.draw()
    
def dh(scale, tilt, file_name):
    split = np.where(height == height[0] + 1)
    reord_height = list(np.flipud(height[split[0][0]:])) + \
                       list(height[:split[0][0]])
    reord_diam = list(np.flipud(diam[split[0][0]:])) + list(diam[:split[0][0]])
    reord_height = np.array(reord_height)      
    reord_diam = np.array(reord_diam)          
    D = reord_diam/scale
    H = (reord_height/scale)/np.sin(np.radians(tilt))
    for i in range(len(D)):
        print D[i], '\t', H[i]    
    with open(file_name, 'w') as f:
        for i, j in zip(D, H):
            f.write("{}\t{}".format(i, j))
            f.write('\n')
    return

def subpixel()
    
def click(event):
    #print event.xdata, event.ydata
    #print int(event.xdata), int(event.ydata)
    n_orig =  int(event.xdata)
    m_orig = int(event.ydata)                                                   # set original coords
   
    rightx = []                                                                 # lists to append edges
    leftx = []
    y = []
    
    m = m_orig
    n = n_orig
    value = pic_array[m, n]  
    #print value                                                                 # starting value (should be 0)

    while value == 0:                                                           # check up
        y.append(m)
        while value == 0:                                                       # check right
            n += 1
            value = pic_array[m, n]
        n = maximum(pic_array, m, n, "Right")
        rightx.append(n)
        holdr = n
        n = n_orig
        value = pic_array[m, n]
        
        while value == 0:                                                       # check left
            n -= 1
            value = pic_array[m, n]
        n = maximum(pic_array, m, n, "Left")
        leftx.append(n)
        holdl = n                                                               # centre x-coord (m)
        n = int((holdr + holdl)/2)
        m += 1
        value = pic_array[m, n]
    
    m = m_orig - 1                                                              # reset to initial coords (-1 in y direction)
    n = n_orig
    value = pic_array[m, n]

    while value == 0:                                                           # check up
        y.append(m)
        while value == 0:                                                       # check right
            n += 1
            value = pic_array[m, n]
        n = maximum(pic_array, m, n, "Right")
        rightx.append(n)
        holdr = n
        n = n_orig
        value = pic_array[m, n]
        while value == 0:                                                       # check left
            n -= 1
            value = pic_array[m, n]
        n = maximum(pic_array, m, n, "Left")
        leftx.append(n)
        holdl = n                                                               # centre x-coord (m)
        n = int((holdr + holdl)/2)
        m -= 1
        value = pic_array[m, n]
    print "Done"
        
    global r, l, h, diam, height
    r = np.array(rightx)
    l = np.array(leftx)
    h = np.array(y)    
    diam = r - l
    height = max(h) - h             
    return

### Parameters ----------------------------------------------------------------
folder = "FOLDER_LOCATION"

picture = "\\IMAGE_NAME.ext" 
output = "OUTPUT_FILE_LOCATION"
img_name = folder + picture                                                     # picture name
sigma = 1                                                                      # for Gaussian filter
#theta = -5
tilt = 20
scale = 0.107                                                                         # lower threshold
t = 6
intensity = 1                                                                       # upper threshold
theta = 0
epsilon = 1e-4
### Main ----------------------------------------------------------------------


if __name__ == "__main__":
    pic_array = filter_pic(img_name, sigma, t, scale, tilt)
    pic_array = misc.imrotate(pic_array, theta)
    thresh_pic = plt.figure(1)
    thresh_pic_ax = thresh_pic.add_subplot(111)
    thresh_pic_ax.imshow(pic_array, interpolation = "None")
    thresh_pic_ax.axis("off")

    cid = thresh_pic.canvas.mpl_connect('button_press_event', click)                # grab a point