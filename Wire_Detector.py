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
    pic = ndimage.gaussian_filter(pic, sigma)                                      # apply Gaussian filter
                                 
    pic_grad = np.gradient(pic)                                                    # take gradient          
    Gx = pic_grad[0]
    Gy = pic_grad[1]
    mag = np.sqrt(Gx**2 + Gy**2)                                                   # magnitude of gradient
    #theta = np.arctan(Gy, Gx)                                                     # orientation of gradient
    
    over = (mag > t2).astype(int)
    thresh = mag*over 
    
    return thresh*intensity

def maximum(picture, row, column, L_or_R):
    if L_or_R == "Right":
        while picture[row, column] <= picture[row, column + 1]:
            column += 1
    else:
        while picture[row, column] <= picture[row, column - 1]:
            column -= 1
    return column

def subpixel(picture, row, column, L_or_R, max_col):
    g = []
    d = []
    while picture[row, column] != 0:
        g.append(picture[row, column])
        d.append(column - max_col)
        if L_or_R == "Right":
            column += 1
        else:
            column -= 1
    g = np.array(g)
    #print "g:", g
    d = np.array(d)
    #print "d:", d
    delta = sum(g*d)/sum(g)
    #print "delta:", delta
    return delta   

#def compare
    

def outline(picture):
    if mode == 1:
        outline_pic = plt.figure(3)
        outline_pic_ax = outline_pic.add_subplot(111)
        outline_pic_ax.imshow(picture)
        outline_pic_ax.axis("off")
        outline_pic_ax.scatter(r, h, color = 'r', marker = '.')
        outline_pic_ax.scatter(l, h, color = 'g', marker  = '.')
        plt.draw()
    elif mode == 2:
        outline_pic = plt.figure(3)
        outline_pic_ax = outline_pic.add_subplot(111)
        outline_pic_ax.imshow(picture)
        outline_pic_ax.axis("off")
        if press == 1:
            outline_pic_ax.scatter(upperright, mtip, color = 'y', marker = '.')
            outline_pic_ax.scatter(upperleft, mtip, color = 'c', marker  = '.')
        elif press == 2:
            outline_pic_ax.scatter(lowerright, mbase, color = 'y', marker = '.')
            outline_pic_ax.scatter(lowerleft, mbase, color = 'c', marker  = '.')
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
        
def wire_trace(m_orig, n_orig):
    rightx = []                                                                 # lists to append edges
    leftx = []
    y = []
    
    m = m_orig
    n = n_orig
    value = pic_array[m, n] 
    print value                                                                # starting value (should be 0)

    while value == 0:                                                           # check up
        y.append(m)
        while value == 0:                                                       # check right
            n += 1
            value = pic_array[m, n]
        n_nonzero = n
        n_max = maximum(pic_array, m, n, "Right")
        subpixel_n = subpixel(pic_array, m, n_nonzero, "Right", n_max) + n_max
        rightx.append(subpixel_n)
        holdr = n_nonzero
        n = n_orig
        value = pic_array[m, n]
        while value == 0:                                                       # check left
            n -= 1
            value = pic_array[m, n]
        n_nonzero = n
        n_max = maximum(pic_array, m, n, "Left")
        subpixel_n = subpixel(pic_array, m, n_nonzero, "Left", n_max) + n_max
        #print "Subpixel:", subpixel_n - n_max
        #print "Corrected:", subpixel_n
        leftx.append(subpixel_n)
        holdl = n_nonzero                                                       # centre x-coord (m)
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
        n_nonzero = n
        n_max = maximum(pic_array, m, n, "Right")
        subpixel_n = subpixel(pic_array, m, n_nonzero, "Right", n_max) + n_max
        rightx.append(subpixel_n)
        holdr = n_nonzero
        n = n_orig
        value = pic_array[m, n]
        while value == 0:                                                       # check left
            n -= 1
            value = pic_array[m, n]
        n_nonzero = n
        n_max = maximum(pic_array, m, n, "Left")
        subpixel_n = subpixel(pic_array, m, n_nonzero, "Left", n_max) + n_max
        #print "Subpixel:", subpixel_n - n_max
        #print "Subpixel_n:", subpixel_n
        leftx.append(subpixel_n)
        holdl = n_nonzero                                                       # centre x-coord (m)
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

def wire_singles(m_orig, n_orig, press):
    #global diamtip, diambase, wireheight
    if (press != 1 and press != 2):
        print 'give a valid value for "press"'
    else:   
        m = m_orig
        n = n_orig
        value = pic_array[m, n] 
        print "val:", value                                                                # starting value (should be 0)
        
        while value == 0:                                                       # check right
            n += 1
            value = pic_array[m, n]
        n_nonzero = n
        n_max = maximum(pic_array, m, n, "Right")
        subpixel_n = subpixel(pic_array, m, n_nonzero, "Right", n_max) + n_max
        if press == 1:
            upperright = subpixel_n
            mtip = m
        elif press == 2:
            lowerright = subpixel_n
            print "lowerright", lowerright
            mbase = m
            print "mbase", mbase
            
        n = n_orig
        
        value = pic_array[m, n]
        while value == 0:                                                       # check left
            n -= 1
            value = pic_array[m, n]
        n_nonzero = n
        n_max = maximum(pic_array, m, n, "Left")
        subpixel_n = subpixel(pic_array, m, n_nonzero, "Left", n_max) + n_max
        if press == 1:
            upperleft = subpixel_n
            #print "upperleft", upperleft
            lowerleft = 'Empty'
            lowerright = 'Empty'
            mbase = 'Empty'
            
        elif press == 2:
            lowerleft = subpixel_n
            print "lowerleft", lowerleft
            upperleft = 'Empty'
            upperright = 'Empty'
            mtip = 'Empty'
        
        return (upperleft, upperright, mtip, lowerleft, lowerright, mbase)   
 


       
def click(event):
    global diamtip, diambase, tip_to_base, mtip, mbase, upperright, upperleft,\
            lowerright, lowerleft
    #print event.xdata, event.ydata
    #print int(event.xdata), int(event.ydata)

    n_orig =  int(event.xdata)
    m_orig = int(event.ydata)                                                   # set original coords
    
    if mode == 1:
        wire_trace(m_orig, n_orig)
        outline(unfiltered)
        dh(scale, tilt, output)
    
    elif mode == 2:
        print "press = ", press
        if press == 1:
            results = wire_singles(m_orig, n_orig, press = 1)
            upperleft = results[0]
            #print "check", upperleft
            upperright = results[1]
            #print "check", upperright
            mtip = results[2]
            print "check mtip", mtip
            ###PLOT upperleft, upperright
            diamtip = abs(upperleft - upperright)/scale
            #press +=
            outline(unfiltered)
        elif press == 2:
            results = wire_singles(m_orig, n_orig, press = 2)
            lowerleft = results[3]
            lowerright = results[4]
            mbase = results[5]
            outline(unfiltered)
            #print "check mbase", mbase

            ### PLOT lowerleft, lowerright
            diambase = abs(lowerleft - lowerright)/scale
            #print "check", diambase
            tip_to_base = (abs(mtip - mbase)/scale)/(np.sin(np.radians(tilt)))
            #print "check tip_to_base"
            print 'Diameter, Height'
            print diamtip, tip_to_base
            print diambase           
    else:
        print 'give a valid "mode"'
        
    return
  
### Parameters ----------------------------------------------------------------
folder = "C:\Users\Paige\Documents\MASc McMaster\Samples\Sample_1909\Images\Jul7"
picture = "\\1909_360D1_20deg_20x_2.tif" 

output = "C:\Users\Paige\Documents\MASc McMaster\Samples\Sample_1909\DH.txt"
img_name = folder + picture               
                                      # picture name                                      
unfiltered = ndimage.imread(img_name)
sigma = 2                                                                    # for Gaussian filter
#theta = -5
tilt = 20
scale = 0.2125                                                                  # threshold
t = 5.2
intensity = 1                                                                   
theta = 0
epsilon = 1e-4
### Main ----------------------------------------------------------------------


if __name__ == "__main__":
    mode = 0
    press = 0
    pic_array = filter_pic(img_name, sigma, t, scale, tilt)
    #pic_array = misc.imrotate(pic_array, theta)
    thresh_pic = plt.figure(1)
    thresh_pic_ax = thresh_pic.add_subplot(111)
    thresh_pic_ax.imshow(pic_array, interpolation = "None")
    thresh_pic_ax.axis("off")
    plt.draw()
    cid = thresh_pic.canvas.mpl_connect('button_press_event', click)            # grab a point