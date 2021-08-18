# -*- coding: utf-8 -*-
"""
Given an image draw an approximation using the discrete fourier transform.

@author: Dylan Wattles
"""

import numpy as np
import matplotlib.pyplot as plt
from PIL import Image
from skimage import filters
from skimage.color import rgb2gray
from matplotlib import animation

DIST=4
SOBEL_CONST=.2
FRAMES=1000

#get edge points from an image that has been filtered using the sobel operator
def get_points(sobel_image, points):
    
    rows=sobel_image.shape[0]
    cols=sobel_image.shape[1]
    
    for i in range(0,rows):
        for j in range(0,cols):
            if sobel_image[i][j] > SOBEL_CONST:
                points.append((j,i))   
    return()

#remove points that are DIST away from eachother
def prune_points(x0,y0):
    
    i=0
    j=0
    last_x0=x0[len(x0)-1]
    last_y0=y0[len(y0)-1]
    
    while x0[i]!=last_x0 or y0[i]!=last_y0:
        j=i+1
        while x0[j]!=last_x0 or y0[j]!=last_y0:
            if np.sqrt((x0[i]-x0[j])**2+(y0[i]-y0[j])**2)<DIST:
                x0.pop(j)
                y0.pop(j)        
            else:
                j=j+1
        i=i+1
        
    return()

#nearest neighbor approximate solution to TSP, orders points for path
def order_points(x0,y0):
    
    ordered_x0=[]
    ordered_y0=[]
    ordered_x0.append(x0[0])
    ordered_y0.append(y0[0])
    
    unvis_x0=x0[:]
    unvis_y0=y0[:]
    unvis_x0.pop(0)
    unvis_y0.pop(0)
    
    xpos=x0[0]
    ypos=y0[0]
    closest_index=0
    
    while len(unvis_x0)>0:
        min_dist=1000
        i=0
        while i<len(unvis_x0):
            if np.sqrt((xpos-unvis_x0[i])**2+(ypos-unvis_y0[i])**2)<min_dist:
                min_dist=np.sqrt((xpos-unvis_x0[i])**2+(ypos-unvis_y0[i])**2)
                closest_x=unvis_x0[i]
                closest_y=unvis_y0[i]
                closest_index=i
            i=i+1
            
        ordered_x0.append(closest_x)
        ordered_y0.append(closest_y)
        xpos=closest_x
        ypos=closest_y
        unvis_x0.pop(closest_index)
        unvis_y0.pop(closest_index)
        
    x0.clear()
    y0.clear()
    for i in ordered_x0:
        x0.append(i)
    for i in ordered_y0:
        y0.append(i)
        
    return()

#convert to complex plane
def to_complex(x,y):
    
    re=np.array(x)
    im=np.array(y)
    z=re+1j*im
    
    return(z)

#discrete fourier transform using definition normalized by 1/N
def dft(z):
    
    N=len(z)
    c_n=np.zeros(N, dtype=complex)

    for k in range(0,N):
        for n in range(0,N):
            c_n[k]+=(1/N)*z[n]*np.exp(-1j*2*np.pi*k*n/N)
    
    return(c_n)

#open original image
image = np.array(Image.open("seahorse.jpg"))
#plt.figure(figsize = (12, 8))
#plt.title("Original")
#plt.imshow(image)

#convert rgb to gray so it can be used with sobel operator
gray_im = rgb2gray(image)
#plt.figure(figsize = (12, 8))
#plt.title("Gray")
#plt.imshow(gray_im, cmap=plt.cm.gray)

#use sobel operator for edge detection
sobel_im = filters.sobel(gray_im)
#plt.figure(figsize = (12, 8))
#plt.title("Sobel")
#plt.imshow(sobel_im, cmap=plt.cm.gray)

#get points on edge
points=[]
get_points(sobel_im, points)
x=[]
y=[]
for i in points:
    x.append(i[0])
    y.append(i[1])
#plt.figure(figsize = (12, 8))
#plt.scatter(x,y)
#plt.title("Points")
#plt.show()

#center image at (0,0), also flip image to compensate for indexing convention
x_center = np.average(x)
y_center = np.average(y)
x0=[]
y0=[]
for a in x:
    x0.append(a-x_center)
for b in y:
    y0.append(y_center-b)
#plt.figure(figsize = (12, 8))
#plt.scatter(x0,y0)
#plt.title("Centered_Points")
#plt.show()

#remove points that are too close to eachother
prune_points(x0,y0)
#plt.figure(figsize = (12, 8))
#plt.scatter(x0,y0)
#plt.title("Pruned_Points")
#plt.show()

#approximate tour through all points
order_points(x0,y0)
#plt.figure(figsize = (12, 8))
#plt.plot(x0,y0,marker="o")
#plt.title("Ordered_Points")
#plt.show()

#put points on complex plane and take the dft to find fourier coefficients
z_seq = to_complex(x0,y0)
c_n = dft(z_seq)

#use if need greater efficiency
#c_n = scipy.fft.fft(z_seq)/len(z_seq)

#set figure
t=0
fig = plt.figure()
fig.set_dpi(100)
fig.set_size_inches(7, 6.5)

ax = plt.axes(xlim=(min(x0)-100, max(x0)+100), ylim=(min(y0)-100, max(y0)+100))

line, = ax.plot([], [], lw = 2) 
circles=[]
arrows=[]
xdata, ydata = [], []

#intialize animation
def init():
    
    #init is called twice by anim, dont append/add_patch twice
    if(len(circles)==0):
    
        x_pos=0
        y_pos=0
    
        for k in range(0,len(c_n)):
            circles.append(plt.Circle((x_pos,y_pos),abs(c_n[k]),fill=False))
            if(k<len(c_n)/2):
                circles[k].center=(x_pos,y_pos)
                arrows.append(plt.Arrow(x_pos,y_pos,
                                        (c_n[k]*np.exp(2*np.pi*1j*k*t/FRAMES)).real,
                                        (c_n[k]*np.exp(2*np.pi*1j*k*t/FRAMES)).imag))
                x_pos+=(c_n[k]*np.exp(2*np.pi*1j*k*t/FRAMES)).real
                y_pos+=(c_n[k]*np.exp(2*np.pi*1j*k*t/FRAMES)).imag
                
            
            else:
                circles[k].center=(x_pos,y_pos)
                arrows.append(plt.Arrow(x_pos,y_pos,
                                        (c_n[k]*np.exp(-2*np.pi*1j*k*t/FRAMES)).real,
                                        (c_n[k]*np.exp(-2*np.pi*1j*k*t/FRAMES)).imag))
                x_pos+=(c_n[k]*np.exp(-2*np.pi*1j*(len(c_n)-k)*t/FRAMES)).real
                y_pos+=(c_n[k]*np.exp(-2*np.pi*1j*(len(c_n)-k)*t/FRAMES)).imag
            
        for k in range(0,len(c_n)):
            ax.add_patch(arrows[k])
            ax.add_patch(circles[k])
    
    line.set_data([], [])
        
    return []

#update dynamic epicycles, arrows at frame t
def update_animation(t,circles,arrows):
    
    ax.cla()
    ax.set_xlim(min(min(x0),min(y0))-100, max(max(x0),max(y0))+100)
    ax.set_ylim(min(min(x0),min(y0))-100, max(max(x0),max(y0))+100)
    
    line, = ax.plot([], [], lw = 2)
    
    x_pos=0
    y_pos=0
    for k in range(0,len(c_n)):
            
        if(k<len(c_n)/2):
            circles[k].center=(x_pos,y_pos)
            arrows[k]=plt.Arrow(x_pos,y_pos,
                                        (c_n[k]*np.exp(2*np.pi*1j*k*t/FRAMES)).real,
                                        (c_n[k]*np.exp(2*np.pi*1j*k*t/FRAMES)).imag,
                                        width=2)
            x_pos+=(c_n[k]*np.exp(2*np.pi*1j*k*t/FRAMES)).real
            y_pos+=(c_n[k]*np.exp(2*np.pi*1j*k*t/FRAMES)).imag
            
        else:
            circles[k].center=(x_pos,y_pos)
            arrows[k]=plt.Arrow(x_pos,y_pos,
                                          (c_n[k]*np.exp(-2*np.pi*1j*(len(c_n)-k)*t/FRAMES)).real,
                                          (c_n[k]*np.exp(-2*np.pi*1j*(len(c_n)-k)*t/FRAMES)).imag,
                                          width=2)
            x_pos+=(c_n[k]*np.exp(-2*np.pi*1j*(len(c_n)-k)*t/FRAMES)).real
            y_pos+=(c_n[k]*np.exp(-2*np.pi*1j*(len(c_n)-k)*t/FRAMES)).imag
            
    
    for k in range(0,len(c_n)):
            ax.add_patch(arrows[k])
            ax.add_patch(circles[k])
    
    xdata.append(x_pos) 
    ydata.append(y_pos) 
    line.set_data(xdata, ydata)
        
    return []

#run animation
anim = animation.FuncAnimation(fig, update_animation,
                               init_func=init,
                               frames=FRAMES+10,
                               fargs=(circles,arrows,),
                               interval=15,
                               blit=True,
                               repeat=True)

#save animation as gif
anim.save('seahorse_dynamic_epicycles.gif',writer='pillow')
