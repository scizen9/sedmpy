# -*- coding: utf-8 -*-
"""
Created on Fri Feb 19 17:55:55 2016

@author: nblago
"""
import numpy as np
from matplotlib import pylab as plt
import scipy.optimize as opt


def twoD_Gaussian(xdata_tuple, amplitude, xo, yo, sigma_x, sigma_y, theta, offset):
    '''
    Produces a 2D gaussian centered in xo, yo with the parameters specified.
    xdata_tuple: coordinates of the points where the 2D Gaussian is computed.
    
    '''
    (x, y) = xdata_tuple                                                        
    xo = float(xo)                                                              
    yo = float(yo)                                                              
    a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)   
    b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)    
    c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2)   
    g = offset + amplitude*np.exp( - (a*((x-xo)**2) + 2*b*(x-xo)*(y-yo) + c*((y-yo)**2)))                                   
    return g.ravel()
    

def fit_gauss(img):

    nx, ny = img.shape    
    x = np.linspace(0, nx, nx)
    y = np.linspace(0, ny, ny)
    X, Y = np.meshgrid(x, y)
    
    pix2ang = 0.394
    def_fwhm = 3./pix2ang
    
    def_x = np.argmax(np.sum(img, axis=0))
    def_y = np.argmax(np.sum(img, axis=1))
        
    initial_guess = (10000, def_x, def_y, def_fwhm, def_fwhm, 0, np.percentile(img, 40))
    
    print initial_guess
    
    popt, pcov = opt.curve_fit(twoD_Gaussian, (X, Y), img.flatten(), p0=initial_guess)
    fwhm_x = np.abs(popt[3])*2*np.sqrt(2*np.log(2))
    fwhm_y = np.abs(popt[4])*2*np.sqrt(2*np.log(2))
    amplitude=popt[0]
    background=popt[-1]   
    
    print popt[1], popt[2], nx, ny, def_x, def_y, amplitude, background
    
    return popt[1], popt[2], fwhm_x, fwhm_y, amplitude, background
        
def twoD_Gauss_test(theta=0):
    '''
    Generates a test Gaussian and fits it using the scisoft optimization software.
    '''
    # Create x and y indices
    x = np.linspace(0, 200, 201)
    y = np.linspace(0, 200, 201)
    x, y = np.meshgrid(x, y)
    
    #create data
    data = twoD_Gaussian((x, y), 3, 100, 100, 20, 40, theta, 10)
    

        
    # plot twoD_Gaussian data generated above
    plt.figure()
    plt.imshow(data.reshape(201, 201), origin="bottom", extent=(x.min(), x.max(), y.min(), y.max()))
    plt.colorbar()
    
    # add some noise to the data and try to fit the data generated beforehand
    initial_guess = (3,100,100,20,40,0,10)
    
    data_noisy = data + 0.2*np.random.normal(size=data.shape)
    
    popt, pcov = opt.curve_fit(twoD_Gaussian, (x, y), data_noisy, p0=initial_guess)
    
    data_fitted = twoD_Gaussian((x, y), *popt)
    
    fig, ax = plt.subplots(1, 1)
    ax.hold(True)
    ax.imshow(data_noisy.reshape(201, 201), cmap=plt.cm.jet, origin='bottom',
        extent=(x.min(), x.max(), y.min(), y.max()))
    ax.contour(x, y, data_fitted.reshape(201, 201), 8, colors='w')
    plt.show()