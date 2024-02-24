'''
The module JroAntSetup contains the pre-defined parameters for beam modelling of the Jicamarca ante-
nna. Any new configuration must be added in this module (if the user decides  that) using a specific
ID (pattern value) or it would be read from a file using pattern=None.

MODULES CALLED:
OS, NUMPY
    
MODIFICATION HISTORY:
Created by Ing. Freddy Galindo (frederickgalindo@gmail.com). ROJ Sep 20, 2009.
'''

import os
import numpy

def select_pattern(path=None,filename=None,pattern=0):
    """
        ReturnSetup is a pre-defined list of  Jicamarca antenna configurations which returns a dic-
    tionary giving the configuration  parameters (e.g. transmitted phases). To  choose one, the
    user must define the input "pattern" (See valid values below). 
    
    Parameters:
    -----------
    
    pattern = A integer (>=0) to specify the setup to choose. The default value is zero. If the
      antenna configuration is user-defined pattern must be None.
    
    path = Set this input a string to specifiy  the folder path where the user-defined configu-
      ration file is placed. If this value is not defined ReturnSetup will return None.
      
    file = Set this input a string to specifiy the name of the user-defined  configuration file
      (*.txt). if this value is not defined ReturnSEtup will return None.
    
    Examples:
    ---------
    
    Choosing a pre-defined antenna configuration
    setup = ReturnSetup(pattern=1)
    
    Reading a user-defined antenna configuration
    setup = ReturnSetup(path="/users/users/Progs/Patterns/",file="ExpSep232009.txt")
    """
    
    
    if pattern == 0:
        title = "for module (rx)"
        
        ues = numpy.array([1.,2.,2.,1.])
        phase = numpy.zeros([8,8])
        phase[0:4,:] = 4
        phase[4:8,:] = 5
        
        gaintx = numpy.zeros([8,8])
        gaintx[0,0] = 1

        gainrx = numpy.zeros([8,8])
        gainrx[0,0] = 1

        justrx = 1
    
    elif pattern==1: 
        # Configuration 1/16 on-axis (rx)
        title = " for 1/16 on-axis (rx)"    
        
        ues = numpy.array([1.,2.,2.,1.])
        phase = numpy.zeros([8,8])
        phase[0:4,:] = 4
        phase[4:8,:] = 5
        
        gaintx = numpy.zeros([8,8])
        gaintx[0:2,0:2] = 1
        
        gainrx = numpy.zeros([8,8])
        gainrx[0:2,0:2] = 1
        
        justrx = 1
    
    elif pattern == 2:
        # Configuration for On-Axis
        title = " for 1/4 on-axis (rx)"    
        
        ues = numpy.array([1.,2.,2.,1.])
        phase = numpy.zeros([8,8])
        phase[0:4,:] = 4
        phase[4:8,:] = 5

        gaintx = numpy.zeros([8,8])
        gaintx[0:4,0:4] = 1
        
        gainrx = numpy.zeros([8,8])
        gainrx[0:4,0:4] = 1
        
        justrx = 1
    
    elif pattern == 3:
        # Configuration for On-Axis
        title = " for all on-axis (rx)"    
        
        ues = numpy.array([1.,2.,2.,1.])
        phase = numpy.zeros([8,8])
        phase[0:4,:] = 4
        phase[4:8,:] = 5
        
        gaintx = numpy.zeros([8,8])
        gaintx[:,:] = 1
        
        gainrx = numpy.zeros([8,8])
        gainrx[:,:] = 1
        
        justrx = 0

    elif pattern == 4:
        # Configuration for oblique ISR On-Axis
        title = " for Oblique ISR On-axis"    
        
        ues = numpy.array([1.,2.,2.,1.])
        phase = numpy.zeros([8,8])
        phase[0:4,:] = 4
        phase[4:8,:] = 5
        
        gaintx = numpy.zeros([8,8])
        gaintx[:,:] = 1
        
        gainrx = numpy.zeros([8,8])
        gainrx[:,:] = 1
        
        justrx = 0

    elif pattern == 5:
        # Configuration for oblique ISR "4.5"
        title = " for Oblique ISR '4.5'"    
        
        ues = numpy.array([1.,2.,2.,1.])
        phase = numpy.array([[4,4,5,5,2,2,3,3],
                          [4,5,5,2,2,3,3,4],
                          [5,5,2,2,3,3,4,4],
                          [5,2,2,3,3,4,4,5],
                          [3,3,4,4,5,5,2,2],
                          [3,4,4,5,5,2,2,3],
                          [4,4,5,5,2,2,3,3],
                          [4,5,5,2,2,3,3,4]],dtype=float)
        
        gaintx = numpy.zeros([8,8])
        gaintx[:,:] = 1
        
        gainrx = numpy.zeros([8,8])
        gainrx[:,:] = 1
        
        justrx = 0

    elif pattern == 6:
        # Configuration for oblique ISR "6.0S"
        title = " for Oblique ISR '6.0S'"    

        ues = numpy.array([1.,2.,2.,1.])
        phase = numpy.array([[4,5,2,3,4,5,2,3],
                          [5,2,3,4,5,2,3,4],
                          [2,3,4,5,2,3,4,5],
                          [3,4,5,2,3,4,5,2],
                          [5,2,3,4,5,2,3,4],
                          [2,3,4,5,2,3,4,5],
                          [3,4,5,2,3,4,5,2],
                          [4,5,2,3,4,5,2,3]],dtype=float)
        
        gaintx = numpy.zeros([8,8])
        gaintx[:,:] = 1
        
        gainrx = numpy.zeros([8,8])
        gainrx[:,:] = 1
        
        justrx = 0

    elif pattern == 7:
        # Configuration for oblique ISR "3.0N"
        title = " for Oblique ISR '3.0N'"    
        
        ues = numpy.array([1.,2.,2.,1.])
        phase = numpy.array([[4,3,2,5,4,3,2,5],
                          [3,2,5,4,3,2,5,4],
                          [2,5,4,3,2,5,4,3],
                          [5,4,3,2,5,4,3,2],
                          [5,4,3,2,5,4,3,2],
                          [4,3,2,5,4,3,2,5],
                          [3,2,5,4,3,2,5,4],
                          [2,5,4,3,2,5,4,3]],dtype=float)
        
        gaintx = numpy.zeros([8,8])
        gaintx[:,:] = 1
        
        gainrx = numpy.zeros([8,8])
        gainrx[:,:] = 1
        
        justrx = 0
    
    elif pattern == 8:
        # Configuration for North Fritts"
        title = " for North (Fritts)"    

        ues = numpy.array([2.513, 1.0, 3.0, 0.413])
        phase = numpy.array([[4.29, 3.55, 2.82, 2.08, 4.20, 3.47, 2.73, 2.00],
                          [2.94, 2.20, 5.44, 4.70, 4.32, 3.59, 2.85, 2.12],
                          [5.56, 4.82, 4.09, 3.35, 4.44, 3.71, 2.97, 2.24],
                          [4.20, 3.47, 2.73, 2.00, 4.56, 3.82, 3.09, 2.35],
                          [4.20, 3.47, 2.73, 2.00, 4.56, 3.82, 3.09, 2.35],
                          [4.32, 3.59, 2.85, 2.12, 2.94, 2.20, 5.44, 4.70],
                          [4.44, 3.71, 2.97, 2.24, 5.56, 4.82, 4.09, 3.35],
                          [4.56, 3.82, 3.09, 2.35, 4.20, 3.47, 2.73, 2.00]],dtype=float)
        
        gaintx = numpy.zeros([8,8])
        gaintx[0:4,0:4] = 1
        gaintx[4:8,4:8] = 1
        
        gainrx = numpy.zeros([8,8])
        gainrx[0:4,0:4] = 1
        gainrx[4:8,4:8] = 1
        
        justrx = 0

    elif pattern == 9:
        # Configuration for West Fritts"
        title = " for West (Fritts)"    

        ues = numpy.array([2.513, 1.0, 3.0, 0.413])
        phase = numpy.array([[4.29, 3.55, 2.82, 2.08, 4.20, 3.47, 2.73, 2.00],
                          [2.94, 2.20, 5.44, 4.70, 4.32, 3.59, 2.85, 2.12],
                          [5.56, 4.82, 4.09, 3.35, 4.44, 3.71, 2.97, 2.24],
                          [4.20, 3.47, 2.73, 2.00, 4.56, 3.82, 3.09, 2.35],
                          [4.20, 3.47, 2.73, 2.00, 4.56, 3.82, 3.09, 2.35],
                          [4.32, 3.59, 2.85, 2.12, 2.94, 2.20, 5.44, 4.70],
                          [4.44, 3.71, 2.97, 2.24, 5.56, 4.82, 4.09, 3.35],
                          [4.56, 3.82, 3.09, 2.35, 4.20, 3.47, 2.73, 2.00]],dtype=float)
        
        gaintx = numpy.zeros([8,8])
        gaintx[4:8,0:4] = 1
        gaintx[4:8,0:4] = 1

        gainrx = numpy.zeros([8,8])
        gainrx[4:8,0:4] = 1
        gainrx[4:8,0:4] = 1

        justrx = 0

    elif pattern == 10:
        # Configuration for South Fritts"
        title = " for South (Fritts)"    

        ues = numpy.array([0.413, 2.0, 1.0, 1.513])
        phase = numpy.array([[2.0 , 2.73, 3.47, 4.2 , 2.08, 2.82, 3.55, 4.29],
                          [2.12, 2.85, 3.59, 4.32, 4.7 , 5.44, 2.20, 2.94],
                          [2.24, 2.97, 3.71, 4.44, 3.35, 4.09, 4.82, 5.56],
                          [2.35, 3.09, 3.82, 4.56, 2.0 , 2.73, 3.47, 4.20],
                          [2.08, 2.82, 3.55, 4.29, 2.0 , 2.73, 3.47, 4.20],
                          [4.70, 5.44, 2.20, 2.94, 2.12, 2.85, 3.59, 4.32],
                          [3.35, 4.09, 4.82, 5.56, 2.24, 2.97, 3.71, 4.44],
                          [2.00, 2.73, 3.47, 4.20, 2.35, 3.09, 3.82, 4.56]],dtype=float)

        gaintx = numpy.zeros([8,8])
        gaintx[0:4,0:4] = 1
        gaintx[4:8,4:8] = 1

        gainrx = numpy.zeros([8,8])
        gainrx[0:4,0:4] = 1
        gainrx[4:8,4:8] = 1

        justrx = 0

    elif pattern == 11:
        # Configuration for East Fritts"
        title = " for East (Fritts)"    

        ues = numpy.array([0.413, 2.0, 1.0, 1.513])
        phase = numpy.array([[2.0 , 2.73, 3.47, 4.2 , 2.08, 2.82, 3.55, 4.29],
                          [2.12, 2.85, 3.59, 4.32, 4.7 , 5.44, 2.20, 2.94],
                          [2.24, 2.97, 3.71, 4.44, 3.35, 4.09, 4.82, 5.56],
                          [2.35, 3.09, 3.82, 4.56, 2.0 , 2.73, 3.47, 4.20],
                          [2.08, 2.82, 3.55, 4.29, 2.0 , 2.73, 3.47, 4.20],
                          [4.70, 5.44, 2.20, 2.94, 2.12, 2.85, 3.59, 4.32],
                          [3.35, 4.09, 4.82, 5.56, 2.24, 2.97, 3.71, 4.44],
                          [2.00, 2.73, 3.47, 4.20, 2.35, 3.09, 3.82, 4.56]],dtype=float)
        
        gaintx = numpy.zeros([8,8])
        gaintx[4:8,0:4] = 1
        gaintx[4:8,0:4] = 1

        gainrx = numpy.zeros([8,8])
        gainrx[4:8,0:4] = 1
        gainrx[4:8,0:4] = 1

        justrx = 0

    elif pattern == 12:
        # Configuration for DEWD position (2009)
        title = " for DEWD position (2009) East Beam"    
        
        ues = numpy.array([0.,0.,0.75,0.75])
        phase = numpy.array([[2,3,3,3,3,4,4,4],
                          [5,2,2,2,2,3,3,3],
                          [3,4,4,4,4,5,5,5],
                          [2,3,3,3,3,4,4,4], 
                          [4,5,5,5,5,2,2,2],
                           [3,4,4,4,4,5,5,5],
                            [5,2,2,2,2,3,3,3],
                           [4,5,5,5,5,2,2,2]],dtype=float)
        
        gaintx = numpy.zeros([8,8])
        gaintx[:,:] = 1
        
        gainrx = numpy.zeros([8,8])
        gainrx[:,0:4] = 1
        
        justrx = 0
    
    elif pattern == 13:
        # Configuration for DEWD position (2009)
        title = " for DEWD position (2009) West Beam"    
        
        ues = numpy.array([1.0,0.5,1.5,2.0])
        phase = numpy.array([[5,4,2,5,3,2,4,3],
                          [2,5,3,2,4,3,5,4],
                          [2,5,3,2,4,3,5,4],
                          [3,2,4,3,5,4,2,5],
                          [3,2,4,3,5,4,2,5],
                          [4,3,5,4,2,5,3,2],  
                          [4,3,5,4,2,5,3,2],
                          [5,4,2,5,3,2,4,3]],dtype=float)

        gaintx = numpy.zeros([8,8])
        gaintx[:,:] = 1
        
        gainrx = numpy.zeros([8,8])
        gainrx[4:8,:] = 1
        
        justrx = 0

    elif pattern == 14:
        # Configuration for DVD position (2009)
        title = " for DVD position (2009)"    
        
        ues = numpy.array([1.0,2.0,2.0,1.25])
        phase = numpy.array([[2,2,5,5,4,4,3,3],
                          [2,5,5,4,4,3,3,2],
                          [5,5,4,4,3,3,2,2],
                             [5,4,4,3,3,2,2,5],
                             [5,5,4,4,3,3,2,2],
                          [5,4,4,3,3,2,2,5],  
                          [4,4,3,3,2,2,5,5],
                          [4,3,3,2,2,5,5,4]],dtype=float)

        gaintx = numpy.zeros([8,8])
        gaintx[:,:] = 1
        
        gainrx = numpy.zeros([8,8])
        gainrx[0:4,0:4] = 1

        justrx = 0

    elif pattern == 15:
        # Configuration for Julia CP2
        title = " for Julia CP2 Ew"    
        
        ues = numpy.array([0.0,1.0,1.0,0.0])
        phase = numpy.array([[2,2,5,4,3,3,2,5],
                          [2,5,4,4,3,2,5,5],
                          [5,4,3,3,2,5,4,4],
                          [4,4,3,2,5,5,4,3],
                          [4,4,3,2,5,5,4,3],
                          [4,3,2,2,5,4,3,3],
                          [3,2,5,5,4,3,2,2],
                          [2,2,5,4,3,3,2,5]],dtype=float)

        gaintx = numpy.zeros([8,8])
        gaintx[0:4,4:8] = 1
        gaintx[4:8,0:4] = 1

        gainrx = numpy.zeros([8,8])
        gainrx[0,0] = 1

        justrx = 0
    
    elif pattern == 16:
        # Configuration for Julia CP2
        title = " for Julia CP2 NS"    
        
        ues = numpy.array([1.0,2.0,2.0,1.0])
        phase = numpy.array([[4,4,3,2,5,5,4,3],
                          [4,3,2,2,5,4,3,3],
                          [3,2,5,5,4,3,2,2],
                          [2,2,5,4,3,3,2,5],
                          [2,2,5,4,3,3,2,5],
                          [2,5,4,4,3,2,5,5],
                          [5,4,3,3,2,5,4,4],
                          [4,4,3,2,5,5,4,3]],dtype=float)
        
        gaintx = numpy.zeros([8,8])
        gaintx[0:4,0:4] = 1
        gaintx[4:8,4:8] = 1
        
        gainrx = numpy.zeros([8,8])
        gainrx[0:4,0:4] = 1

        justrx = 0

    elif pattern == 17:
        # Configuration for Julia CP3
        title = " for Julia CP3 NS"    
        
        ues = numpy.array([1.0,1.0,1.0,1.0])
        phase = numpy.array([[4,4,3,2,5,5,4,3],
                          [4,3,2,2,5,4,3,3],
                          [3,2,5,5,4,3,2,2],
                          [2,2,5,4,3,3,2,5],
                          [2,2,5,4,3,3,2,5],
                          [2,5,4,4,3,2,5,5],
                          [5,4,3,3,2,5,4,4],
                          [4,4,3,2,5,5,4,3]],dtype=float)
        
        gaintx = numpy.zeros([8,8])
        gaintx[0:4,0:4] = 1
        gaintx[4:8,4:8] = 1
        
        gainrx = numpy.zeros([8,8])
        gainrx[0:4,0:4] = 1

        justrx = 0

    elif pattern == 18:
        # Configuration for Julia V
        title = " for Julia V"    
        
        ues = (2/3.)*numpy.array([1.5,3.0+0.75,3.0,1.5-0.75])
        phase = numpy.array([[4,4,3,3,2,2,5,5],
                          [4,3,3,2,2,5,5,4],
                          [3,3,2,2,5,5,4,4],
                          [3,2,2,5,5,4,4,3],
                            [3,3,2,2,5,5,4,4],
                             [3,2,2,5,5,4,4,3],
                          [2,2,5,5,4,4,3,3],
                          [2,5,5,4,4,3,3,2]],dtype=float)

        gaintx = numpy.zeros([8,8])
        gaintx[0:4,0:4] = 1
        gaintx[4:8,4:8] = 1
        
        gainrx = numpy.zeros([8,8])
        gainrx[0:4,0:4] = 1

        justrx = 0

    elif pattern == 19:
        # Configuration for Julia V
        title = " for Julia EW 2006-2007 (W)"    
        
        ues = numpy.array([1.0+0.66,2.0+0.66,2.0,1.0])
        phase = numpy.array([[4,3,2,5,4,3,2,5],
                          [4,3,2,5,4,3,2,5],
                          [4,3,2,5,4,3,2,5],
                          [4,3,2,5,4,3,2,5],
                          [5,4,3,2,5,4,3,2],
                          [5,4,3,2,5,4,3,2],
                            [5,4,3,2,5,4,3,2],
                          [5,4,3,2,5,4,3,2]],dtype=float)

        gaintx = numpy.zeros([8,8])
        gaintx[:,:] = 1
        
        gainrx = numpy.zeros([8,8])
        gainrx[:,:] = 1

        justrx = 0

    elif pattern == 20:
        # Configuration for Julia V
        title = " for Julia EW 2006-2007 (E)"    
        
        ues = numpy.array([1.0,1.0,1.0,1.0])
        phase = numpy.array([[4,4,4,4,5,5,5,5],
                          [3,3,3,3,4,4,4,4],
                          [5,5,5,5,2,2,2,2],
                          [4,4,4,4,5,5,5,5],
                          [2,2,2,2,3,3,3,3],
                          [5,5,5,5,2,2,2,2],
                          [3,3,3,3,4,4,4,4],
                          [2,2,2,2,3,3,3,3]],dtype=float)

        gaintx = numpy.zeros([8,8])
        gaintx[:,:] = 1

        gainrx = numpy.zeros([8,8])
        gainrx[0:4,0:4] = 1
        gainrx[4:8,4:8] = 1

        justrx = 0
    
    elif pattern == 21:
        # Configuration for EW Imaging 1996
        title = " for EW Imaging 1996"    
        
        ues = numpy.array([1.0,2.0,2.0,1.0])
        phase = numpy.array([[4,4,3,2,5,5,4,3],
                          [4,3,2,2,5,4,3,3],
                             [3,2,5,5,4,3,2,2],
                          [2,2,5,4,3,3,2,5],
                          [2,2,5,4,3,3,2,5],
                          [2,5,4,4,3,2,5,5],
                          [5,4,3,3,2,5,4,4],
                          [4,4,3,2,5,5,4,3]],dtype=float)

        gaintx = numpy.zeros([8,8])
        gaintx[0:4,0:4] = 1
        gaintx[4:8,4:8] = 1
        
        gainrx = numpy.zeros([8,8])
        gainrx[0,0] = 1

        justrx = 0

    elif pattern == 22:
        # Configuration for EW Imaging 2003
        title = " for EW Imaging 2003"    

        ues = numpy.array([1.0,1.0,1.0,1.0])
        phase = numpy.array([[4,4,3,2,0,0,0,0],
                          [2,3,2,2,0,0,0,0],
                             [5,0,2,5,0,0,0,0],
                          [2,4,3,4,0,0,0,0],
                          [0,0,0,0,3,3,2,5],
                          [0,0,0,0,2,2,5,5],
                          [0,0,0,0,4,3,5,4],
                          [0,0,0,0,5,3,2,3]],dtype=float)

        gaintx = numpy.zeros([8,8])
        gaintx[0:4,0:4] = 1
        gaintx[4:8,4:8] = 1
        
        gainrx = numpy.zeros([8,8])
        gainrx[0,0] = 1

        justrx = 0

    elif pattern == 23:
        # Configuration for EW Imaging 2003
        title = " for EW Imaging 2006-2008"    
        
        ues = numpy.array([1.0,1.0,1.0,2.0])
        phase = numpy.array([[4,4,3,2,0,0,0,0],
                          [2,3,2,2,0,0,0,0],
                             [5,0,2,5,0,0,0,0],
                          [2,4,3,4,0,0,0,0],
                          [0,0,0,0,3,3,2,5],
                          [0,0,0,0,2,2,5,5],
                          [0,0,0,0,4,3,5,4],
                          [0,0,0,0,5,3,2,3]],dtype=float)

        gaintx = numpy.zeros([8,8])
        gaintx[0:4,0:4] = 1
        gaintx[4:8,4:8] = 1
        
        gainrx = numpy.zeros([8,8])
        gainrx[0,0] = 1

        justrx = 0

    elif pattern == 50:
        # Configuration for vertical drift 1996
        title = " for Vertical drift 1996"    
        
        ues = (2/3.)*numpy.array([0.,1.5,1.5,0.])
        phase = numpy.array([[4,4,3,2,5,5,4,3],
                          [4,3,2,2,5,4,3,3],
                          [3,2,5,5,4,3,2,2],
                          [2,2,5,4,3,3,2,5],
                          [2,2,5,4,3,3,2,5],
                          [2,5,4,4,3,2,5,5],
                          [5,4,3,3,2,5,4,4],
                          [4,4,3,2,5,5,4,3]],dtype=float)

        gaintx = numpy.zeros([8,8])
        gaintx[:,:] = 1
        
        gainrx = numpy.zeros([8,8])
        gainrx[:,:] = 1

        justrx = 0
    
    elif pattern == 51:
        # Configuration for vertical drift 1996
        title = " for East-West Drifts 1996 (W beam)"    
        
        ues = numpy.array([0.0,1.0,2.0,1.0])
        phase = numpy.array([[4,3,5,4,2,5,3,2],
                          [4,3,5,4,2,5,3,2],
                          [4,3,5,4,2,5,3,2],
                          [4,3,5,4,2,5,3,2],
                          [5,4,2,5,3,2,4,3],
                          [5,4,2,5,3,2,4,3],
                          [5,4,2,5,3,2,4,3],
                          [5,4,2,5,3,2,4,3]],dtype=float)

        gaintx = numpy.zeros([8,8])
        gaintx[:,:] = 1
        
        gainrx = numpy.zeros([8,8])
        gainrx[4:8,:] = 1

        justrx = 0

    elif pattern == 52:
        # Configuration for vertical drift 1996
        title = " for East-West Drifts 1996 (E Beam)"    
        
        ues = numpy.array([1.0,1.0,0.0,0.0])
        phase = numpy.array([[4,4,4,4,5,5,5,5],
                          [3,3,3,3,4,4,4,4],
                          [5,5,5,5,2,2,2,2],
                          [4,4,4,4,5,5,5,5],
                          [2,2,2,2,3,3,3,3],
                          [5,5,5,5,2,2,2,2],
                          [3,3,3,3,4,4,4,4],
                          [2,2,2,2,3,3,3,3]],dtype=float)

        gaintx = numpy.zeros([8,8])
        gaintx[:,:] = 1
        
        gainrx = numpy.zeros([8,8])
        gainrx[:,0:4] = 1

        justrx = 0

    elif pattern == 53:
        # Configuration for vertical drift 1996
        title = " for DVD position 3 (2006-2008)"    
        
        ues = numpy.array([1.,2,2,1])
        phase = numpy.array([[4,4,3,3,2,2,5,5],
                          [4,3,3,2,2,5,5,4],
                          [3,3,2,2,5,5,4,4],
                          [3,2,2,5,5,4,4,3],
                          [3,3,2,2,5,5,4,4],
                          [3,2,2,5,5,4,4,3],
                          [2,2,5,5,4,4,3,3],
                          [2,5,5,4,4,3,3,2]],dtype=float)

        gaintx = numpy.zeros([8,8])
        gaintx[:,:] = 1
        
        gainrx = numpy.zeros([8,8])
        gainrx[0:4,4:8] = 1

        justrx = 0

    elif pattern == 54:
        # Configuration for vertical drift 1996
        title = " for DEWD (Mar 2005)"    
        
        ues = numpy.array([0.,1.,1/3.,1])
        phase = numpy.array([[4,3,2,5,3,3,3,3],
                          [4,3,2,5,2,2,2,2],
                          [4,3,2,4,5,5,5,5],
                          [4,3,2,4,4,4,3,3],
                          [5,4,3,2,2,2,2,2],
                          [5,4,3,2,5,5,5,5],
                          [5,4,3,5,4,4,4,4],
                          [5,4,3,5,3,3,2,2]],dtype=float)

        gaintx = numpy.zeros([8,8])
        gaintx[:,:] = 1
        
        gainrx = numpy.zeros([8,8])
        gainrx[:,0:4] = 1

        justrx = 0
    
    elif pattern == 55:
        # Configuration for vertical drift 1996
        title = " for DEWD (Mar 2005)"    
        
        ues = numpy.array([0.,1.,1/3.,1])
        phase = numpy.array([[4,3,2,5,3,3,3,3],
                          [4,3,2,5,2,2,2,2],
                          [4,3,2,4,5,5,5,5],
                          [4,3,2,4,4,4,3,3],
                          [5,4,3,2,2,2,2,2],
                          [5,4,3,2,5,5,5,5],
                          [5,4,3,5,4,4,4,4],
                          [5,4,3,5,3,3,2,2]],dtype=float)

        gaintx = numpy.zeros([8,8])
        gaintx[:,:] = 1
        
        gainrx = numpy.zeros([8,8])
        gainrx[0:4:,4:8] = 1

        justrx = 0

    elif pattern ==56:
        # Configuration using antenna compression
        title = " for antenna compression AA*"
        
        ues = numpy.array([0.0,0.0,0.0,0.0])
        phase = numpy.array([[4,4,4,2,4,4,2,4],
                          [4,4,4,2,4,4,2,4],
                          [2,2,2,4,2,2,4,2],
                          [4,4,4,2,4,4,2,4],
                          [2,2,2,4,2,2,4,2],
                          [2,2,2,4,2,2,4,2],
                          [4,4,4,2,4,4,2,4],
                          [2,2,2,4,2,2,4,2]],dtype=float)

        gaintx = numpy.zeros([8,8])
        gaintx[:,:] = 1
        
        gainrx = numpy.zeros([8,8])
        gainrx[0:4,0:4] = 1

        justrx = 0

    elif pattern ==57:
        # Configuration using antenna compression
        title = " for antenna compression AB*"
        
        ues = numpy.array([0.0,0.0,0.0,0.0])
        phase = numpy.array([[4,4,2,4,2,2,4,2],
                          [4,4,2,4,2,2,4,2],
                          [2,2,4,2,4,4,2,4],
                          [4,4,2,4,2,2,4,2],
                          [2,2,4,2,4,4,2,4],
                          [2,2,4,2,4,4,2,4],
                          [4,4,2,4,2,2,4,2],
                          [2,2,4,2,4,4,2,4]],dtype=float)

        gaintx = numpy.zeros([8,8])
        gaintx[:,:] = 1
        
        gainrx = numpy.zeros([8,8])
        gainrx[0:4,0:4] = 1

        justrx = 0

    elif pattern ==58:
        # Configuration using in Oblique ISR 4.5
        title = " for Oblique ISR 4.5"
        
        ues = numpy.array([1.0,2.0,2.0,1.0])
        phase = numpy.array([[4,4,5,5,2,2,3,3],
                          [4,5,5,2,2,3,3,4],
                          [5,5,2,2,3,3,4,4],
                          [5,2,2,3,3,4,4,5],
                          [3,3,4,4,5,5,2,2],
                          [3,4,4,5,5,2,2,3],
                          [4,4,5,5,2,2,3,3],
                          [4,5,5,2,2,3,3,4]],dtype=float)

        gaintx = numpy.zeros([8,8])
        gaintx[:,:] = 1
        
        gainrx = numpy.zeros([8,8])
        gainrx[:,:] = 1

        justrx = 1
    
    elif pattern == 60:
        title=" for Differential phase 2000"
        ues = (2/3.)*numpy.array([0.,1.5-0.5,1.5,0.+0.5])
        
        phase = numpy.array([[4,4,3,2,5,5,4,3],
                          [4,3,2,2,5,4,3,3],
                          [3,2,5,5,4,3,2,2],
                          [2,2,5,4,3,3,2,5],
                          [2,2,5,4,3,3,2,5],
                          [2,5,4,4,3,2,5,5],
                          [5,4,3,3,2,5,4,4],
                          [4,4,3,2,5,5,4,3]],dtype=float)
        
        gaintx = numpy.zeros([8,8])
        gaintx[:,:] = 1
        
        gainrx = numpy.zeros([8,8])
        gainrx[0:4,0:4] = 1
        
        justrx = 0
    
    elif pattern == 61:
        #for East-West 2003 W
        title=" for East-West 2003"
        
        ues = numpy.array([1.+0.66,2.+0.66,2.,1.]) 
        
        phase = numpy.array([[4,3,2,5,4,3,2,5],
                            [4,3,2,5,4,3,2,5],
                            [4,3,2,5,4,3,2,5],
                            [4,3,2,5,4,3,2,5],
                            [5,4,3,2,5,4,3,2],
                            [5,4,3,2,5,4,3,2],
                            [5,4,3,2,5,4,3,2],
                            [5,4,3,2,5,4,3,2]],dtype=float)
        
        gaintx = numpy.zeros([8,8])
        gaintx[:,:] = 1
        
        gainrx = numpy.zeros([8,8])
        gainrx[4:8,:] = 1
        
        justrx = 0
    
    elif pattern == 62:
        #for East-West 2003 E
        title=" for East-West 2003"
        
        ues = numpy.array([1.,1.,0.+1.0,0.+1.0])
        
        phase = numpy.array([[4,4,4,4,5,5,5,5],
                             [3,3,3,3,4,4,4,4],
                             [5,5,5,5,2,2,2,2],
                             [4,4,4,4,5,5,5,5],
                             [2,2,2,2,3,3,3,3],
                             [5,5,5,5,2,2,2,2],
                             [3,3,3,3,4,4,4,4],
                             [2,2,2,2,3,3,3,3]],dtype=float)

        gaintx = numpy.zeros([8,8])
        gaintx[:,:] = 1
        
        gainrx = numpy.zeros([8,8])
        gainrx[:,0:4] = 1
        
        justrx = 0
    
    elif pattern == 63:
        
        title=" for Differential phase 2004 High Alt."
        
        ues = (2/3.)*numpy.array([0.,1.5-1.0,1.5,0.+1.0])
        
        phase = numpy.array([[4,4,3,2,5,5,4,3],
                             [4,3,2,2,5,4,3,3],
                             [3,2,5,5,4,3,2,2],
                             [2,2,5,4,3,3,2,5],
                             [2,2,5,4,3,3,2,5],
                             [2,5,4,4,3,2,5,5],
                             [5,4,3,3,2,5,4,4],
                             [4,4,3,2,5,5,4,3]],dtype=float)

        gaintx = numpy.zeros([8,8])
        gaintx[:,:] = 1
        
        gainrx = numpy.zeros([8,8])
        gainrx[:,:] = 1
        
        justrx = 0
    
    elif pattern == 64:
        
        title=" for Differential Phase Perp to B 2005-2006"
        
        ues = (2/3.)*numpy.array([1.5,3.0+0.75,3.0,1.5-0.75])
        
        phase = numpy.array([[4,4,3,3,2,2,5,5],
                             [4,3,3,2,2,5,5,4],
                             [3,3,2,2,5,5,4,4],
                             [3,2,2,5,5,4,4,3],
                             [3,3,2,2,5,5,4,4],
                             [3,2,2,5,5,4,4,3],
                             [2,2,5,5,4,4,3,3],
                             [2,5,5,4,4,3,3,2]],dtype=float)

        gaintx = numpy.zeros([8,8])
        gaintx[:,:] = 1
        
        gainrx = numpy.zeros([8,8])
        gainrx[0:4,4:8] = 1
        
        justrx = 0
    
    elif pattern == 65:
        #for JULIA EW 2003 W
        title=" for JULIA EW 2003"
        
        ues = numpy.array([1+0.66,2+0.66,2.,1.])
        
        phase = numpy.array([[4,3,2,5,4,3,2,5],
                             [4,3,2,5,4,3,2,5],
                             [4,3,2,5,4,3,2,5],
                             [4,3,2,5,4,3,2,5],
                             [5,4,3,2,5,4,3,2],
                             [5,4,3,2,5,4,3,2],
                             [5,4,3,2,5,4,3,2],
                             [5,4,3,2,5,4,3,2]],dtype=float)

        gaintx = numpy.zeros([8,8])
        gaintx[:,:] = 1
        
        gainrx = numpy.zeros([8,8])
        gainrx[4:8,:] = 1
        
        justrx = 0
    
    elif pattern == 66:
        #for JULIA EW 2003 E
        title=" for JULIA EW 2003"
        
        ues = numpy.array([1.,1.,0.,0.])
        
        phase = numpy.array([[4,4,4,4,5,5,5,5],
                             [3,3,3,3,4,4,4,4],
                             [5,5,5,5,2,2,2,2],
                             [4,4,4,4,5,5,5,5],
                             [2,2,2,2,3,3,3,3],
                             [5,5,5,5,2,2,2,2],
                             [3,3,3,3,4,4,4,4],
                             [2,2,2,2,3,3,3,3]],dtype=float)

        gaintx = numpy.zeros([8,8])
        gaintx[:,:] = 1
        
        gainrx = numpy.zeros([8,8])
        gainrx[:,0:4] = 1
        
        justrx = 0
    
    elif pattern == 67:
        
        title=" for Vertical (Yellow Cables)"
        
        ues = numpy.array([0.25, 0.25, 0.25, 0.25])
        
        phase = numpy.array([[3.41, 3.41, 3.41, 3.41, 3.41, 3.41, 3.41, 3.41],
                             [2.78, 2.78, 2.78, 2.78, 2.78, 2.78, 2.78, 2.78],
                             [2.15, 2.15, 2.15, 2.15, 2.15, 2.15, 2.15, 2.15],
                             [5.52, 5.52, 5.52, 5.52, 5.52, 5.52, 5.52, 5.52],
                             [4.89, 4.89, 4.89, 4.89, 4.89, 4.89, 4.89, 4.89],
                             [4.26, 4.26, 4.26, 4.26, 4.26, 4.26, 4.26, 4.26],
                             [3.63, 3.63, 3.63, 3.63, 3.63, 3.63, 3.63, 3.63],
                             [3.00, 3.00, 3.00, 3.00, 3.00, 3.00, 3.00, 3.00]],dtype=float)

        gaintx = numpy.zeros([8,8])
        gaintx[:,:] = 1
        
        gainrx = numpy.zeros([8,8])
        gainrx[:,:] = 1
        
        justrx = 0
    elif pattern == 68:
        # MST_ISR_EEJ
        title=" RX_NS_UP:CH3"
        
        ues = numpy.array([0.8,0.0,1.6,0.0])*2/3
        
        phase = numpy.array([[4,5,2,3,2,3,3,3],
                             [5,2,3,4,5,2,2,2],
                             [2,3,4,5,3,4,4,4],
                             [3,4,5,2,2,3,3,3],
                             [4,5,5,5,4,5,2,3],
                             [3,4,4,4,5,2,3,4],
                             [5,2,2,2,2,3,4,5],
                             [4,5,5,5,3,4,5,2]],dtype=float)

        gaintx = numpy.zeros([8,8])
        gaintx[:,:] = 1
        
        gainrx = numpy.zeros([8,8])
        gainrx[0:4,0:4] = 1
        gainrx[4:8,4:8] = 1
        
        justrx = 0

    elif pattern == 69:
        # MST_ISR_EEJ
        title=" RX_WE_UP:CH1"
        
        ues = numpy.array([0.8,0.0,1.6,0.0])*2/3
        
        phase = numpy.array([[4,5,2,3,2,3,3,3],
                             [5,2,3,4,5,2,2,2],
                             [2,3,4,5,3,4,4,4],
                             [3,4,5,2,2,3,3,3],
                             [4,5,5,5,4,5,2,3],
                             [3,4,4,4,5,2,3,4],
                             [5,2,2,2,2,3,4,5],
                             [4,5,5,5,3,4,5,2]],dtype=float)

        gaintx = numpy.zeros([8,8])
        gaintx[:,:] = 1
        
        gainrx = numpy.zeros([8,8])
        gainrx[4:8,0:4] = 1
        gainrx[0:4,4:8] = 1
        
        justrx = 0

    elif pattern == 70:
        # MST_ISR_EEJ
        title=" RX_NS_DW:CH4"
        
        ues = numpy.array([3.4,1.5,0.8,0.25])*2/3
        
        phase = numpy.array([[3.41,3.41,3.41,3.41,2,5,3,2],
                             [2.78,2.78,2.78,2.78,3,2,4,3],
                             [2.15,2.15,2.15,2.15,3,2,4,3],
                             [5.52,5.52,5.52,5.52,4,3,5,4],
                             [4,3,5,4,4.89,4.89,4.89,4.89],
                             [5,4,2,5,4.26,4.26,4.26,4.26],
                             [5,4,2,5,3.63,3.63,3.63,3.63],
                             [2,5,3,2,3,3,3,3]],dtype=float)

        gaintx = numpy.zeros([8,8])
        gaintx[:,:] = 1
        
        gainrx = numpy.zeros([8,8])
        gainrx[0:4,0:4] = 1
        gainrx[4:8,4:8] = 1
        
        justrx = 0

    elif pattern == 71:
        # MST_ISR_EEJ
        title=" RX_WE_DW:CH2"
        
        ues = numpy.array([3.4,1.5,0.8,0.25])*2/3
        
        phase = numpy.array([[3.41,3.41,3.41,3.41,2,5,3,2],
                             [2.78,2.78,2.78,2.78,3,2,4,3],
                             [2.15,2.15,2.15,2.15,3,2,4,3],
                             [5.52,5.52,5.52,5.52,4,3,5,4],
                             [4,3,5,4,4.89,4.89,4.89,4.89],
                             [5,4,2,5,4.26,4.26,4.26,4.26],
                             [5,4,2,5,3.63,3.63,3.63,3.63],
                             [2,5,3,2,3,3,3,3]],dtype=float)

        gaintx = numpy.zeros([8,8])
        gaintx[:,:] = 1
        
        gainrx = numpy.zeros([8,8])
        gainrx[4:8,0:4] = 1
        gainrx[0:4,4:8] = 1
        
        justrx = 0
                    
    elif pattern == 72:
        # MST_ISR_EEJ
        title=" RX_MST_ISR_EEJ_WW"
        
        ues = numpy.array([3.4,1.5,0.8,0.25])*2/3
        
        phase = numpy.array([[3.41,3.41,3.41,3.41,2,5,3,2],
                             [2.78,2.78,2.78,2.78,3,2,4,3],
                             [2.15,2.15,2.15,2.15,3,2,4,3],
                             [5.52,5.52,5.52,5.52,4,3,5,4],
                             [4,3,5,4,4.89,4.89,4.89,4.89],
                             [5,4,2,5,4.26,4.26,4.26,4.26],
                             [5,4,2,5,3.63,3.63,3.63,3.63],
                             [2,5,3,2,3,3,3,3]],dtype=float)

        gaintx = numpy.zeros([8,8])
        gaintx[:,:] = 1
        
        gainrx = numpy.zeros([8,8])
        gainrx[4:8,0:4] = 1
        gainrx[0:4,4:8] = 1
        
        justrx = 0
    
    elif pattern == 73:
        # EW_DRIFT
        title=" chA : wUP"
        
        ues = numpy.array([5,0,3,4.5])*2/3
        
        phase = numpy.array([[4,4,5,5,2,5,3,2],
                             [4,5,5,2,3,2,4,3],
                             [5,5,2,2,4,3,5,4],
                             [5,2,2,3,5,4,2,5],
                             [2,5,3,2,5,5,2,2],
                             [3,2,4,3,5,2,2,3],
                             [4,3,5,4,2,2,3,3],
                             [5,4,2,5,2,3,3,4]],dtype=float)

        gaintx = numpy.zeros([8,8])
        gaintx[0:4,4:8] = 1
        gaintx[4:8,0:4] = 1
        
        gainrx = numpy.zeros([8,8])
        gainrx[4:8,0:4] = 1
        
        justrx = 0
    
    elif pattern == 74:
        # EW_DRIFT
        title=" chB : eUP"
        
        ues = numpy.array([5,0,3,4.5])*2/3
        
        phase = numpy.array([[4,4,5,5,2,5,3,2],
                             [4,5,5,2,3,2,4,3],
                             [5,5,2,2,4,3,5,4],
                             [5,2,2,3,5,4,2,5],
                             [2,5,3,2,5,5,2,2],
                             [3,2,4,3,5,2,2,3],
                             [4,3,5,4,2,2,3,3],
                             [5,4,2,5,2,3,3,4]],dtype=float)

        gaintx = numpy.zeros([8,8])
        gaintx[0:4,4:8] = 1
        gaintx[4:8,0:4] = 1
        
        gainrx = numpy.zeros([8,8])
        gainrx[0:4,4:8] = 1
        
        justrx = 0

    elif pattern == 75:
        # EW_DRIFT
        title=" chC : wDW"
        
        ues = numpy.array([3,4.5,4.5,3])*2/3
        
        phase = numpy.array([[2,2,3,3,2,2,3,3],
                             [2,3,3,4,5,5,2,2],
                             [3,3,4,4,3,3,4,4],
                             [3,4,4,5,2,2,3,3],
                             [4,4,5,5,3,3,4,4],
                             [3,3,4,4,3,4,4,5],
                             [5,5,2,2,4,4,5,5],
                             [4,4,5,5,4,5,5,2]],dtype=float)

        gaintx = numpy.zeros([8,8])
        gaintx[0:4,4:8] = 1
        gaintx[4:8,0:4] = 1
        
        gainrx = numpy.zeros([8,8])
        gainrx[4:8,0:4] = 1
        
        justrx = 0
    
    elif pattern == 76:
        # EW_DRIFT
        title=" chD : eDW"
        
        ues = numpy.array([3,4.5,4.5,3])*2/3
        
        phase = numpy.array([[2,2,3,3,2,2,3,3],
                             [2,3,3,4,5,5,2,2],
                             [3,3,4,4,3,3,4,4],
                             [3,4,4,5,2,2,3,3],
                             [4,4,5,5,3,3,4,4],
                             [3,3,4,4,3,4,4,5],
                             [5,5,2,2,4,4,5,5],
                             [4,4,5,5,4,5,5,2]],dtype=float)

        gaintx = numpy.zeros([8,8])
        gaintx[0:4,4:8] = 1
        gaintx[4:8,0:4] = 1
        
        gainrx = numpy.zeros([8,8])
        gainrx[0:4,4:8] = 1
        
        justrx = 0

    elif pattern == 77:
        # FARADAY
        title=" chE : (4.5) nsUP"
        
        ues = numpy.array([5,0,3,4.5])*2/3
        
        phase = numpy.array([[4,4,5,5,2,5,3,2],
                             [4,5,5,2,3,2,4,3],
                             [5,5,2,2,4,3,5,4],
                             [5,2,2,3,5,4,2,5],
                             [2,5,3,2,5,5,2,2],
                             [3,2,4,3,5,2,2,3],
                             [4,3,5,4,2,2,3,3],
                             [5,4,2,5,2,3,3,4]],dtype=float)

        gaintx = numpy.zeros([8,8])
        gaintx[0:4,0:4] = 1
        gaintx[4:8,4:8] = 1
        
        gainrx = numpy.zeros([8,8])
        gainrx[0:4,0:4] = 1
        gainrx[4:8,4:8] = 1
        
        justrx = 0

    elif pattern == 78:
        # FARADAY
        title=" chF : (4.5) nsDW"
        
        ues = numpy.array([3,4.5,4.5,3])*2/3
        
        phase = numpy.array([[2,2,3,3,2,2,3,3],
                             [2,3,3,4,5,5,2,2],
                             [3,3,4,4,3,3,4,4],
                             [3,4,4,5,2,2,3,3],
                             [4,4,5,5,3,3,4,4],
                             [3,3,4,4,3,4,4,5],
                             [5,5,2,2,4,4,5,5],
                             [4,4,5,5,4,5,5,2]],dtype=float)

        gaintx = numpy.zeros([8,8])
        gaintx[0:4,0:4] = 1
        gaintx[4:8,4:8] = 1
        
        gainrx = numpy.zeros([8,8])
        gainrx[0:4,0:4] = 1
        gainrx[4:8,4:8] = 1
        
        justrx = 0
    
    elif pattern == 79:
        # ON AXIS
        title=" wUP "
        
#        ues = numpy.array([0,0,0,0])*2/3
        ues = numpy.array([0,1.5,1.5,0])*2/3
        
        phase = numpy.array([[4,4,4,4,4,4,4,4],
                             [4,4,4,4,4,4,4,4],
                             [4,4,4,4,4,4,4,4],
                             [4,4,4,4,4,4,4,4],
                             [5,5,5,5,5,5,5,5],
                             [5,5,5,5,5,5,5,5],
                             [5,5,5,5,5,5,5,5],
                             [5,5,5,5,5,5,5,5]],dtype=float)

        gaintx = numpy.zeros([8,8])
        #gaintx[0:4,0:4] = 1
        #gaintx[4:8,4:8] = 1
        
        gainrx = numpy.zeros([8,8])
        gainrx[4:8,0:4] = 1
        
        justrx = 1               

    elif pattern == 80:
        # ON AXIS
        title=" CHD : nUP (ON_AXIS)"
        
        ues = numpy.array([0,0,0,0])*2/3
        
        phase = numpy.array([[4,4,4,4,4,4,4,4],
                             [4,4,4,4,4,4,4,4],
                             [4,4,4,4,4,4,4,4],
                             [4,4,4,4,4,4,4,4],
                             [5,5,5,5,5,5,5,5],
                             [5,5,5,5,5,5,5,5],
                             [5,5,5,5,5,5,5,5],
                             [5,5,5,5,5,5,5,5]],dtype=float)

        gaintx = numpy.zeros([8,8])
        gaintx[0:4,0:4] = 1
        gaintx[4:8,4:8] = 1
        
        gainrx = numpy.zeros([8,8])
        gainrx[0:4,0:4] = 1
        
        justrx = 1
                        
    elif pattern == 81:
        # ON AXIS
        title=" CHE : eUP (ON_AXIS)"
        
        ues = numpy.array([0,0,0,0])*2/3
        
        phase = numpy.array([[4,4,4,4,4,4,4,4],
                             [4,4,4,4,4,4,4,4],
                             [4,4,4,4,4,4,4,4],
                             [4,4,4,4,4,4,4,4],
                             [5,5,5,5,5,5,5,5],
                             [5,5,5,5,5,5,5,5],
                             [5,5,5,5,5,5,5,5],
                             [5,5,5,5,5,5,5,5]],dtype=float)

        gaintx = numpy.zeros([8,8])
        gaintx[0:4,0:4] = 1
        gaintx[4:8,4:8] = 1
        
        gainrx = numpy.zeros([8,8])
        gainrx[0:4,4:8] = 1
        
        justrx = 1

    elif pattern == 82:
        # ON AXIS
        title=" CHH : sUP (ON_AXIS)"
        
        ues = numpy.array([0,0,0,0])*2/3
        
        phase = numpy.array([[4,4,4,4,4,4,4,4],
                             [4,4,4,4,4,4,4,4],
                             [4,4,4,4,4,4,4,4],
                             [4,4,4,4,4,4,4,4],
                             [5,5,5,5,5,5,5,5],
                             [5,5,5,5,5,5,5,5],
                             [5,5,5,5,5,5,5,5],
                             [5,5,5,5,5,5,5,5]],dtype=float)

        gaintx = numpy.zeros([8,8])
        gaintx[0:4,0:4] = 1
        gaintx[4:8,4:8] = 1
        
        gainrx = numpy.zeros([8,8])
        gainrx[4:8,4:8] = 1
        
        justrx = 1

    elif pattern == 83:
        # ON AXIS
        title=" wDW "
        
        ues = numpy.array([1.5,3.0,3.0,1.5])*2/3
        
        phase = numpy.array([[2,2,2,2,2,2,2,2],
                             [2,2,2,2,2,2,2,2],
                             [2,2,2,2,2,2,2,2],
                             [2,2,2,2,2,2,2,2],
                             [3,3,3,3,3,3,3,3],
                             [3,3,3,3,3,3,3,3],
                             [3,3,3,3,3,3,3,3],
                             [3,3,3,3,3,3,3,3]],dtype=float)

        gaintx = numpy.zeros([8,8])
        gaintx[0:4,0:4] = 1
        gaintx[4:8,4:8] = 1
        
        gainrx = numpy.zeros([8,8])
        gainrx[4:8,0:4] = 1
        
        justrx = 1

    elif pattern == 84:
        # ON AXIS
        title=" CHC : nDW (ON_AXIS)"
        
        ues = numpy.array([0,0,0,0])*2/3
        
        phase = numpy.array([[2,2,2,2,2,2,2,2],
                             [2,2,2,2,2,2,2,2],
                             [2,2,2,2,2,2,2,2],
                             [2,2,2,2,2,2,2,2],
                             [3,3,3,3,3,3,3,3],
                             [3,3,3,3,3,3,3,3],
                             [3,3,3,3,3,3,3,3],
                             [3,3,3,3,3,3,3,3]],dtype=float)

        gaintx = numpy.zeros([8,8])
        gaintx[0:4,0:4] = 1
        gaintx[4:8,4:8] = 1
        
        gainrx = numpy.zeros([8,8])
        gainrx[0:4,0:4] = 1
        
        justrx = 1

    elif pattern == 85:
        # ON AXIS
        title=" CHF : eDW (ON_AXIS)"
        
        ues = numpy.array([0,0,0,0])*2/3
        
        phase = numpy.array([[2,2,2,2,2,2,2,2],
                             [2,2,2,2,2,2,2,2],
                             [2,2,2,2,2,2,2,2],
                             [2,2,2,2,2,2,2,2],
                             [3,3,3,3,3,3,3,3],
                             [3,3,3,3,3,3,3,3],
                             [3,3,3,3,3,3,3,3],
                             [3,3,3,3,3,3,3,3]],dtype=float)

        gaintx = numpy.zeros([8,8])
        gaintx[0:4,0:4] = 1
        gaintx[4:8,4:8] = 1
        
        gainrx = numpy.zeros([8,8])
        gainrx[0:4,4:8] = 1
        
        justrx = 1

    elif pattern == 86:
        # ON AXIS
        title=" CHG : sDW (ON_AXIS)"
        
        ues = numpy.array([0,0,0,0])*2/3
        
        phase = numpy.array([[2,2,2,2,2,2,2,2],
                             [2,2,2,2,2,2,2,2],
                             [2,2,2,2,2,2,2,2],
                             [2,2,2,2,2,2,2,2],
                             [3,3,3,3,3,3,3,3],
                             [3,3,3,3,3,3,3,3],
                             [3,3,3,3,3,3,3,3],
                             [3,3,3,3,3,3,3,3]],dtype=float)

        gaintx = numpy.zeros([8,8])
        gaintx[0:4,0:4] = 1
        gaintx[4:8,4:8] = 1
        
        gainrx = numpy.zeros([8,8])
        gainrx[4:8,4:8] = 1
        
        justrx = 1
                         
    elif pattern == 87:
        # JULIA 10/2009
        title=" JULIA D: n_s_UP"
        
        ues = numpy.array([0,0,0,0.25])*2/3
        ues = numpy.array([0,0,0,0.25])
        ues = numpy.array([0,0,0,6*0.25])
        
        phase = numpy.array([[4,4,4,4,4,3,4,4],
                             [3,3,3,3,3,3,4,4],
                             [2,2,2,2,4,3,3,2],
                             [5,5,5,5,5,4,2,2],
                             [4,4,4,4,4,4,4,4],
                             [4,4,4,4,3,3,3,3],
                             [4,4,4,4,2,2,2,2],
                             [4,4,4,4,5,5,5,5]],dtype=float)

        gaintx = numpy.zeros([8,8])
        gaintx[0:4,0:4] = 1
        gaintx[4:8,4:8] = 1
        gaintx[4,6] = 0
        
        gainrx = numpy.zeros([8,8])
        gainrx[0:4,0:4] = 1
        gainrx[4:8,4:8] = 1
        gainrx[4,6] = 0
        
        justrx = 0

    elif pattern == 88:
        # JULIA 10/2009
        title=" JULIA CHC : n_s_DW"
        
        ues = numpy.array([0,0,0,0.75])*2/3
        ues = numpy.array([0,0,0,3*0.75])
        
        phase = numpy.array([[5,4,2,5,3,2,4,3],
                             [2,5,3,2,4,3,5,4],
                             [2,5,3,2,4,3,5,4],
                             [3,2,4,3,5,4,2,5],
                             [3,2,4,3,5,4,2,5],
                             [4,3,5,4,2,5,3,2],
                             [4,3,5,4,3,5,3,2],
                             [5,4,2,5,3,2,4,3]],dtype=float)

        gaintx = numpy.zeros([8,8])
        gaintx[0:4,0:4] = 1
        gaintx[4:8,4:8] = 1
        
        gainrx = numpy.zeros([8,8])
        gainrx[0:4,0:4] = 1
        gainrx[4:8,4:8] = 1
        justrx = 0

    elif pattern == 89:
        # JULIA
        title=" JULIA CHA : wDW"
        
        ues = numpy.array([0,0,0,0.75])*2/3
        ues = numpy.array([0,0,0,3*0.75])
        
        phase = numpy.array([[5,4,3,2,5,4,3,2],
                             [2,5,4,3,2,5,4,3],
                             [2,5,4,3,2,5,4,3],
                             [3,2,5,4,3,2,5,4],
                             [5,4,3,2,3,2,5,4],
                             [2,5,4,3,4,3,2,5],
                             [2,5,4,3,4,3,2,5],
                             [3,2,5,4,5,4,3,2]],dtype=float)

        gaintx = numpy.zeros([8,8])
        gaintx[0:4,0:4] = 1
        gaintx[4:8,4:8] = 1
        
        gainrx = numpy.zeros([8,8])
        gainrx[4:8,0:4] = 1
        
        justrx = 0
        
    elif pattern == 90:
        # JULIA
        title=" JULIA CHB : eDW"
        
        ues = numpy.array([0,0,0,0.75])*2/3
        ues = numpy.array([0,0,0,0.75])
        
        phase = numpy.array([[5,4,3,2,5,4,3,2],
                             [2,5,4,3,2,5,4,3],
                             [2,5,4,3,2,5,4,3],
                             [3,2,5,4,3,2,5,4],
                             [5,4,3,2,3,2,5,4],
                             [2,5,4,3,4,3,2,5],
                             [2,5,4,3,4,3,2,5],
                             [3,2,5,4,5,4,3,2]],dtype=float)

        gaintx = numpy.zeros([8,8])
        gaintx[0:4,0:4] = 1
        gaintx[4:8,4:8] = 1
        
        gainrx = numpy.zeros([8,8])
        gainrx[0:4,4:8] = 1
        
        justrx = 0
        
    elif pattern == 91:
        # MST_ISR_EEJ (Mayo 2014)
        title=" RX_NS_UP:CH3"
        
        ues = numpy.array([5.0,0.0,3.0,0.0])*2/3
        
        phase = numpy.array([[4,5,2,3,2,5,3,2],
                             [5,2,3,4,3,2,4,3],
                             [2,3,4,5,4,3,5,4],
                             [3,4,5,2,5,4,2,5],
                             [2,5,3,2,4,5,2,3],
                             [3,2,4,3,5,2,3,4],
                             [4,3,5,4,2,3,4,5],
                             [5,4,2,5,3,4,5,2]],dtype=float)

        gaintx = numpy.zeros([8,8])
        gaintx[:,:] = 1
        
        gainrx = numpy.zeros([8,8])
        gainrx[0:4,0:4] = 1
        gainrx[4:8,4:8] = 1
        
        justrx = 0

    elif pattern == 92:
        # MST_ISR_EEJ (Mayo 2014)
        title=" RX_WE_UP:CH2"
        
        ues = numpy.array([5.0,0.0,3.0,0.0])*2/3
        
        phase = numpy.array([[4,5,2,3,2,5,3,2],
                             [5,2,3,4,3,2,4,3],
                             [2,3,4,5,4,3,5,4],
                             [3,4,5,2,5,4,2,5],
                             [2,5,3,2,4,5,2,3],
                             [3,2,4,3,5,2,3,4],
                             [4,3,5,4,2,3,4,5],
                             [5,4,2,5,3,4,5,2]],dtype=float)

        gaintx = numpy.zeros([8,8])
        gaintx[:,:] = 1
        
        gainrx = numpy.zeros([8,8])
        gainrx[4:8,0:4] = 1
        gainrx[0:4,4:8] = 1
        
        justrx = 0

    elif pattern == 93:
        # MST_ISR_EEJ (Mayo 2014)
        title=" RX_NS_DW:CH4"
        
        ues = numpy.array([3.0,1.5,4.5,0.25])*2/3
        
        phase = numpy.array([[3.41,3.41,3.41,3.41,2,2,3,3],
                             [2.78,2.78,2.78,2.78,5,5,2,2],
                             [2.15,2.15,2.15,2.15,3,3,4,4],
                             [5.52,5.52,5.52,5.52,2,2,3,3],
                             [4,4,5,5,4.89,4.89,4.89,4.89],
                             [3,3,4,4,4.26,4.26,4.26,4.26],
                             [5,5,2,2,3.63,3.63,3.63,3.63],
                             [4,4,5,5,3,3,3,3]],dtype=float)

        gaintx = numpy.zeros([8,8])
        gaintx[:,:] = 1
        
        gainrx = numpy.zeros([8,8])
        gainrx[0:4,0:4] = 1
        gainrx[4:8,4:8] = 1
        
        justrx = 0

    elif pattern == 94:
        # MST_ISR_EEJ (Mayo 2014)
        title=" RX_WE_DW:CH1"
        
        ues = numpy.array([3.0,1.5,4.5,0.25])*2/3
        
        phase = numpy.array([[3.41,3.41,3.41,3.41,2,2,3,3],
                             [2.78,2.78,2.78,2.78,5,5,2,2],
                             [2.15,2.15,2.15,2.15,3,3,4,4],
                             [5.52,5.52,5.52,5.52,2,2,3,3],
                             [4,4,5,5,4.89,4.89,4.89,4.89],
                             [3,3,4,4,4.26,4.26,4.26,4.26],
                             [5,5,2,2,3.63,3.63,3.63,3.63],
                             [4,4,5,5,3,3,3,3]],dtype=float)

        gaintx = numpy.zeros([8,8])
        gaintx[:,:] = 1
        
        gainrx = numpy.zeros([8,8])
        gainrx[4:8,0:4] = 1
        gainrx[0:4,4:8] = 1
        
        justrx = 0

    elif pattern == 95:
        # MST_ISR_EEJ (Setiembre 2014)
        title=" RX_NS_UP:CH3"
        
        ues = numpy.array([5.0,0.0,3.0,0.0])*2/3
        
        phase = numpy.array([[4,5,2,3,2,5,3,2],
                             [5,2,3,4,3,2,4,3],
                             [2,3,4,5,4,3,5,4],
                             [3,4,5,2,5,4,2,5],
                             [2,5,3,2,4,5,2,3],
                             [3,2,4,3,5,2,3,4],
                             [4,3,5,4,2,3,4,5],
                             [5,4,2,5,3,4,5,2]],dtype=float)

        gaintx = numpy.zeros([8,8])
        gaintx[:,:] = 1
        
        gainrx = numpy.zeros([8,8])
        gainrx[0:4,0:4] = 1
        gainrx[4:8,4:8] = 1
        
        justrx = 0

    elif pattern == 96:
        # MST_ISR_EEJ (Setiembre 2014)
        title=" RX_WE_UP:CH2"
        
        ues = numpy.array([5.0,0.0,3.0,0.0])*2/3
        
        phase = numpy.array([[4,5,2,3,2,5,3,2],
                             [5,2,3,4,3,2,4,3],
                             [2,3,4,5,4,3,5,4],
                             [3,4,5,2,5,4,2,5],
                             [2,5,3,2,4,5,2,3],
                             [3,2,4,3,5,2,3,4],
                             [4,3,5,4,2,3,4,5],
                             [5,4,2,5,3,4,5,2]],dtype=float)

        gaintx = numpy.zeros([8,8])
        gaintx[:,:] = 1
        
        gainrx = numpy.zeros([8,8])
        gainrx[4:8,0:4] = 1
        gainrx[0:4,4:8] = 1
        
        justrx = 0

    elif pattern == 97:
        # MST_ISR_EEJ (Setiembre 2014)
        title=" RX_NS_DW:CH4"
        
        #ues = numpy.array([3.0,1.5,4.5,0.25])*2/3
        ues = numpy.array([3.0,0.0,4.5,0.0])*2/3
        
        phase = numpy.array([[3.41,3.41,3.41,3.41,2,2,3,3],
                             [2.78,2.78,2.78,2.78,5,5,2,2],
                             [2.15,2.15,2.15,2.15,3,3,4,4],
                             [5.52,5.52,5.52,5.52,2,2,3,3],
                             [4,4,5,5,4.89,4.89,4.89,4.89],
                             [3,3,4,4,4.26,4.26,4.26,4.26],
                             [5,5,2,2,3.63,3.63,3.63,3.63],
                             [4,4,5,5,3,3,3,3]],dtype=float)

        gaintx = numpy.zeros([8,8])
        gaintx[:,:] = 1
        
        gainrx = numpy.zeros([8,8])
        gainrx[0:4,0:4] = 1
        gainrx[4:8,4:8] = 1
        
        justrx = 0

    elif pattern == 98:
        # MST_ISR_EEJ (Setiembre 2014)
        title=" RX_WE_DW:CH1"
        
        ues = numpy.array([3.0,0.0,4.5,0.0])*2/3
        
        phase = numpy.array([[3.41,3.41,3.41,3.41,2,2,3,3],
                             [2.78,2.78,2.78,2.78,5,5,2,2],
                             [2.15,2.15,2.15,2.15,3,3,4,4],
                             [5.52,5.52,5.52,5.52,2,2,3,3],
                             [4,4,5,5,4.89,4.89,4.89,4.89],
                             [3,3,4,4,4.26,4.26,4.26,4.26],
                             [5,5,2,2,3.63,3.63,3.63,3.63],
                             [4,4,5,5,3,3,3,3]],dtype=float)

        gaintx = numpy.zeros([8,8])
        gaintx[:,:] = 1
        
        gainrx = numpy.zeros([8,8])
        gainrx[4:8,0:4] = 1
        gainrx[0:4,4:8] = 1
        
        justrx = 0
                        
#     elif pattern == 99:
#         # MST_ISR_EEJ (Setiembre 2014)
#         title=" RX_WE_DW:CH1"
#         
#         ues = numpy.array([0.0,0.0,0.0,0.0])*2/3
#         
#         phase = numpy.array([[3.41,3.41,3.41,3.41,2,2,3,3],
#                              [2.78,2.78,2.78,2.78,5,5,2,2],
#                              [2.15,2.15,2.15,2.15,3,3,4,4],
#                              [5.52,5.52,5.52,5.52,2,2,3,3],
#                              [4,4,5,5,4.89,4.89,4.89,4.89],
#                              [3,3,4,4,4.26,4.26,4.26,4.26],
#                              [5,5,2,2,3.63,3.63,3.63,3.63],
#                              [4,4,5,5,3,3,3,3]],dtype=float)
# 
#         gaintx = numpy.zeros([8,8])
#         gaintx[:,:] = 1
#         
#         gainrx = numpy.zeros([8,8])
#         gainrx[1:2,6:7] = 1
#         
#         justrx = 1
        
    elif pattern == 101:
        # ON AXIS
        title=" ON_AXIS_wUP"
        
        ues = numpy.array([1,2,2,0])
        
        phase = numpy.array([[4,4,4,4,4,4,4,4],
                             [4,4,4,4,4,4,4,4],
                             [4,4,4,4,4,4,4,4],
                             [4,4,4,4,4,4,4,4],
                             [5,5,5,5,5,5,5,5],
                             [5,5,5,5,5,5,5,5],
                             [5,5,5,5,5,5,5,5],
                             [5,5,5,5,5,5,5,5]],dtype=float)

        gaintx = numpy.zeros([8,8])
        gaintx[0:4,4:8] = 1
        gaintx[4:8,0:4] = 1
        
        gainrx = numpy.zeros([8,8])
        gainrx[4:8,0:4] = 1
        
        justrx = 0    
    
    elif pattern == 102:
        # ON AXIS
        title=" ON_AXIS_eUP"
        
        ues = numpy.array([1,2,2,0])
        
        phase = numpy.array([[4,4,4,4,4,4,4,4],
                             [4,4,4,4,4,4,4,4],
                             [4,4,4,4,4,4,4,4],
                             [4,4,4,4,4,4,4,4],
                             [5,5,5,5,5,5,5,5],
                             [5,5,5,5,5,5,5,5],
                             [5,5,5,5,5,5,5,5],
                             [5,5,5,5,5,5,5,5]],dtype=float)

        gaintx = numpy.zeros([8,8])
        gaintx[0:4,4:8] = 1
        gaintx[4:8,0:4] = 1
        
        gainrx = numpy.zeros([8,8])
        gainrx[0:4,4:8] = 1
        
        justrx = 0    
        
    elif pattern == 103:
        # ON AXIS
        title=" ON_AXIS_sUP"
        
        ues = numpy.array([1,2,2,0])
        
        phase = numpy.array([[4,4,4,4,4,4,4,4],
                             [4,4,4,4,4,4,4,4],
                             [4,4,4,4,4,4,4,4],
                             [4,4,4,4,4,4,4,4],
                             [5,5,5,5,5,5,5,5],
                             [5,5,5,5,5,5,5,5],
                             [5,5,5,5,5,5,5,5],
                             [5,5,5,5,5,5,5,5]],dtype=float)

        gaintx = numpy.zeros([8,8])
        gaintx[0:4,4:8] = 1
        gaintx[4:8,0:4] = 1
        
        gainrx = numpy.zeros([8,8])
        gainrx[4:8,4:8] = 1
        
        justrx = 0            
    
    
    elif pattern == 104:
        # ON AXIS
        title=" ON_AXIS_sDW"
        
        ues = numpy.array([0,0,0,1])
        
        phase = numpy.array([[2,2,2,2,2,2,2,2],
                             [2,2,2,2,2,2,2,2],
                             [2,2,2,2,2,2,2,2],
                             [2,2,2,2,2,2,2,2],
                             [3,3,3,3,3,3,3,3],
                             [3,3,3,3,3,3,3,3],
                             [3,3,3,3,3,3,3,3],
                             [3,3,3,3,3,3,3,3]],dtype=float)

        gaintx = numpy.zeros([8,8])
        gaintx[0:4,4:8] = 1
        gaintx[4:8,0:4] = 1
        
        gainrx = numpy.zeros([8,8])
        gainrx[4:8,4:8] = 1
        
        justrx = 1

    elif pattern == 105:
        # POSITION 2 :
        # NDW : ON_AXIS
        # WDW : PERPENDICULAR B
        # EDW : ON_AXIS
        # SDW : ON_AXIS
        
        title=" ON_AXIS:TX nDW, RX wDW"
        
        ues = numpy.array([0,0,0,0])
        
        phase = numpy.array([[2,2,2,2,3,3,3,3],
                             [2,2,2,2,3,3,3,3],
                             [2,2,2,2,3,3,3,3],
                             [2,2,2,2,3,3,3,3],
                             [4,4,5,5,3,3,3,3],
                             [3,3,4,4,3,3,3,3],
                             [5,5,2,2,3,3,3,3],
                             [4,4,5,5,3,3,3,3]],dtype=float)

        gaintx = numpy.zeros([8,8])
        gaintx[0:4,0:4] = 1# N
        #gaintx[0:4,4:8] = 1# E
        #gaintx[4:8,0:4] = 1# W
        #gaintx[4:8,4:8] = 1# S
        
        gainrx = numpy.zeros([8,8])
        #gainrx[0:4,0:4] = 1# N
        #gainrx[0:4,4:8] = 1# E
        gainrx[4:8,0:4] = 1# W
        #gainrx[4:8,4:8] = 1# S
        
        justrx = 0

    elif pattern == 106:
        # POSITION 2 :
        # NDW : ON_AXIS
        # WDW : PERPENDICULAR B
        # EDW : ON_AXIS
        # SDW : ON_AXIS
        
        title=" ON_AXIS:TX nDW, RX eDW"
        
        ues = numpy.array([0,0,0,0])
        
        phase = numpy.array([[2,2,2,2,3,3,3,3],
                             [2,2,2,2,3,3,3,3],
                             [2,2,2,2,3,3,3,3],
                             [2,2,2,2,3,3,3,3],
                             [4,4,5,5,3,3,3,3],
                             [3,3,4,4,3,3,3,3],
                             [5,5,2,2,3,3,3,3],
                             [4,4,5,5,3,3,3,3]],dtype=float)

        gaintx = numpy.zeros([8,8])
        gaintx[0:4,0:4] = 1# N
        #gaintx[0:4,4:8] = 1# E
        #gaintx[4:8,0:4] = 1# W
        #gaintx[4:8,4:8] = 1# S
        
        gainrx = numpy.zeros([8,8])
        #gainrx[0:4,0:4] = 1# N
        gainrx[0:4,4:8] = 1# E
        #gainrx[4:8,0:4] = 1# W
        #gainrx[4:8,4:8] = 1# S
        
        justrx = 0

    elif pattern == 107:
        # POSITION 2 :
        # NDW : ON_AXIS
        # WDW : PERPENDICULAR B
        # EDW : ON_AXIS
        # SDW : ON_AXIS
        
        title=" ON_AXIS:TX nDW, RX sDW"
        
        ues = numpy.array([0,0,0,0])
        
        phase = numpy.array([[2,2,2,2,3,3,3,3],
                             [2,2,2,2,3,3,3,3],
                             [2,2,2,2,3,3,3,3],
                             [2,2,2,2,3,3,3,3],
                             [4,4,5,5,3,3,3,3],
                             [3,3,4,4,3,3,3,3],
                             [5,5,2,2,3,3,3,3],
                             [4,4,5,5,3,3,3,3]],dtype=float)

        gaintx = numpy.zeros([8,8])
        gaintx[0:4,0:4] = 1# N
        #gaintx[0:4,4:8] = 1# E
        #gaintx[4:8,0:4] = 1# W
        #gaintx[4:8,4:8] = 1# S
        
        gainrx = numpy.zeros([8,8])
        #gainrx[0:4,0:4] = 1# N
        #gainrx[0:4,4:8] = 1# E
        #gainrx[4:8,0:4] = 1# W
        gainrx[4:8,4:8] = 1# S
        
        justrx = 0

    elif pattern == 108:
        # POSITION 2 :
        # NUP : ON_AXIS
        # WUP : PERPENDICULAR B
        # EUP : ON_AXIS
        # SUP : ON_AXIS
        
        title=" ON_AXIS:TX nUP, RX wUP"
        
        ues = numpy.array([0,0,0,0])
        
        phase = numpy.array([[4,4,4,4,5,5,5,5],
                             [4,4,4,4,5,5,5,5],
                             [4,4,4,4,5,5,5,5],
                             [4,4,4,4,5,5,5,5],
                             [2,5,3,2,5,5,5,5],
                             [3,2,4,3,5,5,5,5],
                             [4,3,5,4,5,5,5,5],
                             [5,4,2,5,5,5,5,5]],dtype=float)

        gaintx = numpy.zeros([8,8])
        gaintx[0:4,0:4] = 1# N
        #gaintx[0:4,4:8] = 1# E
        #gaintx[4:8,0:4] = 1# W
        #gaintx[4:8,4:8] = 1# S
        
        gainrx = numpy.zeros([8,8])
        #gainrx[0:4,0:4] = 1# N
        #gainrx[0:4,4:8] = 1# E
        gainrx[4:8,0:4] = 1# W
        #gainrx[4:8,4:8] = 1# S
        
        justrx = 0

    elif pattern == 109:
        # POSITION 2 :
        # NUP : ON_AXIS
        # WUP : PERPENDICULAR B
        # EUP : ON_AXIS
        # SUP : ON_AXIS
        
        title=" ON_AXIS:TX nUP, RX sUP"
        
        ues = numpy.array([0,0,0,0])
        
        phase = numpy.array([[4,4,4,4,5,5,5,5],
                             [4,4,4,4,5,5,5,5],
                             [4,4,4,4,5,5,5,5],
                             [4,4,4,4,5,5,5,5],
                             [2,5,3,2,5,5,5,5],
                             [3,2,4,3,5,5,5,5],
                             [4,3,5,4,5,5,5,5],
                             [5,4,2,5,5,5,5,5]],dtype=float)

        gaintx = numpy.zeros([8,8])
        gaintx[0:4,0:4] = 1# N
        #gaintx[0:4,4:8] = 1# E
        #gaintx[4:8,0:4] = 1# W
        #gaintx[4:8,4:8] = 1# S
        
        gainrx = numpy.zeros([8,8])
        #gainrx[0:4,0:4] = 1# N
        #gainrx[0:4,4:8] = 1# E
        #gainrx[4:8,0:4] = 1# W
        gainrx[4:8,4:8] = 1# S
        
        justrx = 0

    elif pattern == 110:
        # POSITION 2 :
        # NUP : ON_AXIS
        # WUP : PERPENDICULAR B
        # EUP : ON_AXIS
        # SUP : ON_AXIS
        
        title=" ON_AXIS:TX nUP, RX eUP"
        
        ues = numpy.array([0,0,0,0])
        
        phase = numpy.array([[4,4,4,4,5,5,5,5],
                             [4,4,4,4,5,5,5,5],
                             [4,4,4,4,5,5,5,5],
                             [4,4,4,4,5,5,5,5],
                             [2,5,3,2,5,5,5,5],
                             [3,2,4,3,5,5,5,5],
                             [4,3,5,4,5,5,5,5],
                             [5,4,2,5,5,5,5,5]],dtype=float)

        gaintx = numpy.zeros([8,8])
        gaintx[0:4,0:4] = 1# N
        #gaintx[0:4,4:8] = 1# E
        #gaintx[4:8,0:4] = 1# W
        #gaintx[4:8,4:8] = 1# S
        
        gainrx = numpy.zeros([8,8])
        #gainrx[0:4,0:4] = 1# N
        gainrx[0:4,4:8] = 1# E
        #gainrx[4:8,0:4] = 1# W
        #gainrx[4:8,4:8] = 1# S
        
        justrx = 0
        
    elif pattern == 111:
        # POSITION 1 :
        # NDW : PERPENDICULAR B
        # WDW : PERPENDICULAR B
        # EDW : ON_AXIS
        # SDW : PERPENDICULAR B
        
        title=" PERPENDICULAR B:TX nDW, RX wDW"
        
        ues = numpy.array([0,0,0,0])
        
        phase = numpy.array([[2,2,3,3,3,3,3,3],
                             [5,5,2,2,3,3,3,3],
                             [3,3,4,4,3,3,3,3],
                             [2,2,3,3,3,3,3,3],
                             [4,4,5,5,4,4,5,5],
                             [3,3,4,4,3,3,4,4],
                             [5,5,2,2,5,5,2,2],
                             [4,4,5,5,4,4,5,5]],dtype=float)

        gaintx = numpy.zeros([8,8])
        gaintx[0:4,0:4] = 1# N
        #gaintx[0:4,4:8] = 1# E
        #gaintx[4:8,0:4] = 1# W
        #gaintx[4:8,4:8] = 1# S
        
        gainrx = numpy.zeros([8,8])
        #gainrx[0:4,0:4] = 1# N
        #gainrx[0:4,4:8] = 1# E
        gainrx[4:8,0:4] = 1# W
        #gainrx[4:8,4:8] = 1# S
        
        justrx = 0

    elif pattern == 112:
        # POSITION 1 :
        # NDW : PERPENDICULAR B
        # WDW : PERPENDICULAR B
        # EDW : ON_AXIS
        # SDW : PERPENDICULAR B
        
        title=" PERPENDICULAR B:TX nDW, RX eDW"
        
        ues = numpy.array([0,0,0,0])
        
        phase = numpy.array([[2,2,3,3,3,3,3,3],
                             [5,5,2,2,3,3,3,3],
                             [3,3,4,4,3,3,3,3],
                             [2,2,3,3,3,3,3,3],
                             [4,4,5,5,4,4,5,5],
                             [3,3,4,4,3,3,4,4],
                             [5,5,2,2,5,5,2,2],
                             [4,4,5,5,4,4,5,5]],dtype=float)

        gaintx = numpy.zeros([8,8])
        gaintx[0:4,0:4] = 1# N
        #gaintx[0:4,4:8] = 1# E
        #gaintx[4:8,0:4] = 1# W
        #gaintx[4:8,4:8] = 1# S
        
        gainrx = numpy.zeros([8,8])
        #gainrx[0:4,0:4] = 1# N
        gainrx[0:4,4:8] = 1# E
        #gainrx[4:8,0:4] = 1# W
        #gainrx[4:8,4:8] = 1# S
        
        justrx = 0

    elif pattern == 113:
        # POSITION 1 :
        # NDW : PERPENDICULAR B
        # WDW : PERPENDICULAR B
        # EDW : ON_AXIS
        # SDW : PERPENDICULAR B
        
        title=" PERPENDICULAR B:TX nDW, RX sDW"
        
        ues = numpy.array([0,0,0,0])
        
        phase = numpy.array([[2,2,3,3,3,3,3,3],
                             [5,5,2,2,3,3,3,3],
                             [3,3,4,4,3,3,3,3],
                             [2,2,3,3,3,3,3,3],
                             [4,4,5,5,4,4,5,5],
                             [3,3,4,4,3,3,4,4],
                             [5,5,2,2,5,5,2,2],
                             [4,4,5,5,4,4,5,5]],dtype=float)

        gaintx = numpy.zeros([8,8])
        gaintx[0:4,0:4] = 1# N
        #gaintx[0:4,4:8] = 1# E
        #gaintx[4:8,0:4] = 1# W
        #gaintx[4:8,4:8] = 1# S
        
        gainrx = numpy.zeros([8,8])
        #gainrx[0:4,0:4] = 1# N
        #gainrx[0:4,4:8] = 1# E
        #gainrx[4:8,0:4] = 1# W
        gainrx[4:8,4:8] = 1# S
        
        justrx = 0

    elif pattern == 114:
        # POSITION 1 :
        # NUP : PERPENDICULAR B
        # WUP : PERPENDICULAR B
        # EUP : ON_AXIS
        # SUP : PERPENDICULAR B
        
        title=" PERPENDICULAR B:TX nUP, RX wUP"
        
        ues = numpy.array([0,0,0,0])
        
        phase = numpy.array([[2,5,3,2,5,5,5,5],
                             [3,2,4,3,5,5,5,5],
                             [4,3,5,4,5,5,5,5],
                             [5,4,2,5,5,5,5,5],
                             [2,5,3,2,2,5,3,2],
                             [3,2,4,3,3,2,4,3],
                             [4,3,5,4,4,3,5,4],
                             [5,4,2,5,5,4,2,5]],dtype=float)

        gaintx = numpy.zeros([8,8])
        gaintx[0:4,0:4] = 1# N
        #gaintx[0:4,4:8] = 1# E
        #gaintx[4:8,0:4] = 1# W
        #gaintx[4:8,4:8] = 1# S
        
        gainrx = numpy.zeros([8,8])
        #gainrx[0:4,0:4] = 1# N
        #gainrx[0:4,4:8] = 1# E
        gainrx[4:8,0:4] = 1# W
        #gainrx[4:8,4:8] = 1# S
        
        justrx = 0

    elif pattern == 115:
        # POSITION 1 :
        # NUP : PERPENDICULAR B
        # WUP : PERPENDICULAR B
        # EUP : ON_AXIS
        # SUP : PERPENDICULAR B
        
        title=" PERPENDICULAR B:TX nUP, RX sUP"
        
        ues = numpy.array([0,0,0,0])
        
        phase = numpy.array([[2,5,3,2,5,5,5,5],
                             [3,2,4,3,5,5,5,5],
                             [4,3,5,4,5,5,5,5],
                             [5,4,2,5,5,5,5,5],
                             [2,5,3,2,2,5,3,2],
                             [3,2,4,3,3,2,4,3],
                             [4,3,5,4,4,3,5,4],
                             [5,4,2,5,5,4,2,5]],dtype=float)

        gaintx = numpy.zeros([8,8])
        gaintx[0:4,0:4] = 1# N
        #gaintx[0:4,4:8] = 1# E
        #gaintx[4:8,0:4] = 1# W
        #gaintx[4:8,4:8] = 1# S
        
        gainrx = numpy.zeros([8,8])
        #gainrx[0:4,0:4] = 1# N
        #gainrx[0:4,4:8] = 1# E
        #gainrx[4:8,0:4] = 1# W
        gainrx[4:8,4:8] = 1# S
        
        justrx = 0

    elif pattern == 116:
        # POSITION 1 :
        # NUP : PERPENDICULAR B
        # WUP : PERPENDICULAR B
        # EUP : ON_AXIS
        # SUP : PERPENDICULAR B
        
        title=" PERPENDICULAR B:TX nUP, RX eUP"
        
        ues = numpy.array([0,0,0,0])
        
        phase = numpy.array([[2,5,3,2,5,5,5,5],
                             [3,2,4,3,5,5,5,5],
                             [4,3,5,4,5,5,5,5],
                             [5,4,2,5,5,5,5,5],
                             [2,5,3,2,2,5,3,2],
                             [3,2,4,3,3,2,4,3],
                             [4,3,5,4,4,3,5,4],
                             [5,4,2,5,5,4,2,5]],dtype=float)

        gaintx = numpy.zeros([8,8])
        gaintx[0:4,0:4] = 1# N
        #gaintx[0:4,4:8] = 1# E
        #gaintx[4:8,0:4] = 1# W
        #gaintx[4:8,4:8] = 1# S
        
        gainrx = numpy.zeros([8,8])
        #gainrx[0:4,0:4] = 1# N
        gainrx[0:4,4:8] = 1# E
        #gainrx[4:8,0:4] = 1# W
        #gainrx[4:8,4:8] = 1# S
        
        justrx = 0

    elif pattern == 117:
        # POSITION 3 :
        # NDW : OFF PERPENDICULAR B
        # WDW : PERPENDICULAR B
        # EDW : ON_AXIS
        # SDW : OFF PERPENDICULAR B
        
        title=" OFF PERPENDICULAR B:TX nDW, RX wDW"
        
        ues = numpy.array([0,0,0,0])
        
        phase = numpy.array([[2,3,4,5,3,3,3,3],
                             [3,4,5,2,3,3,3,3],
                             [4,5,2,3,3,3,3,3],
                             [5,2,3,4,3,3,3,3],
                             [4,4,5,5,3,4,5,2],
                             [3,3,4,4,4,5,2,3],
                             [5,5,2,2,5,2,3,4],
                             [4,4,5,5,2,3,4,5]],dtype=float)

        gaintx = numpy.zeros([8,8])
        gaintx[0:4,0:4] = 1# N
        #gaintx[0:4,4:8] = 1# E
        #gaintx[4:8,0:4] = 1# W
        #gaintx[4:8,4:8] = 1# S
        
        gainrx = numpy.zeros([8,8])
        #gainrx[0:4,0:4] = 1# N
        #gainrx[0:4,4:8] = 1# E
        gainrx[4:8,0:4] = 1# W
        #gainrx[4:8,4:8] = 1# S
        
        justrx = 0

    elif pattern == 118:
        # POSITION 3 :
        # NDW : OFF PERPENDICULAR B
        # WDW : PERPENDICULAR B
        # EDW : ON_AXIS
        # SDW : OFF PERPENDICULAR B
        
        title=" OFF PERPENDICULAR B:TX nDW, RX eDW"
        
        ues = numpy.array([0,0,0,0])
        
        phase = numpy.array([[2,3,4,5,3,3,3,3],
                             [3,4,5,2,3,3,3,3],
                             [4,5,2,3,3,3,3,3],
                             [5,2,3,4,3,3,3,3],
                             [4,4,5,5,3,4,5,2],
                             [3,3,4,4,4,5,2,3],
                             [5,5,2,2,5,2,3,4],
                             [4,4,5,5,2,3,4,5]],dtype=float)

        gaintx = numpy.zeros([8,8])
        gaintx[0:4,0:4] = 1# N
        #gaintx[0:4,4:8] = 1# E
        #gaintx[4:8,0:4] = 1# W
        #gaintx[4:8,4:8] = 1# S
        
        gainrx = numpy.zeros([8,8])
        #gainrx[0:4,0:4] = 1# N
        gainrx[0:4,4:8] = 1# E
        #gainrx[4:8,0:4] = 1# W
        #gainrx[4:8,4:8] = 1# S
        
        justrx = 0

    elif pattern == 119:
        # POSITION 3 :
        # NDW : OFF PERPENDICULAR B
        # WDW : PERPENDICULAR B
        # EDW : ON_AXIS
        # SDW : OFF PERPENDICULAR B
        
        title=" OFF PERPENDICULAR B:TX nDW, RX sDW"
        
        ues = numpy.array([0,0,0,0])
        
        phase = numpy.array([[2,3,4,5,3,3,3,3],
                             [3,4,5,2,3,3,3,3],
                             [4,5,2,3,3,3,3,3],
                             [5,2,3,4,3,3,3,3],
                             [4,4,5,5,3,4,5,2],
                             [3,3,4,4,4,5,2,3],
                             [5,5,2,2,5,2,3,4],
                             [4,4,5,5,2,3,4,5]],dtype=float)

        gaintx = numpy.zeros([8,8])
        gaintx[0:4,0:4] = 1# N
        #gaintx[0:4,4:8] = 1# E
        #gaintx[4:8,0:4] = 1# W
        #gaintx[4:8,4:8] = 1# S
        
        gainrx = numpy.zeros([8,8])
        #gainrx[0:4,0:4] = 1# N
        #gainrx[0:4,4:8] = 1# E
        #gainrx[4:8,0:4] = 1# W
        gainrx[4:8,4:8] = 1# S
        
        justrx = 0

    elif pattern == 120:
        # POSITION 3 :
        # NUP : OFF PERPENDICULAR B
        # WUP : PERPENDICULAR B
        # EUP : ON_AXIS
        # SUP : OFF PERPENDICULAR B
        
        title=" OFF PERPENDICULAR B:TX nUP, RX wUP"
        
        ues = numpy.array([0,0,0,0])
        
        phase = numpy.array([[4,5,2,3,5,5,5,5],
                             [5,2,3,4,5,5,5,5],
                             [2,3,4,5,5,5,5,5],
                             [3,4,5,2,5,5,5,5],
                             [2,5,3,2,5,2,3,4],
                             [3,2,4,3,2,3,4,5],
                             [4,3,5,4,3,4,5,2],
                             [5,4,2,5,4,5,2,3]],dtype=float)

        gaintx = numpy.zeros([8,8])
        gaintx[0:4,0:4] = 1# N
        #gaintx[0:4,4:8] = 1# E
        #gaintx[4:8,0:4] = 1# W
        #gaintx[4:8,4:8] = 1# S
        
        gainrx = numpy.zeros([8,8])
        #gainrx[0:4,0:4] = 1# N
        #gainrx[0:4,4:8] = 1# E
        gainrx[4:8,0:4] = 1# W
        #gainrx[4:8,4:8] = 1# S
        
        justrx = 0

    elif pattern == 121:
        # POSITION 3 :
        # NUP : OFF PERPENDICULAR B
        # WUP : PERPENDICULAR B
        # EUP : ON_AXIS
        # SUP : OFF PERPENDICULAR B
        
        title=" OFF PERPENDICULAR B:TX nUP, RX sUP"
        
        ues = numpy.array([0,0,0,0])
        
        phase = numpy.array([[4,5,2,3,5,5,5,5],
                             [5,2,3,4,5,5,5,5],
                             [2,3,4,5,5,5,5,5],
                             [3,4,5,2,5,5,5,5],
                             [2,5,3,2,5,2,3,4],
                             [3,2,4,3,2,3,4,5],
                             [4,3,5,4,3,4,5,2],
                             [5,4,2,5,4,5,2,3]],dtype=float)

        gaintx = numpy.zeros([8,8])
        gaintx[0:4,0:4] = 1# N
        #gaintx[0:4,4:8] = 1# E
        #gaintx[4:8,0:4] = 1# W
        #gaintx[4:8,4:8] = 1# S
        
        gainrx = numpy.zeros([8,8])
        #gainrx[0:4,0:4] = 1# N
        #gainrx[0:4,4:8] = 1# E
        #gainrx[4:8,0:4] = 1# W
        gainrx[4:8,4:8] = 1# S
        
        justrx = 0

    elif pattern == 122:
        # POSITION 3 :
        # NUP : OFF PERPENDICULAR B
        # WUP : PERPENDICULAR B
        # EUP : ON_AXIS
        # SUP : OFF PERPENDICULAR B
        
        title=" OFF PERPENDICULAR B:TX nUP, RX eUP"
        
        ues = numpy.array([0,0,0,0])
        
        phase = numpy.array([[4,5,2,3,5,5,5,5],
                             [5,2,3,4,5,5,5,5],
                             [2,3,4,5,5,5,5,5],
                             [3,4,5,2,5,5,5,5],
                             [2,5,3,2,5,2,3,4],
                             [3,2,4,3,2,3,4,5],
                             [4,3,5,4,3,4,5,2],
                             [5,4,2,5,4,5,2,3]],dtype=float)

        gaintx = numpy.zeros([8,8])
        gaintx[0:4,0:4] = 1# N
        #gaintx[0:4,4:8] = 1# E
        #gaintx[4:8,0:4] = 1# W
        #gaintx[4:8,4:8] = 1# S
        
        gainrx = numpy.zeros([8,8])
        #gainrx[0:4,0:4] = 1# N
        gainrx[0:4,4:8] = 1# E
        #gainrx[4:8,0:4] = 1# W
        #gainrx[4:8,4:8] = 1# S
        
        justrx = 0
                        
    elif pattern == 123:
        # old julia :
        
        title=" 2008 JULIA"
        
        ues = numpy.array([0,0,0,0])
        
        phase = numpy.array([[2,3,3,3,3,4,4,4],
                             [5,2,2,2,2,3,3,3],
                             [3,4,4,4,4,5,5,5],
                             [2,3,3,3,3,4,4,4],
                             [4,5,5,5,5,2,2,2],
                             [3,4,4,4,4,5,5,5],
                             [5,2,2,2,2,3,3,3],
                             [4,5,5,5,5,2,2,2]],dtype=float)

        gaintx = numpy.zeros([8,8])
        gaintx[0:8,0:8] = 1
        #gaintx[0:4,0:4] = 1# N
        #gaintx[0:4,4:8] = 1# E
        #gaintx[4:8,0:4] = 1# W
        #gaintx[4:8,4:8] = 1# S
        
        gainrx = numpy.zeros([8,8])
        gainrx[0:8,0:8] = 1
        #gainrx[0:4,0:4] = 1# N
        #gainrx[0:4,4:8] = 1# E
        #gainrx[4:8,0:4] = 1# W
        #gainrx[4:8,4:8] = 1# S
        
        justrx = 0
                                
    elif pattern == 124:
        # abs_test :
        
        title=" Perpendicular to B, tx:  Wup, rx: Nup"
        
        ues = numpy.array([0,0,0,0])*2/3
        
        phase = numpy.array([[2,5,3,2,2,5,3,2],
                             [3,2,4,3,3,2,4,3],
                             [4,3,5,4,4,3,5,4],
                             [5,4,2,5,5,4,2,5],
                             [2,5,3,2,2,5,3,2],
                             [3,2,4,3,3,2,4,3],
                             [4,3,5,4,4,3,5,4],
                             [5,4,2,5,5,4,2,5]],dtype=float)

        gaintx = numpy.zeros([8,8])
        #gaintx[0:8,0:8] = 1
        #gaintx[0:4,0:4] = 1# N
        #gaintx[0:4,4:8] = 1# E
        gaintx[4:8,0:4] = 1# W
        #gaintx[4:8,4:8] = 1# S
        
        gainrx = numpy.zeros([8,8])
        #gainrx[0:8,0:8] = 1
        gainrx[0:4,0:4] = 1# N
        #gainrx[0:4,4:8] = 1# E
        #gainrx[4:8,0:4] = 1# W
        #gainrx[4:8,4:8] = 1# S
        
        justrx = 0
    
    elif pattern == 125:
        #  :
        
        title=" tx : Wup, rx: Eup"
        
        ues = numpy.array([0,0,0,0])*2/3    
        
        phase = numpy.array([[2,5,3,2,2,5,3,2],
                             [3,2,4,3,3,2,4,3],
                             [4,3,5,4,4,3,5,4],
                             [5,4,2,5,5,4,2,5],
                             [2,5,3,2,2,5,3,2],
                             [3,2,4,3,3,2,4,3],
                             [4,3,5,4,4,3,5,4],
                             [5,4,2,5,5,4,2,5]],dtype=float)

        gaintx = numpy.zeros([8,8])
        #gaintx[0:8,0:8] = 1
        #gaintx[0:4,0:4] = 1# N
        #gaintx[0:4,4:8] = 1# E
        gaintx[4:8,0:4] = 1# W
        #gaintx[4:8,4:8] = 1# S
        
        gainrx = numpy.zeros([8,8])
        #gainrx[0:8,0:8] = 1
        #gainrx[0:4,0:4] = 1# N
        gainrx[0:4,4:8] = 1# E
        #gainrx[4:8,0:4] = 1# W
        #gainrx[4:8,4:8] = 1# S
        
        justrx = 0

    elif pattern == 126:
        #  :
        
        title=" tx : Wup, rx: Sup"
        
        ues = numpy.array([0,0,0,0])*2/3     
        
        phase = numpy.array([[2,5,3,2,2,5,3,2],
                             [3,2,4,3,3,2,4,3],
                             [4,3,5,4,4,3,5,4],
                             [5,4,2,5,5,4,2,5],
                             [2,5,3,2,2,5,3,2],
                             [3,2,4,3,3,2,4,3],
                             [4,3,5,4,4,3,5,4],
                             [5,4,2,5,5,4,2,5]],dtype=float)

        gaintx = numpy.zeros([8,8])
        #gaintx[0:8,0:8] = 1
        #gaintx[0:4,0:4] = 1# N
        #gaintx[0:4,4:8] = 1# E
        gaintx[4:8,0:4] = 1# W
        #gaintx[4:8,4:8] = 1# S
        
        gainrx = numpy.zeros([8,8])
        #gainrx[0:8,0:8] = 1
        #gainrx[0:4,0:4] = 1# N
        #gainrx[0:4,4:8] = 1# E
        #gainrx[4:8,0:4] = 1# W
        gainrx[4:8,4:8] = 1# S
        
        justrx = 0

    elif pattern == 127:
        # abs_test :
        
        title=" tx : Wdw, rx: Ndw"
        
        ues = numpy.array([0,0,0,0])*2/3
        
        phase = numpy.array([[2,2,3,3,2,2,3,3],
                             [5,5,2,2,5,5,2,2],
                             [3,3,4,4,3,3,4,4],
                             [2,2,3,3,2,2,3,3],
                             [4,4,5,5,2,2,3,3],
                             [3,3,4,4,5,5,2,2],
                             [5,5,2,2,3,3,4,4],
                             [4,4,5,5,2,2,3,3]],dtype=float)

        gaintx = numpy.zeros([8,8])
        #gaintx[0:8,0:8] = 1
        #gaintx[0:4,0:4] = 1# N
        #gaintx[0:4,4:8] = 1# E
        gaintx[4:8,0:4] = 1# W
        #gaintx[4:8,4:8] = 1# S
        
        gainrx = numpy.zeros([8,8])
        #gainrx[0:8,0:8] = 1
        gainrx[0:4,0:4] = 1# N
        #gainrx[0:4,4:8] = 1# E
        #gainrx[4:8,0:4] = 1# W
        #gainrx[4:8,4:8] = 1# S
        
        justrx = 0
    
    elif pattern == 128:
        #  :
        
        title=" tx : Wdw, rx: Edw"
        
        ues = numpy.array([0,0,0,0])*2/3    
        
        phase = numpy.array([[2,2,3,3,2,2,3,3],
                             [5,5,2,2,5,5,2,2],
                             [3,3,4,4,3,3,4,4],
                             [2,2,3,3,2,2,3,3],
                             [4,4,5,5,2,2,3,3],
                             [3,3,4,4,5,5,2,2],
                             [5,5,2,2,3,3,4,4],
                             [4,4,5,5,2,2,3,3]],dtype=float)

        gaintx = numpy.zeros([8,8])
        #gaintx[0:8,0:8] = 1
        #gaintx[0:4,0:4] = 1# N
        #gaintx[0:4,4:8] = 1# E
        gaintx[4:8,0:4] = 1# W
        #gaintx[4:8,4:8] = 1# S
        
        gainrx = numpy.zeros([8,8])
        #gainrx[0:8,0:8] = 1
        #gainrx[0:4,0:4] = 1# N
        gainrx[0:4,4:8] = 1# E
        #gainrx[4:8,0:4] = 1# W
        #gainrx[4:8,4:8] = 1# S
        
        justrx = 0

    elif pattern == 129:
        #  :
        
        title=" tx : Wdw, rx: Sdw"
        
        ues = numpy.array([0,0,0,0])*2/3     
        
        phase = numpy.array([[2,2,3,3,2,2,3,3],
                             [5,5,2,2,5,5,2,2],
                             [3,3,4,4,3,3,4,4],
                             [2,2,3,3,2,2,3,3],
                             [4,4,5,5,2,2,3,3],
                             [3,3,4,4,5,5,2,2],
                             [5,5,2,2,3,3,4,4],
                             [4,4,5,5,2,2,3,3]],dtype=float)

        gaintx = numpy.zeros([8,8])
        #gaintx[0:8,0:8] = 1
        #gaintx[0:4,0:4] = 1# N
        #gaintx[0:4,4:8] = 1# E
        gaintx[4:8,0:4] = 1# W
        #gaintx[4:8,4:8] = 1# S
        
        gainrx = numpy.zeros([8,8])
        #gainrx[0:8,0:8] = 1
        #gainrx[0:4,0:4] = 1# N
        #gainrx[0:4,4:8] = 1# E
        #gainrx[4:8,0:4] = 1# W
        gainrx[4:8,4:8] = 1# S
        
        justrx = 0
        
    elif pattern == 130:
        # abs_test :
        
        title=" Off-Perpendicular, tx: Wup, rx: Nup"
        
        ues = numpy.array([0,0,0,0])*2/3
        
        phase = numpy.array([[5,2,3,4,2,5,3,2],
                             [2,3,4,5,3,2,4,3],
                             [3,4,5,2,4,3,5,4],
                             [4,5,2,3,5,4,2,5],
                             [2,5,3,2,5,2,3,4],
                             [3,2,4,3,2,3,4,5],
                             [4,3,5,4,3,4,5,2],
                             [5,4,2,5,4,5,2,3]],dtype=float)

        gaintx = numpy.zeros([8,8])
        #gaintx[0:8,0:8] = 1
        #gaintx[0:4,0:4] = 1# N
        #gaintx[0:4,4:8] = 1# E
        gaintx[4:8,0:4] = 1# W
        #gaintx[4:8,4:8] = 1# S
        
        gainrx = numpy.zeros([8,8])
        #gainrx[0:8,0:8] = 1
        gainrx[0:4,0:4] = 1# N
        #gainrx[0:4,4:8] = 1# E
        #gainrx[4:8,0:4] = 1# W
        #gainrx[4:8,4:8] = 1# S
        
        justrx = 0
    
    elif pattern == 131:
        #  :
        
        title=" tx : Wup, rx: Eup"
        
        ues = numpy.array([0,0,0,0])*2/3    
        
        phase = numpy.array([[5,2,3,4,2,5,3,2],
                             [2,3,4,5,3,2,4,3],
                             [3,4,5,2,4,3,5,4],
                             [4,5,2,3,5,4,2,5],
                             [2,5,3,2,5,2,3,4],
                             [3,2,4,3,2,3,4,5],
                             [4,3,5,4,3,4,5,2],
                             [5,4,2,5,4,5,2,3]],dtype=float)

        gaintx = numpy.zeros([8,8])
        #gaintx[0:8,0:8] = 1
        #gaintx[0:4,0:4] = 1# N
        #gaintx[0:4,4:8] = 1# E
        gaintx[4:8,0:4] = 1# W
        #gaintx[4:8,4:8] = 1# S
        
        gainrx = numpy.zeros([8,8])
        #gainrx[0:8,0:8] = 1
        #gainrx[0:4,0:4] = 1# N
        gainrx[0:4,4:8] = 1# E
        #gainrx[4:8,0:4] = 1# W
        #gainrx[4:8,4:8] = 1# S
        
        justrx = 0

    elif pattern == 132:
        #  :
        
        title=" tx : Wup, rx: Sup"
        
        ues = numpy.array([0,0,0,0])*2/3     
        
#         phase = numpy.array([[5,5,5,5,2,5,3,2],
#                              [5,5,5,5,3,2,4,3],
#                              [5,5,5,5,4,3,5,4],
#                              [5,5,5,5,5,4,2,5],
#                              [2,5,3,2,5,5,5,5],
#                              [3,2,4,3,5,5,5,5],
#                              [4,3,5,4,5,5,5,5],
#                              [5,4,2,5,5,5,5,5]],dtype=float)
        phase = numpy.array([[5,2,3,4,2,5,3,2],
                             [2,3,4,5,3,2,4,3],
                             [3,4,5,2,4,3,5,4],
                             [4,5,2,3,5,4,2,5],
                             [2,5,3,2,5,2,3,4],
                             [3,2,4,3,2,3,4,5],
                             [4,3,5,4,3,4,5,2],
                             [5,4,2,5,4,5,2,3]],dtype=float)

        gaintx = numpy.zeros([8,8])
        #gaintx[0:8,0:8] = 1
        #gaintx[0:4,0:4] = 1# N
        #gaintx[0:4,4:8] = 1# E
        gaintx[4:8,0:4] = 1# W
        #gaintx[4:8,4:8] = 1# S
        
        gainrx = numpy.zeros([8,8])
        #gainrx[0:8,0:8] = 1
        #gainrx[0:4,0:4] = 1# N
        #gainrx[0:4,4:8] = 1# E
        #gainrx[4:8,0:4] = 1# W
        gainrx[4:8,4:8] = 1# S
        
        justrx = 0

    elif pattern == 133:
        # abs_test :
        
        title=" tx : Wdw, rx: Ndw"
        
        ues = numpy.array([0,0,0,0])*2/3
        
#         phase = numpy.array([[3,3,3,3,2,2,3,3],
#                              [3,3,3,3,5,5,2,2],
#                              [3,3,3,3,3,3,4,4],
#                              [3,3,3,3,2,2,3,3],
#                              [4,4,5,5,3,3,3,3],
#                              [3,3,4,4,3,3,3,3],
#                              [5,5,2,2,3,3,3,3],
#                              [4,4,5,5,3,3,3,3]],dtype=float)
        phase = numpy.array([[3,4,5,2,2,2,3,3],
                             [4,5,2,3,5,5,2,2],
                             [5,2,3,4,3,3,4,4],
                             [2,3,4,5,2,2,3,3],
                             [4,4,5,5,3,4,5,2],
                             [3,3,4,4,4,5,2,3],
                             [5,5,2,2,5,2,3,4],
                             [4,4,5,5,2,3,4,5]],dtype=float)

        gaintx = numpy.zeros([8,8])
        #gaintx[0:8,0:8] = 1
        #gaintx[0:4,0:4] = 1# N
        #gaintx[0:4,4:8] = 1# E
        gaintx[4:8,0:4] = 1# W
        #gaintx[4:8,4:8] = 1# S
        
        gainrx = numpy.zeros([8,8])
        #gainrx[0:8,0:8] = 1
        gainrx[0:4,0:4] = 1# N
        #gainrx[0:4,4:8] = 1# E
        #gainrx[4:8,0:4] = 1# W
        #gainrx[4:8,4:8] = 1# S
        
        justrx = 0
    
    elif pattern == 134:
        #  :
        
        title=" tx : Wdw, rx: Edw"
        
        ues = numpy.array([0,0,0,0])*2/3    
        
#         phase = numpy.array([[3,3,3,3,2,2,3,3],
#                              [3,3,3,3,5,5,2,2],
#                              [3,3,3,3,3,3,4,4],
#                              [3,3,3,3,2,2,3,3],
#                              [4,4,5,5,3,3,3,3],
#                              [3,3,4,4,3,3,3,3],
#                              [5,5,2,2,3,3,3,3],
#                              [4,4,5,5,3,3,3,3]],dtype=float)
        phase = numpy.array([[3,4,5,2,2,2,3,3],
                             [4,5,2,3,5,5,2,2],
                             [5,2,3,4,3,3,4,4],
                             [2,3,4,5,2,2,3,3],
                             [4,4,5,5,3,4,5,2],
                             [3,3,4,4,4,5,2,3],
                             [5,5,2,2,5,2,3,4],
                             [4,4,5,5,2,3,4,5]],dtype=float)

        gaintx = numpy.zeros([8,8])
        #gaintx[0:8,0:8] = 1
        #gaintx[0:4,0:4] = 1# N
        #gaintx[0:4,4:8] = 1# E
        gaintx[4:8,0:4] = 1# W
        #gaintx[4:8,4:8] = 1# S
        
        gainrx = numpy.zeros([8,8])
        #gainrx[0:8,0:8] = 1
        #gainrx[0:4,0:4] = 1# N
        gainrx[0:4,4:8] = 1# E
        #gainrx[4:8,0:4] = 1# W
        #gainrx[4:8,4:8] = 1# S
        
        justrx = 0

    elif pattern == 135:
        #  :
        
        title=" tx : Wdw, rx: Sdw"
        
        ues = numpy.array([0,0,0,0])*2/3     
        
#         phase = numpy.array([[3,3,3,3,2,2,3,3],
#                              [3,3,3,3,5,5,2,2],
#                              [3,3,3,3,3,3,4,4],
#                              [3,3,3,3,2,2,3,3],
#                              [4,4,5,5,3,3,3,3],
#                              [3,3,4,4,3,3,3,3],
#                              [5,5,2,2,3,3,3,3],
#                              [4,4,5,5,3,3,3,3]],dtype=float)
        phase = numpy.array([[3,4,5,2,2,2,3,3],
                             [4,5,2,3,5,5,2,2],
                             [5,2,3,4,3,3,4,4],
                             [2,3,4,5,2,2,3,3],
                             [4,4,5,5,3,4,5,2],
                             [3,3,4,4,4,5,2,3],
                             [5,5,2,2,5,2,3,4],
                             [4,4,5,5,2,3,4,5]],dtype=float)

        gaintx = numpy.zeros([8,8])
        #gaintx[0:8,0:8] = 1
        #gaintx[0:4,0:4] = 1# N
        #gaintx[0:4,4:8] = 1# E
        gaintx[4:8,0:4] = 1# W
        #gaintx[4:8,4:8] = 1# S
        
        gainrx = numpy.zeros([8,8])
        #gainrx[0:8,0:8] = 1
        #gainrx[0:4,0:4] = 1# N
        #gainrx[0:4,4:8] = 1# E
        #gainrx[4:8,0:4] = 1# W
        gainrx[4:8,4:8] = 1# S
        
        justrx = 0
        
    elif pattern == 136:
        # abs_test :
        
        title=" tx : Wup, rx: Nup"
        
        ues = numpy.array([0,0,0,0])*2/3
        
        phase = numpy.array([[5,2,3,4,5,5,5,5],
                             [2,3,4,5,5,5,5,5],
                             [3,4,5,2,5,5,5,5],
                             [4,5,2,3,5,5,5,5],
                             [4,4,4,4,5,2,3,4],
                             [4,4,4,4,2,3,4,5],
                             [4,4,4,4,3,4,5,2],
                             [4,4,4,4,4,5,2,3]],dtype=float)

        gaintx = numpy.zeros([8,8])
        #gaintx[0:8,0:8] = 1
        #gaintx[0:4,0:4] = 1# N
        #gaintx[0:4,4:8] = 1# E
        gaintx[4:8,0:4] = 1# W
        #gaintx[4:8,4:8] = 1# S
        
        gainrx = numpy.zeros([8,8])
        #gainrx[0:8,0:8] = 1
        gainrx[0:4,0:4] = 1# N
        #gainrx[0:4,4:8] = 1# E
        #gainrx[4:8,0:4] = 1# W
        #gainrx[4:8,4:8] = 1# S
        
        justrx = 0

    elif pattern == 137:
        #  :
        
        title=" tx : Wup, rx: Eup"
        
        ues = numpy.array([0,0,0,0])*2/3    
        
        phase = numpy.array([[5,2,3,4,5,5,5,5],
                             [2,3,4,5,5,5,5,5],
                             [3,4,5,2,5,5,5,5],
                             [4,5,2,3,5,5,5,5],
                             [4,4,4,4,5,2,3,4],
                             [4,4,4,4,2,3,4,5],
                             [4,4,4,4,3,4,5,2],
                             [4,4,4,4,4,5,2,3]],dtype=float)

        gaintx = numpy.zeros([8,8])
        #gaintx[0:8,0:8] = 1
        #gaintx[0:4,0:4] = 1# N
        #gaintx[0:4,4:8] = 1# E
        gaintx[4:8,0:4] = 1# W
        #gaintx[4:8,4:8] = 1# S
        
        gainrx = numpy.zeros([8,8])
        #gainrx[0:8,0:8] = 1
        #gainrx[0:4,0:4] = 1# N
        gainrx[0:4,4:8] = 1# E
        #gainrx[4:8,0:4] = 1# W
        #gainrx[4:8,4:8] = 1# S
        
        justrx = 0

    elif pattern == 138:
        #  :
        
        title=" tx : Wdw, rx: Sdw"
        
        ues = numpy.array([0,0,0,0])*2/3     
        
        phase = numpy.array([[3,4,5,2,3,3,3,3],
                             [4,5,2,3,3,3,3,3],
                             [5,2,3,4,3,3,3,3],
                             [2,3,4,5,3,3,3,3],
                             [2,2,2,2,3,4,5,2],
                             [2,2,2,2,4,5,2,3],
                             [2,2,2,2,5,2,3,4],
                             [2,2,2,2,2,3,4,5]],dtype=float)

        gaintx = numpy.zeros([8,8])
        #gaintx[0:8,0:8] = 1
        #gaintx[0:4,0:4] = 1# N
        #gaintx[0:4,4:8] = 1# E
        gaintx[4:8,0:4] = 1# W
        #gaintx[4:8,4:8] = 1# S
        
        gainrx = numpy.zeros([8,8])
        #gainrx[0:8,0:8] = 1
        #gainrx[0:4,0:4] = 1# N
        #gainrx[0:4,4:8] = 1# E
        #gainrx[4:8,0:4] = 1# W
        gainrx[4:8,4:8] = 1# S
        
        justrx = 0
        
    elif pattern == 139:
        # abs_test :
        
        title=" On-axis, tx: Wup, rx:Nup"
        
        ues = numpy.array([0,0,0,0])*2/3
        
        phase = numpy.array([[5,5,5,5,5,5,5,5],
                             [2,5,5,5,5,5,5,5],
                             [5,5,5,5,5,5,5,5],
                             [5,5,5,5,5,5,5,5],
                             [4,4,4,4,5,5,5,5],
                             [4,4,4,4,5,5,5,5],
                             [4,4,4,4,5,5,5,5],
                             [4,4,4,4,5,5,5,5]],dtype=float)

        gaintx = numpy.zeros([8,8])
        #gaintx[0:8,0:8] = 1
        #gaintx[0:4,0:4] = 1# N
        #gaintx[0:4,4:8] = 1# E
        gaintx[4:8,0:4] = 1# W
        #gaintx[4:8,4:8] = 1# S
        
        gainrx = numpy.zeros([8,8])
        #gainrx[0:8,0:8] = 1
        gainrx[0:4,0:4] = 1# N
        #gainrx[0:4,4:8] = 1# E
        #gainrx[4:8,0:4] = 1# W
        #gainrx[4:8,4:8] = 1# S
        
        justrx = 0

    elif pattern == 140:
        # abs_test :
        
        title=" Off-Perpendicular, tx: Wup, rx:Nup"
        
        ues = numpy.array([0,0,0,0])*2/3
        
        phase = numpy.array([[5,2,3,4,5,5,5,5],
                             [2,3,4,5,5,5,5,5],
                             [3,4,5,2,5,5,5,5],
                             [4,5,2,3,5,5,5,5],
                             [4,4,4,4,5,2,3,4],
                             [4,4,4,4,2,3,4,5],
                             [4,4,4,4,3,4,5,2],
                             [4,4,4,4,4,5,2,3]],dtype=float)

        gaintx = numpy.zeros([8,8])
        #gaintx[0:8,0:8] = 1
        #gaintx[0:4,0:4] = 1# N
        #gaintx[0:4,4:8] = 1# E
        gaintx[4:8,0:4] = 1# W
        #gaintx[4:8,4:8] = 1# S
        
        gainrx = numpy.zeros([8,8])
        #gainrx[0:8,0:8] = 1
        gainrx[0:4,0:4] = 1# N
        #gainrx[0:4,4:8] = 1# E
        #gainrx[4:8,0:4] = 1# W
        #gainrx[4:8,4:8] = 1# S
        
        justrx = 0

    elif pattern == 141:
        # abs_test :
        
        title="On axis, tx: Wup, rx:Eup"
        
        ues = numpy.array([0,0,0,0])*2/3
        
        phase = numpy.array([[5,5,5,5,4,4,4,4],
                             [5,5,5,5,4,4,4,4],
                             [5,5,5,5,4,4,4,4],
                             [5,5,5,5,4,4,4,4],
                             [4,4,4,4,5,5,5,5],
                             [4,4,4,4,5,5,5,5],
                             [4,4,4,4,5,5,5,5],
                             [4,4,4,4,5,5,5,5]],dtype=float)

        gaintx = numpy.zeros([8,8])
        #gaintx[0:8,0:8] = 1
        #gaintx[0:4,0:4] = 1# N
        #gaintx[0:4,4:8] = 1# E
        gaintx[4:8,0:4] = 1# W
        #gaintx[4:8,4:8] = 1# S
        
        gainrx = numpy.zeros([8,8])
        #gainrx[0:8,0:8] = 1
        #gainrx[0:4,0:4] = 1# N
        gainrx[0:4,4:8] = 1# E
        #gainrx[4:8,0:4] = 1# W
        #gainrx[4:8,4:8] = 1# S
        
        justrx = 0
                    
    elif pattern == 142:
        # abs_test :
        
        title="Off-on axis, tx:nup, rx: sup"
        
        ues = numpy.array([0,0,0,0])*2/3
        
        phase = numpy.array([[4,5,2,3,4,4,4,4],
                             [5,2,3,4,4,4,4,4],
                             [2,3,4,5,4,4,4,4],
                             [3,4,5,2,4,4,4,4],
                             [4,4,4,4,5,2,3,4],
                             [4,4,4,4,2,3,4,5],
                             [4,4,4,4,3,4,5,2],
                             [4,4,4,4,4,5,2,3]],dtype=float)

        gaintx = numpy.zeros([8,8])
        #gaintx[0:8,0:8] = 1
        gaintx[0:4,0:4] = 1# N
        #gaintx[0:4,4:8] = 1# E
        #gaintx[4:8,0:4] = 1# W
        #gaintx[4:8,4:8] = 1# S
        
        gainrx = numpy.zeros([8,8])
        #gainrx[0:8,0:8] = 1
        #gainrx[0:4,0:4] = 1# N
        #gainrx[0:4,4:8] = 1# E
        #gainrx[4:8,0:4] = 1# W
        gainrx[4:8,4:8] = 1# S
        
        justrx = 0

    elif pattern == 143:
        # abs_test :
        
        title="Off-on axis, tx:nup, rx: wup"
        
        ues = numpy.array([0,0,0,0])*2/3
        
        phase = numpy.array([[4,5,2,3,4,4,4,4],
                             [5,2,3,4,4,4,4,4],
                             [2,3,4,5,4,4,4,4],
                             [3,4,5,2,4,4,4,4],
                             [4,4,4,4,5,2,3,4],
                             [4,4,4,4,2,3,4,5],
                             [4,4,4,4,3,4,5,2],
                             [4,4,4,4,4,5,2,3]],dtype=float)

        gaintx = numpy.zeros([8,8])
        #gaintx[0:8,0:8] = 1
        gaintx[0:4,0:4] = 1# N
        #gaintx[0:4,4:8] = 1# E
        #gaintx[4:8,0:4] = 1# W
        #gaintx[4:8,4:8] = 1# S
        
        gainrx = numpy.zeros([8,8])
        #gainrx[0:8,0:8] = 1
        #gainrx[0:4,0:4] = 1# N
        #gainrx[0:4,4:8] = 1# E
        gainrx[4:8,0:4] = 1# W
        #gainrx[4:8,4:8] = 1# S
        
        justrx = 0
        
    elif pattern == 144:
        # abs_test :
        
        title="Perpendicular B, tx:nup, rx: sup"
        
        ues = numpy.array([0,0,0,0])*2/3
        
        phase = numpy.array([[4,4,5,5,4,4,4,4],
                             [3,3,4,4,4,4,4,4],
                             [2,2,3,3,4,4,4,4],
                             [5,5,2,2,4,4,4,4],
                             [4,4,4,4,4,4,5,5],
                             [4,4,4,4,3,3,4,4],
                             [4,4,4,4,2,2,3,3],
                             [4,4,4,4,5,5,2,2]],dtype=float)

        gaintx = numpy.zeros([8,8])
        #gaintx[0:8,0:8] = 1
        gaintx[0:4,0:4] = 1# N
        #gaintx[0:4,4:8] = 1# E
        #gaintx[4:8,0:4] = 1# W
        #gaintx[4:8,4:8] = 1# S
        
        gainrx = numpy.zeros([8,8])
        #gainrx[0:8,0:8] = 1
        #gainrx[0:4,0:4] = 1# N
        #gainrx[0:4,4:8] = 1# E
        #gainrx[4:8,0:4] = 1# W
        gainrx[4:8,4:8] = 1# S
        
        justrx = 0

    elif pattern == 145:
        # JULIA para ABS
        title=" JULIA CHC : sDW"
        
        ues = numpy.array([0,0,0,0.75])*2/3
        
        phase = numpy.array([[5,4,3,2,5,4,3,2],
                             [2,5,4,3,2,5,4,3],
                             [2,5,4,3,2,5,4,3],
                             [3,2,5,4,3,2,5,4],
                             [5,4,3,2,3,2,5,4],
                             [2,5,4,3,4,3,2,5],
                             [2,5,4,3,4,3,2,5],
                             [3,2,5,4,5,4,3,2]],dtype=float)

        gaintx = numpy.zeros([8,8])
        #gaintx[0:8,0:8] = 1
        #gaintx[0:4,0:4] = 1# N
        gaintx[0:4,4:8] = 1# E
        #gaintx[4:8,0:4] = 1# W
        #gaintx[4:8,4:8] = 1# S
        
        gainrx = numpy.zeros([8,8])
        #gainrx[0:8,0:8] = 1
        #gainrx[0:4,0:4] = 1# N
        #gainrx[0:4,4:8] = 1# E
        #gainrx[4:8,0:4] = 1# W
        gainrx[4:8,4:8] = 1# S
        
        justrx = 0

    elif pattern == 146:
        # JULIA para ABS
        title=" JULIA CHF: sUP "
        
        ues = numpy.array([0,0,0,0.25])*2/3
        
#         phase = numpy.array([[4,4,5,5,4,3,4,2],
#                              [3,3,4,4,3,3,3,5],
#                              [2,2,3,3,3,3,3,2],
#                              [5,5,2,2,5,4,2,2],
#                              [4,3,4,2,2,5,2,2],
#                              [3,3,3,5,5,4,5,5],
#                              [3,3,3,2,3,3,4,4],
#                              [5,4,2,2,2,2,3,3]],dtype=float)
        phase = numpy.array([[4,3,4,3,4,3,4,3],
                             [3,3,3,3,3,3,3,3],
                             [5,3,3,2,5,3,3,2],
                             [5,4,2,2,5,4,2,2],
                             [4,3,4,3,4,3,4,3],
                             [3,3,3,3,3,3,3,3],
                             [5,3,3,2,5,3,3,2],
                             [5,4,2,2,5,4,2,2]],dtype=float)

        gaintx = numpy.zeros([8,8])
        #gaintx[0:8,0:8] = 1
        #gaintx[0:4,0:4] = 1# N
        gaintx[0:4,4:8] = 1# E
#         gaintx[0,7] = 0
#         gaintx[1,6] = 0
#         gaintx[1,7] = 0
#         gaintx[2,4] = 0
        #gaintx[4:8,0:4] = 1# W
        #gaintx[4:8,4:8] = 1# S
        
        gainrx = numpy.zeros([8,8])
        #gainrx[0:8,0:8] = 1
        #gainrx[0:4,0:4] = 1# N
        #gainrx[0:4,4:8] = 1# E
        #gainrx[4:8,0:4] = 1# W
        gainrx[4:8,4:8] = 1# S
        #gainrx[0,7] = 1
        
        justrx = 0
        
    elif pattern == 147:
        # LP FARADAY
        title=" chA : (4.5) all UP"
        
        ues = numpy.array([1,2,2,1])
        
        phase = numpy.array([[4,4,5,5,2,2,3,3],
                             [4,5,5,2,2,3,3,4],
                             [5,5,2,2,3,3,4,4],
                             [5,2,2,3,3,4,4,5],
                             [3,3,4,4,5,5,2,2],
                             [3,4,4,5,5,2,2,3],
                             [4,4,5,5,2,2,3,3],
                             [4,5,5,2,2,3,3,4]],dtype=float)

        gaintx = numpy.zeros([8,8])
        gaintx[0:8,0:8] = 1
        
        gainrx = numpy.zeros([8,8])
        gainrx[0:8,0:8] = 1
        
        justrx = 0

    elif pattern == 148:
        # LP FARADAY
        title=" chB : (4.5) all DW"
        
        ues = numpy.array([0,1,1,0])
        
        phase = numpy.array([[2,2,3,3,4,4,5,5],
                             [2,3,3,4,4,5,5,2],
                             [3,3,4,4,5,5,2,2],
                             [3,4,4,5,5,2,2,3],
                             [5,5,2,2,3,3,4,4],
                             [5,2,2,3,3,4,4,5],
                             [2,2,3,3,4,4,5,5],
                             [2,3,3,4,4,5,5,2]],dtype=float)

        gaintx = numpy.zeros([8,8])
        gaintx[0:8,0:8] = 1
        
        gainrx = numpy.zeros([8,8])
        gainrx[0:8,0:8] = 1
        
        justrx = 0
        
    elif pattern == 149:
        # Moon
        # 2015_10
#         title=" chC : Eup"
#         
#         ues = numpy.array([0,0,0,0])
#         
#         phase = numpy.array([[2,2,3,3,3,3,3,3],
#                              [2,3,3,4,3,3,3,3],
#                              [3,3,4,4,3,3,3,3],
#                              [3,4,4,5,3,3,3,3],
#                              [3,3,3,3,3,3,4,4],
#                              [3,3,3,3,3,4,4,5],
#                              [3,3,3,3,4,4,5,5],
#                              [3,3,3,3,4,5,5,2]],dtype=float)
# 
#         gaintx = numpy.zeros([8,8])
#         gaintx[0:4,4:8] = 1
#         
#         gainrx = numpy.zeros([8,8])
#         gainrx[0:4,4:8] = 1
#         
#         justrx = 0

        # 2017_07
#         title = " for Oblique ISR '4.5'"    
#         
# #        ues = numpy.array([1.,2.,2.,1.])
#         ues = numpy.array([1.5,3.0,3.0+0.15,1.5])*2/3
#         phase = numpy.array([[4,4,5,5,2,2,3,3],
#                           [4,5,5,2,2,3,3,4],
#                           [5,5,2,2,3,3,4,4],
#                           [5,2,2,3,3,4,4,5],
#                           [3,3,4,4,5,5,2,2],
#                           [3,4,4,5,5,2,2,3],
#                           [4,4,5,5,2,2,3,3],
#                           [4,5,5,2,2,3,3,4]],dtype=float)
#         
#         gaintx = numpy.zeros([8,8])
#         gaintx[0:4,4:8] = 1
# #        gaintx[:,:] = 1
#         
#         gainrx = numpy.zeros([8,8])
# #        gainrx[:,:] = 1
#         gainrx[0:4,4:8] = 1
# #        gainrx[4:8,0:4] = 1
# #        gainrx[0,0] = 1
# #        gainrx[7,7] = 1
#         
#         justrx = 0        

#         title = " for Oblique ISR '3.0N'"    
#         
#         ues = numpy.array([1.5,3.0,3.0,1.5])*2/3
# #        ues = numpy.array([0.0,0.0,0.0,0.0])*2/3
#         phase = numpy.array([[4,3,2,5,4,4,3,3],
#                           [3,2,5,4,4,3,3,2],
#                           [2,5,4,3,3,3,2,2],
#                           [5,4,3,2,3,2,2,5],
#                           [5,4,3,2,3,2,2,5],
#                           [4,3,2,5,4,3,2,5],
#                           [3,2,5,4,3,2,5,4],
#                           [2,5,4,3,2,5,4,3]],dtype=float)
#         
#         gaintx = numpy.zeros([8,8])
#         gaintx[0:4,4:8] = 1
#         
#         gainrx = numpy.zeros([8,8])
#         gainrx[0:4,4:8] = 1
#         
#         justrx = 0
# 3.66->3.63
# 4.66->4.89
# 4.33-> 4.26
# 3.33 -> 3.41
        title = " for Oblique ISR '3.0N'"    
         
        ues = numpy.array([1.5,3.0,3.0,1.5])*2/3
#        ues = numpy.array([0.0,0.0,0.0,0.0])*2/3
        phase = numpy.array([[5,3,2,5,5.00,4.89,4.26,4.00],
                          [3,2,5,4,4.89,4.26,4.00,3.63],
                          [2,5,4,3,4.26,4.00,3.63,3.41],
                          [5,4,3,2,4.00,3.63,3.41,3.00],
                          [5,4,3,2,3,2,2,5],
                          [4,3,2,5,4,3,2,5],
                          [3,2,5,4,3,2,5,4],
                          [2,5,4,3,2,5,4,3]],dtype=float)
         
        gaintx = numpy.zeros([8,8])
        gaintx[0:4,4:8] = 1
         
        gainrx = numpy.zeros([8,8])
        gainrx[0:4,4:8] = 1
         
        justrx = 0

#         title = " for Oblique ISR '3.0N'"    
#         
#         ues = numpy.array([1.5,3.0,3.0,1.5])*2/3
# #        ues = numpy.array([0.0,0.0,0.0,0.0])*2/3
#         phase = numpy.array([[5,3,2,5,5.00,4.00,4.00,4.00],
#                           [3,2,5,4,4.00,4.00,4.00,3.00],
#                           [2,5,4,3,4.00,4.00,3.00,3.00],
#                           [5,4,3,2,4.00,3.00,3.00,3.00],
#                           [5,4,3,2,3,2,2,5],
#                           [4,3,2,5,4,3,2,5],
#                           [3,2,5,4,3,2,5,4],
#                           [2,5,4,3,2,5,4,3]],dtype=float)
#         
#         gaintx = numpy.zeros([8,8])
#         gaintx[0:4,4:8] = 1
#         
#         gainrx = numpy.zeros([8,8])
#         gainrx[0:4,4:8] = 1
#         
#         justrx = 0
                
    elif pattern == 150:
        # Moon
        # 2015_10
        # 2017_07
        title=" chG : Wup"
        
        ues = numpy.array([0,0,0,0])
        
        phase = numpy.array([[2,2,3,3,3,3,3,3],
                             [2,3,3,4,3,3,3,3],
                             [3,3,4,4,3,3,3,3],
                             [3,4,4,5,3,3,3,3],
                             [3,3,3,3,3,3,4,4],
                             [3,3,3,3,3,4,4,5],
                             [3,3,3,3,4,4,5,5],
                             [3,3,3,3,4,5,5,2]],dtype=float)

        gaintx = numpy.zeros([8,8])
        gaintx[0:4,4:8] = 1
        
        gainrx = numpy.zeros([8,8])
        gainrx[4:8,0:4] = 1
        
        justrx = 0
    
    elif pattern == 151:
        # Moon
        # 2015_10
        # 2017_07
        title=" chH : Wdw"
        
        ues = numpy.array([0,0,0,0])
        
        phase = numpy.array([[2,2,3,3,4,4,4,4],
                             [2,3,3,4,4,4,4,4],
                             [3,3,4,4,4,4,4,4],
                             [3,4,4,5,4,4,4,4],
                             [4,4,4,4,3,3,4,4],
                             [4,4,4,4,3,4,4,5],
                             [4,4,4,4,4,4,5,5],
                             [4,4,4,4,4,5,5,2]],dtype=float)

        gaintx = numpy.zeros([8,8])
        gaintx[0:4,4:8] = 1
        
        gainrx = numpy.zeros([8,8])
        gainrx[4:8,0:4] = 1
        
        justrx = 0
        
    elif pattern == 152:
        # Moon
        # 2015_10
        # 2017_07
        title=" chB : N1dw"
        
        ues = numpy.array([0,0,0,0])
        
        phase = numpy.array([[4,2,3,3,4,4,4,4],
                             [2,3,3,4,4,4,4,4],
                             [3,3,4,4,4,4,4,4],
                             [3,4,4,5,4,4,4,4],
                             [4,4,4,4,3,3,4,4],
                             [4,4,4,4,3,4,4,5],
                             [4,4,4,4,4,4,5,5],
                             [4,4,4,4,4,5,5,2]],dtype=float)

        gaintx = numpy.zeros([8,8])
        gaintx[0:4,4:8] = 1
        
        gainrx = numpy.zeros([8,8])
        gainrx[0,0] = 1
        
        justrx = 0
    
    elif pattern == 153:
        # Moon
        # 2015_10
        # 2017_07
        title=" chE : S16up"
        
        ues = numpy.array([0,0,0,0])
        
        phase = numpy.array([[4,2,3,3,3,3,3,3],
                             [2,3,3,4,3,3,3,3],
                             [3,3,4,4,3,3,3,3],
                             [3,4,4,5,3,3,3,3],
                             [3,3,3,3,3,3,4,4],
                             [3,3,3,3,3,4,4,5],
                             [3,3,3,3,4,4,5,5],
                             [3,3,3,3,4,5,5,3]],dtype=float)

        gaintx = numpy.zeros([8,8])
        gaintx[0:4,4:8] = 1
        
        gainrx = numpy.zeros([8,8])
        gainrx[7,7] = 1
        
        justrx = 0
        
    elif pattern == 154:
        # 
        title=" Tx all up, Rx nup"
        
        ues = numpy.array([1,0.5,0,0.5])
        
        phase = numpy.array([[3,4,5,2,3,4,5,2],
                             [4,5,2,3,4,5,2,3],
                             [5,2,3,4,5,2,3,4],
                             [2,3,4,5,2,3,4,5],
                             [2,3,4,5,2,3,4,5],
                             [3,4,5,2,3,4,5,2],
                             [4,5,2,3,4,5,2,3],
                             [5,2,3,4,5,2,3,4]],dtype=float)

        gaintx = numpy.zeros([8,8])
        gaintx[0:8,0:8] = 1
        
        gainrx = numpy.zeros([8,8])
        gainrx[0:4,0:4] = 1
        
        justrx = 0
    
    elif pattern == 155:
        # 
        title=" Tx all down, Rx wdw "
        
        ues = numpy.array([2,1.5,1,1.5])
        
        phase = numpy.array([[5,2,3,4,5,2,3,4],
                             [2,3,4,5,2,3,4,5],
                             [3,4,5,2,3,4,5,2],
                             [4,5,2,3,4,5,2,3],
                             [4,5,2,3,4,5,2,3],
                             [5,2,3,4,5,2,3,4],
                             [2,3,4,5,2,3,4,5],
                             [3,4,5,2,3,4,5,2]],dtype=float)

        gaintx = numpy.zeros([8,8])
        gaintx[0:8,0:8] = 1
        
        gainrx = numpy.zeros([8,8])
        gainrx[4:8,0:4] = 1
        
        justrx = 0
 
    elif pattern == 156:
        # 
        title=" Tx all up, Rx nup"
        
        ues = numpy.array([1,0.5,0,0.5])
        
        phase = numpy.array([[3,4,5,2,3,4,5,2],
                             [4,5,2,3,4,5,2,3],
                             [5,2,3,4,5,2,3,4],
                             [2,3,4,5,2,3,4,5],
                             [2,3,4,5,2,3,4,5],
                             [3,4,5,2,3,4,5,2],
                             [4,5,2,3,4,5,2,3],
                             [5,2,3,4,5,2,3,4]],dtype=float)

        gaintx = numpy.zeros([8,8])
        gaintx[0:8,0:8] = 1
        
        gainrx = numpy.zeros([8,8])
        gainrx[0:4,0:4] = 1
        
        justrx = 0
    
    elif pattern == 157:
        # 
        title=" Tx all down, Rx all down "
        
        ues = numpy.array([2,1.5,1,1.5])
        ues = numpy.array([2,0.5,1,1.5])

        phase = numpy.array([[5,2,3,4,5,2,3,4],
                             [2,3,4,5,2,3,4,5],
                             [3,4,5,2,3,4,5,2],
                             [4,5,2,3,4,5,2,3],
                             [4,5,2,3,4,5,2,3],
                             [5,2,3,4,5,2,3,4],
                             [2,3,4,5,2,3,4,5],
                             [3,4,5,2,3,4,5,2]],dtype=float)
                
        phase = numpy.array([[2,3,4,5,5,2,3,4],
                             [3,4,5,2,2,3,4,5],
                             [4,5,2,3,3,4,5,2],
                             [5,2,3,4,4,5,2,3],
                             [4,5,2,3,4,5,2,3],
                             [5,2,3,4,5,2,3,4],
                             [2,3,4,5,2,3,4,5],
                             [3,4,5,2,3,4,5,2]],dtype=float)

        gaintx = numpy.zeros([8,8])
        gaintx[0:8,0:8] = 1
        
        gainrx = numpy.zeros([8,8])
        gainrx[0:8,0:8] = 1
        
        justrx = 0
        
    elif pattern == 158:
        # 
        title=" chA : tx down all, rx wdw "
        
        ues = numpy.array([2,1.5,1,1.5])
        
        phase = numpy.array([[4,4,4,4,4,4,4,4],
                             [4,4,4,4,4,4,4,4],
                             [4,4,4,4,4,4,4,4],
                             [4,4,4,4,4,4,4,4],
                             [5,5,5,5,5,5,5,5],
                             [5,5,5,5,5,5,5,5],
                             [5,5,5,5,5,5,5,5],
                             [5,5,5,5,5,5,5,5]],dtype=float)

        gaintx = numpy.zeros([8,8])
        gaintx[0:8,0:8] = 1
        
        gainrx = numpy.zeros([8,8])
        gainrx[0:8,0:8] = 1
        
        justrx = 0

    elif pattern == 159:
        # 
        title=" chB : tx down all, rx wdw "
        
        ues = numpy.array([2,1.5,1,1.5])
        
        phase = numpy.array([[2,2,2,2,2,2,2,2],
                             [2,2,2,2,2,2,2,2],
                             [2,2,2,2,2,2,2,2],
                             [2,2,2,2,2,2,2,2],
                             [3,3,3,3,3,3,3,3],
                             [3,3,3,3,3,3,3,3],
                             [3,3,3,3,3,3,3,3],
                             [3,3,3,3,3,3,3,3]],dtype=float)

        gaintx = numpy.zeros([8,8])
        gaintx[0:8,0:8] = 1
        
        gainrx = numpy.zeros([8,8])
        gainrx[0:8,0:4] = 1
        
        justrx = 0

    elif pattern == 160:
        # MAGNETOPAUSE 2016_02
        title=" Tx N and S down , Rx Ndw "
        
        ues = numpy.array([1.0,2.0,2.0,1.0])
        
        phase = numpy.array([[2.0,3.0,4.0,5.0,0.5,0.5,0.5,0.5],
                             [3.0,4.0,5.0,2.0,0.5,0.5,0.5,0.5],
                             [4.0,5.0,2.0,3.0,0.5,0.5,0.5,0.5],
                             [5.0,2.0,3.0,4.0,0.5,0.5,0.5,0.5],
                             [0.5,0.5,0.5,0.5,3.0,4.0,5.0,2.0],
                             [0.5,0.5,0.5,0.5,4.0,5.0,2.0,3.0],
                             [0.5,0.5,0.5,0.5,5.0,2.0,3.0,4.0],
                             [0.5,0.5,0.5,0.5,2.0,3.0,4.0,5.0]],dtype=float)

        gaintx = numpy.zeros([8,8])
        gaintx[0:4,0:4] = 1
        gaintx[4:8,4:8] = 1
        
        gainrx = numpy.zeros([8,8])
        gainrx[0:4,0:4] = 1
        
        justrx = 0

    elif pattern == 161:
        # MAGNETOPAUSE 2016_02
        title=" Tx N and S up , Rx Nup "
        
        ues = numpy.array([0.0,1.0,1.0,0.0])
        
        phase = numpy.array([[0.0,1.0,2.0,3.0,0.5,0.5,0.5,0.5],
                             [1.0,2.0,3.0,0.0,0.5,0.5,0.5,0.5],
                             [2.0,3.0,0.0,1.0,0.5,0.5,0.5,0.5],
                             [3.0,0.0,1.0,2.0,0.5,0.5,0.5,0.5],
                             [0.5,0.5,0.5,0.5,1.0,2.0,3.0,0.0],
                             [0.5,0.5,0.5,0.5,2.0,3.0,0.0,5.0],
                             [0.5,0.5,0.5,0.5,3.0,0.0,1.0,2.0],
                             [0.5,0.5,0.5,0.5,0.0,1.0,2.0,3.0]],dtype=float)

        gaintx = numpy.zeros([8,8])
        gaintx[0:4,0:4] = 1
        gaintx[4:8,4:8] = 1
        
        gainrx = numpy.zeros([8,8])
        gainrx[0:4,0:4] = 1
        
        justrx = 0

    elif pattern == 162:
        # MAGNETOPAUSE 2016_02
        title=" Beam 1 : Tx N and S down , Rx Sdw "
        
        ues = numpy.array([1.0,2.0,2.0,1.0])
        
        phase = numpy.array([[2.0,3.0,4.0,5.0,0.5,0.5,0.5,0.5],
                             [3.0,4.0,5.0,2.0,0.5,0.5,0.5,0.5],
                             [4.0,5.0,2.0,3.0,0.5,0.5,0.5,0.5],
                             [5.0,2.0,3.0,4.0,0.5,0.5,0.5,0.5],
                             [0.5,0.5,0.5,0.5,3.0,4.0,5.0,2.0],
                             [0.5,0.5,0.5,0.5,4.0,5.0,2.0,3.0],
                             [0.5,0.5,0.5,0.5,5.0,2.0,3.0,4.0],
                             [0.5,0.5,0.5,0.5,2.0,3.0,4.0,5.0]],dtype=float)

        gaintx = numpy.zeros([8,8])
        gaintx[0:4,0:4] = 1
        gaintx[4:8,4:8] = 1
        
        gainrx = numpy.zeros([8,8])
        gainrx[4:8,4:8] = 1
        
        justrx = 0

    elif pattern == 163:
        # MAGNETOPAUSE 2016_02
        title=" Beam 2 : Tx N and S down , Rx Sdw "
        
        ues = numpy.array([1.0,2.0,2.0,1.0])
        
        phase = numpy.array([[2.0,2.0,3.0,3.0,0.5,0.5,0.5,0.5],
                             [2.0,3.0,3.0,0.0,0.5,0.5,0.5,0.5],
                             [3.0,3.0,0.0,0.0,0.5,0.5,0.5,0.5],
                             [3.0,0.0,0.0,1.0,0.5,0.5,0.5,0.5],
                             [0.5,0.5,0.5,0.5,1.0,1.0,2.0,2.0],
                             [0.5,0.5,0.5,0.5,1.0,2.0,2.0,3.0],
                             [0.5,0.5,0.5,0.5,2.0,2.0,3.0,3.0],
                             [0.5,0.5,0.5,0.5,2.0,3.0,3.0,0.0]],dtype=float)

        gaintx = numpy.zeros([8,8])
        gaintx[0:4,0:4] = 1
        gaintx[4:8,4:8] = 1
        
        gainrx = numpy.zeros([8,8])
        gainrx[4:8,4:8] = 1
        
        justrx = 0
        
    elif pattern == 164:
        # MAGNETOPAUSE 2016_02
        title=" Beam 3 : Tx N and S down , Rx Sdw "
        
        ues = numpy.array([1.0,2.0,2.0,1.0])
        
        phase = numpy.array([[0.0,0.0,0.0,0.0,0.5,0.5,0.5,0.5],
                             [0.0,0.0,0.0,0.0,0.5,0.5,0.5,0.5],
                             [0.0,0.0,0.0,0.0,0.5,0.5,0.5,0.5],
                             [0.0,0.0,0.0,0.0,0.5,0.5,0.5,0.5],
                             [0.5,0.5,0.5,0.5,1.0,1.0,1.0,1.0],
                             [0.5,0.5,0.5,0.5,1.0,1.0,1.0,1.0],
                             [0.5,0.5,0.5,0.5,1.0,1.0,1.0,1.0],
                             [0.5,0.5,0.5,0.5,1.0,1.0,1.0,1.0]],dtype=float)

        gaintx = numpy.zeros([8,8])
        gaintx[0:4,0:4] = 1
        gaintx[4:8,4:8] = 1
        
        gainrx = numpy.zeros([8,8])
        gainrx[4:8,4:8] = 1
        
        justrx = 0
        
    elif pattern == 165:
        # MAGNETOPAUSE 2016_02
        title=" Beam 4 : Tx N and S down , Rx Sdw "
        
        ues = numpy.array([1.0,2.0,2.0,1.0])
        
        phase = numpy.array([[2.0,2.0,1.0,1.0,0.5,0.5,0.5,0.5],
                             [2.0,1.0,1.0,0.0,0.5,0.5,0.5,0.5],
                             [1.0,1.0,0.0,0.0,0.5,0.5,0.5,0.5],
                             [1.0,0.0,0.0,3.0,0.5,0.5,0.5,0.5],
                             [0.5,0.5,0.5,0.5,1.0,1.0,0.0,0.0],
                             [0.5,0.5,0.5,0.5,1.0,0.0,0.0,3.0],
                             [0.5,0.5,0.5,0.5,0.0,0.0,3.0,3.0],
                             [0.5,0.5,0.5,0.5,0.0,3.0,3.0,2.0]],dtype=float)

        gaintx = numpy.zeros([8,8])
        gaintx[0:4,0:4] = 1
        gaintx[4:8,4:8] = 1
        
        gainrx = numpy.zeros([8,8])
        gainrx[4:8,4:8] = 1
        
        justrx = 0

    elif pattern == 166:
        # MAGNETOPAUSE 2016_02
        title=" Beam 5 : Tx N and S down , Rx Sdw "
        
        ues = numpy.array([1.0,2.0,2.0,1.0])
        
        phase = numpy.array([[2.0,1.0,0.0,3.0,0.5,0.5,0.5,0.5],
                             [1.0,0.0,3.0,2.0,0.5,0.5,0.5,0.5],
                             [0.0,3.0,2.0,1.0,0.5,0.5,0.5,0.5],
                             [3.0,2.0,1.0,0.0,0.5,0.5,0.5,0.5],
                             [0.5,0.5,0.5,0.5,3.0,2.0,1.0,0.0],
                             [0.5,0.5,0.5,0.5,2.0,1.0,0.0,3.0],
                             [0.5,0.5,0.5,0.5,1.0,0.0,3.0,2.0],
                             [0.5,0.5,0.5,0.5,0.0,3.0,2.0,1.0]],dtype=float)

        gaintx = numpy.zeros([8,8])
        gaintx[0:4,0:4] = 1
        gaintx[4:8,4:8] = 1
        
        gainrx = numpy.zeros([8,8])
        gainrx[4:8,4:8] = 1
        
        justrx = 0

    elif pattern == 167:
        # MST
        title=" MST : Tx All dw , Rx Wdw "
        
        ues = numpy.array([3.4,0.0,0.0,0.0])*2/3
        ues = numpy.array([0.0,0.0,0.0,0.0])*2/3
        ues = numpy.array([3.0,0.0,1.5,0.0])*2/3
        
        phase = numpy.array([[4.0,5.0,2.0,3.0,2.0,3.0,3.0,3.0],
                             [5.0,2.0,3.0,4.0,5.0,2.0,2.0,2.0],
                             [2.0,3.0,4.0,5.0,3.0,4.0,4.0,4.0],
                             [3.0,4.0,5.0,2.0,2.0,3.0,3.0,3.0],
                             [4.0,3.0,5.0,4.0,4.89,4.89,4.89,4.89],
                             [5.0,4.0,2.0,5.0,4.26,4.26,4.26,4.26],
                             [5.0,4.0,2.0,5.0,3.63,3.63,3.63,3.63],
                             [2.0,5.0,3.0,2.0,3.0,3.0,3.0,3.0]],dtype=float)

        gaintx = numpy.zeros([8,8])+1
        
        gainrx = numpy.zeros([8,8])
        gainrx[4:8,0:4] = 1
        
        justrx = 0
                          
    elif pattern == 168:
        # MST
        title=" MST : Tx All dw , Rx Edw "
        
        ues = numpy.array([3.4,0.0,0.0,0.0])*2/3
        #ues = numpy.array([3.4,0.0,1.6,0.0])*2/3
        ues = numpy.array([0.0,0.0,0.0,0.0])*2/3
        ues = numpy.array([0.0,0.0,3.0,0.0])*2/3
        ues = numpy.array([3.0,0.0,1.5,0.0])*2/3
        
        phase = numpy.array([[4.0,5.0,2.0,3.0,2.0,3.0,3.0,3.0],
                             [5.0,2.0,3.0,4.0,5.0,2.0,2.0,2.0],
                             [2.0,3.0,4.0,5.0,3.0,4.0,4.0,4.0],
                             [3.0,4.0,5.0,2.0,2.0,3.0,3.0,3.0],
                             [4.0,3.0,5.0,4.0,4.89,4.89,4.89,4.89],
                             [5.0,4.0,2.0,5.0,4.26,4.26,4.26,4.26],
                             [5.0,4.0,2.0,5.0,3.63,3.63,3.63,3.63],
                             [2.0,5.0,3.0,2.0,3.0,3.0,3.0,3.0]],dtype=float)

        gaintx = numpy.zeros([8,8])+1
        
        gainrx = numpy.zeros([8,8])
        gainrx[0:4,4:8] = 1
        
        justrx = 0

    elif pattern == 169:
        # MST
        title=" MST : Tx All dw , Rx Sdw "
        
        ues = numpy.array([3.4,0.0,1.6,0.0])*2/3
        ues = numpy.array([0.0,0.0,0.0,0.0])*2/3
        ues = numpy.array([0.0,0.0,3.0,0.0])*2/3
        ues = numpy.array([3.0,0.0,1.5,0.0])*2/3
        
        phase = numpy.array([[4.0,5.0,2.0,3.0,2.0,3.0,3.0,3.0],
                             [5.0,2.0,3.0,4.0,5.0,2.0,2.0,2.0],
                             [2.0,3.0,4.0,5.0,3.0,4.0,4.0,4.0],
                             [3.0,4.0,5.0,2.0,2.0,3.0,3.0,3.0],
                             [4.0,3.0,5.0,4.0,4.89,4.89,4.89,4.89],
                             [5.0,4.0,2.0,5.0,4.26,4.26,4.26,4.26],
                             [5.0,4.0,2.0,5.0,3.63,3.63,3.63,3.63],
                             [2.0,5.0,3.0,2.0,3.0,3.0,3.0,3.0]],dtype=float)

        gaintx = numpy.zeros([8,8])+1
        
        gainrx = numpy.zeros([8,8])
        gainrx[4:8,4:8] = 1
        
        justrx = 0
    
    elif pattern == 170:
        # MST
        title=" MST : Tx All dw , Rx Ndw "
        
        ues = numpy.array([3.4,0.0,0.0,0.0])*2/3
        ues = numpy.array([0.0,0.0,3.0,0.0])*2/3
        ues = numpy.array([3.0,0.0,1.5,0.0])*2/3
        
        phase = numpy.array([[4.0,5.0,2.0,3.0,2.0,3.0,3.0,3.0],
                             [5.0,2.0,3.0,4.0,5.0,2.0,2.0,2.0],
                             [2.0,3.0,4.0,5.0,3.0,4.0,4.0,4.0],
                             [3.0,4.0,5.0,2.0,2.0,3.0,3.0,3.0],
                             [4.0,3.0,5.0,4.0,4.89,4.89,4.89,4.89],
                             [5.0,4.0,2.0,5.0,4.26,4.26,4.26,4.26],
                             [5.0,4.0,2.0,5.0,3.63,3.63,3.63,3.63],
                             [2.0,5.0,3.0,2.0,3.0,3.0,3.0,3.0]],dtype=float)

        gaintx = numpy.zeros([8,8])+1
        
        gainrx = numpy.zeros([8,8])
        gainrx[0:4,0:4] = 1
        
        justrx = 0
    
    elif pattern == 171:
        # Meteors
        title=" Meteors : Tx Wup , Rx Nup "
        
        ues = numpy.array([0.0,0.0,0.0,0.0])*2/3
        
        phase = numpy.array([[5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0],
                             [5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0],
                             [5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0],
                             [5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0],
                             [4.0,4.0,4.0,4.0,4.0,4.0,4.0,4.0],
                             [4.0,4.0,4.0,4.0,4.0,4.0,4.0,4.0],
                             [4.0,4.0,4.0,4.0,4.0,4.0,4.0,4.0],
                             [4.0,4.0,4.0,4.0,4.0,4.0,4.0,4.0]],dtype=float)

        gaintx = numpy.zeros([8,8])
        gaintx[4:8,0:4] = 1
        
        gainrx = numpy.zeros([8,8])
        gainrx[0:4,0:4] = 1
        
        justrx = 0
    
    elif pattern == 172:
        # Meteors
        title=" Meteors : Tx Wup , Rx Eup "
        
        ues = numpy.array([0.0,0.0,0.0,0.0])*2/3
        
        phase = numpy.array([[5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0],
                             [5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0],
                             [5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0],
                             [5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0],
                             [4.0,4.0,4.0,4.0,4.0,4.0,4.0,4.0],
                             [4.0,4.0,4.0,4.0,4.0,4.0,4.0,4.0],
                             [4.0,4.0,4.0,4.0,4.0,4.0,4.0,4.0],
                             [4.0,4.0,4.0,4.0,4.0,4.0,4.0,4.0]],dtype=float)

        gaintx = numpy.zeros([8,8])
        gaintx[4:8,0:4] = 1
        
        gainrx = numpy.zeros([8,8])
        gainrx[0:4,4:8] = 1
        
        justrx = 0

    elif pattern == 173:
        # Meteors
        title=" Meteors : Tx Wup , Rx Sup "
        
        ues = numpy.array([0.0,0.0,0.0,0.0])*2/3
        
        phase = numpy.array([[5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0],
                             [5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0],
                             [5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0],
                             [5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0],
                             [4.0,4.0,4.0,4.0,4.0,4.0,4.0,4.0],
                             [4.0,4.0,4.0,4.0,4.0,4.0,4.0,4.0],
                             [4.0,4.0,4.0,4.0,4.0,4.0,4.0,4.0],
                             [4.0,4.0,4.0,4.0,4.0,4.0,4.0,4.0]],dtype=float)

        gaintx = numpy.zeros([8,8])
        gaintx[4:8,0:4] = 1
        
        gainrx = numpy.zeros([8,8])
        gainrx[4:8,4:8] = 1
        
        justrx = 0

    elif pattern == 174:
        # abs_test :
        
        title=" On-axis, tx: All up, rx: All up"
        
        ues = numpy.array([0,1.5,1.5,0])*2/3
        
        phase = numpy.array([[4,4,4,4,4,4,4,4],
                             [4,4,4,4,4,4,4,4],
                             [4,4,4,4,4,4,4,4],
                             [4,4,4,4,4,4,4,4],
                             [5,5,5,5,5,5,5,5],
                             [5,5,5,5,5,5,5,5],
                             [5,5,5,5,5,5,5,5],
                             [5,5,5,5,5,5,5,5]],dtype=float)

        gaintx = numpy.zeros([8,8])+1
        #gaintx[0:8,0:8] = 1
        #gaintx[0:4,0:4] = 1# N
        #gaintx[0:4,4:8] = 1# E
        #gaintx[4:8,0:4] = 1# W
        #gaintx[4:8,4:8] = 1# S
        
        gainrx = numpy.zeros([8,8])+1
        #gainrx[0:8,0:8] = 1
        #gainrx[0:4,0:4] = 1# N
        #gainrx[0:4,4:8] = 1# E
        #gainrx[4:8,0:4] = 1# W
        #gainrx[4:8,4:8] = 1# S
        
        justrx = 0
        
    elif pattern == 175:
        # abs_test :
        
        title=" On-axis, tx: All down, rx: All down"

        ues = numpy.array([1.5,3.0,3.0,1.5])*2/3
        ues = numpy.array([1.5,1.5,3.0,1.5])*2/3
        
        phase = numpy.array([[3,3,3,3,2,2,2,2],
                             [3,3,3,3,2,2,2,2],
                             [3,3,3,3,2,2,2,2],
                             [3,3,3,3,2,2,2,2],
                             [3,3,3,3,3,3,3,3],
                             [3,3,3,3,3,3,3,3],
                             [3,3,3,3,3,3,3,3],
                             [3,3,3,3,3,3,3,3]],dtype=float)

        gaintx = numpy.zeros([8,8])+1
        #gaintx[0:8,0:8] = 1
        #gaintx[0:4,0:4] = 1# N
        #gaintx[0:4,4:8] = 1# E
        #gaintx[4:8,0:4] = 1# W
        #gaintx[4:8,4:8] = 1# S
        
        gainrx = numpy.zeros([8,8])+1
        #gainrx[0:8,0:8] = 1
        #gainrx[0:4,0:4] = 1# N
        #gainrx[0:4,4:8] = 1# E
        #gainrx[4:8,0:4] = 1# W
        #gainrx[4:8,4:8] = 1# S
        
        justrx = 0

    elif pattern == 176:
        # ON AXIS
        title=" All_ON_AXIS_UP"
        
        ues = numpy.array([0.0,1.5,1.5,0.0])*2/3
        
        phase = numpy.array([[4,4,4,4,4,4,4,4],
                             [4,4,4,4,4,4,4,4],
                             [4,4,4,4,4,4,4,4],
                             [4,4,4,4,4,4,4,4],
                             [5,5,5,5,5,5,5,5],
                             [5,5,5,5,5,5,5,5],
                             [5,5,5,5,5,5,5,5],
                             [5,5,5,5,5,5,5,5]],dtype=float)

        gaintx = numpy.zeros([8,8])+1
        
        gainrx = numpy.zeros([8,8])+1
        
        justrx = 0   
        
    elif pattern == 177:
        # ON AXIS
        title=" All_ON_AXIS_DW"
        
        ues = numpy.array([1.5,1.5,3.0,0.0])*2/3
        
        phase = numpy.array([[3,3,3,3,2,2,2,2],
                             [3,3,3,3,2,2,2,2],
                             [3,3,3,3,2,2,2,2],
                             [3,3,3,3,2,2,2,2],
                             [3,3,3,3,4,4,4,4],
                             [3,3,3,3,4,4,4,4],
                             [3,3,3,3,4,4,4,4],
                             [3,3,3,3,4,4,4,4]],dtype=float)

        gaintx = numpy.zeros([8,8])+1
        
        gainrx = numpy.zeros([8,8])+1
        
        justrx = 0

    elif pattern == 178:
        # Meteors
        title=" Meteors : Tx Wup , Rx Nup "
        
        ues = numpy.array([0.0,1.5,0.0,1.5])*2/3
        
        phase = numpy.array([[3.0,3.0,4.0,4.0,3.0,3.0,4.0,4.0],
                             [3.0,4.0,4.0,5.0,3.0,4.0,4.0,5.0],
                             [4.0,4.0,5.0,5.0,4.0,4.0,5.0,5.0],
                             [4.0,5.0,5.0,2.0,4.0,5.0,5.0,2.0],
                             [2.0,2.0,3.0,3.0,4.0,4.0,5.0,5.0],
                             [2.0,3.0,3.0,4.0,4.0,5.0,5.0,2.0],
                             [3.0,3.0,4.0,4.0,5.0,5.0,2.0,2.0],
                             [3.0,4.0,4.0,5.0,5.0,2.0,2.0,3.0]],dtype=float)

        gaintx = numpy.zeros([8,8])
        gaintx[4:8,0:4] = 1
        
        gainrx = numpy.zeros([8,8])
        gainrx[0:4,0:4] = 1
        
        justrx = 0
    
    elif pattern == 179:
        # Meteors
        title=" Meteors : Tx Wup , Rx Eup "
        
        ues = numpy.array([0.0,1.5,0.0,1.5])*2/3
        
        phase = numpy.array([[3.0,3.0,4.0,4.0,3.0,3.0,4.0,4.0],
                             [3.0,4.0,4.0,5.0,3.0,4.0,4.0,5.0],
                             [4.0,4.0,5.0,5.0,4.0,4.0,5.0,5.0],
                             [4.0,5.0,5.0,2.0,4.0,5.0,5.0,2.0],
                             [2.0,2.0,3.0,3.0,4.0,4.0,5.0,5.0],
                             [2.0,3.0,3.0,4.0,4.0,5.0,5.0,2.0],
                             [3.0,3.0,4.0,4.0,5.0,5.0,2.0,2.0],
                             [3.0,4.0,4.0,5.0,5.0,2.0,2.0,3.0]],dtype=float)

        gaintx = numpy.zeros([8,8])
        gaintx[4:8,0:4] = 1
        
        gainrx = numpy.zeros([8,8])
        gainrx[0:4,4:8] = 1
        
        justrx = 0

    elif pattern == 180:
        # Meteors
        title=" Meteors : Tx Wup , Rx Sup "
        
        ues = numpy.array([0.0,1.5,0.0,1.5])*2/3
        
        phase = numpy.array([[3.0,3.0,4.0,4.0,3.0,3.0,4.0,4.0],
                             [3.0,4.0,4.0,5.0,3.0,4.0,4.0,5.0],
                             [4.0,4.0,5.0,5.0,4.0,4.0,5.0,5.0],
                             [4.0,5.0,5.0,2.0,4.0,5.0,5.0,2.0],
                             [2.0,2.0,3.0,3.0,4.0,4.0,5.0,5.0],
                             [2.0,3.0,3.0,4.0,4.0,5.0,5.0,2.0],
                             [3.0,3.0,4.0,4.0,5.0,5.0,2.0,2.0],
                             [3.0,4.0,4.0,5.0,5.0,2.0,2.0,3.0]],dtype=float)

        gaintx = numpy.zeros([8,8])
        gaintx[4:8,0:4] = 1
        
        gainrx = numpy.zeros([8,8])
        gainrx[4:8,4:8] = 1
        
        justrx = 0

    elif pattern == 181:
        # JULIA 15/07/16
        title=" for JULIA experiment" + '\n' + " NSdw transmission, Wdw reception "
        
        #ues = numpy.array([0.0,1.25 ,0.0,0.0])
        ues = numpy.array([0.0,0.75 ,0.0,0.0])
        #ues = numpy.array([0.0,0.0 ,0.0,0.0])
        
        phase = numpy.array([[4.0,4.0,3.0,3.0,2.0,2.0,5.0,5.0],
                             [5.0,4.0,4.0,3.0,3.0,2.0,2.0,5.0],
                             [5.0,5.0,4.0,4.0,3.0,3.0,2.0,2.0],
                             [2.0,5.0,5.0,4.0,4.0,3.0,3.0,2.0],
                             [2.0,2.0,5.0,5.0,4.0,4.0,3.0,3.0],
                             [3.0,2.0,2.0,5.0,5.0,4.0,4.0,3.0],
                             [3.0,3.0,2.0,2.0,5.0,5.0,4.0,4.0],
                             [4.0,3.0,3.0,2.0,2.0,5.0,5.0,4.0]],dtype=float)

        gaintx = numpy.zeros([8,8])
        gaintx[0:4,0:4] = 1
        gaintx[4:8,4:8] = 1
        
        gainrx = numpy.zeros([8,8])
        gainrx[4:8,0:4] = 1
        
        justrx = 0

    elif pattern == 182:
        # JULIA 15/07/16
        title=" for JULIA experiment" + '\n' + " NSdw transmission, Edw reception "
        
        ues = numpy.array([0.0,0.75 ,0.0,0.0])
        
        phase = numpy.array([[4.0,4.0,3.0,3.0,2.0,2.0,5.0,5.0],
                             [5.0,4.0,4.0,3.0,3.0,2.0,2.0,5.0],
                             [5.0,5.0,4.0,4.0,3.0,3.0,2.0,2.0],
                             [2.0,5.0,5.0,4.0,4.0,3.0,3.0,2.0],
                             [2.0,2.0,5.0,5.0,4.0,4.0,3.0,3.0],
                             [3.0,2.0,2.0,5.0,5.0,4.0,4.0,3.0],
                             [3.0,3.0,2.0,2.0,5.0,5.0,4.0,4.0],
                             [4.0,3.0,3.0,2.0,2.0,5.0,5.0,4.0]],dtype=float)

        gaintx = numpy.zeros([8,8])
        gaintx[0:4,0:4] = 1
        gaintx[4:8,4:8] = 1
        
        gainrx = numpy.zeros([8,8])
        gainrx[0:4,4:8] = 1
        
        justrx = 0

    elif pattern == 183:
        # JULIA 15/07/16
        title=" for JULIA experiment" + '\n' + " NSdw transmission, NSdw reception "
    
        ues = numpy.array([0.0,0.75 ,0.0,0.0])
        
        phase = numpy.array([[4.0,4.0,3.0,3.0,2.0,2.0,5.0,5.0],
                             [5.0,4.0,4.0,3.0,3.0,2.0,2.0,5.0],
                             [5.0,5.0,4.0,4.0,3.0,3.0,2.0,2.0],
                             [2.0,5.0,5.0,4.0,4.0,3.0,3.0,2.0],
                             [2.0,2.0,5.0,5.0,4.0,4.0,3.0,3.0],
                             [3.0,2.0,2.0,5.0,5.0,4.0,4.0,3.0],
                             [3.0,3.0,2.0,2.0,5.0,5.0,4.0,4.0],
                             [4.0,3.0,3.0,2.0,2.0,5.0,5.0,4.0]],dtype=float)

        gaintx = numpy.zeros([8,8])
        gaintx[0:4,0:4] = 1
        gaintx[4:8,4:8] = 1
        
        gainrx = numpy.zeros([8,8])
        gainrx[0:4,0:4] = 1
        gainrx[4:8,4:8] = 1
        
        justrx = 0
    
    elif pattern == 184:
        # JULIA 15/07/16
        title=" for JULIA experiment" + '\n' + " NSup transmission, NSup reception "
    
        ues = numpy.array([0.0,0.0,0.0,1.5])
        
        phase = numpy.array([[4.0,4.0,5.0,5.0,3.0,3.0,4.0,4.0],
                             [3.0,3.0,4.0,4.0,3.0,4.0,4.0,5.0],
                             [2.0,2.0,3.0,3.0,4.0,4.0,5.0,5.0],
                             [5.0,5.0,2.0,2.0,4.0,5.0,5.0,2.0],
                             [2.0,2.0,3.0,3.0,5.0,5.0,4.0,2.0],
                             [5.0,2.0,2.0,3.0,4.0,4.0,5.0,5.0],
                             [5.0,5.0,2.0,2.0,3.0,3.0,4.0,4.0],
                             [2.0,5.0,5.0,2.0,2.0,2.0,3.0,3.0]],dtype=float)
        gaintx = numpy.zeros([8,8])
        gaintx[0:4,0:4] = 1
        gaintx[4:8,4:8] = 1
        gaintx[4,6] = 0
        
        gainrx = numpy.zeros([8,8])
        gainrx[0:4,0:4] = 1
        gainrx[4:8,4:8] = 1
        
        justrx = 0

#     elif pattern == 185:
#         # JULIA 15/07/16
#         title=" Imaging: CHA "
#         
#         ues = numpy.array([0.0,0.0,0.0,0.0])*2/3
#         #ues = numpy.array([1.5,3.0,3.0,1.5])*2/3
#         #ues = numpy.array([0.0,1.5,1.5,0.0])*2/3
#         
#         phase = numpy.array([[4.0,4.0,3.0,3.0,4.0,3.0,4.0,4.0],
#                              [5.0,4.0,4.0,3.0,3.0,3.0,4.0,4.0],
#                              [5.0,5.0,4.0,4.0,4.0,3.0,3.0,2.0],
#                              [2.0,5.0,5.0,4.0,5.0,4.0,2.0,2.0],
#                              [2.0,2.0,5.0,5.0,4.0,4.0,3.0,3.0],
#                              [3.0,2.0,2.0,5.0,5.0,4.0,4.0,3.0],
#                              [3.0,3.0,2.0,2.0,5.0,5.0,4.0,4.0],
#                              [4.0,3.0,3.0,2.0,2.0,5.0,5.0,4.0]],dtype=float)
# 
#         gaintx = numpy.zeros([8,8])+1
#         gaintx[0,7] = 0
#         gaintx[1,7] = 0
#         gaintx[1,6] = 0
#         gaintx[2,4] = 0
#         
#         gainrx = numpy.zeros([8,8])
#         #gainrx[4:8,0:4] = 1
#         gainrx[0,7] = 1
#         
#         justrx = 0
        
    elif pattern == 186:
        # EW_DRIFT 22/07/16
        title=" Tx Edw + Wdw, Rx Edw "
        
        ues = numpy.array([3.0,0.0,3.0,0.0])*2/3
        
        phase = numpy.array([[4.0,4.0,3.0,3.0,4.0,4.0,3.0,3.0],
                             [5.0,4.0,4.0,3.0,5.0,4.0,4.0,3.0],
                             [5.0,5.0,4.0,4.0,5.0,5.0,4.0,4.0],
                             [2.0,5.0,5.0,4.0,2.0,5.0,5.0,4.0],
                             [4.0,4.0,3.0,3.0,4.0,4.0,3.0,3.0],
                             [5.0,4.0,4.0,3.0,5.0,4.0,4.0,3.0],
                             [5.0,5.0,4.0,4.0,5.0,5.0,4.0,4.0],
                             [2.0,5.0,5.0,4.0,2.0,5.0,5.0,4.0]],dtype=float)

        gaintx = numpy.zeros([8,8])
        gaintx[0:4,4:8] = 1
        gaintx[4:8,0:4] = 1
        
        gainrx = numpy.zeros([8,8])
        gainrx[0:4,4:8] = 1
        
        justrx = 0

    elif pattern == 187:
        # EW_DRIFT 22/07/16
        title=" Tx Edw + Wdw, Rx Wdw "
        
        ues = numpy.array([3.0,0.0,3.0,0.0])*2/3
        
        phase = numpy.array([[4.0,4.0,3.0,3.0,4.0,4.0,3.0,3.0],
                             [5.0,4.0,4.0,3.0,5.0,4.0,4.0,3.0],
                             [5.0,5.0,4.0,4.0,5.0,5.0,4.0,4.0],
                             [2.0,5.0,5.0,4.0,2.0,5.0,5.0,4.0],
                             [4.0,4.0,3.0,3.0,4.0,4.0,3.0,3.0],
                             [5.0,4.0,4.0,3.0,5.0,4.0,4.0,3.0],
                             [5.0,5.0,4.0,4.0,5.0,5.0,4.0,4.0],
                             [2.0,5.0,5.0,4.0,2.0,5.0,5.0,4.0]],dtype=float)

        gaintx = numpy.zeros([8,8])
        gaintx[0:4,4:8] = 1
        gaintx[4:8,0:4] = 1
        
        gainrx = numpy.zeros([8,8])
        gainrx[4:8,0:4] = 1
        
        justrx = 0

    elif pattern == 188:
        # EW_DRIFT 22/07/16
        title=" Tx  Eup + Wup, Rx Wup "
        
        ues = numpy.array([5.0,0.0,3.0,0.0])*2/3

        phase = numpy.array([[4.0,4.0,5.0,5.0,4.0,4.0,5.0,5.0],
                             [3.0,3.0,4.0,4.0,3.0,3.0,4.0,4.0],
                             [2.0,2.0,3.0,3.0,2.0,2.0,3.0,3.0],
                             [5.0,5.0,2.0,2.0,5.0,5.0,2.0,2.0],
                             [5.0,5.0,2.0,2.0,5.0,5.0,4.0,2.0],
                             [4.0,4.0,5.0,5.0,4.0,4.0,5.0,5.0],
                             [3.0,3.0,4.0,4.0,3.0,3.0,4.0,4.0],
                             [2.0,2.0,3.0,3.0,2.0,2.0,3.0,3.0]],dtype=float)

        gaintx = numpy.zeros([8,8])
        gaintx[0:4,4:8] = 1
        gaintx[4:8,0:4] = 1
        
        gainrx = numpy.zeros([8,8])
        gainrx[4:8,0:4] = 1
        
        justrx = 0
        
    elif pattern == 189:
        # EW_DRIFT 22/07/16
        title=" Tx Eup + Wup, Rx Eup "
        
        ues = numpy.array([5.0,0.0,3.0,0.0])*2/3
        
        phase = numpy.array([[4.0,4.0,5.0,5.0,4.0,4.0,5.0,5.0],
                             [3.0,3.0,4.0,4.0,3.0,3.0,4.0,4.0],
                             [2.0,2.0,3.0,3.0,2.0,2.0,3.0,3.0],
                             [5.0,5.0,2.0,2.0,5.0,5.0,2.0,2.0],
                             [5.0,5.0,2.0,2.0,5.0,5.0,4.0,2.0],
                             [4.0,4.0,5.0,5.0,4.0,4.0,5.0,5.0],
                             [3.0,3.0,4.0,4.0,3.0,3.0,4.0,4.0],
                             [2.0,2.0,3.0,3.0,2.0,2.0,3.0,3.0]],dtype=float)

        gaintx = numpy.zeros([8,8])
        gaintx[0:4,4:8] = 1
        gaintx[4:8,0:4] = 1
        
        gainrx = numpy.zeros([8,8])
        gainrx[0:4,4:8] = 1
        
        justrx = 0

    elif pattern == 190:
        # ON AXIS
        title=" wUP "
        
        ues = numpy.array([0,0,0,0])*2/3
#        ues = numpy.array([0,1.5,1.5,0])*2/3
        
        phase = numpy.array([[4,4,4,4,4,4,4,4],
                             [4,4,4,4,4,4,4,4],
                             [4,4,4,4,4,4,4,4],
                             [4,4,4,4,4,4,4,4],
                             [5,5,5,5,5,5,5,5],
                             [5,5,5,5,5,5,5,5],
                             [5,5,5,5,5,5,5,5],
                             [5,5,5,5,5,5,5,5]],dtype=float)

        gaintx = numpy.zeros([8,8])
        gaintx[0:4,0:4] = 1
        #gaintx[4:8,4:8] = 1
        
        gainrx = numpy.zeros([8,8])
        gainrx[4:8,4:8] = 1
        
        justrx = 0
               

    elif pattern == 191:
        # ON AXIS
        title=" wUP "
        
        ues = numpy.array([0,0,0,0])*2/3
#        ues = numpy.array([0,1.5,1.5,0])*2/3
        
        phase = numpy.array([[4,4,4,4,4,4,4,4],
                             [4,4,4,4,4,4,4,4],
                             [4,4,4,4,4,4,4,4],
                             [4,4,4,4,4,4,4,4],
                             [5,5,5,5,5,5,5,5],
                             [5,5,5,5,5,5,5,5],
                             [5,5,5,5,5,5,5,5],
                             [5,5,5,5,5,5,5,5]],dtype=float)

        gaintx = numpy.zeros([8,8])
        gaintx[0:4,0:4] = 1
        #gaintx[4:8,4:8] = 1
        
        gainrx = numpy.zeros([8,8])
        gainrx[4:8,4:8] = 1
        
        justrx = 0
        
    elif pattern == 192:
        # LP FARADAY
        title=" chB : (4.5) all DW"
        
        ues = numpy.array([0,1,1,0])
        
        phase = numpy.array([[2,2,3,3,4,4,5,5],
                             [2,3,3,4,4,5,5,2],
                             [3,3,4,4,5,5,2,2],
                             [3,4,4,5,5,2,2,3],
                             [5,5,2,2,3,3,4,4],
                             [5,2,2,3,3,4,4,5],
                             [2,2,3,3,4,4,5,5],
                             [2,3,3,4,4,5,5,2]],dtype=float)

        gaintx = numpy.zeros([8,8])
        gaintx[0:8,0:8] = 1
        
        gainrx = numpy.zeros([8,8])
        gainrx[0:8,0:8] = 1
        
        justrx = 0

    elif pattern == 193:
        # LP FARADAY
        title=" JULIA 2008 dw"
        
        ues = numpy.array([1,0,0,1.4])
        ues = numpy.array([0,0,0,0])
        
        phase = numpy.array([[2,3,3,3,3,4,4,4],
                             [5,2,2,2,2,3,3,3],
                             [3,4,4,4,4,5,5,5],
                             [2,3,3,3,3,4,4,4],
                             [4,5,5,5,5,2,2,2],
                             [3,4,4,4,4,5,5,5],
                             [5,2,2,2,2,3,3,3],
                             [4,5,5,5,5,2,2,2]],dtype=float)

        gaintx = numpy.zeros([8,8])+1
        
        gainrx = numpy.zeros([8,8])+1
        
        justrx = 0        
    
    elif pattern == 194:
        # LP FARADAY
        title=" JULIA 2008 up"
        
        ues = numpy.array([0,0.4,0,0])
        
        phase = numpy.array([[2,5,3,2,5,4,2,5],
                             [3,2,4,3,2,5,3,2],
                             [3,2,4,3,2,5,3,2],
                             [4,3,5,4,3,2,4,3],
                             [4,3,5,4,3,2,4,3],
                             [5,4,2,5,4,3,5,4],
                             [5,4,2,5,4,3,5,4],
                             [2,5,3,2,5,4,2,5]],dtype=float)

        gaintx = numpy.zeros([8,8])+1
        
        gainrx = numpy.zeros([8,8])
        gainrx[0:4,0:4] = 1
        gainrx[4:8,4:8] = 1
        
        justrx = 0    

    elif pattern == 195:
        # JULIA, octubre 2009
        title=" JULIA D: n_s_UP"
        
        ues = numpy.array([0,0,0,0.25])
        
        phase = numpy.array([[4,4,4,4,4,3,4,4],
                             [3,3,3,3,3,3,4,4],
                             [2,2,2,2,4,3,3,2],
                             [5,5,5,5,5,4,2,2],
                             [4,4,4,4,4,4,4,4],
                             [4,4,4,4,3,3,3,3],
                             [4,4,4,4,2,2,2,2],
                             [4,4,4,4,5,5,5,5]],dtype=float)

        gaintx = numpy.zeros([8,8])
        gaintx[0:4,0:4] = 1
        gaintx[4:8,4:8] = 1
        gaintx[4,6] = 0
        
        gainrx = numpy.zeros([8,8])
        gainrx[0:4,0:4] = 1
        gainrx[4:8,4:8] = 1
        gainrx[4,6] = 0
        
        justrx = 0

    elif pattern == 196:
        # JULIA, octubre 2009
        title=" JULIA CHC : n_s_DW"
        
        ues = numpy.array([0,0,0,0.75])
        
        phase = numpy.array([[5,4,2,5,5,4,3,2],
                             [2,5,3,2,2,5,4,3],
                             [2,5,3,2,2,5,4,3],
                             [3,2,4,3,3,2,5,4],
                             [5,4,3,2,5,4,2,5],
                             [2,5,4,3,2,5,3,2],
                             [2,5,4,3,2,5,3,2],
                             [3,2,5,4,3,2,4,3]],dtype=float)

        gaintx = numpy.zeros([8,8])+1
        #gaintx[0:4,0:4] = 1
        #gaintx[4:8,4:8] = 1
        
        gainrx = numpy.zeros([8,8])+1
        #gainrx[0:4,0:4] = 1
        #gainrx[4:8,4:8] = 1
        justrx = 0
        
    elif pattern == 197:
        # JULIA, 26 junio 2008
        title="A: tx = N&S_up, rx = N&S_up"
        
        ues = numpy.array([0.4,0,1,1.4])*2/3
        
        phase = numpy.array([[5,4,2,5,3,2,4,3],
                             [2,5,3,2,4,3,5,4],
                             [2,5,3,2,4,3,5,4],
                             [3,2,4,3,5,4,2,5],
                             [3,2,4,3,5,4,2,5],
                             [4,3,5,4,2,5,3,2],
                             [4,3,5,4,2,5,3,2],
                             [5,4,2,5,3,2,4,3]],dtype=float)

        gaintx = numpy.zeros([8,8])
        gaintx[0:4,0:4] = 1
        gaintx[4:8,4:8] = 1
        
        gainrx = numpy.zeros([8,8])
        gainrx[0:4,0:4] = 1
        gainrx[4:8,4:8] = 1
        #gainrx[0:4,4:8] = 1
        #gainrx[4:8,0:4] = 1
        
        justrx = 0

    elif pattern == 198:
        # JULIA, 26 junio 2008
        title="C : tx = all_dw, rx = all_dw"
        
        ues = numpy.array([0,0,0,0.0])
        
        phase = numpy.array([[2,3,3,3,3,4,4,4],
                             [5,2,2,2,2,3,3,3],
                             [3,4,4,4,4,5,5,5],
                             [2,3,3,3,3,4,4,4],
                             [4,5,5,5,5,2,2,2],
                             [3,4,4,4,4,5,5,5],
                             [5,2,2,2,2,3,3,3],
                             [4,5,5,5,5,2,2,2]],dtype=float)

        gaintx = numpy.zeros([8,8])+1
        
        gainrx = numpy.zeros([8,8])+1

        justrx = 0
    
    elif pattern == 199:
        # JULIA, junio 2013
        title=" JULIA D: n_s_UP"
        
        ues = numpy.array([0,0,0,0.25])
        
        phase = numpy.array([[4,4,5,5,4,3,4,4],
                             [3,3,4,4,3,3,4,4],
                             [2,2,3,3,4,3,3,2],
                             [5,5,2,2,5,4,2,2],
                             [4,4,4,4,5,5,4,2],
                             [4,4,4,4,4,4,5,5],
                             [4,4,4,4,3,3,4,4],
                             [4,4,4,4,2,2,3,3]],dtype=float)

        gaintx = numpy.zeros([8,8])
        gaintx[0:4,0:4] = 1
        gaintx[4:8,4:8] = 1
        gaintx[4,6] = 0
        
        gainrx = numpy.zeros([8,8])
        gainrx[0:4,0:4] = 1
        gainrx[4:8,4:8] = 1
        gainrx[4,6] = 0
        
        justrx = 0
    
    elif pattern == 200:
        # JULIA, junio 2013
        title=" JULIA CHC : All_DW"
        
        ues = numpy.array([0,0,0,0.75])
        
        phase = numpy.array([[5,4,3,2,5,4,3,2],
                             [2,5,4,3,2,5,4,3],
                             [2,5,4,3,2,5,4,3],
                             [3,2,5,4,3,2,5,4],
                             [5,4,3,2,3,2,5,4],
                             [2,5,4,3,4,3,2,5],
                             [2,5,4,3,4,3,2,5],
                             [3,2,5,4,5,4,3,2]],dtype=float)

        gaintx = numpy.zeros([8,8])
        gaintx[0:4,0:4] = 1
        gaintx[4:8,4:8] = 1
                
        gainrx = numpy.zeros([8,8])
        gainrx[0:4,0:4] = 1
        gainrx[4:8,4:8] = 1
        
        justrx = 0

        
    elif pattern == 201:
        # EW_DRIFT 05/10/16
        title=" Tx Eup + Wup, Rx Eup (CHB) "
        
        ues = numpy.array([3.0,0.0,3.0,0.0])*2/3
        
        phase = numpy.array([[4.0,4.0,3.0,3.0,4.0,4.0,3.0,3.0],
                             [5.0,4.0,4.0,3.0,5.0,4.0,4.0,3.0],
                             [5.0,5.0,4.0,4.0,5.0,5.0,4.0,4.0],
                             [2.0,5.0,5.0,4.0,2.0,5.0,5.0,4.0],
                             [4.0,4.0,3.0,3.0,4.0,4.0,3.0,3.0],
                             [5.0,4.0,4.0,3.0,5.0,4.0,4.0,3.0],
                             [5.0,5.0,4.0,4.0,5.0,5.0,4.0,4.0],
                             [2.0,5.0,5.0,4.0,2.0,5.0,5.0,4.0]],dtype=float)

        gaintx = numpy.zeros([8,8])
        gaintx[0:4,4:8] = 1
        gaintx[4:8,0:4] = 1
        
        gainrx = numpy.zeros([8,8])
        gainrx[0:4,4:8] = 1
        
        justrx = 0

    elif pattern == 202:
        # EW_DRIFT 05/10/16
        title=" Tx Eup + Wup, Rx Wup (CHA)"
        
        ues = numpy.array([3.0,0.0,3.0,0.0])*2/3
        
        phase = numpy.array([[4.0,4.0,3.0,3.0,4.0,4.0,3.0,3.0],
                             [5.0,4.0,4.0,3.0,5.0,4.0,4.0,3.0],
                             [5.0,5.0,4.0,4.0,5.0,5.0,4.0,4.0],
                             [2.0,5.0,5.0,4.0,2.0,5.0,5.0,4.0],
                             [4.0,4.0,3.0,3.0,4.0,4.0,3.0,3.0],
                             [5.0,4.0,4.0,3.0,5.0,4.0,4.0,3.0],
                             [5.0,5.0,4.0,4.0,5.0,5.0,4.0,4.0],
                             [2.0,5.0,5.0,4.0,2.0,5.0,5.0,4.0]],dtype=float)

        gaintx = numpy.zeros([8,8])
        gaintx[0:4,4:8] = 1
        gaintx[4:8,0:4] = 1
        
        gainrx = numpy.zeros([8,8])
        gainrx[4:8,0:4] = 1
        
        justrx = 0

    elif pattern == 203:
        # EW_DRIFT 05/10/16
        title=" Tx  Edw + Wdw, Rx Wdw (CHC)"
        
        ues = numpy.array([5.0,0.0,3.0,0.0])*2/3

        phase = numpy.array([[4.0,4.0,5.0,5.0,4.0,4.0,5.0,5.0],
                             [3.0,3.0,4.0,4.0,3.0,3.0,4.0,4.0],
                             [2.0,2.0,3.0,3.0,2.0,2.0,3.0,3.0],
                             [5.0,5.0,2.0,2.0,5.0,5.0,2.0,2.0],
                             [5.0,5.0,2.0,2.0,5.0,5.0,4.0,2.0],
                             [4.0,4.0,5.0,5.0,4.0,4.0,5.0,5.0],
                             [3.0,3.0,4.0,4.0,3.0,3.0,4.0,4.0],
                             [2.0,2.0,3.0,3.0,2.0,2.0,3.0,3.0]],dtype=float)

        gaintx = numpy.zeros([8,8])
        gaintx[0:4,4:8] = 1
        gaintx[4:8,0:4] = 1
        
        gainrx = numpy.zeros([8,8])
        gainrx[4:8,0:4] = 1
        
        justrx = 0
        
    elif pattern == 204:
        # EW_DRIFT 05/10/16
        title=" Tx Edw + Wdw, Rx Edw (CHD)"
        
        ues = numpy.array([5.0,0.0,3.0,0.0])*2/3
        
        phase = numpy.array([[4.0,4.0,5.0,5.0,4.0,4.0,5.0,5.0],
                             [3.0,3.0,4.0,4.0,3.0,3.0,4.0,4.0],
                             [2.0,2.0,3.0,3.0,2.0,2.0,3.0,3.0],
                             [5.0,5.0,2.0,2.0,5.0,5.0,2.0,2.0],
                             [5.0,5.0,2.0,2.0,5.0,5.0,4.0,2.0],
                             [4.0,4.0,5.0,5.0,4.0,4.0,5.0,5.0],
                             [3.0,3.0,4.0,4.0,3.0,3.0,4.0,4.0],
                             [2.0,2.0,3.0,3.0,2.0,2.0,3.0,3.0]],dtype=float)

        gaintx = numpy.zeros([8,8])
        gaintx[0:4,4:8] = 1
        gaintx[4:8,0:4] = 1
        
        gainrx = numpy.zeros([8,8])
        gainrx[0:4,4:8] = 1
        
        justrx = 0

    elif pattern == 205:
        # ON AXIS
        title="ON_AXIS"
        
        ues = numpy.array([0,0,0,0])*2/3
        
        phase = numpy.array([[4,4,4,4,4,4,4,4],
                             [4,4,4,4,4,4,4,4],
                             [4,4,4,4,4,4,4,4],
                             [4,4,4,4,4,4,4,4],
                             [5,5,5,5,5,5,5,5],
                             [5,5,5,5,5,5,5,5],
                             [5,5,5,5,5,5,5,5],
                             [5,5,5,5,5,5,5,5]],dtype=float)

        gaintx = numpy.zeros([8,8])+1
        
        gainrx = numpy.zeros([8,8])
        gainrx[0:4,0:4] = 1
        
        justrx = 1
####
    elif pattern == 206:
        # EW_DRIFT 28/11/16
        title=" Tx Edw + Wdw, Rx Edw "
        
        ues = numpy.array([3.0,4.5,3.0,3.0])*2/3
        
        phase = numpy.array([[2.0,2.0,3.0,3.0,4.0,4.0,3.0,3.0],
                             [2.0,3.0,3.0,4.0,5.0,4.0,4.0,3.0],
                             [3.0,3.0,4.0,4.0,5.0,5.0,4.0,4.0],
                             [3.0,4.0,4.0,5.0,2.0,5.0,5.0,4.0],
                             [4.0,4.0,3.0,3.0,3.0,3.0,4.0,4.0],
                             [5.0,4.0,4.0,3.0,3.0,4.0,4.0,5.0],
                             [5.0,5.0,4.0,4.0,4.0,4.0,5.0,5.0],
                             [2.0,5.0,5.0,4.0,4.0,5.0,5.0,2.0]],dtype=float)

        gaintx = numpy.zeros([8,8])
        gaintx[0:4,4:8] = 1
        gaintx[4:8,0:4] = 1
        
        gainrx = numpy.zeros([8,8])
        gainrx[0:4,4:8] = 1
        
        justrx = 0

    elif pattern == 207:
        # EW_DRIFT 28/11/16
        title=" Tx Edw + Wdw, Rx Wdw "
        
        ues = numpy.array([3.0,4.5,3.0,3.0])*2/3
        
        phase = numpy.array([[2.0,2.0,3.0,3.0,4.0,4.0,3.0,3.0],
                             [2.0,3.0,3.0,4.0,5.0,4.0,4.0,3.0],
                             [3.0,3.0,4.0,4.0,5.0,5.0,4.0,4.0],
                             [3.0,4.0,4.0,5.0,2.0,5.0,5.0,4.0],
                             [4.0,4.0,3.0,3.0,3.0,3.0,4.0,4.0],
                             [5.0,4.0,4.0,3.0,3.0,4.0,4.0,5.0],
                             [5.0,5.0,4.0,4.0,4.0,4.0,5.0,5.0],
                             [2.0,5.0,5.0,4.0,4.0,5.0,5.0,2.0]],dtype=float)

        gaintx = numpy.zeros([8,8])
        gaintx[0:4,4:8] = 1
        gaintx[4:8,0:4] = 1
        
        gainrx = numpy.zeros([8,8])
        gainrx[4:8,0:4] = 1
        
        justrx = 0

    elif pattern == 208:
        # EW_DRIFT 28/11/16
        title=" Tx  Eup + Wup, Rx Wup "
        
        ues = numpy.array([5.0,0.0,3.0,4.5])*2/3

        phase = numpy.array([[4.0,4.0,5.0,5.0,4.0,4.0,5.0,5.0],
                             [4.0,5.0,5.0,2.0,3.0,3.0,4.0,4.0],
                             [5.0,5.0,2.0,2.0,2.0,2.0,3.0,3.0],
                             [5.0,2.0,2.0,3.0,5.0,5.0,2.0,2.0],
                             [5.0,5.0,2.0,2.0,5.0,5.0,2.0,2.0],
                             [4.0,4.0,5.0,5.0,5.0,2.0,2.0,3.0],
                             [3.0,3.0,4.0,4.0,2.0,2.0,3.0,3.0],
                             [2.0,2.0,3.0,3.0,2.0,3.0,3.0,4.0]],dtype=float)

        gaintx = numpy.zeros([8,8])
        gaintx[0:4,4:8] = 1
        gaintx[4:8,0:4] = 1
        
        gainrx = numpy.zeros([8,8])
        gainrx[4:8,0:4] = 1
        
        justrx = 0
        
    elif pattern == 209:
        # EW_DRIFT 28/11/16
        title=" Tx Eup + Wup, Rx Eup "
        
        ues = numpy.array([5.0,0.0,3.0,4.5])*2/3

        phase = numpy.array([[4.0,4.0,5.0,5.0,4.0,4.0,5.0,5.0],
                             [4.0,5.0,5.0,2.0,3.0,3.0,4.0,4.0],
                             [5.0,5.0,2.0,2.0,2.0,2.0,3.0,3.0],
                             [5.0,2.0,2.0,3.0,5.0,5.0,2.0,2.0],
                             [5.0,5.0,2.0,2.0,5.0,5.0,2.0,2.0],
                             [4.0,4.0,5.0,5.0,5.0,2.0,2.0,3.0],
                             [3.0,3.0,4.0,4.0,2.0,2.0,3.0,3.0],
                             [2.0,2.0,3.0,3.0,2.0,3.0,3.0,4.0]],dtype=float)

        gaintx = numpy.zeros([8,8])
        gaintx[0:4,4:8] = 1
        gaintx[4:8,0:4] = 1
        
        gainrx = numpy.zeros([8,8])
        gainrx[0:4,4:8] = 1
        
        justrx = 0
        
    elif pattern == 210:
        # FARADAY 28/11/16
        title=" Tx Nup + Sup, Rx Nup + Sup "
        
        ues = numpy.array([5.0,0.0,3.0,4.5])*2/3

        phase = numpy.array([[4.0,4.0,5.0,5.0,4.0,4.0,5.0,5.0],
                             [4.0,5.0,5.0,2.0,3.0,3.0,4.0,4.0],
                             [5.0,5.0,2.0,2.0,2.0,2.0,3.0,3.0],
                             [5.0,2.0,2.0,3.0,5.0,5.0,2.0,2.0],
                             [5.0,5.0,2.0,2.0,5.0,5.0,2.0,2.0],
                             [4.0,4.0,5.0,5.0,5.0,2.0,2.0,3.0],
                             [3.0,3.0,4.0,4.0,2.0,2.0,3.0,3.0],
                             [2.0,2.0,3.0,3.0,2.0,3.0,3.0,4.0]],dtype=float)

        gaintx = numpy.zeros([8,8])
        gaintx[0:4,0:4] = 1
        gaintx[4:8,4:8] = 1
        
        gainrx = numpy.zeros([8,8])
        gainrx[0:4,0:4] = 1
        gainrx[4:8,4:8] = 1
        
        justrx = 0

    elif pattern == 211:
        # FARADAY 28/11/16
        title=" Tx Ndw + Sdw, Rx Ndw + Sdw "
        
        ues = numpy.array([3.0,4.5,3.0,3.0])*2/3
        
        phase = numpy.array([[2.0,2.0,3.0,3.0,4.0,4.0,3.0,3.0],
                             [2.0,3.0,3.0,4.0,5.0,4.0,4.0,3.0],
                             [3.0,3.0,4.0,4.0,5.0,5.0,4.0,4.0],
                             [3.0,4.0,4.0,5.0,2.0,5.0,5.0,4.0],
                             [4.0,4.0,3.0,3.0,3.0,3.0,4.0,4.0],
                             [5.0,4.0,4.0,3.0,3.0,4.0,4.0,5.0],
                             [5.0,5.0,4.0,4.0,4.0,4.0,5.0,5.0],
                             [2.0,5.0,5.0,4.0,4.0,5.0,5.0,2.0]],dtype=float)

        gaintx = numpy.zeros([8,8])
        gaintx[0:4,0:4] = 1
        gaintx[4:8,4:8] = 1
        
        gainrx = numpy.zeros([8,8])
        gainrx[0:4,0:4] = 1
        gainrx[4:8,4:8] = 1
        
        justrx = 0

    elif pattern == 212:
        # ON AXIS
        title="UP_ON_AXIS"
        
        ues = numpy.array([0.0,0.0,0.0,0.0])*2/3
        
        phase = numpy.array([[4,4,4,4,4,3,2,5],
                             [4,4,4,4,5,4,3,2],
                             [4,4,4,4,2,5,4,3],
                             [4,4,4,4,3,2,5,4],
                             [5,5,5,5,5,4,3,2],
                             [5,5,5,5,2,5,4,3],
                             [5,5,5,5,3,2,5,4],
                             [5,5,5,5,4,3,2,5]],dtype=float)

        gaintx = numpy.zeros([8,8])
        gaintx[0:8,0:4] = 1
        
        gainrx = numpy.zeros([8,8])+1
        gainrx[0:4,0:4] = 1
              
        justrx = 0   
        
    elif pattern == 213:
        # 3.38 west magnetic
        title=" DW_3.38_OFF_AXIS"
        
        ues = numpy.array([0.0,0.0,0.0,0.0])*2/3
        
        phase = numpy.array([[2,2,2,2,2,5,4,3],
                             [2,2,2,2,3,2,5,4],
                             [2,2,2,2,4,3,2,5],
                             [2,2,2,2,5,4,3,2],
                             [3,3,3,3,3,2,5,4],
                             [3,3,3,3,4,3,2,5],
                             [3,3,3,3,5,4,3,2],
                             [3,3,3,3,2,5,4,3]],dtype=float)

        gaintx = numpy.zeros([8,8])
        gaintx[0:8,4:8] = 1
        
        gainrx = numpy.zeros([8,8])
        gainrx[0:4,4:8] = 1
                
        justrx = 1
        
    elif pattern == 214:
        # Pase de sol 2017_02
        title=" ES_TX_RX"
        
#         ues = numpy.array([1.5,0.0,0.0,0])*2/3
#         ues = numpy.array([1.5,0.0,7*(1.5/8.0),0.0])*2/3
#         ues = numpy.array([1.5,0.0,2.25,0.0])*2/3
#         ues = numpy.array([1.5,0.0,1.3,0.0])*2/3
        ues = numpy.array([0.0,0.0,0.0,0.0])*2/3
#         ues = numpy.array([0.0,0.0,1.3,0.0])*2/3
        
        
        phase = numpy.array([[2,5,5,4,4,4,5,5],
                             [3,2,5,5,3,4,4,5],
                             [4,3,2,5,3,3,4,4],
                             [4,4,3,2,2,3,3,4],
                             [4,3,3,2,2,2,3,3],
                             [5,4,3,3,5,2,2,3],
                             [2,5,4,3,5,5,2,2],
                             [2,2,5,4,4,5,5,2]],dtype=float)

        gaintx = numpy.zeros([8,8])
        gaintx[0:8,4:8] = 1
        
        gainrx = numpy.zeros([8,8])
        gainrx[0:8,4:8] = 1
                
        justrx = 0
        
    elif pattern == 215:
        # Pase de sol 2017_02
        title=" NW_JUST_RX"
        
#         ues = numpy.array([1.5,0.0,1.3,0.0])*2/3
#         ues = numpy.array([1.5,0.0,2.25,0.0])*2/3
#         ues = numpy.array([0.0,0.0,0.0,0.0])*2/3
#         ues = numpy.array([0.0,0.0,1.3,0.0])*2/3
        ues = numpy.array([0.0,0.0,0.0,0.0])*2/3
        
        phase = numpy.array([[2,5,5,4,4,4,5,5],
                             [3,2,5,5,3,4,4,5],
                             [4,3,2,5,3,3,4,4],
                             [4,4,3,2,2,3,3,4],
                             [4,3,3,2,2,2,3,3],
                             [5,4,3,3,5,2,2,3],
                             [2,5,4,3,5,5,2,2],
                             [2,2,5,4,4,5,5,2]],dtype=float)

        
        gaintx = numpy.zeros([8,8])
        gaintx[0:8,0:4] = 1
        
        gainrx = numpy.zeros([8,8])
        gainrx[0:8,0:4] = 1
                
        justrx = 1
 
    elif pattern == 216:
        # JULIA 28/02/17
        title=" Tx NSdw , Rx NSdw "
    
#         ues = numpy.array([0.0,0.75 ,0.0,0.0])
#         
#         phase = numpy.array([[4.0,4.0,3.0,3.0,2.0,2.0,5.0,5.0],
#                              [5.0,4.0,4.0,3.0,3.0,2.0,2.0,5.0],
#                              [5.0,5.0,4.0,4.0,3.0,3.0,2.0,2.0],
#                              [2.0,5.0,5.0,4.0,4.0,3.0,3.0,2.0],
#                              [2.0,2.0,5.0,5.0,4.0,4.0,3.0,3.0],
#                              [3.0,2.0,2.0,5.0,5.0,4.0,4.0,3.0],
#                              [3.0,3.0,2.0,2.0,5.0,5.0,4.0,4.0],
#                              [4.0,3.0,3.0,2.0,2.0,5.0,5.0,4.0]],dtype=float)
        ues = numpy.array([0.0,0.0,0.0,1.5])
        ues = numpy.array([0.0,0.0,0.0,0.75])
        ues = numpy.array([0.0,0.0,0.0,1.25])
        ues = numpy.array([0.0,0.0,0.0,1.75])
        
        phase = numpy.array([[4.0,4.0,5.0,5.0,3.0,3.0,4.0,4.0],
                             [3.0,3.0,4.0,4.0,3.0,4.0,4.0,5.0],
                             [2.0,2.0,3.0,3.0,4.0,4.0,5.0,5.0],
                             [5.0,5.0,2.0,2.0,4.0,5.0,5.0,2.0],
                             [2.0,2.0,3.0,3.0,5.0,5.0,4.0,2.0],
                             [5.0,2.0,2.0,3.0,4.0,4.0,5.0,5.0],
                             [5.0,5.0,2.0,2.0,3.0,3.0,4.0,4.0],
                             [2.0,5.0,5.0,2.0,2.0,2.0,3.0,3.0]],dtype=float)

        gaintx = numpy.zeros([8,8])
        gaintx[0:4,0:4] = 1
        gaintx[4:8,4:8] = 1
        gaintx[4,6] = 0
        
        gainrx = numpy.zeros([8,8])
        gainrx[0:4,0:4] = 1
        gainrx[4:8,4:8] = 1
        gainrx[4,6] = 0
        
        justrx = 0
    
    elif pattern == 217:
        # JULIA 28/02/17
        title=" Tx NSup , Rx NSup "
    
        ues = numpy.array([0.0,0.0,0.0,1.75])
         
        phase = numpy.array([[4.0,4.0,5.0,5.0,3.0,3.0,4.0,4.0],
                             [3.0,3.0,4.0,4.0,3.0,4.0,4.0,5.0],
                             [2.0,2.0,3.0,3.0,4.0,4.0,5.0,5.0],
                             [5.0,5.0,2.0,2.0,4.0,5.0,5.0,2.0],
                             [2.0,2.0,3.0,3.0,5.0,5.0,2.0,2.0],
                             [5.0,2.0,2.0,3.0,4.0,4.0,5.0,5.0],
                             [5.0,5.0,2.0,2.0,3.0,3.0,4.0,4.0],
                             [2.0,5.0,5.0,2.0,2.0,2.0,3.0,3.0]],dtype=float)

#         ues = numpy.array([0.0,0.75 ,0.0,0.0])
#          
#         phase = numpy.array([[4.0,4.0,3.0,3.0,2.0,2.0,5.0,5.0],
#                              [5.0,4.0,4.0,3.0,3.0,2.0,2.0,5.0],
#                              [5.0,5.0,4.0,4.0,3.0,3.0,2.0,2.0],
#                              [2.0,5.0,5.0,4.0,4.0,3.0,3.0,2.0],
#                              [2.0,2.0,5.0,5.0,4.0,4.0,3.0,3.0],
#                              [3.0,2.0,2.0,5.0,5.0,4.0,4.0,3.0],
#                              [3.0,3.0,2.0,2.0,5.0,5.0,4.0,4.0],
#                              [4.0,3.0,3.0,2.0,2.0,5.0,5.0,4.0]],dtype=float)

        gaintx = numpy.zeros([8,8])
        gaintx[0:4,0:4] = 1
        gaintx[4:8,4:8] = 1
#         gaintx[4,6] = 0
        
        gainrx = numpy.zeros([8,8])
        gainrx[0:4,0:4] = 1
        gainrx[4:8,4:8] = 1
#         gainrx[4,6] = 0
        
        justrx = 0

    elif pattern == 218:
        # MST_ISR_EEJ (Abril 2017)
        title=" RX_NS_UP:CH3"
        
        ues = numpy.array([4.5,0.0,4.5,0.0])*2/3
#        ues = numpy.array([3.0,0.0,3.0,0.0])*1
        
        phase = numpy.array([[4,5,2,3,2,5,4,3],
                             [5,2,3,4,3,2,5,4],
                             [2,3,4,5,4,3,2,5],
                             [3,4,5,2,4,4,3,2],
                             [2,5,4,3,4,5,2,3],
                             [3,2,5,4,5,2,3,4],
                             [4,3,2,5,2,3,4,5],
                             [4,4,3,2,3,4,5,2]],dtype=float)

        gaintx = numpy.zeros([8,8])
        gaintx[:,:] = 1
        
        gainrx = numpy.zeros([8,8])
        gainrx[0:4,0:4] = 1
        gainrx[4:8,4:8] = 1
        
        justrx = 0

    elif pattern == 219:
        # MST_ISR_EEJ (Abril 2017)
        # MST_ISR_EEJ (Nov 2017)
        title=" RX_WE_UP:CH2"
        
        ues = numpy.array([4.5,0.0,4.5,0.0])*2/3
#        ues = numpy.array([3.0,0.0,3.0,0.0])*1
        
        phase = numpy.array([[4,5,2,3,2,5,4,3],
                             [5,2,3,4,3,2,5,4],
                             [2,3,4,5,4,3,2,5],
                             [3,4,5,2,4,4,3,2],
                             [2,5,4,3,4,5,2,3],
                             [3,2,5,4,5,2,3,4],
                             [4,3,2,5,2,3,4,5],
                             [4,4,3,2,3,4,5,2]],dtype=float)

        gaintx = numpy.zeros([8,8])
        gaintx[:,:] = 1
        
        gainrx = numpy.zeros([8,8])
        gainrx[4:8,0:4] = 1
        gainrx[0:4,4:8] = 1
        
        justrx = 0

    elif pattern == 220:
        # MST_ISR_EEJ (Abril 2017)
        title=" RX_NS_DW:CH4 (ABS)"
        
#        ues = numpy.array([0.0,0.4,4.5,0.75])*2/3
        ues = numpy.array([0.0,0.4,4.5,0.0])*2/3
        ues = numpy.array([0.75,0.4,4.5,0.0])*2/3

        
#         phase = numpy.array([[3.41,3.41,3.41,3.41,3,4,5,5],
#                              [2.78,2.78,2.78,2.78,2,3,3,4],
#                              [2.15,2.15,2.15,2.15,5,5,2,3],
#                              [5.52,5.52,5.52,5.52,3,4,5,2],
#                              [2,3,4,4,4.89,4.89,4.89,4.89],
#                              [5,2,2,3,4.26,4.26,4.26,4.26],
#                              [4,4,5,2,3.63,3.63,3.63,3.63],
#                              [2,3,4,5,3,3,3,3]],dtype=float)

#         phase = numpy.array([[3.50,3.50,3.50,3.00,3,4,5,5],
#                              [3.00,3.00,3.00,2.50,2,3,3,4],
#                              [2.00,2.00,2.00,2.00,5,5,2,3],
#                              [5.50,5.50,5.50,5.50,3,4,5,2],
#                              [2,3,4,4,5.00,5.00,5.00,5.00],
#                              [5,2,2,3,4.50,4.50,4.00,4.00],
#                              [4,4,5,2,3.50,3.50,3.50,3.50],
#                              [2,3,4,5,3,3,3,3]],dtype=float)

        phase = numpy.array([[3.00,3.00,3.00,3.00,3,4,5,5],
                             [2.50,2.50,2.50,2.50,2,3,3,4],
                             [2.00,2.00,1.50,1.50,5,5,2,3],
                             [5.50,5.50,5.50,5.00,3,4,5,2],
                             [2,3,4,4,4.50,4.50,4.50,4.50],
                             [5,2,2,3,4.00,4.00,4.00,4.00],
                             [4,4,5,2,3.50,3.50,3.00,3.00],
                             [2,3,4,5,3,3,3,2]],dtype=float)

        gaintx = numpy.zeros([8,8])
        gaintx[:,:] = 1
        
        gainrx = numpy.zeros([8,8])
        gainrx[0:4,0:4] = 1
        gainrx[4:8,4:8] = 1
        
        justrx = 0

    elif pattern == 221:
        # MST_ISR_EEJ (Abril 2017)
        # MST_ISR_EEJ (Nov 2017)
        title=" RX_WE_DW:CH1"
        
#         ues = numpy.array([4.5,0.0,4.5,0.0])*2/3
        #ues = numpy.array([0.0,0.4,4.5,0.75])*2/3
        ues = numpy.array([0.75,0.4,4.5,0.0])*2/3
        
        phase = numpy.array([[3.0,3.0,3.0,3.0,3,4,5,5],
                             [2.5,2.5,2.5,2.5,2,3,3,4],
                             [2.0,2.0,1.5,1.5,5,5,2,3],
                             [5.5,5.5,5.5,5.0,3,4,5,2],
                             [2,3,4,4,4.5,4.5,4.5,4.5],
                             [5,2,2,3,4.0,4.0,4.0,4.0],
                             [4,4,5,2,3.5,3.5,3.0,3.0],
                             [2,3,4,5,3,3,3,2]],dtype=float)

        
        gaintx = numpy.zeros([8,8])
        gaintx[:,:] = 1
        
        gainrx = numpy.zeros([8,8])
        gainrx[4:8,0:4] = 1
        gainrx[0:4,4:8] = 1
        
        justrx = 0
                            
    elif pattern == 222:
        # JULIA 19/04/17
        # Sup quarter +1m
        title=" Tx NSup , Rx NSup "
    
        ues = numpy.array([0.0,0.0,0.0,0.75])
        
        phase = numpy.array([[4.0,4.0,5.0,5.0,4.0,3.0,4.0,4.0],
                             [3.0,3.0,4.0,4.0,3.0,3.0,4.0,4.0],
                             [2.0,2.0,3.0,3.0,4.0,3.0,3.0,2.0],
                             [5.0,5.0,2.0,2.0,5.0,4.0,2.0,2.0],
                             [4.0,4.0,4.0,4.0,2.0,2.0,5.0,3.0],
                             [4.0,4.0,4.0,4.0,5.0,5.0,2.0,2.0],
                             [4.0,4.0,4.0,4.0,4.0,4.0,5.0,5.0],
                             [4.0,4.0,4.0,4.0,3.0,3.0,4.0,4.0]],dtype=float)
        gaintx = numpy.zeros([8,8])
        gaintx[0:4,0:4] = 1
        gaintx[4:8,4:8] = 1
        gaintx[4,6] = 0
        
        gainrx = numpy.zeros([8,8])
        gainrx[0:4,0:4] = 1
        gainrx[4:8,4:8] = 1
        gainrx[4,6] = 0
        
        justrx = 0

    elif pattern == 223:
        # ON AXIS
 #       title=" TX:EW_UP RX: W_UP"
#        title=" TX:EW_UP RX: E_UP"
#        title=" TX:EW_UP RX: N_UP"
        title=" TX:EW_UP RX: S_UP"
        
        ues = numpy.array([2.0,3.0,3.0,2.0])*1
        ues = numpy.array([0.0,1.0,1.0,0.0])*1
        
        phase = numpy.array([[4,4,4,4,4,4,4,4],
                             [4,4,4,4,4,4,4,4],
                             [4,4,4,4,4,4,4,4],
                             [4,4,4,4,4,4,4,4],
                             [5,5,5,5,5,5,5,5],
                             [5,5,5,5,5,5,5,5],
                             [5,5,5,5,5,5,5,5],
                             [5,5,5,5,5,5,5,5]],dtype=float)

        gaintx = numpy.zeros([8,8])
#        gaintx[4:8,0:4] = 1
        gaintx[0:4,4:8] = 1
        gaintx[4:8,4:8] = 1
#        gaintx[0:4,0:4] = 1
        
        gainrx = numpy.zeros([8,8])
#        gainrx[4:8,0:4] = 1
        gainrx[0:4,4:8] = 1
#        gainrx[0:4,0:4] = 1
        gainrx[4:8,4:8] = 1        
                
        justrx = 0   
        
    elif pattern == 224:
        # ON AXIS
#        title=" TX:E_UP RX: W_UP"
#        title=" TX:E_UP RX: E_UP"
#        title=" TX:E_UP RX: N_UP"
#        title=" TX:W_UP RX: W_UP, N_UP, E_UP, S_UP"
        title=" TX:E_UP RX: W_UP, N_UP, E_UP, S_UP"
        
        ues = numpy.array([3.0,4.0,4.0,3.0])*2/3
        ues = numpy.array([2.0,3.0,3.0,2.0])*1
        
        phase = numpy.array([[4,4,4,4,4,4,4,4],
                             [4,4,4,4,4,4,4,4],
                             [4,4,4,4,4,4,4,4],
                             [4,4,4,4,4,4,4,4],
                             [5,5,5,5,5,5,5,5],
                             [5,5,5,5,5,5,5,5],
                             [5,5,5,5,5,5,5,5],
                             [5,5,5,5,5,5,5,5]],dtype=float)

        gaintx = numpy.zeros([8,8])
#        gaintx[4:8,0:4] = 1
        gaintx[0:4,4:8] = 1
#        gaintx[4:8,4:8] = 1
        
        gainrx = numpy.zeros([8,8])
#        gainrx[4:8,0:4] = 1
#        gainrx[0:4,4:8] = 1
        gainrx[0:4,0:4] = 1
#        gainrx[4:8,4:8] = 1  
                
        justrx = 0   
    
    elif pattern == 225:
        # ON AXIS
        title=" TX:W_DW RX: W_DW"
#        title=" TX:W_DW RX: E_DW"
#        title=" TX:W_DW RX: N_DW"
#        title=" TX:W_DW RX: S_DW"
        
        ues = numpy.array([2.0,3.0,3.0,2.0])*2/3
        
        phase = numpy.array([[2,2,2,2,2,2,2,2],
                             [2,2,2,2,2,2,2,2],
                             [2,2,2,2,2,2,2,2],
                             [2,2,2,2,2,2,2,2],
                             [3,3,3,3,3,3,3,3],
                             [3,3,3,3,3,3,3,3],
                             [3,3,3,3,3,3,3,3],
                             [3,3,3,3,3,3,3,3]],dtype=float)

        gaintx = numpy.zeros([8,8])
        gaintx[4:8,0:4] = 1
        
        gainrx = numpy.zeros([8,8])
        gainrx[4:8,0:4] = 1
#        gainrx[0:4,4:8] = 1
#        gainrx[0:4,0:4] = 1
#        gainrx[4:8,4:8] = 1
                
        justrx = 0

    elif pattern == 226:
        # HM driver
#        title=" TX:ALL_UP RX: W_UP"
#        title=" TX:ALL_UP RX: E_UP"
#        title=" TX:ALL_UP RX: N_UP"
        title=" TX:ALL_UP RX: S_UP"
#        title=" TX:ALL UP RX: ALL_UP"
        
        ues = numpy.array([0.0,1.5,1.5,0.0])*2/3
        
        phase = numpy.array([[4,4,4,4,4,4,4,4],
                             [4,4,4,4,4,4,4,4],
                             [4,4,4,4,4,4,4,4],
                             [4,4,4,4,4,4,4,4],
                             [5,5,5,5,5,5,5,5],
                             [5,5,5,5,5,5,5,5],
                             [5,5,5,5,5,5,5,5],
                             [5,5,5,5,5,5,5,5]],dtype=float)

        gaintx = numpy.zeros([8,8])+1
        
        gainrx = numpy.zeros([8,8])
#        gainrx[4:8,0:4] = 1
#        gainrx[0:4,4:8] = 1
#        gainrx[0:4,0:4] = 1
        gainrx[4:8,4:8] = 1  
                
        justrx = 0   
    
    elif pattern == 227:
        # HM driver
        title=" TX:W_DW RX: W_DW"
#        title=" TX:W_DW RX: E_DW"
#        title=" TX:W_DW RX: N_DW"
#        title=" TX:W_DW RX: S_DW"
        
        ues = numpy.array([0.0,1.5,1.5,0.0])*2/3
        
        phase = numpy.array([[2,2,2,2,2,2,2,2],
                             [2,2,2,2,2,2,2,2],
                             [2,2,2,2,2,2,2,2],
                             [2,2,2,2,2,2,2,2],
                             [3,3,3,3,3,3,3,3],
                             [3,3,3,3,3,3,3,3],
                             [3,3,3,3,3,3,3,3],
                             [3,3,3,3,3,3,3,3]],dtype=float)

        gaintx = numpy.zeros([8,8])
        gaintx[4:8,0:4] = 1
        
        gainrx = numpy.zeros([8,8])
        gainrx[4:8,0:4] = 1
#        gainrx[0:4,4:8] = 1
#        gainrx[0:4,0:4] = 1
#        gainrx[4:8,4:8] = 1
                
        justrx = 0

    elif pattern == 300:
        # EW_DRIFT 28/11/16
        title=" Tx all down, Rx E down "
        
        ues = numpy.array([3.0,4.5,4.5,3.0])*2/3
        
        phase = numpy.array([[4.0,4.0,3.0,3.0,4.0,4.0,3.0,3.0],
                             [5.0,4.0,4.0,3.0,5.0,4.0,4.0,3.0],
                             [5.0,5.0,4.0,4.0,5.0,5.0,4.0,4.0],
                             [2.0,5.0,5.0,4.0,2.0,5.0,5.0,4.0],
                             [4.0,4.0,3.0,3.0,4.0,4.0,3.0,3.0],
                             [5.0,4.0,4.0,3.0,5.0,4.0,4.0,3.0],
                             [5.0,5.0,4.0,4.0,5.0,5.0,4.0,4.0],
                             [2.0,5.0,5.0,4.0,2.0,5.0,5.0,4.0]],dtype=float)

        gaintx = numpy.zeros([8,8])
        gaintx[0:4,4:8] = 1
        gaintx[4:8,0:4] = 1
        
        gainrx = numpy.zeros([8,8])
        gainrx[0:4,4:8] = 1
        
        justrx = 0

    elif pattern == 301:
        # JULIA Sup is not ok

        title=" Tx NSup , Rx NSup "
    
        ues = numpy.array([0.75,0.0,0.0,0.0])
        
        phase = numpy.array([[4.0,4.0,5.0,5.0,4.0,3.0,4.0,4.0],
                             [3.0,3.0,4.0,4.0,3.0,3.0,4.0,4.0],
                             [2.0,2.0,3.0,3.0,4.0,3.0,3.0,2.0],
                             [5.0,5.0,2.0,2.0,5.0,4.0,2.0,2.0],
                             [2.0,2.0,5.0,3.0,2.0,2.0,5.0,3.0],
                             [5.0,5.0,2.0,2.0,5.0,5.0,2.0,2.0],
                             [4.0,4.0,5.0,5.0,4.0,4.0,5.0,5.0],
                             [3.0,3.0,4.0,4.0,3.0,3.0,4.0,4.0]],dtype=float)
        gaintx = numpy.zeros([8,8])
        gaintx[0:4,0:4] = 1
        gaintx[4:8,0:4] = 1
#        gaintx[4,6] = 0
        
        gainrx = numpy.zeros([8,8])
        gainrx[0:4,0:4] = 1
        gainrx[4:8,0:4] = 1
#        gainrx[4,6] = 0
        
        justrx = 0

    elif pattern == 302:
        
        title = " for Oblique ISR '3.0N'"    
         
        ues = numpy.array([1.5,3.0,3.0,1.5])*2/3
#        ues = numpy.array([0.0,0.0,0.0,0.0])*2/3
        phase = numpy.array([[5,3,2,5,5.00,4.66,4.33,4.00],
                          [3,2,5,4,4.66,4.33,4.00,3.66],
                          [2,5,4,3,4.33,4.00,3.66,3.33],
                          [5,4,3,2,4.00,3.66,3.33,3.00],
                          [5,4,3,2,3,2,2,5],
                          [4,3,2,5,4,3,2,5],
                          [3,2,5,4,3,2,5,4],
                          [2,5,4,3,2,5,4,3]],dtype=float)
         
        gaintx = numpy.zeros([8,8])
        gaintx[0:4,4:8] = 1
         
        gainrx = numpy.zeros([8,8])
        gainrx[0:4,4:8] = 1
         
        justrx = 0

    elif pattern == 303:
                            
        # 2017_07
        title = " TX: E, RX: N16"    
         
#        ues = numpy.array([1.,2.,2.,1.])
        ues = numpy.array([1.5,3.0,3.0,1.5])*2/3
        phase = numpy.array([[4,4,4,4,2,2,3,3],
                          [4,4,4,2,2,3,3,4],
                          [4,4,4,4,3,3,4,4],
                          [4,4,4,4,3,4,4,5],
                          [2,2,3,3,4,4,4,4],
                          [2,3,3,4,4,4,4,4],
                          [3,3,4,4,4,4,4,4],
                          [3,4,4,5,4,4,4,4]],dtype=float)
         
        gaintx = numpy.zeros([8,8])
        gaintx[0:4,4:8] = 1
#        gaintx[:,:] = 1
         
        gainrx = numpy.zeros([8,8])
#        gainrx[:,:] = 1
#        gainrx[0:4,4:8] = 1
#        gainrx[4:8,0:4] = 1
#        gainrx[0,0] = 1
        gainrx[7,7] = 1
         
        justrx = 0    

    elif pattern == 304:
# 4 jul DW
        # 2017_07
        #title = "Tx : Ndw, Rx : N1dw"    
        #title = "Tx : Ndw, Rx : Edw"
        title = "Tx : Ndw, Rx : S16dw"
        #title = "Tx : Ndw, Rx : Wdw"
         
        ues = numpy.array([0.0,0.0,0.0,0.0])*2/3
        
#         phase = numpy.array([[2,2,4,4,3,3,4,4],
#                           [2,4,4,5,3,4,4,5],
#                           [4,4,5,5,4,4,5,5],
#                           [4,5,5,2,4,5,5,2],
#                           [3,3,4,4,4,4,4,4],
#                           [3,4,4,5,4,4,4,4],
#                           [4,4,5,5,4,4,4,4],
#                           [4,5,5,2,4,4,4,2]],dtype=float)

        phase = numpy.array([[4,4,5,5,4,4,5,5],
                          [4,5,5,2,4,5,5,2],
                          [5,5,2,3,5,5,2,2],
                          [5,2,3,4,5,2,2,3],
                          [4,4,5,5,4,4,4,4],
                          [4,5,5,2,4,4,4,4],
                          [5,5,2,2,4,4,4,4],
                          [5,2,2,3,4,4,4,4]],dtype=float)
         
        gaintx = numpy.zeros([8,8])
        #gaintx[0:4,4:8] = 1
        gaintx[0:4,0:4] = 1
#        gaintx[:,:] = 1
         
        gainrx = numpy.zeros([8,8])
        #gainrx[0,0] = 1
        #gainrx[0:4,4:8] = 1
        gainrx[7,7] = 1
        #gainrx[4:8,0:4] = 1         
        justrx = 0  

    elif pattern == 305:
# 4jul UP                            
        # 2017_07
        title = "Tx : Nup, Rx : N1up"    
#        title = "Tx : Nup, RX : Eup"
#        title = "Tx : Nup, Rx : S16up"
#        title = "Tx : Nup, Rx : Wup"
         
        ues = numpy.array([0.0,0.0,0.0,0.0])*2/3
        
#         phase = numpy.array([[2,2,3,3,2,2,3,3],
#                           [2,3,3,4,2,3,3,4],
#                           [3,3,4,4,3,3,4,4],
#                           [3,4,4,5,3,4,4,5],
#                           [2,2,3,3,4,4,4,4],
#                           [2,3,3,4,4,4,4,4],
#                           [3,3,4,4,4,4,4,4],
#                           [3,4,4,5,4,4,4,5]],dtype=float)

        phase = numpy.array([[2,2,3,3,2,2,3,3],
                          [2,3,3,4,2,3,3,4],
                          [3,3,4,5,3,3,4,4],
                          [3,4,5,2,3,4,4,5],
                          [2,2,3,3,4,4,4,4],
                          [2,3,3,4,4,4,4,4],
                          [3,3,4,4,4,4,4,4],
                          [3,4,4,5,4,4,4,2]],dtype=float)
                 
        gaintx = numpy.zeros([8,8])
#        gaintx[0:4,4:8] = 1
        gaintx[0:4,0:4] = 1
#        gaintx[:,:] = 1
         
        gainrx = numpy.zeros([8,8])
        gainrx[0,0] = 1
#        gainrx[0:4,4:8] = 1
#        gainrx[7,7] = 1
#        gainrx[4:8,0:4] = 1
         
        justrx = 0    

    elif pattern == 306:
# 3jul UP
        
# 5.66->5.52
# 2.66->2.78
# 2.33-> 2.15
# 5.33 -> 4.89
#        title = "Tx : Nup, Rx : Eup"    
#        title = "Tx : Nup, Rx : Wup"
#        title = "Tx : Nup, Rx : N1up"
        title = "Tx : Nup, Rx : S16up"
         
        ues = numpy.array([0.0,0.0,0.0,0.0])*2/3
# 2.66->2.78
# 3.66->3.63
# 4.66->4.89
# 4.33-> 4.26
# 3.33 -> 3.41
        phase = numpy.array([[5.00,4.00,4.00,3.00,4.89,4.26,4.00,3.63],
                          [4.00,4.00,3.00,3.00,4.26,4.00,3.63,3.41],
                          [4.00,3.00,3.00,3.00,4.00,3.63,3.41,3.00],
                          [3.00,3.00,3.00,3.00,3.63,3.41,3.00,2.78],
                          [4.89,4.26,4.00,3.63,4,4,4,4],
                          [4.26,4.00,3.63,3.41,4,4,4,4],
                          [4.00,3.63,3.41,3.00,4,4,4,4],
                          [3.63,3.41,3.00,2.78,4,4,4,3.00]],dtype=float)
         
        gaintx = numpy.zeros([8,8])
        #gaintx[0:4,4:8] = 1
        gaintx[0:4,0:4] = 1
         
        gainrx = numpy.zeros([8,8])
#        gainrx[0:4,4:8] = 1
#        gainrx[4:8,0:4] = 1
#        gainrx[0,0] = 1
        gainrx[7,7] = 1
         
        justrx = 0
        
    elif pattern == 307:
# 3jul DW        
        
#        title = "Tx : Ndw, Rx : Edw"    
#        title = "Tx : Ndw, Rx : Wdw"
#        title = "Tx : Ndw, Rx : N1dw"
        title = "Tx : Ndw, Rx : S16dw"        
         
        ues = numpy.array([0.0,0.0,0.0,0.0])*2/3
# 2.33->2.15
# 2.66->2.78
# 3.66->3.63
# 4.66->4.89
# 4.33-> 4.26
# 3.33 -> 3.41
# 5.66 -> 5.52
# 5.33 -> 5.00
# 5.00 ->4.89

        phase = numpy.array([[3.00,2.00,2.00,2.00,2.78,2.15,2.00,5.52],
                          [2.00,2.00,2.00,5.00,2.15,2.00,5.52,5.00],
                          [2.00,2.00,5.00,5.00,2.00,5.52,5.00,4.89],
                          [2.00,5.00,5.00,4.00,5.52,5.00,4.89,4.26],
                          [2.78,2.15,2.00,5.52,4,4,4,4],
                          [2.15,2.00,5.52,5.00,4,4,4,4],
                          [2.00,5.52,5.00,4.89,4,4,4,4],
                          [5.52,5.00,4.89,4.26,4,4,4,4.00]],dtype=float)
         
        gaintx = numpy.zeros([8,8])
#        gaintx[0:4,4:8] = 1
        gaintx[0:4,0:4] = 1
         
        gainrx = numpy.zeros([8,8])
#        gainrx[0:4,4:8] = 1
#        gainrx[4:8,0:4] = 1
#        gainrx[0,0] = 1
        gainrx[7,7] = 1
         
        justrx = 0        
        
    elif pattern == 308:
        
#        title = " chD : Edw"    
#        title = " chH : Wdw"
#        title = " chB : N1dw"
        title = " chF : S16dw"        
         
        ues = numpy.array([0.0,0.0,0.0,0.0])*2/3
# 2.33->2.15
# 2.66->2.78
# 3.66->3.63
# 4.66->4.89
# 4.33-> 4.26
# 3.33 -> 3.41
# 5.66 -> 5.52
# 5.33 -> 5.00
# 5.00 ->4.89

        phase = numpy.array([[2.00,4,4,4,2.00,2.78,3.41,4.00],
                          [4,4,4,4,2.78,3.41,4.00,4.89],
                          [4,4,4,4,3.41,4.00,4.89,5.00],
                          [4,4,4,4,4.00,4.89,5.00,2.00],
                          [2.00,2.78,3.41,4.00,4,4,4,4],
                          [2.78,3.41,4.00,4.89,4,4,4,4],
                          [3.41,4.00,4.89,5.00,4,4,4,4],
                          [4.00,4.89,5.00,2.00,4,4,4,2.00]],dtype=float)
         
        gaintx = numpy.zeros([8,8])
        gaintx[0:4,4:8] = 1
         
        gainrx = numpy.zeros([8,8])
#        gainrx[0:4,4:8] = 1
#        gainrx[4:8,0:4] = 1
#        gainrx[0,0] = 1
        gainrx[7,7] = 1
         
        justrx = 0      

    elif pattern == 309:
        
        title = " chC : Eup"    
#        title = " chG : Wup"
#        title = " chA : N1up"
#        title = " chE : S16up"        
         
        ues = numpy.array([0.0,0.0,0.0,0.0])*2/3
# 2.33->2.15
# 2.66->2.78
# 3.66->3.63
# 4.66->4.89
# 4.33-> 4.26
# 3.33 -> 3.41
# 5.66 -> 5.52
# 5.33 -> 5.00
# 5.00 ->4.89

#         phase = numpy.array([[4.00,4,4,4,4.00,4.89,5.00,2.00],
#                           [4,4,4,4,4.89,5.00,2.00,2.78],
#                           [4,4,4,4,5.00,2.00,2.78,3.41],
#                           [4,4,4,4,2.00,2.78,3.41,4.00],
#                           [4.00,4.89,5.00,2.00,4,4,4,4],
#                           [4.89,5.00,2.00,2.78,4,4,4,4],
#                           [5.00,2.00,2.78,3.41,4,4,4,4],
#                           [2.00,2.78,3.41,4.00,4,4,4,4.00]],dtype=float)

        phase = numpy.array([[4.00,4,4,4,3.63,4.26,5.00,5.52],
                          [4,4,4,4,4.26,5.00,5.52,2.15],
                          [4,4,4,4,5.00,5.52,2.15,3.00],
                          [4,4,4,4,5.52,2.15,3.00,3.63],
                          [4.00,4.89,5.00,2.00,4,4,4,4],
                          [4.89,5.00,2.00,2.78,4,4,4,4],
                          [5.00,2.00,2.78,3.41,4,4,4,4],
                          [2.00,2.78,3.41,4.00,4,4,4,4.00]],dtype=float)
                 
        gaintx = numpy.zeros([8,8])
        gaintx[0:4,4:8] = 1
         
        gainrx = numpy.zeros([8,8])
        gainrx[0:4,4:8] = 1
#        gainrx[4:8,0:4] = 1
#        gainrx[0,0] = 1
#        gainrx[7,7] = 1
         
        justrx = 0      

    elif pattern == 310:
                    
        title = " Tx at Wup, Rx at Nup, Eup and Sup"
                 
        ues = numpy.array([0.0,0.0,0.0,0.0])*2/3

        phase = numpy.array([[4,3,2,5,4,3,2,5],
                          [3,2,5,4,3,2,5,4],
                          [2,5,4,3,2,5,4,3],
                          [5,4,3,2,5,4,3,2],
                          [4,3,2,5,4,3,2,5],
                          [3,2,5,4,3,2,5,4],
                          [2,5,4,3,2,5,4,3],
                          [5,4,3,2,5,4,3,2]],dtype=float)
                 
        gaintx = numpy.zeros([8,8])
        #gaintx[:,:] = 1
        #gaintx[0:4,0:4] = 1
        #gaintx[4:8,4:8] = 1
        gaintx[4:8,0:4] = 1
         
        gainrx = numpy.zeros([8,8])
        gainrx[:,:] = 1
#        gainrx[0:4,0:4] = 1
#        gainrx[4:8,4:8] = 1
#        gainrx[0:4,4:8] = 1
         
        justrx = 0  
        
    elif pattern == 311:
        
        title = " for all down"
        title = " for half antenna down"        
                 
        ues = numpy.array([0.0,0.0,0.0,0.0])*2/3

        phase = numpy.array([[3,2,5,4,3,2,5,4],
                          [2,5,4,3,2,5,4,3],
                          [5,4,3,2,5,4,3,2],
                          [4,3,2,5,4,3,2,5],
                          [3,2,5,4,3,2,5,4],
                          [2,5,4,3,2,5,4,3],
                          [5,4,3,2,5,4,3,2],
                          [4,3,2,5,4,3,2,5]],dtype=float)
                 
        gaintx = numpy.zeros([8,8])
        #gaintx[:,:] = 1
        gaintx[0:4,0:4] = 1
        gaintx[4:8,4:8] = 1
         
        gainrx = numpy.zeros([8,8])
        #gainrx[:,:] = 1
        gainrx[0:4,0:4] = 1
        gainrx[4:8,4:8] = 1
         
        justrx = 0  
                                   
    elif pattern == 312:
        
        title = " for Tx & Rx : Nup + Sup"    
                 
        #ues = numpy.array([5.0,0.0,3.0,4.5])*2/3
        ues = numpy.array([6.5,0.0,0.0,4.5])*2/3

#         phase = numpy.array([[4,3,2,5,4,4,5,5],
#                           [3,2,5,4,3,3,4,4],
#                           [2,5,4,3,2,2,3,3],
#                           [5,4,3,2,5,5,2,2],
#                           [5,5,2,2,5,4,3,2],
#                           [4,4,5,5,4,3,2,5],
#                           [3,3,4,4,3,2,5,4],
#                           [2,2,3,3,2,5,4,3]],dtype=float)
        phase = numpy.array([[4,3,2,5,4,5,2,3],
                          [3,2,5,4,3,4,5,2],
                          [2,5,4,3,2,3,4,5],
                          [5,4,3,2,5,2,3,3],
                          [4,5,2,3,5,4,3,2],
                          [3,4,5,2,4,3,2,5],
                          [2,3,4,5,3,2,5,4],
                          [5,2,3,3,2,5,4,3]],dtype=float)
                 
        gaintx = numpy.zeros([8,8])
        gaintx[0:4,0:4] = 1
        gaintx[4:8,4:8] = 1
         
        gainrx = numpy.zeros([8,8])
        gainrx[0:4,0:4] = 1
        gainrx[4:8,4:8] = 1
         
        justrx = 0  
        
    elif pattern == 313:
        
        title = " for Tx & Rx : Ndw + Sdw"        
                 
        ues = numpy.array([3.0,0.0,3.0,3.0])*2/3

        phase = numpy.array([[3,2,5,4,4,4,3,3],
                          [2,5,4,3,5,4,4,3],
                          [5,4,3,2,5,5,4,4],
                          [4,3,2,5,2,5,5,4],
                          [4,4,3,3,5,4,3,2],
                          [5,4,4,3,4,3,2,5],
                          [5,5,4,4,3,2,5,4],
                          [2,5,5,4,2,5,4,3]],dtype=float)        
                                 
        gaintx = numpy.zeros([8,8])
        gaintx[0:4,0:4] = 1
        gaintx[4:8,4:8] = 1
         
        gainrx = numpy.zeros([8,8])
        gainrx[0:4,0:4] = 1
        gainrx[4:8,4:8] = 1
         
        justrx = 0  

    elif pattern == 314:
        
        title = " for Tx : Eup + Wup, Rx : Eup"    
        
#        ues = numpy.array([5.0,0.0,3.0,4.5])*2/3
        ues = numpy.array([6.5,0.0,0.0,4.5])*2/3        

#         phase = numpy.array([[4,3,2,5,4,4,5,5],
#                           [3,2,5,4,3,3,4,4],
#                           [2,5,4,3,2,2,3,3],
#                           [5,4,3,2,5,5,2,2],
#                           [5,5,2,2,5,4,3,2],
#                           [4,4,5,5,4,3,2,5],
#                           [3,3,4,4,3,2,5,4],
#                           [2,2,3,3,2,5,4,3]],dtype=float)

        phase = numpy.array([[4,3,2,5,4,5,2,3],
                          [3,2,5,4,3,4,5,2],
                          [2,5,4,3,2,3,4,5],
                          [5,4,3,2,5,2,3,3],
                          [4,5,2,3,5,4,3,2],
                          [3,4,5,2,4,3,2,5],
                          [2,3,4,5,3,2,5,4],
                          [5,2,3,3,2,5,4,3]],dtype=float)
                        
        gaintx = numpy.zeros([8,8])
        gaintx[0:4,4:8] = 1
        gaintx[4:8,0:4] = 1
         
        gainrx = numpy.zeros([8,8])
        gainrx[0:4,4:8] = 1

        justrx = 0  

    elif pattern == 315:
        
        title = " for Tx : Eup + Wup, Rx : Wup"    
                 
#        ues = numpy.array([5.0,0.0,3.0,4.5])*2/3
        ues = numpy.array([6.5,0.0,0.0,4.5])*2/3

#         phase = numpy.array([[4,3,2,5,4,4,5,5],
#                           [3,2,5,4,3,3,4,4],
#                           [2,5,4,3,2,2,3,3],
#                           [5,4,3,2,5,5,2,2],
#                           [5,5,2,2,5,4,3,2],
#                           [4,4,5,5,4,3,2,5],
#                           [3,3,4,4,3,2,5,4],
#                           [2,2,3,3,2,5,4,3]],dtype=float)
        phase = numpy.array([[4,3,2,5,4,5,2,3],
                          [3,2,5,4,3,4,5,2],
                          [2,5,4,3,2,3,4,5],
                          [5,4,3,2,5,2,3,3],
                          [4,5,2,3,5,4,3,2],
                          [3,4,5,2,4,3,2,5],
                          [2,3,4,5,3,2,5,4],
                          [5,2,3,3,2,5,4,3]],dtype=float)              
        gaintx = numpy.zeros([8,8])
        gaintx[0:4,4:8] = 1
        gaintx[4:8,0:4] = 1
         
        gainrx = numpy.zeros([8,8])
        gainrx[4:8,0:4] = 1
         
        justrx = 0  
 
    elif pattern == 316:
        
        title = " for Tx : Edw + Wdw, Rx : Edw"        
                 
        ues = numpy.array([3.0,0.0,3.0,3.0])*2/3

        phase = numpy.array([[3,2,5,4,4,4,3,3],
                          [2,5,4,3,5,4,4,3],
                          [5,4,3,2,5,5,4,4],
                          [4,3,2,5,2,5,5,4],
                          [4,4,3,3,5,4,3,2],
                          [5,4,4,3,4,3,2,5],
                          [5,5,4,4,3,2,5,4],
                          [2,5,5,4,2,5,4,3]],dtype=float)                 
        gaintx = numpy.zeros([8,8])
        gaintx[0:4,4:8] = 1
        gaintx[4:8,0:4] = 1
         
        gainrx = numpy.zeros([8,8])
        gainrx[0:4,4:8] = 1
         
        justrx = 0  
        
    elif pattern == 317:
        
        title = " for Tx : Edw + Wdw, Rx : Wdw"        
                 
        ues = numpy.array([3.0,0.0,3.0,3.0])*2/3

        phase = numpy.array([[3,2,5,4,4,4,3,3],
                          [2,5,4,3,5,4,4,3],
                          [5,4,3,2,5,5,4,4],
                          [4,3,2,5,2,5,5,4],
                          [4,4,3,3,5,4,3,2],
                          [5,4,4,3,4,3,2,5],
                          [5,5,4,4,3,2,5,4],
                          [2,5,5,4,2,5,4,3]],dtype=float)
                 
        gaintx = numpy.zeros([8,8])
        gaintx[0:4,4:8] = 1
        gaintx[4:8,0:4] = 1
         
        gainrx = numpy.zeros([8,8])
        gainrx[4:8,0:4] = 1
         
        justrx = 0     

    elif pattern == 318:
        
        title = " for half antenna (N&S)"
        title = " for all antenna"        
                 
        ues = numpy.array([0.0,0.0,0.0,0.0])*2/3

        phase = numpy.array([[4,3,2,5,4,3,2,5],
                          [4,3,2,5,4,3,2,5],
                          [4,3,2,5,4,3,2,5],
                          [4,3,2,5,4,3,2,5],
                          [4,3,2,5,4,3,2,5],
                          [4,3,2,5,4,3,2,5],
                          [4,3,2,5,4,3,2,5],
                          [4,3,2,5,4,3,2,5]],dtype=float)                 
        gaintx = numpy.zeros([8,8])
        gaintx[:,:] = 1
        #gaintx[0:4,0:4] = 1
        #gaintx[4:8,4:8] = 1
        
         
        gainrx = numpy.zeros([8,8])
        gainrx[:,:] = 1
        #gainrx[0:4,0:4] = 1
        #gainrx[4:8,4:8] = 1
         
        justrx = 0  
 
    elif pattern == 319:
        
        title = " for half antenna (N&S)"        
                 
        ues = numpy.array([0.0,0.0,0.0,0.0])*2/3

        phase = numpy.array([[4,4,4,4,4,3,2,5],
                          [3,3,3,3,4,3,2,5],
                          [2,2,2,2,4,3,2,5],
                          [5,5,5,5,4,3,2,5],
                          [4,3,2,5,4,4,4,4],
                          [4,3,2,5,3,3,3,3],
                          [4,3,2,5,2,2,2,2],
                          [4,3,2,5,5,5,5,5]],dtype=float)        
        gaintx = numpy.zeros([8,8])
        #gaintx[:,:] = 1
        gaintx[0:4,0:4] = 1
        gaintx[4:8,4:8] = 1
        
         
        gainrx = numpy.zeros([8,8])
        #gainrx[:,:] = 1
        gainrx[0:4,0:4] = 1
        gainrx[4:8,4:8] = 1
         
        justrx = 0  

    elif pattern == 320:
        # MP 2017_09
#        title=" TX:ALL_UP RX: W_UP"
#        title=" TX:ALL_UP RX: E_UP"
#        title=" TX:ALL_UP RX: N_UP"
        title=" TX: ALL_DW, RX: QUARTER_Dw"
#        title=" TX:ALL UP RX: ALL_UP"
        
        ues = numpy.array([0.0,1.5,1.5,0.0])*2/3
        
        phase = numpy.array([[4,4,4,4,4,4,4,4],
                             [4,4,4,4,4,4,4,4],
                             [4,4,4,4,4,4,4,4],
                             [4,4,4,4,4,4,4,4],
                             [5,5,5,5,5,5,5,5],
                             [5,5,5,5,5,5,5,5],
                             [5,5,5,5,5,5,5,5],
                             [5,5,5,5,5,5,5,5]],dtype=float)

        gaintx = numpy.zeros([8,8])+1
        
        gainrx = numpy.zeros([8,8])
#        gainrx[4:8,0:4] = 1
#        gainrx[0:4,4:8] = 1
#        gainrx[0:4,0:4] = 1
        gainrx[4:8,4:8] = 1  
                
        justrx = 0   

    elif pattern == 321:
        # MP 2017_09
        #title=" for TX: ALL [DW] and RX: QUARTER [Dw]"
        #title=" for TX: ALL [DW] and RX: ALL [DW]"
        title=" for Medium Power" + '\n'+ "ALL DW transmission / ALL DW reception"
        title=" for Medium Power" + '\n'+ "ALL DW transmission / A Quarter DW reception"
        
        ues = numpy.array([3.0,1.5,1.5,3.0])*2/3
        
        phase = numpy.array([[5,5,5,5,5,5,5,5],
                             [5,5,5,5,5,5,5,5],
                             [5,5,5,5,5,5,5,5],
                             [5,5,5,5,5,5,5,5],
                             [4,4,4,4,4,4,4,4],
                             [4,4,4,4,4,4,4,4],
                             [4,4,4,4,4,4,4,4],
                             [4,4,4,4,4,4,4,4]],dtype=float)

        gaintx = numpy.zeros([8,8])+1
        
        gainrx = numpy.zeros([8,8])
        gainrx[4:8,0:4] = 1
#        gainrx[0:4,4:8] = 1
#        gainrx[0:4,0:4] = 1
#        gainrx[4:8,4:8] = 1  
                
        justrx = 0
 
    elif pattern == 322:
        # MP 2017_09-Imaging

#        title=" TX: E_UP, RX: E_UP_4"
#        title=" TX: E_UP, RX: E_UP_7"
        title=" TX: E_UP, RX: E_UP_8"
        
        ues = numpy.array([3.0,1.5,1.5,3.0])*2/3
        
        phase = numpy.array([[3,3,3,3,4,3,4,4],
                             [3,3,3,3,3,3,4,4],
                             [3,3,3,3,4,3,3,2],
                             [3,3,3,3,5,4,2,2],
                             [4,4,4,4,2,2,2,2],
                             [4,4,4,4,2,2,2,2],
                             [4,4,4,4,2,2,2,2],
                             [4,4,4,4,2,2,2,2]],dtype=float)

        gaintx = numpy.zeros([8,8])
        gaintx[0:4,4:8] = 1
        gaintx[0,7] = 0
        gaintx[1,6:8] = 0
        gaintx[2,4] = 0
        
        gainrx = numpy.zeros([8,8])
#        gainrx[0,7] = 1
#        gainrx[1,6] = 1
        gainrx[1,7] = 1
#        gainrx[4:8,0:4] = 1
#        gainrx[0:4,4:8] = 1
#        gainrx[0:4,0:4] = 1
#        gainrx[4:8,4:8] = 1  
                
        justrx = 0

    elif pattern == 323:
        # MP 2017_09-ESF

#        title=" TX: N_UP, RX: W_UP"
#        title=" TX: N_UP, RX: S_UP"
        title=" TX: N_UP, RX: E_UP"
        
        ues = numpy.array([3.0,1.5,1.5,3.0])*2/3
        
        phase = numpy.array([[3,3,3,3,4,3,4,4],
                             [3,3,3,3,3,3,4,4],
                             [3,3,3,3,4,3,3,2],
                             [3,3,3,3,5,4,2,2],
                             [2,2,2,2,2,2,2,2],
                             [2,2,2,2,2,2,2,2],
                             [2,2,2,2,2,2,2,2],
                             [2,2,2,2,2,2,2,2]],dtype=float)

        gaintx = numpy.zeros([8,8])
        gaintx[0:4,0:4] = 1
        
        gainrx = numpy.zeros([8,8])
#        gainrx[4:8,0:4] = 1
#         gainrx[7,0] = 0
#         gainrx[7,2] = 0

#        gainrx[4:8,4:8] = 1        
#        gainrx[4,6] = 0
  
        gainrx[0:4,4:8] = 1        
        gainrx[0,7] = 0
        gainrx[1,6] = 0
        gainrx[1,7] = 0
        gainrx[2,4] = 0
                
        justrx = 0

    elif pattern == 324:
        # MP 2017_09-Imaging
        #title=" 324, TX: N_UP, RX: E_UP_4"
        title=" TX: E_UP_4, RX: E_UP_4"
        
        ues = numpy.array([0.0,0.0,0.0,0.0])*2/3
        
        phase = numpy.array([[4,3,4,5,3,3,3,3],
                             [3,3,3,4,3,3,3,3],
                             [4,3,4,3,3,3,3,3],
                             [5,4,3,3,3,3,3,3],
                             [2,2,2,2,2,2,2,2],
                             [2,2,2,2,2,2,2,2],
                             [2,2,2,2,2,2,2,2],
                             [2,2,2,2,2,2,2,2]],dtype=float)

        gaintx = numpy.zeros([8,8])
#        gaintx[0:4,0:4] = 1
        gaintx[0,7] = 1
        
        gainrx = numpy.zeros([8,8])        
        gainrx[0,7] = 1
                
        justrx = 0

    elif pattern == 325:
        # Imaging-JULIA
        
        title=" for Imaging experiment" + '\n'+ "A quarter transmission, Antenna module reception"
        
        ues = numpy.array([0.0,0.0,0.0,0.0])*2/3

#         phase = numpy.array([[3,3,3,3,4,3,4,5],
#                              [3,3,3,4,3,3,3,4],
#                              [3,3,3,3,4,3,3,2],
#                              [3,3,3,3,5,4,2,2],
#                              [2,2,2,2,2,2,2,2],
#                              [2,2,2,2,2,2,2,2],
#                              [2,2,2,2,2,2,2,2],
#                              [2,2,2,2,2,2,2,2]],dtype=float)
                
        phase = numpy.array([[3,3,3,3,4,3,5,2],
                             [3,3,3,4,3,3,3,4],
                             [3,3,3,3,5,4,3,3],
                             [3,3,3,3,2,4,3,3],
                             [2,2,2,2,2,2,2,2],
                             [2,2,2,2,2,2,2,2],
                             [2,2,2,2,2,2,2,2],
                             [2,2,2,2,2,2,2,2]],dtype=float)


        gaintx = numpy.zeros([8,8])
        gaintx[0:4,4:8] = 1
        gaintx[0,7] = 0
        gaintx[1,7] = 0
        gaintx[1,6] = 0
        gaintx[2,4] = 0
        
        gainrx = numpy.zeros([8,8])        
        gainrx[0,7] = 1
                
        justrx = 0

    elif pattern == 326:
        # MP 2017_09
        #title=" 324, TX: N_UP, RX: E_UP_4"
        title="  TX: E_UP_4, RX: W_UP_15"
        title=" for Medium Power" + '\n'+ "H module transmission / F module reception"
        title=" for Imaging experiment" + '\n'+ "Antenna module transmission, Antenna module reception"
        
        ues = numpy.array([0.0,0.0,0.0,0.0])*2/3
        
        phase = numpy.array([[3,3,3,3,3,3,3,3],
                             [3,3,3,3,3,3,3,3],
                             [3,3,3,3,3,3,3,3],
                             [3,3,3,3,3,3,3,3],
                             [2,2,2,2,2,2,2,2],
                             [2,2,2,2,2,2,2,2],
                             [2,2,2,2,2,2,2,2],
                             [2,2,2,2,2,2,2,2]],dtype=float)

        gaintx = numpy.zeros([8,8])
#        gaintx[0:4,0:4] = 1
        gaintx[0,7] = 1
        
        gainrx = numpy.zeros([8,8])        
        gainrx[7,2] = 1
                
        justrx = 0

    elif pattern == 327:
        # MP 2017_09

        title=" for Medium Power" + '\n'+ "All UP transmission / All UP reception"
        title=" for Medium Power" + '\n'+ "All transmission / All reception"
        
        ues = numpy.array([0.0,0.0,0.0,0.0])*2/3
        
        phase = numpy.array([[4,3,2,5,4,3,2,5],
                          [5,4,3,2,5,4,3,2],
                          [2,5,4,3,2,5,4,3],
                          [3,2,5,4,3,2,5,4],
                          [4,3,2,5,4,3,2,5],
                          [5,4,3,2,5,4,3,2],
                          [2,5,4,3,2,5,4,3],
                          [3,2,5,4,3,2,5,4]],dtype=float)

#         ues = numpy.array([0.0,0.0,0.0,0.0])*2/3

#         phase = numpy.array([[4,3,2,5,4,3,2,5],
#                           [4,3,2,5,4,3,2,5],
#                           [4,3,2,5,4,3,2,5],
#                           [4,3,2,5,4,3,2,5],
#                           [4,3,2,5,4,3,2,5],
#                           [4,3,2,5,4,3,2,5],
#                           [4,3,2,5,4,3,2,5],
#                           [4,3,2,5,4,3,2,5]],dtype=float)  

        gaintx = numpy.zeros([8,8])+1
#        gaintx[0:4,0:4] = 1
#        gaintx[0,7] = 1
        
        gainrx = numpy.zeros([8,8])+1        
#        gainrx[7,2] = 1
        #gainrx[0:4,0:4] = 1
                
        justrx = 0
 
    elif pattern == 328:
        # MP 2017_09

        title=" for Medium Power" + '\n'+ "All DW transmission / All DW reception"
#        title=" for Medium Power" + '\n'+ "All transmission / A quarter reception"
        
        ues = numpy.array([0.0,0.0,0.0,0.0])*2/3
        
        phase = numpy.array([[4,4,3,3,2,2,5,5],
                           [5,4,4,3,3,2,2,5],
                           [5,5,4,4,3,3,2,2],
                           [2,5,5,4,4,3,3,2],
                           [2,2,5,5,4,4,3,3],
                           [3,2,2,5,5,4,4,3],
                           [3,3,2,2,5,5,4,4],
                           [4,3,3,2,2,5,5,4]],dtype=float)

        gaintx = numpy.zeros([8,8])+1
        #gaintx[0:4,0:8] = 1
        #gaintx[4:8,0:8] = 1
        #gaintx[4:8,0:4]=1
        #gaintx[4:8,4:8]=1
        
        gainrx = numpy.zeros([8,8])+1    
        #gainrx[0:4,0:8] = 1
        #gainrx[4:8,0:8] = 1
        #gainrx[4:8,0:4]=1
        #gainrx[4:8,4:8]=1
                
        justrx = 0

    elif pattern == 329:
        # MP 2017_09

        title=" for Medium Power" + '\n'+ "All UP transmission / All UP reception"
        title=" for Medium Power" + '\n'+ "All transmission / All reception"
        
        ues = numpy.array([1.5,3.0,3.0,1.5])*2/3
        
        phase = numpy.array([[4,5,2,3,4,5,2,3],
                             [3,4,5,2,3,4,5,2],
                             [2,3,4,5,2,3,4,5],
                             [5,2,3,4,5,2,3,4],
                             [4,5,2,3,4,5,2,3],
                             [3,4,5,2,3,4,5,2],
                             [2,3,4,5,2,3,4,5],
                             [5,2,3,4,5,2,3,4]],dtype=float)

#         ues = numpy.array([0.0,0.0,0.0,0.0])*2/3

#         phase = numpy.array([[4,3,2,5,4,3,2,5],
#                           [4,3,2,5,4,3,2,5],
#                           [4,3,2,5,4,3,2,5],
#                           [4,3,2,5,4,3,2,5],
#                           [4,3,2,5,4,3,2,5],
#                           [4,3,2,5,4,3,2,5],
#                           [4,3,2,5,4,3,2,5],
#                           [4,3,2,5,4,3,2,5]],dtype=float)  

        gaintx = numpy.zeros([8,8])+1
#        gaintx[0:4,0:4] = 1
#        gaintx[0,7] = 1
        
        gainrx = numpy.zeros([8,8])+1        
#        gainrx[7,2] = 1
        #gainrx[0:4,0:4] = 1
                
        justrx = 0
        
        
    elif pattern == 330:
        # MP 2017_09

        title=" for Medium Power" + '\n'+ "All UP transmission / All UP reception"
        title=" for Medium Power" + '\n'+ "All UP transmission / All UP reception"
        
        ues = numpy.array([0.0,0.0,0.0,0.0])*2/3
        
#         phase = numpy.array([[4,4,5,5,2,2,3,3],
#                              [3,4,4,5,5,2,2,3],
#                              [3,3,4,4,5,5,2,2],
#                              [4,3,3,4,2,5,5,2],
#                              [5,5,4,4,3,3,2,2],
#                              [2,5,5,4,4,3,3,2],
#                              [2,2,5,5,4,4,3,3],
#                              [3,2,2,5,5,4,4,3]],dtype=float)

#         ues = numpy.array([0.0,0.0,0.0,0.0])*2/3

        phase = numpy.array([[5,2,2,2,2,3,3,4],
                           [5,5,5,5,2,2,3,3],
                           [4,4,4,5,5,2,2,3],
                           [3,3,4,4,5,5,2,2],
                           [2,3,3,4,4,5,5,2],
                           [2,2,3,3,4,4,5,5],
                           [5,2,2,3,3,4,4,5],
                           [5,5,2,2,3,3,4,4]],dtype=float)

#         phase = numpy.array([[4,4,3,3,2,2,5,5],
#                              [5,4,4,3,3,2,2,5],
#                              [5,5,4,4,3,3,2,2],
#                              [2,5,5,4,4,3,3,2],
#                              [2,2,5,5,4,4,3,3],
#                              [3,2,2,5,5,4,4,3],
#                              [3,3,2,2,5,5,4,4],
#                              [4,3,3,2,2,5,5,4]],dtype=float)

        gaintx = numpy.zeros([8,8])+1
 #       gaintx[4:8,0:4] = 1
#        gaintx[0,7] = 1
        
        gainrx = numpy.zeros([8,8])+1       
#        gainrx[7,2] = 1
#        gainrx[4:8,0:4] = 1
                
        justrx = 0

    elif pattern == 331:
        # JULIA 15/07/16
        title=" for Medium Power" + '\n' + " All dw transmission, All dw reception "
        title=" for JULIA experiment" + '\n' + " NSdw transmission, NSdw reception "
        #ues = numpy.array([0.0,1.75 ,0.0,0.0])
        #ues = numpy.array([0.0,1.5 ,0.0,0.0])
        #ues = numpy.array([0.0,0.0 ,0.0,0.0])
        #ues = numpy.array([0.0,0.25 ,0.0,0.0])
        ues = numpy.array([0.0,1.0,0.0,0.0])
        
#         phase = numpy.array([[4,4,3,3,2,2,5,5],
#                              [5,4,4,3,3,2,2,5],
#                              [5,5,4,4,3,3,2,2],
#                              [2,5,5,4,4,3,3,2],
#                              [2,2,5,5,4,4,3,3],
#                              [3,2,2,5,5,4,4,3],
#                              [3,3,2,2,5,5,4,4],
#                              [4,3,3,2,2,5,5,4]],dtype=float)
#         phase = numpy.array([[4,4,3,3,2,2,5,5],
#                              [5,4,4,3,3,2,2,5],
#                              [5,5,4,4,3,3,2,2],
#                              [2,5,5,4,4,3,3,2],
#                              [2,2,5,5,4,4,3,3],
#                              [3,2,2,5,5,4,4,3],
#                              [3,3,2,2,5,5,4,4],
#                              [4,3,3,2,2,5,5,4]],dtype=float)
        phase = numpy.array([[3.0,3.0,2.0,2.0,2.0,2.0,5.0,5.0],
                             [4.0,3.0,3.0,2.0,3.0,2.0,2.0,5.0],
                             [4.0,4.0,3.0,3.0,3.0,3.0,2.0,2.0],
                             [5.0,4.0,4.0,3.0,4.0,3.0,3.0,2.0],
                             [2.0,2.0,5.0,5.0,4.0,4.0,3.0,3.0],
                             [3.0,2.0,2.0,5.0,5.0,4.0,4.0,3.0],
                             [3.0,3.0,2.0,2.0,5.0,5.0,4.0,4.0],
                             [4.0,3.0,3.0,2.0,2.0,5.0,5.0,4.0]],dtype=float)

        #gaintx = numpy.zeros([8,8])+1
        gaintx = numpy.zeros([8,8])
        gaintx[0:4,0:4] = 1
        gaintx[4:8,4:8] = 1
        
        #gainrx = numpy.zeros([8,8])+1
        gainrx = numpy.zeros([8,8])
        gainrx[0:4,0:4] = 1
        gainrx[4:8,4:8] = 1
        
        justrx = 0
    
    elif pattern == 332:
        # JULIA 15/07/16
        #title=" for JULIA experiment" + '\n' + " NSup transmission, NSup reception "
        #title=" for Medium Power" + '\n' + " All up transmi2,2,3,3ssion, All up reception "
        title=" for JULIA experiment" + '\n' + " NSup transmission, NSup reception "
        
        #ues = numpy.array([0.8,0.0,0.0,0.8])*2/3
        #ues = numpy.array([0.0,0.0,0.0,0.0])*2/3
        ues = numpy.array([0.0,0.0,0.0,1.5])
        
        phase = numpy.array([[4,4,5,5,2,2,3,3],
                             [3,3,4,4,5,5,2,2],
                             [2,2,3,3,4,4,5,5],
                             [5,5,2,2,3,3,4,4],
                             [3,3,4,4,5,5,2,2],
                             [2,2,3,3,4,4,5,5],
                             [5,5,2,2,3,3,4,4],
                             [4,4,5,5,2,2,3,3]],dtype=float)
#         phase = numpy.array([[4,4,5,5,2,2,3,3],
#                              [3,3,4,4,5,5,2,3],
#                              [2,2,3,3,4,4,5,5],
#                              [5,5,2,3,3,3,4,4],
#                              [4,4,5,5,2,2,3,3],
#                              [3,3,4,4,5,5,2,3],
#                              [2,2,3,3,4,4,5,5],
#                              [5,5,2,3,3,3,4,4]],dtype=float)
#         phase = numpy.array([[3,4,5,5,2,2,3,3],
#                              [3,3,4,4,5,5,2,3],
#                              [2,2,3,3,3,4,5,5],
#                              [5,5,2,2,3,3,4,4],
#                              [3,4,5,5,2,2,3,3],
#                              [3,3,4,4,5,5,2,3],
#                              [2,2,3,3,3,4,5,5],
#                              [5,5,2,3,3,3,4,4]],dtype=float)
#         phase = numpy.array([[4,4,5,5,2,2,3,3],
#                              [3,3,4,4,5,5,2,2],
#                              [2,2,3,3,4,4,5,5],
#                              [5,5,2,2,3,3,4,4],
#                              [4,4,5,5,2,2,3,3],
#                              [3,3,4,4,5,5,2,2],
#                              [2,2,3,3,4,4,5,5],
#                              [5,5,2,2,3,3,4,4]],dtype=float)
        
        
#         gaintx = numpy.zeros([8,8])+1
        gaintx = numpy.zeros([8,8])
        gaintx[0:4,0:4] = 1
        #gaintx[4:8,0:4] = 1
        gaintx[4:8,4:8] = 1
        
#         gainrx = numpy.zeros([8,8])+1
        gainrx = numpy.zeros([8,8])
        gainrx[0:4,0:4] = 1
        #gainrx[4:8,0:4] = 1
        gainrx[4,6] = 0
        gainrx[4:8,4:8] = 1
        
        justrx = 0

    elif pattern == 333:
        # JULIA 15/07/16
        title=" for JULIA experiment" + '\n' + " NSdw transmission, Edw reception "
        
        #ues = numpy.array([0.0,0.75 ,0.0,0.0])
        ues = numpy.array([0.0,1-0.0,0.0,0.0])
        
#         phase = numpy.array([[4.0,4.0,3.0,3.0,2.0,2.0,5.0,5.0],
#                              [5.0,4.0,4.0,3.0,3.0,2.0,2.0,5.0],
#                              [5.0,5.0,4.0,4.0,3.0,3.0,2.0,2.0],
#                              [2.0,5.0,5.0,4.0,4.0,3.0,3.0,2.0],
#                              [2.0,2.0,5.0,5.0,4.0,4.0,3.0,3.0],
#                              [3.0,2.0,2.0,5.0,5.0,4.0,4.0,3.0],
#                              [3.0,3.0,2.0,2.0,5.0,5.0,4.0,4.0],
#                              [4.0,3.0,3.0,2.0,2.0,5.0,5.0,4.0]],dtype=float)

        phase = numpy.array([[3.0,3.0,2.0,2.0,2.0,2.0,5.0,5.0],
                             [4.0,3.0,3.0,2.0,3.0,2.0,2.0,5.0],
                             [4.0,4.0,3.0,3.0,3.0,3.0,2.0,2.0],
                             [5.0,4.0,4.0,3.0,4.0,3.0,3.0,2.0],
                             [2.0,2.0,5.0,5.0,4.0,4.0,3.0,3.0],
                             [3.0,2.0,2.0,5.0,5.0,4.0,4.0,3.0],
                             [3.0,3.0,2.0,2.0,5.0,5.0,4.0,4.0],
                             [4.0,3.0,3.0,2.0,2.0,5.0,5.0,4.0]],dtype=float)

        gaintx = numpy.zeros([8,8])
        gaintx[0:4,0:4] = 1
        gaintx[4:8,4:8] = 1
        
        gainrx = numpy.zeros([8,8])
        gainrx[0:4,4:8] = 1
        
        justrx = 0

    elif pattern == 334:
        # JULIA 15/07/16
        title=" for JULIA experiment" + '\n' + " NSdw transmission, Wdw reception "
        
        #ues = numpy.array([0.0,1.25 ,0.0,0.0])
        #ues = numpy.array([0.0,0.75 ,0.0,0.0])
        #ues = numpy.array([0.0,0.0 ,0.0,0.0])
        ues = numpy.array([0.0,1-0.0,0.0,0.0])
        
#         phase = numpy.array([[4.0,4.0,3.0,3.0,2.0,2.0,5.0,5.0],
#                              [5.0,4.0,4.0,3.0,3.0,2.0,2.0,5.0],
#                              [5.0,5.0,4.0,4.0,3.0,3.0,2.0,2.0],
#                              [2.0,5.0,5.0,4.0,4.0,3.0,3.0,2.0],
#                              [2.0,2.0,5.0,5.0,4.0,4.0,3.0,3.0],
#                              [3.0,2.0,2.0,5.0,5.0,4.0,4.0,3.0],
#                              [3.0,3.0,2.0,2.0,5.0,5.0,4.0,4.0],
#                              [4.0,3.0,3.0,2.0,2.0,5.0,5.0,4.0]],dtype=float)
        phase = numpy.array([[3.0,3.0,2.0,2.0,2.0,2.0,5.0,5.0],
                             [4.0,3.0,3.0,2.0,3.0,2.0,2.0,5.0],
                             [4.0,4.0,3.0,3.0,3.0,3.0,2.0,2.0],
                             [5.0,4.0,4.0,3.0,4.0,3.0,3.0,2.0],
                             [2.0,2.0,5.0,5.0,4.0,4.0,3.0,3.0],
                             [3.0,2.0,2.0,5.0,5.0,4.0,4.0,3.0],
                             [3.0,3.0,2.0,2.0,5.0,5.0,4.0,4.0],
                             [4.0,3.0,3.0,2.0,2.0,5.0,5.0,4.0]],dtype=float)
#         phase = numpy.array([[4.0,4.0,3.0,3.0,2.0,2.0,5.0,5.0],
#                              [5.0,4.0,4.0,3.0,3.0,2.0,2.0,5.0],
#                              [5.0,5.0,4.0,4.0,3.0,3.0,2.0,2.0],
#                              [2.0,5.0,5.0,4.0,4.0,3.0,3.0,2.0],
#                              [5.0,5.0,4.0,4.0,4.0,4.0,3.0,3.0],
#                              [2.0,5.0,5.0,4.0,5.0,4.0,4.0,3.0],
#                              [2.0,2.0,5.0,5.0,5.0,5.0,4.0,4.0],
#                              [3.0,2.0,2.0,5.0,2.0,5.0,5.0,4.0]],dtype=float)

        gaintx = numpy.zeros([8,8])
        gaintx[0:4,0:4] = 1
        gaintx[4:8,4:8] = 1
        
        gainrx = numpy.zeros([8,8])
        gainrx[4:8,0:4] = 1
        
        justrx = 0

    elif pattern == 335:
        # Pase de sol 2017_10
        title=" for Sun experiment" + '\n' + " NSup transmission, NSup reception "
        title=" for N, E, S transmission / A quarter reception "
        
        #27/10
#         ues = numpy.array([0.0,0.0,0.0,0.0])*2/3
        #26/10
        ues = numpy.array([0.0,1.5,0.75,0.0])*2/3
        #27/10
#         phase = numpy.array([[2,2,3,3,4,4,5,5],
#                              [5,2,2,3,3,4,4,5],
#                              [5,5,2,2,3,3,4,4],
#                              [4,5,5,2,2,3,3,4],
#                              [4,3,3,2,2,2,3,3],
#                              [5,4,3,3,5,2,2,3],
#                              [2,5,4,3,5,5,2,2],
#                              [2,2,5,4,4,5,5,2]],dtype=float)
        #26/10 con cambio de ues + antena
#         phase = numpy.array([[2,2,3,3,4,5,5,5],
#                              [5,2,2,3,4,5,4,5],
#                              [5,5,2,2,4,4,5,4],
#                              [4,5,5,2,3,4,4,4],
#                              [4,3,3,2,3,3,3,3],
#                              [5,4,3,3,5,2,2,3],
#                              [2,5,4,3,5,5,2,2],
#                              [2,2,5,4,4,5,5,2]],dtype=float)
        #26/10 con cambio de ues
        phase = numpy.array([[2,2,3,3,4,4,5,5],
                             [5,2,2,3,3,4,4,5],
                             [5,5,2,2,3,3,4,4],
                             [4,5,5,2,2,3,3,4],
                             [4,3,3,2,2,2,3,3],
                             [5,4,3,3,5,2,2,3],
                             [2,5,4,3,5,5,2,2],
                             [2,2,5,4,4,5,5,2]],dtype=float)

        gaintx = numpy.zeros([8,8])
        gaintx[0:8,4:8] = 1
        gaintx[0:4,0:4] = 1
        
        gainrx = numpy.zeros([8,8])
        #gainrx[0:4,0:4] = 1
        #gainrx[0:4,4:8] = 1
        gainrx[4:8,4:8] = 1
                
        justrx = 0
        
    elif pattern == 336:
        # Pase de sol 2017_10
        title=" for W reception"
        #26/10
#         ues = numpy.array([0.0,0.0,0.0,0.0])*2/3
        #27/10
        ues = numpy.array([0.0,0.0,0.0,0.0])*2/3
        #26/10
#         phase = numpy.array([[2,5,5,4,4,4,5,5],
#                              [3,2,5,5,3,4,4,5],
#                              [4,3,2,5,3,3,4,4],
#                              [4,4,3,2,2,3,3,4],
#                              [4,3,3,2,2,2,3,3],
#                              [5,4,3,3,5,2,2,3],
#                              [2,5,4,3,5,5,2,2],
#                              [2,2,5,4,4,5,5,2]],dtype=float)
        #27/10
        phase = numpy.array([[2,5,5,4,4,4,5,5],
                             [3,2,5,5,3,4,4,5],
                             [4,3,2,5,3,3,4,4],
                             [4,4,3,2,2,3,3,4],
                             [4,3,3,2,2,2,3,3],
                             [5,4,4,3,5,2,2,3],
                             [2,5,4,4,5,5,2,2],
                             [2,2,5,5,4,5,5,2]],dtype=float)
        
        gaintx = numpy.zeros([8,8])
        gaintx[4:8,0:4] = 1
        
        gainrx = numpy.zeros([8,8])
        gainrx[4:8,0:4] = 1
#         gainrx[5,1] = 0
#         gainrx[6,1] = 0
#         gainrx[6,2] = 0
                
        justrx = 1

    elif pattern == 337:
        # MIMO 2017_11

        #title=" for Imaging experiment" + '\n'+ "W15 transmission,  E4 reception"
        title=" for Imaging experiment" + '\n'+ "W15 transmission,  W15 reception"
        
        ues = numpy.array([0.0,0.0,0.0,0.0])*2/3
        
        phase = numpy.array([[3,3,3,3,3,3,3,3],
                             [3,3,3,3,3,3,3,3],
                             [3,3,3,3,3,3,3,3],
                             [3,3,3,3,3,3,3,3],
                             [2,2,2,2,2,2,2,2],
                             [2,2,2,2,2,2,2,2],
                             [2,2,2,2,2,2,2,2],
                             [2,2,2,2,2,2,2,2]],dtype=float)

        gaintx = numpy.zeros([8,8])
#        gaintx[0:4,0:4] = 1 # N quarter
        #gaintx[0,7] = 1 # E4
        gaintx[7,2] = 1 # W15 
        
        gainrx = numpy.zeros([8,8])        
        gainrx[7,2] = 1 # W15
        #gainrx[0,7] = 1 # E4
                
        justrx = 0

    elif pattern == 339:
        # Pase de sol 2018_02
        title=" for Sun experiment" + '\n' + " NSup transmission, NSup reception "
        title=" for N, E, S transmission / A quarter reception "
        
        #12/02
        ues = numpy.array([0.0,0.0,0.75,1.5])*2/3

        #13/02
#        ues = numpy.array([0.0,0.0,0.0,0.0])*2/3
        #14/02
#        ues = numpy.array([0.0,1.5,0.75,0.0])*2/3
        #13/02
#         phase = numpy.array([[2,2,3,3,4,4,5,5],
#                              [5,2,2,3,3,4,4,5],
#                              [5,5,2,2,3,3,4,4],
#                              [4,5,5,2,2,3,3,4],
#                              [4,3,3,2,2,2,3,3],
#                              [5,4,4,3,5,2,2,3],
#                              [2,5,4,4,5,5,2,2],
#                              [2,2,5,5,4,5,5,2]],dtype=float)
        #12/02    
        phase = numpy.array([[2,2,3,3,4,4,5,5],
                             [5,2,2,3,3,4,4,5],
                             [5,5,2,2,3,3,4,4],
                             [4,5,5,2,2,3,3,4],
                             [4,3,3,2,2,2,3,3],
                             [5,4,4,3,5,2,2,3],
                             [2,5,4,4,5,5,2,2],
                             [2,2,5,5,4,5,5,2]],dtype=float)            
#         phase = numpy.array([[2.00,2.33,2.33,2.66,4.00,4.33,4.33,4.66],
#                              [2.33,2.33,2.66,2.66,4.33,4.33,4.66,4.66],
#                              [2.33,2.66,2.66,3.00,4.33,4.66,4.66,5.00],
#                              [2.66,2.66,3.00,3.00,4.66,4.66,5.00,5.00],
#                              [0.00,0.00,0.00,0.00,3.66,4.00,4.00,4.33],
#                              [0.00,0.00,0.00,0.00,4.00,4.00,4.33,4.33],
#                              [0.00,0.00,0.00,0.00,8.00,8.33,8.33,8.66],
#                              [0.00,0.00,0.00,0.00,8.33,8.33,8.66,8.66]],dtype=float)
#         phase = numpy.array([[4,4,4,4,4,4,4,4],
#                              [4,4,4,4,4,4,4,4],
#                              [5,5,5,5,5,5,5,5],
#                              [5,5,5,5,5,5,5,5],
#                              [4,3,3,2,3,3,3,3],
#                              [5,4,3,3,3,3,3,3],
#                              [2,5,4,3,4,4,4,4],
#                              [2,2,5,4,4,4,4,4]],dtype=float)        
#         phase = numpy.array([[2.00,2.33,2.33,2.66,4.00,4.33,4.33,4.66],
#                              [2.33,2.33,2.66,2.66,4.33,4.33,4.66,4.66],
#                              [2.33,2.66,2.66,3.00,4.33,4.66,4.66,5.00],
#                              [2.66,2.66,3.00,3.00,4.66,4.66,5.00,5.00],
#                              [0.00,0.00,0.00,0.00,3.66,4.00,4.00,4.33],
#                              [0.00,0.00,0.00,0.00,4.00,4.00,4.33,4.33],
#                              [0.00,0.00,0.00,0.00,4.00,4.33,4.33,4.66],
#                              [0.00,0.00,0.00,0.00,4.33,4.33,8.66,4.66]],dtype=float)
#         phase = numpy.array([[2.00,2.00,2.00,3.00,4.00,4.00,4.00,5.00],
#                              [2.00,2.00,3.00,3.00,4.00,4.00,5.00,5.00],
#                              [3.00,3.00,3.00,3.00,4.00,5.00,5.00,5.00],
#                              [3.00,3.00,3.00,3.00,5.00,5.00,5.00,5.00],
#                              [0.00,0.00,0.00,0.00,4.00,4.00,4.00,4.00],
#                              [0.00,0.00,0.00,0.00,4.00,4.00,4.00,4.00],
#                              [0.00,0.00,0.00,0.00,8.00,8.00,8.00,9.00],
#                              [0.00,0.00,0.00,0.00,8.00,8.00,9.00,9.00]],dtype=float)
#         phase = numpy.array([[4,4,5,5,4,4,5,5],
#                              [4,4,5,5,4,4,5,5],
#                              [4,4,5,5,4,4,5,5],
#                              [4,4,5,5,4,4,5,5],
#                              [4,3,3,2,5,5,2,2],
#                              [5,4,3,3,5,5,2,2],
#                              [2,5,4,3,5,5,2,2],
#                              [2,2,5,4,5,5,2,2]],dtype=float)

#         phase = numpy.array([[2,2,3,3,5,5,2,2],
#                              [2,2,3,3,5,5,2,2],
#                              [2,2,3,3,5,5,2,2],
#                              [2,2,3,3,5,5,2,2],
#                              [4,3,3,2,2,2,3,3],
#                              [5,4,3,3,2,2,3,3],
#                              [2,5,4,3,2,2,3,3],
#                              [2,2,5,4,2,2,3,3]],dtype=float)
#         phase = numpy.array([[3,3,3,3,4,4,4,4],
#                              [3,3,3,3,4,4,4,4],
#                              [4,4,4,4,5,5,5,5],
#                              [4,4,4,4,5,5,5,5],
#                              [4,3,3,2,3,3,3,3],
#                              [5,4,3,3,3,3,3,3],
#                              [2,5,4,3,4,4,4,4],
#                              [2,2,5,4,4,4,4,4]],dtype=float)        
        #14/02 con cambio de ues + antena
#         phase = numpy.array([[2,2,3,3,4,5,5,5],
#                              [5,2,2,3,4,5,4,5],
#                              [5,5,2,2,4,4,5,4],
#                              [4,5,5,2,3,4,4,4],
#                              [4,3,3,2,3,3,3,3],
#                              [5,4,3,3,5,2,2,3],
#                              [2,5,4,3,5,5,2,2],
#                              [2,2,5,4,4,5,5,2]],dtype=float)
        #14/02 con cambio de ues
#         phase = numpy.array([[2,2,3,3,4,4,5,5],
#                              [5,2,2,3,3,4,4,5],
#                              [5,5,2,2,3,3,4,4],
#                              [4,5,5,2,2,3,3,4],
#                              [4,3,3,2,2,2,3,3],
#                              [5,4,3,3,5,2,2,3],
#                              [2,5,4,3,5,5,2,2],
#                              [2,2,5,4,4,5,5,2]],dtype=float)

        gaintx = numpy.zeros([8,8])
        gaintx[0:8,4:8] = 1
        gaintx[0:4,0:4] = 1
        gaintx[6,5] = 0
        gaintx[6,7] = 0
        
        gainrx = numpy.zeros([8,8])
        #gainrx[0:8,4:8] = 1
        #gainrx[0:4,4:8] = 1
        gainrx[4:8,4:8] = 1
        #gainrx[0:4,0:4] = 1
#         gainrx[4:8,0:4] = 1
                
        justrx = 0

    elif pattern == 340:
        # Configuration for DVD position (2013)
        title = " for DVD position (2013)"    
        
        #ues = numpy.array([1.0,2.0,2.0,1.25])
        ues = numpy.array([3.5,3.8,0.2,4.0])*(2/3.0)
        phase = numpy.array([[2,2,5,5,4,4,3,3],
                             [2,5,5,4,4,3,3,2],
                             [5,5,4,4,3,3,2,2],
                             [5,4,4,3,3,2,2,5],
                             [5,5,4,4,3,3,2,2],
                             [5,4,4,3,3,2,2,5],  
                             [4,4,3,3,2,2,5,5],
                             [4,3,3,2,2,5,5,4]],dtype=float)

        gaintx = numpy.zeros([8,8])
        gaintx[:,:] = 1
        
        gainrx = numpy.zeros([8,8])
        gainrx[0:4,0:4] = 1

        justrx = 0

    elif pattern == 341:
        # 06/2018
        title=" for Medium Power" + '\n' + " All dw transmission, All dw reception "
        ues = numpy.array([0.0,0.0 ,0.0,0.0])
        
        phase = numpy.array([[4,4,3,3,2,2,5,5],
                             [5,4,4,3,3,2,2,5],
                             [5,5,4,4,3,3,2,2],
                             [2,5,5,4,4,3,3,2],
                             [2,2,5,5,4,4,3,3],
                             [3,2,2,5,5,4,4,3],
                             [3,3,2,2,5,5,4,4],
                             [4,3,3,2,2,5,5,4]],dtype=float)

        gaintx = numpy.zeros([8,8])+1
        
        gainrx = numpy.zeros([8,8])+1
        
        justrx = 0
    
    elif pattern == 342:
        # 06/2018
        title=" for Medium Power" + '\n' + " All up transmission, All up reception "

        ues = numpy.array([0.26,0.0,0.0,0.89])
        ues = numpy.array([0.53,0.0,0.0,0.89])
        
        phase = numpy.array([[3,4,5,5,2,2,3,3],
                             [3,3,4,4,5,5,2,3],
                             [2,2,3,3,3,4,5,5],
                             [5,5,2,2,3,3,4,4],
                             [3,4,5,5,2,2,3,3],
                             [3,3,4,4,5,5,2,3],
                             [2,2,3,3,3,4,5,5],
                             [5,5,2,3,3,3,4,4]],dtype=float)
          
        gaintx = numpy.zeros([8,8])+1
        
        gainrx = numpy.zeros([8,8])+1
        
        justrx = 0

    elif pattern == 343:
        # 06/2018
        title="  TX: E_UP_4, RX: W_UP_15"
        title=" for Medium Power" + '\n'+ "H module transmission / F module reception"
        #title=" for Imaging experiment" + '\n'+ "Antenna module transmission, Antenna module reception"
        
        ues = numpy.array([0.0,0.0,0.0,0.0])*2/3
        
        phase = numpy.array([[3,3,3,3,3,3,3,3],
                             [3,3,3,3,3,3,3,3],
                             [3,3,3,3,3,3,3,3],
                             [3,3,3,3,3,3,3,3],
                             [2,2,2,2,2,2,2,2],
                             [2,2,2,2,2,2,2,2],
                             [2,2,2,2,2,2,2,2],
                             [2,2,2,2,2,2,2,2]],dtype=float)

        gaintx = numpy.zeros([8,8])
#        gaintx[0:4,0:4] = 1
        gaintx[0,7] = 1
        
        gainrx = numpy.zeros([8,8])        
        gainrx[7,2] = 1
                
        justrx = 0
 
    elif pattern == 344:
        # 06/2018, 09/2018
        title="  TX: E_UP_4, RX: W_UP_15"
        title=" for Medium Power" + '\n' + " All up transmission, All up reception "
        title=" for all up"
        
        ues = numpy.array([0.0,0.0,0.0,0.0])*2/3
        
        phase = numpy.array([[4,3,2,5,4,3,2,5],
                             [3,2,5,4,3,2,5,4],
                             [2,5,4,3,2,5,4,3],
                             [5,4,3,2,5,4,3,2],
                             [4,3,2,5,4,3,2,5],
                             [3,2,5,4,3,2,5,4],
                             [2,5,4,3,2,5,4,3],
                             [5,4,3,2,5,4,3,2]],dtype=float)

        gaintx = numpy.zeros([8,8])
#        gaintx[0:4,0:4] = 1
        gaintx[0:8,0:8] = 1
        
        gainrx = numpy.zeros([8,8])        
        gainrx[0:8,0:8] = 1
                
        justrx = 0
                                                                                       
    elif pattern == 345:
        # 06/2018, 09/2018
        title="  TX: E_UP_4, RX: W_UP_15"
        title=" for Medium Power" + '\n' + " All dw transmission, All dw reception "
        title=" for all down"
        
        ues = numpy.array([0.0,0.0,0.0,0.0])*2/3
        
        phase = numpy.array([[3,2,5,4,3,2,5,4],
                             [2,5,4,3,2,5,4,3],
                             [5,4,3,2,5,4,3,2],
                             [4,3,2,5,4,3,2,5],
                             [3,2,5,4,3,2,5,4],
                             [2,5,4,3,2,5,4,3],
                             [5,4,3,2,5,4,3,2],
                             [4,3,2,5,4,3,2,5]],dtype=float)

        gaintx = numpy.zeros([8,8])
#        gaintx[0:4,0:4] = 1
        gaintx[0:8,0:8] = 1
        
        gainrx = numpy.zeros([8,8])        
        gainrx[0:8,0:8] = 1
                
        justrx = 0    

    elif pattern == 346:
        # valley 8 (Agosto 2018)

        title=" for Valley experiment" + '\n'+ "Tx : All up, Rx : wdw"
        
        ues = numpy.array([3.0,1.5,1.5,3.0])*2.0/3
        
        phase = numpy.array([[5,5,5,5,5,5,5,5],
                             [5,5,5,5,5,5,5,5],
                             [5,5,5,5,5,5,5,5],
                             [5,5,5,5,5,5,5,5],
                             [4,4,4,4,4,4,4,4],
                             [4,4,4,4,4,4,4,4],
                             [4,4,4,4,4,4,4,4],
                             [4,4,4,4,4,4,4,4]],dtype=float)

        gaintx = numpy.zeros([8,8])+1
        
        gainrx = numpy.zeros([8,8])
        gainrx[4:8,0:4] = 1
                
        justrx = 0                                                                                   

    elif pattern == 347:
        # valley 8 (Agosto 2018)

        title=" for Valley experiment" + '\n'+ "Tx : All up, Rx : nup"
        
        ues = numpy.array([1.5,0.0,0.0,1.5])*2/3
        
        phase = numpy.array([[3,3,3,3,3,3,3,3],
                             [3,3,3,3,3,3,3,3],
                             [3,3,3,3,3,3,3,3],
                             [3,3,3,3,3,3,3,3],
                             [2,2,2,2,2,2,2,2],
                             [2,2,2,2,2,2,2,2],
                             [2,2,2,2,2,2,2,2],
                             [2,2,2,2,2,2,2,2]],dtype=float)

        gaintx = numpy.zeros([8,8])+1
        
        gainrx = numpy.zeros([8,8])
        gainrx[0:4,0:4] = 1
                
        justrx = 0
        
    elif pattern == 348:
        # Configuration for oblique ISR "6.0S"
        # Septiembre 2018
        title = " for Up Oblique ISR '6.0S'"    

#        ues = numpy.array([1.,2.,2.,1.])
        ues = numpy.array([1.5,3.0,3.0,1.5])*2/3
        
        phase = numpy.array([[4,5,2,3,4,5,2,3],
                          [5,2,3,4,5,2,3,4],
                          [2,3,4,5,2,3,4,5],
                          [3,4,5,2,3,4,5,2],
                          [5,2,3,4,5,2,3,4],
                          [2,3,4,5,2,3,4,5],
                          [3,4,5,2,3,4,5,2],
                          [4,5,2,3,4,5,2,3]],dtype=float)
        
        gaintx = numpy.zeros([8,8])
        gaintx[4:8,0:4] = 1
        gaintx[0:4,4:8] = 1
        
        gainrx = numpy.ones([8,8])
        #gainrx[:,:] = 1
        
        justrx = 0

    elif pattern == 349:
        # Configuration for oblique ISR "6.0S"
        # Septiembre 2018
        title = " for Down Oblique ISR '6.0S'"    

#        ues = numpy.array([1.,2.,2.,1.])
#        ues = numpy.array([1.5,3.0,3.0,1.5])*2/3
        ues = numpy.array([1.5,3.0,3.0,1.5])*2/3
        
        phase = numpy.array([[3,4,5,2,3,4,5,2],
                          [4,5,2,3,4,5,2,3],
                          [5,2,3,4,5,2,3,4],
                          [2,3,4,5,2,3,4,5],
                          [4,5,2,3,4,5,2,3],
                          [5,2,3,4,5,2,3,4],
                          [2,3,4,5,2,3,4,5],
                          [3,4,5,2,3,4,5,2]],dtype=float)
        
        gaintx = numpy.zeros([8,8])
        gaintx[4:8,0:4] = 1
        gaintx[0:4,4:8] = 1
        
        gainrx = numpy.ones([8,8])
        #gainrx[:,:] = 1
        
        justrx = 0

    elif pattern == 350:
        # Pase de sol 2018_02
        title=" for Sun experiment" + '\n' + " NSup transmission, NSup reception "
        title=" for N, E, S transmission / A quarter reception "
        
        #26/10
        ues = numpy.array([0.0,0.0,0.0,0.0])*2/3
        
         #26/10    
        phase = numpy.array([[3,3,4,4,5,5,2,2],
                             [2,3,3,4,4,5,5,2],
                             [5,5,2,2,3,3,4,4],
                             [4,5,5,2,2,3,3,4],
                             [4,4,5,5,2,2,3,3],
                             [3,4,4,5,5,2,2,3],
                             [3,3,4,4,5,5,2,2],
                             [2,3,3,4,4,5,5,2]],dtype=float)
               
 
        gaintx = numpy.zeros([8,8])+1
        #gaintx[0:8,4:8] = 1
        #gaintx[0:4,0:4] = 1
        #gaintx[6,5] = 0
        #gaintx[6,7] = 0
        
        gainrx = numpy.zeros([8,8])
        #gainrx[4:8,0:4] = 1
        #gainrx[0:4,0:4] = 1
        #gainrx[0:4,4:8] = 1
        gainrx[4:8,4:8] = 1
                
        justrx = 0

    elif pattern == 351:
        # Pase de sol 2018_02
        title=" for Sun experiment" + '\n' + " NSup transmission, NSup reception "
        title=" for N, E, S transmission / A quarter reception "
        
        #26/10
        ues = numpy.array([0.0,0.0,0.0,0.0])*2/3
        
         #26/10    
        phase = numpy.array([[2,2,3,3,4,4,5,5],
                             [5,2,2,3,3,4,4,5],
                             [5,5,2,2,3,3,4,4],
                             [4,5,5,2,2,3,3,4],
                             [4,4,5,5,2,2,3,3],
                             [3,4,4,5,5,2,2,3],
                             [3,3,4,4,5,5,2,2],
                             [2,3,3,4,4,5,5,2]],dtype=float)
               
 
        gaintx = numpy.zeros([8,8])+1
        #gaintx[0:8,4:8] = 1
        #gaintx[0:4,0:4] = 1
        #gaintx[6,5] = 0
        #gaintx[6,7] = 0
        
        gainrx = numpy.zeros([8,8])
        #gainrx[4:8,0:4] = 1
        #gainrx[0:4,0:4] = 1
        #gainrx[0:4,4:8] = 1
        gainrx[4:8,4:8] = 1
                
        justrx = 0

    elif pattern == 352:
        # Pase de sol 2018_02
        title=" for Sun experiment" + '\n' + " NSup transmission, NSup reception "
        title=" for N, E, S transmission / A quarter reception "
        
        #26/10
        ues = numpy.array([0.0,0.0,0.0,0.0])*2/3
        
         #26/10    
        phase = numpy.array([[4,3,3,2,4,4,5,5],
                             [5,4,4,3,3,4,4,5],
                             [3,5,4,4,3,3,4,4],
                             [2,3,5,5,2,3,3,4],
                             [4,4,5,5,4,3,3,2],
                             [3,4,4,5,5,4,4,3],
                             [3,3,4,4,3,5,4,4],
                             [2,3,3,4,2,3,5,5]],dtype=float)
               
 
        gaintx = numpy.zeros([8,8])+1
        #gaintx[0:8,4:8] = 1
        #gaintx[0:4,0:4] = 1
        #gaintx[6,5] = 0
        #gaintx[6,7] = 0
        
        gainrx = numpy.zeros([8,8])
        gainrx[4:8,0:4] = 1
        #gainrx[0:4,0:4] = 1
        #gainrx[0:4,4:8] = 1
        #gainrx[4:8,4:8] = 1
                
        justrx = 1

    elif pattern == 353:
        # Pase de sol 2018_02
        title=" for Sun experiment" + '\n' + " NSup transmission, NSup reception "
        title=" for N, E, S transmission / A quarter reception "
        
        #26/10
        ues = numpy.array([0.0,0.0,0.0,0.0])*2/3
        
         #26/10    
        phase = numpy.array([[2,2,3,3,4,5,2,5],
                             [5,2,2,3,3,5,5,5],
                             [5,5,2,2,3,4,5,4],
                             [4,5,5,2,2,4,4,4],
                             [4,4,5,5,2,3,4,3],
                             [3,4,4,5,5,3,3,3],
                             [3,3,4,4,5,2,3,2],
                             [2,3,3,4,4,2,2,2]],dtype=float)
               
 
        gaintx = numpy.zeros([8,8])+1
        #gaintx[0:8,4:8] = 1
        #gaintx[0:4,0:4] = 1
        #gaintx[6,5] = 0
        #gaintx[6,7] = 0
        
        gainrx = numpy.zeros([8,8])
        #gainrx[4:8,0:4] = 1
        gainrx[0:4,0:4] = 1
        #gainrx[0:4,4:8] = 1
        #gainrx[4:8,4:8] = 1
                
        justrx = 0

    elif pattern == 354:
        # Pase de sol 2018_02
        title=" for Sun experiment" + '\n' + " NSup transmission, NSup reception "
        title=" for N, E, S transmission / A quarter reception "
        
        #26/10
        ues = numpy.array([0.0,0.0,0.0,0.0])*2/3
        
         #26/10    
        phase = numpy.array([[4,3,3,2,4,5,2,5],
                             [5,4,4,3,3,5,5,5],
                             [3,5,4,4,3,4,5,4],
                             [2,3,5,5,2,4,4,4],
                             [4,4,5,5,4,3,3,2],
                             [3,4,4,5,5,4,4,3],
                             [3,3,4,4,3,5,4,4],
                             [2,3,3,4,2,3,5,5]],dtype=float)
               
 
        gaintx = numpy.zeros([8,8])+1
        #gaintx[0:8,4:8] = 1
        #gaintx[0:4,0:4] = 1
        #gaintx[6,5] = 0
        #gaintx[6,7] = 0
        
        gainrx = numpy.zeros([8,8])
        #gainrx[4:8,0:4] = 1
        gainrx[0:4,0:4] = 1
        #gainrx[0:4,4:8] = 1
        #gainrx[4:8,4:8] = 1
                
        justrx = 1

    elif pattern == 355:
        # Configuration for North Fritts"
        title = " for North (Fritts)"    

        ues = numpy.array([2.513, 1.0, 3.0, 0.413])
        phase = numpy.array([[4.29, 3.55, 2.82, 2.08, 4.20, 3.47, 2.73, 2.00],
                             [2.94, 2.20, 5.44, 4.70, 4.32, 3.59, 2.85, 2.12],
                             [5.56, 4.82, 4.09, 3.35, 4.44, 3.71, 2.97, 2.24],
                             [4.20, 3.47, 2.73, 2.00, 4.56, 3.82, 3.09, 2.35],
                             [4.20, 3.47, 2.73, 2.00, 4.56, 3.82, 3.09, 2.35],
                             [4.32, 3.59, 2.85, 2.12, 2.94, 2.20, 5.44, 4.70],
                             [4.44, 3.71, 2.97, 2.24, 5.56, 4.82, 4.09, 3.35],
                             [4.56, 3.82, 3.09, 2.35, 4.20, 3.47, 2.73, 2.00]],dtype=float)
        
        gaintx = numpy.zeros([8,8])
        gaintx[0:4,0:4] = 1
        gaintx[4:8,4:8] = 1
        
        gainrx = numpy.zeros([8,8])
        gainrx[0:4,0:4] = 1
        gainrx[4:8,4:8] = 1
        
        justrx = 0

    elif pattern == 356:
        # Configuration for West Fritts"
        title = " for West (Fritts)"    

        ues = numpy.array([2.513, 1.0, 3.0, 0.413])
        phase = numpy.array([[4.29, 3.55, 2.82, 2.08, 4.20, 3.47, 2.73, 2.00],
                             [2.94, 2.20, 5.44, 4.70, 4.32, 3.59, 2.85, 2.12],
                             [5.56, 4.82, 4.09, 3.35, 4.44, 3.71, 2.97, 2.24],
                             [4.20, 3.47, 2.73, 2.00, 4.56, 3.82, 3.09, 2.35],
                             [4.20, 3.47, 2.73, 2.00, 4.56, 3.82, 3.09, 2.35],
                             [4.32, 3.59, 2.85, 2.12, 2.94, 2.20, 5.44, 4.70],
                             [4.44, 3.71, 2.97, 2.24, 5.56, 4.82, 4.09, 3.35],
                             [4.56, 3.82, 3.09, 2.35, 4.20, 3.47, 2.73, 2.00]],dtype=float)
        
        gaintx = numpy.zeros([8,8])
        gaintx[4:8,0:4] = 1
#        gaintx[4:8,0:4] = 1
        gaintx[0:4,4:8] = 1

        gainrx = numpy.zeros([8,8])
        gainrx[4:8,0:4] = 1
        #gainrx[4:8,0:4] = 1
        gainrx[0:4,4:8] = 1

        justrx = 0

    elif pattern == 357:
        # Configuration for South Fritts"
        title = " for South (Fritts)"    

        ues = numpy.array([0.413, 2.0, 1.0, 1.513])
        phase = numpy.array([[2.0 , 2.73, 3.47, 4.2 , 2.08, 2.82, 3.55, 4.29],
                             [2.12, 2.85, 3.59, 4.32, 4.7 , 5.44, 2.20, 2.94],
                             [2.24, 2.97, 3.71, 4.44, 3.35, 4.09, 4.82, 5.56],
                             [2.35, 3.09, 3.82, 4.56, 2.0 , 2.73, 3.47, 4.20],
                             [2.08, 2.82, 3.55, 4.29, 2.0 , 2.73, 3.47, 4.20],
                             [4.70, 5.44, 2.20, 2.94, 2.12, 2.85, 3.59, 4.32],
                             [3.35, 4.09, 4.82, 5.56, 2.24, 2.97, 3.71, 4.44],
                             [2.00, 2.73, 3.47, 4.20, 2.35, 3.09, 3.82, 4.56]],dtype=float)

        gaintx = numpy.zeros([8,8])
        gaintx[0:4,0:4] = 1
        gaintx[4:8,4:8] = 1

        gainrx = numpy.zeros([8,8])
        gainrx[0:4,0:4] = 1
        gainrx[4:8,4:8] = 1

        justrx = 0

    elif pattern == 358:
        # Configuration for East Fritts"
        title = " for East (Fritts)"    

        ues = numpy.array([0.413, 2.0, 1.0, 1.513])
        phase = numpy.array([[2.0 , 2.73, 3.47, 4.2 , 2.08, 2.82, 3.55, 4.29],
                             [2.12, 2.85, 3.59, 4.32, 4.7 , 5.44, 2.20, 2.94],
                             [2.24, 2.97, 3.71, 4.44, 3.35, 4.09, 4.82, 5.56],
                             [2.35, 3.09, 3.82, 4.56, 2.0 , 2.73, 3.47, 4.20],
                             [2.08, 2.82, 3.55, 4.29, 2.0 , 2.73, 3.47, 4.20],
                             [4.70, 5.44, 2.20, 2.94, 2.12, 2.85, 3.59, 4.32],
                             [3.35, 4.09, 4.82, 5.56, 2.24, 2.97, 3.71, 4.44],
                             [2.00, 2.73, 3.47, 4.20, 2.35, 3.09, 3.82, 4.56]],dtype=float)
        
        gaintx = numpy.zeros([8,8])
        gaintx[4:8,0:4] = 1
        gaintx[0:4,4:8] = 1

        gainrx = numpy.zeros([8,8])
        gainrx[4:8,0:4] = 1
        gainrx[0:4,4:8] = 1

        justrx = 0

    elif pattern == 359:
        # Configuration for North Fritts"
        title = " for North (Fritts)"    

        #ues = numpy.array([2.513, 1.0, 3.0, 0.413])*2/3
        ues = numpy.array([2.5, 1.0, 3.0, 0.4])*2/3
        #ues = numpy.array([3.0, 1.0, 3.0, 0.5])*2/3

        phase = numpy.array([[4.5,3.5,3.0,2.0,4.0,3.5,2.5,2.0],
                            [3.0,2.0,5.5,4.5,4.5,3.5,3.0,2.0],
                            [5.5,5.0,4.0,3.5,4.5,3.5,3.0,2.0],
                            [4.0,3.5,2.5,2.0,4.5,4.0,3.0,2.5],
                            [4.0,3.5,2.5,2.0,4.5,4.0,3.0,2.5],
                            [4.5,3.5,3.0,2.0,3.0,2.0,5.5,4.5],
                            [4.5,3.5,3.0,2.0,5.5,5.0,4.0,3.5],
                            [4.5,4.0,3.0,2.5,4.0,3.5,2.5,2.0]],dtype=float)
        
        gaintx = numpy.zeros([8,8])
        gaintx[0:4,0:4] = 1
        gaintx[4:8,4:8] = 1
        
        gainrx = numpy.zeros([8,8])
        gainrx[0:4,0:4] = 1
        gainrx[4:8,4:8] = 1
        
        justrx = 0

    elif pattern == 360:
        # Configuration for West Fritts"
        title = " for West (Fritts)"    

        ues = numpy.array([2.513, 1.0, 3.0, 0.413])

        phase = numpy.array([[4.5,3.5,3.0,2.0,4.0,3.5,2.5,2.0],
                            [3.0,2.0,5.5,4.5,4.5,3.5,3.0,2.0],
                            [5.5,5.0,4.0,3.5,4.5,3.5,3.0,2.0],
                            [4.0,3.5,2.5,2.0,4.5,4.0,3.0,2.5],
                            [4.0,3.5,2.5,2.0,4.5,4.0,3.0,2.5],
                            [4.5,3.5,3.0,2.0,3.0,2.0,5.5,4.5],
                            [4.5,3.5,3.0,2.0,5.5,5.0,4.0,3.5],
                            [4.5,4.0,3.0,2.5,4.0,3.5,2.5,2.0]],dtype=float)
        
        gaintx = numpy.zeros([8,8])
        gaintx[4:8,0:4] = 1
        gaintx[4:8,0:4] = 1

        gainrx = numpy.zeros([8,8])
        gainrx[4:8,0:4] = 1
        gainrx[4:8,0:4] = 1

        justrx = 0

    elif pattern == 361:
        # Configuration for South Fritts"
        title = " for South (Fritts)"    

        #ues = numpy.array([0.413, 2.0, 1.0, 1.513])*2/3
        ues = numpy.array([0.5, 2.0, 1.0, 1.5])*2/3

        phase = numpy.array([[2.0,2.5,3.5,4.0,2.0,3.0,3.5,4.5],
                            [2.0,3.0,3.5,4.5,4.5,5.5,2.0,3.0],
                            [2.0,3.0,3.5,4.5,3.5,4.0,5.0,5.5],
                            [2.5,3.0,4.0,4.5,2.0,2.5,3.5,4.0],
                            [2.0,3.0,3.5,4.5,2.0,2.5,3.5,4.0],
                            [4.5,5.5,2.0,3.0,2.0,3.0,3.5,4.5],
                            [3.5,4.0,5.0,5.5,2.0,3.0,3.5,4.5],
                            [2.0,2.5,3.5,4.0,2.5,3.0,4.0,4.5]],dtype=float)

        gaintx = numpy.zeros([8,8])
        gaintx[0:4,0:4] = 1
        gaintx[4:8,4:8] = 1

        gainrx = numpy.zeros([8,8])
        gainrx[0:4,0:4] = 1
        gainrx[4:8,4:8] = 1

        justrx = 0

    elif pattern == 362:
        # Configuration for East Fritts"
        title = " for East (Fritts)"    

        ues = numpy.array([0.413, 2.0, 1.0, 1.513])

        phase = numpy.array([[2.0,2.5,3.5,4.0,2.0,3.0,3.5,4.5],
                            [2.0,3.0,3.5,4.5,4.5,5.5,2.0,3.0],
                            [2.0,3.0,3.5,4.5,3.5,4.0,5.0,5.5],
                            [2.5,3.0,4.0,4.5,2.0,2.5,3.5,4.0],
                            [2.0,3.0,3.5,4.5,2.0,2.5,3.5,4.0],
                            [4.5,5.5,2.0,3.0,2.0,3.0,3.5,4.5],
                            [3.5,4.0,5.0,5.5,2.0,3.0,3.5,4.5],
                            [2.0,2.5,3.5,4.0,2.5,3.0,4.0,4.5]],dtype=float)        
        gaintx = numpy.zeros([8,8])
        gaintx[4:8,0:4] = 1
        gaintx[4:8,0:4] = 1

        gainrx = numpy.zeros([8,8])
        gainrx[4:8,0:4] = 1
        gainrx[4:8,0:4] = 1

        justrx = 0

    elif pattern == 363:
        # Configuration for North Fritts"
        title = " for North (Fritts)"    

        ues = numpy.array([2.513, 1.0, 3.0, 0.413])

        phase = numpy.array([[4.0,4.0,3.0,2.0,4.0,3.0,3.0,2.0],
                            [3.0,2.0,5.0,5.0,4.0,4.0,3.0,2.0],
                            [6.0,5.0,4.0,3.0,4.0,4.0,3.0,2.0],
                            [4.0,3.0,3.0,2.0,5.0,4.0,3.0,2.0],
                            [4.0,3.0,3.0,2.0,5.0,4.0,3.0,2.0],
                            [4.0,4.0,3.0,2.0,3.0,2.0,5.0,5.0],
                            [4.0,4.0,3.0,2.0,6.0,5.0,4.0,3.0],
                            [5.0,4.0,3.0,2.0,4.0,3.0,3.0,2.0]],dtype=float)        
        gaintx = numpy.zeros([8,8])
        gaintx[0:4,0:4] = 1
        gaintx[4:8,4:8] = 1
        
        gainrx = numpy.zeros([8,8])
        gainrx[0:4,0:4] = 1
        gainrx[4:8,4:8] = 1
        
        justrx = 0

    elif pattern == 364:
        # Configuration for West Fritts"
        title = " for West (Fritts)"    

        ues = numpy.array([2.513, 1.0, 3.0, 0.413])

        phase = numpy.array([[4.0,4.0,3.0,2.0,4.0,3.0,3.0,2.0],
                            [3.0,2.0,5.0,5.0,4.0,4.0,3.0,2.0],
                            [6.0,5.0,4.0,3.0,4.0,4.0,3.0,2.0],
                            [4.0,3.0,3.0,2.0,5.0,4.0,3.0,2.0],
                            [4.0,3.0,3.0,2.0,5.0,4.0,3.0,2.0],
                            [4.0,4.0,3.0,2.0,3.0,2.0,5.0,5.0],
                            [4.0,4.0,3.0,2.0,6.0,5.0,4.0,3.0],
                            [5.0,4.0,3.0,2.0,4.0,3.0,3.0,2.0]],dtype=float)
        gaintx = numpy.zeros([8,8])
        gaintx[4:8,0:4] = 1
        gaintx[4:8,0:4] = 1

        gainrx = numpy.zeros([8,8])
        gainrx[4:8,0:4] = 1
        gainrx[4:8,0:4] = 1

        justrx = 0

    elif pattern == 365:
        # Configuration for South Fritts"
        title = " for South (Fritts)"    

        ues = numpy.array([0.413, 2.0, 1.0, 1.513])

        phase = numpy.array([[2.0,3.0,3.0,4.0,2.0,3.0,4.0,4.0],
                        [2.0,3.0,4.0,4.0,5.0,5.0,2.0,3.0],
                        [2.0,3.0,4.0,4.0,3.0,4.0,5.0,6.0],
                        [2.0,3.0,4.0,5.0,2.0,3.0,3.0,4.0],
                        [2.0,3.0,4.0,4.0,2.0,3.0,3.0,4.0],
                        [5.0,5.0,2.0,3.0,2.0,3.0,4.0,4.0],
                        [3.0,4.0,5.0,6.0,2.0,3.0,4.0,4.0],
                        [2.0,3.0,3.0,4.0,2.0,3.0,4.0,5.0]],dtype=float)
        
        gaintx = numpy.zeros([8,8])
        gaintx[0:4,0:4] = 1
        gaintx[4:8,4:8] = 1

        gainrx = numpy.zeros([8,8])
        gainrx[0:4,0:4] = 1
        gainrx[4:8,4:8] = 1

        justrx = 0

    elif pattern == 366:
        # Configuration for East Fritts"
        title = " for East (Fritts)"    

        ues = numpy.array([0.413, 2.0, 1.0, 1.513])

        phase = numpy.array([[2.0,3.0,3.0,4.0,2.0,3.0,4.0,4.0],
                            [2.0,3.0,4.0,4.0,5.0,5.0,2.0,3.0],
                            [2.0,3.0,4.0,4.0,3.0,4.0,5.0,6.0],
                            [2.0,3.0,4.0,5.0,2.0,3.0,3.0,4.0],
                            [2.0,3.0,4.0,4.0,2.0,3.0,3.0,4.0],
                            [5.0,5.0,2.0,3.0,2.0,3.0,4.0,4.0],
                            [3.0,4.0,5.0,6.0,2.0,3.0,4.0,4.0],
                            [2.0,3.0,3.0,4.0,2.0,3.0,4.0,5.0]],dtype=float)

        gaintx = numpy.zeros([8,8])
        gaintx[4:8,0:4] = 1
        gaintx[4:8,0:4] = 1

        gainrx = numpy.zeros([8,8])
        gainrx[4:8,0:4] = 1
        gainrx[4:8,0:4] = 1

        justrx = 0

    elif pattern == 367:
        # Configuration for North Fritts"
        title = " for North (Fritts)"    

        ues = numpy.array([2.513, 1.0, 3.0, 0.413])

        phase = numpy.array([[4.5,3.5,3.0,2.0,4.0,3.0,3.0,2.0],
                            [3.0,2.0,5.5,4.5,4.0,4.0,3.0,2.0],
                            [5.5,5.0,4.0,3.5,4.0,4.0,3.0,2.0],
                            [4.0,3.5,2.5,2.0,5.0,4.0,3.0,2.0],
                            [4.0,3.0,3.0,2.0,4.5,4.0,3.0,2.5],
                            [4.0,4.0,3.0,2.0,3.0,2.0,5.5,4.5],
                            [4.0,4.0,3.0,2.0,5.5,5.0,4.0,3.5],
                            [5.0,4.0,3.0,2.0,4.0,3.5,2.5,2.0]],dtype=float)
        
        gaintx = numpy.zeros([8,8])
        gaintx[0:4,0:4] = 1
        gaintx[4:8,4:8] = 1
        
        gainrx = numpy.zeros([8,8])
        gainrx[0:4,0:4] = 1
        gainrx[4:8,4:8] = 1
        
        justrx = 0

    elif pattern == 368:
        # Configuration for West Fritts"
        title = " for West (Fritts)"    

        #ues = numpy.array([2.513, 1.0, 3.0, 0.413])*2/3
        ues = numpy.array([2.5, 1.0, 3.0, 0.4])*2/3

        phase = numpy.array([[4.50,3.50,3.00,2.00,4.26,3.41,2.78,2.00],
                         [3.00,2.00,5.50,4.50,4.26,3.63,2.78,2.15],
                         [5.50,5.00,4.00,3.50,4.26,3.63,3.00,2.15],
                         [4.00,3.50,2.50,2.00,4.89,4.00,3.00,2.15],
                         [4.26,3.41,2.78,2.00,4.50,4.00,3.00,2.50],
                         [4.26,3.63,2.78,2.15,3.00,2.00,5.50,4.50],
                         [4.26,3.63,3.00,2.15,5.50,5.00,4.00,3.50],
                        [4.89,4.00,3.00,2.15,4.00,3.50,2.5,2.00]],dtype=float)
         
        gaintx = numpy.zeros([8,8])
        gaintx[4:8,0:4] = 1
        gaintx[0:4,4:8] = 1

        gainrx = numpy.zeros([8,8])
        gainrx[4:8,0:4] = 1
        gainrx[0:4,4:8] = 1

        justrx = 0

    elif pattern == 369:
        # Configuration for South Fritts"
        title = " for South (Fritts)"    

        ues = numpy.array([0.413, 2.0, 1.0, 1.513])

        phase = numpy.array([[2.0,2.5,3.5,4.0,2.0,3.0,4.0,4.0],
                            [2.0,3.0,3.5,4.5,5.0,5.0,2.0,3.0],
                            [2.0,3.0,3.5,4.5,3.0,4.0,5.0,6.0],
                            [2.5,3.0,4.0,4.5,2.0,3.0,3.0,4.0],
                            [2.0,3.0,4.0,4.0,2.0,2.5,3.5,4.0],
                            [5.0,5.0,2.0,3.0,2.0,3.0,3.5,4.5],
                            [3.0,4.0,5.0,6.0,2.0,3.0,3.5,4.5],
                            [2.0,3.0,3.0,4.0,2.5,3.0,4.0,4.5]],dtype=float)

        gaintx = numpy.zeros([8,8])
        gaintx[0:4,0:4] = 1
        gaintx[4:8,4:8] = 1

        gainrx = numpy.zeros([8,8])
        gainrx[0:4,0:4] = 1
        gainrx[4:8,4:8] = 1

        justrx = 0

    elif pattern == 370:
        # Configuration for East Fritts"
        title = " for East (Fritts)"    

        #ues = numpy.array([0.413, 2.0, 1.0, 1.513])*2/3
        ues = numpy.array([0.5, 2.0, 1.0, 1.5])*2/3

        phase = numpy.array([[2.0 , 2.73, 3.47, 4.2 ,2.00, 2.78, 3.63, 4.26],
                            [2.12, 2.85, 3.59, 4.32, 4.89, 5.52, 2.15, 3.00],
                            [2.24, 2.97, 3.71, 4.44, 3.41, 4.00, 4.89, 5.52],
                            [2.35, 3.09, 3.82, 4.56, 2.00, 2.78, 3.41, 4.00],
                            [2.00, 2.78, 3.63, 4.26, 2.0 , 2.73, 3.47, 4.20],
                            [4.89, 5.52, 2.15, 3.00, 2.12, 2.85, 3.59, 4.32],
                            [3.41, 4.00, 4.89, 5.52, 2.24, 2.97, 3.71, 4.44],
                            [2.00, 2.78, 3.41, 4.00, 2.35, 3.09, 3.82, 4.56]],dtype=float)
        
        gaintx = numpy.zeros([8,8])
        gaintx[4:8,0:4] = 1
        gaintx[0:4,4:8] = 1

        gainrx = numpy.zeros([8,8])
        gainrx[4:8,0:4] = 1
        gainrx[0:4,4:8] = 1

        justrx = 0
    
    elif pattern == 371:
        # Pase de sol 2019_02, variando solo ues
        title=" for All transmission / A quarter reception "
        
        ues = numpy.array([0.3,0.0,0.0,0.3])*2/3  #13/02
        #ues = numpy.array([0.0,0.0,0.0,0.0])*2/3  #14/02
        #ues = numpy.array([0.0,1.5,1.5,0.0])*2/3  #15/02
        #ues = numpy.array([0.0,0.0,1.5,1.5])*2/3
        
        # 13,14 / 02
        phase = numpy.array([[2,2,3,3,4,4,5,5],
                             [5,2,2,3,3,4,4,5],
                             [5,5,2,2,3,3,4,4],
                             [4,5,5,2,2,3,3,4],
                             [4,4,5,5,2,2,3,3],
                             [3,4,4,5,5,2,2,3],
                             [3,3,4,4,5,5,2,2],
                             [2,3,3,4,4,5,5,2]],dtype=float)

#         # 15 / 02
#         phase = numpy.array([[3,3,3,3,5,5,5,5],
#                              [5,2,2,3,3,4,4,5],
#                              [5,5,2,2,3,3,4,4],
#                              [4,5,5,2,2,3,3,4],
#                              [5,5,5,5,3,3,3,3],
#                              [3,4,4,5,5,2,2,3],
#                              [3,3,4,4,5,5,2,2],
#                              [2,3,3,4,4,5,5,2]],dtype=float)

        # 15 / 02
#         phase = numpy.array([[4,4,4,4,4,4,4,4],
#                              [4,4,4,4,4,4,4,4],
#                              [3,3,3,3,3,3,3,3],
#                              [3,3,3,3,3,3,3,3],
#                              [2,2,2,2,2,2,2,2],
#                              [2,2,2,2,2,2,2,2],
#                              [5,5,5,5,5,5,5,5],
#                              [5,5,5,5,5,5,5,5]],dtype=float)
      
 
        gaintx = numpy.zeros([8,8])
        #gaintx[0:8,4:8] = 1 
        gaintx[0:8,0:8] = 1
        #gaintx[6,5] = 0
        #gaintx[6,7] = 0
        
        gainrx = numpy.zeros([8,8])
        #gainrx[4:8,0:4] = 1 # W
        #gainrx[0:4,0:4] = 1 # N
        #gainrx[0:4,4:8] = 1 # E
        gainrx[4:8,4:8] = 1 # S
                
        justrx = 0
 
    elif pattern == 372:
        # Pase de sol 2019_02
        # [-1.67406905793, 14.5668994145] , 16.24096847243 
        # [ -2.12381895409, 14.6668438359] , 16.79066278999 
        
        title=" DTime = [-1.67406905793, 14.5668994145] , 16.24096847243 "        
        title=" DTime = [ -1.87395790067, 14.6668438359] , 16.54080173657 "        
        title=" DTime = [ -3.02331874642, 13.6174274115] , 16.64074615792 "
        
        title=" DTime = [ -3.22320758915, 13.6174274115] , 16.84063500065 "
        
        title=" DTime = [ -3.22320758915, 13.6174274115] , 16.84063500065 "
        title=" for All transmission / A quarter reception "
        
        ues = numpy.array([0.3,0.0,0.0,0.3])*2/3    #13/02
        #ues = numpy.array([0.0,0.0,0.0,0.0])*2/3  #14/02
        #ues = numpy.array([0.0,1.5,1.5,0.0])*2/3  #15/02
        #ues = numpy.array([0.0,0.0,0.0,0.0])*2/3  #14/02
        
         #13,14 / 02    
        phase = numpy.array([[4,3,3,2,5,5,2,2],
                             [5,4,4,3,4,5,5,2],
                             [3,5,4,4,3,3,4,4],
                             [2,3,5,5,2,3,3,4],
                             [4,4,5,5,4,3,3,2],
                             [3,4,4,5,5,4,4,3],
                             [3,3,4,4,3,5,4,4],
                             [2,3,3,4,2,3,5,5]],dtype=float)      
         #15 / 02    
#         phase = numpy.array([[5,4,4,2,5,5,2,2],
#                              [5,4,4,3,4,5,5,2],
#                              [3,5,4,4,3,3,4,4],
#                              [2,3,5,5,2,3,3,4],
#                              [4,4,5,5,5,4,4,2],
#                              [3,4,4,5,5,4,4,3],
#                              [3,3,4,4,3,5,4,4],
#                              [2,3,3,4,2,3,5,5]],dtype=float)            
 
        gaintx = numpy.zeros([8,8])
        #gaintx[0:8,4:8] = 1
        #gaintx[0:4,0:4] = 1
        #gaintx[6,5] = 0
        #gaintx[6,7] = 0
        
        gainrx = numpy.zeros([8,8])
        #gainrx[4:8,0:4] = 1
        #gainrx[0:4,0:4] = 1 # N
        #gainrx[0:4,4:8] = 1
        gainrx[4:8,4:8] = 1 # S
                
        justrx = 1 

    elif pattern == 373:
        # Pase de sol 2019_02, variando solo ues
        title=" for All transmission / A quarter reception "
        
        ues = numpy.array([0.3,0.0,0.0,0.3])*2/3  #13/02
        #ues = numpy.array([0.0,0.0,0.0,0.0])*2/3  #14/02
        ues = numpy.array([0.0,1.5,1.5,0.0])*2/3  #15/02
        #ues = numpy.array([0.0,0.0,0.0,0.0])*2/3
        
        # 13,14 / 02
#         phase = numpy.array([[2,2,3,3,4,4,5,5],
#                              [5,2,2,3,3,4,4,5],
#                              [5,5,2,2,3,3,4,4],
#                              [4,5,5,2,2,3,3,4],
#                              [4,4,5,5,2,2,3,3],
#                              [3,4,4,5,5,2,2,3],
#                              [3,3,4,4,5,5,2,2],
#                              [2,3,3,4,4,5,5,2]],dtype=float)
      
        # 15 / 02
        phase = numpy.array([[3,3,3,3,5,5,5,5],
                             [5,2,2,3,3,4,4,5],
                             [5,5,2,2,3,3,4,4],
                             [4,5,5,2,2,3,3,4],
                             [5,5,5,5,3,3,3,3],
                             [3,4,4,5,5,2,2,3],
                             [3,3,4,4,5,5,2,2],
                             [2,3,3,4,4,5,5,2]],dtype=float)
         
        gaintx = numpy.zeros([8,8])+1
        #gaintx[0:8,4:8] = 1
        #gaintx[0:4,0:4] = 1
        #gaintx[6,5] = 0
        #gaintx[6,7] = 0
        gaintx[0,0] = 0
        gaintx[2,0] = 0
        gaintx[5,1] = 0
        gaintx[6,1] = 0
        gaintx[6,2] = 0
        gaintx[6,5] = 0
        gaintx[6,7] = 0
        
        gainrx = numpy.zeros([8,8])
        gainrx[4:8,0:4] = 1  # W
        gainrx[5,1] = 0
        gainrx[6,1] = 0        
        gainrx[6,2] = 0        
#         gainrx[0:4,0:4] = 1 # N
#         gainrx[0,0] = 0
#         gainrx[2,0] = 0
#         gainrx[4:8,4:8] = 1 # S
#         gainrx[6,5] = 0
#         gainrx[6,7] = 0
         
        justrx = 0

    elif pattern == 374:
        # Pase de sol 2019_02
        # [-1.67406905793, 14.5668994145] , 16.24096847243 
        # [ -2.12381895409, 14.6668438359] , 16.79066278999 
        
        title=" DTime = [-1.67406905793, 14.5668994145] , 16.24096847243 "
        title=" DTime = [ -1.87395790067, 14.6668438359] , 16.54080173657 "
        title=" for All transmission / A quarter reception "
        title=" DTime = [ -3.02331874642, 13.6174274115] , 16.64074615792 "
        title=" DTime = [-1.67406905793, 14.0671773077] , 15.74124636563 "
        title=" DTime = [-1.67406905793, 14.7667882573] , 16.44085731523 "        
        
        
        title=" DTime = [-3.02331874642, 12.8678442513] , 15.891162997719999 "
        title=" DTime = [-3.17323537847, 13.3675663581] , 16.540801736577 "
        title=" for All transmission / A quarter reception "
        title=" DTime = [-3.17323537847, 12.8678442513] , 16.04107962977 "
        title=" DTime = [-1.87395790067, 14.0671773077] , 15.94113520837 "
        title=" for All transmission / A quarter reception "
        
        ues = numpy.array([0.3,0.0,0.0,0.3])*2/3    #13/02
        ues = numpy.array([3.3,3.0,3.0,3.3])*2/3    #13/02
        #ues = numpy.array([0.0,0.0,0.0,0.0])*2/3  #14/02
        ues = numpy.array([0,1.5,1.5,0.0])*2/3  #15/02
        #ues = numpy.array([0.0,0.0,0.0,0.0])*2/3  #14/02
        
         #13,14 / 02    
#         phase = numpy.array([[4,3,3,2,4,4,5,5],
#                              [5,4,4,3,3,4,4,5],
#                              [3,5,4,4,3,3,4,4],
#                              [2,3,5,5,2,3,3,4],
#                              [4,4,5,5,4,3,3,2],
#                              [3,4,4,5,5,4,4,3],
#                              [3,3,4,4,3,5,4,4],
#                              [2,3,3,4,2,3,5,5]],dtype=float)      
#          #15 / 02    
        phase = numpy.array([[5,4,4,2,5,5,2,2],
                             [5,4,4,3,4,5,5,2],
                             [3,5,4,4,3,3,4,4],
                             [2,3,5,5,2,3,3,4],
                             [4,4,5,5,5,4,4,2],
                             [3,4,4,5,5,4,4,3],
                             [3,3,4,4,2,5,5,4],
                             [2,3,3,4,2,3,5,5]],dtype=float)            
 
        gaintx = numpy.zeros([8,8])
        #gaintx[0:8,4:8] = 1
        #gaintx[0:4,0:4] = 1
        #gaintx[6,5] = 0
        #gaintx[6,7] = 0
        
        gainrx = numpy.zeros([8,8])
        #gainrx[4:8,0:4] = 1
#         gainrx[0:4,0:4] = 1 # N
#         gainrx[0,0] = 0
#         gainrx[2,0] = 0
        #gainrx[0:4,4:8] = 1
        gainrx[4:8,4:8] = 1 # S
        gainrx[6,5] = 0
        gainrx[6,7] = 0
                        
        justrx = 1 

    elif pattern == 375:
        # Pase de sol 2019_02, variando solo ues
        title=" for All down transmission / All down reception "
        
        ues = numpy.array([6.0,6.75,0.0,0.0])*2/3  #13/02

#         phase = numpy.array([[4,4,4,4,4,4,4,4],
#                              [4,4,4,4,4,4,4,4],
#                              [4,4,4,4,4,4,4,4],
#                              [4,4,4,4,4,4,4,4],
#                              [4,4,4,4,4,4,4,4],
#                              [4,4,4,4,4,4,4,4],
#                              [4,4,4,4,4,4,4,4],
#                              [4,4,4,4,4,4,4,4]],dtype=float)

#         phase = numpy.array([[4.00,4.33,4.33,4.66,4,4,4,4],
#                              [4.33,4.33,4.66,4.66,4,4,4,4],
#                              [4.33,4.66,4.66,5.00,4,4,4,4],
#                              [4.66,4.66,5.00,5.00,4,4,4,4],
#                              [4,4,4,4,4,4,4,4],
#                              [4,4,4,4,4,4,4,4],
#                              [4,4,4,4,4,4,4,4],
#                              [4,4,4,4,4,4,4,4]],dtype=float)
#         phase = numpy.array([[4.00,4.00,4.50,4.50,5.00,5.00,5.50,5.50],
#                              [4.00,4.50,4.50,5.00,5.00,5.50,5.50,2.00],
#                              [4.50,4.50,5.00,5.00,5.50,5.50,2.00,2.00],
#                              [4.50,5.00,5.00,5.50,5.50,2.00,2.00,2.50],
#                              [5.00,5.00,5.50,5.50,2.00,2.00,2.50,2.50],
#                              [5.00,5.50,5.50,2.00,2.00,2.50,2.50,3.00],
#                              [5.50,5.50,2.00,2.00,2.50,2.50,3.00,3.00],
#                              [5.50,2.00,2.00,2.50,2.50,3.00,3.00,3.50]],dtype=float)

        phase = numpy.array([[4.00,4.00,4.50,4.50,5.00,5.00,2.00,2.00],
                             [4.00,4.50,4.50,5.00,5.00,2.00,2.00,2.00],
                             [4.50,4.50,5.00,5.00,2.00,2.00,2.00,2.00],
                             [4.50,5.00,5.00,5.50,2.00,2.00,2.00,2.00],
                             [5.00,5.00,2.00,2.00,2.00,2.00,2.50,2.50],
                             [5.00,2.00,2.00,2.00,2.00,2.50,2.50,3.00],
                             [2.00,2.00,2.00,2.00,2.50,2.50,3.00,3.00],
                             [2.00,2.00,2.00,2.00,2.50,3.00,3.00,3.50]],dtype=float)
         
        gaintx = numpy.zeros([8,8])+1
        #gaintx[0:4,0:4] = 1  # N
        
        gainrx = numpy.zeros([8,8])+1
        #gainrx[0:4,0:4] = 1  # N

         
        justrx = 0
 
    elif pattern == 376:
        # Pase de sol 2019_02, variando solo ues
        title=" for All up reception "
        
        ues = numpy.array([6.0,6.75,0.0,0.0])*2/3  #13/02

        phase = numpy.array([[2.00,2.00,2.50,2.50,3.00,3.00,4.00,4.00],
                             [2.00,2.50,2.50,3.00,3.00,4.00,4.00,4.00],
                             [2.50,2.50,3.00,3.00,4.00,4.00,4.00,4.00],
                             [2.50,3.00,3.00,3.50,4.00,4.00,4.00,4.00],
                             [3.00,3.00,4.00,4.00,4.00,4.00,4.50,4.50],
                             [3.00,4.00,4.00,4.00,4.00,4.50,4.50,5.00],
                             [4.00,4.00,4.00,4.00,4.50,4.50,5.00,5.00],
                             [4.00,4.00,4.00,4.00,4.50,5.00,5.00,5.50]],dtype=float)
         
        gaintx = numpy.zeros([8,8])+1
        #gaintx[0:4,0:4] = 1  # N
        
        gainrx = numpy.zeros([8,8])+1
        #gainrx[0:4,0:4] = 1  # N

         
        justrx = 1
                                                                                                     
    elif pattern==None:
        
        inputs = numpy.array(["title","ues_tx","phase_tx","gain_tx","gain_rx","just_rx"])

        # Reading user-defined configuration.
        if path==None:path = os.getcwd() + os.sep + "patterns" + os.sep
        if filename==None:filename = "jropattern.txt"
        
        ff = open(os.path.join(path,filename),'r')
        
        while  1:
            # Checking EOF.
            init = ff.tell()
            if not ff.readline():break
            else:ff.seek(init)

            line = ff.readline().lstrip()
            if line.__len__()!=0:
                if line[0]!='#':                    
                    keys = line.split("=")
                    key = keys[0].lstrip().rstrip().lower()
                    vv = numpy.where(inputs==key)
                    if vv[0][0]==0:
                        title = keys[1].lstrip().rstrip()
                    elif vv[0][0]==1:
                        ues = (keys[1].lstrip().rstrip())
                        ues = numpy.float32(ues[1:-1].split(","))
                    elif vv[0][0]==2:
                        phase = numpy.zeros([8,8])
                        tx = (keys[1].lstrip().rstrip())
                        tx = numpy.float32(tx[2:-3].split(","))
                        phase[0,:] = tx
                        for ii in numpy.arange(7):
                            tx = ff.readline().lstrip().rstrip()
                            tx = numpy.float32(tx[1:-3+(ii==6)].split(","))
                            phase[ii+1,:] = tx
                    elif vv[0][0]==3:
                        gaintx = numpy.zeros([8,8])
                        gg = (keys[1].lstrip().rstrip())
                        gg = numpy.float32(gg[2:-3].split(","))
                        gaintx[0,:] = gg
                        for ii in numpy.arange(7):
                            gg = ff.readline().lstrip().rstrip()
                            gg = numpy.float32(gg[1:-3+(ii==6)].split(","))
                            gaintx[ii+1,:] = gg
                    elif vv[0][0]==4:
                        gainrx = numpy.zeros([8,8])
                        gg = (keys[1].lstrip().rstrip())
                        gg = numpy.float32(gg[2:-3].split(","))
                        gainrx[0,:] = gg
                        for ii in numpy.arange(7):
                            gg = ff.readline().lstrip().rstrip()
                            gg = numpy.float32(gg[1:-3+(ii==6)].split(","))
                            gainrx[ii+1,:] = gg
                    elif vv[0][0]==5:
                        justrx = numpy.float(keys[1].lstrip().rstrip())
        
        ff.close()
    
    
    
    setup = {"ues":ues, "phase":phase, "gaintx":gaintx, "gainrx":gainrx, "justrx":justrx, \
     "title":title}
    
    return setup
    
    
        
    
    