
import os
import time
import numpy
import scipy
import base64
import datetime
import scipy.interpolate
import math


#from .muf import *
from io import BytesIO
import gaussfit
import TimeTools
import Astro_Coords
import Misc_Routines
from matplotlib.figure import Figure
from patterns import select_pattern
import matplotlib.pyplot as plt
from PIL import Image


attenuation = numpy.array([[[-21.25,-15.25,-9.25,-3.25,3.25,9.25,15.25,21.25],
		         [-21.25,-15.25,-9.25,-3.25,03.25,09.25,15.25,21.25],
   		         [-21.25,-15.25,-9.25,-3.25,03.25,09.25,15.25,21.25],
		         [-21.25,-15.25,-9.25,-3.25,03.25,09.25,15.25,21.25],
		         [-21.25,-15.25,-9.25,-3.25,03.25,09.25,15.25,21.25],
		         [-21.25,-15.25,-9.25,-3.25,03.25,09.25,15.25,21.25],
		         [-21.25,-15.25,-9.25,-3.25,03.25,09.25,15.25,21.25],
		         [-21.25,-15.25,-9.25,-3.25,03.25,09.25,15.25,21.25]],
       		         [[21.25,21.25,21.25,21.25,21.25,21.25,21.25,21.25],
		         [15.25,15.25,15.25,15.25,15.25,15.25,15.25,15.25],
		         [09.25,09.25,09.25,09.25,09.25,09.25,09.25,09.25],
		         [03.25,03.25,03.25,03.25,03.25,03.25,03.25,03.25],
		         [-03.25,-03.25,-03.25,-03.25,-03.25,-03.25,-03.25,-03.25],
		         [-09.25,-09.25,-09.25,-09.25,-09.25,-09.25,-09.25,-09.25],
		         [-15.25,-15.25,-15.25,-15.25,-15.25,-15.25,-15.25,-15.25],
    		         [-21.25,-21.25,-21.25,-21.25,-21.25,-21.25,-21.25,-21.25]]])

class BField():
    def __init__(self,year=None,doy=None,site=1,heights=None,alpha_i=90):
        """
        BField class creates an object to get  the Magnetic  field for a  specific date and
        height(s).

        Parameters
        ----------
        year = A scalar giving the desired year.  If the value is None (default value) then
          the current year will be used.
        doy = A scalar giving the desired day of the year. If the value is None (default va-
          lue) then the current doy will be used.
        site = An integer to choose the geographic coordinates of the place where the magne-
          tic field will be computed. The default value is over Jicamarca (site=1)
        heights = An array giving the heights (km) where the magnetic field will be modeled By default the magnetic field will be computed at 100, 500 and 1000km.
        alpha_i = Angle to interpolate the magnetic field.

        Modification History
        --------------------
        Converted to Object-oriented Programming by Freddy Galindo, ROJ, 07 October 2009.
        """

        tmp = time.localtime()
        if year==None: year = tmp[0]
        if doy==None:  doy = tmp[7]
        self.year = year
        self.doy = doy
        self.site = site
        if heights is None:
            heights = numpy.array([100,500,1000])
        else:
            heights = numpy.array(heights)
        self.heights = heights
        self.alpha_i = alpha_i
        

    def getBField(self,maglimits=numpy.array([-7,-7,7,7])):
        """
        getBField models the magnetic field for a different heights in a specific date.

        Parameters
        ----------
        maglimits = An 4-elements array giving ..... The default value is [-7,-7,7,7].

        Return
        ------
        dcos = An 4-dimensional array giving the directional cosines of the magnetic field
          over the desired place.
        alpha = An 3-dimensional array giving the angle of the magnetic field over the desi-
          red place.

        Modification History
        --------------------
        Converted to Python by Freddy R. Galindo, ROJ, 07 October 2009.
        Updated to include AMISR by Joab Apaza, ROJ, July 2023.
        """

        x_ant = numpy.array([1,0,0])
        y_ant = numpy.array([0,1,0])
        z_ant = numpy.array([0,0,1])
        
        if self.site==0:
            title_site = "Magnetic equator"
            coord_site = numpy.array([-76+52./60.,-11+57/60.,0.5])
            grid_res = 0.5
        elif self.site==1:
            title_site = 'Jicamarca'
            coord_site = [-76-52./60.,-11-57/60.,0.5]
            theta = (45+5.35)*numpy.pi/180.      # (50.35 and 1.46 from Fleish Thesis)
            delta = -1.46*numpy.pi/180
            grid_res = 0.5
        elif self.site==2:
            title_site = 'AMISR 14'
            coord_site = [-76.874913, -11.953371,0.52984]
            theta = (0.0977)*numpy.pi/180.       #  0.0977
            delta = 0.110*numpy.pi/180        # 0.11
            maglimits=numpy.array([-45,-45,45,45])
            grid_res = 2.5
        else:
#            print "No defined Site. Skip..."
            return None
        
        x_ant1 = numpy.roll(self.rotvector(self.rotvector(x_ant,1,delta),3,theta),1)
        y_ant1 = numpy.roll(self.rotvector(self.rotvector(y_ant,1,delta),3,theta),1)
        z_ant1 = numpy.roll(self.rotvector(self.rotvector(z_ant,1,delta),3,theta),1)

        ang0 = -1*coord_site[0]*numpy.pi/180.
        ang1 = coord_site[1]*numpy.pi/180.
        x_ant = self.rotvector(self.rotvector(x_ant1,2,ang1),3,ang0)
        y_ant = self.rotvector(self.rotvector(y_ant1,2,ang1),3,ang0)
        z_ant = self.rotvector(self.rotvector(z_ant1,2,ang1),3,ang0)
        
        nhei = self.heights.size
        pt_intercep = numpy.zeros((nhei,2))
        nfields = 1

        
        nlon = int(numpy.int(maglimits[2] - maglimits[0])/grid_res + 1)
        nlat = int(numpy.int(maglimits[3] - maglimits[1])/grid_res + 1)

        location = numpy.zeros((nlon,nlat,2))
        mlon = numpy.atleast_2d(numpy.arange(nlon)*grid_res + maglimits[0])
        mrep = numpy.atleast_2d(numpy.zeros(nlat) + 1)
        location0 = numpy.dot(mlon.transpose(),mrep)

        mlat = numpy.atleast_2d(numpy.arange(nlat)*grid_res + maglimits[1])
        mrep = numpy.atleast_2d(numpy.zeros(nlon) + 1)
        location1 = numpy.dot(mrep.transpose(),mlat)

        location[:,:,0] = location0
        location[:,:,1] = location1

        alpha = numpy.zeros((nlon,nlat,nhei))
        rr = numpy.zeros((nlon,nlat,nhei,3))
        dcos = numpy.zeros((nlon,nlat,nhei,2))

        global first_time

        first_time = None
        for ilon in numpy.arange(nlon):
            for ilat in numpy.arange(nlat):
                outs = self.__bdotk(self.heights,
                                    self.year + self.doy/366.,
                                    coord_site[1],
                                    coord_site[0],
                                    coord_site[2],
                                    coord_site[1]+location[ilon,ilat,1],
                                    location[ilon,ilat,0]*720./180.)

                alpha[ilon, ilat,:] = outs[1]
                rr[ilon, ilat,:,:] = outs[3]

                mrep = numpy.atleast_2d((numpy.zeros(nhei)+1)).transpose()
                tmp = outs[3]*numpy.dot(mrep,numpy.atleast_2d(x_ant))
                tmp = tmp.sum(axis=1)
                dcos[ilon,ilat,:,0] = tmp/numpy.sqrt((outs[3]**2).sum(axis=1))

                mrep = numpy.atleast_2d((numpy.zeros(nhei)+1)).transpose()
                tmp = outs[3]*numpy.dot(mrep,numpy.atleast_2d(y_ant))
                tmp = tmp.sum(axis=1)
                dcos[ilon,ilat,:,1] = tmp/numpy.sqrt((outs[3]**2).sum(axis=1))

        return dcos, alpha, nlon, nlat


    def __bdotk(self,heights,tm,gdlat=-11.95,gdlon=-76.8667,gdalt=0.0,decd=-12.88, ham=-4.61666667):

        global first_time
        # Mean Earth radius in Km WGS 84
        a_igrf = 6371.2

        bk = numpy.zeros(heights.size)
        alpha = numpy.zeros(heights.size)
        bfm = numpy.zeros(heights.size)
        rr = numpy.zeros((heights.size,3))
        rgc = numpy.zeros((heights.size,3))

        ObjGeodetic = Astro_Coords.Geodetic(gdlat,gdalt)
        [gclat,gcalt] = ObjGeodetic.change2geocentric()

        gclat = gclat*numpy.pi/180.
        gclon = gdlon*numpy.pi/180.

        # Antenna position from center of Earth
        ca_vector = [numpy.cos(gclat)*numpy.cos(gclon),numpy.cos(gclat)*numpy.sin(gclon),numpy.sin(gclat)]
        ca_vector = gcalt*numpy.array(ca_vector)

        dec = decd*numpy.pi/180.

        # K  vector respect to the center of earth.
        klon = gclon + ham*numpy.pi/720.
        k_vector = [numpy.cos(dec)*numpy.cos(klon),numpy.cos(dec)*numpy.sin(klon),numpy.sin(dec)]
        k_vector = numpy.array(k_vector)
        
        for ih in numpy.arange(heights.size):
            # Vector from Earth's center to volume of interest
            rr[ih,:] = k_vector*heights[ih]
            cv_vector = numpy.squeeze(ca_vector) + rr[ih,:]

            cv_gcalt = numpy.sqrt(numpy.sum(cv_vector**2.))
            cvxy = numpy.sqrt(numpy.sum(cv_vector[0:2]**2.))

            radial = cv_vector/cv_gcalt
            east = numpy.array([-1*cv_vector[1],cv_vector[0],0])/cvxy
            comp1 = east[1]*radial[2] - radial[1]*east[2]
            comp2 = east[2]*radial[0] - radial[2]*east[0]
            comp3 = east[0]*radial[1] - radial[0]*east[1]
            north = -1*numpy.array([comp1, comp2, comp3])

            rr_k = cv_vector - numpy.squeeze(ca_vector)
            u_rr = rr_k/numpy.sqrt(numpy.sum(rr_k**2.))

            cv_gclat = numpy.arctan2(cv_vector[2],cvxy)
            cv_gclon = numpy.arctan2(cv_vector[1],cv_vector[0])

            bhei = cv_gcalt-a_igrf
            blat = cv_gclat*180./numpy.pi
            blon = cv_gclon*180./numpy.pi
            bfield = self.__igrfkudeki(bhei,tm,blat,blon)

            B = (bfield[0]*north + bfield[1]*east - bfield[2]*radial)*1.0e-5

            bfm[ih] = numpy.sqrt(numpy.sum(B**2.)) #module
            bk[ih] = numpy.sum(u_rr*B)
            alpha[ih] = numpy.arccos(bk[ih]/bfm[ih])*180/numpy.pi
            rgc[ih,:] = numpy.array([cv_gclon, cv_gclat, cv_gcalt])

        return bk, alpha, bfm, rr, rgc


    def __igrfkudeki(self,heights,time,latitude,longitude,ae=6371.2):
        """
        __igrfkudeki calculates the International Geomagnetic Reference Field for given in-
        put conditions based on IGRF2005 coefficients.

        Parameters
        ----------
        heights = Scalar or vector giving the height above the Earth  of the point in ques-
          tion in kilometers.
        time = Scalar or vector giving the decimal year of time in question (e.g. 1991.2).
        latitude = Latitude of point in question in decimal degrees. Scalar or vector.
        longitude = Longitude of point in question in decimal degrees. Scalar or vector.
        ae =
        first_time =

        Return
        ------
        bn =
        be =
        bd =
        bmod =
        balpha =
        first_time =

        Modification History
        --------------------
        Converted to Python by Freddy R. Galindo, ROJ, 03 October 2009.
        """

        global first_time
        global gs, hs, nvec, mvec, maxcoef

        heights = numpy.atleast_1d(heights)
        time = numpy.atleast_1d(time)
        latitude = numpy.atleast_1d(latitude)
        longitude = numpy.atleast_1d(longitude)

        if numpy.max(latitude)==90:
#            print "Field calculations are not supported at geographic poles"
            pass

        # output arrays
        bn = numpy.zeros(heights.size)
        be = numpy.zeros(heights.size)
        bd = numpy.zeros(heights.size)

        if first_time==None:first_time=0

        time0 = time[0]
        if time!=first_time:
            #print "Getting coefficients for", time0
            [periods,g,h ] = self.__readIGRFcoeff()
            top_year = numpy.max(periods)
            nperiod = (top_year - 1900)/5 + 1

            maxcoef = 10
            if time0>=2000:maxcoef = 12


            # Normalization array for Schmidt fucntions
            multer = numpy.zeros((2+maxcoef,1+maxcoef)) + 1
            for cn in (numpy.arange(maxcoef)+1):
                for rm in (numpy.arange(cn)+1):
                    tmp = numpy.arange(2*rm) + cn - rm + 1.
                    multer[rm+1,cn] = ((-1.)**rm)*numpy.sqrt(2./tmp.prod())

            schmidt = multer[1:,1:].transpose()

            # n and m arrays
            nvec = numpy.atleast_2d(numpy.arange(maxcoef)+2)
            mvec = numpy.atleast_2d(numpy.arange(maxcoef+1)).transpose()

            # Time adjusted igrf g and h with Schmidt normalization
            # IGRF coefficient arrays: g0(n,m), n=1, maxcoeff,m=0, maxcoeff, ...
            if time0<top_year:
                dtime = (time0 - 1900) % 5
                ntime = int((time0 - 1900 - dtime)/5)
            else:
                # Estimating coefficients for times > top_year
                dtime = (time0 - top_year) + 5
                ntime = int(g[:,0,0].size - 2)

            
            g0 = g[ntime,1:maxcoef+1,:maxcoef+1]
            h0 = h[ntime,1:maxcoef+1,:maxcoef+1]
            gdot = g[ntime+1,1:maxcoef+1,:maxcoef+1]-g[ntime,1:maxcoef+1,:maxcoef+1]
            hdot = h[ntime+1,1:maxcoef+1,:maxcoef+1]-h[ntime,1:maxcoef+1,:maxcoef+1]
            gs = (g0 + dtime*(gdot/5.))*schmidt[:maxcoef,0:maxcoef+1]
            hs = (h0 + dtime*(hdot/5.))*schmidt[:maxcoef,0:maxcoef+1]

            first_time = time0

        for ii in numpy.arange(heights.size):
            # Height dependence array rad = (ae/(ae+height))**(n+3)
            rad = numpy.atleast_2d((ae/(ae + heights[ii]))**(nvec+1))

            # Sin and Cos of m times longitude phi arrays
            mphi = mvec*longitude[ii]*numpy.pi/180.
            cosmphi = numpy.atleast_2d(numpy.cos(mphi))
            sinmphi = numpy.atleast_2d(numpy.sin(mphi))

            # Cos of colatitude theta
            c = numpy.cos((90 - latitude[ii])*numpy.pi/180.)

            # Legendre functions p(n,m|c)
            [p,dp]= scipy.special.lpmn(maxcoef+1,maxcoef+1,c)
            p = p[:,:-1].transpose()
            s = numpy.sqrt((1. - c)*(1 + c))

            # Generate derivative array dpdtheta = -s*dpdc
            dpdtheta = c*p/s
            for m in numpy.arange(maxcoef+2):    dpdtheta[:,m] = m*dpdtheta[:,m]
            dpdtheta = dpdtheta + numpy.roll(p,-1,axis=1)

            # Extracting arrays required for field calculations
            p = p[1:maxcoef+1,:maxcoef+1]
            dpdtheta = dpdtheta[1:maxcoef+1,:maxcoef+1]

            # Weigh p and dpdtheta with gs and hs coefficients.
            gp = gs*p
            hp = hs*p
            gdpdtheta = gs*dpdtheta
            hdpdtheta = hs*dpdtheta
            # Calcultate field components
            matrix0 = numpy.dot(gdpdtheta,cosmphi)
            matrix1 = numpy.dot(hdpdtheta,sinmphi)
            bn[ii] = numpy.dot(rad,(matrix0 + matrix1))
            matrix0 = numpy.dot(hp,(mvec*cosmphi))
            matrix1 = numpy.dot(gp,(mvec*sinmphi))
            be[ii] = numpy.dot((-1*rad),((matrix0 - matrix1)/s))
            matrix0 = numpy.dot(gp,cosmphi)
            matrix1 = numpy.dot(hp,sinmphi)
            bd[ii] = numpy.dot((-1*nvec*rad),(matrix0 + matrix1))

        bmod = numpy.sqrt(bn**2. + be**2. + bd**2.)
        btheta = numpy.arctan(bd/numpy.sqrt(be**2. + bn**2.))*180/numpy.pi
        balpha = numpy.arctan(be/bn)*180./numpy.pi

        #bn : north
        #be : east
        #bn : radial
        #bmod : module


        return bn, be, bd, bmod, btheta, balpha

    def str2num(self, datum):
        try:
            return int(datum)
        except:
            try:
                return float(datum)
            except:
                return datum

    def __readIGRFfile(self, filename):
        list_years=[]
        for i in range(1,26):
            list_years.append(1895.0 + i*5)

        epochs=list_years
        epochs.append(epochs[-1]+5)
        nepochs = numpy.shape(epochs)

        gg = numpy.zeros((13,14,nepochs[0]),dtype=float)
        hh = numpy.zeros((13,14,nepochs[0]),dtype=float)

        coeffs_file=open(filename)
        lines=coeffs_file.readlines()

        coeffs_file.close()

        for line in lines:
            items = line.split()
            g_h = items[0]
            n = self.str2num(items[1])
            m = self.str2num(items[2])

            coeffs = items[3:]

            for i in range(len(coeffs)-1):
                coeffs[i] = self.str2num(coeffs[i])

            #coeffs = numpy.array(coeffs)
            ncoeffs = numpy.shape(coeffs)[0]

            if g_h == 'g':
    #            print n," g ",m
                gg[n-1,m,:]=coeffs
            elif g_h=='h':
    #            print n," h ",m
                hh[n-1,m,:]=coeffs
    #        else :
    #            continue

    #    Ultimo Reordenamiento para  almacenar .
        gg[:,:,nepochs[0]-1] = gg[:,:,nepochs[0]-2] + 5*gg[:,:,nepochs[0]-1]
        hh[:,:,nepochs[0]-1] = hh[:,:,nepochs[0]-2] + 5*hh[:,:,nepochs[0]-1]

#        return numpy.array([gg,hh])
        periods = numpy.array(epochs)
        g = gg
        h = hh
        return periods, g, h


    def __readIGRFcoeff(self,filename="igrf10coeffs.dat"):
        """
        __readIGRFcoeff reads the coefficients from  a binary file which is located  in the
        folder "resource."

        Parameter
        ---------
        filename = A string to specify the name of the file which contains thec coeffs. The
        default value is "igrf10coeffs.dat"

        Return
        ------
        periods = A lineal array giving...
        g1 =
        h1 =

        Modification History
        --------------------
        Converted to Python by Freddy R. Galindo, ROJ, 03 October 2009.
        """

        base_path = os.path.dirname(os.path.abspath(__file__))
        filename = os.path.join(base_path, "igrf13coeffs.txt")

        period_v, g_v, h_v = self.__readIGRFfile(filename)
        g2 = numpy.zeros((14,14,26))
        h2 = numpy.zeros((14,14,26))
        g2[1:14,:,:] = g_v
        h2[1:14,:,:] = h_v

        g = numpy.transpose(g2, (2,0,1))
        h = numpy.transpose(h2, (2,0,1))
        periods = period_v.copy()

        return periods, g, h

    def rotvector(self,vector,axis=1,ang=0):
        """
        rotvector function returns the new vector generated rotating the rectagular coords.

        Parameters
        ----------
        vector = A lineal 3-elements array (x,y,z).
        axis = A integer to specify  the axis used to rotate the coord systems. The default
          value is 1.
            axis = 1 -> Around "x"
            axis = 2 -> Around "y"
            axis = 3 -> Around "z"
        ang = Angle of rotation (in radians). The default value is zero.

        Return
        ------
        rotvector = A lineal array of 3 elements giving the new coordinates.

        Modification History
        --------------------
        Converted to Python by Freddy R. Galindo, ROJ, 01 October 2009.
        """

        if axis==1:
            t = [[1,0,0],[0,numpy.cos(ang),numpy.sin(ang)],[0,-numpy.sin(ang),numpy.cos(ang)]]
        elif axis==2:
            t = [[numpy.cos(ang),0,-numpy.sin(ang)],[0,1,0],[numpy.sin(ang),0,numpy.cos(ang)]]
        elif axis==3:
            t = [[numpy.cos(ang),numpy.sin(ang),0],[-numpy.sin(ang),numpy.cos(ang),0],[0,0,1]]

        rotvector = numpy.array(numpy.dot(numpy.array(t),numpy.array(vector)))

        return rotvector


class AntPatternPlot:

    def __init__(self,ploteo=0):
        """
        AntPatternPlot creates an object to call methods to plot the antenna pattern.

        Modification History
        --------------------
        Created by Freddy Galindo, ROJ, 06 October 2009.
        """
        ######################## CAMBIO
        #self.fig = Figure(figsize=(8,8), facecolor='white')
        #self.ax = self.fig.add_subplot(111)        
        ########################
        if ploteo == 1:
            self.fig,self.ax = plt.subplots()
            
        else:
            self.fig = Figure(figsize=(8,8), facecolor='white')
            self.ax  = self.fig.add_subplot(111)

    def contPattern(self,site=1, iplot=0,gpath='',filename='',mesg='',amp=None ,x=None ,y=None ,getCut=None,title='', save=False):
        """
        contPattern plots a contour map of the antenna pattern.

        Parameters
        ----------
        iplot = A integer to specify if the plot is  the first, second, ...  The default va-
          lue is 0.

        Examples
        --------
        >> Over_Jro.JroPattern(pattern=2).contPattern()

        Modification history
        --------------------
        Converted to Python by Freddy R. Galindo, ROJ, 06 October 2009.
        Updated to include AMISR by Joab Apaza, ROJ, July 2023.
        """

        if getCut == 1:
            return

        xmax = numpy.max(x)
        xmin = numpy.min(x)
        ymax = numpy.max(y)
        ymin = numpy.min(y)

        levels = numpy.array([1e-3,1e-2,1e-1,0.5,1.0])
        tmp = numpy.round(10*numpy.log10(levels),decimals=1)
        labels = list(range(5))
        for i in numpy.arange(5):
            labels[i] = str(numpy.int(tmp[i]))


        colors = ((0,0,1.),(0,170/255.,0),(127/255.,1.,0),(1.,109/255.,0),(128/255.,0,0))
        if site== 1:
            CS = self.ax.contour(x,y,amp.transpose(),levels,colors=colors)
        else:
            CS = self.ax.contour(x,y,amp,levels,colors=colors)
        fmt = {}
        for l,s in zip(CS.levels,labels):
            fmt[l] = s

        
        self.ax.annotate(mesg,xy=(0,0),xytext=(0.01,0.01),xycoords='figure fraction')
        self.ax.clabel(CS,CS.levels,inline=True,fmt=fmt,fontsize=10)
        self.ax.set_xlim(xmin,xmax)
        self.ax.set_ylim(ymin,ymax)
        self.ax.set_title("Total Pattern: " + title)
        if site==1:
            self.ax.annotate('Ng',xy=(-0.05,1.04),xytext=(0.01,0.962),xycoords='axes fraction',arrowprops=dict(facecolor='black', width=1.,shrink=0.2),fontsize=15.)
            self.ax.set_xlabel("West to South")
            self.ax.set_ylabel("West to North")
        else:
            self.ax.set_xlabel("West  to  East")
            self.ax.set_ylabel("South  to  North")
        self.ax.grid(True)      


    def plotRaDec(self,site=1, gpath=None,filename=None,jd=2452640.5,ra_obs=None,xg=None,yg=None,x=None,y=None,allAmisr_x=[], allAmisr_y=[], save=True):
        """
        plotRaDec draws right ascension and declination lines on a JRO plane. This function
        must call after conPattern.

        Parameters
        ----------
        jd = A scalar giving the Julian date.
        ra_obs = Scalar giving the right ascension of the observatory.
        xg = A 3-element array to specify ..
        yg = A 3-element array to specify ..

        Examples
        --------
        >> Over_Jro.JroPattern(pattern=2).contPattern()
        >> Over_Jro.JroPattern(pattern=2).plotRaDec(jd=jd,ra_obs=ra_obs,xg=xg,yg=yg)

        Modification history
        --------------------
        Converted to Python by Freddy R. Galindo, ROJ, 06 October 2009.
        """

        # Finding RA of observatory for a specific date
        if ra_obs is None:ra_obs = numpy.array([23.37060849])
        if xg is None:xg = numpy.array([0.62918474,-0.77725579,0.])
        if yg is None:yg = numpy.array([0.77700346,0.62898048,0.02547905])

        # Getting HA and DEC axes
        if site==1:
            mindec = -28; maxdec = 4; incdec = 2.
            minha = -20; maxha = 20; incha = 2.
        else:
            mindec = -45; maxdec = 25; incdec = 4
            minha = -80; maxha = 80; incha = 6

        ndec = numpy.int((maxdec - mindec)/incdec) + 1

        
        nha = numpy.int((maxha - minha)/incha) + 1

        #mcosx = numpy.zeros((nha,ndec))
        #mcosy = numpy.zeros((nha,ndec))

        ha_axes = numpy.reshape(numpy.arange(nha)*incha + minha,(nha,1))
        ones_dec = numpy.reshape(numpy.zeros(ndec) + 1,(ndec,1))
        ha_axes = numpy.dot(ha_axes,ones_dec.transpose())
        ha_axes2 = numpy.array(ra_obs - ha_axes)

        dec_axes = numpy.reshape(numpy.arange(ndec)*incdec + mindec,(ndec,1))
        ones_ra = numpy.reshape(numpy.zeros(nha) + 1,(nha,1))
        dec_axes = numpy.dot(ones_ra,dec_axes.transpose())
        dec_axes2 = numpy.array(dec_axes)

        ObjHor = Astro_Coords.Equatorial(ha_axes2,dec_axes2,jd)
        [alt,az,ha] = ObjHor.change2AltAz()

        z = numpy.transpose(alt)*Misc_Routines.CoFactors.d2r  ; z = z.flatten()
        az = numpy.transpose(az)*Misc_Routines.CoFactors.d2r  ; az = az.flatten()

        vect = numpy.array([numpy.cos(z)*numpy.sin(az),numpy.cos(z)*numpy.cos(az),numpy.sin(z)])

        xg = numpy.atleast_2d(xg)
        dcosx = numpy.array(numpy.dot(xg,vect))
        yg = numpy.atleast_2d(yg)
        dcosy = numpy.array(numpy.dot(yg,vect))

        mcosx = dcosx.reshape(ndec,nha)
        mcosy = dcosy.reshape(ndec,nha)

        # Defining NAN for points outof limits.
        xmax = numpy.max(x)
        xmin = numpy.min(x)
        ymax = numpy.max(y)
        ymin = numpy.min(y)

        factor = 1.3

        noval = numpy.where((mcosx>(xmax*factor)) | (mcosx<(xmin*factor)))
        if noval[0].size>0:mcosx[noval] = numpy.nan
        noval = numpy.where((mcosy>(ymax*factor)) | (mcosy<(ymin*factor)))
        if noval[0].size>0:mcosy[noval] = numpy.nan

        # Plotting HA and declination grid.
        iha0 = numpy.int((0 - minha)/incha)
        idec0 = numpy.int((-14 - mindec)/incdec)

        colorgrid = (1.,109/255.,0)
        self.ax.plot(mcosx.transpose(),mcosy.transpose(),color=colorgrid,linestyle='--', lw=0.5)
        for idec in numpy.arange(ndec):
            if idec != idec0:
                valx = (mcosx[idec,iha0]<=xmax) & (mcosx[idec,iha0]>=xmin)
                valy = (mcosy[idec,iha0]<=ymax) & (mcosy[idec,iha0]>=ymin)
                if valx & valy:
                    text = str(numpy.int(mindec + incdec*idec))+'$^o$'
                    self.ax.text(mcosx[idec,iha0],mcosy[idec,iha0],text)

        self.ax.plot(mcosx,mcosy,color=colorgrid,linestyle='--',lw=0.5)
        for iha in numpy.arange(nha):
            if iha != iha0:
                valx = (mcosx[idec0,iha]<=xmax) & (mcosx[idec0,iha]>=xmin)
                valy = (mcosy[idec0,iha]<=ymax) & (mcosy[idec0,iha]>=ymin)
                if valx & valy:
                    text = str(numpy.int(minha + incha*iha))+"'"
                    self.ax.text(mcosx[idec0,iha],mcosy[idec0,iha],text)
        
        if len(allAmisr_x) > 0 and len(allAmisr_y) > 0:
            self.ax.scatter(allAmisr_x, allAmisr_y, c='yellow', marker='x', s = 40)
        if save:
            save_fig = os.path.join(gpath,filename)
            self.fig.savefig(save_fig,format='png')
    
    def ploteo(self):
        plt.show()



    def plotBField(self,gpath,filename,dcos,alpha, nlon, nlat, dcosxrange, dcosyrange, heights, alpha_i, plot=True):
        """
        plotBField draws the magnetic field in a directional cosines plot.

        Parameters
        ----------
        dcos = An 4-dimensional array giving the directional cosines of the magnetic field
          over the desired place.
        alpha = An 3-dimensional array giving the angle of the magnetic field over the desi-
          red place.
        nlon = An integer to specify the number of elements per longitude.
        nlat = An integer to specify the number of elements per latitude.
        dcosxrange = A 2-element array giving the range of the  directional cosines  in the
          "x" axis.
        dcosyrange = A 2-element array giving the range of the  directional cosines  in the
          "y" axis.
        heights = An array giving the heights (km) where the magnetic field will be modeled               By default the magnetic field will be computed at 100, 500 and 1000km.
        alpha_i = Angle to interpolate the magnetic field.
        Modification History
        --------------------
        Converted to Python by Freddy R. Galindo, ROJ, 07 October 2009.
        """

        handles = []
        objects = []
        colors = ['k','m','c','b','g','r','y']
        marker = ['-+','-*','-D','-x','-s','->','-o','-^']

        alpha_location = numpy.zeros((nlon,2,heights.size))

        for ih in numpy.arange(heights.size):
            alpha_location[:,0,ih] = dcos[:,0,ih,0]
            for ilon in numpy.arange(nlon):
                myx = (alpha[ilon,:,ih])[::-1]
                myy = (dcos[ilon,:,ih,0])[::-1]
                tck = scipy.interpolate.splrep(myx,myy,s=0)
                mydcosx = scipy.interpolate.splev(alpha_i,tck,der=0)

                myx = (alpha[ilon,:,ih])[::-1]
                myy = (dcos[ilon,:,ih,1])[::-1]
                tck = scipy.interpolate.splrep(myx,myy,s=0)
                mydcosy = scipy.interpolate.splev(alpha_i,tck,der=0)
                alpha_location[ilon,:,ih] = numpy.array([mydcosx, mydcosy])

            if plot:
                ObjFig, = self.ax.plot(alpha_location[:,0,ih],alpha_location[:,1,ih],
                    marker[ih % 8],color=colors[numpy.int(ih/8)],ms=5.5,lw=0.5)
                handles.append(ObjFig)
                objects.append(numpy.str(heights[ih]) + ' km')        

        if plot:
            legend = self.ax.legend(handles, objects,loc="lower right", numpoints=1, handlelength=0.3,
                        handletextpad=0.02, borderpad=0.3, labelspacing=0.1)
            self.ax.add_artist(legend)

        self.alpha_location = alpha_location

    def plotCelestial(self, jd, main_dec, tod, maxha_min, objects, glat, glon, xg, yg, dcosxrange, dcosyrange, plot=True):

        self.tod = tod
                
        self.dcosx_sun = 1
        self.dcosy_sun = 1
        self.ha_sun = 1
        self.time_sun = 1
        
        self.dcosx_moon = 1
        self.dcosy_moon = 1
        self.ha_moon = 1
        self.time_moon = 1
        
        self.dcosx_hydra = 1
        self.dcosy_hydra = 1
        self.ha_hydra = 1
        self.time_hydra = 1
        
        self.dcosx_galaxy = 1
        self.dcosy_galaxy = 1
        self.ha_galaxy = 1
        self.time_galaxy = 1
        
        tod = self.tod
        
        maxlev = 24; minlev = 0; maxcol = 39; mincol = 10        
        handles = []
        titles = ['','$Sun$','$Moon$','$Hydra$','$Galaxy$']
        marker = ['','--^','--s','--*','--o']
        
        # Getting RGB table to plot celestial object over Jicamarca
        colortable = ['', 'olive', 'indigo', 'cyan', 'red']
        labels = [] 
        for io in objects:
            ObjBodies = Astro_Coords.CelestialBodies()
            if io==0:
                continue
            elif io==1:
                [ra,dec,sunlon,sunobliq] = ObjBodies.sunpos(jd)
            elif io==2:                    
                [ra,dec,dist,moonlon,moonlat] = ObjBodies.moonpos(jd)
            elif io==3:                    
                [ra,dec] = ObjBodies.hydrapos()
            elif io==4:
                [maxra,ra] = ObjBodies.skynoise_jro(dec_cut=main_dec)
                ra = maxra*15.
                dec = main_dec

            ObjEq = Astro_Coords.Equatorial(ra, dec, jd, lat=glat, lon=glon)
            [alt, az, ha] = ObjEq.change2AltAz()
            vect = numpy.array([az,alt]).transpose()
            vect = Misc_Routines.Vector(vect,direction=0).Polar2Rect()
            
            dcosx = numpy.array(numpy.dot(vect,xg))
            dcosy = numpy.array(numpy.dot(vect,yg))
            wrap = numpy.where(ha>=180.)
            
            if wrap[0].size>0:
                ha[wrap] = ha[wrap] - 360.
            
            val = numpy.where((numpy.abs(ha))<=(maxha_min*0.25))
            
            if val[0].size>2:
                tod_1 = tod*1.
                shift_1 = numpy.where(tod>12.)
                tod_1[shift_1] = tod_1[shift_1] - 24.
                tod_2 = tod*1.
                shift_2 = numpy.where(tod<12.)
                tod_2[shift_2] = tod_2[shift_2] + 24.
                
                diff0 = numpy.nanmax(tod[val])  -  numpy.nanmin(tod[val])
                diff1 = numpy.nanmax(tod_1[val]) - numpy.nanmin(tod_1[val])
                diff2 = numpy.nanmax(tod_2[val]) - numpy.nanmin(tod_2[val])
                
                if ((diff0<=diff1) & (diff0<=diff2)):
                    tod_0 = tod
                elif ((diff1<diff0) & (diff1<diff2)):
                    tod_0 = tod_1
                else: 
                    tod_0 = tod_2

                if io==1:
                    self.dcosx_sun = dcosx[val]
                    self.dcosy_sun = dcosy[val]    
                    self.ha_sun = ha[val]
                    self.time_sun = numpy.median(tod_0[val])
                elif io==2:                    
                    self.dcosx_moon = dcosx[val]
                    self.dcosy_moon = dcosy[val]    
                    self.ha_moon = ha[val]
                    self.time_moon = numpy.median(tod_0[val])
                elif io==3:    
                    self.dcosx_hydra = dcosx[val]
                    self.dcosy_hydra = dcosy[val]    
                    self.ha_hydra = ha[val]
                    self.time_hydra = numpy.mean(tod_0[val])
                elif io==4:
                    self.dcosx_galaxy = dcosx[val]
                    self.dcosy_galaxy = dcosy[val]    
                    self.ha_galaxy = ha[val]
                    self.time_galaxy = numpy.mean(tod_0[val])

                #index = numpy.mean(tod_0[val]) - minlev
                #index = (index*(maxcol - mincol)/(maxlev - minlev)) + mincol
                #index = numpy.int(index)
                if plot:
                    figobjects, = self.ax.plot(dcosx[val], dcosy[val], marker[io], lw=1, ms=7, mew=0, color=colortable[io])    
                    handles.append(figobjects)
                    labels.append(titles[io])
         
        if plot:
            legend = self.ax.legend(handles,labels,loc="lower left", numpoints=1, handlelength=0.5, \
                                    borderpad=0.3, handletextpad=0.1,labelspacing=0.2,fontsize='small')
            
            self.ax.add_artist(legend)
        
class PatternCutPlot:

    def __init__(self, nsubplots):
        self.nsubplots = nsubplots

        self.fig = None

        self.__plot_width  = 8

        if self.nsubplots == 5:
            self.__plot_height = 11

        if self.nsubplots == 4:
            self.__plot_height = 9

        if self.nsubplots == 3:
            self.__plot_height = 7

        if self.nsubplots == 2:
            self.__plot_height = 5

        if self.nsubplots == 1:
            self.__plot_height = 3

        self.fig = Figure(figsize=(self.__plot_width, self.__plot_height), facecolor='white')

        if self.nsubplots < 5:
            self.__height_inch = 1.1 #altura de los subplots (pulgadas)
            top_inch = 1.5/2.7 #espacio entre el primer subplot y el limite superior del plot
            self.__vspace_plot_inch = 1.0#1.5/2 # espacio vertical entre subplots
            self.__left = 0.1
        else:
            self.__height_inch = 1.1 #altura de los subplots (pulgadas)
            top_inch = 1.5/2.7 #espacio entre el primer subplot y el limite superior del plot
            self.__vspace_plot_inch = 1.0 # espacio vertical entre subplots
            self.__left = 0.1

        self.__bottom_inch = self.__plot_height - (self.__height_inch + top_inch)
        self.__height = self.__height_inch/self.__plot_height

        self.__width = 0.8


    def drawCut(self,io,patterns,npatterns,ha,otitle,subtitle,ptitle):

        t_cuts = ['B','Sun','Moon','Hydra','Galaxy']
        self.__bottom = self.__bottom_inch/self.__plot_height


        subp = self.fig.add_axes([self.__left,self.__bottom,self.__width,self.__height])

        on_axis_angle = -4.65562
        for icut in numpy.arange(npatterns):
            # Getting Antenna cut.
            pattern = patterns[icut]
            power = numpy.abs(pattern/numpy.nanmax(pattern))
            max_power_db = numpy.round(10.*numpy.log10(numpy.nanmax(pattern)),2)

            bval = numpy.where(power[:,0]==numpy.nanmax(power))
            beta = -0.25*(ha[bval[0]] + on_axis_angle)
#            print 'Angle (deg): '+"%f"%(beta)

            subp.plot(ha,power)


        xmax = numpy.max(numpy.nanmin(ha))
        xmin = numpy.min(numpy.nanmax(ha))
        ymax = numpy.max(1)
        ymin = numpy.min(0)


        subp.set_xlim(xmin, xmax)

        subp.set_ylim(ymin, ymax)

        subp.set_title(otitle + ' ' + ptitle,size="medium")

        subp.text(0.5, 1.26,subtitle[0],
         horizontalalignment='center',
         verticalalignment='center',
         transform = subp.transAxes)

        xlabels = subp.get_xticks()

        subp.set_xticklabels(xlabels,size="small")

        ylabels = subp.get_yticks()

        subp.set_yticklabels(ylabels,size="small")

        subp.set_xlabel('Hour angle (min) (+ve to West)',size="small")

        subp.set_ylabel("Power [Max: " + str(max_power_db) + ' dB]',size="small")

        subp.grid()


        self.__bottom_inch = self.__bottom_inch - (self.__height_inch + self.__vspace_plot_inch)


class JroPattern():
    def __init__(self,pattern=0,path=None,filename=None,nptsx=101,nptsy=101,maxphi=5,fftopt=0, \
        getcut=0,dcosx=[],dcosy=[],eomwl=6,airwl=4, **kwargs):
        """
        JroPattern class creates an object to represent the useful parameters for beam mode-
        lling of the Jicamarca VHF radar.

        Parameters
        ----------
        pattern = An integer (See JroAntSetup to know the available values) to load a prede-
          fined configuration. The default value  is 0. To use a  user-defined configuration
          pattern  must be None.
        path = A string giving the directory that contains the user-configuration file. PATH
          will work if pattern is None.
        filename = A string giving the  name of the  user-configuration  file. FILENAME will
          work if pattern is None.
        nptsx = A scalar to specify the number of points  used to define the angular resolu-
          tion in the "x" axis. The default value is 101.
        nptsy = A scalar to specify the number of points  used to define the angular resolu-
          tion in the "x" axis. The default value is 101.
        maxphi = A scalar giving the maximum (absolute) angle (in degree) to model the ante-
          nna pattern. The default value is 5 degrees.
        fftopt = Set this input  to 1 to model  the beam using FFT.  To model  using antenna
          theory set to 0 (default value).
        getcut = Set to 1 to show an antenna cut instead of a contour plot of itself (set to
          0). The defautl value is 0.
        dcosx = An array giving the directional cosines for the  x-axis. DCOSX  will work if
          getcut is actived.
        dcosy = An array giving the directional cosines for the  y-axis. DCOSY  will work if
          getcut is actived.
        eomwl = A scalar giving the radar wavelength. The default value is 6m (50 MHZ).
        airwl = Set this input to float (or intger) to specify the wavelength (in meters) of
          the transmitted EOM wave in the air. The default value is 4m.

        Modification History
        --------------------
        Converted to Object-oriented Programming by Freddy Galindo, ROJ, 20 September 2009.
        """



        # Getting antenna configuration.
        if filename:
            setup = JroAntSetup.ReturnSetup(path=path,filename=filename,pattern=pattern)
    
            ues = setup["ues"]
            phase = setup["phase"]
            gaintx = setup["gaintx"]
            gainrx = setup["gainrx"]
            justrx = setup["justrx"]
            self.title = setup["title"]
        else:
            ues = kwargs["ues"]
            phase = kwargs["phases"]
            gaintx = kwargs["gain_tx"]
            gainrx = kwargs["gain_rx"]
            justrx = kwargs["just_rx"]
            self.title = kwargs.get("title", "JRO Pattern")

        # Defining attributes for JroPattern class.
        # Antenna configuration
        
        self.uestx = ues
        self.phasetx = phase
        self.gaintx = gaintx
        self.uesrx = ues
        self.phaserx = phase
        self.gainrx = gainrx
        self.justrx = justrx

        # Pattern resolution & method to model
        self.maxphi = maxphi
        self.nptsx = nptsx
        self.nptsy = nptsy
        self.fftopt = fftopt

        # To get a cut of the pattern.
        self.getcut = getcut

        maxdcos = numpy.sin(maxphi*Misc_Routines.CoFactors.d2r)
        if len(dcosx) == 0:
            dcosx = ((numpy.arange(nptsx,dtype=float)/(nptsx-1))-0.5)*2*maxdcos
        if len(dcosy) == 0:
            dcosy = ((numpy.arange(nptsy,dtype=float)/(nptsy-1))-0.5)*2*maxdcos
        self.dcosx = dcosx
        self.dcosy = dcosy
        self.nx = dcosx.size
        self.ny = dcosy.size*(getcut==0) + (getcut==1)

        self.eomwl = eomwl
        self.airwl = airwl

        self.kk = 2.*numpy.pi/eomwl

        self.pattern = None
        self.meanpos = None
        self.norpattern = None
        self.maxpattern = None

        

        self.getPattern()

    def getPattern(self):
        """
        getpattern method returns the modeled total antenna pattern and its mean position.

        Return
        ------
        pattern = An array giving the Modelled antenna pattern.
        mean_pos = A 2-elements array giving the mean position of the main beam.

        Examples
        --------
        >> [pattern, mean_pos] = JroPattern(pattern=2).getPattern()
        >> print meanpos
        [  8.08728085e-14  -4.78193873e-14]

        Modification history
        --------------------
        Developed by Jorge L. Chau.
        Converted to Python by Freddy R. Galindo, ROJ, 20 September 2009.
        """

        if (self.fftopt>0) and (self.getcut>0):
            #print  "Conflict bewteen fftopt and getcut"
            #print  "To get a cut of the antenna pattern uses ffopt=0"
            return None, None

        if (self.fftopt==0):
            # Getting antenna pattern using the array method
            self.pattern = self.__usingArray(rx=1)
            if (self.justrx==0):
                self.pattern = self.pattern*self.__usingArray(rx=0)

        elif (self.fftopt>0):
            # Getting antenna pattern using FFT method
            self.pattern = self.__usingFFT(rx=1)
            if (self.justrx==0):
                self.pattern = self.pattern*self.__usingFFT(rx=0)

        self.maxpattern = numpy.nanmax(self.pattern)
        self.norpattern = self.pattern/self.maxpattern
        if self.getcut==0:
            self.__getBeamPars()

    def __usingArray(self,rx):
        """
        __usingArray method returns the Jicamarca antenna pattern computed using array model

        pattern = dipolepattern x modulepattern

        Parameters
        ----------
        rx = Set to 1 to use the Rx information. Otherwise set to 0 for Tx.

        Return
        ------
        pattern = An array giving the modelled antenna pattern using the array model.

        Modification history
        --------------------
        Developed by Jorge L. Chau.
        Converted to Python by Freddy R. Galindo, ROJ, 20 September 2009.
        """

        if rx==1:
            ues = self.uesrx
            phase = self.phaserx
            gain = self.gainrx
        elif rx==0:
            ues = self.uestx
            phase = self.phasetx
            gain = self.gaintx

        ues = ues*360./self.airwl
        phase = phase*360./self.airwl

        for ii in range(4):
            if ii==0:dim = numpy.array([4,0,8,4])   # WEST
            elif ii==1:dim = numpy.array([0,0,4,4]) # NORTH
            elif ii==2:dim = numpy.array([0,4,4,8]) # EAST
            elif ii==3:dim = numpy.array([4,4,8,8]) # SOUTH
            xi = dim[0]; xf = dim[2]; yi = dim[1]; yf = dim[3]
            phase[xi:xf,yi:yf] = phase[xi:xf,yi:yf] + ues[ii]

        phase = -phase

        ar = self.eomwl*numpy.array([[0.5,6., 24.5],[0.5,6.,24.5]])
        nr = numpy.array([[12.,4.,2.],[12.,4.,2.]])
        lr = 0.25*self.eomwl*numpy.array([[0,0.,0],[0.,0,0]])

        # Computing module and dipole patterns.
        pattern = (numpy.abs(self.__dipPattern(ar,nr,lr)*self.__modPattern(phase,gain)))**2

        return pattern

    def __usingFFT(self,rx):
        """
        __usingFFT method returns the Jicamarca antenna pattern computed using The Fast Fou-
        rier Transform.

        pattern = iFFT(FFT(gain*EXP(j*phase)))

        Parameters
        ----------
        rx = Set to 1 to use the Rx information. Otherwise set to 0 for Tx.

        Return
        ------
        pattern = An array giving the modelled antenna pattern using the array model.

        Modification history
        --------------------
        Developed by Jorge L. Chau.
        Converted to Python by Freddy R. Galindo, ROJ, 20 September 2009.
        """

        if rx==1:
            ues = self.uesrx
            phase = self.phaserx
            gain = self.gainrx
        elif rx==0:
            ues = self.uestx
            phase = self.phasetx
            gain = self.gaintx

        ues = ues*360./self.airwl
        phase = phase*360./self.airwl

        for ii in range(4):
            if ii==0:dim = numpy.array([4,0,8,4])   # WEST
            elif ii==1:dim = numpy.array([0,0,4,4]) # NORTH
            elif ii==2:dim = numpy.array([0,4,4,8]) # EAST
            elif ii==3:dim = numpy.array([4,4,8,8]) # SOUTH
            xi = dim[0]; xf = dim[2]; yi = dim[1]; yf = dim[3]
            phase[xi:xf,yi:yf] = phase[xi:xf,yi:yf] + ues[ii]

        phase = -phase

        delta_x = self.eomwl/2.
        delta_y = self.eomwl/2.

        nxfft = 2048
        nyfft = 2048
        dcosx = (numpy.arange(nxfft) - (0.5*nxfft))/(nxfft*delta_x)*self.eomwl
        dcosy = (numpy.arange(nyfft) - (0.5*nyfft))/(nyfft*delta_y)*self.eomwl

        fft_gain = numpy.zeros((nxfft,nyfft))
        fft_phase = numpy.zeros((nxfft,nyfft))

        nx = 8
        ny = 8
        ndx =12
        ndy =12
        for iy in numpy.arange(ny):
            for ix in numpy.arange(nx):
                ix1 = nxfft/2-self.nx/2*ndx+ix*ndx
                if ix<(nx/2):ix1 = ix1 - 1
                if ix>=(nx/2):ix1 = ix1 + 1

                iy1 = nyfft/2-ny/2*ndx+iy*ndy
                if iy<(ny/2):iy1 = iy1 - 1
                if iy>=(ny/2):iy1 = iy1 + 1

                fft_gain[ix1:ix1+ndx-1,iy1:iy1+ndy-1] = gain[ix,ny-1-iy]
                fft_phase[ix1:ix1+ndx-1,iy1:iy1+ndy-1] = phase[ix,ny-1-iy]


        fft_phase = fft_phase*Misc_Routines.CoFactors.d2r

        pattern = numpy.abs(numpy.fft.fft2(fft_gain*numpy.exp(numpy.complex(0,1)*fft_phase)))**2
        pattern = numpy.fft.fftshift(pattern)

        xvals = numpy.where((dcosx>=(numpy.min(self.dcosx))) & (dcosx<=(numpy.max(self.dcosx))))
        yvals = numpy.where((dcosy>=(numpy.min(self.dcosy))) & (dcosy<=(numpy.max(self.dcosy))))

        pattern = pattern[xvals[0][0]:xvals[0][-1],yvals[0][0]:yvals[0][-1]]

        return pattern


    def __dipPattern(self,ar,nr,lr):
        """
        _dipPattern function computes the dipole's pattern to  the Jicamarca radar. The next
        equation defines the pattern as a function of the mainlobe direction:

        sincx = SIN(k/2*n0x*(a0x*SIN(phi)*COS(alpha)))/SIN(k/2*(a0x*SIN(phi)*COS(alpha)))
        sincy = SIN(k/2*n0y*(a0y*SIN(phi)*SIN(alpha)))/SIN(k/2*(a0y*SIN(phi)*SIN(alpha)))
        A0(phi,alpha) = sincx*sincy
        Parameters
        ----------
        ar = ?
        nr = ?
        lr = ?

        Return
        ------
        dipole = An array giving antenna pattern from the dipole point of view..

        Modification history
        --------------------
        Developed by Jorge L. Chau.
        Converted to Python by Freddy R. Galindo, ROJ, 20 September 2009.
        """

        dipole = numpy.zeros((self.nx,self.ny),dtype=complex)
        for iy in range(self.ny):
            for ix in range(self.nx):
                yindex = iy*(self.getcut==0) + ix*(self.getcut==1)

                argx = ar[0,0]*self.dcosx[ix] - lr[0,0]
                if argx == 0.0:    
                    junkx = nr[0,0]
                else:
                    junkx = numpy.sin(0.5*self.kk*nr[0,0]*argx)/numpy.sin(0.5*self.kk*argx)
                

                argy = ar[1,0]*self.dcosy[yindex] - lr[1,0]
                if argy == 0.0: 
                    junky = nr[1,0]
                else:
                    junky = numpy.sin(0.5*self.kk*nr[1,0]*argy)/numpy.sin(0.5*self.kk*argy)
                

                dipole[ix,iy] = junkx*junky

        return dipole

    def __modPattern(self,phase,gain):
        """
        ModPattern computes the module's pattern to the Jicamarca radar.  The next equation
        defines    the pattern as a function mainlobe direction:

        phasex = pos(x)*SIN(phi)*COS(alpha)
        phasey = pos(y)*SIN(phi)*SIN(alpha)

        A1(phi,alpha) = TOTAL(gain*EXP(COMPLEX(0,k*(phasex+phasey)+phase)))

        Parameters
        ----------
        phase = Bidimensional array (8x8) giving the phase (in meters) of each module.
        gain  = Bidimensional array (8x8) giving to  define modules  will be active  (ones)
          and which will not (zeros).

        Return
        ------
        module = An array giving antenna pattern from the module point of view..

        Modification history
        --------------------
        Developed by Jorge L. Chau.
        Converted to Python by Freddy R. Galindo, ROJ, 20 September 2009.
        """

        pos = self.eomwl*attenuation
        posx = pos[0,:,:]
        posy = pos[1,:,:]

        phase = phase*Misc_Routines.CoFactors.d2r
        module = numpy.zeros((self.nx,self.ny),dtype=complex)
        for iy in range(self.ny):
            for ix in range(self.nx):
                yindex = iy*(self.getcut==0) + ix*(self.getcut==1)
                phasex = posx*self.dcosx[ix]
                phasey = posy*self.dcosy[yindex]
                tmp = gain*numpy.exp(numpy.complex(0,1.)*(self.kk*(phasex+phasey)+phase))
                module[ix,iy] = tmp.sum()

        return module

    def __getBeamPars(self):
        """
        _getBeamPars computes the main-beam parameters of the antenna.

        Modification history
        --------------------
        Developed by Jorge L. Chau.
        Converted to Python by Freddy R. Galindo, ROJ, 20 September 2009.
        """

        dx = self.dcosx[1] - self.dcosx[0]
        dy = self.dcosy[1] - self.dcosy[0]

        amp = self.norpattern

        xx =  numpy.resize(self.dcosx,(self.nx,self.nx)).transpose()
        yy =  numpy.resize(self.dcosy,(self.ny,self.ny))

        mm0 = amp[numpy.where(amp > 0.5)]
        xx0 = xx[numpy.where(amp > 0.5)]
        yy0 = yy[numpy.where(amp > 0.5)]

        xc = numpy.sum(mm0*xx0)/numpy.sum(mm0)
        yc = numpy.sum(mm0*yy0)/numpy.sum(mm0)
        rc = numpy.sqrt(mm0.size*dx*dy/numpy.pi)

        nnx = numpy.where(numpy.abs(self.dcosx - xc) < rc)
        nny = numpy.where(numpy.abs(self.dcosy - yc) < rc)

        mm1 = amp[numpy.min(nnx):numpy.max(nnx)+1,numpy.min(nny):numpy.max(nny)+1]
        xx1 = self.dcosx[numpy.min(nnx):numpy.max(nnx)+1]
        yy1 = self.dcosy[numpy.min(nny):numpy.max(nny)+1]

        # fitting data into the main beam.
        
        params = gaussfit.fitgaussian(mm1)

        # Tranforming from indexes to axis' values
        xcenter = xx1[0] + (((xx1[xx1.size-1] - xx1[0])/(xx1.size -1))*(params[1]))
        ycenter = yy1[0] + (((yy1[yy1.size-1] - yy1[0])/(yy1.size -1))*(params[2]))
        xwidth  = ((xx1[xx1.size-1] - xx1[0])/(xx1.size-1))*(params[3])*(1/Misc_Routines.CoFactors.d2r)
        ywidth  = ((yy1[yy1.size-1] - yy1[0])/(yy1.size-1))*(params[4])*(1/Misc_Routines.CoFactors.d2r)
        meanwx = (xwidth*ywidth)
        meanpos = numpy.array([xcenter,ycenter])

        #print  'Position: %f %f' %(xcenter,ycenter)
        #print  'Widths:   %f %f' %(xwidth, ywidth)
        #print  'BWHP:     %f' %(2*numpy.sqrt(2*meanwx)*numpy.sqrt(-numpy.log(0.5)))

        self.meanpos = meanpos

class AmisrPattern():

    def __init__(self,azimuth=0, elevation=90,filename=None,nptsx=101,nptsy=101,maxphi=40,\
                dcosx=[],dcosy=[], fc=445, just_rx=False):
        """
        AMISR-14 class creates an object to represent the useful parameters for beam mode-
        lling of the Jicamarca VHF radar.

        Modification History
        --------------------
        Writen by Joab Apaza, ROJ, July 2023.
        """


        # Pattern resolution & method to model
        self.maxphi = maxphi
        self.nptsx = nptsx
        self.nptsy = nptsy
        
        self.just_rx = just_rx

        self.getcut = 0

        maxdcos = numpy.sin(maxphi*Misc_Routines.CoFactors.d2r)
        if len(dcosx) == 0:
            dcosx = ((numpy.arange(nptsx,dtype=float)/(nptsx-1))-0.5)*2*maxdcos
        if len(dcosy) == 0:
            dcosy = ((numpy.arange(nptsy,dtype=float)/(nptsy-1))-0.5)*2*maxdcos/3
        self.dcosx = dcosx
        self.dcosy = dcosy
        self.nx = dcosx.size
        self.ny = dcosy.size*(self.getcut==0) + (self.getcut==1)

        self.Cx, self.Cy = numpy.meshgrid(dcosx, dcosy)
        self.Cx0 = numpy.cos(numpy.radians(elevation))*numpy.sin(numpy.radians(azimuth))
        self.Cy0 = numpy.cos(numpy.radians(elevation))*numpy.cos(numpy.radians(azimuth))

        self.eomwl = 300/fc

        self.kk = 2.*numpy.pi/self.eomwl

        self.pattern = None
        self.meanpos = None
        self.norpattern = None
        self.maxpattern = None

        

        self.getPattern()

    def getPattern(self, xy_panelPos=[], wgts=[]):
        """
        getpattern method returns the modeled total antenna pattern and its mean position.

        Return
        ------
        pattern = An array giving the Modelled antenna pattern.
        xy_panelPos = A 2-elements array giving the mean position of panels.

        Modification history
        --------------------
        Writen by Joab Apaza, ROJ, July 2023.
        """
        panelPos = numpy.array([[-9.90000000e-01,  9.90000000e-01, -9.90000000e-01,
         9.90000000e-01, -9.90000000e-01,  9.90000000e-01,
        -9.90000000e-01,  9.90000000e-01, -9.90000000e-01,
         9.90000000e-01, -9.90000000e-01,  9.90000000e-01,
        -9.90000000e-01,  9.90000000e-01],
       [-1.04100000e+01, -1.04100000e+01, -6.94000000e+00,
        -6.94000000e+00, -3.47000000e+00, -3.47000000e+00,
        -1.77635684e-15, -1.77635684e-15,  3.47000000e+00,
         3.47000000e+00,  6.94000000e+00,  6.94000000e+00,
         1.04100000e+01,  1.04100000e+01]])
        
        if len(xy_panelPos)>0:
            xpos = xy_panelPos[0]
            ypos = xy_panelPos[1]

        else:
            xpos = panelPos[0]
            ypos = panelPos[1]
        
        if len(wgts)<1:
            wgts = numpy.ones(len(xpos))

        panel = self.__modPattern()

        farr = 0
        for i in range(len(xpos)):

            farr += wgts[i]*numpy.exp( 1j*self.kk*(xpos[i]*(self.Cx-self.Cx0) + ypos[i]*(self.Cy-self.Cy0)))

        #antenna pattern
        self.pattern = abs(farr*panel)**2
        

        if not self.just_rx:
            self.pattern *=self.pattern

        self.maxpattern = numpy.nanmax(self.pattern)
        self.norpattern = self.pattern/self.maxpattern
        if self.getcut==0:
            self.__getBeamPars()

    def __getBeamPars(self):
        """
        _getBeamPars computes the main-beam parameters of the antenna.

        Modification history
        --------------------
        Developed by Jorge L. Chau.
        Converted to Python by Freddy R. Galindo, ROJ, 20 September 2009.
        """

        dx = self.dcosx[1] - self.dcosx[0]
        dy = self.dcosy[1] - self.dcosy[0]

        amp = self.norpattern

        xx =  numpy.resize(self.dcosx,(self.nx,self.nx)).transpose()
        yy =  numpy.resize(self.dcosy,(self.ny,self.ny))

        mm0 = amp[numpy.where(amp > 0.5)]
        xx0 = xx[numpy.where(amp > 0.5)]
        yy0 = yy[numpy.where(amp > 0.5)]

        xc = numpy.sum(mm0*xx0)/numpy.sum(mm0)
        yc = numpy.sum(mm0*yy0)/numpy.sum(mm0)
        rc = numpy.sqrt(mm0.size*dx*dy/numpy.pi)

        nnx = numpy.where(numpy.abs(self.dcosx - xc) < rc)
        nny = numpy.where(numpy.abs(self.dcosy - yc) < rc)

        mm1 = amp[numpy.min(nnx):numpy.max(nnx)+1,numpy.min(nny):numpy.max(nny)+1]
        xx1 = self.dcosx[numpy.min(nnx):numpy.max(nnx)+1]
        yy1 = self.dcosy[numpy.min(nny):numpy.max(nny)+1]

        # fitting data into the main beam.
        
        params = gaussfit.fitgaussian(mm1)

        # Tranforming from indexes to axis' values
        xcenter = xx1[0] + (((xx1[xx1.size-1] - xx1[0])/(xx1.size -1))*(params[1]))
        ycenter = yy1[0] + (((yy1[yy1.size-1] - yy1[0])/(yy1.size -1))*(params[2]))
        xwidth  = ((xx1[xx1.size-1] - xx1[0])/(xx1.size-1))*(params[3])*(1/Misc_Routines.CoFactors.d2r)
        ywidth  = ((yy1[yy1.size-1] - yy1[0])/(yy1.size-1))*(params[4])*(1/Misc_Routines.CoFactors.d2r)
        meanwx = (xwidth*ywidth)
        meanpos = numpy.array([xcenter,ycenter])

        #print  'Position: %f %f' %(xcenter,ycenter)
        #print  'Widths:   %f %f' %(xwidth, ywidth)
        #print  'BWHP:     %f' %(2*numpy.sqrt(2*meanwx)*numpy.sqrt(-numpy.log(0.5)))

        self.meanpos = meanpos


    def __dipPattern(self):
        """
        _dipPattern function computes the dipole's pattern to  the AMISR radar. The next
        equation defines the pattern as a function of the mainlobe direction:


        Return dipole pattern
        ------
        dipole = An array giving antenna pattern from the dipole point of view..

        Modification history
        --------------------
        Developed by M. Milla.
        Writen by Joab Apaza, ROJ, July 2023.
        """

        Cz2 = 1-(self.Cx**2+self.Cy**2)
        Cz2[Cz2<=0] = numpy.nan
        Cz = numpy.sqrt(Cz2)
        Cy2 = 1 - self.Cy**2
        Cy2[Cy2<=0] = numpy.nan
        sinth = numpy.sqrt(Cy2)
        DIP = numpy.cos(numpy.pi/2*self.Cy)/sinth
        Gnd= 2j*numpy.sin(numpy.pi/2*Cz)

        return DIP*Gnd


    def __modPattern(self):
        """
        ModPattern computes the module's pattern to the amisr radar. We consider each panel as
        a module.


        Return
        ------
        module = An array giving antenna pattern from the module point of view.

        Modification history
        --------------------
        Writen by Joab Apaza, ROJ, July 2023.
        """
        _xypos = numpy.asarray([[ 0.        ,  0.        ],
            [-0.15802247,  0.43434   ],
            [ 0.        ,  0.86868   ],
            [-0.15802247,  1.30302   ],
            [ 0.        ,  1.73736   ],
            [-0.15802247,  2.1717    ],
            [ 0.        ,  2.60604   ],
            [-0.15802247,  3.04038   ],
            [ 0.495808  ,  0.        ],
            [ 0.33778553,  0.43434   ],
            [ 0.495808  ,  0.86868   ],
            [ 0.33778553,  1.30302   ],
            [ 0.495808  ,  1.73736   ],
            [ 0.33778553,  2.1717    ],
            [ 0.495808  ,  2.60604   ],
            [ 0.33778553,  3.04038   ],
            [ 0.991616  ,  0.        ],
            [ 0.83359353,  0.43434   ],
            [ 0.991616  ,  0.86868   ],
            [ 0.83359353,  1.30302   ],
            [ 0.991616  ,  1.73736   ],
            [ 0.83359353,  2.1717    ],
            [ 0.991616  ,  2.60604   ],
            [ 0.83359353,  3.04038   ],
            [ 1.487424  ,  0.        ],
            [ 1.32940153,  0.43434   ],
            [ 1.487424  ,  0.86868   ],
            [ 1.32940153,  1.30302   ],
            [ 1.487424  ,  1.73736   ],
            [ 1.32940153,  2.1717    ],
            [ 1.487424  ,  2.60604   ],
            [ 1.32940153,  3.04038   ]])
        
         
        FPanel = 0
        for i in range(32):
            xa = _xypos[i, 0]/self.eomwl
            ya = _xypos[i, 1]/self.eomwl
            FPanel += numpy.exp(1j*2*numpy.pi*xa*(self.Cx-self.Cx0) + 1j*2*numpy.pi*ya*(self.Cy-self.Cy0))
        
        dipole = self.__dipPattern()
        
        return FPanel*dipole

    


class overJroShow:

    __serverdocspath = ''
    __tmpDir = ''

    def __init__(self, site=1, title='', heights=None, maxphi=None,ploteo=0):
        self.year = None
        self.month = None
        self.dom = None
        self.pattern = None
        self.maxphi = maxphi
        self.heights = None
        self.filename = None
        self.showType = None
        self.path = None
        self.objects = None
        self.nptsx = 101
        self.nptsy = 101
        self.fftopt = 0
        self.dcosx = 1
        self.dcosy = 1
        self.dcosxrange = None
        self.dcosyrange = None
        self.maxha_min= 0.
        self.show_object = None
        self.dcosx_mag = None
        self.dcosy_mag = None
        self.ha_mag = None
        self.time_mag = None
        self.main_dec = None
        self.ObjC = None
        self.ptitle = title
        self.path4plotname = None
        self.plotname0 = None
        self.plotname1 = None
        self.plotname2 = None
        self.scriptHeaders = 0
        self.ploteo=ploteo
        if site==1:
            self.glat = -11.95
            self.glon = -76.8667
        else:
            self.glat = -11.953371
            self.glon = -76.874913

        self.UT = 5 #timezone

        self.site=site

        if heights is None:
            self.heights = numpy.array([100.,500.,1000.])
        else:
            self.heights = numpy.array(heights)
        
        self.pattern_plot = AntPatternPlot(self.ploteo)

    def initParameters(self, dt):

        self.year = dt.year
        self.month = dt.month
        self.dom = dt.day
        self.doy = datetime.datetime(dt.year,dt.month,dt.day).timetuple().tm_yday
        # Defining plot filenames
        self.path4plotname = os.path.join(self.__serverdocspath,self.__tmpDir)
        self.plotname0 = 'over_jro_0_%i.png'% (time.time()) #plot pattern & objects
        self.plotname1 = 'over_jro_1_%i.png'% (time.time()) #plot antenna cuts
        self.plotname2 = 'over_jro_2_%i.png'% (time.time()) #plot sky noise

        if self.site==1:
            # Defining antenna axes respect to geographic coordinates (See Ochs report).
    #        alfa = 1.46*Misc_Routines.CoFactors.d2r
    #        theta = 51.01*Misc_Routines.CoFactors.d2r

            alfa = 1.488312*Misc_Routines.CoFactors.d2r
            th = 6.166710 + 45.0
            theta = th*Misc_Routines.CoFactors.d2r
        else:
            #AMISR
            alfa = 0.1*Misc_Routines.CoFactors.d2r
            th = 0.0977
            theta = th*Misc_Routines.CoFactors.d2r

        self.maxha_min = 4*self.maxphi*numpy.sqrt(2)*1.25


        sina = numpy.sin(alfa)
        cosa = numpy.cos(alfa)
        MT1 = numpy.array([[1,0,0],[0,cosa,-sina],[0,sina,cosa]])
        sinb = numpy.sin(theta)
        cosb = numpy.cos(theta)
        MT2 = numpy.array([[cosb,sinb,0],[-sinb,cosb,0],[0,0,1]])
        self.MT3 = numpy.array(numpy.dot(MT2, MT1)).transpose()

        self.xg = numpy.dot(self.MT3.transpose(),numpy.array([1,0,0]))
        self.yg = numpy.dot(self.MT3.transpose(),numpy.array([0,1,0]))
        self.zg = numpy.dot(self.MT3.transpose(),numpy.array([0,0,1]))    

    def plotPattern(self, site, azimuth, elevation, date, phases, gain_tx, gain_rx, ues, just_rx, angle, plot=True):
        fullDCOSX = []
        fullDCOSY = []
        # Plotting Antenna patterns.
        self.maxphi = angle
        self.initParameters(date)
        self.junkjd = TimeTools.Time(self.year,self.month,self.dom).change2julday()
        self.junklst = TimeTools.Julian(self.junkjd).change2lst(longitude=self.glon)
        self.ra_obs = self.junklst*Misc_Routines.CoFactors.h2d
        
        date = TimeTools.Time(date.year,date.month,date.day).change2strdate(mode=2)
        if site==1:
            mesg = 'Over Jicamarca: ' + date[0]
            ObjAnt = JroPattern(pattern=0,
                            filename=None,
                            path=None,
                            nptsx=self.nptsx,
                            nptsy=self.nptsy,
                            maxphi=angle,
                            fftopt=self.fftopt,
                            phases=phases,
                            gain_tx=gain_tx,
                            gain_rx=gain_rx,
                            ues=ues,
                            just_rx=just_rx
                            )
        else:
            mesg = 'Over AMISR-14: ' + date[0]
            
            ObjAnt = AmisrPattern(azimuth,
                            elevation,
                            maxphi=angle,
                            nptsx=self.nptsx,
                            nptsy=self.nptsy,
                            just_rx=False
                            )
            cwd = os.getcwd()
            pointings = numpy.genfromtxt(cwd+'/utils/UMET_beamcodes.csv', delimiter=',')
            fullDCOSX = numpy.cos(numpy.radians(pointings[:,2]))*numpy.sin(numpy.radians(pointings[:,1]))
            fullDCOSY = numpy.cos(numpy.radians(pointings[:,2]))*numpy.cos(numpy.radians(pointings[:,1]))
        
        
        if plot:
            self.pattern_plot.contPattern(site=site,
                            iplot=0,
                            gpath=self.path4plotname,
                            filename=self.plotname0,
                            mesg=mesg,
                            amp=ObjAnt.norpattern,
                            x=ObjAnt.dcosx,
                            y=ObjAnt.dcosy,
                            getCut=ObjAnt.getcut,
                            title=self.ptitle,
                            save=False)

            
            self.pattern_plot.plotRaDec(site=site, 
                        gpath=self.path4plotname,
                        filename=self.plotname0,
                        jd=self.junkjd, 
                        ra_obs=self.ra_obs, 
                        xg=self.xg, 
                        yg=self.yg, 
                        x=ObjAnt.dcosx, 
                        y=ObjAnt.dcosy,
                        allAmisr_x=fullDCOSX,
                        allAmisr_y=fullDCOSY,
                        save=False)
            
        vect_ant = numpy.array([ObjAnt.meanpos[0],ObjAnt.meanpos[1],numpy.sqrt(1-numpy.sum(ObjAnt.meanpos**2.))])
        vect_geo = numpy.dot(scipy.linalg.inv(self.MT3),vect_ant)
        vect_polar = Misc_Routines.Vector(numpy.array(vect_geo),direction=1).Polar2Rect()
        [ra,dec,ha] = Astro_Coords.AltAz(vect_polar[1],vect_polar[0],self.junkjd).change2equatorial()
        self.main_dec = dec

    def plotBfield(self, date, plot=True):

        self.initParameters(date)
        ObjB = BField(self.year,self.doy,self.site,self.heights)
        [dcos, alpha, nlon, nlat] = ObjB.getBField()
        
        self.pattern_plot.plotBField('', '', dcos, alpha, nlon, nlat, self.dcosxrange, self.dcosyrange, ObjB.heights, ObjB.alpha_i, plot=plot)
        
        if not plot:
            self.junkjd = TimeTools.Time(self.year,self.month,self.dom).change2julday()
            Bhei = 0
            dcosx = self.pattern_plot.alpha_location[:,0,Bhei]
            dcosy = self.pattern_plot.alpha_location[:,1,Bhei]
            vect_ant = [dcosx,dcosy,numpy.sqrt(1.-(dcosx**2. + dcosy**2.))]
            vect_ant = numpy.array(vect_ant)
            vect_geo = numpy.dot(scipy.linalg.inv(self.MT3),vect_ant)
            vect_geo = numpy.array(vect_geo).transpose()
            vect_polar = Misc_Routines.Vector(vect_geo,direction=1).Polar2Rect()
            [ra,dec,ha] = Astro_Coords.AltAz(vect_polar[1,:],vect_polar[0,:],self.junkjd).change2equatorial()
            val = numpy.where(ha>=180)
            
            if val[0].size>0:
                ha[val] = ha[val] -360.
            
            val = numpy.where(numpy.abs(ha)<=self.maxphi)
            
            if val[0].size>2:
                self.dcosx_mag = dcosx[val]
                self.dcosy_mag = dcosy[val]
                self.ha_mag = ha[val]
                self.time_mag = 0
    


    def plotCelestial(self, objects, plot=True):

        ntod = 24.*16.

        tod = numpy.arange(ntod)/ntod*24.

        [month,dom] = TimeTools.Doy2Date(self.year, self.doy).change2date()

        jd = TimeTools.Time(self.year, month, dom, tod+self.UT).change2julday()

        self.pattern_plot.plotCelestial(
            jd,
            self.main_dec,
            tod,
            self.maxha_min,
            objects,
            self.glat,
            self.glon,
            self.xg,
            self.yg,
            self.dcosxrange,
            self.dcosyrange,
            plot=plot
            )

        
    def plotAntennaCuts(self, objects, phases, gain_tx, gain_rx, ues, just_rx):

        incha = 0.05 # min
        nha = numpy.int32(2*self.maxha_min/incha) + 1.
        newha = numpy.arange(nha)/nha*2.*self.maxha_min - self.maxha_min
        nha_star = numpy.int32(200./incha)
        star_ha = (numpy.arange(nha_star) - (nha_star/2))*nha_star

        #Init ObjCut for PatternCutPlot()
        subplots = len(objects)
        self.ObjCut = PatternCutPlot(subplots)

        for io in objects:
            if io==0:
                if self.dcosx_mag.size!=0:
                    dcosx = self.dcosx_mag
                    dcosy = self.dcosy_mag
                    dcosz = 1 - numpy.sqrt(dcosx**2. + dcosy**2.)
                    # Finding rotation of B respec to antenna coords.
                    [mm,bb] = scipy.polyfit(dcosx,dcosy,1)
                    alfa = 0.0
                    theta = -1.*numpy.arctan(mm)
                    sina = numpy.sin(alfa); cosa = numpy.cos(alfa)
                    MT1 = [[1,0,0],[0,cosa,-sina],[0,sina,cosa]]
                    MT1 = numpy.array(MT1)
                    sinb = numpy.sin(theta); cosb = numpy.cos(theta)
                    MT2 = [[cosb,sinb,0],[-sinb,cosb,0],[0,0,1]]
                    MT2 = numpy.array(MT2)
                    MT3_mag = numpy.dot(MT2, MT1)
                    MT3_mag = numpy.array(MT3_mag).transpose()
                    # Getting dcos respec to B coords
                    vector = numpy.array([dcosx,dcosy,dcosz])
                    nvector = numpy.dot(MT3_mag,vector)
                    nvector = numpy.array(nvector).transpose()
                    yoffset = numpy.sum(nvector[:,1])/nvector[:,1].size
                    ha = self.ha_mag*4.
                    time = self.time_mag
                    width_star = 0.1 # half width in minutes
                    otitle = 'B Perp. cut'
            elif io==1:
                if self.pattern_plot.dcosx_sun.size!=0:
                    dcosx = self.pattern_plot.dcosx_sun
                    dcosy = self.pattern_plot.dcosy_sun
                    ha = self.pattern_plot.ha_sun*4.0
                    time = self.pattern_plot.time_sun
                    width_star = 2. # half width in minutes
                    otitle = 'Sun cut'
            elif io==2:
                if self.pattern_plot.dcosx_moon.size!=0:
                    dcosx = self.pattern_plot.dcosx_moon
                    dcosy = self.pattern_plot.dcosy_moon
                    ha = self.pattern_plot.ha_moon*4
                    time = self.pattern_plot.time_moon
                    m_distance = 404114.6  # distance to the Earth in km
                    m_diameter = 1734.4    # diameter in km.
                    width_star = numpy.arctan(m_distance/m_diameter)
                    width_star = width_star/2./Misc_Routines.CoFactors.d2r*4.
                    otitle = 'Moon cut'
            elif io==3:
                if self.pattern_plot.dcosx_hydra.size!=0:
                    dcosx = self.pattern_plot.dcosx_hydra
                    dcosy = self.pattern_plot.dcosy_hydra
                    ha = self.pattern_plot.ha_hydra*4.
                    time = self.pattern_plot.time_hydra
                    width_star = 0.25 # half width in minutes
                    otitle = 'Hydra cut'
            elif io==4:
                if self.pattern_plot.dcosx_galaxy.size!=0:
                    dcosx = self.pattern_plot.dcosx_galaxy
                    dcosy = self.pattern_plot.dcosy_galaxy
                    ha = self.pattern_plot.ha_galaxy*4.
                    time = self.pattern_plot.time_galaxy
                    width_star = 25. # half width in minutes
                    otitle = 'Galaxy cut'
            hour = numpy.int32(time)
            mins = numpy.int32((time - hour)*60.)
            secs = numpy.int32(((time - hour)*60. - mins)*60.)

            ObjT = TimeTools.Time(self.year,self.month,self.dom,hour,mins,secs)
            subtitle = ObjT.change2strdate()

            star_cut = numpy.exp(-(star_ha/width_star)**2./2.)

            pol = scipy.polyfit(ha,dcosx,3.)
            polx = numpy.poly1d(pol)
            newdcosx = polx(newha)
            pol = scipy.polyfit(ha,dcosy,3.)
            poly = numpy.poly1d(pol)
            newdcosy = poly(newha)

            patterns = []
            for icut in numpy.arange(1):
                # Getting Antenna cut.
                Obj = JroPattern(dcosx=newdcosx,
                                    dcosy=newdcosy,
                                    getcut=1,
                                    path=self.path,
                                    phases=phases,
                                    gain_tx=gain_tx,
                                    gain_rx=gain_rx,
                                    ues=ues,
                                    just_rx=just_rx
                                    )

                patterns.append(Obj.pattern)


            self.ObjCut.drawCut(io,
                            patterns,
                            1,
                            newha,
                            otitle,
                            subtitle,
                            self.ptitle)

    

def skyNoise(jd, ut=-5.0, longitude=-76.87, filename='/app/utils/galaxy.txt'):
    """
    hydrapos returns RA and Dec provided by Bill Coles (Oct 2003).
    
    Parameters
    ----------
    jd = The julian date of the day (and time), scalar or vector.
    
    dec_cut = A scalar giving  the declination  to get a cut of the skynoise over Jica-
        marca. The default value is -11.95.
    filename = A string to specify name the skynoise file. The default value is skynoi-
        se_jro.dat
    
    Return
    ------
    maxra = The maximum right ascension to the declination used to get a cut.
    ra = The right ascension.
    
    Examples
        --------
    >> [maxra,ra] = skynoise_jro()
    >> print maxra, ra
    139.45 -12.0833333333
    
    Modification history
    --------------------             
    Converted to Python by Freddy R. Galindo, ROJ, 06 October 2009.
    """

    # Defining date to compute SkyNoise.
    [year, month, dom, hour, mis, secs] = TimeTools.Julian(jd).change2time()
    is_dom = (month==9) & (dom==21)
    if is_dom:
        tmp = jd
        jd = TimeTools.Time(year,9,22).change2julian()
        dom = 22

    # Reading SkyNoise        
    
    g = os.getcwd()

    f = open(filename)
    
    lines = f.read()
    f.close()
    
    nlines = 99
    lines = lines.split('\n')
    data = numpy.zeros((2,nlines))*numpy.float32(0.)
    for ii in numpy.arange(nlines):
        line = numpy.array([lines[ii][0:6],lines[ii][6:]])
        data[:,ii] = numpy.float32(line)

    # Getting SkyNoise to the date desired.
    otime = data[0,:]*60.0
    opowr = data[1,:]
    
    hour = numpy.array([0,23])
    mins = numpy.array([0,59])
    secs = numpy.array([0,59])
    LTrange = TimeTools.Time(year,month,dom,hour,mins,secs).change2julday()
    LTtime  = LTrange[0] + numpy.arange(1440)*((LTrange[1] - LTrange[0])/(1440.-1))
    lst = TimeTools.Julian(LTtime + (-3600.*ut/86400.)).change2lst()

    ipowr = lst*0.0
    # Interpolating using scipy (inside max and min "x")
    otime = otime/3600.
    val = numpy.where((lst>numpy.min(otime)) & (lst<numpy.max(otime))); val = val[0]
    tck = scipy.interpolate.interp1d(otime,opowr)
    ipowr[val] = tck(lst[val])

    # Extrapolating above maximum time data (23.75).
    uval = numpy.where(lst>numpy.max(otime))
    if uval[0].size>0:
        ii = numpy.min(uval[0])
        m = (ipowr[ii-1] - ipowr[ii-2])/(lst[ii-1] - lst[ii-2])
        b = ipowr[ii-1] - m*lst[ii-1]            
        ipowr[uval] = m*lst[uval] + b

    if is_dom:
        lst = numpy.roll(lst,4)        
        ipowr = numpy.roll(ipowr,4)

    new_lst = numpy.int32(lst*3600.)
    new_pow = ipowr

    return ipowr, LTtime, lst

def skynoise_plot(year, month, day):
    """
    getPlot method creates a skynoise map over Jicamarca for a desired day. Additionally
    save a PNG file of this plot.
    
    Examples
    --------
    >> SkyNoisePlot(skypower,skytime,skytime_lst).getPlot()
    
    Modification history
    --------------------             
    Created by Freddy R. Galindo, ROJ, 18 October 2009.
    """

    # Working with the time before to plot the SkyNoise

    julian = TimeTools.Time(year, month, day).change2julday()
    power, tm, time_lst = skyNoise(julian)

    secs = TimeTools.Julian(tm).change2secs()    
    secs_lst = time_lst*1.
    if time_lst.argmin()>0:
        secs_lst[time_lst.argmin():] = secs_lst[time_lst.argmin():] + 24.
        secs_lst = secs_lst*3600. 

    label_secs = time.localtime(secs[power.argmax()] + time.timezone)
    label_secs_lst = time.localtime(secs_lst[power.argmax()] + time.timezone)
    
    # Creating labels for main x-labelticks (Local Time):
    snow = TimeTools.Time(year, month, day, 0, 0, 0).change2secs()
    today = secs - snow
    xticks_dpos = []
    xticks_dval = []
    for ix in [0,120,240,360,480,600,720,840,960,1080,1200,1320,1439]:             
        xticks_dpos.append(today[ix])
        time_tuple = time.localtime(today[ix] + time.timezone)
        xticks_dval.append(time.strftime('%H:%M',time_tuple))
        if ix==1439:xticks_dval[12] = ''

    # Creating labels for secondary x-labelticks (Sidereal Time):
    xticks_upos = [secs_lst[0],secs_lst[-1]]
    xticks_uval = ['','']
    for ix in [0,120,240,360,480,600,720,840,960,1080,1200,1320]:               
        index = numpy.argmin(numpy.abs(time_lst - (ix/60. + 1.)))
        xticks_upos.append(secs_lst[index])
        time_tuple = time.localtime(secs_lst[index] + time.timezone)
        if time_tuple[4]==59:
            tmp_list = list(time_tuple)
            tmp_list[4] = 0
            tmp_list[3] = tmp_list[3] + 1
            time_tuple = tuple(tmp_list)            
        xticks_uval.append(time.strftime('%H:%M',time_tuple))

    # Creating SkyNoise Map.
    fig = Figure()
    fig.subplots_adjust(top=0.80)
    
    main_mesg = "SKY BRIGHTNESS AT 50Mhz - Date: " 
    main_mesg = main_mesg + time.strftime("%d-%b-%Y (%j)",label_secs)
    main_mesg = main_mesg + "\nGalaxy Pass at " + time.strftime("%H:%M:%S LT",label_secs)
    main_mesg = main_mesg + time.strftime("  (%H:%M:%S LST)",label_secs_lst)
    
    fig.suptitle(main_mesg, fontsize=12)
    axl = fig.add_subplot(1,1,1)
    
    axl.plot(today,power,'b-')       
    axr = axl.twinx()        
    axl.set_xlabel('Local Time (Hours)')
    axl.set_ylabel('Signal Strengh (mA)')
    axr.set_ylabel('Temperature [K]x10^3')
    axl.set_xlim(0,24)
    axl.set_ylim(0,115)
    axr.set_ylim(0,50)
    axl.set_xticks(xticks_dpos)
    axl.set_xticklabels(xticks_dval, fontsize=8)
    axl.grid()

    axt = axl.twiny()
    axt.set_xlabel('Local Sidereal Time (Hours)')
    axt.set_xlim(secs_lst[0],secs_lst[-1])
    axt.set_xticks(xticks_upos)
    axt.set_xticklabels(xticks_uval, fontsize=8)
    
    buf = BytesIO()
    fig.savefig(buf, format="png")
    
    return buf

def overjro_plot(site, pattern_id, date, angle, height, bodys,  azimuth, elevation):
    
    if site!=1:
        pattern_id = 0

    pattern = select_pattern(pattern = pattern_id)

    if site==1: #Jicamarca
        ob = overJroShow(site, pattern['title'], heights=height)

    else:   #AMISR-14
        if azimuth > 180 or azimuth < -180:
            azimuth = 180
        if elevation > 90:
            elevation = 90
        if elevation < 49.93:
            elevation = 49.93
        beam,azimuth,elevation = findRealBeam(azimuth,elevation)
        tamisr = 'AMISR-14 at {:2.2f}° Az - {:2.2f}° El'.format(azimuth, elevation)
        ob = overJroShow(site, title=tamisr, heights=height, maxphi=10)

    if site==2 and angle==5:
        angle=40

    ob.plotPattern(
        site,
        azimuth,
        elevation,
        date, 
        pattern['phase'], 
        pattern['gaintx'], 
        pattern['gainrx'], 
        pattern['ues'], 
        pattern['justrx'], 
        angle
        )
    
    if 'bfield' in bodys:
        ob.plotBfield(date)
    
    objects = []
    for body in bodys:
        if body=='sun':
            objects.append(1)
        elif body=='moon':
            objects.append(2)
        elif body=='hydra':
            objects.append(3)
        elif body=='galaxy':
            objects.append(4)
    
    if objects:
        ob.plotCelestial(objects)

    buf = BytesIO()
    ob.pattern_plot.fig.savefig(buf, format="png")
    
    return buf

def antennacuts_plot(pattern_id, date, angle, height, bodys):

    pattern = select_pattern(pattern = pattern_id)

    ob = overJroShow(pattern['title'], heights=height, maxphi=angle)
    
    if 'bfield' in bodys:
        ob.plotBfield(date, plot=False)

    
    objects = []
    for body in bodys:
        if body=='bfield':
            objects.append(0)
        if body=='sun':
            objects.append(1)
        elif body=='moon':
            objects.append(2)
        elif body=='hydra':
            objects.append(3)
        elif body=='galaxy':
            ob.plotPattern(date, pattern['phase'], pattern['gaintx'], pattern['gainrx'], pattern['ues'], pattern['justrx'], angle, plot=False)
            objects.append(4)
    
    if objects:
        ob.plotCelestial(objects, plot=False)
        ob.plotAntennaCuts(
            objects,
            pattern['phase'], 
            pattern['gaintx'], 
            pattern['gainrx'], 
            pattern['ues'], 
            pattern['justrx'],
        )

    buf = BytesIO()
    ob.ObjCut.fig.savefig(buf, format="png")
    
    return buf


def muf_plot(lat,lon,version,date=None):
    now=datetime.datetime.now()
    obj=plot_IRI(latitude=lat,
              longitude=lon,
              version_IRI=version)
    
    if (date!=None ):
        if ((now.day==date.day) and (now.month==date.month) and (now.year==date.year)):
            
            pass
        
        else:
            obj.set_date(year=date.year,
                     month=date.month,
                     day=date.day)
        
    obj.run()
    buf=obj.get_buf()

    del obj 
    
    return buf

def findRealBeam(azimuth, elevation,path=''):
    cwd = os.getcwd()
    
    if len(path) != 0:
        validpointings =  numpy.genfromtxt(path, delimiter=',')
    else:
        validpointings =  numpy.genfromtxt(cwd+'/utils/UMET_beamcodes.csv', delimiter=',')
    nbhr = 10 #angle margin
    nbhr = 6
    beam1min = numpy.where(validpointings[:,1] > (azimuth-nbhr))
    beam1max = numpy.where(validpointings[:,1] < (azimuth+nbhr))
    beam1 = numpy.intersect1d(beam1min, beam1max)

    errs=[]
    i = 0
    DCOSX = numpy.cos(numpy.radians(elevation))*numpy.sin(numpy.radians(azimuth))
    DCOSY = numpy.cos(numpy.radians(elevation))*numpy.cos(numpy.radians(azimuth))
    for beam in beam1:
        el = validpointings[beam,2]
        az = validpointings[beam,1]
        dcosx = numpy.cos(numpy.radians(el))*numpy.sin(numpy.radians(az))
        dcosy = numpy.cos(numpy.radians(el))*numpy.cos(numpy.radians(az))
        err = math.sqrt( 0.5*((DCOSX-dcosx)**2) + 0.5*((DCOSY-dcosy)**2))
        errs.append(err)
        i+=1
    
    errsArr = numpy.asarray(errs)
    sorted = numpy.argsort(errsArr)
    nbeam = beam1[sorted[0]]

    return validpointings[nbeam,0:3]
