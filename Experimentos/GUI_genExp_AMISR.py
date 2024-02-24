import datetime
import os
import numpy as np
import scipy.io
from matplotlib import pyplot as plt

def insert_code(tuf=None,tx_start=0,sCode=None, nbaud=0, ncode=0, tx_baud=0):
 # Missing ncode, now is just for 1 code input, MISSING for multiple codes
    print("inserting code")
    for i in range(len(sCode)):
        if sCode[i]=='0':
            tuf.write("5\t\t%d" %(tx_start+(tx_baud*(i))) + "\t\t%d"%(tx_baud) + "\t\t0\t\t0\t\tphase %d \n"%i)


def rangecorrect(rangepath=None, sampletime=0, off=0):
    st = sampletime
    if st<10:
        file = rangepath+'blackman_0%2.2f'%st+'usec.rxc.mat'
    else:
        file = rangepath+'blackman_%2.2f'%st+'usec.rxc.mat'
    print(file)
    f = scipy.io.loadmat(file) 

    if st >= 1:
        mul = np.ceil(150/st)
        modu = off%st
        o = mul * st + modu
        i = f['offset']
        ind = np.where(i==o)
        c = f['correction']
        corr = c[ind]
    else:
        print("no corr")
        corr = 0
    print("corr = {}".format(corr))
    return corr

# variables a modificar
# rutas y nombres de archivos
''' 
mainpath = os.getcwd()
rangepath = os.path.join(mainpath,"offsets/" )
print (mainpath,rangepath)
#sExpName = 'winds'

sExpName = 'ISR_1Beam_Faraday_11_23_v5'

# rxc_filename = 'blackman_01.00usec_020505_20.0.rxc'     #Meteors
rxc_filename = 'blackman_10.00usec_051010_20.0.rxc'     #ESF, ISR_5

#rxc_filename = 'blackman_05.00usec_020525_20.0.rxc'    #EEJ, ISR

a = datetime.date.today()
date = a.strftime("%Y%m%d")
print(date)
path = mainpath+'/' + sExpName +'_' + date
if not (os.path.exists(path)):
    print("creating path...",path)
    os.mkdir(path)
'''
# ##################################################################
# ##################################################################
class Parametros():

    def __init__(self,Title,pps,ipp,TX,sCode,nCode,nBaud,ndh,nsa,nProfBlock,nBlockFile,sbeam,file,user):

        self.Title      = Title
        self.pps        = pps
        self.ipp        = ipp
        self.TX         = TX
        self.sCode      = sCode
        self.ncode      = nCode
        self.nbaud      = nBaud
        self.ndh        = ndh
        self.nsa        = nsa
        self.nProfBlock = nProfBlock
        self.nBlockFile = nBlockFile
        self.sbeam      = sbeam
        self.file       = file
        self.user       = user

class Generate_Experiments():
    
    def __init__(self,param):

        #Nombre del Experimento
        sExpName          =  param.Title
        #Configuracion de experimento
        self.SYNC         = param.pps # us =  1 second
        #Pulsos ISR 2023-ISR 
        self.IPP          =  param.ipp        # 10000
        self.nTxA         =  param.TX         # 600
        self.sCode        =  param.sCode      #None
        self.ncode        =  param.ncode      #0
        self.nbaud        =  param.nbaud      #1
        self.ndh          =  param.ndh        #10        #en us
        self.nsa          =  param.nsa        #980
        self.nProfBlock   =  param.nProfBlock #600  #multiplo de los canales e ipps
        self.nBlockFile   =  param.nBlockFile #10 #
        self.file         =  param.file

        sbeam             =  param.sbeam #[ '0xF923'] 
        nbeam             =  len(sbeam)
        user              =  param.user
        '''
        Pointing : 	 Azimuth 	 Elevation 	    Hex
        Beam  1 :  63779 	180.0	70.27 	0xF923
        '''
        bsync = 0            #boolean 0:ext (PPS) 1: interno
        RxBand = 445
        TxBand = 445

        # Consideraciones para procesamiento
        bBeamProf = 0     #boolean 0:beam cambia cada IPP 1: beams seguidos hasta completar nFFT IPPs
        nFFT =    1          #si puntos de FFT, entonces perfiles seguidos con mismo apunte ( 16 for EEJ;  32 for ESF;16 for ESF_EWDrift)

        nTX = 0      #IPPs transmitidos en 1 intervalo de SYNC, si se nTx = 0 se llena todo el SYNC

        # ##################################################################
        # ##################################################################
        # tiempos predeterminados de AMISR
        samp_start = 200 #wind=400, else = 200,  meteors = 810
        aeu_start = 401
        uc_start = 408
        dc_start = 408
        tx_start = 410
        rx_reset = 15

        # ##################################################################
        # otros parametros calculados
        baud = self.nbaud
        tx_baud = self.nTxA/self.nbaud
        osamp = tx_baud*1.0/self.ndh # oversampling
        print("oversampling: ",osamp)
        if self.nTxA > 0:
            tx_enable = 1

        mainpath = os.getcwd()
        mainpath = mainpath.split('/')
        mainpath="/".join(mainpath[:-1])+"/Experimentos"
        rangepath = os.path.join(mainpath,"offsets/" )
        print (mainpath,rangepath)

        a = datetime.date.today()
        date = a.strftime("%Y%m%d")
        print(date)
        path = mainpath+'/' + sExpName +'_' + date
        
        if not (os.path.exists(path)):
            print("creating path...",path)
            os.mkdir(path)

        # NOTA: hay que corregir rangecorrection y firstrange dependiendo del experimento
        if self.sCode != None:
            phase_code_enable = 1
            aeu_length = tx_start - aeu_start + tx_baud*baud + 1
            uc_length = tx_start - uc_start + tx_baud*baud + 1
            dc_length = tx_start - dc_start + tx_baud*baud + 1
            firstrange = ((samp_start-(tx_start+tx_baud/2))*150)
        else:
            phase_code_enable = 0
            aeu_length = tx_start - aeu_start + tx_baud*baud + 1
            uc_length = tx_start - uc_start + tx_baud*baud + 1
            dc_length = tx_start - dc_start + tx_baud*baud + 1
            firstrange = ((samp_start-(tx_start+tx_baud/2))*150)

        #def rangecorrect(rangepath=None, sampletime=0, offset=0):
        rangecorrection = rangecorrect(rangepath, self.ndh, (samp_start - rx_reset))
        print("range correction")
        if bsync == 0: #boolean 0:ext (PPS) 1: interno
            if (self.SYNC%self.IPP) != 0:
                print("WARNING: If PPS is used, then IPP = %d"%self.IPP +"us should be submultiple of 1 sec")
                var = input('Do you want to proceed? (y/n): ')
                if var == 'n':
                    exit()

            if (not float(nTX).is_integer()) and (nTX != 0):
                print("WARNING: the nTx = {%d} must be and integer".format(nTX))
                var = input('Do you want to proceed? (y/n): ')
                if var == 'n':
                    exit()

            max_nIPPs = self.SYNC/self.IPP # max NTX for SYNC range
            if nTX == 0:
                nIPPs = max_nIPPs   #se llena todo el SYNC
            else:
                nIPPs = nTX*nbeam*nFFT

            if max_nIPPs < nIPPs :
                print("WARNING:  The number of nTx = {}".format(nTX) +" exced = {}".format(max_nIPPs))
            print("nIPPs = {}".format(nIPPs))

        else:
            #nIPPs =  round(SYNC/(nbeam*IPP))
            nIPPs = round(self.SYNC/self.IPP)

        if (self.nProfBlock%nIPPs)!=0:
            print("WARNING: Number of pulses per block  = {}".format(self.nProfBlock)+" must be multiple of nIPPs = {}".format(nIPPs))
            var = input("Do you want to proceed? (y/n): ")
            if var == 'n':
                exit()

        if bBeamProf == 1: #repetir beam dependiendo de nFFT
            if (nIPPs*1.0%(nFFT*nbeam))!=0:
                print("WARNING: If nFFT is used, then ProfBlock = {}".format(nIPPs) +" should be multiple of number of FFT {}".format(nFFT) + " * #beams {}".format(nbeam) + " = {}".format(nFFT*nbeam) )
                repeat = nIPPs/nFFT*nbeam
                var = input('Do you want to proceed with repeat = {}? (y/n): '.format(repeat))
                if var == 'n':
                    exit()
            else:
                repeat =  nIPPs*1.0/(nFFT*nbeam)
        else:
            repeat =  nIPPs*1.0/(nbeam)

        ## ##################################################################
        pulse_pps = np.zeros(int(self.SYNC))
        pulse_pps[0:500000] = 1
        #print(pulse_pps)
        pulse_line = np.zeros(int(self.SYNC))

        for i in range(int(nIPPs)):
            tx_init = int(tx_start)+(i*int(self.IPP))
            tx_end = tx_init + int(self.nTxA)
            pulse_line[tx_init:tx_end] = 1

        fig,ax = plt.subplots(2,1, figsize=(15, 8))

        fig.suptitle('EXPERIMENT PULSES OVER {} us'.format(self.SYNC), fontsize=14)
        ax[0].plot(pulse_pps)
        ax[0].set_xlabel("Microseconds")
        ax[0].set_ylabel("Pulses pps")
        ax[0].grid()
        ax[1].plot(pulse_line)
        ax[1].set_xlabel("Microseconds")
        ax[1].set_ylabel("Pulses RF")
        ax[1].grid()
        #plt.show()

        # ##################################################################
        #Copy file
        if self.file == "01":
            fco = "blackman_%s.00usec_020505.fco"%(self.file)
            rxc = "blackman_%s.00usec_020505_20.0.rxc"%(self.file)
        elif self.file == "05":
            fco = "blackman_%s.00usec_020525.fco"%(self.file)
            rxc = "blackman_%s.00usec_020525_20.0.rxc"%(self.file)
        elif self.file == "10":
            fco = "blackman_%s.00usec_051010.fco"%(self.file)
            rxc = "blackman_%s.00usec_051010_20.0.rxc"(self.file)

        command_fco = "cp %s/Filtros/%s %s/"%(mainpath,fco,path)
        command_rxc = "cp %s/Filtros/%s %s/"%(mainpath,rxc,path)
        print("FILES *.FCO ",command_fco)
        print("FILES *.RXC ",command_rxc)
        os.system(command_fco)
        os.system(command_rxc)
        print("")

        # ##################################################################
        # Creacion de EXP file
        EXPfile= path + '/' + sExpName + '.exp'
        print("Creating %s" %sExpName + "File" )
        exp = open(EXPfile,'w')
        exp.write("[File] \n")
        exp.write(" Name=%s" %sExpName + ".v01\n")
        exp.write(" Description=Experiment, %d" %self.IPP + " us IPP with %d" %(baud*tx_baud)  +" us pulse and %d repeats per SYNC\n"%repeat)
        exp.write("\n")
        exp.write("[Log Info]\n")
        exp.write("%s\n"%(user))
        exp.write("Purpose=Pulse to pulse coded pulse\n")
        exp.write("Comment=Version v01\n")
        exp.write("\n")
        exp.write("[Common Parameters]\n")
        # In an integration mode is the number of pulses to integrate, in other modes is the number of pulses to write to 1 data record.
        '''
        1. Integrating mode
        You decide how much online time integration you want. Then you compute your frametime (total time before we start over)
        in the tufile basically adding all the IPP times up. Now you divide the frametime by the total number of pulses transmitted
        in a frame. This yields average time / IPP (pulse).
        You can now divide the integration time with that number and find the total number of pulses to integrate to achieve that
        online integration time. Then you set npulsesint to that number. If you are doing complex timing unit files like alternating
        codes etc you adjust the npulsesint to get a multiple of the number of pulses needed to complete the code and so on. This is
        why we don't just specify integration time to begin with. The system would have no way of knowing about complex tufiles so
        we might truncate the last pulse sequence in the code.

        2. Non-integrating mode
        Here you would have the same considerations as above concerning complex codes but other than that you are simply setting the
        npulsesint to the number of pulses to write to 1 data record in the file. I you use a big number you are using more memory
        in the computer to hold all the pulses before a disk write. If you set the number too small you end up not being able to keep
        up writing to disk.
        '''
        exp.write("npulsesint=%d\n" %self.nProfBlock)
        exp.write("recordsperfile=%d\n" %self.nBlockFile )
        exp.write("\n")
        exp.write("[copy files]\n")
        exp.write("\n")
        exp.write("[include data]\n")
        exp.write("\n")
        exp.write("[Hosts]\n")
        exp.write("DTC0=Minimum peak sidelobe code\n")
        exp.write("APS=Array Proxy\n")
        exp.write("\n")
        exp.write("[Modes]\n")
        exp.write("0=raw\n")
        exp.write("\n")
        exp.write("[Common Mode:0]\n")
        exp.write("beamcodefile=" + sExpName + ".v01.bco\n")
        exp.write("nbeamcodes=%d\n" %nbeam)
        exp.write("\n")
        exp.write("[DTC0 Mode:0]\n")
        exp.write("writedisplayrecord=0\n")
        exp.write("modes=raw\n") # This parameter determine which kind of data is going to be stored. If this is just "data" then the RadacHeader is going to be missing.
        exp.write("RxAttenuation=0\n")
        exp.write("txenabled=1\n")
        exp.write("tufile= " + sExpName + ".v01.dtc0.tuf\n")
        exp.write("rxconfig= %s\n" %rxc)
        exp.write("headerenable= 1\n")
        exp.write("maxsamples=%d\n" %self.nsa)
        exp.write("PulseWidth=%d\n" %(baud*tx_baud))
        exp.write("baud=%d\n" %baud)
        exp.write("txbaud=%d\n" %tx_baud)
        if self.sCode!= None:
            #exp.write("codelength=%d\n" %(2*baud)) # ?? length of total time of code or number of bauds in code?
            exp.write("codelength=%d\n" %(self.nTxA/self.ndh)) # ?? ?
        else:
            exp.write("codelength=%d\n" %(1)) # ?? length of total time of code or number of bauds in code?
        exp.write("rangecorrection=%f\n" %rangecorrection)
        exp.write("CollectStatistics=0\n")
        exp.write("RxBand=%d\n" %RxBand)
        exp.write("TxBand=%d\n" %TxBand)
        exp.write("txfrequency=tx0frequency1\n")
        exp.write("\n")
        exp.write("[DTC0 Mode:0,raw]\n")
        exp.write("displayexclude=1\n")
        exp.write("name=Data\n")
        exp.write("mode=raw11\n")
        exp.write("datatype=1\n")
        exp.write("modegroup=1\n")
        exp.write("indexsample=0\n")
        exp.write("ngates=%d\n" %self.nsa) # The number of samples you want to store. This will most often = "the number of samples in the data sample window" but could be set to less if one wanted to.
        exp.write("firstrange= %d\n" %firstrange)
        exp.write("\n")
        if bsync!=None:
            exp.write("[dtc0] \n")
            # External trigger (1 PPS signal)
            exp.write("internaltrig=%d"%bsync + "\n")
            exp.write("\n")
        exp.close()

        # ##################################################################
        # Creacion de TUF file
        TUFfile= path + '/' + sExpName + ".v01.dtc0.tuf"
        print("Creating " +sExpName + " TUFFile")
        tuf = open(TUFfile,'w')
        for i in range(int(nIPPs)):
            tuf.write("*\n")
            tuf.write("-1\t\t0\t\t%d"%self.IPP+"\t\t0\t\t0\t\tIPP %d"%i+"\n")
            tuf.write("-2\t\t1\t\t0\t\t2\t\t1\t\tDatatype\n")
            tuf.write("-2\t\t2\t\t0\t\t3\t\t1\t\tModegroup\n")
            tuf.write("-2\t\t3\t\t0\t\t4\t\t0\t\tNsPulse high word\n")
            tuf.write("-2\t\t4\t\t0\t\t5\t\t%d" %self.nsa+ "\t\t  NsPulse low word\n")
            # FALTA PONER EL CODIGO EN HEX PARA PONERLO AQUI
            if ((phase_code_enable == 1) and (tx_enable == 1)):
                tuf.write("-2\t\t5\t\t0\t\t6\t\t0x00C0\t\tCode MSW\n") ## ??
                tuf.write("-2\t\t6\t\t0\t\t7\t\t0xFF03\t\tCode \n") ## ??
                tuf.write("-2\t\t7\t\t0\t\t8\t\t0x0303\t\tCode \n") ## ??
                tuf.write("-2\t\t8\t\t0\t\t9\t\t0x0CF3\t\tCode LSW\n") ## ??

            tuf.write("12\t\t1\t\t5\t\t0\t\t0\t\tBeam code trigger\n")
            #tuf.write("11\t\t15\t\t5\t\t0\t\t0\t\tLoad header/sync Rx\n")
            tuf.write("11\t\t%d"%rx_reset + "\t\t5\t\t0\t\t0\t\tLoad header/sync Rx\n")
            tuf.write("7\t\t%d" %samp_start + "\t\t%d" %(self.nsa*self.ndh) + "\t\t0\t\t0\t\tData window, %d" %self.nsa + " samples\n")
            tuf.write("0\t\t%d"%aeu_start + "\t\t%d" %(aeu_length) + "\t\t0\t\t0\t\tAEU T/R\n")
            tuf.write("10\t\t%d" %(uc_start) + "\t\t%d" %(uc_length) + "\t\t0\t\t0\t\tUC T/R\n")
            tuf.write("13\t\t%d" %(dc_start) + "\t\t%d" %(dc_length) + "\t\t0\t\t0\t\tDC T/R\n")
            if (tx_enable == 1):
                tuf.write("-2\t\t%d"%(tx_start) + "\t\t0\t\t0\t\t1\t\tTx profile 1\n")
                if (phase_code_enable==1):
                    # si multiples codigos debe depender de repeticiones por beam
                    insert_code(tuf, tx_start,self.sCode, self.nbaud, self.ncode, tx_baud)
                tuf.write("-2\t\t%d"%(baud*tx_baud+tx_start) + "\t\t0\t\t0\t\t0\t\tTx profile 0, zero phase\n")
            tuf.write("-9\t\t0\t\t0\t\t0\t\t0\t\tend of ipp\n")

        tuf.close()


        # ##################################################################
        # Creacion de BCO file
        BCOfile= path + '/' +sExpName + ".v01.bco"
        print("Creating %s" %sExpName + " BCOFile")
        bco = open(BCOfile,'w')
        beamlist = []

        #reordena beams
        if (bBeamProf==1):
            print("con nFFT")
            for i in range(len(sbeam)):
                temp = [sbeam[i]]*nFFT
                beamlist.extend(temp)
        else:
            beamlist = sbeam*(int(nIPPs/nbeam))
            print("sin nFFT")

        # repite beams de ser necesario
        beamlist = beamlist*int(repeat)

        for i in range(int(nIPPs)):
            bco.write("%s"%beamlist[i] + "\n")

        bco.close()

'''

Title     = "ISR_Experiment"
pps       = 1000000
ipp       =   10000
nTxA      =     600
sCode     =    None
ncode     =       0 
nbaud     =       1
nh0       =       0
ndh       =      10
nsa       =     980
nProfBlock=     600
nBlockFile=      10
sbeam     = ['0xF923']
files     =     "05"   #05   #10
 
op =Parametros(Title,pps,ipp,nTxA,sCode,ncode,nbaud,nh0,ndh,nsa,nProfBlock,nBlockFile,sbeam,files)

Generate_Experiments(op)
'''