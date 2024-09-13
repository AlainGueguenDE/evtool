#!/usr/bin/python3
"""! Copyright (C) 2019-2023:
!   Max Planck Institut fuer extraterrestrische Physik (MPE)
!   Leibniz Institute for Astrophysics Potsdam (AIP)
!
! This file is part of eROSITA's Science Analysis Software System
! (eSASS) (https://erosita.mpe.mpg.de/). It is free software; you can
! redistribute it and/or modify it under the terms of the GNU General
! Public License as published by the Free Software Foundation; either
! version 2 of the License, or (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
! General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this library; if not, write to the Free Software
! Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307,
! USA, or see https://www.gnu.org/licenses/.
!-----------------------------------------------------------------------
!
"""
import sys
import os
from datetime import datetime
from astropy.io import fits
#import astropy.units as u
#import astropy.wcs as wcs
from astropy import wcs
from astropy.coordinates import SkyCoord
import numpy as np

"""
!-----------------------------------------------------------------------------
! AGueguen 20230315
! inspired by evtool*.f90 from eSASS
!
! Program read a source file, filter it and extract event falling in a given
! region, passed in parameters.
! save the result in a new file in the download folder.
! This python version is optimised for the DATool backend, it offers only a
! subset of the evtool official features
!Input parameters:
!
! file          FITS file name
! ra0           RA [deg] of projection centre
! dec0          DEC [deg] of projection centre
! Radius        Radius of the circle in degree , decimal values
! resultfolder  folder to save the extracted region
!Output :
! produce a new file
!
"""
# evtool eventfiles='%s' outfile='%s' region='%s' image='yes' size='%s' rebin='80'
# overlap='2.0'
# python evtool.py sourcefile=./eb01_123051_020_EventList_c020.fits.gz ra0=122 dec0=37
# Radius=0.25 result=./

neededkeys=['sourcefile','ra0','dec0','Radius','outfolder']

def extractonregion(sourcefilename,resultfilename,ra_cen,dec_cen,radius,image,dimimg,rebin,overlap):
    """ main extraction function , produce a new file """
    headermain = fits.getheader(sourcefilename)
    hdul1=fits.open(sourcefilename)
    new_hdul = fits.HDUList()
    new_hdul.append(fits.PrimaryHDU(header=headermain))

    # create the structure of the outpput fits file
    for indexextens in range(len(hdul1)):
        if indexextens>0:
            collist=[]
            for  colloc in hdul1[indexextens].columns:
                #arraytmp=hdul1[indexextens].data[colloc.name][0:20]
                col1 = fits.Column(name=colloc.name, format=colloc.format,array=[])
                collist.append(col1)
            #header1=hdul1[indexextens].header
            coldefs = fits.ColDefs(collist)
            tmpextens = fits.BinTableHDU.from_columns(coldefs)
            tmpextens.name = hdul1[indexextens].name
            new_hdul.append(tmpextens )
    # copy the original header to the new extracted fits
    for indexextens in range(len(hdul1)):
        if indexextens>0:
            header1=hdul1[indexextens].header
            newhead=new_hdul[indexextens].header
            lsth=list(header1.keys())
            for i in lsth :#range(len(header1)):
                newhead[i]=(header1[i],header1.comments[i])

    # go once again through all the extensions and copy interesting data#
    # NOTE  : On october 23 23 deactivate all filtering on extensions other than EVENTS
    #         (all data are needed by srctool)
    #center = SkyCoord(ra_cen,dec_cen,unit="deg")
    t_min =None
    t_max= None
    for indexextens in range(len(hdul1)):
        if indexextens>0:
            data1=hdul1[indexextens].data
            #newdata=new_hdul[indexextens].data
            if indexextens==1 :#was a filtering on all extensions with   'RA' in lstcolnames:
                print (f"extract extension {hdul1[indexextens].name} on RA DEC")
                print (f"         t_min {t_min} t_max {t_max}")
                #here check for each line RA DEC and if in region add it to
                #the new fitstable extension
                tmpskycoord=SkyCoord(data1['RA'],data1['DEC'],unit="deg")
                #mask=(center.separation(tmpskycoord).degree <=radius)
                mask=(SkyCoord(ra_cen,dec_cen,unit="deg").separation(tmpskycoord).degree <=radius)
                new_hdul[indexextens].data=data1[mask]
            else:
                new_hdul[indexextens].data=data1
            #lstcolnames=data1.names
            #if 'RA' in lstcolnames:
            #    print (f"extract extension {hdul1[indexextens].name} on RA DEC")
            #    print (f"         t_min {t_min} t_max {t_max}")
            #    #here check for each line RA DEC and if in region add it to
            #    #the new fitstable extension
            #    tmpskycoord=SkyCoord(data1['RA'],data1['DEC'],unit="deg")
            #    #mask=(center.separation(tmpskycoord).degree <=radius)
            #    mask=(SkyCoord(ra_cen,dec_cen,unit="deg").separation(tmpskycoord).degree <=radius)
            #    new_hdul[indexextens].data=data1[mask]
            #else:
            #    new_hdul[indexextens].data=data1
            #    if 'TIME' in lstcolnames and t_min is not None:
            #        #print (f"extract extension {hdul1[indexextens].name} on TIME")
            #        #print (f"         t_min {t_min} t_max {t_max}")
            #        mask2=((data1['TIME']>= t_min ) & (data1['TIME']<=t_max ))
            #        ## TODO : keep info from the past
            #        new_hdul[indexextens].data=data1[mask2]

        if indexextens==1:
            #only works if event extension is the 1st extension
            # after the main header
            t_min = min(new_hdul[indexextens].data['TIME'])
            t_max = max(new_hdul[indexextens].data['TIME'])

    #here must compute the image if immage=yes
    # and update wcs in primary header
    if image:
        new_hdul= setimagewithwcs(new_hdul ,ra_cen,dec_cen,radius,dimimg,rebin,overlap)

    #resultfilename='%s/%s_%s.fits'%(os.path.basename(sourcefilename).
    #split('.')[0],timestampe)
    new_hdul.writeto(resultfilename)

    hdul1.close()
    return True


def setimagewithwcs(fitsreference,racen,deccen,radius,dimmig,rebin,overlap):
    """
    create the image in the primary header and update the wcs section of the header
    Input :
        fitsreference : astropy.io.fits.hdu.hdulist.HDUList reference to the new hdulist
        racen         : central ra coordinate of the region to extract in degree
        deccen        : central dec coordinate of the region to extract in degree
        radius        : radius of the region to extract in degree
        dimmig        : dimension of the imag in pixel
        rebin         : TBD
        overlap       : TBD
    """
    #
    #evtool eventfiles='%s' outfile='%s' region='%s' image='yes' size='%s' rebin='80'  overlap='2.0'
    #python=fits.open('~/Erosita/datool_pythonV2/testspython/res12238_75/2
    #/event_python1679390257_1679478586.470288.fits')

    #this%delta(:) = this%delta(:) * rebin
    #this%refpix(1) = this%refpix(1) + resolution(1) * rebin(1) / 2.0_r8
    #this%refpix(2) = this%refpix(2) + resolution(2) * rebin(2) / 2.0_r8
    #! 0.5 / rebin to account for rebin
    #this%refpix(1) = this%refpix(1) / rebin(1) - 0.5_r8 + 0.5_r8 / rebin(1)
    #this%refpix(2) = this%refpix(2) / rebin(2) - 0.5_r8 + 0.5_r8 / rebin(2)
    #this%cdelt1p = rebin(1)
    #this%cdelt2p = rebin(2)
    #fitsreference[0].header
    newimg=np.zeros((dimmig, dimmig,))
    fitsreference[0].data=newimg
    radat     =list(fitsreference[1].data['ra'] )
    decdat    =list(fitsreference[1].data['dec  '])
    locw = wcs.WCS(fitsreference[0].header)

    locw.wcs.ctype=['RA---SIN' , 'DEC--SIN']
    locw.wcs.crval=[racen,deccen]
    locw.wcs.crpix=[dimmig/2+0.5,dimmig/2+0.5]
    cdelt=2*1.15*radius/dimmig
    # 2.*1.15*.75/1552.

    locw.wcs.cdelt=[-1*cdelt,cdelt ]
    ##x, y = w.all_world2pix(radat, decdat, 1)
    locx, locy  = locw.all_world2pix(radat, decdat, 1)
    for i in range(len(locx)):
        #print ("i=%s ,locx[i]=%s,locy[i]=%s, DIMIG=%s" %(i,locx[i],locy[i],dimmig))
        if locx[i]< dimmig and locy [i]<dimmig:
            newimg[int(locx[i]),int(locy[i])]+=1
    fitsreference[0].data=newimg
    fitsreference[0].header.update(locw.to_header())

    return fitsreference

def evtool_main(filename,outfolder,ra_cen,dec_cen,radius,image,dimimg,rebin,overlap):
    """
     Python version of main evtool.
     Used for the region extraciotn in DATOOL:
      -in this version the programme should keep event from the region given in
       parameeter and only them . plus the extension
     - all extension must be filtered on the region .
     - header must probably be updated.
     - Needs some parameters
         - Eventfile name
         - Region
     return : name of the output file extracted on the region
    """

    #define center of region as a skycoord
    #centerreg = SkyCoord(ra_cen*u.deg,dec_cen*u.deg, frame='icrs')

    #open file
    source=fits.open(filename)
    # check the region is IN the file
    secondhdu=source[1].data
    #if ( ra_cen> max(source[1].data['RA']) or ra_cen<   min(source[1].data['RA'])):
    if ( ra_cen> max(secondhdu['RA']) or ra_cen<   min(secondhdu['RA'])):

        print (" RA is out of source file, check the coordinate of the region ")
        print ("RAmax ", max(secondhdu['RA']) ,"min RA", min(secondhdu['RA']))
        #to fit in the source file
        return "ERROR: RA is out of source file, check the coordinate of the region "
        #to fit in the source file
    #if ( dec_cen> max(source[1].data['DEC']) or dec_cen< min(source[1].data['DEC'])):
    if ( dec_cen> max(secondhdu['DEC']) or dec_cen< min(secondhdu['DEC'])):
        print (" DEC is out of source file, check the coordinate of the region ")
        print ("DEC min",max(secondhdu['DEC']) ,"DEC min ",min(secondhdu['DEC']))
        #to fit in the source file
        return "ERROR: DEC is out of source file, check the coordinate of the region"
        # to fit in the source fil
    source.close()
    timestampe=datetime.now().timestamp()

    #resultfilename='%s/%s_%s.fits'%(outfolder,os.path.basename(filename).split('.')[0],timestampe)
    firstpart=os.path.basename(filename).split('.')[0]
    resultfilename=f"{outfolder}/{firstpart}_{timestampe}.fits"
    #june27 remove the timestamp for realscale test
    resultfilename=f"{outfolder}/{firstpart}.fits"

    resultfinal=extractonregion(filename,
                                resultfilename,
                                ra_cen,
                                dec_cen,
                                radius,
                                image,
                                dimimg,
                                rebin,
                                overlap)

    #loop on extension then on data
    # do it in a function

    source.close()
    if resultfinal:
        return "SUCCES"
    return "UNEXPECTED ERROR"



def help_err_message(scriptname):
    """ function to display the help message """
    defaultfilename='./eb01_123051_020_EventList_c020.fits'
    basescriptname=os.path.basename(scriptname)
    print (f"{scriptname} python version of radec2xy taks of the esass")
    print ("by A.Gueguen")
    print ("Current version : 2023")
    print ("This script must be called with 1 or 2 parameters")
    #print (f"{os.path.basename(scriptname)} --help (or -H), print this message ")
    print (f"{basescriptname}  --help (or -H), print this message ")
    print ("")
    print ("\t*\tFirst  parameter, sourcefile, is the input file, with absolute or relative path")
    print ("\t*\tSecond parameter, ra0, is the RA coordinate [deg] of region")
    print ("\t*\tThird  parameter, dec0, is the DEC coordinate [deg] of region")
    print ("\t*\tFourth parameter, Radius, is the Radius of the region [deg] for extracted region")
    print ("\t*\tFifth  parameter, resultfolder, output folder name")
    print ("\t*\t       evtool will create a new file with a new name in the output folder file")
    print ("\t*\tsixth  parameter, Optional, image boolean ( or yes no ) To activate ")
    print ("\t*\t       the creation of an image in the primary header of the result file")
    print ("\t*\tseventh parameter, Optional, size integer mandatory if image is true")

    print ("")
    print (" EXAMPLES")
    print ("")
    #print (f"./{os.path.basename(scriptname)}
    # sourcefile={defaultfilename} ra0=122 dec0=34 result=./")
    print (f"./{basescriptname} sourcefile={defaultfilename} ra0=122 dec0=34 result=./")
    print (f"    This call will run Pyradec2xy on a file {defaultfilename} ")
    print ("    and save the results in the current folder ")
    print ("")
    #print (f"If one or more parameters are missing {os.path.basename(scriptname)}
    #will ask for these info from the command line")
    print (f"If one or more parameters are missing {basescriptname} will ask for these info from the command line")
    print ("This program needs python3")

def check_input(args):#,dictparam):
    """
     check the arguments of thie script (in var args) and fill/ expand
     a dictionnary with the entries
     return none if an error occurs
    """
    dictparam={}
    if len (args)>1:
        if  args[1].lower().strip()=="--help"  or args[1].lower().strip()=="-h":
            help_err_message(args[0])
            return None
        #fill a dictionnary of the arguments
        for elem in args[1:]:
            if '=' in elem :
                #spltparam=elem.split('=')
                #dictparam[spltparam[0]]=spltparam[1]
                dictparam[elem.split('=')[0]]=elem.split('=')[1]
    #check all needed params are given, if not ask for them
    #neededkeys=['sourcefile','ra0','dec0','Radius','outfolder']
    print (dictparam)
    print ("-------")
    print (dictparam.keys())
    print ("-------")
    for crtkey in neededkeys:
        if crtkey not in dictparam.keys():
            #var = input(f"Please enter value for {crtkey}>")
            #dictparam[crtkey]=var
            dictparam[crtkey]= input(f"Please enter value for {crtkey}>")

    infile=dictparam['sourcefile']
    if not os.path.isfile(infile):
        lclstr=f"ERROR: input file {infile}\n\tcan't be found,"
        lclstr=f"{lclstr} please check it exists, or that you have access to it"
        print (lclstr)
        return None
    outfolder=dictparam['outfolder']
    if not os.path.isdir(outfolder ):
        try:
            lclstr=f"WARNING: result folder{outfolder} \n\tcan't be found, try to create it "
            print (lclstr)
            os.mkdir(outfolder)
        except Exception as errortest:
            lclstr = f"ERROR: result folder {outfolder} \n\tdoes not exist and can't be created,"
            lclstr =f'{lclstr} error was {errortest}.'
            lclstr =f'{lclstr} Please check the output folder exist or you have right to create it'
            print (lclstr)
            return  None

    #chek ra-0 is a numerical value between 0 and 360
    try:
        ra0=float(dictparam['ra0'])
    except ValueError:
        print ("the reference RA must be a float value ")
        return None
    #chek dec-0 is a numerical value between -90 and 90

    try:
        dec0=float(dictparam['dec0'])
    except ValueError:
        print ("the reference DEC must be a float value ")
        return None
    try:
        radius=float(dictparam['Radius'])
        #float(dictparam['Radius'])
    except ValueError:
        print ("the Radius must be a float value (in  degree)")
        return None

    if ra0<0 or ra0>360:
        print ("ERROR: ra0 must be greater than 0 and lower than 360 degree, given as a numerical value ")
        return None

    if dec0<-90 or dec0>90:
        print("ERROR: dec0 must be greater than -90 and lower than 90 deg, given as a numerical value")
        return None
    #print ("radius =",radius)
    imageflag=False
    if 'image' in  dictparam.keys():
        #if isinstance(dictparam['image'], (bool)):
        if dictparam['image'].lower() in ('yes','y','true','t',1):
            imageflag=True
            if 'size' not in dictparam.keys() :
                print ("ERROR: when parameter image is given the parameter size MUST be provided")
                return None
        elif dictparam['image'].lower() in ('no','n','false','f',0):
            imageflag=False
            #strflag='false'

        else:
            print ("ERROR: parameter image must be boolean  ")
            return None
    #dimimg=None
    if 'size' in  dictparam.keys():
        try :
            #dimimg=int(dictparam['size'])
            int(dictparam['size'])
        except ValueError :
            print ("ERROR: parameter size when provided must be integer")
            return None
    if not imageflag:
         dictparam['size']=0

#    if 'image' in  dictparam.keys() and 'size' not in dictparam.keys() :
#        print ("ERROR: when parameter image is given the parameter size MUST be provided")
#        return None
    return dictparam




def main(args):
    """
    Main function which call the check argument function and, if succesfull,
         launch the function which create the extracted file.
    """
    starttimesecs=datetime.now()
    #dictparam={}
    dictparam=check_input(args)#,dictparam)
    if dictparam is None:
        return None
    #if len (args)>1:
    #    if  args[1].lower().strip()=="--help"  or args[1].lower().strip()=="-h":
    #        help_err_message(args[0])
    #        return None
    #    #fill a dictionnary of the arguments
    #    for elem in args[1:]:
    #        if '=' in elem :
    #            spltparam=elem.split('=')
    #            dictparam[spltparam[0]]=spltparam[1]
    ##check all needed params are given, if not ask for them
    #for crtkey in neededkeys:
    #    if crtkey not in dictparam.keys():
    #        var = input(f"Please enter value for {crtkey}>")
    #        dictparam[crtkey]=var
    #infile=dictparam['sourcefile']
    #if not os.path.isfile(infile):
    #    lclstr=f"ERROR: input file {infile}\n\tcan't be found, please check it exists,
    # or that you have access to it"
    #    print (lclstr)
    #    return None
    #outfolder=dictparam['outfolder']
    #if not os.path.isdir(outfolder ):
    #    try:
    #        lclstr=f"WARNING: result folder{outfolder} \n\tcan't be found, try to create it "
    #        print (lclstr)
    #        os.mkdir(outfolder)
    #    except:
    #        lclstr = f"ERROR: result folder {outfolder} \n\tdoes not exist and can't be created,
    #Please check the output folder exist or you have right to create it"
    #        print (lclstr)
    #        return  None
    #
    ##chek ra-0 is a numerical value between 0 and 360
    #try:
    #    ra0=float(dictparam['ra0'])
    #except:
    #    print ("the reference RA must be a float value ")
    #    return None
    ##chek dec-0 is a numerical value between -90 and 90
    #
    #try:
    #    dec0=float(dictparam['dec0'])
    #except:
    #    print ("the reference DEC must be a float value ")
    #    return None
    #try:
    #    radius=float(dictparam['Radius'])
    #except:
    #    print ("the Radius must be a float value (in  degree)")
    #    return None
    #
    #if ra0<0 or ra0>360:
    #    print ("ERROR: ra0 must be > 0 and < 360 degree, given as a numerical value ")
    #    return None
    #
    #if dec0<-90 or dec0>90:
    #    print ("ERROR: dec0 must be > -90 and < 90 degree, given as a numerical value ")
    #    return None
    ##print ("radius =",radius)
    imageflag=False
    if 'image' in  dictparam.keys():
    #    #if isinstance(dictparam['image'], (bool)):
        if dictparam['image'].lower() in ('yes','y','true','t',1):
            imageflag=True

    #        if 'size' not in dictparam.keys() :
    #            print ("ERROR: when parameter image is given the parameter size MUST be provided")
    #            return None
    #    elif dictparam['image'].lower() in ('no','n','false','f',0):
    #        imageflag=False
    #    else:
    #        print ("ERROR: parameter image must be boolean  ")
    #        return None
    #dimimg=None
    #if 'size' in  dictparam.keys():
    #    try :
    #        dimimg=int(dictparam['size'])
    #    except :
    #        print ("ERROR: parameter size when provided must be integer")
    #        return None


#    if 'image' in  dictparam.keys() and 'size' not in dictparam.keys() :
#        print ("ERROR: when parameter image is given the parameter size MUST be provided")
#        return None
    rebin=80
    overlap=2
    #evtool_main(infile ,outfolder,ra0,dec0,radius,imageflag,dimimg,rebin,overlap)
    evtool_main(dictparam['sourcefile'],
                dictparam['outfolder'] ,
                float(dictparam['ra0']),
                float(dictparam['dec0']),
                float(dictparam['Radius']),
                imageflag,int(dictparam['size']),
                rebin,overlap)
    #radec2xy(infile,ra0,dec0)

    endtimesec=datetime.now()
    end_time = endtimesec.strftime("%H:%M:%S")
    #starttimesecs=datetime.now()
    #start_time = starttimesecs.strftime("%H:%M:%S")
    lclstr=f"DONE: region extraction finished at {end_time} took: {(endtimesec-starttimesecs)} sec"
    print(lclstr)
    return endtimesec-starttimesecs

if __name__ == "__main__":
    # execute only if run as a script
    main(sys.argv)
