'''
ARCHIVAL - script to cross-check whether a transient source
           coincides with a galaxy, star, or a previously known
           transient source.

Written for Python 3.x

Requirements: pandas, requests, astropy, astroquery,
              yattag, matplotlib, numpy, pillow

Run from the folder with the input_list as:

> python /path_to_script/archival.py -i input_list.txt [-p] [-d]

input_list of the form (Name RA DEC [DATE]):
GWT1	03:55:50.1	-45:06:20.0	2019-01-06.70
GWT2	04:11:34.02	-09:07:21.2
GWT3	03:51:59.81	-08:51:18.5
GWT4	03:49:55.34	-14:01:55.8
GWT5	04:22:17.87	-05:28:05.41

date is necessary for the search of coincident Solar system objects.

-p: parameter file
-d: luminosity distance and error in Mpc (eg. -d 80 20)

Updates:
05. 02. 2019
		- fix several bugs, TNS parsing

11. 01. 2019
		- DSS, terminal output

10. 01. 2019
		- update glade, web design

16. 12. 2018
		- script written
'''

from optparse import OptionParser
import datetime, sys
import urllib.request
import shutil, os, warnings, re
import pandas as pd
from astroquery.simbad import Simbad
from astroquery.vizier import Vizier
import astropy.coordinates as coord
import astropy.units as u
from astropy.coordinates import Angle
from astropy.visualization import PercentileInterval, AsinhStretch
from astropy.io import fits
from astropy.wcs import WCS
from astropy.table import Table
from astropy.cosmology import FlatLambdaCDM
from yattag import Doc
import numpy as np
import requests
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrowPatch
from scipy import interpolate
from PIL import Image
from io import BytesIO
import matplotlib.patheffects as PathEffects

'''
Coordinates to degrees
'''
def transCoord(ra,dec):
	rat = ra.split(':')
	dect = dec.split(':')
	rac = float(rat[0])*15 + float(rat[1])/4 + float(rat[2])/240
	if (float(dect[0]) < 0):
		decc = float(dect[0]) - float(dect[1])/60 - float(dect[2])/3600
	else:
		decc = float(dect[0]) + float(dect[1])/60 + float(dect[2])/3600
	return rac, decc

'''
Get redshift from luminosity distance
'''
def redshift(ld):
	#
	# even though we are at low z, let's use the real model
	#
	cosmo = FlatLambdaCDM(H0=67.3, Om0=0.315)
	zz=np.linspace(0,0.1,100)
	ldd=cosmo.luminosity_distance(zz)
	f=interpolate.interp1d(ldd,zz)
	
	z1=f(ld[0])
	z0=f(ld[0]-ld[1])
	z2=f(ld[0]+ld[1])
	
	ang=cosmo.angular_diameter_distance([z0,z1,z2]).value
	
	return [z0,z1,z2], ang

'''
Read the user-provided parameter file
'''
def readparfile(pat):
	print('Reading user-provided parameters...')
	try:
		par=open(pat).readlines()
	except IOError:
		print('error: There was an error opening the parameter file! Please check if you are in the right folder.')
		sys.exit()
	pard={}
	for lines in par:
		if lines[0]=='#':
			continue
		line=lines.strip('\n')
		if line[:4]=='PAR1':
			pard['radS']=float(line[5:])
		if line[:4]=='PAR2':
			pard['radT']=float(line[5:])
		if line[:4]=='PAR3':
			pard['dateT']=line[5:]
		if line[:4]=='PAR4':
			pard['radM']=float(line[5:])
		if line[:4]=='PAR5':
			pard['magM']=float(line[5:])
		if line[:4]=='PAR6':
			pard['radG']=float(line[5:])
		if line[:4]=='PAR7':
			pard['offG']=float(line[5:])
		if line[:4]=='PAR8':
			pard['web']=line[5:]
		if line[:4]=='PAR9':
			pard['fol']=line[5:]
	
	return pard


'''
This piece of code is taken from: https://ps1images.stsci.edu/ps1image.html
'''
def getimages(ra,dec,size=1200,filters="grizy"):
    
    """Query ps1filenames.py service to get a list of images
    
    ra, dec = position in degrees
    size = image size in pixels (0.25 arcsec/pixel)
    filters = string with filters to include
    Returns a table with the results
    """
    
    service = "https://ps1images.stsci.edu/cgi-bin/ps1filenames.py"
    url = ("{service}?ra={ra}&dec={dec}&size={size}&format=fits"
           "&filters={filters}").format(**locals())
    table = Table.read(url, format='ascii')
    return table

'''
This piece of code is taken from: https://ps1images.stsci.edu/ps1image.html
'''
def geturl(ra, dec, size=1200, output_size=None, filters="grizy", format="jpg", color=False):
    
    """Get URL for images in the table
    
    ra, dec = position in degrees
    size = extracted image size in pixels (0.25 arcsec/pixel)
    output_size = output (display) image size in pixels (default = size).
                  output_size has no effect for fits format images.
    filters = string with filters to include
    format = data format (options are "jpg", "png" or "fits")
    color = if True, creates a color image (only for jpg or png format).
            Default is return a list of URLs for single-filter grayscale images.
    Returns a string with the URL
    """
    
    if color and format == "fits":
        raise ValueError("color images are available only for jpg or png formats")
    if format not in ("jpg","png","fits"):
        raise ValueError("format must be one of jpg, png, fits")
    table = getimages(ra,dec,size=size,filters=filters)
    url = ("https://ps1images.stsci.edu/cgi-bin/fitscut.cgi?"
           "ra={ra}&dec={dec}&size={size}&format={format}").format(**locals())
    if output_size:
        url = url + "&output_size={}".format(output_size)
    # sort filters from red to blue
    flist = ["yzirg".find(x) for x in table['filter']]
    table = table[np.argsort(flist)]
    if color:
        if len(table) > 3:
            # pick 3 filters
            table = table[[0,len(table)//2,len(table)-1]]
        for i, param in enumerate(["red","green","blue"]):
            url = url + "&{}={}".format(param,table['filename'][i])
    else:
        urlbase = url + "&red="
        url = []
        for filename in table['filename']:
            url.append(urlbase+filename)
    return url

'''
Specifically for Paranal.
'''
def visibility(name,ra,dec,fol):
	
	#
	# sexagesimal/degrees transformation
	#
	ra1,dec1=transCoord(ra,dec)
	ra1s=str(round(ra1,6))
	dec1s=str(round(dec1,6))
	
	#
	# transformation to galactic coords
	#
	c_icrs = coord.SkyCoord(ra1*u.degree, dec1*u.degree, frame='icrs')
	gal=c_icrs.galactic.to_string('decimal',precision=6).split()
	
	#
	# current date/time
	#
	now = datetime.datetime.now()
	
	nday=str(now.day)
	if len(nday) == 1:
		nday="0"+nday
	
	nmonth=str(now.month)
	if len(nmonth) == 1:
		nmonth="0"+nmonth
	
	nyear=str(now.year)[2:]
	
	coordlist=name+'+'+ra+'+'+dec
	
	#
	# retrieving and saving the visibility plot
	#
	link="http://catserver.ing.iac.es/staralt/index.php?action=showImage&form%5Bmode%5D=1&form%5Bday%5D="+nday+"&form%5Bmonth%5D="+nmonth+"&form%5Byear%5D="+str(now.year)+"&form%5Bobservatory%5D=Cerro+Paranal+Observatory+(Chile)&form%5Bcoordlist%5D="+coordlist+"&form%5Bcoordfile%5D=&form%5Bparamdist%5D=2&form%5Bformat%5D=gif&submit=+Retrieve+"
	
	outf=fol+'/'+name+'/'+name+'.gif'
	
	try:
		with urllib.request.urlopen(link) as response, open(outf, 'wb') as outf:
			shutil.copyfileobj(response, outf)
		
		print('\t... plot saved to:', fol+'/'+name+'/'+name+'.gif')
		return {"out":name+'/'+name+'.gif',"coord":[ra,dec],"coordD":[ra1s,dec1s],"coordG":[gal[0],gal[1]]}
		
	except Exception as e:
		print(str(e))
		return {"out":"-99","coord":[ra,dec],"coordD":[ra1s,dec1s],"coordG":[gal[0],gal[1]]} 

'''
From NASA/IPAC Extragalactic Database 
'''
def extinction(name,ra,dec,fol):
	
	#
	# prepare coordinates
	#
	ra1=ra.split(':')
	dec1=dec.split(':')
	coordlist='lon='+ra1[0]+'%3A'+ra1[1]+'%3A'+ra1[2]+'&lat='+dec1[0]+'%3A'+dec1[1]+'%3A'+dec1[2]
	
	#
	# retrieving data and saving them to ascii file
	#
	link="https://ned.ipac.caltech.edu/cgi-bin/calc?in_csys=Equatorial&in_equinox=J2000.0&obs_epoch=2000.0&"+coordlist+"&pa=0.0&out_csys=Equatorial&out_equinox=J2000.0"
	
	try:
		tables = pd.read_html(link)
		
		filters=['Landolt','SDSS','UKIRT']
		nt=89
		
		nn=0
		tabe=open(fol+'/'+name+'/'+name+'_extinction.txt','w')
		for i in range(nt):
			if tables[1][0][i] in filters:
				tabe.write(tables[1][0][i]+'\t'+tables[1][1][i]+'\t'+str(tables[1][2][i])+'\t'+str(tables[1][3][i])+'\n')
				if tables[1][0][i] == 'Landolt' and tables[1][1][i] == 'V':
					av=tables[1][3][i]
				nn += 1
		
		tabe.close()
		
		print('\t... Av = ', np.round(av,2))
		print('\t... extinction in additional filters saved to:',fol+'/'+name+'/'+name+'_extinction.txt')
		return av
	except Exception as e:
		print(str(e))
		return -99

'''
Query Simbad database for existing objects
'''
def simbad(name,ra,dec,fol,rad):
	ra1,dec1=transCoord(ra,dec)
	
	pagea = requests.get('http://simbad.u-strasbg.fr/simbad/')
	if pagea.status_code is not 200:
		print('error: Page not available.')
		return -99
	else:
		print('\t... Search radius: ',rad,' arcsec')
		#
		# query simbad and write the results to a dictionary
		#
		with warnings.catch_warnings():
			warnings.simplefilter("ignore")
			result = Simbad.query_region(coord.SkyCoord(ra1, dec1, unit=(u.deg, u.deg), frame='icrs'),radius=Angle(rad,"arcsec"))
		
		try:
			tblA=result['MAIN_ID']
			n=len(tblA)
			raS=result['RA']
			decS=result['DEC']
			dist=[]
			tbla=[]
			for i in range(n):
				ra2,dec2=transCoord(":".join(raS[i].split()),":".join(decS[i].split()))
				c1=coord.SkyCoord(ra1,dec1,unit='deg')
				c2=coord.SkyCoord(ra2,dec2,unit='deg')
				sep=c1.separation(c2)
				dist.append(sep.arcsecond)
				tbla.append(tblA[i].decode("utf-8") )
			
			print('\t... Number of sources: ',str(n))
			tabb={"sname":tbla,"sdist":dist}
		except:
			print('\t... Number of sources: 0')
			tabb={"sname":None,"sdist":None}
		
		return tabb

'''
Query Transient Name Server for existing transients
'''
def weizmann(name,ra,dec,fol,rad,datestart):
	
	now=datetime.datetime.now()
	now=now.strftime("%Y-%m-%d")
	ra1=ra.replace(':','%3A')
	dec1=dec.replace(':','%3A')
	
	pagea = requests.get('https://wis-tns.weizmann.ac.il/search')
	if pagea.status_code is not 200:
		print('error: Page not available.')
		return -99
	
	else:
		link="https://wis-tns.weizmann.ac.il/search?&name=&name_like=0&isTNS_AT=all&public=all&unclassified_at=0&classified_sne=0&ra="+ra1+"&decl="+dec1+ "&radius="+str(rad)+"&coords_unit=arcsec&groupid%5B%5D=null&classifier_groupid%5B%5D=null&type%5B%5D=null&date_start%5Bdate%5D="+datestart+"&date_end%5Bdate%5D="+now+"&discovery_mag_min=&discovery_mag_max=&internal_name=&redshift_min=&redshift_max=&spectra_count=&discoverer=&classifier=&discovery_instrument%5B%5D=&classification_instrument%5B%5D=&hostname=&associated_groups%5B%5D=null&ext_catid=&num_page=50&display%5Bredshift%5D=1&display%5Bhostname%5D=1&display%5Bhost_redshift%5D=1&display%5Bsource_group_name%5D=1&display%5Bclassifying_source_group_name%5D=1&display%5Bdiscovering_instrument_name%5D=0&display%5Bclassifing_instrument_name%5D=0&display%5Bprograms_name%5D=0&display%5Binternal_name%5D=1&display%5BisTNS_AT%5D=0&display%5Bpublic%5D=1&display%5Bend_pop_period%5D=0&display%5Bspectra_count%5D=1&display%5Bdiscoverymag%5D=1&display%5Bdiscmagfilter%5D=1&display%5Bdiscoverydate%5D=1&display%5Bdiscoverer%5D=1&display%5Bsources%5D=0&display%5Bbibcode%5D=0&display%5Bext_catalogs%5D=0"
		
		print('\t... Search radius:',rad,' arcsec')
		print('\t... Start date:',datestart)
		
		#a silly piece of code: do this because the table is not easily parsable :(
		fp = urllib.request.urlopen(link)
		mybytes = fp.read()

		mystr = mybytes.decode("utf8")
		fp.close()
		pos=mystr.find('class=\"count rsults\"')
		if pos>0:
			
			pos=int(mystr.find('out of <em class="placeholder">') + len('out of <em class="placeholder">'))
			mystr1=mystr[pos:]
			pos=mystr1.find('<')
			numobj=int(mystr1[:pos])
		
			print('\t... Number of sources:',str(numobj))
			
			#search for each event by class="cell-name"
			
			tabb={"sname":[],"sra":[],"sdec":[],"sdist":[],"slink":[]}
			
			iii=-1
			for m in re.finditer('class=\"cell-name\">',mystr):
				iii += 1
				if iii==0:	#the first occurrence is skipped
					continue 
				mystr1=mystr[m.end()+9:]
				lin=mystr1[:mystr1.find('\"')]
				nam=mystr1[mystr1.find('\"')+2:mystr1.find('</a>')]
				mystr2=mystr1[mystr1.find('class=\"cell-ra\">')+len('class=\"cell-ra\">'):]
				raS=mystr2[:mystr2.find('</td>')]
				mystr3=mystr2[mystr2.find('class=\"cell-decl\">')+len('class=\"cell-decl\">'):]
				decS=mystr3[:mystr3.find('</td>')]
				
				tabb["sname"].append(nam)
				tabb["sra"].append(raS)
				tabb["sdec"].append(decS)
				
				ra1,dec1=transCoord(ra,dec)
				raS,decS=transCoord(raS,decS)
				c1=coord.SkyCoord(ra1,dec1,unit='deg')
				c2=coord.SkyCoord(raS,decS,unit='deg')
				sep=c1.separation(c2)
				
				tabb["sdist"].append(sep.arcsecond)
				tabb["slink"].append(lin)
		
		else:
			tabb={"sname":[],"sra":[],"sdec":[],"sdist":[],"slink":[]}
			print('\t... Number of sources: 0')

		return tabb

'''
Query for minor planets
'''
def mplanet(name,ra,dec,tim,rad,mag,fol):
	
	#
	# handle time
	#
	nyear=tim.split('-')[0]
	nmonth=tim.split('-')[1]
	nday=tim.split('-')[2]

	#
	# check page
	#
	pagea = requests.get('https://cgi.minorplanetcenter.net/cgi-bin/checkmp.cgi')
	if pagea.status_code is not 200:
		print('error: Page not available.')
		return -99
	
	else:
		print('\t... Search radius:',rad,' arcmin')
		print('\t... Magnitude limit:',mag,'mag')
		print('\t... Epoch:',tim)
		ra1=ra.replace(':','+')
		dec1=dec.replace(':','+')
		link="https://minorplanetcenter.net/cgi-bin/mpcheck.cgi?year="+nyear+"&month="+nmonth+"&day="+nday+"&ra="+ra1+"&decl="+dec1+"&which=pos&TextArea=&radius="+str(rad)+"&limit="+str(mag)+"&oc=500&sort=d&mot=h&tmot=s&pdes=u&needed=f&ps=n&type=p"
	
		htm=urllib.request.urlopen(link).read().decode('utf-8')
		
		
		n1=htm.find("<pre>")
		n2=htm.find("</pre>")
		
		if n1 == -1:
			print('\t... Minor planets found: 0')
			return 'No known minor planets, brighter than V = '+str(mag)+', were found in the '+str(rad)+'-arcminute region around the source on '+tim
		else:
			print('\t... Minor planets found - see html output for details!')
			return htm[n1+5:n2]
		
		
	

'''
Query the GLADE catalog directly from Vizier
'''
def glade(name,ra,dec,fol,radG,ld,offG):
	
	ra1,dec1=transCoord(ra,dec)
	#
	# luminosity distance is not provided
	#
	if ld == -1:
		print('\t... Luminosity distance is not provided - a simple search within a radius of',radG,'arcsec')
		#
		# querry catalog 
		#
		
		v=Vizier()
		v.ROW_LIMIT = -1
		with warnings.catch_warnings():
			warnings.simplefilter("ignore")
			result=v.query_region(coord.SkyCoord(ra1, dec1, unit=(u.deg, u.deg), frame='icrs'), catalog='VII/281/glade2',radius=Angle(radG,"arcsec"))
		
		a=str(result).split()
		
		if a[0]=='Empty':
			print('\t... Number of sources: 0')
			tabb={"sname":[],"pars":[]}
		else:
			a=result['VII/281/glade2']
			n=len(a)
			pos=(6,7,8,9,10,11,12,13,14,15,16)
			nameF=['ra','dec','dist','red','Bmag','MB','Jmag','Hmag','Kmag','flag2','flag3']
			tabG=np.zeros((n,len(pos)+1))
			nameG=[]
			j=0
			with warnings.catch_warnings():
				warnings.simplefilter("ignore")
				for i in range(n):
					if a[i][5]=='Q':
						nameG.append(a[i][4])
					elif a[i][5]=='C':
						nameG.append(a[i][1])
					else:
						if a[i][3]=="---":
							nameG.append(a[i][2])
						else:
							nameG.append('2MASS J'+a[i][3])
					tabG[i,:-1]=[a[i][j] for j in pos]
					c1=coord.SkyCoord(ra1,dec1,unit='deg')
					c2=coord.SkyCoord(tabG[i,0],tabG[i,1],unit='deg')
					sep=c1.separation(c2)
					tabG[i,-1]=sep.arcsecond
				tabG[np.isnan(tabG)]=-99
		
			print('\t... Number of sources:',str(n))
			tabb={"sname":nameG,"pars":tabG}
		
		return tabb
	#
	# luminosity distance is provided
 	#
	else:
		red1,ang=redshift(ld)
		print('\t... Luminosity distance:',str(ld[0]),'+-',str(ld[1]))
		#
		# compute the approximate region assuming a physical offset
		# note that the lower lumd distance is assumed for calculation
		# in order not to underestimate the region.
		#
		reg=int((offG/(1000.*ang[0]))*3600*180./np.pi) # [arcsec]
		print('\t... Search radius:',str(int(reg)),' arcsec (',str(offG),' kpc at the lower luminosity distance limit)')
		
		#
		# query catalog
		#
		v=Vizier()
		v.ROW_LIMIT = -1
		with warnings.catch_warnings():
			warnings.simplefilter("ignore")
			result=v.query_region(coord.SkyCoord(ra1, dec1, unit=(u.deg, u.deg), frame='icrs'), catalog='VII/281/glade2',radius=Angle(reg,"arcsec"))
		
		a=str(result).split()
		
		if a[0]=='Empty':
			print('\t... Number of sources: 0')
			tabb={"sname":[],"pars":[]}
		else:
			a=result['VII/281/glade2']
			n=len(a)
			pos=(6,7,8,9,10,11,12,13,14,15,16)
			nameF=['ra','dec','dist','red','Bmag','MB','Jmag','Hmag','Kmag','flag2','flag3']
			rrr=[2,3,4,5,6,7,8]
			tabG=np.zeros((n,len(pos)+4))
			nameG=[]
			j=0
			with warnings.catch_warnings():
				warnings.simplefilter("ignore")
				for i in range(n):
					if a[i][5]=='Q':
						nameG.append(a[i][4])
					elif a[i][5]=='C':
						nameG.append(a[i][1])
					else:
						if a[i][3]=="---":
							nameG.append(a[i][2])
						else:
							nameG.append('2MASS J'+a[i][3])
					tabG[i,:-4]=[a[i][j] for j in pos]
					for k in range(len(pos)):
						if k in rrr:
							tabG[i,k]=round(tabG[i,k],2)
					c1=coord.SkyCoord(ra1,dec1,unit='deg')
					c2=coord.SkyCoord(tabG[i,0],tabG[i,1],unit='deg')
					sep=c1.separation(c2)
					tabG[i,-4]=np.round(sep.arcsecond,2)	#projected seperation
					tabG[i,-3:]=np.round(sep.arcsecond*(np.pi/(3600*180))*ang*1000,2) #projected seperation in kpc
			
			print('\t... Number of sources:',str(n))
			tabb={"sname":nameG,"pars":tabG}
		
		return tabb

'''
Get a cutout of the region from PanSTARRS survey and create a simple FC
'''
def panscut(name,ra,dec,fol):
	
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')
	
	#
	# basic parameters, can be changed in principle
	#
	size=1200 #==5 arcmin: the size [px] of the retrieved and saved image
	size1=240 #=1 arcmin: the size [px] of the central part of the image for the FC
	fil='r' #filters
	
	#
	# get the cutout image from PanSTARRS
	#
	pix=0.25
	ra1,dec1=transCoord(ra,dec)
	fitsurl = geturl(ra1, dec1, size=size, filters=fil, format="fits")
	
	if not fitsurl:	#no image returned, not in the field
		print('\t... Image not available')
		return -99
	else:
		fh = fits.open(fitsurl[0])
		fim = fh[0].data
		fhe = fh[0].header
		
		#
		# save fits image
		#
		with warnings.catch_warnings(): #to prevent the warnings in case the script is run several times and the files are rewritten
			fits.writeto(fol+'/'+name+'/'+name+'_panstarrs_cutout_'+str(size)+'px_'+fil+'.fits', fim, fhe,overwrite=True)
		
		#
		# make the cut and apply scale
		#
		x1=int(size/2. - size1/2.)
		x2=int(size/2. + size1/2.)
		y1=int(size/2. - size1/2.)
		y2=int(size/2. + size1/2.)
		
		fim=fim[y1:y2,x1:x2]
		fhe["CRPIX1"]=fhe["CRPIX1"] - x1
		fhe["CRPIX2"]=fhe["CRPIX2"] - y1
		
		fim[np.isnan(fim)] = 0.0
		transform = AsinhStretch() + PercentileInterval(99.9)
		bfim = transform(fim)
		
		with warnings.catch_warnings():	#because there are deprecated keywords in the header, no need to write it out
			warnings.simplefilter("ignore")
			wcs = WCS(fhe)
		
		#
		# produce and save the FC
		#
		fig=plt.figure(1, figsize=(12,6))
		fig1=fig.add_subplot(121,aspect='equal', projection=wcs)
		
		plt.imshow(bfim,cmap='gray_r',origin='lower')
		c = Circle((ra1, dec1), 0.00028, edgecolor='k', lw=2,facecolor='none',transform=fig1.get_transform('fk5'))
		fig1.add_patch(c)
		c = Circle((ra1, dec1), 0.00028, edgecolor='yellow', facecolor='none',transform=fig1.get_transform('fk5'))
		fig1.add_patch(c)
		txta=fig.text(0.14,0.8,name,fontsize=23,color='yellow')
		txta.set_path_effects([PathEffects.withStroke(linewidth=2, foreground='k')])
		txtb=fig.text(0.32,0.8,'Pan-STARRS',fontsize=23,color='yellow')
		txtb.set_path_effects([PathEffects.withStroke(linewidth=2, foreground='k')])
		
		fig1.add_patch(FancyArrowPatch((size1-10/0.25-10,30),(size1-10,30),arrowstyle='-',color='k',linewidth=3.5))
		fig1.add_patch(FancyArrowPatch((size1-10/0.25-10,30),(size1-10,30),arrowstyle='-',color='yellow',linewidth=2.0))

		txtc=fig.text(0.42,0.18,'10"',fontsize=20,color='yellow')
		txtc.set_path_effects([PathEffects.withStroke(linewidth=2, foreground='k')])
		
		fname=name+'_panstarrs.png'
		
		#plt.savefig(fol+'/'+fname,dpi=150,format='PNG',bbox_inches='tight')
		#fig.clear()
		
		#
		# get colour image
		#
		url = geturl(ra1,dec1,size=size1,filters='grz',output_size=None,format='png',color=True)
		r = requests.get(url)
		
		im = Image.open(BytesIO(r.content))
		
		#figA=plt.figure(2)
		figA1=fig.add_subplot(122,aspect='equal',projection=wcs)
		
		plt.imshow(im,origin="lower")
		c = Circle((ra1, dec1), 0.00028, lw=2, edgecolor='k', facecolor='none',transform=figA1.get_transform('fk5'))
		figA1.add_patch(c)
		c = Circle((ra1, dec1), 0.00028, edgecolor='yellow', facecolor='none',transform=figA1.get_transform('fk5'))
		figA1.add_patch(c)
		#fname=name+'_pan_'+'color_grz.png'
		
		fig.text(0.5,0.05,r'$\alpha$',fontsize=22,horizontalalignment='center')
		fig.text(0.03,0.5,r'$\delta$',fontsize=22,verticalalignment='center',rotation=90)
		
		plt.savefig(fol+'/'+name+'/'+fname,dpi=120,format='PNG')
		fig.clear()
		
		return fname

'''
Get a cutout of the region from DSS survey and create a simple FC
To be checked in case the image is not available in PanSTARRS
'''
def dsscut(name,ra,dec,fol):
	
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')

	#
	# basic parameters, can be changed in principle
	#
	size=5 # [arcmin] the size of the retrieved and saved image
	size1=2 # [arcmin] the size of the image used for the FC
	pix=1.008 # approx pixel scale

	#
	# get the cutout image from the ESO archive
	#
	ra1,dec1=transCoord(ra,dec)
	raS=ra.replace(':','%3A')
	decS=dec.replace(':','%3A')
	link='http://archive.eso.org/dss/dss/image?ra='+raS+'&dec='+decS+'&equinox=J2000&name=&x='+str(size)+'&y='+str(size)+'&Sky-Survey=DSS2-red&mime-type=download-fits&statsmode=WEBFORM'
	
	outf=fol+'/'+name+'/'+name+'_dss_red.fits'
	
	try:
		with urllib.request.urlopen(link) as response, open(outf, 'wb') as outf:
			shutil.copyfileobj(response, outf)
		
		print('\t... image saved to:', fol+'/'+name+'/'+name+'_dss_red.fits')
		
		#
		# load image and make a FC
		#
		fh=fits.open(fol+'/'+name+'/'+name+'_dss_red.fits')
		fim = fh[0].data
		fhe = fh[0].header
		
		#
		# cut image and apply scale
		# 
		x1=int(30*(size - size1))
		x2=int(30*(size + size1))
		y1=int(30*(size - size1))
		y2=int(30*(size + size1))
		
		fim=fim[y1:y2,x1:x2]
		fhe["CRPIX1"]=fhe["CRPIX1"] - x1
		fhe["CRPIX2"]=fhe["CRPIX2"] - y1
		
		fim[np.isnan(fim)] = 0.0
		transform = AsinhStretch() + PercentileInterval(99.7)
		bfim = transform(fim)
		
		with warnings.catch_warnings():	#because there are deprecated keywords in the header, no need to write it out
			warnings.simplefilter("ignore")
			wcs = WCS(fhe)
		
		#
		# produce and save the FC
		#
		fig=plt.figure(2)
		fig1=fig.add_subplot(111,aspect='equal', projection=wcs)
		
		plt.imshow(bfim,cmap='gray_r',origin='lower')
		c = Circle((size1*30, size1*30), pix*2, edgecolor='k', lw=2,facecolor='none')
		fig1.add_patch(c)
		c = Circle((size1*30, size1*30), pix*2, edgecolor='yellow', facecolor='none')
		fig1.add_patch(c)
		txta=fig.text(0.25,0.8,name,fontsize=23,color='yellow')
		txta.set_path_effects([PathEffects.withStroke(linewidth=2, foreground='k')])
		txtb=fig.text(0.7,0.8,'DSS',fontsize=23,color='yellow')
		txtb.set_path_effects([PathEffects.withStroke(linewidth=2, foreground='k')])
		
		fig1.add_patch(FancyArrowPatch((size1*60-15/pix-10,20),(size1*60-10,20),arrowstyle='-',color='k',linewidth=3.5))
		fig1.add_patch(FancyArrowPatch((size1*60-15/pix-10,20),(size1*60-10,20),arrowstyle='-',color='yellow',linewidth=2.0))

		txtc=fig.text(0.71,0.16,'15"',fontsize=20,color='yellow')
		txtc.set_path_effects([PathEffects.withStroke(linewidth=2, foreground='k')])
		
		plt.xlabel(r'$\alpha$',fontsize=22)
		plt.ylabel(r'$\delta$',fontsize=22)
		
		fname=name+'_dss.png'
		
		plt.savefig(fol+'/'+name+'/'+fname,dpi=120,format='PNG')
		fig.clear()
		
		return fname
		
		
	except Exception as e:
		print(str(e))
		return -99
	

'''
Create additional links/queries that are not imperative but maybe still useful
'''
def additional(name,ra,dec):
	#
	# link for ESO archive search
	#
	linkESO='http://archive.eso.org/wdb/wdb/eso/eso_archive_main/query?ra='+ra+'&dec='+dec+'&deg_or_hour=hours&box=00+10+00&max_rows_returned=2000'

	#
	# link for NED
	#
	linkNED='http://ned.ipac.caltech.edu/cgi-bin/nph-objsearch?in_csys=Equatorial&in_equinox=J2000.0&lon='+ra+'&lat='+dec+'&radius=1.0&search_type=Near+Position+Search'

	tabb={'eso':linkESO, 'ned':linkNED}
	
	return tabb


'''
Create an html summary file. 
'''
def makehtml(name,parin,vis,ext,sim,wei,mpp,gla,ld,pan,dss,add,fol):
	
	n=len(name)
	
	doc, tag, text = Doc().tagtext()
	with tag('html'):
		#
		# head, include simple internal css
		#
		with tag('head'):
			with tag('title'):
				text('GW candidates check')
			doc.asis("<!-- gsfc meta tags --><meta http-equiv=\"Content-Type\" content=\"text/html;charset=UTF-8\">")
			with tag('style'):
				doc.asis('p { font-family: Geneva,sans-serif; font-size: 1em; font-weight: 500; margin: 1em; }')
				doc.asis(' a { color: #ff6600; transition: .5s; -moz-transition: .5s; -webkit-transition: .5s; -o-transition: .5s; font-family: Geneva, sans-serif; }')
				doc.asis(' table { font-family: Geneva,sans-serif; font-size: 1em; margin: 1em; border-bottom: 1px solid black;}')
				doc.asis(' #table1 th {padding-top: 2px; padding-bottom: 2px; text-align: left; border-bottom: 1px solid black; color: black;}')
				doc.asis(' #table1 tr:nth-child(even){background-color: #f2f2f2;}')
				doc.asis(' #table1 td {padding-left: 10px; padding-right: 10px}')
				doc.asis(' #table1 td:nth-child(4) {background: #FFDEAD}')
				doc.asis(' #table1 td:nth-child(12) {background: #FFDEAD}')
		# 
		# body
		#
		with tag('body'):
			for i in range(n):
				with tag('h2'):
					text(name[i])
				#
				# coordinates, basic info
				#
				with tag('p'):
					text('ra = '+vis[i]["coord"][0]+' = '+vis[i]["coordD"][0])
					doc.asis('</br>')
					text('dec = '+vis[i]["coord"][1]+' = '+vis[i]["coordD"][1])
					doc.asis('</br></br>')
					text('Galactic (l,b) = ('+vis[i]["coordG"][0]+','+vis[i]["coordG"][1]+')')
					
				#
				# gal extinction
				#
				with tag('p'):
					with tag('b', style='border-top: 1px solid; border-bottom: 1px solid'):
						text('Galactic extinction (')
						with tag('a', href="http://ned.ipac.caltech.edu/forms/calculator.html", target="_blank"):
									text('NASA/IPAC')
						doc.asis(')</br></br>')
					
					if ext[i] == -99:
						text('Not available - see the output in the terminal for possible errors.')
					else:
						text('Av = '+str(round(ext[i],2))+' (')
						with tag('a', href=name[i]+'/'+name[i]+'_extinction.txt', target="_blank"):
							text('details')
						doc.asis(')')
				#
				# simbad query
				#
				with tag('p'):
					with tag('b', style='border-top: 1px solid; border-bottom: 1px solid'):
						text('Simbad (')
						with tag('a', href="http://simbad.u-strasbg.fr/simbad/", target="_blank"):
									text('Simbad-CDS')
						doc.asis(')</br></br>')
					
					if sim[i] == -99:
						text('Page not available.')
					else:
						text('Search radius: ' + str(parin[0]) + ' arcsec')
						doc.asis('</br>')
						doc.asis('</br>')
						if not sim[i]["sname"]:
							text('Sources found: 0')
						else:
							nn=len(sim[i]["sname"])
							
							text('Sources found: '+str(nn))
							
							tabh=['Name','sep ["]','Link']
							with tag('table', id='table1'):
								with tag('tr'):
									for ii in range(len(tabh)):
										with tag('th'):
											text(tabh[ii])
								for ii in range(nn):
									with tag('tr'):
										with tag('td'):
											text(sim[i]["sname"][ii])
										with tag('td'):
											text(str(round(sim[i]["sdist"][ii],2)))
										with tag('td'):
											with tag('a', href="http://simbad.u-strasbg.fr/simbad/sim-id?Ident="+sim[i]["sname"][ii]+"&NbIdent=1&Radius=2&Radius.unit=arcmin&submit=submit+id", target="_blank"):
												text('more')
																						
				#
				# TNS query
				#
				with tag('p'):
					with tag('b', style='border-top: 1px solid; border-bottom: 1px solid'):
						text('Transient Name Server (')
						with tag('a', href="https://wis-tns.weizmann.ac.il/search", target="_blank"):
									text('TNS')
						doc.asis(')</br></br>')
					
					if wei[i] == -99:
						text('Page not available.')
					else:
						text('Search radius: '+ str(parin[1]) + ' arcsec, since '+str(parin[2]))
						doc.asis('</br>')
						doc.asis('</br>')
						if not wei[i]["sname"]:
							text('Sources found: 0')
						else:
							nn=len(wei[i]["sname"])
							text('Sources found: '+str(nn))
							
							tabh=['Name','sep ["]','Link']
							with tag('table', id='table1'):
								with tag('tr'):
									for ii in range(len(tabh)):
										with tag('th'):
											text(tabh[ii])
								for ii in range(nn):
									with tag('tr'):
										with tag('td'):
											text(wei[i]["sname"][ii])
										with tag('td'):
											text(str(round(wei[i]["sdist"][ii],2)))
										with tag('td'):
											with tag('a', href="https://wis-tns.weizmann.ac.il/"+wei[i]["slink"][ii], target="_blank"):
												text('more')

				#
				# Minor planets querry
				#
				if mpp[i] is not 'noplanet':
					with tag('p'):
						with tag('b', style='border-top: 1px solid; border-bottom: 1px solid'):
							text('Minor planets (')
							with tag('a', href="https://cgi.minorplanetcenter.net/cgi-bin/checkmp.cgi", target="_blank"):
										text('MPChecker')
							doc.asis(')</br></br>')
						
						if mpp[i] == -99:
							text('Page not available.')
						else:
							if mpp[i][:2]=='No':
								text(mpp[i])
							else:
								with tag('pre'):
									text(mpp[i])
				
				#
				# GLADE query
				#
				with tag('p'):
					with tag('b', style='border-top: 1px solid; border-bottom: 1px solid'):
						text('GLADE catalogue (')
						with tag('a', href="http://aquarius.elte.hu/glade/", target="_blank"):
									text('GLADE home')
						doc.asis(', see also ')
						with tag('a', href="https://gxgwtest.herokuapp.com/", target="_blank"):
									text('Flask')
						doc.asis(')</br></br>')
					#
					# no lumd provided
					#
					if ld ==-1:
						if not gla[i]["sname"]:
							text('No galaxies were found within '+str(parin[5])+' arcsec radius.')
						
						else:
							nn=len(gla[i]["sname"])
							text(str(nn)+' galaxies were found within '+str(parin[5])+' arcsec radius.')
							doc.asis('</br>')
							for ii in range(nn):
								doc.asis(gla[i]["sname"][ii]+'&nbsp;&nbsp;&nbsp;&nbsp;'+str(round(gla[i]["pars"][ii][-1],2))+'  arcsec from the candidate,&nbsp;&nbsp;&nbsp;&nbsp; M_B='+str(round(gla[i]["pars"][ii][5]))+',&nbsp;&nbsp;&nbsp;&nbsp; z='+str(round(gla[i]["pars"][ii][3]))+',&nbsp;&nbsp;&nbsp;&nbsp; lumd='+str(round(gla[i]["pars"][ii][2],2))+' Mpc (')
								with tag('a', href="http://simbad.u-strasbg.fr/simbad/sim-id?Ident="+gla[i]["sname"][ii]+"&NbIdent=1&Radius=2&Radius.unit=arcmin&submit=submit+id", target="_blank"):
									text('link')
								doc.asis(')</br>')
					#
					# lumd provided
					#
					else:
						text('Input luminosity distance: ' + str(ld[0]) + ' +- ' + str(ld[1]) + ' Mpc')
						doc.asis('</br>')
						doc.asis('</br>')
						if not gla[i]["sname"]:
							text('Search radius: '+str(parin[6])+' kpc (computed at the lower distance limit)')
							doc.asis('</br>')
							doc.asis('</br>')
							text('Sources found: 0')
						else:
							tabh=['Name','ra','dec','lum dist [Mpc]','red','Bmag','MB','Jmag','Hmag','Kmag','flag2','flag3','sep ["]*','sep1*','sep2*','sep3*']
							exc=[10,11,12]
							exc1=[9,10,11]
							nn=len(gla[i]["sname"])
							text('Search radius: '+str(parin[6])+' kpc (computed at the lower distance limit)')
							doc.asis('</br>')
							doc.asis('</br>')
							text('Sources found: '+str(nn))
							with tag('table', id='table1'):
								with tag('caption', style="caption-side:bottom; text-align:left;"):
									text('*projected seperations [kpc] at the luminosity distances of (1) ' + str(ld[0]-ld[1]) + ', (2) ' + str(ld[0]) + ', (3) ' + str(ld[0] + ld[1])+' Mpc')
								with tag('tr'):
									for ii in range(len(tabh)):
										if ii in exc:
											continue
										with tag('th'):
											text(tabh[ii])
								for ii in range(nn):
									with tag('tr'):
										with tag('td'):
											text(gla[i]["sname"][ii])
										for jj in range(len(tabh)-1):
											if jj in exc1:
												continue
											with tag('td'):
												text(str(gla[i]["pars"][ii][jj]))
								

				with tag('p'):
					with tag('b', style='border-top: 1px solid; border-bottom: 1px solid'):
						text('Visibility (Paranal)')
					
					doc.asis('</br>')
					doc.asis('</br>')
					
					if vis[i]["out"] is not "-99":
						doc.stag('img', src=vis[i]["out"])
				
				with tag('p'):
					with tag('b', style='border-top: 1px solid; border-bottom: 1px solid'):
						text('Pan-STARRS')
					
					doc.asis('</br>')
					doc.asis('</br>')
					
					if pan[i] == -99:
						text('Images not available.')
					
					else:
						text('The region around the source in r and grz image (5\'x5\' image saved to '+name[i]+'/'+name[i]+'_panstarrs_cutout_1200px_r.fits):')
						
						doc.asis('</br>')
						doc.asis('</br>')
						
						doc.stag('img', src=name[i]+'/'+pan[i])
				
				
				if pan[i] == -99:
					with tag('p'):
						with tag('b', style='border-top: 1px solid; border-bottom: 1px solid'):
							text('DSS')
						
						doc.asis('</br>')
						doc.asis('</br>')
				
						if dss[i] == -99:
							text('Images not available.')
						
						else:
							text('The region around the source in r (5\'x5\' image saved to '+name[i]+'/'+name[i]+'_dss_red.fits):')
						
							doc.asis('</br>')
						doc.asis('</br>')
						
						doc.stag('img', src=name[i]+'/'+dss[i])
				
				with tag('p'):
					with tag('b', style='border-top: 1px solid; border-bottom: 1px solid'):
						text('Additional links and information')
					
					doc.asis('</br>')
					doc.asis('</br>')
					
					text('ESO archive ')
					with tag('a', href=add[i]["eso"], target="_blank"):
									text('search')
					text(' (within 10 arcmin box)')
					doc.asis('</br>')
				
					text('NED ')
					with tag('a', href=add[i]["ned"], target="_blank"):
									text('search')
					text(' (within 1 arcmin box)')
				
				doc.stag('hr', style="width: 100%")
					
	
	result = doc.getvalue()
	
	f=open(fol+'/'+parin[-1],'w')
	f.write(result)
	f.close()
	
	
def main():
	usage = "usage: %prog [options] \n\nBasic checks for GW candidates."
	parser = OptionParser(usage)
	parser.add_option("-i", "--input", action="store", nargs=1, type="string", dest="inputfile", help="Input list of sources, each line of the form: NAME RA DEC [TIME_DET]. TIME_DET should be in the form yyyy-mm-dd.dd")
	parser.add_option("-p", "--input1", action="store", nargs=1, type="string", dest="parfile", help="[optional] File with user defined parameters")
	parser.add_option("-d", "--input2", action="store", nargs=2, type="float", dest="lumd", help="[optional] luminosity distance and error in Mpc")

	(options, args) = parser.parse_args()

	print('\nARCHIVAL\n')

	try:
		inp=open(options.inputfile).readlines()
	except IOError:
		print('error: There was an error opening the file with coordinates!')
		sys.exit()
		
	n=len(inp)

	if options.parfile:
		parin=readparfile(options.parfile)
		rad=parin['radS']
		radW=parin['radT']
		datestart=parin['dateT']
		radM=parin['radM']
		magM=parin['magM']
		radG=parin['radG']
		offG=parin['offG']
		webname=parin['web']
		fol=parin['fol']
	else:
		rad=20 # [arcsec] radius in which known objects are searched for in Simbad
		radW=20 # [arcsec] radius in which preexisting transients are searched for in TNS server
		datestart='2001-01-01' # start date for searching in TNS server
		radM=5 # [arcmin] radius in which to search for potential minor planets
		magM=24.0 # limiting V magnitude of minor planets
		radG=15 # [arcsec] radius in which potential host galaxies are searched for in GLADE, if lumd unknown
		offG=100 # [kpc] upper limit for the physical offset of the source and a galaxy in GLADE, if lumd known
		webname='summary.html'
		fol='output'

	
	print('Input file and content:')
	print('\t',options.inputfile)
	for lines in inp:
		print('\t\t',lines.strip('\n'))
	
	print('Start program')
	print('Note that the execution speed depends on the servers of the pages from where the content is parsed.')
	
	if not os.path.exists(fol):
		os.mkdir(fol)
	
	if options.lumd:
		ld=options.lumd
	else:
		ld=-1
	
	nam=[]
	vis=[]
	ext=[]
	sim=[]
	wei=[]
	gla=[]
	pan=[]
	mpp=[]
	add=[]
	dss=[]
	for lines in inp:
		line=lines.strip('\n')
		lin=line.split()
		if not os.path.exists(fol+'/'+lin[0]):
			os.mkdir(fol+'/'+lin[0])
		nam.append(lin[0])
		print(lin[0])
		print('\tVisibility')
		vis.append(visibility(lin[0],lin[1],lin[2],fol))
		print('\tExtinction')
		ext.append(extinction(lin[0],lin[1],lin[2],fol))
		print('\tSimbad')
		sim.append(simbad(lin[0],lin[1],lin[2],fol,rad))
		print('\tTNS')
		wei.append(weizmann(lin[0],lin[1],lin[2],fol,radW,datestart))
		print('\tGLADE')
		gla.append(glade(lin[0],lin[1],lin[2],fol,radG,ld,offG))
		print('\tMinor planets')
		if len(lin)==4:
			mpp.append(mplanet(lin[0],lin[1],lin[2],lin[3],radM,magM,fol))
		else:
			mpp.append('noplanet')
		print('\tRetrieving Pan-STARRS cutout')
		pan.append(panscut(lin[0],lin[1],lin[2],fol))
		if pan[-1]==-99:
			print('\tRetrieving DSS cutout')
			dss.append(dsscut(lin[0],lin[1],lin[2],fol))
		else:
			dss.append(-99)
		print('\tAdditional info')
		add.append(additional(lin[0],lin[1],lin[2]))
		#break

	
	print('\nBuilding an html summary page: ', fol+'/'+webname)
	
	parin=[rad,radW,datestart,radM,magM,radG,offG,webname]
	makehtml(nam,parin,vis,ext,sim,wei,mpp,gla,ld,pan,dss,add,fol)
	
	print('\nScript completed.\n')



if __name__ == "__main__":
    main()
