#Research Techniques in Astronomy
#Coursework 2 Problem 2
#Oonagh-Alice Parker
#URN:6542728

#This code uses the MilkyWaySats.txt data provided to create two aitoff maps
#of all the positions of the satellites in the Milky Way in both RA, Dec and
#Galactic coordinates, with different shapes depending on the boundary 
#magnitude classification. 

#------------------------------------------------------------------------------

#import all the modules needed for the code
import numpy as np
import matplotlib.pyplot as plt
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy import coordinates as coord
import matplotlib as mpl

#------------------------------------------------------------------------------

#Load in all the data from the text file
RA = np.loadtxt('MilkyWaySats.txt', skiprows=1, usecols=1)
Dec = np.loadtxt('MilkyWaySats.txt', skiprows=1, usecols=2)
Mv = np.loadtxt('MilkyWaySats.txt', skiprows=1, usecols=3)
vel = np.loadtxt('MilkyWaySats.txt', skiprows=1, usecols=4)

#------------------------------------------------------------------------------

#Create dictionaries to store the ultra faints and classical sats
uf={'RA':[],'Dec':[],'Mv':[],'vel':[]}
cs={'RA':[],'Dec':[],'Mv':[],'vel':[]}


#Define the boundary magnitude
for i in range(len(Mv)):
	if Mv[i]>-7.5:
		uf['RA'].append(RA[i])
		uf['Dec'].append(Dec[i])
		uf['Mv'].append(Mv[i])
		uf['vel'].append(vel[i])
	if Mv[i]<-7.5:
		cs['RA'].append(RA[i])
		cs['Dec'].append(Dec[i])
		cs['Mv'].append(Mv[i])
		cs['vel'].append(vel[i])

#------------------------------------------------------------------------------

#Store all of the coordinates into a single array with the correct units
ccs =SkyCoord(ra=cs['RA']*u.degree, dec=cs['Dec']*u.degree, frame='icrs')
cuf =SkyCoord(ra=uf['RA']*u.degree, dec=uf['Dec']*u.degree, frame='icrs')

#Convert the coordinates to galactic coordinates
csgal=ccs.galactic
ufgal=cuf.galactic

#------------------------------------------------------------------------------

#Create an aitoff map
def aitoff(coords, fig=None, ax=None, color='k', marker='', label=''):
	sph = coords.spherical
	cs = ax.scatter(-sph.lon.wrap_at(180*u.deg).radian, sph.lat.radian, color=color, marker=marker, label=label)
	return fig, ax
def fmt_func(x, pos):
		val = coord.Angle(-x*u.radian).wrap_at(360*u.deg).degree
		return f'${val:.0f}' + r'^{\circ}$'

#------------------------------------------------------------------------------

#Create two subplots for each of the sets of data
fig,(axis1, axis2)=plt.subplots(2,1, subplot_kw=dict(projection="aitoff"))
ticker = mpl.ticker.FuncFormatter(fmt_func)
axis1.xaxis.set_major_formatter(ticker)
axis1.grid()
axis2.xaxis.set_major_formatter(ticker)
axis2.grid()

#Create RA, Dec Plot
fig, axis1 = aitoff(ccs, fig, axis1, color='b', marker='*', label='Ultra Faints')
fig, axis1 = aitoff(cuf, fig, axis1, color='r', marker='o', label='Classical Satellites')
axis1.legend()
axis1.set_xlabel('RA [degrees]')
axis1.set_ylabel('dec [degrees]')
axis1.set_title('Aitoff Map in RA, Dec Coordinates', loc='right')

#Create Galactic Coordinates Plot
fig, axis2 = aitoff(csgal, fig, axis2, color='b', marker='*', label='Ultra Faints')
fig, axis2 = aitoff(ufgal, fig, axis2, color='r', marker='o',label='Classical Satellites')
axis2.legend()
axis2.set_xlabel('l')
axis2.set_ylabel('b')
axis2.set_title('Aitoff Map in Galactic Coordinates', loc='right')
plt.show()

#------------------------------------------------------------------------------

