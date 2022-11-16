#Research Techniques in Astronomy
#Coursework 2 Problem 1
#Oonagh-Alice Parker
#URN:6542728

#This code uses the LeoV-RT data to find the mass-to-light ratio of the Leo V
#dwarf galaxies. It normalises the spectra of each of the stars and finds their
#velocities, metallicities and signal to noise ratios for each, plotting their
#spectra both normalised and unnormalised with their velocities. It also 
#calculates the mass-to-light ratio for the galaxy.

#------------------------------------------------------------------------------

#Import all the necessary modules
from astropy.io import fits
from astropy import units as u
from astropy.visualization import quantity_support
import matplotlib.pyplot as plt
from specutils import Spectrum1D, SpectralRegion
import numpy as np
from specutils.fitting import fit_generic_continuum, fit_continuum
from astropy.modeling import models, fitting
from astropy.nddata import InverseVariance
from specutils.analysis import snr_derived, equivalent_width
from specutils.manipulation import SplineInterpolatedResampler
from scipy import ndimage as ndi

#------------------------------------------------------------------------------

#Load in all the data

#Create empty lists needed to store all of the data
lam = []
flux = []
ivar = []
header = []
v = []
std = []

#Load in all the data needed from each of the spectra and append them to the 
#empty lists
ExpId = ['001', '002', '003', '004', '005', '006', '007', '008']

for expID in ExpId:
	f = fits.open('LeoV-RT/LeoV_%3s.fits.gz' % expID)
	spec = f[2].data
	hdr = f[2].header
	lamb = spec['LAMBDA'].flatten()
	fl = spec['SPEC'].flatten()
	inv = InverseVariance(spec['IVAR'].flatten())
	header.append(hdr)
	lam.append(lamb)
	flux.append(fl)
	ivar.append(inv)

#Load in all the magnitudes for all of the stars
vmag = np.loadtxt('LeoV-RT/LeoV_magnitudes.txt', skiprows=1, usecols=1)

#------------------------------------------------------------------------------

#Sample Spectra

#Load in the sample spectra and normalise it
with fits.open("lte05000-3.00-1.5.fits") as lte:
	data = lte[0].data 
	h = lte[0].header
	tflux = np.array(data).flatten()
	wav = np.arange(h['CRVAL1'], h['CRVAL1']+h['NAXIS1']*h['CDELT1'], 
	h['CDELT1'])
	#Vacuum correct the data
	vac = 10**4 / wav
	sky = 1. - 0.00008336624212083 - (0.02408926869968 / 
	(130.1065924522 - vac**2)) + (0.0001599740894897/(38.92568793293 - vac**2))
	tlam = wav*sky
	#Chebyshev normalise the spectra
	fit = fitting.LinearLSQFitter()
	cheby = models.Chebyshev1D(7)
	cont = fit(cheby, tlam, tflux)
	tspec = Spectrum1D(spectral_axis = tlam * u.AA, flux = tflux/cont(tlam) * u.dimensionless_unscaled)

#------------------------------------------------------------------------------

#Write all the functions necessary for the code

#Define a fitting function for the data
vr = np.arange(-300,300,1)
def tfit(spec, tspec, vr, iv):
	chisq = np.zeros(len(vr))
	for i in range(len(vr)):
		lt = tspec.spectral_axis * (1 + vr[i]/3e5)
		reso = ndi.filters.gaussian_filter(tspec.flux, sigma=0.32/0.1)
		tempspec = Spectrum1D(spectral_axis = lt , 
		flux = reso *u.dimensionless_unscaled)
		resamp = SplineInterpolatedResampler()
		trspec = resamp(tempspec, spec.spectral_axis)
		chisq[i] = np.sum((spec.flux - trspec.flux)**2*(iv.array))
	return chisq

#Define a function to help determine the metallicities 
def met(met1,met2,met3, vmag):
	eqw = 0.5*met1 + met2 + 0.5*met3
	a = -2.9
	b = 0.187
	c = 0.422
	d = -0.882
	e = 0.0133
	meta = a + b * vmag + c * eqw + d * eqw ** (-1.5) + e * vmag * eqw
	return meta

#------------------------------------------------------------------------------

#Complete all the necessary calculations on each of the stars printing the
#results and plotting the data

#Define the regions where the calcium triplets can be found
cal1 = SpectralRegion(8494 * u.AA, 8551 * u.AA)
cal2 = SpectralRegion(8538 * u.AA, 8546 * u.AA)
cal3 = SpectralRegion(8658 * u.AA, 8666 * u.AA)

#Normalise each of the spectra, find the velocity, metallicity and SNR 
#for each and plot them
for i in range(len(lam)):
	#Define the region of the spectra we are interested in
	reg = (lam[i]<8700) & (lam[i]>8000)
	iv = ivar[i][reg]
	#Add noise to the spectra
	for k in range(80):
		sig = 1/np.sqrt(iv.array)
		noise = np.random.normal(scale = sig)
		nflux = (flux[i][reg]+noise)
	#Chebyshev normalise the data
	fit = fitting.LinearLSQFitter()
	cheby = models.Chebyshev1D(7)
	cont = fit(cheby,lam[i][reg],nflux)
	normf = nflux/cont(lam[i][reg])
	spec = Spectrum1D(spectral_axis = lam[i][reg]*u.AA, 
	flux = normf*u.dimensionless_unscaled, uncertainty = iv)
	#Find the velocities of the stars
	chisq = tfit(spec,tspec, vr, iv)
	velo = chisq.argmin()
	v.append(vr[velo])
	print('The velocity of star:',i+1,'=', vr[velo])
	#Correcting the wavelengths standard deviation in python to 
	#create a new spectra
	newlam = spec.spectral_axis * (1 - vr[velo] /(3e5))
	newspec = Spectrum1D(spectral_axis = newlam, 
	flux = normf*u.dimensionless_unscaled, uncertainty = iv)
	#Find the metallicity of the stars  
	met1 = equivalent_width(newspec, regions=cal1)
	met2 = equivalent_width(newspec, regions=cal2)
	met3 = equivalent_width(newspec, regions=cal3)
	print('The metallicity of star:',i+1,"=",
	met(met1.value,met2.value,met3.value,vmag[i]))
	#Measure the signal to noise ratio
	sn = snr_derived(newspec)
	print('The signal to noise ratio of star',i+1,'=', sn)

#Find the standard deviation of the velocity
dev=np.std(v)

#------------------------------------------------------------------------------

#Find the mass-to-light ratio and its standard deviation of the Leo V Dwarf 
#Galaxy and print it
halfmass=580*70*(dev)**2
mtol=halfmass/(4.9e3)
print('The mass within the half-light radius is:', halfmass)
print('The standard deviation of the half-light mass is:', dev)
print('The mass-to-light ratio is:', mtol)

#------------------------------------------------------------------------------

