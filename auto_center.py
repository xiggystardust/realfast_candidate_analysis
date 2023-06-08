from itertools import islice, chain, tee
import code
import numpy as np

# Enter the imaging output name base here, and where you want to put the outputs
resultname='18B-405'
resultdir='/lorule/scratch/sbs0016/realfast_imaging/calvin/imaging/18B-405/' # Be sure to add the tailing /

# Enter the full directory for the measurement set here.
visname='18B-405.sb36671532.eb36805737.58655.94299563658.ms'

# Enter the source name here (get source name from Realfast plots or from listobs output).
sourcename='R3'

# Enter the phase calibrator name here (get phasecal name from listobs output).
phasename='J0217+7349'

# Enter the phase center of the image (get peak RA and Dec from Realfast plot), enter in hm dm. For no centering, make ''
phase_center = 'J2000 hm dm'

print("Starting to process...")
print(visname, sourcename)

print("Will put results named in the directory and naming scheme of...")
print(resultdir,resultname)

print("\nRunning listobs.")
listobs(vis=visname,listfile=resultdir+visname+'.txt')

# Reset initial flags
print("\nResetting initial flags.")
flagmanager(vis=visname,mode='restore',versionname='main')


# Here is a list of all numbers 2^i*3^j*5^k up to 15000.
# We'll use this later.
ham_list = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 9.0, 10.0, 12.0, 15.0, 16.0, 18.0, 20.0, 24.0, 25.0, 27.0, 30.0, 32.0, 36.0, 40.0, 45.0, 48.0, 50.0, 54.0, 60.0, 64.0, 72.0, 75.0, 80.0, 81.0, 90.0, 96.0, 100.0, 108.0, 120.0, 125.0, 128.0, 135.0, 144.0, 150.0, 160.0, 162.0, 180.0, 192.0, 200.0, 216.0, 225.0, 240.0, 243.0, 250.0, 256.0, 270.0, 288.0, 300.0, 320.0, 324.0, 360.0, 375.0, 384.0, 400.0, 405.0, 432.0, 450.0, 480.0, 486.0, 500.0, 512.0, 540.0, 576.0, 600.0, 625.0, 640.0, 648.0, 675.0, 720.0, 729.0, 750.0, 768.0, 800.0, 810.0, 864.0, 900.0, 960.0, 972.0, 1000.0, 1024.0, 1080.0, 1125.0, 1152.0, 1200.0, 1215.0, 1250.0, 1280.0, 1296.0, 1350.0, 1440.0, 1458.0, 1500.0, 1536.0, 1600.0, 1620.0, 1728.0, 1800.0, 1875.0, 1920.0, 1944.0, 2000.0, 2025.0, 2048.0, 2160.0, 2187.0, 2250.0, 2304.0, 2400.0, 2430.0, 2500.0, 2560.0, 2592.0, 2700.0, 2880.0, 2916.0, 3000.0, 3072.0, 3125.0, 3200.0, 3240.0, 3375.0, 3456.0, 3600.0, 3645.0, 3750.0, 3840.0, 3888.0, 4000.0, 4050.0, 4096.0, 4320.0, 4374.0, 4500.0, 4608.0, 4800.0, 4860.0, 5000.0, 5120.0, 5184.0, 5400.0, 5625.0, 5760.0, 5832.0, 6000.0, 6075.0, 6144.0, 6250.0, 6400.0, 6480.0, 6561.0, 6750.0, 6912.0, 7200.0, 7290.0, 7500.0, 7680.0, 7776.0, 8000.0, 8100.0, 8192.0, 8640.0, 8748.0, 9000.0, 9216.0, 9375.0, 9600.0, 9720.0, 10000.0, 10125.0, 10240.0, 10368.0, 10800.0, 10935.0, 11250.0, 11520.0, 11664.0, 12000.0, 12150.0, 12288.0, 12500.0, 12800.0, 12960.0, 13122.0, 13500.0, 13824.0, 14400.0, 14580.0, 15000.0]
max_image_size = 8192.0 # MUST BE A HAMMING NUMBER FROM THE LIST ABOVE.

# Carry out flagging on full dataset
print("\nRunning flag command. This may take a little while.")
flagdata(vis=visname,mode='rflag',timedevscale=3.0,freqdevscale=3.0,field=sourcename,flagbackup=True)
flagdata(vis=visname,mode='rflag',timedevscale=3.0,freqdevscale=3.0,field=phasename,flagbackup=True)


# Calculate the cell size (stripped from our previous cellsize.py script).
# This ensures that our image has three pixels across the resolution element of the array.
print("\nCalculating imaging parameters.")
ms.open(visname)
uv_range = ms.range(['uvdist'])
uv_max = uv_range['uvdist'][1]
freq_range = ms.range(['ref_frequency'])
cent_freq = np.median(freq_range['ref_frequency'])
ms.close()
c = 2.997925e8
wave = c/cent_freq
cellsize=206265.*wave/uv_max/3.
formatted_cell = [str(cellsize) + 'arcsec']


# Below here, we run an automatic calculation of the required image
# size.  The goal is to image out to twice the primary beam FWHM
# (approximately the first null in PSF sensitivity).

# Primary beam FWHM = lambda/dish size  (in radians).
primarybeam=206265.*wave/25. #in arcsec

# Make number of pixels amount to a field 2 * the primary beam FWHM.
image_size_in_pixels_initial = 2.*primarybeam/cellsize

# The number of pixels your image is across is required to be an
# integer, and also be a number that has prime factors of only 2, 3,
# and/or 5 (it has to be a "hamming number"). This makes the fourier
# transform code run a LOT faster.  The lines below ensures that we
# make the # of pixels obey this rule.
if (image_size_in_pixels_initial >= max_image_size):
    image_size_in_pixels_final = int(max_image_size)
    print("Image size forced from "+str(image_size_in_pixels_initial)+" to max size "+str(image_size_in_pixels_final))
else:
    for ham in ham_list:
        if (image_size_in_pixels_initial <= ham):
            image_size_in_pixels_final = int(ham)
            print("Image size forced from "+str(image_size_in_pixels_initial)+" to hamming number "+str(image_size_in_pixels_final))
            break

print("Proceeding with tclean with image size "+str(image_size_in_pixels_final))
print("Cell size will be "+str(formatted_cell))


# Here we run tclean to deconvolve the telescope's
# point-spread-function (PSF, or sensitivity pattern) from the actual
# sky image. There are a few different runs this does, based on
# different weighting schemes for how we treat measurement noise in the
# data. Looking at images creating from different weighting helps us
# get a good sense of where image artefacts (if any) might be coming from.
#
# We also create a small image of the phase calibrator with both
# weighting schemes to make sure it looks good.
#
# Note that the inputs into tclean have "niter" set to 10000. These
# runs will run automatic, rather than interactive, imaging. We might
# need niter to be larger if the images end up being very complex (or
# the residual images look like they still have a lot of structure).
#
print("\nReady for imaging. You'll know it's working if after a while you start seeing a series of 0%...10%.... type messages.\nThese may repeat several times.\n\nBombs away!")
for weight_scheme in ['natural','briggs']:

    outname=resultdir+resultname+'_'+sourcename+'_'+weight_scheme
    phaseoutname=resultdir+resultname+'_phase_'+phasename
    if weight_scheme == 'natural':
        print("Producing natural-weighted target image\n")
        tclean(vis=visname,imsize=[image_size_in_pixels_final],niter=10000,field=sourcename,phasecenter=phase_center,cell=formatted_cell,weighting=weight_scheme,pblimit=-0.001,imagename=outname)
        print("Producing natural-weighted phasecal image\n")
        tclean(vis=visname,imsize=[256],niter=10000,field=phasename,phasecenter=phase_center,cell=formatted_cell,weighting=weight_scheme,pblimit=-0.001,imagename=phaseoutname)
    else:
        print("Producing briggs-weighted target image\n")
        tclean(vis=visname,imsize=[image_size_in_pixels_final],niter=10000,field=sourcename,phasecenter=phase_center,cell=formatted_cell,weighting=weight_scheme,pblimit=-0.001,robust=0.0,imagename=outname)
        print("Producing briggs-weighted phasecal image\n")
        tclean(vis=visname,imsize=[256],niter=10000,field=phasename,phasecenter=phase_center,cell=formatted_cell,weighting=weight_scheme,robust=0.0,pblimit=-0.001,imagename=phaseoutname)

    # Here is the CASA command to convert the .image file into a FITS
    # file. We also export the target source's residual plot to see
    # "what's left". This can help us judge whether to go back and
    # clean more thoroughly.
    exportfits(imagename=outname+'.image',fitsimage=outname+'.fits',overwrite=True)
    exportfits(imagename=outname+'.residual',fitsimage=outname+'_residual.fits',overwrite=True)
    exportfits(imagename=phaseoutname+'.image',fitsimage=phaseoutname+'.fits',overwrite=True)

    #imstat(imagename=


print('Primary beam [arcsec]: '+str(primarybeam))
print('Initial image size in pixels: '+str(image_size_in_pixels_initial), 'Final image size in pixels: '+str(image_size_in_pixels_final))

