"""
Welcome to the Duct Taped Pipeline!
To get things going, create a .csv file named info.csv with the following format:

observation,targetName,phaseName,phaseCenter,arrConfig
XXX-XXX.sbXXXXXXXX.ebXXXXXXXX.XXXXX.XXXXXXXXXX,XXXXX,XXXXX+XXXX,J2000 XXhXXmXX.XXXX XXdXXmXX.XXXX,X
XXXXXX.sbXXXXXXXX.ebXXXXXXXX.XXXXX.XXXXXXXXXX,XXXXX,XXXXX-XXXX,J2000 XXhXXmXX.XXXX XXdXXmXX.XXXX,X
XXXXXX.X.sbXXXXXXXX.ebXXXXXXXX.XXXXX.XXXXXXXXXX,XXXXX,XXXXX+XXXX,J2000 XXhXXmXX.XXXX XXdXXmXX.XXXX,X
...

The first line is a header, with each line following being a comma separated list of the observation ID, target, phase calibrator (for that target!), the RA and Dec of the detection, and the array configuration.
For the RA and Dec, be sure to include the J2000 hm dm in your syntax.
If you do not have the RA and Dec of the detection, or you don't know the array configuration, that is okay. These sections can be omitted, but be sure to still include the commas.
Keep in mind that the array configuration IS case sensitive. Example: 'A', 'CnB'

When you have your .csv file, place this script and the .csv file in the same directory that contains all your observation subdirectories.
Be sure that your observation subdirectories and .ms files match the name of the observation, which they should if using NRAO's provided wget commands during download.
If you want, you can specify a different output directory below. It is surrounding by two lines of #
Otherwise, the default convention will be used and everything will be output in an 'output' directory found where the script is placed.

Happy imaging!
"""

from itertools import islice, chain, tee
import numpy as np
import casatools.msmetadata
import casaviewer
import os
import math
import csv
import subprocess
import logging


# First thing's first, we grab all the directories in the cwd, read info.csv into a dictionary, and open an output file where we will put in all of our parameters for future review.

print("\nOpening and reading info.csv ...")
dirs = next(os.walk('.'))[1]
if 'output' not in dirs:
	subprocess.run(['mkdir','output'])
file = open('info.csv','r',newline='')
dictReader = csv.DictReader(file)
myDict = {}
for item in dictReader:
	myDict[item['observation']] = item
print("Creating or accessing output file ...")
files = next(os.walk('output/'))[2]
outFieldNames = ['rmsNatural','rmsBriggs','noiseMicroJy','briggsNoiseMicroJy','rmsMicroToNoiseRatioNatural','rmsMicroToNoiseRatioBriggs','SEFD','natConfLevel','briggsConfLevel','antennaNum','arrConfig','repFreqGHz','totBandGHz','timeOnSource']
if 'output.csv' not in files:
	outputFile = open('output/output.csv','w',newline='')
	dictWriter = csv.DictWriter(outputFile,fieldnames=outFieldNames)
	dictWriter.writeheader()
	outputFile.close()



# Now we begin the pipeline proper. The first part checks to see if the key read in from info.csv matches a directory. If so, it enters the directory. If not, it sends a message and continues.

print("Input and output files accessed. Beginning pipeline ...\n")
for key, value in myDict.items():
	if key in dirs:
		os.chdir(key)
		print("\n\nDirectory match! Proceeding with pipeline on observation",key,"...")
	else:
		print("\n\nHmmm... It seems",key,"wasn't a directory. Moving to the next observation.")
		continue
	with open('../output/output.csv','a',newline='') as outputFile:
		try:
			dictWriter = csv.DictWriter(outputFile,fieldnames=outFieldNames)

# Here we assign our most important variables. This is also where the output directory can be changed if desired.
	
			resultName=key
#################################################################################################################################################
			resultDir='../output/' # This is where to change the output! Be sure to add the tailing /
#################################################################################################################################################
			msFile=key+'.ms'
			targetName=value['targetName']
			phaseName=value['phaseName']
			phaseCenter=value['phaseCenter']
			arrConfig=value['arrConfig']
	
# Now we run listobs and reset the flags

			print("\n\nRunning listobs. Resetting flag manager.")
			listobs(vis=msFile,listfile=resultDir+msFile+'/'+msFile+'.txt', overwrite=True)
			flagmanager(vis=msFile,mode='restore',versionname='main')
		
			ham_list = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 9.0, 10.0, 12.0, 15.0, 16.0, 18.0, 20.0, 24.0, 25.0, 27.0, 30.0, 32.0, 36.0, 40.0, 45.0, 48.0, 50.0, 54.0, 60.0, 64.0, 72.0, 75.0, 80.0, 81.0, 90.0, 96.0, 100.0, 108.0, 120.0, 125.0, 128.0, 135.0, 144.0, 150.0, 160.0, 162.0, 180.0, 192.0, 200.0, 216.0, 225.0, 240.0, 243.0, 250.0, 256.0, 270.0, 288.0, 300.0, 320.0, 324.0, 360.0, 375.0, 384.0, 400.0, 405.0, 432.0, 450.0, 480.0, 486.0, 500.0, 512.0, 540.0, 576.0, 600.0, 625.0, 640.0, 648.0, 675.0, 720.0, 729.0, 750.0, 768.0, 800.0, 810.0, 864.0, 900.0, 960.0, 972.0, 1000.0, 1024.0, 1080.0, 1125.0, 1152.0, 1200.0, 1215.0, 1250.0, 1280.0, 1296.0, 1350.0, 1440.0, 1458.0, 1500.0, 1536.0, 1600.0, 1620.0, 1728.0, 1800.0, 1875.0, 1920.0, 1944.0, 2000.0, 2025.0, 2048.0, 2160.0, 2187.0, 2250.0, 2304.0, 2400.0, 2430.0, 2500.0, 2560.0, 2592.0, 2700.0, 2880.0, 2916.0, 3000.0, 3072.0, 3125.0, 3200.0, 3240.0, 3375.0, 3456.0, 3600.0, 3645.0, 3750.0, 3840.0, 3888.0, 4000.0, 4050.0, 4096.0, 4320.0, 4374.0, 4500.0, 4608.0, 4800.0, 4860.0, 5000.0, 5120.0, 5184.0, 5400.0, 5625.0, 5760.0, 5832.0, 6000.0, 6075.0, 6144.0, 6250.0, 6400.0, 6480.0, 6561.0, 6750.0, 6912.0, 7200.0, 7290.0, 7500.0, 7680.0, 7776.0, 8000.0, 8100.0, 8192.0, 8640.0, 8748.0, 9000.0, 9216.0, 9375.0, 9600.0, 9720.0, 10000.0, 10125.0, 10240.0, 10368.0, 10800.0, 10935.0, 11250.0, 11520.0, 11664.0, 12000.0, 12150.0, 12288.0, 12500.0, 12800.0, 12960.0, 13122.0, 13500.0, 13824.0, 14400.0, 14580.0, 15000.0]
			max_image_size = 15000.0 # MUST BE A HAMMING NUMBER FROM THE LIST ABOVE.

# Next we flag the data for both the target and phase calibrator. timedev and freqdev are assinged for when rflag determines something as an outlier to be flagged.

			print("Beginning flagging. This may take a while.")
			flagdata(vis=msFile,mode='rflag',timedevscale=3.0,freqdevscale=3.0,field=targetName,flagbackup=True)
			flagdata(vis=msFile,mode='rflag',timedevscale=3.0,freqdevscale=3.0,field=phaseName,flagbackup=True)
			print("Flagging successful! Continuing pipeline.")


	
# Below we perform some calculations to get values we will need later. These include the size of each pixel, the primary beam size, the representative frequency, and the total image sizes.
# Note that after the total size is found, it is then rounded up to one of the hamming numbers found in the list above. This is to make the FFTs faster and reduce time.
	
			ms.open(msFile)
			msmd.open(msFile)
			
		
			print("\n\nCalculating representative frequency, cellsize, primary beam, and image sizes ...")
			specWins = msmd.spwsforfield(targetName)
			uv_range = ms.range(['uvdist'])
			uv_max = uv_range['uvdist'][1]
			minFreq = msmd.chanfreqs(spw=specWins[0],unit='Hz')[0]
			maxFreq = msmd.chanfreqs(spw=specWins[-1],unit='Hz')[-1]
			maxWidth = msmd.chanwidths(spw=specWins[-1],unit='Hz')[-1]
			maxFreq += maxWidth
			repFreq = (minFreq + maxFreq)/2
			c = 2.997925e8
			wave = c/repFreq
			cellsize=206265.*wave/uv_max/3
			formatted_cell = [str(cellsize) + 'arcsec']
			primarybeam=206265.*wave/25
		
			image_size_in_pixels_initial = 2.*primarybeam/cellsize
			if (image_size_in_pixels_initial >= max_image_size):
				image_size_in_pixels_final = int(max_image_size)
			else:
				for ham in ham_list:
					if (image_size_in_pixels_initial <= ham):
						image_size_in_pixels_final = int(ham)
						break
		
			imWidth = 15
			midpoint = int(image_size_in_pixels_final/2)
			bottomLeft = [midpoint-imWidth,midpoint-imWidth]
			topRight = [midpoint+imWidth,midpoint+imWidth]
			rmsNatural = 0
			rmsBriggs = 0
		
			ms.close()
			msmd.close()



# With our values calculated, we move on to imaging. We create natural and briggs images of the target and phase calibrator, find RMS values, and export raster, contour map, and fits images.
# Note that the target images are phasecentered for the detection and that the contour images are zoomed in to ~20arcsecs from the detection.
# Created rasters and contour maps are in a png format.
# The RMS values are calulated by taking the RMS via imstat of the residual image created in tclean.

			print("Calculations successful. Running cleaning and image creation tasks (tclean and imview). This may take some time.")
			for weight_scheme in ['natural','briggs']:
				outname=resultDir+resultName+'/'+resultName+'_'+targetName+'_'+weight_scheme
				phaseoutname=resultDir+resultName+'/'+resultName+'_phase_'+phaseName
		
				if weight_scheme == 'natural':
					print("\n\nBeginning tclean on natural-weighted target image ...")
					tclean(vis=msFile,imsize=[image_size_in_pixels_final],niter=10000,field=targetName,phasecenter=phaseCenter,cell=formatted_cell,weighting=weight_scheme,pbcor=True,pblimit=-0.001,imagename=outname)
					print("Beginning tclean on natural-weighted phase calibrator ...")
					tclean(vis=msFile,imsize=[256],niter=10000,field=phaseName,cell=formatted_cell,weighting=weight_scheme,pblimit=-0.001,imagename=phaseoutname)
		
					print("Extracting rms value and creating raster and contour map images ...")
					rmsNatural = imstat(imagename=outname+'.residual')['rms'][0].item()
					casaviewer.imview(raster={'file':outname+'.image','scaling':0,'range':[0,rmsNatural*10],'colormap':'Greyscale 2','colorwedge':True},out={'file':outname+'.png','format':'png'})
					casaviewer.imview(contour={'file':outname+'.image','levels':[-3,3,6,9,12,15],'base':0,'unit':rmsNatural},zoom={'blc':bottomLeft,'trc':topRight,'coord':'pixel'},out=outname+'_contour.png')
					subprocess.run(['convert','-negate',outname+'.png',outname+'.png'])
					subprocess.run(['convert','-negate',outname+'_contour.png',outname+'_contour.png'])
					print("Natural weighting complete.")
		
				else:
					print("\n\nBeginning tlcean on briggs-weighted target image ...")
					tclean(vis=msFile,imsize=[image_size_in_pixels_final],niter=10000,field=targetName,phasecenter=phaseCenter,cell=formatted_cell,weighting=weight_scheme,pbcor=True,pblimit=-0.001,robust=0.0,imagename=outname)
					print("Beginning tclean on briggs-weighted phase calibrator ...")
					tclean(vis=msFile,imsize=[256],niter=10000,field=phaseName,cell=formatted_cell,weighting=weight_scheme,robust=0.0,pblimit=-0.001,imagename=phaseoutname)
		
					print("Extracting rms value and creating raster and contour map images ...")
					rmsBriggs = imstat(imagename=outname+'.residual')['rms'][0].item()
					casaviewer.imview(raster={'file':outname+'.image','scaling':0,'range':[0,rmsNatural*10],'colormap':'Greyscale 2','colorwedge':True},out={'file':outname+'.png','format':'png'})
					casaviewer.imview(contour={'file':outname+'.image','levels':[-3,3,6,9,12,15],'base':0,'unit':rmsNatural},zoom={'blc':bottomLeft,'trc':topRight,'coord':'pixel'},out=outname+'_contour.png')
					subprocess.run(['convert','-negate',outname+'.png',outname+'.png'])
					subprocess.run(['convert','-negate',outname+'_contour.png',outname+'_contour.png'])
					print("Briggs weighting complete.")
		
				print("Exporting images to fits files ...")
				exportfits(imagename=outname+'.image',fitsimage=outname+'.fits',overwrite=True)
				exportfits(imagename=outname+'.residual',fitsimage=outname+'_residual.fits',overwrite=True)
				exportfits(imagename=phaseoutname+'.image',fitsimage=phaseoutname+'.fits',overwrite=True)



# Below we calculate the theoretical noise expected from the observation.
# For array configurations C and D, this also includes confusion level with the thermal noise, which are added in quadrature. For array configurations A and B, the confusion component is negligible.
# Note that this calculation is a rough estimate, designed to see if there are extreme problems with an observation. If you wish to check this calculation, you can go to https://obs.vla.nrao.edu/ect/

			msmd.open(msFile)
		
			print("\n\nImaging process complete. Beginning theoretical noise calculations.")
	
# The observation band is determined by checking the bounds of each band with a series of if statements. The total bandwidth is found by summing up all the spectral windows.
# If the total is then multiplied by a factor determined by the band. This factor removes the unusable section of the band due to RFI.
# If the new total is above a certain threshold still, it is cropped down to the maximum amount for the band.

			print("\nDeterming observation band and total bandwidth ...")
			repFreqMHz = repFreq/1000000
			totBand = sum(msmd.bandwidths(specWins))
			if repFreqMHz < 2000:
				totBand = totBand*.6
				if totBand > 600:
					totBand = 600
				band = 'L'
			elif repFreqMHz < 4000:
				totBand = totBand*.75
				if totBand > 1500:
					totBand = 1500
				band = 'S'
			elif repFreqMHz < 12000:
				totBand = totBand*.85
				if totBand > 3400:
					totBand = 3400
				if repFreq < 8000:
					band = 'C'
				else:
					band = 'X'
			else:
				totBand = totBand*.88
				if totBand > 5280:
					totBand = 5280
				band = 'Ku'
			totBandGHz = totBand/1000
			totBandHz = totBand * 1000000

# We assume if the bandwidth is less than 2GHz, then the 8-bit sampler was used. If the bandwidth is greater than 2GHz, the 3-bit sampler is required.
# Based on the sampler, we assign a correlator efficiency.

			print("Assigning correlator efficiency ...")
			if totBandGHz > 2:
				correlatorEff = .80
			else:
				correlatorEff = .93

# Total observation time is found by subtraction the start time from end time in seconds.

			print("Calculating total observation time ...")
			totTime = 0
			sourceScans = msmd.scansforfield(targetName)
			for scn in sourceScans:
				temp = msmd.timesforscan(scn)
				totTime += temp[-1] - temp[0]

# The SEFD (System Equivalent Flux Density) is determined from a list of posted values by the VLA at https://science.nrao.edu/facilities/vla/docs/manuals/oss/performance/sensitivity

			print("Assigning SEFD ...")
			SEFDs = {'L':420,'S':370,'C':310,'X':250,'Ku':320}
			SEFD = SEFDs[band]

# If the array configuration is unknown, the code below runs an algorithm based on the antennas used during the observation to determine the configuration.
# The algorithm is:
	# Search through each station id, identify which arm the antenna is in. If the antenna is further out then the current max for that arm, make that antenna the new max.
	# If the furthest north antenna is the same distance as the furthest east or west, then it is not a hybrid. Assign configuration according to north.
	# Else if the furthest west is equal or longer than the furthest east, and if the furthest north is longest, then this is a hybrid. Assign configuration accordingly to west.
	# Else if the furthest north is longer than the furthest east, then this is a hybrid. Assign configuration according to east.

			if arrConfig == '':
				print("Array configuration not found. Beginning array configuration algorithm ...")
				stations = msmd.antennastations()
				maxN = 0
				maxW = 0
				maxE = 0
				hybrid = False
				for station in stations:
					try:
						branch = station[0]
						id = int(station[1:])
						if branch == 'N' and id > maxN:
							maxN = id
						elif branch == 'W' and id > maxW:
							maxW = id
						elif branch == 'E' and id > maxE:
							maxE = id
					except:
						continue
			
				antIds = {9:'D',18:'C',36:'B',72:'A'}
				if maxN == maxW or maxN == maxE:
					try:
						arrConfig = antIds[maxN]
					except:
						# raise some kind of exception maybe?
						temp = [9,18,36,72]
						for index, value in enumerate(temp):
							if maxN < value:
								maxN = value
								break
						arrConfig = antIds[maxN]
				elif maxW >= maxE:
					if maxN > maxW:
						hybrid = True
					try:
						arrConfig = antIds[maxW]
					except:
						temp = [9,18,36,72]
						for index, value in enumerate(temp):
							if maxW < value:
								maxW = value
								break
						arrConfig = antIds[maxW]
				else:
					if maxN > maxE:
						hybrid = True
					try:
						arrConfig = antIds[maxE]
					except:
						temp = [9,18,36,72]
						for index, value in enumerate(temp):
							if maxE < value:
								maxE = value
								break   
						arrConfig = antIds[maxE]        
				
				hybConfig = {'D':'DnC','C':'CnB','B':'BnA'}
		
			else:
				print("Array configuration found. Proceeding ...")
				for a,b in hybConfig.items():
					if arrConfig == b:
						arrConfig = a
						hybrid = True
						break

# The confusion level is caclulated based on a list from the same source as the SEFDs. The confusion is assigned based on the array configuration and band.

			print("Calculating confusion level ...")
			confusion = False
			confusionLevels = {'C':{'L':[14.8,7.4],'S':[2.4,1.2],'C':[0.4,0.2],'X':[0,0],'Ku':[0,0]},'D':{'L':[148,74],'S':[24,12],'C':[4,2],'X':[0,0],'Ku':[0,0]}}
			if arrConfig == 'C' or arrConfig == 'D':
				confusion = True
				naturalConf = confusionLevels[arrConfig][band][0]
				briggsConf = confusionLevels[arrConfig][band][1]

# We assume that all calculations are with two polarizations, add the confusion and rms calculations in quadrature, and run the values through a final noise calculation.
# After this calculation, a ratio of observed and expected noise is generated.

			print("Performing final noise calculations ...")
			antNum = msmd.antennaids()[-1] + 1
			polNum = 2
			noise = SEFD / ( correlatorEff * math.sqrt( polNum * antNum * ( antNum - 1 ) * totTime * totBandHz ) )
			noiseMicro = noise * 1000000
			briggsNoiseMicro = noiseMicro * 1.2
			if confusion:
				noiseMicro = math.sqrt(noiseMicro**2 + naturalConf**2)
				briggsNoiseMicro = math.sqrt(briggsNoiseMicro**2 + briggsConf**2)
			noiseRatioNatural = (rmsNatural * 1000000) / noiseMicro
			noiseRatioBriggs = (rmsBriggs * 1000000) / briggsNoiseMicro

# Lastly, we write out all of our parameters and calculated values to an output file which can be read with csv.DictReader.

			print("\nCalculations complete. Writing output parameters and data to output.csv.")
			if not hybrid:
				dictWriter.writerow({"rmsNatural":float(rmsNatural),"rmsBriggs":float(rmsBriggs),"noiseMicroJy":float(noiseMicro),"briggsNoiseMicroJy":float(briggsNoiseMicro),"rmsMicroToNoiseRatioNatural":float(noiseRatioNatural),"rmsMicroToNoiseRatioBriggs":float(noiseRatioBriggs),"SEFD":int(SEFD),"natConfLevel":float(naturalConf),"briggsConfLevel":float(briggsConf),"antennaNum":int(antNum),"arrConfig":arrConfig,"repFreqGHz":float(repFreq/1e9),"totBandGHz":float(totBandGHz),"timeOnSource":float(totTime)})
			else:
				dictWriter.writerow({"rmsNatural":float(rmsNatural),"rmsBriggs":float(rmsBriggs),"noiseMicroJy":float(noiseMicro),"briggsNoiseMicroJy":float(briggsNoiseMicro),"rmsMicroToNoiseRatioNatural":float(noiseRatioNatural),"rmsMicroToNoiseRatioBriggs":float(noiseRatioBriggs),"SEFD":int(SEFD),"natConfLevel":float(naturalConf),"briggsConfLevel":float(briggsConf),"antennaNum":int(antNum),"arrConfig":hybConfig[arrConfig],"repFreqGHz":float(repFreq/1e9),"totBandGHz":float(totBandGHz),"timeOnSource":float(totTime)})
		
			print("\n\nObservation",resultName,"complete. Continuing pipeline.")
			msmd.close()
		except Exception as Argument:
			print("\n\noh NO! An error has occurred. Please check info.csv and note where during pipeline the error occurred. Error sent to error.txt")
			with open('../error.txt','a') as er:
				er.write(str(Argument))
		os.chdir("../")

file.close()
print("\n\n\nPipeline complete. Output files can be found in the assigned output directory.")
