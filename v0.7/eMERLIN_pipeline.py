# The e-MERLIN data reduction pipeline
#    Copyright (C) 2013  Megan Argo
#    Bug fixes and modifications also by: Neal Jackson, Rob Beswick, Danielle Fenech, Pierre Emmanuel Belles, Anthony Rushton
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.

####################################################################################
# Data loading, calibration and imaging script for e-MERLIN data processing.
#
# The script uses inputs specified in a file given on the command line.  It loads
# all the data it finds from the specified directories, averages it in time and
# frequency (if requested - but highly recommended!), runs the flagger SERPent on
# each file in the AIPS catalogue, then runs the pipeline carrying out calibration
# and imaging.  All stages are optional and can be run separately.
#
# It is highly unlikely that the resulting images will be of publication quality.
# You will get better images by checking the flags, deleting any residual bad data,
# and re-imaging the data.  It may (if there are significant bad data remaining) be
# wise to also redo the calibration.
#
# You can run the script in segments, loading, averaging and flagging first, then
# checking the flags yourself and adding to them if needed before running the rest
# of the calibration and imaging routine.  See the doall.readme file.
#
# If the script fails, the first thing to do is run 'source /aips/LOGIN.CSH'.
####################################################################################

# Changelog
# 20140117 - Get around source names with spaces (causes fittp issues) by replacing with a dash.
# 20140116 - Fixed the flagging of IF2 after the LO change in December.  Added a forced delete
#	of CL#1 at end of the DBCON stage.
# 20131022 - Added SEFD calculation (from RJB) for QA purposes (v0.7).
# 20130902 - Added the option of applying a flag mask of known RFI.  Added some bits from NJJ.
# 20130822 - Fixed calculation of appropriate cellsize for any array.  Added flux calibration.
# 20130104 - Tidied up and fixed calibration section (v0.6).
# 20121217 - Fixed specindx issue.
# 20121123 - Bug-fixing merged script.
# 20121119 - Merged orginal script (v0.3) with version by RJB -> v0.4.
# 20120618 - Changed the averaging function to use SPLAT in one step instead of AVSPC 
#	and UVAVG seperately.  Added MULTI following SPLIT to avoid DBCON issues after flagging.
# 20120614 - Added the flagger and defined several program entry points depending on the
#	state of your data (v0.2).
# 20120611 - merged together the load and dbcon scripts from the original seperate scripts
# 	and parts of the pipeline script to do with reading the inputs file (v0.1).


# Import useful things

import os, re, time, datetime, sys, math, fnmatch
from os.path import join, getsize
from datetime import date
from collections import deque
import Utilities
#from multiprocessing import Process	# No longer needed now SERPent is parallel
#from multiprocessing import Pool
from AIPS import AIPS, AIPSDisk
from AIPSTask import AIPSTask, AIPSList
from AIPSData import AIPSUVData, AIPSImage, AIPSCat
from Wizardry.AIPSData import AIPSUVData as WizAIPSUVData
from eMERLIN_tasks import *
import eMERLIN_tasks
import math, time, datetime
from numpy import *
import itertools
from time import gmtime, strftime, localtime
ti = time.time()    # To time the script


#################################################################################################################
# Define functions needed in this file
#################################################################################################################

def time2hms(seconds):
	''' Function to convert seconds to hh:mm:ss.ss format, returns a string'''
	h=int(seconds/3600)
	m=int(seconds % 3600)/60
	s=seconds-(h*3600)-(m*60)
	h=`h`
	m=`m`
	s="%4.2f" % s
	hms=h.zfill(2)+":"+m.zfill(2)+":"+s.zfill(4)
	return hms

def parse_inp(filename):
	''' Parse the list of inputs given in the specified file. (Modified from evn_funcs.py)'''
	INPUTFILE = open(filename, "r")
	control = dict()

	# a few useful regular expressions
	newline = re.compile(r'\n')
	space = re.compile(r'\s')
	char = re.compile(r'\w')
	comment = re.compile(r'#.*')

	# parse the input file assuming '=' is used to separate names from values
	for line in INPUTFILE:
		if char.match(line):
			line = comment.sub(r'', line)
			line = line.replace("'", '')
			(param, value) = line.split('=')

			param = newline.sub(r'', param)
			param = param.strip()
			param = space.sub(r'', param)

			value = newline.sub(r'', value)
			value = value.strip()
			valuelist = value.split(', ')
			control[param] = valuelist

	return control


def checkin(control):
	''' convert the control hash to a list of global variables '''

	global userno, indisk, fitsdir, fitsfil, toload, plotdir, fittpdir, msgkill
	global targets, phsrefs, fluxcals, pointcals, pointflux, bpasscals
	global refant, refantnames, band, imsize, doback
	global solint, snver, fring_snr, setjy_fluxes
	global nchan, bchan, echan, nif, stokes, plotstokes, polarizations, dotherest
	global doload, doavg, doflag, doflagmask, doconcat, docalib, dotidy, doback	# process management
	global dosort, doselfcal, tint, dosortfirst, dodiagnostic1, dodiagnostic2, dodiagnostic3, dosefd
	global sbandl, sbandu, smootha, smoothb, smoothc, dodelete, dofringfit, dosetflux, dofluxcal, dobandpass
	global aggr1, max1, aggr2, max2, rho, ncpu, kickoutsigma			# flagger options
	global 	Sfield, beta, kappa, nu, tau, SEFDbif, SEFDeif, SEFDbchan, SEFDechan	# SEFD options

	AIPS.userno = int(control.get('userno',[0])[0])
	indisk = int(control.get('indisk', [0])[0])
	dosort = int(control.get('dosort',[1])[0])
	doselfcal = int(control.get('doselfcal',[0])[0])

	# default is no averaging
	nchan = int(control.get('nchan',[-1])[0])
	tint  = int(control.get('tint',[-1])[0])

	# default is to run all procedures
	doload = int(control.get('doload',[1])[0])
	doavg  = int(control.get('doavg',[1])[0])
	doflag = int(control.get('doflag',[1])[0])
	doflagmask = int(control.get('doflagmask',[1])[0])
	doconcat = int(control.get('doconcat',[1])[0])
	docalib = int(control.get('docalib',[1])[0])
	dobandpass = int(control.get('dobandpass',[1])[0])
	dotidy = int(control.get('dotidy',[1])[0])
	dosortfirst = int(control.get('dosortfirst',[-1])[0])
	dodiagnostic1 = int(control.get('dodiagnostic1',[1])[0])
	dodiagnostic2 = int(control.get('dodiagnostic2',[1])[0])
	dofringfit = int(control.get('dofringfit',[1])[0])
	dosetflux = int(control.get('dosetflux',[1])[0])
	dofluxcal = int(control.get('dofluxcal',[1])[0])
	dodiagnostic3 = int(control.get('dodiagnostic3',[1])[0])
	doback = int(control.get('doback',[1])[0])
	sbandl = int(control.get('sbandl',[0])[0])
	sbandu = int(control.get('sbandu',[0])[0])
	smootha = int(control.get('smootha',[0])[0])
	smoothb = int(control.get('smoothb',[0])[0])
	smoothc = int(control.get('smoothc',[0])[0])
	dotherest = int(control.get('dotherest',[1])[0])
	dodelete = int(control.get('dodelete',[-1])[0])
	dosefd = int(control.get('dosefd',[-1])[0])

	# defalut options for SERPent
	aggr1 = int(control.get('aggr1',[25])[0])
	max1 = int(control.get('max1',[32])[0])
	aggr2 = int(control.get('aggr2',[25])[0])
	max2 = int(control.get('max2',[256])[0])
	rho = float(control.get('rho',[1.5])[0])
	ncpu = float(control.get('ncpu',[8])[0])
	kickoutsigma = float(control.get('kickoutsigma',[1.5])[0])

	refant = control.get('refant',['none'])[0]

	# Options relating to loading from disk
	fitsdir = control.get('fitsdir', [])
	toload = 0

	# plotdir defaults to the cwd if not specified
	plotdir = control.get('plotdir', [False])[0]
	if (not plotdir):
		plotdir = os.getcwd()
	if (not os.path.isdir(plotdir)):
		print "Error:", plotdir, "does not exist. Check your inputs."
		sys.exit()
	if (not os.access(plotdir, os.W_OK) ):
		print "Error:", plotdir, "is not writable by you. Check your inputs."

	# fittpdir defaults to the cwd if not specified
	fittpdir = control.get('fittpdir', [False])[0]
	if (not fittpdir):
		fittpdir = os.getcwd()
	if (not os.path.isdir(fittpdir)):
		print "Error:", fittpdir, "does not exist. Check your inputs."
		sys.exit()
	if (not os.access(fittpdir, os.W_OK) ):
		print "Error:", fittpdir, "is not writable by you. Check your inputs."

	solint = float(control.get('solint', [10])[0])

	bpasscals = control.get('bpasscals', [])
	for i in range(len(bpasscals)):
		bpasscals[i] = bpasscals[i].upper()
	nbp_table = 0
	if bpasscals:
		nbp_table = 1

	phsrefs = control.get('phsrefs', [])
	for i in range(len(phsrefs)):
		phsrefs[i] = phsrefs[i].upper()

	targets = control.get('targets', [])
	for i in range(len(targets)):
		targets[i] = targets[i].upper()

	pointcals = control.get('pointcals', [])
	for i in range(len(pointcals)):
		pointcals[i] = pointcals[i].upper()

	fluxcals = control.get('fluxcals', [])
	for i in range(len(fluxcals)):
		fluxcals[i] = fluxcals[i].upper()

	pointflux = control.get('pointflux', [])

	fring_snr = float(control.get('fring_snr', [9])[0])

	interactive = int(control.get('interactive', [False])[0])

	imsize = int(control.get('imsize', [1024])[0])
	if imsize:
		imsize = imsize,imsize

#	cellsize = control.get('cellsize',[0])[0]

	msgkill = int(control.get('msgkill', [-5])[0])

	assert (len(phsrefs) == len(targets)), """unmatched number of phase reference and target sources! """

	if (not AIPS.userno) or (not indisk):
		print 'Error: You must set both your AIPS userno and indisk.'
		sys.exit()
	if (not targets) or (not phsrefs) or (not pointcals) or (not bpasscals):
		print 'Error: You must specify targets, phaserefs and bpasscals.'
		sys.exit()
	if not pointflux:
		print 'Error: You must specify the fluxes for your pointcal!'
		sys.exit()

	# SEFD control parameters -- RJB
	SEFDbchan = int(control.get('SEFDbchan', [1])[0])
	SEFDechan = int(control.get('SEFDechan', [0])[0])
	Sfield= float(control.get('Sfield', [0])[0])
	beta= float(control.get('beta', [0])[0])
	tau= float(control.get('tau', [0])[0])
	nu= float(control.get('nu', [0])[0])
	kappa= float(control.get('kappa', [0])[0])
	SEFDbif = int(control.get('SEFDbif', [1])[0])
	SEFDeif = int(control.get('SEFDeif', [1])[0])


	# Get the user to check the inputs before we actually begin.
	print "***********************************************************"
	print "Key parts of your input file is read as follows"
	print "AIPS number = " + format(AIPS.userno) + " on disk " + format(indisk) + "."
	print "Targets = " + format(targets)
	print "Phase cals = " + format(phsrefs)
	print "Flux density cals = "+ format(fluxcals)
	print "Bandpass cals = "+ format(bpasscals)
	print "Point flux cals = "+ format(pointcals)
	print ""
	print "Options choosen in:", sys.argv[1]
	print "doload = "+ format(doload)
	if doload == 1 :
		print "With:"  
		print "             **With dosortfirst =" + format(dosortfirst)
	print "doavg = "+ format(doavg)
	if doavg == 1 :
		print "With:"
		print "             **Spectral averaging to a final number of sub-band channels =" + format(nchan)
		print "             **Time averaging to =" + format(tint)+"sec"
		print "             **Lowest sub-band extracted =" + format(sbandl)
		print "             **Highest sub-band extracted =" + format(sbandu)
		print "             **Spectral smoothing parameters (smooth parms) =" + format(smootha)+","+ format(smoothb)+","+ format(smoothc)
		print "             **UV sorting averaged data dosort =" + format(dosort)
		print "dodiagnostic1 = "+ format(dodiagnostic1)
	print "doflagmask = "+ format(doflagmask)
	print "doflag = "+ format(doflag)
	if doflag == 1 :
		print "With:"
		print "             Auto-flagger settings:"
		print "             **aggr1 = " +format(aggr1)
		print "             **max1 = " +format(max1)
		print "             **aggr2 = " +format(aggr2)
		print "             **max2 = " +format(max2)
		print "             **rho = " +format(rho)
		print "             **kickoutsigma = " +format(kickoutsigma)
		print "             **ncpu = " +format(ncpu) #MAX= number of baselines or CORES available
		print "             See SERPent cookbook for meanings - default values set"
	print "dodiagnostic2 = "+ format(dodiagnostic2)
	print "doback = "+ format(doback)
	print "doconcat = "+ format(doconcat)
	if doconcat == 1 :
		print "With:"
		print "             Deleting work files after combining(1=yes) =" + format(dodelete)
	print "docalib = "+ format(docalib)
	if docalib == 1 :
		print "With sub-sections:"
		print "             **Initial pre-cal diagnostic plots = " + format(dodiagnostic3)
		print "             **fring fitting of data = " + format(dofringfit)
		print "             **Setting flux density of 3C286 = " + format(dosetflux)
		print "             **Flux calibration (bootstrap from 3C286 values in input file) = " + format(dofluxcal)
		print "             **Do a bandpass calibration using " + format(bpasscals) + " dobandpass " + format(dobandpass)
		print "             **run other untested stuff..  = " + format(dotherest)
	print "dosefd = "+ format(dosefd)
	if dosefd == 1 :
		print "With SEFD calculation settings:"
		print "              intergration time (s) = " + format(tau)
		print "              Final spec chan width (Hz)= " + format(beta)
		print "              widar corrlator efficieny (esitimate)= " + format(nu)
		print "              Hanning smoothing SNR improvement= " + format(kappa)

	print "***********************************************************"
	print AIPSCat(indisk)
	print "***********************************************************"




#def setup(uvdata, dosort, cellsize):
def setup(uvdata, dosort):
	''' extract useful information and set up a few useful parameters '''

	# Create a list of all cals and a list of all sources
	fringlist=[]
	for i in range(len(phsrefs)):
		fringlist.append(phsrefs[i])
	for i in range(len(fluxcals)):
		fringlist.append(fluxcals[i])
	for i in range(len(pointcals)):
		fringlist.append(pointcals[i])
	for i in range(len(bpasscals)):
		fringlist.append(bpasscals[i])
	# Remove duplicate source entries
	tmpfringlist = set(fringlist)
	fringlist=[]
	for source in tmpfringlist:
		fringlist.append(source)
	sourcelist=[]
	for i in range(len(fringlist)):
		sourcelist.append(fringlist[i])
	for i in range(len(targets)):
		sourcelist.append(targets[i])
	caliblist=[]
	for i in range(len(phsrefs)):
		caliblist.append(phsrefs[i])
	for i in range(len(pointcals)):
		caliblist.append(pointcals[i])
	bpasslist=[]
	for i in range(len(bpasscals)):
		bpasslist.append(bpasscals[i])
	allsources = uvdata.sources

	# Set up uvdata and make sure it's in a fit state for use
	if dosort == 1 :
		print 'Tidying up your data before we start:'
		runmsort(uvdata)
		runindxr(uvdata)
	if dotidy == 1:
		print 'Deleting old PL, SN, BP and CL tables.'
		uvdata.zap_table('SN', -1)
		uvdata.zap_table('PL', -1)
		uvdata.zap_table('BP', -1)
		for table in uvdata.tables :
			if 'CL' in table[1] :
				if table[0]>1 :
					uvdata.zap_table('CL', table[0])

	stokes = uvdata.stokes

	global plotstokes
	assert(len(stokes) > 0)
	if len(stokes) == 1:
		plotstokes = stokes[0]
	else:
		plotstokes = 'HALF'

	polarizations = uvdata.polarizations
	polarizations.sort()
	assert(len(polarizations) > 0)

	chanwidth = uvdata.header.cdelt[2]  # in Hz
	nchan = uvdata.header.naxis[2]
	# ignore the outer 10% of channels for certain tasks
	bchan = int(math.floor(nchan * 0.1))
	echan = nchan - bchan
	nif = uvdata.header.naxis[3]
	ifwidth = chanwidth*nchan

	freq = uvdata.header.crval[2] / 1000000000
	band = ''
	if (freq > 1.3) and (freq < 1.7):
		band = 'L'
	if (freq > 4) and (freq < 8):
		band = 'C'
	if (freq > 22) and (freq < 24):
		band = 'X'	# some day

	return (sourcelist, fringlist, caliblist, bpasslist, nchan, nif, stokes, plotstokes, polarizations, bchan, echan, freq, band)


def dovplot(uvdata, docalib, gainuse, flagver, doband, bpver, step, mode):
	''' set up parameters and run vplot diagnostics '''
	if (mode == 'a&p'):
		# Plot amp on cals (UVPLT, VPLOT)
		aparm = 0,0
		bparm = 12,1,0
		for pol in stokes:
			runvplot(uvdata, fringlist, pol, timer, refantlist, refantlist, 1, nif, bchan, echan, -1, 0, 0, -1, 0, aparm, bparm, refantn, dotv, 9)

		# LWPLA to disk
		filename = 'VPLOT_' + str(step) + '_amp.ps'
		outfile = os.path.join(plotdir, filename)
		print 'Plotting VPLOT diagnostics:' + outfile
		runlwpla(uvdata, outfile)

	# Plot phase on cals (UVPLT, VPLOT)
	aparm = 0,0
	bparm = 12,2,0
	runvplot(uvdata, fringlist, 'I', timer, refantlist, refantlist, 1, nif, bchan, echan, -1, 0, 0, -1, 0, aparm, bparm, refantn, dotv, 9)

	# LWPLA to disk
	filename = 'VPLOT_' + str(step) + '_phase.ps'
	outfile = os.path.join(plotdir, filename)
	print 'Plotting VPLOT diagnostics:' + outfile
	runlwpla(uvdata, outfile)


def flagger(name, klass, seq, indisk):
	''' create directory, write control file and run SERPent '''
	mkdir = 'serpent_' + name + '_' + klass + '_' + str(seq)
	os.system('mkdir ' + mkdir)
	path2folder = os.path.join(os.getcwd(), mkdir)
	path2file = os.path.join(path2folder, 'SERPent_input.py')
	print "Writing ", path2file, "for file: ", name, klass, seq
	# write ctrl file
	inputfile = open(path2file, 'w')
	print >> inputfile, 'AIPS_user_number = ', AIPS.userno
	print >> inputfile, 'Name = \'' + name + '\''
	print >> inputfile, 'Klass = \'' + klass + '\''
	print >> inputfile, 'Disk = ' + str(indisk)
	print >> inputfile, 'Seq = ' + str(seq)
	print >> inputfile, 'path2folder = \'' + path2folder + '/\''
	print >> inputfile, 'which_baselines = \'all\''
	print >> inputfile, 'flagging_options = \'choose\''
	print >> inputfile, 'aggressiveness_first_run = ', aggr1
	print >> inputfile, 'max_subset_first_run = ', max1
	print >> inputfile, 'aggressiveness_second_run = ', aggr2
	print >> inputfile, 'max_subset_second_run = ', max2
	print >> inputfile, 'rho = ', rho
	print >> inputfile, 'kickout_sigma_level = ', kickoutsigma
	print >> inputfile, 'NCPU = ', int(ncpu)
	print >> inputfile, 'zero_level = \'yes\'' 
	print >> inputfile, 'phasecal = \'no\'' 
	print >> inputfile, 'write2log = \'no\''
	inputfile.close()

	# run the flagger
	os.chdir(path2folder)
	os.system('cp ../SERPent.py .') 
	os.system('ParselTongue SERPent.py')
	os.chdir('../')

#print 'MEAN SEFDLT:',numpy.mean(SEFDLT)
def nanmean(data,**args):
	return numpy.ma.filled(numpy.ma.masked_array(data,numpy.isnan(data)).mean(**args),fill_value=numpy.nan)



today = date.today()
now = today.strftime("%Y%m%d")
mylist=[]


# AIPS is chatty. Supress most messages on the terminal but log everything.
AIPSTask.msgkill = -4





#################################################################################################################
# End of function definitions - start doing useful things
#################################################################################################################


# Sanity checks before we start

print "-----------------------------------------------------------"
print "********** Data reduction pipeline for e-MERLIN ***********"
print "-----------------------------------------------------------"
print "  Your data will be loaded, averaged, flagged and DBCON'd  "
print "  or some selection of these tasks depending on the inputs "
print "  before being calibrated and imaged. CHECK YOUR OUTPUTS!!!"

# Check for an inputs file
if len(sys.argv)==2 :
	print "Using inputs specified in", sys.argv[1]
	afile = sys.argv[1]
	if os.path.isfile(afile) :
		control = parse_inp(afile)
	else :
		print "Error:" + afile + "does not exist, quitting."
		sys.exit()
else :
	print "Error: no parameter file specified, quitting."
	print "Usage: parseltongue pipeline.py inputs.txt"
	sys.exit()


# Check the inputs - must do this to set AIPS.userno *BEFORE* task definitions!
checkin(control)


inputinput('Continue? ([y]/n): ')

if (doavg != 1) or ((doavg == 1) and ((nchan == -1) and (tint == -1))):
	inputinput("You requested NO averaging. Are you ABSOLUTELY SURE? ([y]/n): ")

inputinput('If you observed in mixed-mode, you MUST separate your data before running the pipeline! Proceed? ([y]/n): ')


dotv = -1
AIPSTask.msgkill = msgkill
logfile = os.path.join(plotdir, 'eMERLIN_pipeline.log')
AIPS.log = open(logfile, 'a')
print "Logging AIPS messages to", logfile

outdisk = indisk


#################################################################################################################
# Start loading data
#################################################################################################################

if doload == 1 :

	print "------------------------------------------------------------"
	print " Loading e-MERLIN data from disk and doing some houekeeping"
	print "------------------------------------------------------------"


	# Check we can see the data directory, otherwise quit
	loadlist=[]
	print "Checking we can see your data directories.... "
	for directory in fitsdir:
		if os.path.exists(directory[:-1]):
			loadlist.append(directory,)
		else :
			print "Error: the specified directory doesn't seem to exist!"
			print "Problem:", directory
			print "Problem:", fitsdir
			print "Check your inputs."
			sys.exit()

	print "OK, proceeding with loading."
	loadlist.sort()				# Sort the list before we load (makes DBCON more efficient later)

	thisdir = 0
	nfiles = 0
	fitslist = []
	for dir in loadlist:						# for each sub-directory...
		thisdir = thisdir + 1
		numdirs = len(dir)
		fitslist = [f for f in filter(nodot, os.listdir(dir))	# make a list of files to load
		if os.path.isfile(os.path.join(dir, f))]
		for fitsfile in fitslist:				# for each file...
			if fnmatch.fnmatch(fitsfile, '*.fits'):
				descriptor = dir + "/" + fitsfile
				print "Loading " + descriptor
				# njj: Use softlink to get around long datafile name problems
				os.system('rm EMTEMP.FITS')
				os.system('ln -s '+descriptor+' EMTEMP.FITS')
				runfitld('./EMTEMP.FITS', indisk, thisdir)

				# njj: set sort order to T* (as it is written this way) - this means that
				#      indxr will not bomb out and we can sort to TB later if needed
				#      I think we need WizAIPSUVData to do the header update here
				print 'Setting sortorder to T*'
				uvdata = WizAIPSUVData('TMP','UVDATA',outdisk,thisdir)
				uvdata.header.sortord='T*'
				uvdata.header.update()

				# reset the filename based on the source name in the SU table
				uvdata = AIPSUVData('TMP','UVDATA',outdisk,thisdir)
				if len(uvdata.sources) > 1 :
					print "Error: Too many sources in the file! Run on single-source files only."
					sys.exit()

				sourcelist = uvdata.table('SU',1)
				outn = re.sub(r"\s+", '-', uvdata.sources[0])	# catch source names with spaces
				if len(outn)>12 :
					outn = outn[0:12]			# catch long source names
#				uvdata.rename(outn,'UVDATA',0)
				uvdata.rename(outn,'UVDATA',thisdir)
				inname = outn

				# uvfix each file here, only if not averaging later
				if doavg != 1:
					uvdata = AIPSUVData(inname,'UVDATA',indisk,thisdir)
					runuvfix(uvdata, inname, indisk, thisdir)

				if dosortfirst == 1 :
					# sort each file, if requested
					uvdata = AIPSUVData(inname,'UVDATA',indisk,thisdir)
					runuvsrt(uvdata, inname, indisk, thisdir)
				else :
					print "You didn't request Time-Baseline sorting on full (unaveraged) data, skipping."

				# index each file
				uvdata = AIPSUVData(inname,'UVDATA',indisk,thisdir)
				runindxr(uvdata)

				# flag the Lovell - Mk2 baseline (SERPent does it anyway, but not all users will make use of it).
				flagLOMK(uvdata)



#################################################################################################################
# Flag the data using the flag mask before averaging.
#################################################################################################################

if doflagmask == 1 :

	print "------------------------------------------------------------"
	print "                  Applying the flag mask                    "
	print "------------------------------------------------------------"

	pca = AIPSCat(indisk)
	for fitsfil in pca[indisk]:
		uvdata = AIPSUVData(fitsfil.name, fitsfil.klass, indisk, fitsfil.seq)
		flagmask(uvdata)

	print "------------------------------------------------------------"
	print "           Flag mask applied - check your data!             "
	print "------------------------------------------------------------"



#################################################################################################################
# Average the data
#################################################################################################################

if doavg == 1 :

	pca = AIPSCat(indisk)

	# average only if requested
	if tint != -1 or nchan != -1 :
		print "Averaging in time and/or frequency with SPLAT."
		aparm1 = 0
		for fitsfil in pca[indisk]:
			uvdata = AIPSUVData(fitsfil.name, fitsfil.klass, indisk, fitsfil.seq)
			if len(uvdata.sources) > 1 :
				print "Error: Too many sources in the file! Run on single-source files only."
				sys.exit()

			runsplat(uvdata, nchan, tint, sbandl, sbandu, smootha, smoothb, smoothc, indisk)

			# Remove unaveraged file, then rename the averaged file
			uvdata.zap()
			uvdata = AIPSUVData(fitsfil.name,'SPLAT',indisk,fitsfil.seq)
			uvdata.rename(fitsfil.name,'UVDATA',fitsfil.seq)

			uvdata = AIPSUVData(fitsfil.name,'UVDATA',indisk,fitsfil.seq)
			runuvfix(uvdata, fitsfil.name, indisk, fitsfil.seq)

			# njj: think the outputs are also T* sorted, relabel them
			print 'Setting sortorder to T*'
			uvdata = WizAIPSUVData(fitsfil.name,'UVDATA',indisk,fitsfil.seq)
			uvdata.header.sortord='T*'
			uvdata.header.update()

		print "Your data are now averaged as requested. Next the data should be flagged with SERPent."

	else :
		print "You requested no averaging, skipping."


#################################################################################################################
# Optional dosort (probably superfluous - esp if we can get some kind of check working)
#################################################################################################################

if dosort == 1 :
	pca = AIPSCat(indisk)

	for fitsfil in pca[indisk]:
		uvdata = AIPSUVData(fitsfil.name, fitsfil.klass, indisk, fitsfil.seq)
		print "UVSRT: SORTING DATA  " + fitsfil.name + '.' + fitsfil.klass + '.' + format(fitsfil.seq)
		runuvsrt(uvdata, fitsfil.name, indisk, fitsfil.seq)
		runindxr(uvdata)


#################################################################################################################
# Diagnostic plotting Section 1 = before auto-flagging
#################################################################################################################

if dodiagnostic1 == 1 :
	print "------------------------------------------------------------"
	print " STAGE 1 diagnostics for e-MERLIN data (before Auto-flagger)"
	print "                     doiagnostic1 = 1                       "
	print "------------------------------------------------------------"

	pca = AIPSCat(indisk)
	for fitsfil in pca[indisk]:
			uvdata = AIPSUVData(fitsfil.name, fitsfil.klass, indisk, fitsfil.seq)
			if len(uvdata.sources) > 1 :
				print "Error: Too many sources in the file! Run on single-source files only."
				sys.exit()

			print "POSSM: Making plots of DATA  " + fitsfil.name + '.' + fitsfil.klass + '.' + format(fitsfil.seq)
			source = uvdata.sources[0],
			runquickpossm(uvdata, source)
			plotfile =  fitsfil.name + '.' + '.' + str(fitsfil.seq) + '_PREAOF_POSSM_' + strftime("%y%m%d_%H%M") + '.ps'
			outfile = os.path.join(plotdir,plotfile)
			print 'Plotting POSSM diagnostics:' + outfile
			runlwpla(uvdata, outfile)

			runquickvplot(uvdata)
			plotfile =  fitsfil.name + '.' + '.' + str(fitsfil.seq) + '_PREAOF_VPLOT_' + strftime("%y%m%d_%H%M") + '.ps'
			outfile = os.path.join(plotdir,plotfile)
			print 'Plotting VPLOT diagnostics:' + outfile
			runlwpla(uvdata, outfile)


	print '********************************************************************'
	print '    Finished plotting STAGE 1 pre-flag/calibration diagnostic plots.'
	print '         It would be wise to inspect them closely.                  '
	print '********************************************************************'



#################################################################################################################
# Flag the data with SERPent after averaging (if averaging).
#################################################################################################################
	

if doflag == 1 :

	print "------------------------------------------------------------"
	print "               Autoflagging with SERPent                    "
	print "------------------------------------------------------------"

	pca = AIPSCat(indisk)

	# Run parallelised version of the flagger on all files in catalog
	for fitsfil in pca[indisk]:
		uvdata = AIPSUVData(fitsfil.name, fitsfil.klass, indisk, fitsfil.seq)
		if len(uvdata.sources) > 1 :
			print "Error: Too many sources in the file! Run on single-source files only."
			sys.exit()

		apply(flagger, (uvdata.name, uvdata.klass, uvdata.seq, indisk,),)
	print "------------------------------------------------------------"
	print "         Autoflagging complete - check your data!           "
	print "------------------------------------------------------------"

	# TASAV the FG tables and write to disk
	for fitsfil in pca[indisk]:
		uvdata = AIPSUVData(fitsfil.name, fitsfil.klass, indisk, fitsfil.seq)
		srcname = re.sub(r"\s+", '-', fitsfil.name)
		fittpfile =  srcname + '.TASAV.' + str(fitsfil.seq) + '_AOF_' + strftime("%y%m%d_%H%M") + '.fits'
		runtasav(uvdata, fittpdir, fittpfile, indisk)




#################################################################################################################
# Diagnostic plotting 2 = after auto-flagging
#################################################################################################################

if dodiagnostic2 == 1 :
	print "------------------------------------------------------------"
	print " STAGE 2 diagnostics for e-MERLIN data (AFTER Auto-flagger)"
	print "                        dodiagnostic2 = 1                   " 
	print "------------------------------------------------------------"

	pca = AIPSCat(indisk)
	for fitsfil in pca[indisk]:
			uvdata = AIPSUVData(fitsfil.name, fitsfil.klass, indisk, fitsfil.seq)
			if len(uvdata.sources) > 1 :
				print "Error: Too many sources in the file! Run on single-source files only."
				sys.exit()


			print "POSSM: Making plots of DATA  " + fitsfil.name + '.' + fitsfil.klass + '.' + format(fitsfil.seq)
			source = uvdata.sources[0],
			runquickpossm(uvdata, source)
			plotfile =  uvdata.name + '.' + str(uvdata.seq) + '_POSTAOF_POSSM_' + strftime("%y%m%d_%H%M") + '.ps'
			outfile = os.path.join(plotdir,plotfile)
			print 'Plotting POSSM diagnostics:' + outfile
			runlwpla(uvdata, outfile)

			runquickvplot(uvdata)
			plotfile =  uvdata.name + '.' + str(uvdata.seq) + '_POSTAOF_VPLOT_' + strftime("%y%m%d_%H%M") + '.ps'
			outfile = os.path.join(plotdir,plotfile)
			print 'Plotting VPLOT diagnostics:' + outfile
			runlwpla(uvdata, outfile)


	print '***********************************************************'
	print 'Finished plotting STAGE 2 pre-calibration diagnostic plots.'
	print '         It would be wise to inspect them closely.         '
	print '***********************************************************'


#################################################################################################################
# Backup the data before we munge it
#################################################################################################################

if doback == 1 :
	print "------------------------------------------------------------"
	print " BACKUP individual (AVERAGED) files including tables        "
	print "------------------------------------------------------------"

	pca = AIPSCat(indisk)
	for fitsfil in pca[indisk]:
		uvdata = AIPSUVData(fitsfil.name, fitsfil.klass, indisk, fitsfil.seq)
		fittpfile =  fitsfil.name + '.DATA.' + str(fitsfil.seq) + '_PRE_CONCAT_' + strftime("%y%m%d_%H%M") + '.fits'
		runfittp(uvdata, fittpdir, fittpfile)


#################################################################################################################
# DBCON the data to create one dataset
#################################################################################################################

if doconcat == 1 :

	pca = AIPSCat(indisk)
	print "Combining your datasets to create one uv database for calibration."
	for fitsfil in pca[indisk]:
		uvdata = AIPSUVData(fitsfil.name, fitsfil.klass, indisk, fitsfil.seq)
		if len(uvdata.sources) > 1 :
			print "Error: Too many sources in the file! Run on single-source files only."
			sys.exit()


	# Backup tables in case the user did any extra flagging
	for fitsfil in pca[indisk]:
		fittpfile =  fitsfil.name + '.TASAV.' + str(fitsfil.seq) + strftime("%y%m%d_%H%M") + '.fits'
		uvdata = AIPSUVData(fitsfil.name, fitsfil.klass, indisk, fitsfil.seq)
		runtasav(uvdata, fittpdir, fittpfile, indisk)

	print "First apply the HIGHEST FG table using UVCOP"
	for fitsfil in pca[indisk]:
		uvdata = AIPSUVData(fitsfil.name, fitsfil.klass, indisk, fitsfil.seq)
		uvcop = AIPSTask('UVCOP')
		uvcop.indata = uvdata
		uvcop.outclass = 'UVCOP'
		uvcop.outdisk = indisk
		uvcop.outseq = fitsfil.seq
		uvcop.flagver = get_tab(uvdata, 'FG')
		print "Apply top flag table number", get_tab(uvdata, 'FG')
		uvcop.go()
		uvdata.zap()
		# Delete any remaining 'old' FG tables so DBCON is not confused
		uvdata = AIPSUVData(fitsfil.name, 'UVCOP', indisk, fitsfil.seq)
		uvdata.zap_table('FG', -1)



	# Combine the data with DBCON
	pca = AIPSCat(indisk)
	queueNAME = deque([])
	queueKLAS = deque([])
	queueSEQ  = deque([])

	for fitsfil in pca[indisk]:
		queueNAME.append(fitsfil.name)
		queueKLAS.append(fitsfil.klass)
		queueSEQ.append(fitsfil.seq)

	print "Total files to combine: ", len(queueNAME)

	i = 0
	while len(queueNAME)>1 :
		print "Combining " + queueNAME[0] + '.' + queueKLAS[0] + '.' + format(queueSEQ[0]) + " and " + queueNAME[1] + '.' + queueKLAS[1] + '.' + format(queueSEQ[1])
		uvdata1 = AIPSUVData(queueNAME.popleft(), queueKLAS.popleft(), indisk, queueSEQ.popleft())
		uvdata2 = AIPSUVData(queueNAME.popleft(), queueKLAS.popleft(), indisk, queueSEQ.popleft())
		dbcon = AIPSTask('DBCON')
		dbcon.indata = uvdata1
		dbcon.in2data = uvdata2
		dbcon.outname = format(i)
		dbcon.outclass = 'DBCON'
		dbcon.outdisk = indisk
		dbcon.go()

		# in the interests of disk efficiency:
#		uvdata1.zap()
#		uvdata2.zap()

		indxr = AIPSTask('INDXR')
		indxr.indata = AIPSUVData(dbcon.outname, dbcon.outclass, indisk, dbcon.outseq)
		indxr.go()

		queueNAME.append(format(i))
		queueKLAS.append('DBCON')
		queueSEQ.append(1)

		i = i + 1


	# Run MSORT on the final DBCONd file
	print "Final sort and index...."
	msort = AIPSTask('MSORT')
	msort.indata = AIPSUVData(queueNAME[0], queueKLAS[0], indisk, queueSEQ[0])
	uvdata = AIPSUVData(queueNAME[0], 'MSORT', indisk, queueSEQ[0])
	msort.outdata = uvdata
	msort.go()
	####################################################################################################################
	### Delete CL#1 and remake it here!
	####################################################################################################################
	uvdata.zap_table('CL', -1)	
	indxr.indata = uvdata
	indxr.go()
	uvdata.rename('ALL','DBCON',0)

	# QUACK dBCONd data file 
	print "Then quack the data!..quack-quack-OOPS..FG#1 created.."
	quack = AIPSTask('QUACK')
	quack.indata = uvdata
	quack.opcode ='BEG'
	quack.reason = 'quack1min'
	quack.aparm[2] = 1.0
	quack.go()

	# LISTR DBCONd data file for records
	print "Create a scan list of data in file"
	listr = AIPSTask('LISTR')
	listr.indata = uvdata
	listr.optype = 'SCAN'
	textfile = 'DBCONNED_FILE_Listr_SCAN_' + strftime("%y%m%d_%H%M") + '.txt'
	listr.outprint = os.path.join(plotdir,textfile)
	listr.docrt = -1
	listr.go()



	# Tidy up: delete UVDATA/SPLIT and most of DBCON files
	print "Tidying up your AIPS catalogue."
	pca = AIPSCat(indisk)
	for fitsfil in pca[indisk]:
		if fitsfil.klass == 'DBCON' :
			if fitsfil.name.isdigit() :
				print "Zapping: " + fitsfil.name + '.' + fitsfil.klass + '.' + format(fitsfil.seq)
				uvdata1 = AIPSUVData(fitsfil.name, fitsfil.klass, indisk, fitsfil.seq)
				uvdata1.zap()

	if dodelete != -1 :

		# Tidy up: delete UVDATA/UVCOP and most of DBCON files
		print "Tidying up your AIPS catalogue."
		pca = AIPSCat(indisk)
		for fitsfil in pca[indisk]:
			if fitsfil.klass == 'UVDATA' or fitsfil.klass == 'UVCOP' :
				print "Zapping: " + fitsfil.name + '.' + fitsfil.klass + '.' + format(fitsfil.seq)
				uvdata1 = AIPSUVData(fitsfil.name, fitsfil.klass, indisk, fitsfil.seq)
				uvdata1.zap()
			if fitsfil.klass == 'DBCON' :
				if fitsfil.name.isdigit() :
					print "Zapping: " + fitsfil.name + '.' + fitsfil.klass + '.' + format(fitsfil.seq)
					uvdata1 = AIPSUVData(fitsfil.name, fitsfil.klass, indisk, fitsfil.seq)
					uvdata1.zap()
	else :
		print "Not deleteting ORIGINAL UV files."
		print "You will want to do this but obviously you want to check first, skipping.. for now..."
		print "NOTE FLAGS have been applied!!!"




#################################################################################################################
# Run the calibration and imaging pipeline
#################################################################################################################

# Need to have specified the datafile above, or only have one file in the aips disk

if docalib == 1 :

	pca = AIPSCat(indisk)
	if len(pca[indisk]) == 1 :
		uvdata = AIPSUVData(pca[indisk][0]["name"], pca[indisk][0]["klass"], indisk, pca[indisk][0]["seq"])
		print "Calibrating:", uvdata
	else :
		uvdata = AIPSUVData('ALL', 'DBCON', indisk, 1)
		if uvdata.exists() :
			pass
		else:
			print "################################################################"
			print "Error: ALL.DBCON.1 does not exist! Don't know what to calibrate."
			print "################################################################"
			sys.exit()

#	(sourcelist, fringlist, caliblist, bpasslist, nchan, nif, stokes, plotstokes, polarizations, bchan, echan, freq, band, cellsize) = setup(uvdata, dosort)
	(sourcelist, fringlist, caliblist, bpasslist, nchan, nif, stokes, plotstokes, polarizations, bchan, echan, freq, band) = setup(uvdata, dosort)

	(refant, refantn, refantlist) = set_refant(uvdata, refant)

	cellsize = findmaxb(uvdata)

	snver = 0
	clver = 0
	flagver = get_tab(uvdata, 'FG')

	if dodiagnostic3 == 1 :
		print "------------------------------------------------------------"
		print "          STAGE 3 diagnostics for e-MERLIN data"
		print "            dodiagnostic3 = 1 (part of docalib)                   " 
		print "------------------------------------------------------------"

		timer = 0,0,0,0,0,0,0,0
		aparm = 0,0,0,0,0,0,0,0,0,1
		bparm = 0,0,0,0,0,0,0,0,0,0
		for i in range(len(fringlist)) :
			sourceID = fringlist[i],
			print 'Running POSSM on', fringlist[i]
			#runpossm(uvdata, sourceID, timer, refantlist, refantlist, aparm, bparm, 0, nif, -1, 0, flagver, plotstokes, -1, 0, 'A&P', 30, 9, -1, 1)
			runquickpossm(uvdata, sourceID)

			# LWPLA to disk
			plotfile = 'POSSM_NOCAL_' + str(fringlist[i]) + '_' + strftime("%y%m%d_%H%M") + '.ps'
			outfile = os.path.join(plotdir, plotfile)
			print 'Plotting POSSM diagnostics:' + outfile
			runlwpla(uvdata, outfile)

		# Plot VPLOT diagnostics
#		for i in range(len(fringlist)) :
#			sourceID = fringlist[i],
#			print 'Running VPLOT on', fringlist[i]
#			timer = 0,0,0,0,0,0,0,0
#			runvplot(uvdata, fringlist, 'HALF', timer, refantlist, refantlist, 1, nif, bchan, echan, -1, 0, 0, -1, 0, aparm, bparm, refantn, dotv, 9)
#			# LWPLA to disk
#			outfile = os.path.join(plotdir, 'VPLOT_NOCAL_' + str(fringlist[i]) + '_' + strftime("%d%m%y_%H%M") + '.ps')
#			print 'Plotting VPLOT diagnostics:' + outfile
#			runlwpla(uvdata, outfile)

		print 'Running VPLOT on', fringlist
		timer = 0,0,0,0,0,0,0,0
		runvplot(uvdata, fringlist, 'HALF', timer, refantlist, refantlist, 1, nif, bchan, echan, -1, 0, 0, -1, 0, aparm, bparm, refantn, dotv, 9)
		# LWPLA to disk
		outfile = os.path.join(plotdir, 'VPLOT_NOCAL_' + strftime("%y%m%d_%H%M") + '.ps')
		print 'Plotting VPLOT diagnostics:' + outfile
		runlwpla(uvdata, outfile)




		print '***********************************************************'
		print '    Finished plotting pre-calibration diagnostic plots.'
		print '         It would be wise to inspect them closely.'
		print '***********************************************************'


	else:
		print "*** Plot pre-calibration diagnostics not selected.. skipping ***"



	if dofringfit == 1:
		print "------------------------------------------------------------"
		print "            Fringe fitting e-MERLIN data                    "
		print "           dofringfit = 1 (part of docalib)                 " 
		print "------------------------------------------------------------"

		# Initial fringe fit on everything *but* the target(s)
		aparm = 3, 0, 0, 0, 2, 1, 2.5, 0, 0, 0
		dparm = 3, 450, 0, 0, 0, 0, 0, 1, 1, 0 
		timer = 0,0,0,0,0,0,0,0
		solint = 10
		print 'Running initial fringe fit (SN#1)'
		newlist = list(set(fluxcals+bpasscals+pointcals+phsrefs))
		runfring(uvdata, newlist, timer, -1, 0, 0, -1, -1, refantn, refantlist, solint, aparm, dparm, 0, fring_snr, bchan, echan)

		# Plot Delay solutions (SNPLT)
		optype = 'DELA'
		snver=0
		snver = get_tab(uvdata, 'SN')
		print 'SN version' + format(snver) + 'contains fringe fit solutions for:', fringlist
		runsnplt(uvdata, 'SN', 1, fringlist, 9, optype, dotv)

		# LWPLA to disk
		plotfile = 'SNPLT_POSTFRING_' + strftime("%y%m%d_%H%M") + '.ps'
		outfile = os.path.join(plotdir, plotfile)
		print 'Plotting SNPLT diagnostics:' + outfile
		runlwpla(uvdata, outfile)



		# Apply solutions with CLCAL
		runclcal(uvdata, sourcelist, fringlist, 'CALI', 'AMBG', snver, snver, 0, 0, refantn)



		for i in range(len(fringlist)) :
			sourceID = fringlist[i],
			bparm = 0,0
			print 'Running POSSM on', fringlist[i]
			runpossm(uvdata, sourceID, timer, refantlist, refantlist, aparm, bparm, 0, nif, 100, 0, flagver, plotstokes, -1, 0, 'A&P', 30, 9, -1, 1)

			# LWPLA to disk
			plotfile = 'POSSM_POSTFRING_' + str(fringlist[i]) + '_' + strftime("%y%m%d_%H%M") + '.ps'
			outfile = os.path.join(plotdir, plotfile)
			print 'Plotting POSSM diagnostics:' + outfile
			runlwpla(uvdata, outfile)
			print 'CHECK YOUR INITIAL FRINGE FIT SOLUTIONS. You may need to edit SN table and/or re-run fring.'

	else :
		print "****************************************************************"
		print "*** Fringe fitting was not selected................ skipping.***"
		print "****************************************************************"




	if dosetflux == 1 :          # njj
		print "------------------------------------------------------------"
		print "            Setting flux density scale for e-MERLIN data    "
		print "                 dosetflux = 1 (part of docalib)            " 
		print "------------------------------------------------------------"          

		new = list(set(fluxcals+bpasscals+pointcals))
		clver = get_tab(uvdata, 'CL')
		print "Assuming fring fitting solutions are in CLver", clver, "(or have been already applied with SPLAT etc)"
		for i in range(len(fluxcals)):
			print 'Running SETJY on', fluxcals[i]

			# Now we calculate them properly:
			# First extract IF frequencies from the header and overwrite subbandfreq
			fluxcalflux = []
			subbandfreq = []
			fqtab=uvdata.table('FQ',1)
			for ifstep in fqtab[0].if_freq :
				setfreq = (uvdata.header.crval[2] + ifstep) / 1000000
				subbandfreq.append(setfreq)
				# Then calculate the appropriate flux and overwite fluxcalflux
				fluxcalflux.append(dfluxpy(setfreq,uvdata))
			print fluxcalflux, 'Jy per Sub-band [calculated]'
			print subbandfreq, 'Central frequency of each sub-band [calculated]'





			for j in range(len(fluxcalflux)):
				thing = float(fluxcalflux[j])
				thingfreq = float(subbandfreq[j])
				print "Setting flux density for source " + str(fluxcals[i]) + ", Sub-band", str(j+1),'(',thingfreq,'MHz) =', thing, 'Jy'
				runsetjy(uvdata, fluxcals, j+1, j+1, thing, '')
				i = 0
				j = 0
		for i in range(len(pointcals)):
			print 'Running SETJY on', pointcals[i]
			for j in range(len(pointflux)):
				thing = float(pointflux[j])
				thingfreq = float(subbandfreq[j])
				print "Setting flux for source " + str(pointcals[i]) + ", Sub-band", str(j+1),'(',thingfreq,'MHz) =', thing, 'Jy'
				runsetjy(uvdata, pointcals, j+1, j+1, thing, '')

	if dofluxcal == 1 :          # njj
		print "------------------------------------------------------------"
		print "            Running CALIB on fluxcal     "
		print "                 dofluxcal = 1 (part of docalib)            " 
		print "------------------------------------------------------------"

		# Run calib - P&A on FLUX cal with limited uvrange
		aparm = 4, 0, 0, 0, 0, 3, 5, 0, 0, 0
		cparm = 10, 0, 0, 0, 0, 0, 0
		timer = 0,0,0,0,0,0,0,0
		solint=10
		x = nif / 2
		uvrang = 0, float(subbandfreq[x])*0.07*10
		antwt = 0, 0
		clver = get_tab(uvdata, 'CL')
		snver = get_tab(uvdata, 'SN') + 1
		print "Running A&P calib on", new, "with CLver", clver
		print "Using limited UVrange of " + str(uvrang) + " klambda, derived from your central frequency"
		print "Writing SN table #", snver
		runcalib(uvdata, new, timer, uvrang, 100, clver, 0, -1, 0, 'DFT', refantn, solint, aparm, 1, 'L1', 'A&P', 10, 10, cparm, snver, antwt, 1)
		print "Re-run A&P calib using all baselines on", caliblist, "with CLver", clver
		print "Appending solutions to SN table #", snver
		uvrang = 0, 0
		runcalib(uvdata, caliblist, timer, uvrang, 100, clver, 0, -1, 0, 'DFT', refantn, solint, aparm, 1, 'L1', 'A&P', 10, 10, cparm, snver, antwt, 1)

		snver = get_tab(uvdata, 'SN')
		#for phscal in phsrefs :
		for i in range(len(caliblist)) :
			sourceID = caliblist[i],
			print "Running GETJY on caliblist", sourceID[0]
			rungetjy(uvdata, sourceID, fluxcals, 1, nif, snver)

		# Run clcal
		snver = get_tab(uvdata, 'SN')
		print "Running clcal on snver", snver, "with", caliblist
		runclcal(uvdata, '', caliblist, 'CALI', 'AMBG', snver, snver, 0, 0, refantn)

		print "Attempted to set FLUX density scale using 3C286 and bootstrapping to other calibrators - CHECK YOUR RESULTS in your SU TABLE"

	else :
		print "You opted not to derive flux scale... assuming you want to do this manually"



	if dobandpass == 1 :
		print "------------------------------------------------------------"
		print "            Doing a bandpass of your e-MERLIN data "
		print "                 dobandpass = 1 (part of docalib)                   "
		print "------------------------------------------------------------"   

		#First run SOUSP to derive a source spectral index for bpass
		print 'Running SOUSP.'
		bpasssrc = bpasscals[0],
		specindx = runsousp(uvdata, bpasssrc)

		# Bandpass solutions
		soltyp = 'L1'
		solint = -1
		docal = 100
		bpassprm = 0, 0, 0, 0, 1, 0, 0, 0, 1, 3, 1 
		print 'Running BPASS on', uvdata
		bpasssrc = bpasscals[0],
		runbpass(uvdata, bpasssrc, refantn, bpassprm, soltyp, solint, docal, specindx)

		timer = 0,0,0,0,0,0,0,0
		aparm = 0,0,0,0,0,0,0,0,0,1
		bparm = 0,0,0,0,0,0,0,0,0,0
		for i in range(len(caliblist)) :
			sourceID = caliblist[i],
			print 'Running POSSM on', caliblist[i]
			runpossm(uvdata, sourceID, timer, refantlist, refantlist, aparm, bparm, 0, nif, 1, 0, 0, plotstokes, 1, 1, 'A&P', 30, 9, -1, 1)

		# LWPLA to disk
		plotfile = 'POSSM_BANDPASS_' + strftime("%y%m%d_%H%M") + '.ps'
		outfile = os.path.join(plotdir, plotfile)
		print 'Plotting POSSM diagnostics:' + outfile
		runlwpla(uvdata, outfile)



	if dotherest == 1 :
		step = 0
		# Run calib - phase-only on phase and point source cals
		aparm = 3, 0, 0, 0, 0, 0, 5, 0, 0, 0
		cparm = 10, 0, 0, 0, 0, 0, 0
		timer = 0,0,0,0,0,0,0,0
		solint=10
		uvrang = 0, 0
		antwt = 0, 0
		clver = get_tab(uvdata, 'CL')
		print "Attempting phase-only calibration with CALIB", caliblist, "and clver", clver
		runcalib(uvdata, caliblist, timer, uvrang, 100, clver, 0, -1, 0, 'DFT', refantn, solint, aparm, 1, 'L1', 'P', 10, 10, cparm, 0, antwt, 1)

		# Run clcal
		snver = get_tab(uvdata, 'SN')
		print "Running clcal on snver", snver, "with", caliblist
		runclcal(uvdata, '', caliblist, 'CALI', 'AMBG', snver, snver, 0, 0, refantn)

		# Run calib - P&A on phase cal and point source cal
		aparm = 3, 0, 0, 0, 0, 3, 5, 0, 0, 0
		cparm = 10, 0, 0, 0, 0, 0, 0
		timer = 0,0,0,0,0,0,0,0
		solint=10
		uvrang = 0, 0
		antwt = 0, 0
		clver = get_tab(uvdata, 'CL')
		print "Running A&P calib on", caliblist, "with clver", clver
		runcalib(uvdata, caliblist, timer, uvrang, 100, clver, 0, -1, 0, 'DFT', refantn, solint, aparm, 1, 'L1', 'A&P', 10, 10, cparm, 0, antwt, 1)


		# runpossm() using solint ~30mins on all cals.
		timer = 0,0,0,0,0,0,0,0
		aparm = 0,0,0,0,0,0,0,0,0,1
		bparm = 0,0,0,0,0,0,0,0,0,0
		for i in range(len(caliblist)) :
			sourceID = caliblist[i],
			print 'Running POSSM on', caliblist[i]
			runpossm(uvdata, sourceID, timer, refantlist, refantlist, aparm, bparm, 0, nif, 1, 0, 0, plotstokes, -1, 0, 'A&P', 30, 9, -1, 1)

		# LWPLA to disk
		plotfile = 'POSSM_CALIB_' + strftime("%y%m%d_%H%M") + '.ps'
		outfile = os.path.join(plotdir, plotfile)
		print 'Plotting POSSM diagnostics:' + outfile
		runlwpla(uvdata, outfile)


		# Plot SNPLT diagnostics
		clver = get_tab(uvdata, 'CL')
		runsnplt(uvdata, 'CL', clver, caliblist, nplots=9, optype='PHAS', dotv=dotv)
		plotfile = 'SNPLT_CALIB_' + strftime("%y%m%d_%H%M") + '.ps'
		outfile = os.path.join(plotdir, plotfile)
		print 'Plotting SNPLT diagnostics:' + outfile
		runlwpla(uvdata, outfile)




		# Do selfcal if requested
		# For now, just assume we have only one phase cal and one source.
		imver = 1
		if doselfcal :
			mode = 'P'
			print "Self-calibration requested. Repeats:", doselfcal
			# Setup the base table for the start of each calibrator
			clverstart = get_tab(uvdata, 'CL')

			for i in range(doselfcal):
				j = i
				step = step + 1
				if i == (doselfcal - 1):
					mode = 'a&p'

				# image phase ref
				print "Mapping phase cal", phsrefs[0]
				calimsize = 256,256
				sourceID = phsrefs[0],
				runimagr(uvdata, sourceID, 1, 0, 0, 1, 0, bchan, echan, nchan, 1, cellsize, calimsize, 500, -1, indisk)


				aparm = 3, 0, 0, 0, 0, 0, 5, 0, 0, 0
				cparm = 10, 0, 0, 0, 0, 0, 0
				timer = 0,0,0,0,0,0,0,0
				solint=10
				uvrang = 0, 0
				antwt = 0, 0
				clver = get_tab(uvdata, 'CL')
				loop = i + 1
				source = str(sourceID[0])
				if len(source)>12 :
					source = source[0:12]
				uvdata2 = AIPSImage(source,'ICL001',indisk,loop)
				print "Attempting self-calibration with", source, "clver", clver, " and selfcal mode:", mode
				runselfcalib(uvdata, uvdata2, sourceID, timer, uvrang, 100, clver, 0, 1, 0, 'DFT', refantn, solint, aparm, 1, 'L1', mode, 10, 10, cparm, 0, antwt, 1)

				# Run clcal
				snver = get_tab(uvdata, 'SN')
				print "Running clcal on snver", snver, "with", sourceID[0]
				runclcal(uvdata, sourceID, sourceID, 'CALI', 'AMBG', snver, snver, 0, 0, refantn)

				# Plot VPLOT diagnostics
				gainuse = get_tab(uvdata, 'CL')
				runsnplt(uvdata, 'CL', gainuse, caliblist, nplots=9, optype='PHAS', dotv=dotv)
				plotfile = 'SNPLT_SELFCAL_' + str(step) + strftime("%y%m%d_%H%M") + '.ps'
				outfile = os.path.join(plotdir, plotfile)
				print 'Plotting SNPLT diagnostics:' + outfile
				runlwpla(uvdata, outfile)


			imver = j




#################################################################################################################
# Calculate SEFDs for each antenna and IF 
#################################################################################################################

if dosefd == 1 :

	print "------------------------------------------------------------"
	print " STAGE 3 SEFD diagnostics for e-MERLIN data (CALIB)"
	print "                        dosefd = 1                   " 
	print "------------------------------------------------------------"

	pca = AIPSCat(indisk)
	if len(pca[indisk]) == 1 :
		uvdata = AIPSUVData(pca[indisk][0]["name"], pca[indisk][0]["klass"], indisk, pca[indisk][0]["seq"])
		print "Calibrating:", uvdata
	else :
		uvdata = AIPSUVData('ALL', 'DBCON', indisk, 1)
		if uvdata.exists() :
			pass
		else:
			print "################################################################"
			print "Error: ALL.DBCON.1 does not exist! Don't know what to calibrate."
			print "################################################################"
			sys.exit()

	print "------------------------------------------------------------"
	print "      Running UVHGM averaging all IFs between BIF="+ format(SEFDbif) + " and EIF=" + format(SEFDeif)   
	print "                       mode = 1                            "
	print "------------------------------------------------------------"

	print "Running UVHGM on " + format(uvdata.name)+  '.' +format(uvdata.klass)+ '.' + format(uvdata.seq)

	antlist = []
	(refant, refantn, antlist) = set_refant(uvdata, refant)
	print "Antenna list : " +format(antlist)
	print "Number of Visibilities in data = " + format(len(uvdata))
	#print "Observation date = " + format(AIPSUVData.header.date-obs)
	#print uvdata.header.OBSRA

	logfileSEFD = os.path.join(plotdir, 'eMERLIN_SEFDSANTENNA.log')
	sefdfile = os.path.join(plotdir, 'SEFDs.log')
	sefdtxt = open(sefdfile, 'w')
	AIPS.log = open(logfileSEFD, 'w')
	sefdtxt.write("##########################################################\n")
	thisif = SEFDbif
	if SEFDeif > uvdata.header.naxis[3] :
		SEFDeif = uvdata.header.naxis[3]
		print "WARNING: you requested more IFs than I can find.  Setting SEFDeif to", SEFDeif

	uvdata.zap_table('PL', -1)
	finalsefds = zeros((10,4))
	while thisif <= SEFDeif :
		results = zeros((21,4))	# For individual baseline results
		sefds = zeros((7,15))	# For final SEFD results per antenna
		index = 0
		# Loop run of UVHGM over all baselines
		for i in range(len(antlist)):
			for j in range(i, len(antlist)-1):	# -1 here, and j+1 below, to avoid trying to calculate for not-baselines
				ant = antlist[i],0
				basel = antlist[j+1],0
				runuvhgm(uvdata, targets, ant, basel, thisif, thisif, SEFDbchan, SEFDechan)
				INPUTFILE = open(logfileSEFD, "rw").readlines()
				#RMSant0 = []
				for line in INPUTFILE:
					# Ugly fix to get around "CLRMSG" lines printed to the log (contains the string "RMS"!)
					if "UVHGM" in line:
						if "RMS" in line:
							GauRMS = float(re.findall(r'\S+',line)[4])
					elif "HISTOGRAM" in line:
						GauRMS = numpy.nan

				ijsefd = math.sqrt(GauRMS**2 - Sfield**2)*2*nu*kappa*math.sqrt(tau*beta)
				results[index] = (antlist[i],antlist[j+1],GauRMS,ijsefd)

				index = index + 1

		# Loop over triangles
		thing = []
		antlist = []
		(refant, refantn, antlist) = set_refant(uvdata, refant)
		antlist.sort()
		master = tuple(antlist)
		for i in range(len(master)):
			n=0
			antlist = []
			(refant, refantn, antlist) = set_refant(uvdata, refant)
			antlist.sort()
			ant1 = master[i]
			antlist.remove(ant1)
			bslist = []
			bslist.extend(itertools.combinations(antlist,2))
			sefd = []
			for j in range(len(bslist)):
				ant2,ant3 = bslist[j]
				bsl1 = 0
				bsl2 = 0
				bsl3 = 0
				for index in range(len(results)):
					if not numpy.isnan(results[index,3]):
						if (results[index,0] == ant1) or (results[index,1] == ant1):
							if (results[index,0] == ant2) or (results[index,1] == ant2):
								bsl1 = results[index,3]
						if (results[index,0] == ant1) or (results[index,1] == ant1):
							if (results[index,0] == ant3) or (results[index,1] == ant3):
								bsl2 = results[index,3]
						if (results[index,0] == ant2) or (results[index,1] == ant2):
							if (results[index,0] == ant3) or (results[index,1] == ant3):
								bsl3 = results[index,3]
				if bsl1 > 0:
					if bsl2>0:
						if bsl3>0:
							n = n + 1
							nargh = (bsl1 * bsl2) / bsl3
							sefd.append(nargh)

			if n>0:
				sefd = numpy.array(sefd)
				print "average sefd for antenna", ant1,"in IF",thisif,"is",nanmean(sefd)
				mystring = "average sefd for antenna " + format(ant1) + " in IF " + format(thisif) + " is " + format(nanmean(sefd)) + "\n"
				sefdtxt.write(mystring)
				finalsefds[ant1,0] = ant1
				finalsefds[ant1,1] = finalsefds[ant1,1] + nanmean(sefd)
				finalsefds[ant1,2] = finalsefds[ant1,2] + 1


		thisif = thisif + 1
		sefdtxt.write("##########################################################\n")

	for index in range(len(master)):
		finalsefds[master[index],3] = finalsefds[master[index],1] / finalsefds[master[index],2]

	print finalsefds

	for index in range(len(master)):
		print "Overall average SEFD for antenna ", format(int(finalsefds[master[index],0])), " is ", format(finalsefds[master[index],3]), "\n"
		mystring = "Overall average SEFD for antenna " + format(int(finalsefds[master[index],0])) + " is " + format(float(finalsefds[master[index],3])) + "\n"
		sefdtxt.write(mystring)
		
	sefdtxt.write("##########################################################\n")

	for index in range(len(master)):
		antsefd = finalsefds[refantn,3] / finalsefds[master[index],3]
		print "Overall scaled SEFD for antenna ", format(int(finalsefds[master[index],0])), " is ", format(antsefd), "\n"
		mystring = "Overall scaled SEFD for antenna " + format(int(finalsefds[master[index],0])) + " is " + format(float(antsefd)) + "\n"
		sefdtxt.write(mystring)

	plotfile = 'SEFD_UVHGM_' + strftime("%y%m%d_%H%M") + '.ps'
	outfile = os.path.join(plotdir, plotfile)
	print 'Plotting SEFD histogram diagnostics:' + outfile
	runlwpla(uvdata, outfile)



	print '***********************************************************'
	print 'Finished plotting STAGE 3 SEFD diagnostic.                 '
	print '      Your SEFDs are stored in: ', sefdfile
	print '      It would be wise to inspect them closely.            '
	print '      Note SEFD numbers are only indcative.                '
	print '      Values should be squared for use in WTMOD            '
	print '***********************************************************'






#################################################################################################################
# Final imaging run
#################################################################################################################

if docalib == 1 :

	print '***********************************************************'
	print 'End of calibration procedures. Mapping.'
	print '***********************************************************'

	pca = AIPSCat(indisk)
	if len(pca[indisk]) == 1 :
		uvdata = AIPSUVData(pca[indisk][0]["name"], pca[indisk][0]["klass"], indisk, pca[indisk][0]["seq"])
		print "Calibrating:", uvdata
	else :
		uvdata = AIPSUVData('ALL', 'DBCON', indisk, 1)
		if uvdata.exists() :
			pass
		else:
			print "################################################################"
			print "Error: ALL.DBCON.1 does not exist! Don't know what to calibrate."
			print "################################################################"
			sys.exit()


	print "Mapping phase cal(s)", phsrefs
	calimsize = 512,512
	for i in range(len(phsrefs)) :
		sourceID = phsrefs[i],
		runimagr(uvdata, sourceID, 1, 0, 0, 1, 0, bchan, echan, nchan, 1, cellsize, calimsize, 500, -1, indisk)
		source = str(sourceID[0])
		if len(source)>12 :
			source = source[0:12]
		image = AIPSImage(source,'ICL001',indisk,imver)
		runkntr(image,1,AIPS.userno)
		plotfile = sourceID[0] + str(imver) + '_' + strftime("%y%m%d_%H%M") + '.ps'
		outfile = os.path.join(plotdir, plotfile)
		runlwpla(image,outfile)

	print "Mapping target(s)", targets
	for i in range(len(targets)) :
		sourceID = targets[i],
		runimagr(uvdata, sourceID, 1, 0, 0, 1, 0, bchan, echan, nchan, 1, cellsize, imsize, 500, -1, indisk)
		source = str(sourceID[0])
		if len(source)>12 :
			source = source[0:12]
		image = AIPSImage(source,'ICL001',indisk,1)
		runkntr(image,5,AIPS.userno)
		plotfile = sourceID[0] + '_' + strftime("%y%m%d_%H%M") + '.ps'
		outfile = os.path.join(plotdir, plotfile)
		runlwpla(image,outfile)



	print '************************************************************'
	print "Finished. Now admire the awesomeness of your e-MERLIN maps."
	print "Your diagnostic plots and pipelined images can be found in:"
	print plotdir
	print '------------------------------------------------------------'
	print 'Note that it is highly unlikely that your images will be of'
	print 'publication quality. It is recommended that you refine the '
	print 'flagging and calibration by hand and re-image your target(s)'
	print '************************************************************'



print "Your AIPS catalog now looks like this:"
print AIPSCat(indisk)


# Finished whatever was requested, print the time taken.
print "Completed the requested procedures.  Check your results carefully!"
print "Time taken (hh:mm:ss):", time2hms(time.time()-ti)
print 'Finished on %s' % strftime("%d-%b-%y"), 'at:', strftime("%H:%M:%S", localtime()),'\n'

