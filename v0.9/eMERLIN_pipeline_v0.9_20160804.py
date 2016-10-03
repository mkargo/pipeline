# The e-MERLIN data reduction pipeline
#    Copyright (C) 2013  Megan Argo
#    Bug fixes, modifications and extra bits of code also by: 
#       Neal Jackson, Rob Beswick, Nick Wrigley, Javier Moldon, Danielle Fenech, 
#       Pierre Emmanuel Belles, Anthony Rushton, Adam Avison
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
# 20160504 - Merging in new QA code (JM), new widefield code (NW), cillitbang peeler (NJJ) model-based amplitude calibration and time-order checks in the doload section (MKA).
# 20160311 - Added new inputs for updated SERPent package, moved to run after doconcat section.
# 20150902 - Bug-fixed new doload section.
# 20150810 - Added fix for the Cambridge axis offset problem (CLCOR and SPLAT).
# 20150729 - Major re-write of the doload section to solve mixed-mode issues.
# 20140813 - Added a check in the DBCON stage to make sure there are files to combine.
# 20140513 - Added pol line in SERPent inputs, NJJ's fix for multiple files per directory, and
#	Nick Wrigley's fast wide-field imaging routines.
# 20140506 - Updated dfluxpy() to properly calculate projected baselines, and to use 2012 coefficients.
#	Also added in Neal Jackson's snedit routine and associated things.
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
from AIPS import AIPS, AIPSDisk
from AIPSTask import AIPSTask, AIPSList
from AIPSData import AIPSUVData, AIPSImage, AIPSCat
from Wizardry.AIPSData import AIPSUVData as WizAIPSUVData
from eMERLIN_tasks import *
import eMERLIN_tasks
import math, time, datetime
from numpy import *
import itertools, socket, glob
from time import gmtime, strftime, localtime
ti = time.time()

INDE = 3140.892822265625	# value corresponding to aips INDE

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
	global project, targets, phsrefs, fluxcals, pointcals, pointflux, bpasscals
	global refant, refantnames, band, imsize, doback
	global solint, snver, fring_snr, setjy_fluxes
	global nchan, bchan, echan, nif, stokes, plotstokes, polarizations, dotherest
	global doload, dosingleload, doavg, doflag, doflagmask, doconcat, docalib, dotidy, doback	# process management
	global dosort, doselfcal, tint, dosortfirst, dodiagnostic1, dodiagnostic3, dosefd
	global sbandl, sbandu, smootha, smoothb, smoothc, dodelete, dofringfit, dosetflux, dofluxcal, dobandpass, dosnedit, doimage, doqa, qaanten, qabasel
	global 	Sfield, beta, kappa, nu, tau, SEFDbif, SEFDeif, SEFDbchan, SEFDechan	# SEFD options
	global dowide, dosplit, dochop, domap, doflatn, dopbcor, dosad, fovradius, noiter, widetargets, imgsize, dopeel	# wide-field imaging variables
	global mixed

	AIPS.userno = int(control.get('userno',[0])[0])
	indisk = int(control.get('indisk', [0])[0])
	dosort = int(control.get('dosort',[1])[0])
	doselfcal = int(control.get('doselfcal',[0])[0])
	project = control.get('project',['none'])[0]

	# default is no averaging
	nchan = int(control.get('nchan',[-1])[0])
	tint  = int(control.get('tint',[-1])[0])

	# default is to run all procedures
	doload = int(control.get('doload',[-1])[0])
	dosingleload = int(control.get('dosingleload',[-1])[0])
	doavg  = int(control.get('doavg',[1])[0])
	doflag = int(control.get('doflag',[1])[0])
	doflagmask = int(control.get('doflagmask',[1])[0])
	doconcat = int(control.get('doconcat',[1])[0])
	docalib = int(control.get('docalib',[1])[0])
	dobandpass = int(control.get('dobandpass',[1])[0])
	dotidy = int(control.get('dotidy',[1])[0])
	dosortfirst = int(control.get('dosortfirst',[-1])[0])
	dodiagnostic1 = int(control.get('dodiagnostic1',[-1])[0])
	dosnedit = int(control.get('dosnedit',[-1])[0])
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
	dowide = int(control.get('dowide',[-1])[0])
	doimage = int(control.get('doimage',[-1])[0])
	doqa = int(control.get('doqa',[-1])[0])
	qaanten = control.get('qaanten',[])
	qabasel = control.get('qabasel',[])
	mixed = False


	# wide-field imaging options
	dosplit = int(control.get('dosplit',[1])[0])
	dochop = int(control.get('dochop',[1])[0])
	dopeel = int(control.get('dopeel',[-1])[0])
	domap = int(control.get('domap',[1])[0])
	doflatn = int(control.get('doflatn',[1])[0])
	dopbcor = int(control.get('dopbcor',[1])[0])
	dosad=int(control.get('dosad',[1])[0])
	fovradius = float(control.get('fovradius',[7.5])[0])
	noiter = int(control.get('noiter',[-1])[0])
	widetargets= control.get('widetargets', [])
	for i in range(len(widetargets)):
		widetargets[i] = widetargets[i].upper()
	

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
	print "doback = "+ format(doback)
	print "doconcat = "+ format(doconcat)
	if doconcat == 1 :
		print "With:"
		print "             Deleting work files after combining(1=yes) =" + format(dodelete)
	print "doflag = "+ format(doflag)
	if doflag == 1 :
		print "             OK: will run the SERPent autoflagger on your files"
	print "doqa = "+ format(doqa)
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
		print "              Final spec chan width (Hz) = " + format(beta)
		print "              widar corrlator efficieny (esitimate) = " + format(nu)
		print "              Hanning smoothing SNR improvement = " + format(kappa)
	print "dowide = "+ format(dowide)
	if dowide == 1 :
		print "With sub-sections:"
		print "             **Targets to map = " + format(widetargets)
		print "             **Split each target from the main file = " + format(dosplit)
		print "             **Chop and average each target into chessboard = " + format(dochop)
		print "             **Run clean on each facet = " + format(domap)
		print "             **Trim and flatn the 64 files into one wide field = " + format(doflatn)
		print "             **Beam correct the images to typical parameters = " + format(dopbcor)
		print "             **Resample, create noise map, and search using SAD = " + format(dosad)
		print "             **With field-of-view radius in arcminutes = " + format(fovradius)
		print "             **Number of clean iterations per facet = " + format(noiter)



	print "***********************************************************"
	print AIPSCat(indisk)
	print "***********************************************************"





def setup(uvdata, dosort, dotidy):
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

	chanwidth = uvdata.header.cdelt[2]	# in Hz
	nchan = uvdata.header.naxis[2]		# assumes continuum ONLY!
	# ignore the outer 10% of channels for certain tasks
	bchan = int(math.floor(nchan * 0.1))
	echan = nchan - bchan
	nif = uvdata.header.naxis[3]
	ifwidth = chanwidth*nchan

	freq = uvdata.header.crval[2] / 1000000000	# This may be wrong if the user has run DBCON themselves with FQCENTER=1 (but I'm not sure it matters) 20140522
	band = ''
	if (freq > 1.3) and (freq < 1.7):
		band = 'L'
	if (freq > 4) and (freq < 8):
		band = 'C'
	if (freq > 18) and (freq < 26):
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
		runlwpla(uvdata, outfile, 0)

	# Plot phase on cals (UVPLT, VPLOT)
	aparm = 0,0
	bparm = 12,2,0
	runvplot(uvdata, fringlist, 'I', timer, refantlist, refantlist, 1, nif, bchan, echan, -1, 0, 0, -1, 0, aparm, bparm, refantn, dotv, 9)

	# LWPLA to disk
	filename = 'VPLOT_' + str(step) + '_phase.ps'
	outfile = os.path.join(plotdir, filename)
	print 'Plotting VPLOT diagnostics:' + outfile
	runlwpla(uvdata, outfile, 0)


def flagger(name, klass, seq, indisk):
	''' create directory, write control file and run SERPent '''
	mkdir = 'serpent_' + name + '_' + klass + '_' + str(seq)
	os.system('mkdir ' + mkdir)
	path2folder = os.path.join(os.getcwd(), mkdir)

	# run the flagger
	os.chdir(path2folder)
	try: os.system('cp ../SERPent.py .')
	except IOError: print "Error: can't find SERPent!"
	try: os.system('cp ../SERPent_input.py .')
	except IOError: print "Error: can't find SERPent inputs file!"
	# Work-around for Jodrell numpy problems
	host = socket.gethostname()
	if 'compute-0' in host:
		print 'WARNING: I think you are on the JBO PHOENIX cluster, running parseltongue.broken'
		try: os.system('parseltongue.broken SERPent.py')
		except: print "Error: unknown SERPent problem!"
	else:
		try: os.system('parseltongue SERPent.py')
		except: print "Error: unknown SERPent problem!"

	os.chdir('../')


def nanmean(data,**args):
	return numpy.ma.filled(numpy.ma.masked_array(data,numpy.isnan(data)).mean(**args),fill_value=numpy.nan)


# ----  njj: functions for the automatic editor ------------------------
#  NB need Wizardry.AIPSData for this since it changes tables
#  NB  the way to update elements of a table sn[0] is:
#      i=sn[0]; i['delay_1']=1; i.update()
#      AND NOT sn[0]['delay_1']=1; sn[0].update() which DOES NOT WORK!!
#      (if anyone knows why, please write in here....)
#      This means there is no way of dealing with the case of one IF except
#      an explicit if block, which is what is done.

def ifedit (ain,clip):
	a=numpy.copy(ain)
	aval = a[(a!=0.)&(~numpy.isnan(a))]
	if len(aval)<3:
		return a
	# remove discrepant values
	medscat = numpy.median(abs(numpy.diff(aval)))
	fval = numpy.array([False]*len(aval))
	if abs(aval[0]-aval[1])>clip*medscat:
		fval[0]=True
	if abs(aval[-1]-aval[-2])>clip*medscat:
		fval[-1]=True
	for i in range(1,len(aval)-1):
		if abs(aval[i]-aval[i-1])>clip*medscat and \
		   abs(aval[i]-aval[i+1])>clip*medscat:
			fval[i]=True
	aval=aval[fval]
	for i in range(len(a)):
		if a[i] in aval:
			a[i] = numpy.nan
	# interpolate over to remove nans/zeroes, except where this would
	# interpolate over delay jumps
	bad = numpy.where((a==0)|(numpy.isnan(a)))[0]
	badgroup = numpy.array_split(bad,numpy.where(numpy.diff(bad)!=1)[0]+1)
	#print badgroup
	for i in badgroup:
		if len(i) and i[0] and i[-1]!=a.shape[0]-1:
			bdiff = -a[i[0]-1]+a[i[-1]+1]
			if abs(bdiff) < 2.0*clip*medscat:
				for j in i:
					a[j] = a[i[0]-1]+bdiff*float(j-i[0]+1)/float(len(i)+1)
	return a

# edit the data. s is in nanoseconds with zero or INDE recorded as numpy.nan
def delay_edit(s,clip=3.0):
	nant,npol,nif,nt = s.shape
	so = numpy.zeros_like(s)
	for iant in range(nant):
		aflat = numpy.ravel(s[iant])
		if len(aflat) == len(aflat[aflat==0])+len(aflat[numpy.isnan(aflat)]):
			print 'Detected reference antenna: ',iant,': setting delay=0'
			continue
		for ipol in range(npol):
			for iif in range(nif):
				ain = numpy.copy(s[iant,ipol,iif])
				if not len(ain[~numpy.isnan(ain)]):
					print 'Failed: ant=',iant,'pol=',ipol,'if=',iif
					so[iant,ipol,iif]=numpy.copy(ain)
					continue
				aout = ifedit (ain,clip)
				so[iant,ipol,iif]=numpy.copy(aout)
	numpy.save('s',s)
	numpy.save('so',so)
	return so
	
def ifphase(a,clip):
	b=a[(a!=0.0)&(~numpy.isnan(a))]; c=[]
	# identify and set to nan all outlying points
	medscat = numpy.median(abs(numpy.diff(b)))
	for i in range(len(b)):
		if (i==0 and abs(b[0]-b[1])>clip*medscat) or (i==len(b)-1 and abs(b[-1]-b[-2])>clip*medscat):
			c=numpy.append(c,b[i])
		if i in [0,len(b)-1]:
			continue
		if abs(b[i]-b[i+1])>clip*medscat and abs(b[i]-b[i-1])>clip*medscat:
			c=numpy.append(c,b[i])
	for i in range(len(a)):
		if a[i] in c:
			a[i]=numpy.nan
	bad = numpy.where((a==0)|(numpy.isnan(a)))[0]
	badgroup = numpy.array_split(bad,numpy.where(numpy.diff(bad)!=1)[0]+1)
	for i in badgroup:
		if len(i) and i[0] and i[-1]!=a.shape[0]-1:
			bdiff = -a[i[0]-1]+a[i[-1]+1]
			if abs(bdiff) < 2.0*clip*medscat:
				for j in i:
					a[j] = a[i[0]-1]+bdiff*float(j-i[0]+1)/float(len(i)+1)
	return a

def phase_edit(s,clip):
	nant,npol,nif,nt = s.shape
	so = numpy.zeros_like(s)
	for iant in range(nant):
		for ipol in range(npol):
			for iif in range(nif):
				a=numpy.copy(s[iant,ipol,iif])
				if len(a) > 1:
					so[iant,ipol,iif]=ifphase(a,clip)
	return so

def edit_table(inna, incl, itype, inver=0, inseq=1, indist=1, indisk=1, clip=3.0):
	uvdata = AIPSUVData (inna,incl,indisk,inseq)
	tacop = AIPSTask ('tacop')
	tacop.indata = uvdata
	tacop.inext = 'SN'
	tacop.inver = inver
	tacop.outna = 'ALL'
	tacop.outcl = 'DBCON'
	tacop.outver = 0
	tacop.outdisk = indisk
	tacop.go()
	inver = 0	# have now copied to the last SN table, which we can edit

	uvdata = WizAIPSUVData(inna, incl, indisk, inseq)
	sntables = numpy.asarray(uvdata.tables,dtype='S')
	try:
		asntables = sntables[sntables[:,1]=='AIPS SN']
		nsntables = int(max(numpy.asarray(asntables[:,0],dtype='int')))
	except:
		print 'There are no SN tables - cannot edit'
		return

	if inver>nsntables:
		print 'The requested SN table does not exist - cannot edit'
		return

	inver = nsntables if not inver else inver
	print 'Amending table',inver
	sn = uvdata.table('SN',inver)

	try:
		nif = len(sn[0]['delay_1'])
	except:
		nif = 1
	antkey = numpy.array([]); tkey = numpy.array([])
	for i in sn:
		antkey = numpy.append(antkey,i['antenna_no'])
		tkey = numpy.append(tkey,i['time'])

	ant = numpy.sort(numpy.unique(antkey)); nant = len(ant)
	t = numpy.sort(numpy.unique(tkey)); nt = len(t)
	src = numpy.zeros((nant,nt))
	if itype=='delay':
		s = numpy.zeros((nant,2,nif,nt))
	elif itype=='phase':
		s = numpy.zeros((nant,2,nif,nt),dtype='complex')

	for i in sn:
		xant = numpy.argwhere(ant==i['antenna_no'])[0,0]
		xt = numpy.argwhere(t==i['time'])[0,0]
		if itype=='delay':
			s[xant,0,:,xt]=numpy.copy(i['delay_1'])
			s[xant,1,:,xt]=numpy.copy(i['delay_2'])
		elif itype=='phase':
			for j in range(s.shape[2]):
				if nif>1:
					s[xant,0,j,xt]=complex(i['real1'][j],i['imag1'][j])
					s[xant,1,j,xt]=complex(i['real2'][j],i['imag2'][j])
				else:
					s[xant,0,0,xt]=complex(i['real1'],i['imag1'])
					s[xant,1,0,xt]=complex(i['real2'],i['imag2'])
		src[xant,xt] = i['source_id']

	if itype=='delay':
		numpy.putmask(s,s>1.0,numpy.nan)
		so=1.0E-9*delay_edit(s*1.0E9,clip)
		numpy.putmask(so,numpy.isnan(so),INDE)
		for i in sn:
			xant = numpy.argwhere(ant==i['antenna_no'])[0,0]
			xt = numpy.argwhere(t==i['time'])[0,0]
			if nif > 1:
				for iif in range(nif):
					if i['delay_1'][iif] != INDE:
						i['delay_1'][iif]=so[xant,0,iif,xt]
						i.update()
					if i['delay_2'][iif] != INDE:
						i['delay_2'][iif]=so[xant,1,iif,xt]
						i.update()
			else:
				if i['delay_1'] != INDE:
					i['delay_1']=so[xant,0,0,xt]
					i.update()
				if i['delay_2'] != INDE:
					i['delay_2']=so[xant,1,0,xt]
					i.update()

	elif itype=='phase':
		numpy.putmask(s,s==complex(INDE,INDE),numpy.nan)
		so = phase_edit(s,clip)
		numpy.putmask(so,numpy.isnan(so),complex(INDE,INDE))
		for i in sn:
			xant = numpy.argwhere(ant==i['antenna_no'])[0,0]
			xt = numpy.argwhere(t==i['time'])[0,0]
			if nif>1:
				for iif in range(nif):
					i['real1'][iif],i['imag1'][iif] = \
						 so[xant,0,iif,xt].real,so[xant,0,iif,xt].imag
					i.update()
					i['real2'][iif],i['imag2'][iif] = \
						 so[xant,1,iif,xt].real,so[xant,1,iif,xt].imag
					i.update()
			else:
				i['real1'],i['imag1'] = \
					 so[xant,0,0,xt].real,so[xant,0,0,xt].imag
				i.update()
				i['real2'],i['imag2'] = \
					 so[xant,1,0,xt].real,so[xant,1,0,xt].imag
				i.update()
  
#   -------- njj: end functions for automatic editor ------



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
print "-----------------------------------------------------------"
print "     Please read and respond to the following prompts      "
print "-----------------------------------------------------------"

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

try: os.system('mkdir plots')
except: pass
try: os.system('mkdir fits')
except: pass

# Check the inputs - must do this to set AIPS.userno *BEFORE* task definitions!
checkin(control)

# Changed the defaults of inputinput() so that users have to type something! 20150804
inputinput('Continue? (y/[n]): ')

if (doavg != 1) or ((doavg == 1) and ((nchan == -1) and (tint == -1))):
	inputinput("You requested NO averaging. Are you ABSOLUTELY SURE? (y/[n]): ")

#if docalib != -1:
#	inputinput('If you observed in mixed-mode, you MUST separate your data before calibrating your data! Proceed? (y/[n]): ')


if (dowide != -1):
	inputinput("You requested some or all of the widefield imaging procedures. Is your AIPS catalogue in the expected state? (y/[n]): ")


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
	print " Loading e-MERLIN data from disk and doing some houekeeping "
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
				# Use softlink to get around long datafile name problems
				os.system('rm EMTEMP.FITS')
				os.system('ln -s '+descriptor+' EMTEMP.FITS')
				runfitld('./EMTEMP.FITS', indisk, thisdir)

				# Force sort order to T* (as it is written this way)
				print 'Setting sortorder to T*'
				uvdata = WizAIPSUVData('TMP','UVDATA',outdisk,thisdir)
				uvdata.header.sortord='T*'
				uvdata.header.update()

				uvdata = AIPSUVData('TMP','UVDATA',outdisk,thisdir)
				# reset the filename based on the source name in the SU table
				if len(uvdata.sources) > 1 :
					print "Error: Too many sources in the file! Run on single-source files only."
					sys.exit()

				# index each file
				runindxr(uvdata)

				sourcelist = uvdata.table('SU',1)
				outn = re.sub(r"\s+", '-', uvdata.sources[0])	# catch source names with spaces
				if len(outn)>12 :
					outn = outn[0:12]			# catch long source names

				# to get around multiple files of a single source in one subdirectory
				seqno = 1
				pca = AIPSCat()
				for ipca in pca[indisk]:
					if ipca['name']==outn:
						seqno = ipca['seq']+1
				uvdata.rename(outn,'UVDATA',seqno)
				inname = outn


	# check each file, make sure they are all in the same mode
	mixed = False
	pca = AIPSCat(indisk)
	for fitsfil in pca[indisk]:
		uvdata1 = AIPSUVData(fitsfil.name, fitsfil.klass, indisk, fitsfil.seq)
		fqtab1 = uvdata1.table('FQ',1)
		for row in fqtab1:
			# Are we dealing with mixed mode
			mixed = checkEqual(row.ch_width)
		try:
			uvdata2
			fqtab2 = uvdata2.table('FQ',1)
			if fqtab1[0].ch_width == fqtab2[0].ch_width:
				ifstruc = fqtab1[0].ch_width
			else:
				print "Error: your datasets have different frequency setups, quitting."
				print "Consult the e-MERLIN staff for advice on how to deal with this."
				sys.exit()
			uvdata2 = AIPSUVData(fitsfil.name, fitsfil.klass, indisk, fitsfil.seq)
		except NameError:
			uvdata2 = AIPSUVData(fitsfil.name, fitsfil.klass, indisk, fitsfil.seq)

	print "All datasets have the same setup, continuing."


	# Deal with mixed-mode datasets by splat-ing and re-combining.
	if mixed:
		print "Your data are in mixed mode (sub bands of different widths)."
		print "Will now split your data and create a sane continuum database."
		pca = AIPSCat(indisk)
		for fitsfil in pca[indisk]:
			filelist = deque([])
			uvdata = AIPSUVData(fitsfil.name, fitsfil.klass, indisk, fitsfil.seq)
			fqtab = uvdata.table('FQ',1)
			# because on loading, AIPS gets the total bandwidths wrong!
			ifstruc = fqtab[0].ch_width
			ifs_tot = len(ifstruc)
			modes = set(ifstruc)
			contf = max(modes)
			print "The continuum width in your dataset is", contf
			for freq in modes:
				ns = ifstruc.count(freq)
				print "there are", ns, "instances of", freq
			start = 0
			stop = 0
			iflist = list(enumerate(ifstruc))
			# Make list of files to operate on, and then splat out IFs.
			for eachif in iflist:
				outcl = 'SPLT' + str(eachif[0]+1)
				runsplatmixed(uvdata, eachif[0]+1, eachif[0]+1, outcl, indisk)
				# test for line IF, move if so
				if eachif[1] == contf:
					fitout = (fitsfil.name, outcl, fitsfil.seq, eachif[1]) 
					filelist.append(fitout)
				else:
					lineuv = AIPSUVData(fitsfil.name, outcl, indisk, fitsfil.seq)
					runmove(lineuv, indisk)

			# move the original data
			runmove(uvdata, indisk)

			print "Files to combine:"
			print filelist
			vbgluloop(filelist, indisk, contf)





	pca = AIPSCat(indisk)
	for fitsfil in pca[indisk]:
		uvdata = AIPSUVData(fitsfil.name, fitsfil.klass, indisk, fitsfil.seq)

		# Check for time-ordering problems
		needsort = checkTime(uvdata)
		print "SORT CHECK: needsort = ", needsort
		if needsort:
			print "Warning: file not time-ordered, running UVSRT!"
			runuvsrt(uvdata, fitsfil.name, indisk, fitsfil.seq)
			runindxr(uvdata)

		# uvfix each file here, only if not averaging later
		if doavg != 1:
			runuvfix(uvdata, fitsfil.name, indisk, fitsfil.seq)

#		# index each file
		runindxr(uvdata)

		# flag the Lovell - Mk2 baseline (SERPent does it anyway, but not all users will make use of it).
		flagLOMK(uvdata)


if dosingleload == 1:
	print 'LOADING single fits file from: '+ fitsdir[0]
	#fitslist = [f for f in filter(nodot, os.listdir(fitsdir)) if os.path.isfile(os.path.join(fitsdir, f))]
	fitslist = glob.glob(fitsdir[0]+'/*fits')
	for fitsfile in fitslist:  # for each file...
		uvdata = AIPSUVData('MULTI', 'UVDATA', indisk, 1)
		runfitld2(outdata=uvdata, infile=fitsfile)
		runuvsrt(uvdata, uvdata.name, indisk, uvdata.seq)
		runindxr(uvdata)
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
#			if len(uvdata.sources) > 1 :
#				print "Error: Too many sources in the file! Run on single-source files only."
#				sys.exit()

			print 'Running SPLAT on ', fitsfil.name + '.' + fitsfil.klass + '.' + format(fitsfil.seq)
			runsplat(uvdata, nchan, tint, sbandl, sbandu, smootha, smoothb, smoothc, indisk)

			# Remove unaveraged file, then rename the averaged file
			uvdata.zap()
			uvdata = AIPSUVData(fitsfil.name,'SPLAT',indisk,fitsfil.seq)
			uvdata.rename(fitsfil.name,'UVDATA',fitsfil.seq)

			uvdata = AIPSUVData(fitsfil.name,'UVDATA',indisk,fitsfil.seq)
			print 'Running UVFIX on ', fitsfil.name + '.' + fitsfil.klass + '.' + format(fitsfil.seq)
			runuvfix(uvdata, fitsfil.name, indisk, fitsfil.seq)

			# think the outputs are also T* sorted, relabel them
			print 'Setting sortorder to T*'
			uvdata = WizAIPSUVData(fitsfil.name,'UVDATA',indisk,fitsfil.seq)
			uvdata.header.sortord='T*'
			uvdata.header.update()

		print "Your data are now averaged as requested. Next the data should be flagged with SERPent."

	else :
		print "You requested no averaging, skipping."


#################################################################################################################
# Optional dosort (not required now, but left in for completeness)
#################################################################################################################

if dosort == 1 :
	pca = AIPSCat(indisk)

	for fitsfil in pca[indisk]:
		uvdata = AIPSUVData(fitsfil.name, fitsfil.klass, indisk, fitsfil.seq)
		# Check for time-ordering problems
		needsort = checkTime(uvdata)
		print "In dosort() SORT CHECK: needsort = ", needsort
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
			runlwpla(uvdata, outfile, 0)

			runquickvplot(uvdata)
			plotfile =  fitsfil.name + '.' + '.' + str(fitsfil.seq) + '_PREAOF_VPLOT_' + strftime("%y%m%d_%H%M") + '.ps'
			outfile = os.path.join(plotdir,plotfile)
			print 'Plotting VPLOT diagnostics:' + outfile
			runlwpla(uvdata, outfile, 0)


	print '********************************************************************'
	print '    Finished plotting STAGE 1 pre-flag/calibration diagnostic plots.'
	print '         It would be wise to inspect them closely.                  '
	print '********************************************************************'





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
		print 'Running FITTP on ', fitsfil.name + '.' + fitsfil.klass + '.' + format(fitsfil.seq)
		runfittp(uvdata, fittpdir, fittpfile)


#################################################################################################################
# DBCON the data to create one dataset
#################################################################################################################

if doconcat == 1 :

	pca = AIPSCat(indisk)
	print "------------------------------------------------------------------"
	print "Combining your datasets to create one uv database for calibration."
	print "------------------------------------------------------------------"

	# Sanity check: do we have files to combine?
	uvcount = 0
	for fitsfil in pca[indisk]:
		if fitsfil.type == 'UV':
			uvcount = uvcount + 1

	if uvcount > 1:

		for fitsfil in pca[indisk]:
			uvdata = AIPSUVData(fitsfil.name, fitsfil.klass, indisk, fitsfil.seq)
			if len(uvdata.sources) > 1 :
				print "Error: Too many sources in the file! Run on single-source files only."
				sys.exit()


#		# Backup tables in case the user did any extra flagging
#		for fitsfil in pca[indisk]:
#			fittpfile =  fitsfil.name + '.TASAV.' + str(fitsfil.seq) + strftime("%y%m%d_%H%M") + '.fits'
#			uvdata = AIPSUVData(fitsfil.name, fitsfil.klass, indisk, fitsfil.seq)
#			runtasav(uvdata, fittpdir, fittpfile, indisk)

		print "First apply the HIGHEST FG table using UVCOP"
		for fitsfil in pca[indisk]:
			uvdata = AIPSUVData(fitsfil.name, fitsfil.klass, indisk, fitsfil.seq)
			if get_tab(uvdata, 'FG') != 0:
				fittpfile =  fitsfil.name + '.TASAV.' + str(fitsfil.seq) + strftime("%y%m%d_%H%M") + '.fits'
				uvdata = AIPSUVData(fitsfil.name, fitsfil.klass, indisk, fitsfil.seq)
				runtasav(uvdata, fittpdir, fittpfile, indisk)

				uvcop = AIPSTask('UVCOP')
				uvcop.indata = uvdata
				uvcop.outclass = 'UVCOP'
				uvcop.outdisk = indisk
				uvcop.outseq = fitsfil.seq
				uvcop.flagver = get_tab(uvdata, 'FG')
				print "Apply top flag table number", get_tab(uvdata, 'FG')
				try:
					uvcop.fqcenter = -1		# Added 20150918
				except AttributeError as detail:
					print "Your version of AIPS does not have this parameter,"
				uvcop.go()
				uvdata.zap()
				# Delete any remaining 'old' FG tables so DBCON is not confused
				uvdata = AIPSUVData(fitsfil.name, 'UVCOP', indisk, fitsfil.seq)
				uvdata.zap_table('FG', -1)
			else:
				uvdata.rename(fitsfil.name,'UVCOP',fitsfil.seq)


		# Combine the data with DBCON
		pca = AIPSCat(indisk)
		queueNAME = deque([])

		for fitsfil in pca[indisk]:
			tmp = (fitsfil.name, fitsfil.klass, indisk, fitsfil.seq)
			queueNAME.append(tmp)

		print "Total files to combine: ", len(queueNAME)

		i = 0
		while len(queueNAME)>1 :
			tmp1 = queueNAME.popleft()
			tmp2 = queueNAME.popleft()
			print "Combining", tmp1, "and", tmp2
			uvdata1 = AIPSUVData(tmp1[0], tmp1[1], indisk, tmp1[3])
			uvdata2 = AIPSUVData(tmp2[0], tmp2[1], indisk, tmp2[3])
			dbcon = AIPSTask('DBCON')
			dbcon.indata = uvdata1
			dbcon.in2data = uvdata2
			dbcon.outname = format(i)
			dbcon.outclass = 'DBCON'
			dbcon.outdisk = indisk
			dbcon.dopos = [None, [None, -1, -1], [None, -1, -1], [None, -1, -1], [None, -1, -1],]
			try:
				dbcon.fqcenter = -1		# Added 20140521 due to a change in AIPS
			except AttributeError as detail:
				print "Your version of AIPS does not have the fqcenter parameter, so this does not matter"
			dbcon.go()

			# in the interests of disk efficiency:
#			uvdata1.zap()
#			uvdata2.zap()		# If there is a problem during this stage of the pipeline, deleting files here means starting from scratch!

			indxr = AIPSTask('INDXR')
			indxr.indata = AIPSUVData(dbcon.outname, dbcon.outclass, indisk, dbcon.outseq)
			indxr.go()

			tmp = (format(i), 'DBCON', indisk, 1)
			queueNAME.append(tmp)

			i = i + 1


		# Run MSORT on the final DBCONd file
		print "Final sort and index...."
		msort = AIPSTask('MSORT')
		msort.indata = AIPSUVData(queueNAME[0][0], queueNAME[0][1], indisk, queueNAME[0][3])
		uvdata = AIPSUVData(queueNAME[0][0], 'MSORT', indisk, queueNAME[0][3])
		msort.outdata = uvdata
		msort.go()

		print "Reset all fluxes to zero"
		runsetjy(uvdata, '', 0, 0, 0, 'REJY')

		# Delete CL#1 and remake it here, since it is almost certainly incomplete after DBCON process.
		uvdata.zap_table('CL', -1)
		indxr.indata = uvdata
		indxr.go()
		uvdata.rename('ALL','DBCON',0)

		# Run CLCOR if needed to correct the Cambridge axis offset problem in data
		# before mid-2015, then splat data to apply correction.
		antab = uvdata.table('AN',1)
		for row in antab:
			if 'Cm' in row.anname:
				if row.staxof == 0:
					print "Doing Cambridge axis offset correction."
					runclcor(uvdata, row.nosta, indisk)

		# QUACK dBCONd data file 
		print "Then quack the data!..quack-quack-OOPS..FG#1 created.."
		quack = AIPSTask('QUACK')
		quack.indata = uvdata
		quack.opcode ='BEG'
		quack.reason = 'quack1min'
		quack.aparm[2] = 0.2
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
			print "You will probably want to do this later to save space."
			print "NOTE that any FLAGS have been applied."


	else:
		print "WARNING: no files to combine! DBCON stage skipped."



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
		try: apply(flagger, (uvdata.name, uvdata.klass, uvdata.seq, indisk,),)
		except NotImplementedError: print "ERROR: there appears to be a SERPent problem, quitting."
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
# Make plots for QA
#################################################################################################################


pid_list = []
solint_vplot = -1.5
#spplot_path  = '/home/jmoldon/software/SPPlot_230615_alt_jm.py'
spplot_path  = './SPPlot_230615_alt_jm.py'



if doqa == 1 :

	print "------------------------------------------------------------------"
	print "                   Producing summary and QA plots                 "
	print "------------------------------------------------------------------"

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

	# Data channel information
	sourcelist, fringlist, caliblist, bpasslist, nchan, nif, stokes, plotstokes, polarizations, bchan, echan, freq, band = setup(uvdata, 0, 0)
	(refant, refantn, refantlist) = set_refant(uvdata, refant)

	chanwidth = uvdata.header.cdelt[2]	  # in Hz
	nchan = uvdata.header.naxis[2]		  # assumes continuum ONLY!
	ifwidth = chanwidth*nchan

  

	# Run INDXR (duplicated if whole pipeline is run, but needed if not)
	print 'Running INDXR...'
	runindxr(uvdata)

	## Run UVFIX (duplicated if whole pipeline is run, but needed if not. Needed when exporting data from console)
	print 'Running UVFIX...'
	runuvfix(uvdata, uvdata.name, uvdata.disk, uvdata.seq)

	# scan list
	print("Producing scan list")
	outfile = plotdir + '/'+project+'.SCAN.txt'
	delfile(outfile)
	runlistr(indata=uvdata, optype='SCAN', outprint=outfile)

	# list of stations and number of vis.
	print("Producing data summary")
	outfile = plotdir + '/'+project+ '.DTSUM.txt'
	delfile(outfile)
	rundtsum(indata=uvdata, outprint=outfile)

	# possm plot of the cross-correlations
	print("Ploting POSSM uncalibrated")
	aparm = set_default_aparms('possm')
	plot_file = plotdir + '/'+project+ '.POSSM_CALIBRATORS_UNCAL.PS'
	uvdata.zap_table('AIPS PL', -1)
	src = numpy.unique(numpy.hstack([pointcals, fluxcals,bpasscals]))
	runpossm2(indata=uvdata, sources = fringlist, aparm=aparm, docalib = -1)
	plot(uvdata, plot_file, dopng=False)
	
	# Now only the target source (probably noise)
	aparm = set_default_aparms('possm')
	plot_file = plotdir + '/'+project+ '.POSSM_TARGET_UNCAL.PS'
	uvdata.zap_table('AIPS PL', -1)
	runpossm2(indata=uvdata, sources = targets, aparm=aparm, docalib = -1)
	plot(uvdata, plot_file, dopng=False)

	# vplot of the raw data
	print("Ploting VPLOT LOGamp uncalibrated")
	bparm = set_default_bparms('vplot')
	bparm[2] = 22
	bparm[3] = 1 
	uvdata.zap_table('AIPS PL', -1)
	plot_file = plotdir + '/'+project+ '.VPLOT_LOGAMP_UNCAL.PS'
	runvplot2(indata=uvdata, bparm=bparm, docalib = -1, stokes='HALF', solint=0, nplots = nif, do3col = 1, factor = 0.001, xyratio = 1.5)
	plot(uvdata, plot_file, dopng=False)

	# vplot of the raw data
	print("Ploting VPLOT uncalibrated")
	bparm = set_default_bparms('vplot')
	bparm[1:] = [0, -22, 1, 0,0,0,0,0,0]
	for stk_i in ['RR', 'LL']:
		uvdata.zap_table('AIPS PL', -1)
		plot_file = plotdir + '/'+project+ '.VPLOT_UNCAL_'+stk_i+'.PS'
		runvplot2(indata=uvdata, bparm=bparm, docalib = -1, stokes=stk_i, solint=0, nplots = 4, do3col=1, factor = 0.001, xyratio = 1.5)
		plot(uvdata, plot_file, dopng=False)

	# Source elevation vs time
	print("Ploting source elevation")
	uvdata.zap_table('AIPS PL', -1)
	bparm = set_default_bparms('uvplt')
	bparm[1:] = [11, 15, 1, 0, 0, 0, 90, 0]
	runuvplt2(indata=uvdata, bparm=bparm, docalib = -1, nchav = nchan, stokes = 'I', refant = refantn)
	plot_file = plotdir + '/'+project+ '.ELEV.PS'
	plot(uvdata, plot_file, dopng=True)

	# uvplt the u,v coverage
	print("Ploting UV coverage")
	uvdata.zap_table('AIPS PL', -1)
	bparm = set_default_bparms('uvplt')
	bparm[1:] = [6, 7, 2, 0]
	allsources = uvdata.sources
	for s in allsources:
		runuvplt2(indata=uvdata, bparm=bparm, sources=[s], docalib = -1, xinc = 1, nchav = nchan, do3col = 1, stokes = 'I')
	plot_file = plotdir + '/'+project+ '.UVCOV.PS'
	plot(uvdata, plot_file, dopng=True)

	# Antenna info
	an = uvdata.table('AN',0)
	anten_info = []
	for a in an:
		print a['nosta'], a['anname'], a['diameter']
		anten_info.append([a['nosta'], a['anname'], a['diameter']])


	f = open(plotdir + '/'+project+ '.DTSUM.txt', 'rb').readlines()
	for fi in f:
		if fi.split()[0] == 'Timerange:':
			timerange_str = fi[10:]
		if fi.split()[0] == 'Data':
			int_time = float(fi.split()[5])
  
	###  Start weblog  ###

	try: os.mkdir('weblog')
	except: pass
	######  Home page
	wlog = open("weblog/index.html","w")
	weblog_header(wlog, 'Home')
   #------------------------------------------
	wlog.write('<table bgcolor="#eeeeee" border="3px" cellspacing = "0" cellpadding = "4px" style="width:30%">\n')
	wlog.write('<tr><td>Project </td> 	<td> '+project+'</td>\n')
	wlog.write('<tr><td>Date observed </td>  <td> '+uvdata.header['date_obs']+'</td>\n')
	wlog.write('<tr><td>Timerange </td>      <td> '+timerange_str+'</td>\n')
	wlog.write('<tr><td>Band </td>           <td> '+band+'</td>\n')
	wlog.write('<tr><td>Frequency </td>      <td> {0:5.2f} - {1:5.2f} GHz'.format(freq, freq+nif*ifwidth/1e9)+'</td>\n')
	wlog.write('<tr><td>IF </td>             <td> '+str(nif)+'</td>\n')
	wlog.write('<tr><td><abbr title="Number of channels per IF after averaging. The raw data may provide more channels">Channels/IF</abbr></td><td> '+str(nchan)+'</td>\n')
	wlog.write('<tr><td>IF bandwidth </td>   <td> {0:6.0f} MHz</td>\n'.format(ifwidth/1.0e6))
	wlog.write('<tr><td>Total bandwidth </td><td> {0:6.0f} MHz</td>\n'.format(ifwidth*nif/1.0e6))
	wlog.write('<tr><td><abbr title="Integration time after averaging. The raw data may provide higher time resolution">Integration time</abbr></td><td> {0:3.1f} s</td>\n'.format(int_time))
	wlog.write('<tr><td>Polarizations </td>  <td> '+', '.join(str(x) for x in polarizations)+'</td>\n')
	wlog.write('<tr><td>Stokes </td>         <td> '+', '.join(str(x) for x in stokes)+'</td>\n')
	wlog.write('<tr><td>Observer </td>       <td> '+uvdata.header['observer']+'</td>\n')
	wlog.write('<tr><td>Telescope </td>     <td> '+uvdata.header['telescop']+'</td>\n')
	wlog.write('</table><br>\n')
	write_linksize('txt','./'+project+'.notes.txt', 'Notes and observing comments', wlog)
	#------------------------------------------
	weblog_foot(wlog)
	wlog.close()

	###### Observation summary page
	wlog = open("weblog/obs_summary.html","w")
	weblog_header(wlog, 'Observation summary')
	#------------------------------------------
	wlog.write('<h3>Summary:</h3>\n')
	write_linksize('txt','../plots/'+project+'.DTSUM.txt', 'Brief data summary (DTSUM)', wlog)
	write_linksize('txt','../plots/'+project+'.SCAN.txt', 'Scan list (LISTR)', wlog)
	wlog.write('<h3>Sources:</h3>\n')
	wlog.write("Targets = {0}</br>\n".format(", ".join(targets)))
	wlog.write("Phase cals = {0}</br>\n".format(", ".join(phsrefs)))
	wlog.write("Flux density cals = {0}</br>\n".format(", ".join(fluxcals)))
	wlog.write("Bandpass cals = {0}</br>\n".format(", ".join(bpasscals)))
	wlog.write("Point flux cals = {0}</br>\n".format(", ".join(pointcals)))
	wlog.write('<h3>Antennas:</h3>\n')
	wlog.write('<pre>\n')
	for a in anten_info:
		wlog.write("{0:2.0f} {1:3s} {2:3.0f} m\n".format(a[0],a[1],a[2]))
	wlog.write('</pre>\n')
	#elev = glob.glob('./plots/'+project+ '.ELEV_*.png')[0]
	uvcov = glob.glob('./plots/'+project+ '.UVCOV_*.png')
	wlog.write('<h3>Source elevation:</h3>\n')
	wlog.write('<img src=".{0}"><br>\n'.format('./plots/'+project+ '.ELEV_001.png'))
	wlog.write('<h3>UV coverage:</h3>\n')
	for u in uvcov:
		wlog.write('<img src=".{0}"><br>\n'.format(u))
	#------------------------------------------
	weblog_foot(wlog)
	wlog.close()

	###### Plots page
	wlog = open("weblog/plots.html","w")
	weblog_header(wlog, 'Plots')
	#------------------------------------------

	wlog.write('<h3>Observation</h3>\n')
	wlog.write('<table border="1px" cellspacing = "0" cellpadding = "6px" style="width:100%">\n<col width="40%">\n<col width="10%">\n<col width="50%">')
	write_linksize('pdf', '../plots/'+project+'.ELEV.pdf', 'Source elevation vs time', 'Source elevation for reference antenna: {0} ({1}) in this observation.'.format(refant, refantn), '../plots/'+project+'.ELEV_001.png')
	write_linksize('pdf','../plots/'+project+'.UVCOV.pdf', 'UV coverage for all sources', 'One source per page. Each IF shown with a different color. Channels averaged by IF')
	wlog.write('</table>\n')

	wlog.write('<h3>Uncalibrated data</h3>\n')
	wlog.write('<table border="1px" cellspacing = "0" cellpadding = "6px" style="width:100%">\n<col width="40%">\n<col width="10%">\n<col width="50%">')
	write_linksize('pdf', '../plots/'+project+'.POSSM_CALIBRATORS_UNCAL.pdf', 'Uncalibrated amplitude and phase against frequency channel (calibrators)', 'Cross-power spectra per baseline per polarization. Each scan plotted separately. No calibration applied yet. All calibrators plotted.')
	write_linksize('pdf', '../plots/'+project+'.POSSM_TARGET_UNCAL.pdf', 'Uncalibrated amplitude and phase against frequency channel (targets)', 'Cross-power spectra per baseline per polarization. Each scan plotted separately. No calibration applied yet. All calibrators plotted. Probably noise for faint sources.')
	write_linksize('pdf', '../plots/'+project+'.VPLOT_LOGAMP_UNCAL.pdf', 'Uncalibrated amplitude against time (log scale)', 'Logarithm of visibility amplitude for RR (red) and LL (blue) polarization.')
	write_linksize('pdf', '../plots/'+project+'.VPLOT_UNCAL_RR.pdf', 'Uncalibrated amplitude and phase against time (RR)', 'Visibility amplitude and phase for RR polarization.')
	write_linksize('pdf', '../plots/'+project+'.VPLOT_UNCAL_LL.pdf', 'Uncalibrated amplitude and phase against time (LL)', 'Visibility amplitude and phase for LL polarization.')
	wlog.write('</table>\n')


	wlog.write('<h3>Dynamic spectrum</h3>\n')
	sources_spplot = numpy.unique(numpy.hstack([pointcals, fluxcals, bpasscals, phsrefs]))
	#sources_spplot = numpy.unique(numpy.hstack([pointcals, fluxcals, bpasscals]))
	for source_spplot in sources_spplot:
		wlog.write(source_spplot+'<br>\n')
		wlog.write('<table border="1px" cellspacing = "0" cellpadding = "6px" style="width:100%">\n<col width="40%">\n<col width="10%">\n<col width="50%">')
		write_linksize('pdf', '../plots/'+project+'.SPPLOT.PHS.UNCAL_'+source_spplot+'.pdf', 'Uncalibrated dynamic spectrum (phase)', 'Visibility phase as a function of time (increases towards the top) and frequency (increases towards the right. One baseline per page.')
		write_linksize('pdf', '../plots/'+project+'.SPPLOT.PHS.UNCAL_'+source_spplot+'.combined.pdf', 'Uncalibrated combined dynamic spectrum (phase)','Combined phase dynamic spectrum with one baseline per subplot.', '../plots/'+project+'.SPPLOT.PHS.UNCAL_'+source_spplot+'.combined.png')
		write_linksize('pdf', '../plots/'+project+'.SPPLOT.AMP.UNCAL_'+source_spplot+'.pdf', 'Uncalibrated dynamic spectrum (amplitude)', 'Visibility amplitude as a function of time (increases towards the top) and frequency (increases towards the right. One baseline per page.')
		write_linksize('pdf', '../plots/'+project+'.SPPLOT.AMP.UNCAL_'+source_spplot+'.combined.pdf', 'Uncalibrated combined dynamic spectrum (amplitude)', 'Combined amplitude dynamic spectrum with one baseline per subplot.', '../plots/'+project+'.SPPLOT.AMP.UNCAL_'+source_spplot+'.combined.png')
		wlog.write('</table><br>\n')

	#wlog.write('<h3>Calibration solution tables</h3><br>\n')
	#wlog.write('<h3>Calibrated data</h3><br>\n')
	#wlog.write('<h3>Images</h3><br>\n')
	#------------------------------------------
	weblog_foot(wlog)
	wlog.close()

	###### FITS files
	wlog = open("weblog/files.html","w")
	weblog_header(wlog, 'FITS files')
	#------------------------------------------
	wlog.write('<table border="1px" cellspacing = "0" cellpadding = "6px" style="width:100%">\n<col width="40%">\n<col width="10%">\n<col width="50%">')
	write_linksize('fits', '../fits/'+project+'.RAW.AVG.UNCAL.FITS', 'Raw averaged uncalibrated data', 'Raw data concatenated, averaged, and TB sorted.')
	wlog.write('</table><br>\n')

   #------------------------------------------
	weblog_foot(wlog)
	wlog.close()






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

	if uvdata.header.crpix[2] != 1:
		print "Error: FREQ reference pixel is not channel 1, quitting."


	# Run CLCOR if needed to correct the Cambridge axis offset problem in data
	# before mid-2015, then splat data to apply correction.
	antab = uvdata.table('AN',1)
	for row in antab:
		if 'Cm' in row.anname:
			if row.staxof == 0:
				print "Doing Cambridge axis offset correction."
				runclcor(uvdata, row.nosta, indisk)

	(sourcelist, fringlist, caliblist, bpasslist, nchan, nif, stokes, plotstokes, polarizations, bchan, echan, freq, band) = setup(uvdata, dosort, 1)

	(refant, refantn, refantlist) = set_refant(uvdata, refant)

	cellsize = findmaxb(uvdata)

	snver = 0
	clver = 0
	flagver = get_tab(uvdata, 'FG')

	if dodiagnostic3 == 1 :
		print "------------------------------------------------------------"
		print "          STAGE 3 diagnostics for e-MERLIN data"
		print "            dodiagnostic3 = 1 (part of docalib)             " 
		print "------------------------------------------------------------"

		timer = 0,0,0,0,0,0,0,0
		aparm = 0,0,0,0,0,0,0,0,1,0
		bparm = 0,0,0,0,0,0,0,0,0,0
		for i in range(len(fringlist)) :
			sourceID = fringlist[i],
			print 'Running POSSM on', fringlist[i]
			runquickpossm(uvdata, sourceID)

			plotfile = 'POSSM_NOCAL_' + str(fringlist[i]) + '_' + strftime("%y%m%d_%H%M") + '.ps'
			outfile = os.path.join(plotdir, plotfile)
			print 'Plotting POSSM diagnostics:' + outfile
			runlwpla(uvdata, outfile, 0)


		print 'Running VPLOT on', fringlist
		timer = 0,0,0,0,0,0,0,0
		runvplot(uvdata, fringlist, 'HALF', timer, refantlist, refantlist, 1, nif, bchan, echan, -1, 0, 0, -1, 0, aparm, bparm, refantn, dotv, 9)
		outfile = os.path.join(plotdir, 'VPLOT_NOCAL_' + strftime("%y%m%d_%H%M") + '.ps')
		print 'Plotting VPLOT diagnostics:' + outfile
		runlwpla(uvdata, outfile, 0)




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
		aparm = 3, 0, 0, 0, 0, 0, 2.5, 0, 0, 0
		dparm = 2, 450, 0, 0, 0, 0, 0, 1, 1, 0 
		timer = 0,0,0,0,0,0,0,0
		solint = 10
		print 'Running initial fringe fit (SN#1)'
		newlist = list(set(fluxcals+bpasscals+pointcals+phsrefs))
		runfring(uvdata, newlist, timer, -1, 0, 0, -1, -1, refantn, refantlist, solint, aparm, dparm, 0, fring_snr, bchan, echan)

		# Edit solutions
		if dosnedit == 1:
			snver = get_tab (uvdata, 'SN')
			wuvdata = WizAIPSUVData ('ALL','DBCON', indisk, 1)
			try: edit_table ('ALL','DBCON','delay')
			except: pass

		# Plot Delay solutions (SNPLT)
		optype = 'DELA'
		snver=0
		snver = get_tab(uvdata, 'SN')
		print 'SN version' + format(snver) + 'contains fringe fit solutions for:', fringlist
		runsnplt(uvdata, 'SN', 1, fringlist, 9, optype, dotv)
		plotfile = 'SNPLT_POSTFRING_' + strftime("%y%m%d_%H%M") + '.ps'
		outfile = os.path.join(plotdir, plotfile)
		print 'Plotting SNPLT diagnostics:' + outfile
		runlwpla(uvdata, outfile, 0)

		# Apply solutions with CLCAL
		runclcal(uvdata, sourcelist, fringlist, 'CALI', 'AMBG', snver, snver, 0, 0, refantn)

		for i in range(len(fringlist)) :
			sourceID = fringlist[i],
			bparm = 0,0
			aparm = 0,0,0,0,0,0,0,0,1,0
			print 'Running POSSM on', fringlist[i]
			runpossm(uvdata, sourceID, timer, refantlist, refantlist, aparm, bparm, 0, nif, 0, 0, 100, 0, flagver, plotstokes, -1, 0, 'A&P', 30, 9, -1, 1)
			plotfile = 'POSSM_POSTFRING_' + str(fringlist[i]) + '_' + strftime("%y%m%d_%H%M") + '.ps'
			outfile = os.path.join(plotdir, plotfile)
			print 'Plotting POSSM diagnostics:' + outfile
			runlwpla(uvdata, outfile, 0)
			print 'CHECK YOUR INITIAL FRINGE FIT SOLUTIONS. You may need to edit SN table and/or re-run fring.'

	else :
		print "****************************************************************"
		print "*** Fringe fitting was not selected................ skipping.***"
		print "****************************************************************"




	if dosetflux == 1 :
		print "------------------------------------------------------------"
		print "            Setting flux density scale for e-MERLIN data    "
		print "                 dosetflux = 1 (part of docalib)            " 
		print "------------------------------------------------------------"

		print "Reset all fluxes to zero"
		runsetjy(uvdata, '', 0, 0, 0, 'REJY')


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
			index=0
			reffreq = uvdata.header.crval[2]
			if uvdata.header.crpix[2] > 1:
				reffreq = reffreq - (uvdata.crpix[2] * ch_width)
			for ifstep in fqtab[0].if_freq :		# gives freq at middle of each IF
				setfreq = (reffreq + ifstep + (fqtab[0].total_bandwidth[index]/2) - (fqtab[0].ch_width[index]/2) ) / 1000000
				subbandfreq.append(setfreq)
				# Then calculate the appropriate flux and overwite fluxcalflux
				fluxcalflux.append(dfluxpy(setfreq,uvdata))

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

	if dofluxcal == 1 :
		print "------------------------------------------------------------"
		print "            Running CALIB on fluxcal     "
		print "                 dofluxcal = 1 (part of docalib)            " 
		print "------------------------------------------------------------"

		# Load appropriate 3C286 model, if not there already
		model = AIPSImage('3C286','MODEL',indisk,1)
		if not model.exists():
			if band == 'L':
				modfil = os.path.join(os.getcwd(), '3C286_CWT_LBAND.ICL001')
				try: runimlod(modfil, indisk, 1)
				except: "Error: can't load L-band 3C286 model!"
			elif band == 'C':
				modfil = os.path.join(os.getcwd(), '3C286_C.ICL001')
				try: runimlod(modfil, indisk, 1)
				except: "Error: can't load C-band 3C286 model!"
			else:
				print "Error: I don't know how to calibrate data in this frequency range!"
				sys.exit()
			model = AIPSImage('3C286','MODEL',indisk,1)
			print "Loading 3C286 model."
		else:
			print "WARNING: I found a 3C286 model, assuming it is correct."

		# Run P-only CALIB on 3C286
		aparm = 3, 0, 0, 0, 0, 3, 5, 0, 0, 0
		cparm = 10, 0, 0, 0, 0, 0, 0
		timer = 0,0,0,0,0,0,0,0
		solint=0.5
		uvrang = 0, 0
		antwt = 0, 0
		clver = get_tab(uvdata, 'CL')
		snver = get_tab(uvdata, 'SN') + 1
		ncomp = 1000,0
		weightit = 1
		calsour = '1331+305',''
		runselfcalib(uvdata, model, calsour, timer, uvrang, 100, clver, 0, -1, -1, 'DFT', refantn, solint, aparm, -1, 'L1', 'P', 0, 0, cparm, snver, antwt, 1, 1, ncomp)

		# Run P-only CALIB on other cals
		model = False
		ncomp = 0,0
		print 'CALIBLIST = ', caliblist
		runselfcalib(uvdata, model, caliblist, timer, uvrang, 100, clver, 0, -1, -1, 'DFT', refantn, solint, aparm, -1, 'L1', 'P', 0, 0, cparm, snver, antwt, weightit, 0, ncomp)

		# Edit solutions
		if dosnedit == 1:
			try: edit_table ('ALL','DBCON','phase')
			except: pass
			snver = get_tab(uvdata, 'SN')

		# Run SNPLT
		# Run CLCAL
		snver = get_tab(uvdata, 'SN')
		print "Running clcal on snver", snver, "with", caliblist
		runclcal(uvdata, '', '', 'CALI', '', snver, snver, 0, 0, refantn)

		# Run POSSM
		# Run first BPASS
		bpassprm = 0, 0, 0, 0, 1, 0, 0, 0, 1, 3, 1 
		print 'Running BPASS on', uvdata
		bpasssrc = bpasscals[0],
		runbpass(uvdata, bpasssrc, refantn, bpassprm, 'L1', -1, 100, 0)

		# Run A&P CALIB on 3C286
		aparm = 3, 0, 0, 0, 0, 3, 5, 0, 0, 0
		cparm = 10, 0, 0, 0, 0, 0, 0
		timer = 0,0,0,0,0,0,0,0
		solint=2
		uvrang = 0, 0
		antwt = 0, 0
		clver = get_tab(uvdata, 'CL')
		snver = get_tab(uvdata, 'SN') + 1
		calsour = '1331+305',''
		ncomp = 1000,0
		model = AIPSImage('3C286','MODEL',indisk,1)
		runselfcalib(uvdata, model, calsour, timer, uvrang, 100, 0, 0, 4, 1, 'DFT', refantn, solint, aparm, -1, 'L1', 'A&P', 10, 10, cparm, snver, antwt, 1, 1, ncomp)

		# Run A&P CALIB on other cals
		model = False
		ncomp = 0,0
		runselfcalib(uvdata, model, caliblist, timer, uvrang, 100, 0, 0, 4, 1, 'DFT', refantn, solint, aparm, -1, 'L1', 'A&P', 10, 10, cparm, snver, antwt, 1, 0, ncomp)

		# Edit solutions
		if dosnedit == 1:
			try: edit_table ('ALL','DBCON','phase')
			except: pass
			snver = get_tab (uvdata, 'SN')

		# Run GETJY
		for i in range(len(caliblist)) :
			sourceID = caliblist[i],
			print "Running GETJY on caliblist", sourceID[0]
			rungetjy(uvdata, sourceID, fluxcals, 1, nif, snver)

		print "Attempted to set FLUX density scale using 3C286 model - CHECK YOUR RESULTS in your SU TABLE"

	else :
		print "You opted not to derive flux scale... assuming you want to do this manually"



	if dobandpass == 1 :
		print "------------------------------------------------------------"
		print "            Doing a bandpass of your e-MERLIN data          "
		print "                 dobandpass = 1 (part of docalib)           "
		print "------------------------------------------------------------"   
		# Run SOUSP on point-cal
		print 'Running SOUSP.'
		bpasssrc = bpasscals[0],
		bpspecindx = runsousp(uvdata, bpasssrc)

		# Run SOUSP on phase cal(s)
		for i in range(len(phsrefs)):
			specindx = runsousp(uvdata, bpasssrc)

		# Remove BP#1 and top SN table
		uvdata.zap_table('BP', -1)
		uvdata.zap_table('SN', snver)

		# Run BPASS
		soltyp = 'L1'
		solint = 0
		docal = 100
		bpassprm = 0, 0, 0, 0, 1, 0, 0, 0, 1, 3, 1 
		print 'Running BPASS on', uvdata
		runbpass(uvdata, bpasssrc, refantn, bpassprm, soltyp, solint, docal, bpspecindx)


	if dotherest == 1 :
		step = 0

		# Run A&P CALIB on 3C286
		aparm = 3, 0, 0, 0, 0, 3, 5, 0, 0, 0
		cparm = 10, 0, 0, 0, 0, 0, 0
		timer = 0,0,0,0,0,0,0,0
		solint=2
		uvrang = 0, 0
		antwt = 0, 0
		clver = get_tab(uvdata, 'CL')
		snver = get_tab(uvdata, 'SN') + 1
		calsour = '1331+305',''
		ncomp = 1000,0
		model = AIPSImage('3C286','MODEL',indisk,1)
		runselfcalib(uvdata, model, calsour, timer, uvrang, 100, 0, 0, 4, 1, 'DFT', refantn, solint, aparm, -1, 'L1', 'A&P', 10, 10, cparm, snver, antwt, 1, 1, ncomp)

		# Run A&P CALIB on other cals
		model = False
		ncomp = 0,0
		runselfcalib(uvdata, model, caliblist, timer, uvrang, 100, 0, 0, 4, 1, 'DFT', refantn, solint, aparm, -1, 'L1', 'A&P', 10, 10, cparm, snver, antwt, 1, 0, ncomp)

		# Edit solutions
		if dosnedit == 1:
			try: edit_table('ALL','DBCON','phase')
			except: pass
			snver = get_tab (uvdata, 'SN')

		# Run CLCAL
		snver = get_tab(uvdata, 'SN')
		print "Running final clcal on snver", snver, "with all sources"
		runclcal(uvdata, '', '', 'CALI', '', snver, snver, 0, 0, refantn)

		# Run POSSM using solint ~30mins on all cals.
		timer = 0,0,0,0,0,0,0,0
		aparm = 0,0,0,0,0,0,0,0,1,0
		bparm = 0,0,0,0,0,0,0,0,0,0
		for i in range(len(caliblist)) :
			sourceID = caliblist[i],
			print 'Running POSSM on', caliblist[i]
			runpossm(uvdata, sourceID, timer, refantlist, refantlist, aparm, bparm, 0, nif, 0, 0, 1, 0, 0, plotstokes, -1, 0, 'A&P', 30, 9, -1, 1)

		plotfile = 'POSSM_CALIB_' + strftime("%y%m%d_%H%M") + '.ps'
		outfile = os.path.join(plotdir, plotfile)
		print 'Plotting POSSM diagnostics:' + outfile
		runlwpla(uvdata, outfile, 0)

		# Plot SNPLT diagnostics
		clver = get_tab(uvdata, 'CL')
		runsnplt(uvdata, 'CL', clver, caliblist, nplots=9, optype='PHAS', dotv=dotv)
		plotfile = 'SNPLT_CALIB_' + strftime("%y%m%d_%H%M") + '.ps'
		outfile = os.path.join(plotdir, plotfile)
		print 'Plotting SNPLT diagnostics:' + outfile
		runlwpla(uvdata, outfile, 0)


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
				nmaps = 1
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
				ccneg = 0
				cctab = uvdata2.table('CC',1)
				for row in cctab:
					if row.flux < 0:
						break
					ccneg = ccneg + 1
				ccneg = ccneg,0
				print "Attempting self-calibration with", source, "clver", clver, " and selfcal mode:", mode
				runselfcalib(uvdata, uvdata2, sourceID, timer, uvrang, 100, clver, 0, 1, 0, 'DFT', refantn, solint, aparm, 1, 'L1', mode, 10, 10, cparm, 0, antwt, 1, nmaps, ncomp)

				# Edit solutions
				if dosnedit == 1:
					try: edit_table ('ALL','DBCON','phase',inver=0)
					except: pass

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
				runlwpla(uvdata, outfile, 0)


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
	runlwpla(uvdata, outfile, 0)



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

if (doimage == 1) and (dowide != 1) :

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
			print "Error: ALL.DBCON.1 does not exist! Don't know what to map."
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
		runlwpla(image,outfile, 0)

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
		runlwpla(image,outfile, 0)


	print '************************************************************'
	print "Finished. Now admire the awesomeness of your e-MERLIN maps."
	print "Your diagnostic plots and pipelined images can be found in:"
	print plotdir
	print '------------------------------------------------------------'
	print 'Note that it is highly unlikely that your images will be of'
	print 'publication quality. It is recommended that you refine the '
	print 'flagging and calibration by hand and re-image your target(s)'
	print '************************************************************'




#################################################################################################################
# Fast wide field imaging routines
#################################################################################################################

if dowide == 1:

	# split the target sources off the multisource file
	if dosplit == 1:

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

		# Spilt off each target from the main dbconned uvdata and apply calibration
		for i in range(len(widetargets)):
			runsplit(uvdata,widetargets[i])


	pca = AIPSCat(indisk)
	cont = 0
	for i in range(len(widetargets)):
		widetarget=widetargets[i]
		if len(widetarget)>12:		# check that the target name wasn't too long
			widetarget=widetarget[0:12]
		for fitsfil in pca[indisk]:
			if (widetarget in fitsfil.name) and (fitsfil.klass == 'SPLIT'):
				cont = 1
			
	if cont == 1:
		pass
	else:
		print "###############################################################################"
		print "Error: Can't find one or more widetarget SPLIT files. Don't know what to image."
		print "###############################################################################"
		sys.exit()



	#chop up the field
	if dochop == 1:

		pca = AIPSCat(indisk)

		#execute the chessboard chopper with a maximum radius in arcminutes
		for i in range(len(widetargets)):
			widetarget=widetargets[i]
			if len(widetarget)>12:		# check that the target name wasn't too long
				widetarget=widetarget[0:12]
			for fitsfil in pca[indisk]:
				if (widetarget in fitsfil.name) and (fitsfil.klass == 'SPLIT'):
					uvdata = AIPSUVData(fitsfil.name,fitsfil.klass,indisk,fitsfil.seq)
					chessboard64wide(uvdata,fovradius)

	if dopeel == 1:
		pca = AIPSCat(indisk)

		#run the 'peeler'
		for i in range(len(widetargets)):
			widetarget=widetargets[i]
			if len(widetarget)>12:		# check that the target name wasn't too long
				widetarget=widetarget[0:12]
			for fitsfil in pca[indisk]:
				if (widetarget in fitsfil.name) and ('UV' in fitsfil.klass):
					uvdata = AIPSUVData(fitsfil.name,fitsfil.klass,indisk,fitsfil.seq)
					print 'Renaming ' + fitsfil.name + '.' + fitsfil.klass + ' to ' + 'OR'+ fitsfil.klass[2:]
					uvdata.rename(fitsfil.name,'OR'+ fitsfil.klass[2:],fitsfil.seq)
					print 'Executing the peeler on ' + 'OR'+ fitsfil.klass[2:]				
					cbang(fitsfil.name,'OR'+ fitsfil.klass[2:],indisk,fitsfil.seq,cellsize=0.045,imsize=128,niter=100,zap=True,infile='')

	#clean the facets
	if domap == 1:
		# now image all the uvfiles, base the imaging size larger than each square

		pca = AIPSCat(indisk)

		for i in range(len(widetargets)):
			widetarget=widetargets[i]
			#check that the target name wasn't too long
			if len(widetarget)>12:
				widetarget=widetarget[0:12]
			for fitsfil in pca[indisk]:
				if (widetarget in fitsfil.name) and ('ICL' in fitsfil.klass):
					print "Warning: found existing image data - this *WILL* cause problems!"
					inputinput("Are you sure you want to continue? ([y]/n): ")
				if (widetarget in fitsfil.name) and ('UV' in fitsfil.klass):
					uvdata = AIPSUVData(fitsfil.name,fitsfil.klass,indisk,fitsfil.seq)
					widebandmaps(uvdata,indisk,fovradius,widetarget,noiter) #temporary change here to force frequency dependent imaging nwrigley 20141105
					#widemaps(uvdata,indisk,fovradius,widetarget,noiter) #temporary change here to force standard  imaging nwrigley 20150219

		#now use APCLN  added nwrigley 20141106
		pca = AIPSCat(indisk)

		for i in range(len(widetargets)):
			widetarget=widetargets[i]
			if len(widetarget)>12:			# check that the target name wasn't too long
				widetarget=widetarget[0:12]		
		for fitsfil in pca[indisk]:					
			if (widetarget in fitsfil.name) and ('IIM' in fitsfil.klass):
				imdata = AIPSImage(fitsfil.name,fitsfil.klass,indisk,fitsfil.seq)
				apclean(imdata,widetargets[i]) #widebandmaps
				#apclean2(imdata,widetargets[i]) #widemaps





	# trim and flatten the images
	if doflatn == 1:

		pca = AIPSCat(indisk)

		for i in range(len(widetargets)):
			widetarget=widetargets[i]
			if len(widetarget)>12:			# check that the target name wasn't too long
				widetarget=widetarget[0:12]
			for fitsfil in pca[indisk]:
				if (widetarget in fitsfil.name) and ('ICL' in fitsfil.klass):
					imdata = AIPSImage(fitsfil.name,fitsfil.klass,indisk,fitsfil.seq)
					flat(imdata,fovradius,widetargets[i])
			cellpitch = (imdata.header.bmin * 3600) / 3

			imgsize = int((fovradius*2*60/(cellpitch*16)))
			thing = imgsize
			count=2
			while thing > 1.0 :
				thing = thing / 2
				count = count + 1

			imgsize = 2**count

			cellpitch = abs(imdata.header.cdelt[0] * 3600)

			#now use flatn
			print 'Flattening facets'
			flatn=AIPSTask('FLATN')
			suboutname='SUB' + widetarget[-4:-1]+widetarget[-1]
			flatn.inname=suboutname
			flatn.inseq=0
			flatn.outdisk=indisk
			flatn.outname=widetarget
			flatn.outclass='FLATN'
			flatn.imsize[1]=256+(fovradius*60*2/cellpitch)
			flatn.imsize[2]=256+(fovradius*60*2/cellpitch)
			flatn.nfield=256
			flatn.inseq=0
			flatn.nmaps=1
			flatn.go()



	#beam correction
	if dopbcor == 1:

		pca = AIPSCat(indisk)

		for i in range(len(widetargets)):
			widetarget=widetargets[i]
			if len(widetarget)>12:
				widetarget=widetarget[0:12]
			for fitsfil in pca[indisk]:
				if (widetarget in fitsfil.name) and ('FLATN' in fitsfil.klass):
					imdata = AIPSImage(fitsfil.name,fitsfil.klass,indisk,fitsfil.seq)
					print "Initialising primary beam correction.."
					pbcorr(imdata)




	#Catalogue
	if dosad == 1:
		for i in range(len(widetargets)):
			widetarget=widetargets[i]
			#check that the target name wasn't too long
			if len(widetarget)>12:
				widetarget=widetarget[0:12]
			for fitsfil in pca[indisk]:
				if (widetarget in fitsfil.name) and ('PBCOR' in fitsfil.klass):
					imdata = AIPSImage(fitsfil.name,fitsfil.klass,indisk,fitsfil.seq)
					print "Searching for sources.."
					SAD(imdata)
			#create a contour plot of the resampled image
			#image = AIPSImage(widetarget,'REGRD',uvdata[2],0)
			#runkntr(image,3,AIPS.userno)
			#plotfile = widetarget + '_' + strftime("%y%m%d_%H%M") + '.ps'
			#outfile = os.path.join(plotdir, plotfile)
			#runlwpla(image,outfile)








print "Your AIPS catalog now looks like this:"
print AIPSCat(indisk)

if mixed:
	print "WARNING: your data were mixed mode (including narrow spectral IFs)"
	print "WARNING: the continuum data have been split and reassembled"
	print "WARNING: you can find your line data in AIPS number", AIPS.userno + 1
else:
	print "WARNING: if your data were observed in mixed mode (mixed IF widths),"
	print "WARNING: then you WILL need to seperate your data before calibrating."
	print "WARNING: If you are unsure, contact e-MERLIN staff for advice."



# Finished whatever was requested, print the time taken.
print "Completed the requested procedures.  Check your results carefully!"
print "Time taken (hh:mm:ss):", time2hms(time.time()-ti)
print 'Finished on %s' % strftime("%d-%b-%y"), 'at:', strftime("%H:%M:%S", localtime()),'\n'


