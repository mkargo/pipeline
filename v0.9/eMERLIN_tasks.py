# Define the AIPS tasks to use in the e-MERLIN pipeline.
# mkargo 2011
# Last updated mkargo 20160504, amended 20160801 added SAD routine output nw

import os, sys, math, numpy, re
from AIPS import AIPS, AIPSDisk
from AIPSTask import AIPSTask, AIPSList
from AIPSData import AIPSUVData, AIPSImage, AIPSCat
from Wizardry.AIPSData import AIPSUVData as WizAIPSUVData
import astLib;from astLib import astCoords
from collections import deque
import numpy as np


def runfitld(datain, indisk, thisdir):
	fitld = AIPSTask('FITLD')
	fitld.datain = datain		# "datain" is correct here (use "indata" in most other places)
	fitld.outdisk = indisk
	fitld.digicor = -1
	fitld.douvcomp = -1
	fitld.clint = 8/60
	fitld.bif = 1
	fitld.eif = 0
	fitld.outdata = AIPSUVData('TMP','UVDATA',indisk,thisdir)
	fitld.go()

def runfitld2(outdata, infile, nfits=0):
	"""Must set outdata, infile"""
	#assert (outdata != None, infile != None), "set indata, infile in runfitld" 
	fitld = AIPSTask('fitld')
	fitld.outdata = outdata
	fitld.ncount = nfits
	fitld.douvcomp = 1
	fitld.doconcat = 1
	fitld.clint = .25
	try:
		fitld.datain = infile
	except AttributeError:
		fitld.infile = infile
	# trick to ignore data with '0' time
	fitld.timerang[1:] = [0, 0, 0, 0, 0, 0, 0, 0]
	#print "fits file = " + fitld.infile
	fitld()

def runimlod(datain, indisk, thisdir):
	imlod = AIPSTask('IMLOD')
	imlod.outdata = AIPSUVData('3C286','MODEL',indisk,thisdir)
	imlod.datain = datain
	imlod.go()

def runmsort(indata):
	print 'Running MSORT.'
	msort = AIPSTask('MSORT')
	msort.indata = indata
	msort.outdata = indata
	msort.sort = 'TB'
	msort.go()

def runindxr(indata):
	print 'Running INDXR.'
	indxr = AIPSTask('INDXR')
	indxr.cparm[3] = 8/60
	indxr.indata = indata
	indxr.go()

def runfring(indata, calsour, timer, docalib, gainuse, flagver, doband, bpver, refant, refantlist, solint, aparm, dparm, snver, snr, bchan, echan):
	print 'Running FRING on ', format(calsour)
	fring = AIPSTask('FRING')
	fring.indata = indata
	fring.calsour[1:] = calsour
	fring.timer[1:] = timer
	fring.docalib = docalib
	fring.gainuse = gainuse
	fring.flagver = flagver
	fring.doband = doband
	fring.bpver = bpver
	fring.refant = refant
	fring.search[1:] = refantlist
	fring.solint = solint
	fring.cmethod = 'DFT'
	fring.aparm[1:] = aparm
	fring.aparm[7] = snr
	fring.dparm[1:] = dparm
	fring.snver = snver
	fring.bchan = bchan
	fring.echan = echan
	fring.go()

def runsnplt(indata, inext, invers, sources, nplots, optype, dotv):
	snplt = AIPSTask('SNPLT')
	snplt.indata = indata
	snplt.inext= inext
	snplt.invers = invers
	snplt.sources[1:] = sources
	snplt.nplots = nplots
	snplt.optype = optype
	snplt.dotv = dotv
	snplt.go()

def runlwpla(indata, outfile, qa):
	plver=0
	outfile = re.sub(r"\s+", '-', outfile)
	for tab in indata.tables :
		if 'PL' in tab[1] :
			plver = plver + 1
	if os.path.exists(outfile) :
		os.remove(outfile)
	lwpla = AIPSTask('LWPLA')
	lwpla.indata = indata
	lwpla.plver = 1
	lwpla.inver = plver
	lwpla.outfile = outfile
	lwpla.msgkill = -10
	if qa == 1:
		lwpla.dparm[5] = 2
		lwpla.docolor = 1
		lwpla.plcolors[10][1:] = 1,1,1
	lwpla.go()
	# Remove PL tables
	indata.zap_table('PL', -1)

def runclcal(indata, sources, calsour, opcode, interpol, snver, inver, gainver, gainuse, refant):
	clcal = AIPSTask('CLCAL')
	clcal.indata = indata
	clcal.sources[1:] = sources
	clcal.calsour[1:] = calsour
	clcal.opcode = opcode
	clcal.interpol = interpol
	clcal.snver = snver
	clcal.inver = inver
	clcal.gainver = gainver
	clcal.gainuse = gainuse
	clcal.refant = refant
	clcal.go()


def runpossm(indata, sources, timer, anten, basel, aparm, bparm, bif, eif, bchan, echan, docalib, gainuse, flagver, stokes, doband, bpver, codetype, solint, nplots, dotv, freqid):
	possm = AIPSTask('POSSM')
	possm.indata = indata
	possm.sources[1:] = sources
	possm.timer[1:] = timer
	possm.aparm[1:] = aparm
	possm.bparm[1:] = bparm
	i = 1
	for ant in anten:
		possm.antennas[i] = int(ant)
		i = i + 1
	i = 1
	for ant in basel:
		possm.baseline[i] = int(ant)
		i = i + 1
	possm.bif = bif
	possm.eif = eif
	possm.bchan = bchan
	possm.echan = echan
	possm.docalib = docalib
	possm.gainuse = gainuse
	possm.flagver = flagver
	possm.stokes = stokes
	possm.doband = doband
	possm.bpver = bpver
	possm.codetype = codetype
	possm.solint = solint
	possm.nplots = nplots
	possm.dotv = dotv
#	possm.freqid = 1
	possm.go()

def runquickpossm(indata, source):
	possm = AIPSTask('POSSM')
	possm.indata = indata
	possm.source[1:] = source
	possm.aparm[9] = 1
	possm.nplots = 6
	possm.solint = 30 # plot every 30min
	possm.stokes = 'HALF'
	possm.flagver = 0
	possm.go()


def runbpass(indata, calsour, refant, bpassprm, soltyp, solint, docal, specindx):
#	bpass = AIPSTask('BPASS', version='OLD')
	bpass = AIPSTask('BPASS')
#	bpass.version='OLD'	# Only if on the JODRELL system!
#	bpass.uvrange[1:] = 0,0
	bpass.indata = indata
	bpass.calsour[1:] = calsour
	bpass.refant = refant
	bpass.bpassprm[1:] = bpassprm
#	bpass.bif = bif
#	bpass.eif = eif
	bpass.soltype = soltyp
	bpass.solint = solint
	bpass.docal = docal
	bpass.specindx = specindx
	bpass.go()

def runsetjy(indata, sources, bif, eif, zerosp, optype):
	setjy = AIPSTask('SETJY')
	setjy.indata = indata
	setjy.sources[1:] = sources
	setjy.bif = bif
	setjy.eif = eif
	setjy.zerosp[1:] = zerosp, 0, 0, 0
	setjy.optype = optype
	setjy.go()

def runcalib(indata, calsour, timer, uvrang, docalib, gainuse, flagver, doband, bpver, cmethod, refant, solint, aparm, doflag, soltype, solmode, minamper, minphser, cparm, snver, antwt, weightit):
	calib = AIPSTask('CALIB')
	calib.indata = indata
	calib.calsour[1:] = calsour
	calib.timer[1:] = timer
	calib.uvrange[1:] = uvrang
	calib.docalib = docalib
	calib.gainuse = gainuse
	calib.flagver = flagver
	calib.doband = doband
	calib.bpver = bpver
	calib.cmethod = cmethod
	calib.refant = refant
	calib.solint = solint
	calib.cparm[1:] = cparm
	if(solmode == 'P'):
		calib.aparm[1] = 3
	else:
		calib.aparm[1] = 4
	calib.aparm[6] = 3
	calib.aparm[7] = 5
	calib.cparm[1] = 10
	calib.doflag = doflag
	calib.soltype = soltype
	calib.solmode = solmode
	calib.minamper = minamper
	calib.minphser = minphser
	calib.snver = snver
	calib.antwt[1:] = antwt
	calib.weightit = weightit
	calib.go()



def runselfcalib(indata, model, calsour, timer, uvrang, docalib, gainuse, flagver, doband, bpver, cmethod, refant, solint, aparm, doflag, soltype, solmode, minamper, minphser, cparm, snver, antwt, weightit, nmaps, ncomp):
	# CALIB procedure for self-calibration, or for doing flux calibration using the 3C286 model
	calib = AIPSTask('CALIB')
	calib.indata = indata
	try:
		model.exists()
		calib.in2name = model.name
		calib.in2class = model.klass
		calib.in2seq = model.seq
		calib.in2disk = model.disk
		calib.invers = 1
		#calib.in2data = in2data
	except:
		print "no model!"
	calib.calsour[1:] = calsour
	calib.timer[1:] = timer[1:]
	calib.uvrange[1:] = uvrang[1:]
	calib.docalib = docalib
	calib.gainuse = gainuse
	calib.flagver = flagver
	calib.doband = doband
	calib.bpver = bpver
	calib.cmethod = cmethod
	calib.refant = refant
	calib.solint = solint
#	calib.aparm[1:] = aparm[1:]
	if(solmode == 'P'):
		calib.aparm[1] = 3
	else:
		calib.aparm[1] = 4
	calib.aparm[6] = 3
	calib.aparm[7] = 5
	calib.cparm[1] = 10
	calib.doflag = doflag
	calib.soltype = soltype
	calib.solmode = solmode
	calib.minamper = minamper
	calib.minphser = minphser
	calib.cparm[1:] = cparm[1:]
	calib.snver = snver
	calib.antwt[1:] = antwt[1:]
	calib.weightit = weightit
	calib.nmaps = nmaps
	calib.ncomp[1:] = ncomp[1:]
	calib.go()



def runuvplt(indata, sources, stokes, timer, anten, basel, bif, eif, docalib, gainuse, flagver, doband, bpver, aparm, bparm, doweight, refant, dotv):
	uvplt = AIPSTask('UVPLT')
	uvplt.indata = indata
	uvplt.sources = sources
	uvplt.stokes = stokes
	uvplt.timerang = timer
	uvplt.antennas[1:] = anten
	uvplt.baseline[1:] = basel
	uvplt.bif = bif
	uvplt.eif = eif
	uvplt.docalib = docalib
	uvplt.gainuse = gainuse
	uvplt.flagver = flagver
	uvplt.doband = doband
	uvplt.bpver = bpver
	uvplt.aparm[1:] = aparm[1:]
	uvplt.bparm[1:] = bparm[1:]
	uvplt.doweight = doweight
	uvplt.refant = refant
	uvplt.dotv = dotv
	uvplt.go()



def runvplot(indata, sources, stokes, timer, anten, basel, bif, eif, bchan, echan, docalib, gainuse, flagver, doband, bpver, aparm, bparm, refant, dotv, nplots):
	vplot = AIPSTask('VPLOT')
	vplot.indata = indata
	vplot.sources[1:] = sources
	vplot.stokes = stokes
	vplot.timerang[1:] = timer
	vplot.antennas[1:] = anten
	vplot.baseline[1:] = basel
	vplot.bchan = bchan
	vplot.echan = echan
	vplot.avgchan = 1
	vplot.bif = bif
	vplot.eif = eif
	vplot.docalib = docalib
	vplot.gainuse = gainuse
	vplot.flagver = flagver
	vplot.doband = doband
	vplot.bpver = bpver
	vplot.aparm[1:] = aparm
	vplot.bparm[1:] = bparm
	vplot.crowded = 0
	vplot.refant = refant
	vplot.dotv = dotv
	vplot.xinc = 1
	vplot.nplots = nplots
	vplot.go()

def runquickvplot(uvdata):
	vplot = AIPSTask('VPLOT')
	vplot.indata = uvdata
	vplot.bchan = 3*uvdata.header['naxis'][2] /8 
	vplot.bchan = 5*uvdata.header['naxis'][2] /8
	vplot.solint = 0.5
	vplot.crowded = 3
	vplot.do3col = 1
	vplot.bparm[1] = 0
	vplot.bparm[2]=-2
	vplot.nplot = 6
	vplot.go()


def rungetjy(indata, sources, calsour, bif, eif, snver):
	getjy = AIPSTask('GETJY')
	getjy.indata = indata
	getjy.sources[1:] = sources
	getjy.calsour[1:] = calsour
	getjy.bif = bif
	getjy.eif = eif
	getjy.snver = snver
	getjy.go()


def runimagr(indata, sources, docalib, gainuse, flagver, doband, bpver, bchan, echan, nchav, chinc, cellsiz, imsiz, niter, dotv, outdisk):
	imagr = AIPSTask('IMAGR')
	imagr.indata = indata
	imagr.sources[1:] = sources
	source = str(sources[0])
	if len(source)>12 :
		source = source[0:12]
	imagr.outname = source
	imagr.outdisk = outdisk
	imagr.docalib = docalib
	imagr.gainuse = gainuse
	imagr.flagver = flagver
	imagr.doband = doband
	imagr.bpver = bpver
	imagr.bchan = bchan
	imagr.echan = echan
	imagr.nchav = nchav
	imagr.chinc = chinc
	imagr.cellsize[1:] = cellsiz
	imagr.imsize[1:] = imsiz
	imagr.niter = niter
	imagr.dotv = dotv
	imagr.go()


def runkntr(image, factor, userno):
	imean = AIPSTask('IMEAN')
	imean.userid = userno
	imean.indata = image
	imean()
	kntr = AIPSTask('KNTR')
	kntr.dogrey = -1
	kntr.dovect = -1
	kntr.docont = 1
	kntr.indata = image
	kntr.blc[1:] = 0,0
	kntr.trc[1:] = 0,0
#	kntr.clev = 3*imean.pixstd
#	kntr.levs[1:] = -1,1,2,4,6,8,16,32,64,128
	kntr.cbpl = 16
	kntr.dotv = -1
	kntr.plev = 1
	rms = imean.pixstd
	peak = image.header.datamax
	# Set the kntr levels based on dynamic range
	firstlev = factor*abs(rms/peak)*300.
	kntr.levs[1] = firstlev * (-1.)
	kntr.levs[2] = firstlev
	i = 3
	while (kntr.levs[i-1]*2. < 100. and i<=30):
		kntr.levs[i] = 2 * kntr.levs[i - 1]
		i += 1
	kntr.go()


def aipsuvname(aipsdata):
	return aipsdata.name + '.' + aipsdata.klass + '.' + str(aipsdata.seq)


def runfittp(uvdata, fittpdir, fittpfile):
	fittp = AIPSTask('FITTP')
	fittp.indata = uvdata
	fittp.doall = -1
	fittp.intype = ''
	fittp.outtape = 1
	srcname = re.sub(r"\s+", '-', fittpfile)
	fittp.dataout = os.path.join(fittpdir, srcname)
	print "Saving to disk: ", uvdata.name, uvdata.klass, uvdata.seq
	fittp.go()


def runtasav(uvdata, fittpdir, fittpfile, indisk):
	# Backup the FG tables after flagging
	tasav = AIPSTask('TASAV')
	tasavfil = AIPSUVData(uvdata.name, 'TASAV', indisk, uvdata.seq)
	tasav.indata = uvdata
	tasav.outdata = tasavfil
	tasav.outseq = uvdata.seq
	print "Saving tables for ", uvdata.name, uvdata.klass, uvdata.seq
	tasav.go()
	runfittp(tasavfil, fittpdir, fittpfile)
	tasavfil.zap()

def runsousp(indata, sources):
	sousp = AIPSTask('SOUSP')
	sousp.indata = indata
	sousp.sources[1:] = sources
	sousp.order = 1
	sousp.dotv = -1
	sousp.go()
	return sousp.specindx

def runuvfix(datafile, inname, indisk, thisdir):
	uvfix = AIPSTask('UVFIX')
	klass = datafile.klass
	uvfix.indata = datafile
	uvfix.outdata = AIPSUVData(datafile.name,'UVFIX',indisk,thisdir)
	try:
		uvfix.fqcenter = -1		# Added 20140521 due to a change in AIPS, this resets to old behaviour uvdata.header.crval[2]
	except AttributeError as detail:
		print "Your version of AIPS does not have this parameter,"
		print "so this does not matter"
	uvfix.go()
	datafile.zap()
	datafile = AIPSUVData(inname,'UVFIX',indisk,thisdir)
	datafile.rename(inname, klass, thisdir)

def runuvsrt(datafile, inname, indisk, thisdir):
	uvsrt = AIPSTask('UVSRT')
	klass = datafile.klass
	uvsrt.indata = datafile
	uvsrt.outdata = AIPSUVData(datafile.name,'UVSRT',indisk,datafile.seq)
	uvsrt.sort = 'TB'
	try:
		uvsrt.fqcenter = -1		# Added 20150918
	except AttributeError as detail:
		print "Your version of AIPS does not have this parameter,"
		print "so this does not matter"

	uvsrt.go()
	datafile.zap()
	datafile = AIPSUVData(inname,'UVSRT',indisk,thisdir)
	datafile.rename(inname,klass,thisdir)

def flagLOMK(datafile):
	# WRITE A FLAG FOR MK2-LT SPACING
	uvflg = AIPSTask('UVFLG')
	print "UVFLG: FLAGGING MK2-LT Baseline " + datafile.name + '.' + datafile.klass + '.' + format(datafile.seq)
	uvflg.indata = datafile
	uvflg.outfgver = 1
	uvflg.bif = 1
	uvflg.eif = 0
	uvflg.opcode = 'FLAG'
	uvflg.reason = 'Lovell - Mk2 baseline'
	uvflg.antennas[1] = 1 #Lovell
	uvflg.baseline[1] = 2 #MK2
	uvflg.go()

def runuvflg(datafile, bchan, echan, bif, eif, opcod, reason, antennas, baseline, sources, timer):
	"Run UVFLG with the specified parameters"
	uvflg = AIPSTask('UVFLG')
	uvflg.indata = datafile
	uvflg.outfgver = 1
	uvflg.bchan = bchan
	uvflg.echan = echan
	uvflg.bif = bif
	uvflg.eif = eif
	uvflg.opcode = opcod
	uvflg.reason = reason
	uvflg.antennas[1:] = antennas
	uvflg.baseline[1:] = baseline
	uvflg.sources[1:] = sources
	uvflg.timer[1:] = timer
	uvflg.go()

def runsplat(uvdata, outchan, tint, sbandl, sbandu, smootha, smoothb, smoothc, indisk):
	splat = AIPSTask('SPLAT')
	splat.indata = uvdata
	splat.outname = uvdata.name
	splat.outdisk = indisk
	splat.outclass = 'SPLAT'
	splat.outseq = uvdata.seq
	splat.solint = 0
	if outchan != -1 :
		aparm1 = 3
		splat.channel = uvdata.header['naxis'][2] / outchan
		splat.chinc = uvdata.header['naxis'][2] / outchan
	if tint != -1 :
		splat.solint = tint / 60
	splat.bif = sbandl # SBLANL -- Lower sub-band for extraction
	splat.eif = sbandu # SBANDU - Upper sub-band for extraction
	splat.douvcomp = -1
	splat.aparm[1] = 3
	splat.aparm[2] = 1
	# Smoothing options!?!?
	splat.smooth[1] = smootha # 1 = hanning default
	splat.smooth[2] = smoothb # 4 = hanning default
	splat.smooth[3] = smoothc # 1 = hanning default
	splat.flagver = 0
	try:
		splat.fqcenter = -1		# Added 20150918
	except AttributeError as detail:
		print "Your version of AIPS does not have this parameter,"
		print "so this does not matter"
	splat.go()

def runsplatmixed(uvdata, sbandl, sbandu, outcl, indisk):
	splat = AIPSTask('SPLAT')
	splat.indata = uvdata
	splat.outname = uvdata.name
	splat.outdisk = indisk
	splat.outclass = outcl
	splat.outseq = uvdata.seq
	splat.solint = 0
	splat.channel = -1
	splat.solint = -1
	splat.bif = sbandl # SBLANL -- Lower sub-band for extraction
	splat.eif = sbandu # SBANDU - Upper sub-band for extraction
	splat.douvcomp = -1
	splat.flagver = 0
	try:
		splat.fqcenter = -1		# Added 20150918
	except AttributeError as detail:
		print "Your version of AIPS does not have this parameter,"
		print "so this does not matter"
	splat.go()

def runprtmsg(prtask, outprint):
	prtmsg = AIPSTask('PRTMSG')
	clrmsg = AIPSTask('CLRMSG')
	prtmsg.prtask = prtask
	prtmsg.docrt = -1
	prtmsg.outprint = outprint
	prtmsg.go()
	clrmsg()

def runuvhgm(indata, sources, anten, basel, SEFDbif, SEFDeif, SEFDbchan, SEFDechan):
	uvhgm = AIPSTask('UVHGM')
	uvhgm.indata = indata
	uvhgm.sources[1:] = sources
	uvhgm.antennas[1:] = anten
	uvhgm.baseline[1:] = basel
	uvhgm.bif = SEFDbif
	uvhgm.eif = SEFDeif
	uvhgm.flagver = 0
	uvhgm.stokes ='HALF'
	uvhgm.bchan = SEFDbchan
	uvhgm.echan = SEFDechan
	uvhgm.doall = 1
	uvhgm.axtype ='H'
	uvhgm.timerang[1:] = 0, 0, 0, 0, 0, 0, 0, 0
#	uvhgm.pixrange[1:] = -20, 20
	uvhgm.dotv = -1
	uvhgm.docal = 1
	uvhgm.gainuse = 0
	uvhgm.pixrange[1:] = -10, 10
	uvhgm.doband = -1
	uvhgm.nboxes = 1000
	uvhgm.go()

def runlistr(indata,sources,optype,docrt,outprint):
	listr = AIPSTask('LISTR')
	listr.indata = indata
	listr.sources[1:] = sources
	listr.optype = optype
	listr.docrt = docrt
	listr.outprint = outprint
	listr.go()


def runclcor(indata, antn, indisk):
	"Run correction for incorrect Cambridge axis offset in data prior to mid-2015."
	clcor = AIPSTask('CLCOR')
	clcor.indata = indata
	clcor.antennas[1] = antn
	clcor.opcod = 'ANAX'
	clcor.clcorprm[1:] = 1.503,0
	clcor.inp()
	clcor.go()
	splat = AIPSTask('SPLAT')
	splat.indata = indata
	splat.docalib = 100
	splat.outname = indata.name
	splat.outdisk = indisk
	splat.outclass = 'SPLAT'
	splat.outseq = indata.seq
	try:
		splat.fqcenter = -1		# Added 20150918
	except AttributeError as detail:
		print "Your version of AIPS does not have this parameter,"
		print "so this does not matter"
	splat.inp()
	splat.go()
	uvdata = AIPSUVData(indata.name, 'SPLAT', indisk, indata.seq)
	indata.zap()
	uvdata.rename(uvdata.name,splat.inclass,int(splat.outseq))
	print "Cambridge axis offset correction applied."

def runmove(indata, indisk):
	move = AIPSTask('MOVE')
	move.indata = indata
	move.outdata = indata
	move.userid = AIPS.userno + 1
	move.opcode = 'MOVE'
	move.go()


################################################################################
# Mixed-mode stuff
################################################################################

def vbgluloop(filelist, indisk, contf):
	vbseq = 0
	vbglu = AIPSTask('VBGLU')
	while len(filelist) >1:
		vbglu.default()
		vbglu.indata  = AIPSUVData(filelist[0][0], filelist[0][1], indisk, filelist[0][2])
		vbglu.in2data = AIPSUVData(filelist[1][0], filelist[1][1], indisk, filelist[1][2])
		try:
			vbglu.fqcenter = -1
		except:
			pass
		try:
			vbglu.in3data = AIPSUVData(filelist[2][0], filelist[2][1], indisk, filelist[2][2])
			vbglu.in4data = AIPSUVData(filelist[3][0], filelist[3][1], indisk, filelist[3][2])
		except:
			pass
		vbglu.outname = filelist[0][0]
		vbglu.outclass = 'VBG' + str(int(vbseq))
		print "TYPE OUTSEQ", type(filelist[0][2])
		vbglu.outseq = math.floor(filelist[0][2])
		print "TYPE OUTSEQ", vbglu.outseq
		vbglu.outdisk = indisk
		vbglu.inp()
		vbglu.go()
		print "before we try to remove files, filelist =", filelist
		try:
			tmp = AIPSUVData(filelist[0][0], filelist[0][1], indisk, int(filelist[0][2]))
			print "removing ", tmp.name, tmp.klass, tmp.seq
			tmp.zap()
			filelist.popleft()
			tmp = AIPSUVData(filelist[0][0], filelist[0][1], indisk, int(filelist[0][2]))
			print "removing ", tmp.name, tmp.klass, tmp.seq
			tmp.zap()
			filelist.popleft()
			tmp = AIPSUVData(filelist[0][0], filelist[0][1], indisk, int(filelist[0][2]))
			print "removing ", tmp.name, tmp.klass, tmp.seq
			tmp.zap()
			filelist.popleft()
			tmp = AIPSUVData(filelist[0][0], filelist[0][1], indisk, int(filelist[0][2]))
			print "removing ", tmp.name, tmp.klass, tmp.seq
			tmp.zap()
			filelist.popleft()
		except:
			pass
		vbout = (vbglu.outname, vbglu.outclass, vbglu.outseq, contf)
		filelist.append(vbout)
		#print "after we try to remove files, filelist =", filelist
		vbseq += 1



################################################################################
# Other (non-AIPS) tasks
################################################################################

def get_tab(uvdata, table):
	# find the number of tables of a certain type
	ver = 0
	for i in range(len(uvdata.tables)) :
		if table in uvdata.tables[i][1] :
			ver = uvdata.tables[i][0]
	print "HIGHEST TABLE OF TYPE", table, "is", ver
	return ver


def get_ant_num(uvdata, refant_name):
	# convert antenna name to number
	antab = uvdata.table('AN',1)
	for row in antab :
		if refant_name in row.anname :
			return row.nosta


def set_refant(uvdata, refant):
	# Check the reference antenna is present, set to sensible default if not.
	# Returns referance antenna name, refant number, and a list of antennas present in the data.
	searchants = ['Mk2', 'Pi', 'Da', 'Kn', 'Cm', 'Lo', 'De']
	antab = uvdata.table('AN',1)
	refantn = 0
	refantn = get_ant_num(uvdata, refant)
	if not refantn :
		print "Warning: No refant specified or requested reference antenna not present."
		for item in searchants:
			if refantn:
				break
			for row in antab :
				if item in row.anname :
					print "Warning: Using reference antenna", item
					refant = item
					refantn = get_ant_num(uvdata, refant)
					break

	refantlist = []
	for antenna in searchants:
		antnum = get_ant_num(uvdata, antenna)
		if antnum:
			refantlist.append(antnum)

	return (refant, refantn, refantlist)


def nodot(item):
	"""Filter out hidden files."""
	return item[0] != '.'



def time2hms(seconds):
	# Function to convert seconds to hh:mm:ss.ss format, returns a string
	h=int(seconds/3600)
	m=int(seconds % 3600)/60
	s=seconds-(h*3600)-(m*60)
	h=`h`
	m=`m`
	s="%4.2f" % s
	hms=h.zfill(2)+":"+m.zfill(2)+":"+s.zfill(4)
	return hms




def inputinput(prompt):
	retries=1
	complaint='Yes or no, please!'
	while True :
		ok = raw_input(prompt)
		if ok == '' :
			#print "You wanted to continue but were too lazy to say so."
			print "OK, quitting."
			sys.exit()
		if ok in ('y', 'ye', 'yes','Y'):
			print "Thank you, continuing."
			break
		if ok in ('n', 'N', 'no', 'nop', 'nope'):
			print "OK, quitting."
			sys.exit()
		retries = retries - 1
		if retries < 0:
			print "Don't be daft."
			sys.exit()
		print complaint
	return


def findmaxb(uvdata):
	maxbaseline = 0
	antab = uvdata.table('AN',1)
	for row in antab :
		for row2 in antab :
			xsep = row.stabxyz[0] - row2.stabxyz[0]
			ysep = row.stabxyz[1] - row2.stabxyz[1]
			zsep = row.stabxyz[2] - row2.stabxyz[2]
			hypxy = math.sqrt((xsep * xsep) + (ysep * ysep))
			hypxyz = math.sqrt((zsep * zsep) + (hypxy * hypxy))
			if hypxyz > maxbaseline :
				maxbaseline = hypxyz
	cellsize = (1.22 * (300000000 / uvdata.header.crval[2]) / maxbaseline) / 3.141592 * 180 * 3600 / 5
	print "maxbaseline = ", maxbaseline, "cellsize = ", cellsize
	return cellsize,cellsize



def checkEqual(lst):
	"Check whether data are in mixed mode (True) or not (False)"
	return not len(set(lst)) == 1


def checkTime(uvdata):
	"Check the data to see if it genuinely out of time order."
	data=WizAIPSUVData(uvdata.name, uvdata.klass, uvdata.disk, uvdata.seq)
	itertime = 0
	sort = False
	for vis in data:
		if vis.time < itertime:
			sort = True
			break
		itertime = vis.time
	return sort



########################################################################################################################################
# Flux-calibration calculations
########################################################################################################################################


def dfluxpy(freq,uvdata):

	"Function to calculate 3C286 values, modified from Danielle's dfluxpy.py."
	"Updated 20140506 to include new coefficients and a calculation for projected"
	"baseline length from inspecting the u's and v's."


	data=WizAIPSUVData(uvdata.name, uvdata.klass, uvdata.disk, uvdata.seq)

	antab = data.table('AN',1)
	numLO = 10
	numM2 = 10
	for row in antab:			# Ignore the Lovell baselines
		if 'Lo' in row.anname:
			numLO = row.nosta
		if 'Mk2' in row.anname:
			numM2 = row.nosta

	u=[]
	v=[]
	proj=[]
	basel = []
	baseline = 10000000
	for visibility in data:
		if ( visibility.baseline != [numLO, numM2] ) and ( visibility.baseline != [numM2, numLO] ):
			u.append(visibility.uvw[0])
			v.append(visibility.uvw[1])
			newbasel = (299792458.0/data.header.crval[2]) * math.sqrt((visibility.uvw[0] ** 2) + (visibility.uvw[1] ** 2))
			proj.append( newbasel )
			if newbasel < baseline:
				baseline = newbasel
				basel = visibility.baseline


	antab = data.table('AN',1)
	for row in antab:
		if row.nosta == basel[0]:
			ant1=row.anname

	for row in antab:
		if row.nosta == basel[1]:
			ant2=row.anname

	print "For projected baseline", baseline/1000, "km, between", ant1, "and", ant2


	# Perley & Butler 2012 values
	A = 1.2515
	B = -0.4605
	C = -0.1715
	D = 0.0336


	log10f = (math.log(freq)/2.3025851) - 3.0; # Why the -3? Because freq has to be GHz for the formula to work.
	log_flux = A + B*log10f + C*log10f*log10f + D*log10f*log10f*log10f
	vlaflux = math.pow(10.0,log_flux)

	ref_bl_length = 11236.79 # MK-TA separation in metres.
	ref_freq = 5000.0
	ref_rho = 0.04

	bl_length = baseline

	frac = (freq / ref_freq) * (bl_length / ref_bl_length)
	rho = frac * frac * ref_rho
	merlinflux = vlaflux / (1.0 + rho)

	print "\tfor IF with freq =", freq, ", e-MERLIN flux =", merlinflux

	return merlinflux



########################################################################################################################################
# Flagging tasks
########################################################################################################################################


def flagmask(uvdata):

	"Function to apply flags from a mask of known bad frequencies, derived by the e-MERLIN team."

	afile = 'flagmask512.fg'
	if os.path.isfile(afile) :
		flagfile = open(afile, "r")
	else :
		print "Error:" + afile + "does not exist, don't know what to flag."
		print "Error: If you do not have this file, ask the e-MERLIN team."
		sys.exit()

	start = 0

	space = re.compile(r'\s')
	newline = re.compile(r'\n')

	sources = ['','']
	timer = 0,0,0,0,999,23,59,59

	fqtab=uvdata.table('FQ',1)
	band_strt = uvdata.header.crval[2]
	nchan = uvdata.header.naxis[2]
	if uvdata.header.naxis[3] > 1:
		ifstep = fqtab[0].total_bandwidth[0]
		band_stop = uvdata.header.crval[2] + (len(fqtab[0].total_bandwidth) * ifstep)
		chan_wdth = fqtab[0].total_bandwidth[0] / nchan
	else:
		ifstep = fqtab[0].total_bandwidth
		band_stop = uvdata.header.crval[2] + ifstep
		chan_wdth = fqtab[0].total_bandwidth / nchan


	for line in flagfile:
		line = newline.sub(r'', line)
		if start == 0:
			print line
		if line == "***BEGIN*PASS***":
			start = 1
		elif line == "***END*PASS***":
			print line
			start = 0
		elif start == 1:
			freq1 = float(space.sub(r'', line[87:98]))
			#print freq1, band_strt, band_stop, chan_wdth, ifstep
			if (freq1 >= band_strt) and (freq1 <= band_stop):
				ant1 = int(space.sub(r'', line[50:52]))
				if uvdata.header.naxis[3] > 1:
					ifs1 = int( 1 + ((freq1 - uvdata.header.crval[2]) / ifstep) )
					chans1 = int( (freq1 - (band_strt + (ifstep * (ifs1 - 1)))) / chan_wdth)
				else:
					ifs1 = 1
					chans1 = int( (freq1 - band_strt) / chan_wdth)

				start = 2
				print line[0:49], str(ant1).rjust(2), line[53:76], str(ifs1).rjust(2), line[81:86], str(chans1).rjust(12), line[99:140]
				pass
		elif start == 2:
			freq2 = float(space.sub(r'', line[87:98]))
			if (freq2 >= band_strt) and (freq2 <= band_stop):
				ant2 = int(space.sub(r'', line[50:52]))
				if uvdata.header.naxis[3] > 1:
					ifs2 = int( 1 + ((freq2 - uvdata.header.crval[2]) / ifstep) )
					chans2 = int( (freq2 - (band_strt + (ifstep * (ifs2 - 1)))) / chan_wdth)
					if chans2 == 0:
						chans2 = nchan
				else:
					ifs2 = 1
					chans2 = int( (freq2 - band_strt) / chan_wdth)
					
				# Deal with commands which cross a band edge
				if chans2 < chans1 :
					
					print line[0:49], str(ant2).rjust(2), line[53:76], str(ifs1).rjust(2), line[81:86], str(nchan).rjust(12), line[99:140]
					print line[0:49], str(ant1).rjust(2), line[53:76], str(ifs2).rjust(2), line[81:86], str(1).rjust(12), line[99:140]
					print line[0:49], str(ant2).rjust(2), line[53:76], str(ifs2).rjust(2), line[81:86], str(chans2).rjust(12), line[99:140]

					runuvflg(uvdata, chans1, nchan, ifs1, ifs1, 'FLAG', 'e-MERLIN flag mask', [ant1,0], [ant2,0], sources, timer)
					
					runuvflg(uvdata, 1, chans2, ifs2, ifs2, 'FLAG', 'e-MERLIN flag mask', [ant1,0], [ant2,0], sources, timer)
					
				else:
					print line[0:49], str(ant2).rjust(2), line[53:76], str(ifs2).rjust(2), line[81:86], str(chans2).rjust(12), line[99:140]
					runuvflg(uvdata, chans1, chans2, ifs1, ifs2, 'FLAG', 'e-MERLIN flag mask', [ant1,0], [ant2,0], sources, timer)

				start = 1

	return


########################################################################################################################################
# Wide-field tasks
########################################################################################################################################

def runsplit(uvdata,widetarget):
	# This task splits off each target from the final uvdataset
	# first need to SPLIT off each target and give it the name of the source
	print "Splitting sources off muilt-source file."
	split = AIPSTask('SPLIT')
	split.indata = uvdata
	split.docalib=100
	split.outclass='SPLIT' # set the outclass to 'SPLIT'
	split.outdisk=uvdata.disk
	split.outseq=0
	split.gainuse=0
	split.flagver=0
	split.doband=-1
	split.bchan=1
	split.echan=0
	split.bif=1
	split.eif=0
	split.bpver=1
	split.sources[1] = widetarget	#choose the source from widetarget
	try:
		split.fqcenter = -1		# Added 20150918
	except AttributeError as detail:
		print "Your version of AIPS does not have this parameter,"
		print "so this does not matter"
	split.go()


def chessboard64wide(uvdata,radius):
	# This function creates a chessboardlike set of smaller uvfiles.
	# Can be used for rapid imaging.	
	# radius is the half the half power beamwidth of the array.
	# uvdata is a list containing name,class,disk,sequence
	print "Chopping up uv data into chessboard-like uv datasets. This will take some time..."

	# set main input parameters
	uvfix = AIPSTask('UVFIX')
	uvfix.indata=uvdata
	try:
		uvfix.fqcenter = -1		# Added 20150918
	except AttributeError as detail:
		print "Your version of AIPS does not have this parameter,"
		print "so this does not matter"

	#now shift the phase centre to make 4 quadrants using uvfix
	print "Shifting the phase centre to make quadrant 1."
	uvfix.shift[1]=radius*60/2
	uvfix.shift[2]=radius*60/2 
	uvfix.outdata=AIPSUVData(uvdata.name,'UV1000',uvdata.disk,1)
	uvfix.go()
	print "Shifting the phase centre to make quadrant 2."
	uvfix.shift[1]=-radius*60/2
	uvfix.shift[2]=radius*60/2
	uvfix.outdata=AIPSUVData(uvdata.name,'UV2000',uvdata.disk,1)
	uvfix.go()
	print "Shifting the phase centre to make quadrant 3."
	uvfix.shift[1]=radius*60/2
	uvfix.shift[2]=-radius*60/2
	uvfix.outdata=AIPSUVData(uvdata.name,'UV3000',uvdata.disk,1)
	uvfix.go()
	print "Shifting the phase centre to make quadrant 4."
	uvfix.shift[1]=-radius*60/2
	uvfix.shift[2]=-radius*60/2
	uvfix.outdata=AIPSUVData(uvdata.name,'UV4000',uvdata.disk,1)
	uvfix.go()

	# now average each quadrant in time
	uvavg=AIPSTask('UVAVG')
	uvavg.yinc=2
	for i in range(1,5):
		klass="UV"+ str(i) +"000"
		uvavg.indata=AIPSUVData(uvdata.name,klass,uvdata.disk,1)
		uvavg.outdata=AIPSUVData(uvdata.name,klass,uvdata.disk,2)
		print "Time averaging quadrant "+str(i)
		try:
			uvavg.fqcenter = -1		# Added 20150918
		except AttributeError as detail:
			print "Your version of AIPS does not have this parameter,"
			print "so this does not matter"
		uvavg.go()
		uvdata = AIPSUVData(uvdata.name,klass,uvdata.disk,1)	# zap the input
		uvdata.zap()
		uvdata = AIPSUVData(uvdata.name,klass,uvdata.disk,2)
		uvdata.rename(uvdata.name,klass,uvdata.disk)

	# now average each quadrant in frequency
	avspc=AIPSTask('AVSPC')
	avspc.avoption='SUBS'
	avspc.channel=2 #averaging
	try:
		avspc.fqcenter = -1		# Added 20150918
	except AttributeError as detail:
		print "Your version of AIPS does not have this parameter,"
		print "so this does not matter"
	'''
	for i in range(1,5):
		klass="UV"+ str(i) +"000"
		avspc.indata=AIPSUVData(uvdata.name,klass,uvdata.disk,1)
		avspc.outdata=AIPSUVData(uvdata.name,klass,uvdata.disk,2)
		avspc.go()
		uvdata = AIPSUVData(uvdata.name,klass,uvdata.disk,1)	# zap the input
		uvdata.zap()
		uvdata = AIPSUVData(uvdata.name,klass,uvdata.disk,2)
		uvdata.rename(uvdata.name,klass,uvdata.disk)
	'''


	# **********************************************************************
	# Now shift the phase centre for each quadrant again making 16 files
	# **********************************************************************
	
	for i in range(1,5):
		klass='UV'+str(i)+'000'
		outkl='UV'+str(i)+'100'
		try:
			uvfix.fqcenter = -1		# Added 20150918
		except AttributeError as detail:
			print "Your version of AIPS does not have this parameter,"
			print "so this does not matter"

		aipsuvdata=AIPSUVData(uvdata.name,klass,uvdata.disk,1)
		uvfix.indata=aipsuvdata
		uvfix.outdata=AIPSUVData(uvdata.name,'UV'+str(i)+'100',uvdata.disk,1)
		uvfix.shift[1]=radius*60/4
		uvfix.shift[2]=radius*60/4
		uvfix.go()

		uvfix.outdata=AIPSUVData(uvdata.name,'UV'+str(i)+'200',uvdata.disk,1)
		uvfix.shift[1]=-radius*60/4
		uvfix.shift[2]=radius*60/4
		uvfix.go()

		uvfix.outdata=AIPSUVData(uvdata.name,'UV'+str(i)+'300',uvdata.disk,1)
		uvfix.shift[1]=radius*60/4
		uvfix.shift[2]=-radius*60/4
		uvfix.go()

		uvfix.outdata=AIPSUVData(uvdata.name,'UV'+str(i)+'400',uvdata.disk,1)
		uvfix.shift[1]=-radius*60/4
		uvfix.shift[2]=-radius*60/4
		uvfix.go()
		uvdata = AIPSUVData(uvdata.name,klass,uvdata.disk,1)	# zap the input
		uvdata.zap()


	# **********************************************************************	
	# Now average these 16 files in time to 4 seconds
	uvavg.yinc=4 #desired average in seconds
	uvavg.zinc=2 #previous average in seconds
	try:
		uvavg.fqcenter = -1		# Added 20150918
	except AttributeError as detail:
		print "Your version of AIPS does not have this parameter,"
		print "so this does not matter"

	for i in range(1,5):
		for j in range(1,5):
			uvavg.indata=AIPSUVData(uvdata.name,"UV"+ str(i) + str(j)+"00",uvdata.disk,1)
			uvavg.outdata=AIPSUVData(uvdata.name,"UV"+ str(i) + str(j)+"00",uvdata.disk,2)
			print "Time averaging subquadrant "+str(i)
			uvavg.go()
			uvdata=AIPSUVData(uvdata.name,"UV"+ str(i) + str(j)+"00",uvdata.disk,1)
			uvdata.zap()
			uvdata=AIPSUVData(uvdata.name,"UV"+ str(i) + str(j)+"00",uvdata.disk,2)
			uvdata.rename(uvdata.name,"UV"+ str(i) + str(j)+"00",uvdata.disk)

	# **********************************************************************
	# Now channel average by a factor of 2 again
	avspc.channel=2
	try:
		avspc.fqcenter = -1		# Added 20150918
	except AttributeError as detail:
		print "Your version of AIPS does not have this parameter,"
		print "so this does not matter"
	for i in range(1,5):
		for j in range(1,5):
			avspc.indata=AIPSUVData(uvdata.name,"UV"+ str(i) + str(j)+"00",uvdata.disk,1)
			avspc.outdata=AIPSUVData(uvdata.name,"UV"+ str(i) + str(j)+"00",uvdata.disk,2)
			print "Averaging channels in subquadrant"
			avspc.go()
			uvdata=AIPSUVData(uvdata.name,"UV"+ str(i) + str(j)+"00",uvdata.disk,1)
			uvdata.zap()
			uvdata=AIPSUVData(uvdata.name,"UV"+ str(i) + str(j)+"00",uvdata.disk,2)
			uvdata.rename(uvdata.name,"UV"+ str(i) + str(j)+"00",uvdata.disk)



	# **********************************************************************
	# Now shift the phase centre for each quadrant again making 64 files
	# **********************************************************************
	
	for i in range(1,5):
		for j in range(1,5):
			aipsuvdata=AIPSUVData(uvdata.name,'UV'+str(i)+str(j)+'00',uvdata.disk,1)
			uvfix.indata=aipsuvdata
			uvfix.outdata=AIPSUVData(uvdata.name,'UV'+str(i)+str(j)+'10',uvdata.disk,1)
			uvfix.shift[1]=radius*60/8
			uvfix.shift[2]=radius*60/8
			uvfix.go()

			uvfix.outdata=AIPSUVData(uvdata.name,'UV'+str(i)+str(j)+'20',uvdata.disk,1)
			uvfix.shift[1]=-radius*60/8
			uvfix.shift[2]=radius*60/8
			uvfix.go()

			uvfix.outdata=AIPSUVData(uvdata.name,'UV'+str(i)+str(j)+'30',uvdata.disk,1)
			uvfix.shift[1]=radius*60/8
			uvfix.shift[2]=-radius*60/8
			uvfix.go()

			uvfix.outdata=AIPSUVData(uvdata.name,'UV'+str(i)+str(j)+'40',uvdata.disk,1)
			uvfix.shift[1]=-radius*60/8
			uvfix.shift[2]=-radius*60/8
			uvfix.go()

			uvdata = AIPSUVData(uvdata.name,'UV'+str(i)+str(j)+'00',uvdata.disk,1)
			uvdata.zap()


	# **********************************************************************	
	# Now average these 64 files in time to 8 seconds
	uvavg.yinc=8 #desired average in seconds
	uvavg.zinc=4 #previous average in seconds
	for i in range(1,5):
		for j in range(1,5):
			for k in range(1,5):
				uvavg.indata=AIPSUVData(uvdata.name,"UV"+ str(i) + str(j)+str(k)+"0",uvdata.disk,1)
				uvavg.outdata=AIPSUVData(uvdata.name,"UV"+ str(i) + str(j)+str(k)+"0",uvdata.disk,2)
				print "Time averaging subsubquadrant "+str(i)+str(j)+str(k)
				uvavg.go()
				uvdata=AIPSUVData(uvdata.name,"UV"+ str(i) + str(j)+str(k)+"0",uvdata.disk,1)
				uvdata.zap()
				uvdata=AIPSUVData(uvdata.name,"UV"+ str(i) + str(j)+str(k)+"0",uvdata.disk,2)
				uvdata.rename(uvdata.name,"UV"+ str(i) + str(j) + str(k) + "0", uvdata.disk)


	# **********************************************************************
	# Now channel average by a factor of 2 again
	avspc.channel=2
	for i in range(1,5):
		for j in range(1,5):
			for k in range(1,5):
				avspc.indata=AIPSUVData(uvdata.name,"UV"+ str(i) + str(j)+str(k)+"0",uvdata.disk,1)
				avspc.outdata=AIPSUVData(uvdata.name,"UV"+ str(i) + str(j)+str(k)+"0",uvdata.disk,2)
				print "Averaging channels in subsubquadrant"
				avspc.go()
				uvdata = AIPSUVData(uvdata.name,"UV"+ str(i) + str(j)+str(k)+"0",uvdata.disk,1)
				uvdata.zap()
				uvdata = AIPSUVData(uvdata.name,"UV"+ str(i) + str(j)+str(k)+"0",uvdata.disk,2)
				uvdata.rename(uvdata.name,"UV"+ str(i) + str(j)+str(k)+"0",uvdata.disk)

	print "Completed averaging of " + uvdata.name + '.' + uvdata.klass + '.' + format(uvdata.seq) + " into 64 uvfiles."


	# **********************************************************************
	# Now shift the phase centre for each quadrant again making 256 files
	# **********************************************************************
	
	for i in range(1,5):
		for j in range(1,5):
			for k in range(1,5):
				aipsuvdata=AIPSUVData(uvdata.name,'UV'+str(i)+str(j)+str(k)+'0',uvdata.disk,1)
				uvfix.indata=aipsuvdata
				uvfix.outdata=AIPSUVData(uvdata.name,'UV'+str(i)+str(j)+str(k)+'1',uvdata.disk,1)
				uvfix.shift[1]=radius*60/16
				uvfix.shift[2]=radius*60/16
				uvfix.go()

				uvfix.outdata=AIPSUVData(uvdata.name,'UV'+str(i)+str(j)+str(k)+'2',uvdata.disk,1)
				uvfix.shift[1]=-radius*60/16
				uvfix.shift[2]=radius*60/16
				uvfix.go()

				uvfix.outdata=AIPSUVData(uvdata.name,'UV'+str(i)+str(j)+str(k)+'3',uvdata.disk,1)
				uvfix.shift[1]=radius*60/16
				uvfix.shift[2]=-radius*60/16
				uvfix.go()

				uvfix.outdata=AIPSUVData(uvdata.name,'UV'+str(i)+str(j)+str(k)+'4',uvdata.disk,1)
				uvfix.shift[1]=-radius*60/16
				uvfix.shift[2]=-radius*60/16
				uvfix.go()

				uvdata = AIPSUVData(uvdata.name,'UV'+str(i)+str(j)+str(k)+'0',uvdata.disk,1)
				uvdata.zap()


	# **********************************************************************	
	# Now average these 256 files in time to 16 seconds
	uvavg.yinc=16 #desired average in seconds
	uvavg.zinc=8 #previous average in seconds
	for i in range(1,5):
		for j in range(1,5):
			for k in range(1,5):
				for l in range(1,5):
					uvavg.indata=AIPSUVData(uvdata.name,"UV"+ str(i) + str(j)+str(k)+str(l),uvdata.disk,1)
					uvavg.outdata=AIPSUVData(uvdata.name,"UV"+ str(i) + str(j)+str(k)+str(l),uvdata.disk,2)
					print "Time averaging subsubquadrant "+str(i)+str(j)+str(k)+str(l)
					uvavg.go()
					uvdata=AIPSUVData(uvdata.name,"UV"+ str(i) + str(j)+str(k)+str(l),uvdata.disk,1)
					uvdata.zap()
					uvdata=AIPSUVData(uvdata.name,"UV"+ str(i) + str(j)+str(k)+str(l),uvdata.disk,2)
					uvdata.rename(uvdata.name,"UV"+ str(i) + str(j) + str(k) + str(l), uvdata.disk)


	# **********************************************************************
	# Now channel average by a factor of 2 again
	avspc.channel=2
	for i in range(1,5):
		for j in range(1,5):
			for k in range(1,5):
				for l in range(1,5):
					avspc.indata=AIPSUVData(uvdata.name,"UV"+ str(i) + str(j)+str(k)+str(l),uvdata.disk,1)
					avspc.outdata=AIPSUVData(uvdata.name,"UV"+ str(i) + str(j)+str(k)+str(l),uvdata.disk,2)
					print "Averaging channels in subsubquadrant "+str(i)+str(j)+str(k)+str(l)
					avspc.go()
					uvdata = AIPSUVData(uvdata.name,"UV"+ str(i) + str(j)+str(k)+str(l),uvdata.disk,1)
					uvdata.zap()
					uvdata = AIPSUVData(uvdata.name,"UV"+ str(i) + str(j)+str(k)+str(l),uvdata.disk,2)
					uvdata.rename(uvdata.name,"UV"+ str(i) + str(j)+str(k)+str(l),uvdata.disk)

	print "Completed averaging of " + uvdata.name + '.' + uvdata.klass + '.' + format(uvdata.seq) + " into 256 uvfiles."


def widemaps(uvdata,indisk,radius,target,noiter):
	# This images each uv file created in chessboard
	# First set up imaging defaults	
	# Enter permanent parameters here:

	# The next 10 lines moved here from the main code.  Added findmaxb() call to calculate correct pixel size for any freq.  mkargo 20140515
	cellsiz = findmaxb(uvdata)
	cellpitch = cellsiz[0]
	cellpitch=0.1	#this overrides the previous line

	imgsize = int((radius*2*60/(cellpitch*8)))
	thing = imgsize
	count=2
	#print "imgsize =", imgsize, "i =", count, "radius =", radius, "cellpitch =", cellpitch
	while thing > 1.0 :
		thing = thing / 2
		count = count + 1
	imgsize = 2**count
	#print "imgsize =", imgsize, "i =", count, "radius =", radius, "cellpitch =", cellpitch

	docalib=-1
	gainuse=0
	flagver=0
	doband=-1
	bpver=0
	bchan=1
	echan=0
	nchav=16
	chinc=nchav
	dotv=-1

	imagr = AIPSTask('IMAGR')
	imagr.sources[1] = target
	imagr.outname = uvdata.name
	imagr.outdisk = indisk
	imagr.docalib = docalib
	imagr.gainuse = gainuse
	imagr.flagver = flagver
	imagr.doband = doband
	imagr.bpver = bpver
	imagr.bchan = bchan
	imagr.echan = echan
	imagr.nchav = nchav
	imagr.chinc = chinc
	imagr.cellsize[1] = cellpitch
	imagr.cellsize[2] = cellpitch
	imagr.imsize[1] = imgsize
	imagr.imsize[2] = imgsize
#	imagr.imagrprm[1] = 48      #effective array antenna diameter
	imagr.imagrprm[10:]= 1,0,0,0 #make the dirty beams the same size as the maps
	#imagr.uvrang=
	imagr.niter = 0#noiter
	imagr.dotv = dotv
	imagr.robust = 0
#	imagr.flux = 1.5e-5
#	imagr.bmaj = 0.4
#	imagr.bmin = 0.4
	imagr.bmaj = 0.4 #make sure that the convolution beam size is the same for all IFs
	imagr.bmin = 0.4
	indata=uvdata
	print "Now imaging field from UV file " + uvdata.name + '.' + uvdata.klass + '.' + format(uvdata.seq)
	imagr.indata = indata
	imagr.go()


def beamgain(freq,throw,pbparms):
	#this calculates the attenuation at an angle throw armin for freq in GHz using pbparms
	p=1+pbparms[0]*1e-3*(throw*freq)**2 + pbparms[1]*1e-7*(throw*freq)**4 + pbparms[2]*1e-10*(throw*freq)**6 + pbparms[3]*1e-13*(throw*freq)**8 + pbparms[4]*1e-16*(throw*freq)**10
	return p



def widebandmaps(uvdata,indisk,radius,target,noiter):
	# This generates multiple dirty beams and images using APCLN for each uv file created in chessboard
	# It creates accurate Point Spread Functions accounting for induced spectral slope caused by primary beam
	# First set up imaging defaults	
	# Enter permanent parameters here:

	# The next 10 lines moved here from the main code.  Added findmaxb() call to calculate correct pixel size for any freq.  mkargo 20140515
	cellsiz = findmaxb(uvdata)
	cellpitch = cellsiz[0]
	cellpitch=0.045	#this overrides the previous line which actually gets the size wrong
	facetnumber=255
	
	
	imgsize = int((radius*2*60/(cellpitch*(facetnumber)**(0.5))))
	#now modify the image size
	thing = imgsize
	count=2
	#print "imgsize =", imgsize, "i =", count, "radius =", radius, "cellpitch =", cellpitch
	while thing > 1.0 :
		thing = thing / 2
		count = count + 1
	imgsize = 2**count
	#print "imgsize =", imgsize, "i =", count, "radius =", radius, "cellpitch =", cellpitch

	#find out relevant info about the UV file
	nifs=uvdata.header['naxis'][3]		#number of IFs or spectral windows in dataset
	basefreq=uvdata.header['crval'][2]	#frequency of IF 1 in Hz
	atpixel=uvdata.header['crpix'][2]	#at pixel (channel) number
	RAobs=uvdata.header['obsra']		#RA of pointing centre in degrees
	DECobs=uvdata.header['obsdec']		#DEC of pointing centre in degrees
	RAphase=uvdata.header['crval'][4]	#RA of phase centre in degrees
	DECphase=uvdata.header['crval'][5]	#DEC of phase centre in degrees
	channels=uvdata.header['naxis'][2]	#number of channels per IF
	chanwidth=uvdata.header['cdelt'][2]	#bandwidth of each channel in Hz
	centrefreq=float((basefreq-atpixel*chanwidth) + (nifs/2 * chanwidth*channels)+(chanwidth*atpixel))
	#calculate approximate throw from pointing centre to use for finding primary beam attenuation
	deltaRA=numpy.abs((RAobs-RAphase)*numpy.cos(DECobs*3.14159265359/180))
	deltaDEC=numpy.abs(DECobs-DECphase)
	deltadist=60*((numpy.abs(deltaRA**2+deltaDEC**2)**0.5))	#distance from pc in arcmin
	print "Throw from pointing centre of " + uvdata.name + " is " + format(deltadist) + " arcmin."


	#create dirty beams & maps *************************
	docalib=-1
	gainuse=0
	flagver=0
	doband=-1
	bpver=0
	bchan=1
	echan=0
	nchav=channels
	chinc=nchav
	dotv=-1

	imagr = AIPSTask('IMAGR')
	imagr.sources[1] = target
	imagr.outname = uvdata.name
	imagr.outdisk = indisk
	imagr.docalib = docalib
	imagr.gainuse = gainuse
	imagr.flagver = flagver
	imagr.doband = doband
	imagr.bpver = bpver
	imagr.bchan = bchan
	imagr.echan = echan
	imagr.nchav = nchav
	imagr.chinc = chinc
	imagr.cellsize[1] = cellpitch
	imagr.cellsize[2] = cellpitch
	imagr.imsize[1] = imgsize
	imagr.imsize[2] = imgsize
	#imagr.antennas[1]=-1	#omit the Lovell
	#imagr.antennas[2]=0
#	imagr.imagrprm[1] = 48      #effective array antenna diameter only used for built-in (unused)
	imagr.imagrprm[10:]= 1,0,0,0 #make the dirty beams the same size as the maps

	imagr.niter = 0	#this forces dirty maps
	imagr.dotv = dotv
	imagr.robust = 7
#	imagr.flux = 1.5e-5
	imagr.bmaj = 0.2 #make sure that the restoration convolution beam size is the same for all IFs
	imagr.bmin = 0.2
	indata=uvdata
	print "Now creating dirty beams from UV file " + uvdata.name + '.' + uvdata.klass + '.' + format(uvdata.seq)
	imagr.indata = indata

	#temp set pbparms here for the amount of synth beam scaling

	#pbparms=[-5.6712,191.83,-378.306,303.408,0.186]#e-MERLIN LO=30
	pbparms=[-5.8722,195.198,-390.066,333.9,1.098]#e-MERLIN LO=50
	#pbparms=[-6.3174,242.64,-600.162,597.27,100.026]#e-MERLIN LO=70
	#pbparms=[-1.343,6.579,-1.186,0,0]#JVLA 
	
	scaledown=[]
	#loop all IFs to create dirty beams and maps per IF
	for i in range (nifs): #from zero to number of IFs - 1
		#create scaling factors for primary beam
		#centre frequency of IF i+1
		frequency=1e-9*(basefreq+(channels*chanwidth*(i)))	#calculate the frequency of this IF in GHz
		scaledown.append(beamgain(frequency,deltadist,pbparms))	#calculate the scaling factor by calling for primary beamshape
		print "Primary beam gain in IF" + format(i+1) + " (" + format(frequency) + "GHz) at " + format(deltadist) + " arcmin is " + format(scaledown[i]) + " of " + format(nifs) + " IFs"
		imagr.bif = i+1
		imagr.eif = i+1
		imagr.outseq=i+1
		imagr.go()
	

	#now scale and combine the dirty beams *********************************
	comb = AIPSTask('COMB')
	comb.inname=uvdata.name
	comb.in2name=uvdata.name
	comb.indisk=indisk
	comb.outdisk=indisk
	comb.opcode='SUM'
	ifs=float(nifs)	#convert to a number

	#this is the new way of doing it using mcube and sqash
	for i in range(nifs): #scale all the beams first
		comb.inclass='IBM001'
		comb.in2class='IBM001'
		comb.inseq=i+1
		comb.in2seq=i+1
		comb.aparm[1]=0.5*(1/float(nifs))*float(scaledown[i])
		comb.aparm[2]=0.5*(1/float(nifs))*float(scaledown[i])
		comb.outseq=i+1
		comb.outclass='IBMSCL'
		print 'Scaling beam in IF ' + format(i+1) + ' for facet ' + format(uvdata.klass)
		comb.go()

	#now put all the scaled beams into a cube, they should all have class ibm002
	mcube=AIPSTask('MCUBE')
	mcube.inclass='IBMSCL'
	mcube.inname=uvdata.name
	#mcube.in2name=uvdata.name
	mcube.indisk=indisk
	mcube.outdisk=indisk
	mcube.outclass='IBMCUB'
	mcube.outname=uvdata.name
	mcube.inseq=1 #start of sequence
	mcube.outseq=0
	mcube.in2seq=int(float(nifs))
	mcube.npoints=int(float(nifs))
	mcube.axref=0
	mcube.ax2ref=0
	print 'Building Beam Cube for facet ' + format(uvdata.klass)
	mcube.go()

	#now collapse the cube
	sqash=AIPSTask('SQASH')
	sqash.inclass='IBMCUB'
	sqash.outclass='IBMSQA'
	sqash.inname=uvdata.name
	sqash.outname=uvdata.name
	sqash.indisk=indisk
	sqash.outdisk=indisk
	sqash.inseq=0
	sqash.dparm[1]=4 #sum
	sqash.bdrop=3#collapse in frequency
	sqash.dparm[2]=1 #ignore blank pixels
	sqash.outseq=0
	print 'Collapsing Beam Cube for facet ' + format(uvdata.klass)
	sqash.go()


	#now normalise the beam	
	#ncomb.inclass='IBM001'
	ima=AIPSImage(uvdata.name,'IBMSQA',indisk,1)
	peakvalue=ima.header.datamax
	comb.inname=uvdata.name
	comb.in2name=uvdata.name
	comb.inclass='IBMSQA'
	comb.in2class='IBMSQA'
	comb.inseq=0
	comb.in2seq=0
	comb.aparm[1]=0.5*(1/peakvalue)#(1/float(nifs))*float(scaledown[i])
	comb.aparm[2]=0.5*(1/peakvalue)#(1/float(nifs))*float(scaledown[i])
	comb.outseq=0
	comb.outclass='IBMNOR'
	print 'Normalising Beam for facet ' + format(uvdata.klass)
	comb.go()

	#repair the frequency in the header
	#puthead=AIPSTask('PUTHEAD')
	from Wizardry.AIPSData import AIPSImage as WizAIPSImage
	imagebeam=WizAIPSImage(uvdata.name,'IBMNOR',indisk,1)#+ format(int(numpy.log(ifs)/numpy.log(2))+1),indisk,1)
	imagebeam.header.crval[2]=centrefreq
	imagebeam.header.update()
	#imagecrval=imagebeam.header['crval']
	#imagecrval[2]=centrefreq
	#basefreq=uvdata.header['crval']
	#puthead.indisk=indisk
	#puthead.inname=uvdata.name
	#puthead.inclass='IBM00'+ format(int(numpy.log(ifs)/numpy.log(2))+1)
	#puthead.inseq=1
	#puthead.keyword='CRVAL3'
	#puthead.keyvalue[1]=centrefreq
	#puthead.keyvalue[2]=0
	#puthead.inp()
	#puthead.go()
	
		
	#now zap all but the IBMNOR beam
	for i in range(nifs):
		image = AIPSImage(uvdata.name,'IBM001',indisk,i+1) #zapping the original dirty beams
		image.zap()
		image = AIPSImage(uvdata.name,'IBMSCL',indisk,i+1) #zapping the scaled dirty beams
		image.zap()

	image = AIPSImage(uvdata.name,'IBMCUB',indisk,1) #zapping the cube of the dirty beams
	image.zap()
	image = AIPSImage(uvdata.name,'IBMSQA',indisk,1) #zapping the squashed cube of the dirty beams
	image.zap()
	

	#now do the dirty maps ****************************************************

	#now put all the scaled beams into a cube, they should all have class ibm002
	mcube=AIPSTask('MCUBE')
	mcube.inclass='IIM001'
	mcube.inname=uvdata.name
	#mcube.in2name=uvdata.name
	mcube.indisk=indisk
	mcube.outdisk=indisk
	mcube.outclass='IIMCUB'
	mcube.outname=uvdata.name
	mcube.inseq=1 #start of sequence
	mcube.outseq=0
	mcube.in2seq=int(float(nifs))
	mcube.npoints=int(float(nifs))
	mcube.axref=0
	mcube.ax2ref=0
	print 'Building Image Cube for facet ' + format(uvdata.klass)
	mcube.go()

	#now collapse the cube
	sqash=AIPSTask('SQASH')
	sqash.inclass='IIMCUB'
	sqash.outclass='IIMNOR'
	sqash.inname=uvdata.name
	sqash.outname=uvdata.name
	sqash.indisk=indisk
	sqash.outdisk=indisk
	sqash.inseq=0
	sqash.dparm[1]=2 #average
	sqash.bdrop=3#collapse in frequency
	sqash.dparm[2]=1 #ignore blank pixels
	sqash.outseq=0
	print 'Collapsing Image Cube for facet ' + format(uvdata.klass)
	sqash.go()



	#################



	'''


	comb = AIPSTask('COMB')
	comb.inname=uvdata.name
	comb.in2name=uvdata.name
	comb.indisk=indisk
	comb.outdisk=indisk
	comb.opcode='SUM'
	ifs=float(nifs)	#convert to a number




	#start the main loop
	for j in range(1,int(numpy.log(ifs)/numpy.log(2))+2):

		comb.inclass='IIM00' + format(j)
		comb.in2class='IIM00' + format(j)
		print "Creating combination dirty maps tier " + format(j)
		#start the inner loop
		for i in range (1,int(nifs)/(2**(j-1)),2): #from 1 to nifs/ step 2

			comb.inseq=i+1
			comb.in2seq=i
			comb.outclass='IIM00' + format(j+1)	#the output name
			comb.outseq=0
			if j==1:#i.e. the first loop so scale according to pbeam
				comb.aparm[1]=0.5#(1/float(nifs))
				comb.aparm[2]=0.5#(1/float(nifs))
			else:#not the first loop so combine beams without scaling again!
				comb.aparm[1]=0.5#1
				comb.aparm[2]=0.5#1
			comb.go()
	'''

	imageimage=WizAIPSImage(uvdata.name,'IIMNOR',indisk,1)
	imageimage.header.crval[2]=centrefreq
	imageimage.header.update()



	#now zap all but the IIMNOR map
	for i in range(nifs):
		image = AIPSImage(uvdata.name,'IIM001',indisk,i+1) #zapping the original dirty beams
		image.zap()

	image = AIPSImage(uvdata.name,'IIMCUB',indisk,1) #zapping the cube of the dirty beams
	image.zap()

#	puthead.indisk=indisk
#	puthead.inclass='IIM00'+ format(int(numpy.log(ifs)/numpy.log(2))+1)
#	puthead.inname=uvdata.name
#	puthead.inseq=1
#	puthead.keyword='CRVAL3'
#	puthead.keyvalue[1]=centrefreq
#	puthead.keyvalue[2]=0
#	puthead.go()	
	'''
	if ifs>1: #zap dirty maps - only do this if there was more than one combination tier	
		#zap the files that have already been combined

		#start the main zap loop
		for j in range(1,(int(numpy.log(ifs)/numpy.log(2))+1)):# 1 less than the total tiers

			#start the inner loop
			for i in range (1,int(nifs)/(2**(j-1))+1):
				
				image = AIPSImage(uvdata.name,'IIM00' + format(j),indisk,i)
				
				print 'zapping ' + format(uvdata.name) + "." + format('IIM00' + format(j)) + "." + format(i) + " on disk " + format(indisk)
				image.zap()
	
	
	'''
	'''
	#normalise the final dirty beam
	#divide by the total scaling factor
	sumscale=0.0
	for i in range (1,int(nifs)+1):
		sumscale=sumscale+float(scaledown[i-1])

	#now get the last file and boost it by 1/sumscale
	#use comb to do this with the same 'in and in2'
	comb.inname=uvdata.name
	comb.in2name=uvdata.name
	comb.indisk=indisk
	comb.outdisk=indisk
	comb.opcode='SUM'
	comb.inseq=1
	comb.in2seq=1
	comb.outseq=0
	comb.inclass ='IBM00'+format(int(numpy.log(ifs)/numpy.log(2))+1)
	comb.in2class ='IBM00'+format(int(numpy.log(ifs)/numpy.log(2))+1)
	comb.aparm[1]=0.5*(float(nifs))/(sumscale)
	comb.aparm[2]=0.5*(float(nifs))/(sumscale)
	comb.outclass='IBMNOR'
	
	comb.go()
	
	#now zap the original dirty beam
	image = AIPSImage(uvdata.name,'IBM00' + format(int(numpy.log(ifs)/numpy.log(2))+1),indisk,1)
	image.zap()
	'''
	'''

	#rename the dirty map
	#image = AIPSImage(uvdata.name,'IIM00' + format(int(numpy.log(ifs)/numpy.log(2))+1),indisk,1)	
	#image.rename(uvdata.name,'IIMNOR',0)
	comb.inname=uvdata.name
	comb.in2name=uvdata.name
	comb.indisk=indisk
	comb.outdisk=indisk
	comb.opcode='SUM'
	comb.inseq=1
	comb.in2seq=1
	comb.outseq=0
	comb.inclass ='IIM00'+format(int(numpy.log(ifs)/numpy.log(2))+1)
	comb.in2class ='IIM00'+format(int(numpy.log(ifs)/numpy.log(2))+1)
	comb.aparm[1]=0.5
	comb.aparm[2]=0.5
	comb.outclass='IIMNOR'
	comb.go()

	#now zap the original dirty map
	image = AIPSImage(uvdata.name,'IIM00' + format(int(numpy.log(ifs)/numpy.log(2))+1),indisk,1)
	image.zap()
	'''
def apclean(imdata,target):

	#now clean the dirty map with the dirty beam using APCLN
	apcln=AIPSTask('APCLN')
	apcln.inname=imdata.name
	apcln.in2name=imdata.name
	apcln.inclass='IIMNOR'
	apcln.in2class='IBMNOR'
	apcln.indisk=imdata.disk
	apcln.inseq=imdata.seq
	apcln.outname=imdata.name
	apcln.outclass='ICL001'
	apcln.outseq=0
	apcln.outdisk=imdata.disk
	apcln.niter=15625 #per facet
	apcln.minpatch=255
	apcln.dotv=-1
	print "Cleaning using APCLN"
	apcln.go()
	
def apclean2(imdata,target):

	#now clean the dirty map with the dirty beam using APCLN
	apcln=AIPSTask('APCLN')
	apcln.inname=imdata.name
	apcln.in2name=imdata.name
	apcln.inclass='IIM001'
	apcln.in2class='IBM001'
	apcln.indisk=imdata.disk
	apcln.inseq=imdata.seq
	apcln.outname=imdata.name
	apcln.outclass='ICL001'
	apcln.outseq=0
	apcln.outdisk=imdata.disk
	apcln.niter=0#15625 #per facet
	apcln.minpatch=255
	apcln.dotv=-1
	print "Cleaning using APCLN"
	apcln.go()

def trimmer(imdata,fovradius,target):

	# Determine the cellsize used from the image header
	cellpitch = abs(imdata.header.cdelt[0] * 3600)

	imgsize = int((fovradius*2*60/(cellpitch*8)))
	thing = imgsize
	count=2
	#print "imgsize =", imgsize, "i =", count, "radius =", fovradius, "cellpitch =", cellpitch
	while thing > 1.0 :
		thing = thing / 2
		count = count + 1
	imgsize = (2**count)#/2
	#print "imgsize =", imgsize, "i =", count, "radius =", fovradius, "cellpitch =", cellpitch

	
 	# First trim the images
	print 'Trimming images'
 	subim=AIPSTask('SUBIM')
	subsize=int(fovradius*60*2/(cellpitch*8))
	subblc=int((imgsize-subsize)/2)
	subtrc=int(imgsize-subblc)
	#print "imsize, subsize, subblc, subtrc =", imgsize, subsize, subblc, subtrc
	buffe=24 #overlap in pixels to account for projection effects
 	subim.blc[1]=subblc-buffe
	subim.blc[2]=subblc-buffe
	subim.trc[1]=subtrc+buffe
	subim.trc[2]=subtrc+buffe
	subim.xinc=1
	subim.yinc=1
	subim.indata = imdata
	suboutname='SUB' + target[-4:-1]+target[-1]
	subim.outname=suboutname
	subim.outseq=0
	subim.outdisk=imdata.disk
	if imdata.seq < 10 :
		subim.outclass='ICL00' + format(imdata.seq)
	else:
		subim.outclass='ICL0' + format(imdata.seq)
	print 'trimming facet '+ format(imdata.seq)
	subim.go()




def flat(imdata,fovradius,target):

	# Determine the cellsize used from the image header
	cellpitch = abs(imdata.header.cdelt[0] * 3600)

	imgsize = int((fovradius*2*60/(cellpitch*16)))
	thing = imgsize
	count=2
	#print "imgsize =", imgsize, "i =", count, "radius =", fovradius, "cellpitch =", cellpitch
	while thing > 1.0 :
		thing = thing / 2
		count = count + 1
	imgsize = 2**count
	#print "imgsize =", imgsize, "i =", count, "radius =", fovradius, "cellpitch =", cellpitch

	
 	# First trim the images
	print 'Trimming images'
 	subim=AIPSTask('SUBIM')
	subsize=int(fovradius*60*2/(cellpitch*16))
	subblc=int((imgsize-subsize)/2)
	subtrc=int(imgsize-subblc)
	#print "imsize, subsize, subblc, subtrc =", imgsize, subsize, subblc, subtrc
	buffe=24 #overlap in pixels to account for projection effects
 	subim.blc[1]=subblc-buffe
	subim.blc[2]=subblc-buffe
	subim.trc[1]=subtrc+buffe
	subim.trc[2]=subtrc+buffe
	subim.xinc=1
	subim.yinc=1
	subim.indata = imdata
	suboutname='SUB' + target[-4:-1]+target[-1]
	subim.outname=suboutname
	subim.outseq=0
	if imdata.seq < 10 :
		subim.outclass='ICL00' + format(imdata.seq)
	else:
		if imdata.seq < 100 :
			subim.outclass='ICL0' + format(imdata.seq)
		else:
			if imdata.seq > 99 :
				subim.outclass='ICL' + format(imdata.seq)
	print 'trimming facet '+ format(imdata.seq)
	subim.go()




def pbcorr(imdata):
 	# beam correct using pbcor
	print 'Beam correction based typical WTMOD values...'
 	pbcor=AIPSTask('PBCOR')
	pbcor.indata=imdata
	pbcor.inclass='FLATN'
	pbcor.pbparm[1]=0.3 # minimum extent of beam to believe
	pbcor.pbparm[2]=1 # use the parameters below
	pbcor.pbparm[3]=-5.679
	pbcor.pbparm[4]=168.018			# these are for the emerge field with the lo weighted at 50
	pbcor.pbparm[5]=-283.092			# eventually we want to change these dynamically, right?
	pbcor.pbparm[6]=206.940
	pbcor.pbparm[7]=-2.742
	pbcor.outclass='PBCOR'
	#pbcor.outdisk=uvdata.disk
	pbcor.outname=imdata.name
	pbcor.outseq=0
	pbcor.outdisk=imdata.disk
	pbcor.go()

def SAD(imdata):
	# Catalogue the bright sources in the field.
	print "Adjusting dimensions to create an RMS map..."
	regrd=AIPSTask('REGRD')
	regrd.imsize[1]=16384 #maximum size for RMSD to work.
	regrd.imsize[2]=16384
	regrd.inname=imdata.name
	regrd.outname=imdata.name
	regrd.inclass='PBCOR'
	regrd.outclass='REGRD'
	regrd.outseq=0
	regrd.outdisk=imdata.disk
	regrd.go()

	print "Creating RMS map..."
	rmsd=AIPSTask('RMSD')
	rmsd.imsize[1]=255 #size of area to get the rms
	rmsd.imsize[2]=255
	rmsd.inname=imdata.name
	rmsd.outname=imdata.name
	rmsd.inclass='REGRD'
	rmsd.outclass='RMSD'
	rmsd.outseq=0
	rmsd.outdisk=imdata.disk
	rmsd.go()

	print "Searching field for bright sources..."
	sad=AIPSTask('SAD')
	sad.inname=imdata.name
	sad.outname=imdata.name
	sad.inclass='REGRD' #the regridded image
	sad.in2name=imdata.name
	sad.in2class='RMSD' #the rms map
	sad.in2disk=imdata.disk
	sad.inseq=0
	sad.in2disk=imdata.disk
	sad.dparm[9]=2 #use the rms map
	sad.cparm=5 #s/n ratio
	sad.ngauss=500 #maximum number of components to search for
	sad.stvers=0 #create an ST table
	sad.go()

	tbout=AIPSTask('TBOUT')
	tbout.inname=imdata.name
	tbout.outname=imdata.name
	tbout.inclass='REGRD' #the regridded image
	tbout.inext='ST'
	tbout.invers=0
	tbout.bcount=1
	tbout.ecount=0
	tbout.outtext='PWD:WIDEFIELD_SOURCE_COORDS.TXT'
	tbout.go()

#######################################
#cillit bang peeling rountines
#######################################

def cbang_imagr (uvdata,rashift,decshift,bif,eif,cellsize,imsize,niter,nchan):
	imagr = AIPSTask('imagr')
	imagr.indata = uvdata
   	imagr.outname = 'CBANG_IF%d'%bif
	imagr.docalib = -1
	imagr.cellsize[1:] = [cellsize,cellsize]
	imagr.imsize[1:] = [imsize,imsize]
	imagr.bif = bif
	imagr.eif = eif
	imagr.nfield = len(rashift)
	for i in range(len(rashift)):
		imagr.rashift[i+1] = float(rashift[i])
		imagr.decshift[i+1] = float(decshift[i])
	imagr.dotv = -1
	imagr.niter = niter
	imagr.nchav = nchan
	imagr.go()

def cbang_ccmrg (image):
	ccmrg = AIPSTask('ccmrg')
	ccmrg.indata = image
	ccmrg.go()

def cbang_dbcon (inna,incl,in2na,in2cl,outna,outcl):
	dbcon = AIPSTask('dbcon')
	dbcon.inname = inna
	dbcon.inclass = incl
	dbcon.in2name = in2na
	dbcon.in2class = in2cl
	dbcon.outname = outna
	dbcon.outclass = outcl
	dbcon.dopos[1][1:] = [-1,-1]
	dbcon.dopos[2][1:] = [-1,-1]
	dbcon.doarray = 1
	dbcon.go()

def cbang_uvsub (ifno, nmaps, disk):
	uvsub = AIPSTask('uvsub')
	uvsub.inname = 'CBANG_IF%d' % ifno
	uvsub.inclass = 'UVCOP'
	uvsub.in2name = 'CBANG_IF%d' % ifno
	uvsub.in2cl = 'ICL001'
	uvsub.indisk = disk
	uvsub.nmaps = nmaps
	uvsub.outname = 'CBANG_IF%d' % ifno
	uvsub.outcl = 'UVSUB'
	uvsub.outdisk = disk
	uvsub.inp()
	uvsub.go()

def cbang_uvcop (uvdata, i):
	uvcop = AIPSTask('uvcop')
	uvcop.indata = uvdata
	uvcop.outcl = 'UVCOP'
	uvcop.uvcopprm[1:] = [1,0,0]
	uvcop.outname = 'CBANG_IF%d' % i
	uvcop.bif = i
	uvcop.eif = i
	uvcop.go()

def cbang_getsrc (ra,dec):
	first, isfirst = np.load('first_2008.simple.npy'), True
	corr = correlate(np.array([[ra,dec]]),0,1,first,0,1,0.2)
	if not corr:
		first, isfirst = np.load('NVSS.npy'), False
		corr = correlate(np.array([[ra,dec]]),0,1,first,0,2,0.2)
	sra,sdec,sdist,sflux = [],[],[],[]
	for i in corr:
		sra.append(first[int(i[1]),0])
		sdec.append(first[int(i[1]),1 if isfirst else 2])
		sflux.append(first[int(i[1]),3 if isfirst else 5])
		sdist.append(i[2])
	print 'Central RA, dec:',astCoords.decimal2hms(ra,' '),astCoords.decimal2dms(dec,' ')
	print 'Processing the following sources:\n  RA         DEC   flux   distance/deg'
	for i in range(len(sra)):
		print '%s %s %7.2f %7.2f'%(astCoords.decimal2hms(sra[i],' '),\
		astCoords.decimal2dms(sdec[i],' '),sflux[i],sdist[i])
	src = np.column_stack((sra,sdec,sflux,sdist))
	src = np.column_stack((src,3600.*np.cos(np.deg2rad(dec))*(src[:,0]-ra)))
	src = np.column_stack((src,3600.*(src[:,1]-dec)))
	return src, isfirst

def cbang_parse (ra,dec,infile):
	a = np.loadtxt(infile,dtype='S')
	src = np.zeros((len(a),6))
	src[:,2] = 1.0
	for i in range(len(a)):
		src[i,0] = astCoords.hms2decimal(a[i,0]+' '+a[i,1]+' '+a[i,2],' ')
		src[i,1] = astCoords.dms2decimal(a[i,3]+' '+a[i,4]+' '+a[i,5],' ')
		src[i,4] = 3600.*np.cos(np.deg2rad(dec))*(src[i,0]-ra)
		src[i,5] = 3600.*(src[i,1]-dec)
		src[i,3] = np.hypot(src[i,4],src[i,5])
	print 'Central RA, dec:',astCoords.decimal2hms(ra,' '),astCoords.decimal2dms(dec,' ')
	print 'Processing the following sources from file %s:\n  RA         DEC   flux   distance/deg'%infile
	for i in range(len(src)):
		print '%s %s %7.2f %7.2f %7.2f'%(astCoords.decimal2hms(src[i,0],' '),\
 		astCoords.decimal2dms(src[i,1],' '),src[i,3]/3600.,src[i,4]/3600.,src[i,5]/3600.)
	return src, True

def cbang_clean (n_if,nfield):
	for i in range(1,1+n_if):
		uvdata = AIPSUVData('CBANG_IF%d'%i,'UVSUB',1,1)
		try:
			uvdata.zap()
		except:
			pass
		uvdata = AIPSUVData('CBANG_IF%d'%i,'UVCOP',1,1)
		try:
	  		uvdata.zap()
		except:
			 pass
		if i<n_if-1:
			uvdata = AIPSUVData('CBANG','PROC%d'%i,1,1)
			uvdata.zap()
		for j in range(1,1+nfield):
			imdata = AIPSImage('CBANG_IF%d'%i,'ICL%03d'%j,1,1)
			try:
				imdata.zap()
			except:
				pass
			imdata = AIPSImage('CBANG_IF%d'%i,'IBM%03d'%j,1,1)
			try:
				imdata.zap()
			except:
				pass


def cbang (inna,incl,indi,inseq,cellsize,imsize=128,niter=100,zap=True,infile=''):
	uvdata = AIPSUVData (inna,incl,indi,inseq)
	ra,dec = uvdata.header['obsra'], uvdata.header['obsdec']
	ctypes=np.asarray(uvdata.header['ctype'],dtype='S')
	n_if = uvdata.header['naxis'][np.argwhere(ctypes=='IF')[0][0]]
	nchan = uvdata.header['naxis'][np.argwhere(ctypes=='FREQ')[0][0]]
	print 'Found ',n_if,' IFs each with ',nchan,' channels'
	if infile:
		src, isfirst = cbang_parse(ra,dec,infile)
	else:
		src, isfirst = cbang_getsrc(ra,dec)
	f = np.zeros((len(src),n_if))
	if not isfirst:
 		imsize = 512       # only NVSS, positions less accurate
	for i in range(1,1+n_if):
  		cbang_uvcop (uvdata, i)
		try:
			cbang_imagr (uvdata,src[:,-2],src[:,-1],i,i,\
				cellsize,imsize,niter,nchan)
			for j in range(1,1+len(src)):
				imdata = AIPSImage ('CBANG_IF%d'%i,'ICL%03d'%j,indi,inseq)
				cbang_ccmrg(imdata)
				imtab = imdata.table('CC',2)
				for k in imtab:
					if k['flux']<0.0:
						break
					f[j-1,i-1]+=k['flux']
			cbang_uvsub (i,len(src),indi)
		except:
			print 'No data for IF:',i,', continuing'
			tmpdata = AIPSUVData('CBANG_IF%d'%i,'UVCOP',1,1)
			tmpdata.rename ('CBANG_IF%d'%i, 'UVSUB')

	cbang_dbcon ('CBANG_IF1','UVSUB','CBANG_IF2','UVSUB','CBANG','PROC1')
	for i in range(n_if-2):
		d_out = inna if i==n_if-3 else 'CBANG'
		d_cl = 'UV' + incl[2:] if i==n_if-3 else 'PROC%d'%(i+2)
		cbang_dbcon ('CBANG','PROC%d'%(i+1),'CBANG_IF%d'%(i+3),'UVSUB',d_out,d_cl)

	if zap:
		cbang_clean (n_if, len(src))


# end of cillitbang peeler

# correlate two arrays sorted in ra. array2 is the bigger one.

def correlate (array1, ra1, dec1, array2, ra2, dec2, dist, \
               mindist=0.0, isabs=False, noisy=True):
    fstart=nfstart=0
    fend=array2.shape[0]
    icou=0
    correl=np.array([])
    decfac=1.0 if isabs else min(1./np.cos(np.deg2rad(array1[:,dec1].max())),\
               1./np.cos(np.deg2rad(array2[:,dec2].max())))
    
    for i in range(array1.shape[0]):
        i10 = np.linspace(0,array1.shape[0],10,dtype='int')
        i100 = np.linspace(0,array1.shape[0],100,dtype='int')
        if i in i10 and noisy:
            sys.stdout.write('*')
            sys.stdout.flush()
        elif i in i100 and noisy:
            sys.stdout.write('.')
            sys.stdout.flush()
        else:
            pass
        fstart=nfstart
        for j in range(fstart,fend):
            radiff=array2[j,ra2]-array1[i,ra1]
            if radiff<-decfac*dist:
                nfstart=j
            if radiff>decfac*dist:
                break
            if abs(array2[j,dec2]-array1[i,dec1])>dist:
                continue
            adist = np.hypot(array1[i,ra1]-array2[j,ra2],\
                             array1[i,dec1]-array2[j,dec2]) if isabs \
              else astCoords.calcAngSepDeg(array1[i,ra1],array1[i,dec1],\
                                          array2[j,ra2],array2[j,dec2])
            if adist<dist and abs(radiff)<90.0 and adist>=mindist:
                try:
                    correl=np.vstack((correl,np.array([i,j,adist])))
                except:
                    correl=np.array([[i,j,adist]])

    return correl
            



########################################################################################################################################
# QA plotting functions
########################################################################################################################################

def plot(data, plotname, dopng=False):
	# standard plotting routine

	# Rename any old plot files before writing the new one
	#save_old_file(plotname)		
	##print >>FuncLog, 'Creating plot file: ', plotname
	delfile(plotname)
	runlwpla2(indata=data, outfile=plotname)

	#plot_base = os.path.splitext(plotname)[0] 
	##print >>FuncLog, 'Converting ', plotname, ' to pdf.'
	##if dopng:
		##print >>FuncLog, 'Converting ', plotname, ' to png.'

	# convert from PS. Conversion can be quite slow so do in parallel while the
	# pipeline continues.
	global pid_list
	pid_list = pid_cleanup(pid_list)
	pid = os.fork()
	if pid == 0:
		plot_convert(plotname, dopng)
		os._exit(0)

	pid_list.append(pid)
	#print >>FuncLog, 'pid_list = ', pid_list

def pid_cleanup(pid_list):
	# clean up child processes, and return list of child processes still
	# running
	pid_notdone = []
	for pid in pid_list:
		#print >>FuncLog, 'pid=', pid
		wait_pid, status = os.waitpid(pid, os.WNOHANG)
		#print >>FuncLog, 'wait_pid, status=', wait_pid, status
		if wait_pid != pid:
			pid_notdone.append(pid)

	return pid_notdone

def plot_convert(plotname, dopng):
	# convert the postscript files to pdf (and png in some cases). Could use
	# 'convert', but ps2pdf seems to produce smaller files with the same
	# quality.	 

	# Convert to pdf
	print 'Converting ...', plotname
	pdfname = os.path.splitext(plotname)[0] + '.pdf'
	##save_old_file(pdfname)
	os.system('ps2pdf  -sPAPERSIZE=a4 ' +  plotname + ' ' + pdfname)
	#pid = os.spawnv(os.P_NOWAIT, 'ps2pdf', ['plotname', 'pdfname'])

	# Convert to png if requested. 
	if dopng:
		pngname = os.path.splitext(plotname)[0] + '_%03d.png'
		##save_old_file(pngname)
		#print >>FuncLog, 'Converting plot file to png: ', pngname
		os.system('gs -dNOPAUSE -dQUIET -sDEVICE=png16m -r80 -dTextAlphaBits=4 -dGraphicsAlphaBits=4  -sOutputFile={1} -dBATCH {0}'.format(pdfname, pngname))
		#os.system('convert ' + plotname + ' ' + pngname)
		#print 'Finished converting to png: ', pngname
		#pid = os.spawnv(os.P_NOWAIT, 'convert', ['plotname', 'pngname'])

	# delete the postscript
	##print 'Deleting the postscript file: ', plotname
	print 'Finish converting plot', plotname
	delfile(plotname)

def delfile(file_to_del):
	#print 'Deleting...', file_to_del
	if os.path.isfile(file_to_del):
		os.system('rm -r '+file_to_del)

def filesize(file):
	""" Returns file size in MB"""
	try:
		filesize = os.stat(file).st_size/1e6
	except:
		filesize = 0.0
	return filesize

def set_default_aparms(task):
	# simple function to set default aparms for various tasks. Returns default
	# values as a list, which can then be changed before calling the function
	# that runs the task. Note that the bparm elements match those in AIPS - do
	# not pass element zero to the tasks (e.g. use possm.aparm[1:] = aparm[1:])
	aparm = []
	aparm[1:] = [0 for i in range(10)]
	if task == 'possm':
		aparm[1] = 0
		aparm[2] = 1
		aparm[5] = -180
		aparm[6] = 180
		aparm[9] = 1  # IF together
		#aparm[9] = 3  # IF and polarizations together


	if task == 'fring':
		aparm[1] = 3
		aparm[6] = 2
		# aparm[7] (the snr cutoff) is set in the call to fring
		#aparm[7] = 11
		aparm[9] = 1

	return aparm

def set_default_bparms(task):
	# cf aparms above
	bparm = []
	bparm[1:] = [0 for i in range(10)]
	if task == 'vplot':
		bparm[2] = -1
		bparm[3] = 1
		bparm[8] = -180
		bparm[9] = 180

	return bparm

def runlwpla2(indata, outfile, inver=None, plver=1):

	'''Must set indata, outfile'''
	#assert (indata != None, outfile != None)

	lwpla = AIPSTask('lwpla')
	lwpla.indata = indata
	lwpla.outfile = outfile 
	lwpla.lpen = 1.
	lwpla.dparm[5] = 2
	lwpla.dparm[6] = 4
	lwpla.dparm[7] = 31  # Helvetica
	lwpla.dparm[8] = 10   # Font size
	lwpla.plver = plver
	if (inver == None):
		inver = indata.table_highver('AIPS PL')
	lwpla.inver = inver
	#lwpla.inputs()
	lwpla()

def runlistr(indata, outprint, optype='', inver=0, freqid=1,
		docalib=0, inext='SN', stokes='HALF', sources=[], dparm=[]):

	"""Must set indata, outprint"""
	#assert (indata != None, outprint != None)

	listr = AIPSTask('listr')
	listr.indata = indata
	listr.outprint = outprint
	listr.optype = optype
	listr.inext = inext
	listr.inver = inver
	listr.stokes = stokes
	listr.freqid = freqid
	listr.sources[1:] = sources
	listr.bif = 0
	listr.eif = 0
	listr.bchan = 1
	listr.echan = 0
	listr.docalib = docalib
	listr.gainuse = 0
	listr.flagver = 0
	listr.dparm[1:] = dparm[1:]
	listr.dparm[3] = 6
	listr.docrt = -3
	listr()

def rundtsum(indata, outprint=''):
	'''Must set indata'''

	#assert (indata != None)
	dtsum = AIPSTask('dtsum')
	dtsum.indata = indata
	dtsum.aparm[1] = 1.0
	dtsum.outprint = outprint
	dtsum.docrt = -3
	dtsum()

def runpossm2(indata, aparm, solint=-1, freqid=1, sources=[],
		timerang=[0,0,0,1], codetype='A&P',
		nplots=9, stokes='HALF', antennas=[], docalib=0, gainuse=0,
		doband=0, bpver=0, flagver=1, dotv=-1, factor = 0):

	"""Must set indata, aparm"""

	#assert (aparm != None, indata != None), "set indata, aparm in runpossm!"

	#print 'possm antennas=', antennas
	#raw_input('press return')

	possm = AIPSTask('possm')

	possm.stokes = stokes
	possm.indata = indata
	possm.solint = solint
	possm.freqid = freqid
	possm.sources[1:] = sources
	possm.timerang[1:] = timerang
	possm.antennas[1:] = antennas
	possm.docalib = docalib
	possm.gainuse = gainuse
	possm.flagver = flagver
	possm.doband = doband
	possm.bpver = bpver
	possm.aparm[1:] = aparm[1:]
	possm.codetype = codetype
	possm.nplots = nplots
	possm.factor = factor
	possm.dotv = dotv
	#possm.inputs()
	possm()


def is_aipsdata(aipsdata):
	'''Check whether the passed object has the valid attributes for an AIPS
	data object'''

	got_attr = False
	if (hasattr(aipsdata, 'name') and hasattr(aipsdata, 'disk') and
			hasattr(aipsdata, 'seq') and hasattr(aipsdata, 'klass') ):
		got_attr = True

	return got_attr

def runvplot2(indata, bparm, antenna=[], in2data=None, freqid=1, bchan=1, 
		echan=4096, bif=1, eif=0, sources=[], docalib=1, gainuse=0, flagver=0,
		doband=-1, bpver=0, solint=0, nplots=4, ncomp=0, nmaps=0, stokes='HALF',
		dotv=-1, do3col = 0, symbol = 0, factor = 0, xyratio = 0):

	"""Must set indata, bparm"""
	#assert (indata != None, bparm != None  and antenna !=
	#		None), """bparms or antenna not set in runvplot!"""
	vplot = AIPSTask('vplot')
	vplot.indata = indata
	if is_aipsdata(in2data):
		vplot.in2data = in2data
	vplot.ncomp[1] = ncomp
	vplot.nmaps = nmaps
	vplot.sources[1:] = sources
	vplot.freqid = freqid 
	vplot.bchan = bchan
	vplot.echan = echan 
	vplot.bif = bif
	vplot.eif = eif 
	vplot.stokes = stokes
	vplot.antenna[1:] = antenna
	vplot.docalib = docalib
	vplot.gainuse = gainuse
	vplot.flagver = flagver
	vplot.doband = doband
	vplot.doebar = -1
	vplot.solint = solint 
	vplot.bparm[1:] = bparm[1:]
	vplot.nplots = nplots 
	vplot.ltype = 3 
	vplot.dotv = dotv
	vplot.grchan = 1
	vplot.do3col = do3col
	vplot.symbol = symbol
	vplot.factor = factor
	vplot.xyratio = xyratio
	#vplot.inputs()
	error = 0
	try:
		vplot()
		error = 0
	except:
		error = 1
	return error


def runuvplt2(indata, bparm, freqid=1, docalib=2, gainuse=0, sources=[], 
		dotv=-1, antennas=[], baseline=[], grchan=0, bchan=1, echan=0, nchav = 0, refant=0, xinc = 1, do3col = 0, symbol = 0, factor = 0, stokes = 'HALF'):

	"""Must set indata, bparm"""
	#assert (indata != None, bparm != None), 'set indata, bparm in runuvplt'

	uvplt = AIPSTask('uvplt')
	uvplt.indata = indata
	uvplt.sources[1:] = sources
	uvplt.qual = -1
	uvplt.calcod = ''
	uvplt.stokes = stokes
	uvplt.selband = -1
	uvplt.selfreq = -1
	uvplt.freqid = freqid
	uvplt.timerang[1:] = [0]
	uvplt.anten[1:] = antennas
	uvplt.basel[1:] = baseline
	uvplt.uvrange[1:] = [0]
	uvplt.subarray = 0
	uvplt.bchan = bchan
	uvplt.echan = echan
	uvplt.nchav = nchav
	uvplt.chinc = 1
	uvplt.bif = 1
	uvplt.eif = 0
	uvplt.docalib = docalib
	uvplt.gainuse = gainuse
	uvplt.dopol = -1
	uvplt.blver = -1
	uvplt.flagver = 0
	uvplt.doband = -1
	uvplt.bpver = 0
	uvplt.smooth[1:] = [0]
	uvplt.xinc = xinc
	uvplt.bparm[1:] = bparm[1:]
	uvplt.doweight = 1
	uvplt.refant = refant
	uvplt.rotate = 0
	uvplt.factor = 0
	uvplt.do3col = do3col
	uvplt.ltype = 3
	uvplt.dotv = dotv
	uvplt.grchan = grchan

	uvplt()

def runfittp2(indata, outfile):

	"""Must set indata, outfile"""
	#assert (outfile != None), 'set outfile in runfittp'

	fittp = AIPSTask('fittp')
	fittp.indata = indata
	fittp.doall = -1
	fittp.intype = ''
	fittp.outtape = 1
	try:
		fittp.dataout = outfile
	except AttributeError:
		fittp.outfile = outfile
	fittp()

def save_old_file(filename, save_old=False):
	# check if an external file exists and delete/rename if it does
	if os.path.isfile(filename):
		if save_old:
			oldfilename = filename
			uniq = 0
			while os.path.isfile(oldfilename):
				oldfilename = filename + '.' + str(uniq)
				uniq += 1
			os.rename(filename, oldfilename)
			output = "\nrenaming " + str(filename) + " to " + str(oldfilename)
		else:
			os.remove(filename)
			output = "\ndeleting " + str(filename) 
		#print >>FuncLog, output

def write_spplot_input(outname, conf_file, uvdata, scale, amporphas, spplot_dir = './plots/', baselines = [], sources_spplot = ''):
	if len(numpy.atleast_1d(sources_spplot)) == 1:
		if sources_spplot == '':
			choose = 'all'
			src_choose = ''
		else:
			choose = 'choose'
			src_choose = sources_spplot
	else:
		choose = 'choose'
		src_choose = "', '".join(sources_spplot)
	f = open(conf_file, 'wb')
	f.write("AIPS.userno = {0}\n".format(AIPS.userno))
	f.write("Name = '{0}'\n".format(uvdata.name))
	f.write("Klass = '{0}'\n".format(uvdata.klass))
	f.write("Disk = {0}\n".format(uvdata.disk))
	f.write("Seq = {0}\n".format(uvdata.seq))
	f.write("path2folder = '{0}'\n".format(spplot_dir))
	f.write("picklepath = '{0}'\n".format(spplot_dir+'picklefiles/'))
	f.write("choosesources = '{0}'\n".format(choose))
	f.write("specifysources = ['{0}']\n".format(src_choose))
	f.write("choosebaselines = '{0}'\n".format('all'))
	f.write("baselines = {0}\n".format(baselines))
	f.write("outfilename = '{0}'\n".format(outname))
	f.write("stokes = [{0}]\n".format("'RR','LL'"))
	f.write("scale = '{0}'\n".format(scale))
	f.write("scale_over_all_IFs = '{0}'\n".format('yes'))
	f.write("colourscheme = '{0}'\n".format('CMRmap'))
	f.write("amporphas = '{0}'\n".format(amporphas)) #'A' or 'P'
	f.write("timeperpage = {0}\n".format(100000))
	f.write("IF = {0}\n".format(nif))
	f.write("IF_start = {0}\n".format(1))
	f.write("IF_end = {0}\n".format(nif))
	f.close()
 
def run_spplot(conf_file, spplot_dir, sources_spplot):
	try:
		os.mkdir(spplot_dir)
	except:
		pass
	if len(numpy.atleast_1d(sources_spplot)) == 1:
		if sources_spplot == '':
			outname = 'all'
		else:
			outname = sources_spplot
	else:
		outname = 'diverse'
	print 'SP Running SPPLOT: ', 'parseltongue '+spplot_path+' '+ conf_file
	os.system('parseltongue '+spplot_path+' '+ conf_file)
	print 'SP Moving spplot pdf: ', 'mv '+spplot_dir+project+'.SPPLOT.*.pdf '+ plotdir
	os.system('mv '+spplot_dir+project+'.SPPLOT.*.pdf '+ plotdir)
	print 'SP Converting 1:', 'gs -dNOPAUSE -dQUIET -sDEVICE=png16m -r150 -dTextAlphaBits=4 -dGraphicsAlphaBits=4  -sOutputFile={1} -dBATCH {0}'.format(plotdir+project+'.SPPLOT.AMP.UNCAL_'+outname+'.combined.pdf', plotdir+project+'.SPPLOT.AMP.UNCAL_'+outname+'.combined.png')
	os.system('gs -dNOPAUSE -dQUIET -sDEVICE=png16m -r150 -dTextAlphaBits=4 -dGraphicsAlphaBits=4  -sOutputFile={1} -dBATCH {0}'.format(plotdir+project+'.SPPLOT.AMP.UNCAL_'+outname+'.combined.pdf', plotdir+project+'.SPPLOT.AMP.UNCAL_'+outname+'.combined.png'))
	print 'SP Converting 2: ', 'gs -dNOPAUSE -dQUIET -sDEVICE=png16m -r150 -dTextAlphaBits=4 -dGraphicsAlphaBits=4  -sOutputFile={1} -dBATCH {0}'.format(plotdir+project+'.SPPLOT.PHS.UNCAL_'+outname+'.combined.pdf', plotdir+project+'.SPPLOT.PHS.UNCAL_'+outname+'.combined.png')
	os.system('gs -dNOPAUSE -dQUIET -sDEVICE=png16m -r150 -dTextAlphaBits=4 -dGraphicsAlphaBits=4  -sOutputFile={1} -dBATCH {0}'.format(plotdir+project+'.SPPLOT.PHS.UNCAL_'+outname+'.combined.pdf', plotdir+project+'.SPPLOT.PHS.UNCAL_'+outname+'.combined.png'))
	print 'SP Removing: ', 'rm -r '+spplot_dir
	os.system('rm -r '+spplot_dir)


def write_linksize(filetype, infile, intext, explanation = '', inpng=False):
	if filetype == 'txt':
		wlog.write('{0}: [<a href="{1}" target="_blank">txt</a>]<br>\n'.format(intext, infile))
	elif filetype == 'fits':
		wlog.write('<tr><td valign="top">{0}</td><td valign="top">[<a href="{1}" target="_blank">fits</a>] <small>{2:5.1f} GB</small></td><td valign="top"><small>{3}</small></td>\n'.format(intext, infile, filesize(infile[1:])/1000., explanation))
	elif filetype == 'pdf':
		if inpng == False:
			wlog.write('<tr><td valign="top">{0}</td><td valign="top">[<a href="{1}" target="_blank">pdf</a>] <small>{2:5.1f} MB</small></td><td valign="top"><small>{3}</small></td>\n'.format(intext, infile, filesize(infile[1:]), explanation))
		else:
			wlog.write('<tr><td valign="top">{0}</td><td valign="top">[<a href="{1}" target="_blank">pdf</a>] <small>{2:5.1f} MB</small> [<a href="{3}" target="_blank">png</a>] <small>{4:5.1f} MB</small></td><td valign="top"><small>{5}</small></td>\n'.format(intext, infile, filesize(infile[1:]), inpng, filesize(inpng[1:]), explanation))
	else:
		wlog.write('{0}: [<a href="{1}" target="_blank">file</a>] {2}<br>\n'.format(intext, infile, explanation)) 

def weblog_header(wlog, section):
	wlog.write('<html>\n')
	wlog.write('<head>\n')
	wlog.write('<title>eMERLIN Pipeline Web Log</title>\n')
	wlog.write('</head>\n')
	wlog.write('<body>\n')
	wlog.write('<h5>eMERLIN Pipeline Web Log</h5>\n')
	wlog.write('<form>\n')
	wlog.write('<input type="button" value="Home"  onclick="window.location.href=\'./index.html\'">\n')
	wlog.write('<input type="button" value="Observation summary"  onclick="window.location.href=\'./obs_summary.html\'">\n')
	wlog.write('<input type="button" value="Plots"  onclick="window.location.href=\'./plots.html\'">\n')
	wlog.write('<input type="button" value="Files"  onclick="window.location.href=\'./files.html\'">\n')
	wlog.write('</form>\n')
	#wlog.write('<br>\n')
	wlog.write('<h1>{0}</h1>\n'.format(section))
	wlog.write('<hr>\n')

def weblog_foot(wlog):
	wlog.write('<hr><br>\n')
	wlog.write('<small><a href="http://www.e-merlin.ac.uk/observe/pipeline/">e-MERLIN Pipeline</a><br>\n')
	wlog.write('<a href="https://github.com/mkargo/pipeline">e-MERLIN Pipeline Github</a><br>\n')
	wlog.write('<a href="http://www.e-merlin.ac.uk/">e-MERLIN</a><br>\n')
	wlog.write('<a href="http://www.e-merlin.ac.uk/data_red/">User support and data reduction for e-MERLIN</a><br>\n')
	wlog.write('<a href="http://www.jive.nl/jivewiki/doku.php?id=parseltongue:parseltongue">Parseltongue</a></small><br>\n')
	wlog.write('<hr>\n')
	wlog.write('</body>\n')
	wlog.write('</html>\n')


print "Specific e-MERLIN procedures loaded! Running pipeline."
