import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.axes_grid1.inset_locator import BboxPatch, BboxConnector, mark_inset, zoomed_inset_axes
from matplotlib.transforms import Bbox, TransformedBbox, IdentityTransform
from matplotlib.patches import Rectangle
from matplotlib.path import Path
import matplotlib.patches as patches
import matplotlib
import math
import numpy as np

#as recorded in original fits file
delta_t = 0.00102399999741

#MCMC estimated dm
dm = 622.8

def run(base_path= '/Users/alex/burst_data/', data_file='filtered_short',aspct=1.0,
	f_rebin=0,t_rebin=0,ave_bins=(1,1),dist=False,zoom=2.0,std_show=[2.0,2.0],
	x_anch=1.0,band_sample=False,rad=75,lineshade=0.5,save=False,dpi=300.0):
	#set mpl params
	matplotlib.rcParams.update({
	    'pdf.fonttype': 42,
	    'ps.fonttype': 42,
	    'font.serif' : ['Times New Roman'],
	    'font.sans-serif' : ['Helvetica', 'Arial'],
	    })
	matplotlib.rcParams['savefig.dpi'] = dpi
	ft = np.load(base_path + data_file + '.npy')

	#select intensity data
	ft = ft[:,0,:] # + ft[:,-1,:]

	#there is a small error here;
	#t_len should be 2400, but it
	#is not used anywhere substantial,
	#only to place the zoomed plot
	#region on the full data. Correcting
	#this error would only change the
	#zoomed region scale factors, but to
	#ensure that plotting is repeatable, I
	#will leave this as-is.
	t_len = 3000
	f_len = 4096
	spec_tmin , spec_tmax = 400*delta_t, 2800*delta_t

	cs = ft[:,spec_tmin/delta_t:spec_tmax/delta_t]

	#convert NaN to zero
	cs = np.nan_to_num(cs)
	gs = cm.gray_r

	dat = cs[:,:]

	#mean subtract each frequency channel
	mean_sub(dat)

	#diagnostic plots
	if(dist):
		plt.figure()
		hist, bedges = np.histogram(dat,bins=1000)
		plt.plot(bedges[0: len(bedges) - 1], hist)

		plt.figure()
		for p in integrated_power(dat):
			plt.plot(p[0,:])

	print 'original shape', dat.shape

	#unbinned should be called 'less binned'
	#it is used in the zoomed inset to provide
	#additional resolution
	unbinned = simple_rebin(dat,f_rebin/2,t_rebin/2)
	dat = simple_rebin(dat,f_rebin,t_rebin)
	normalize(unbinned)

	ub_min = np.mean(unbinned) - std_show[0]*np.std(unbinned)
	ub_max = np.mean(unbinned) + std_show[1]*np.std(unbinned)

	print 'inset min/max', ub_min, ub_max

	#apply computed statistical cuts
	unbinned = np.clip(unbinned,ub_min,ub_max)

	#compute mean and stdev of data
	#prior to clipping
	disp_mean = np.mean(dat)
	disp_stdev = np.std(dat)
	mn, mx = normalize(dat)
	mean = np.mean(dat)
	stdev = np.std(dat)

	cmin = mean - std_show[0]*stdev
	cmax = mean + std_show[1]*stdev

	print 'main panel min/max', cmin, cmax

	#apply computed cuts. Note that the
	#range (in stdev) is the same as used for
	#the inset data (unbinned)
	dat = np.clip(dat,cmin,cmax)

	dat_min = cmin
	dat_max = cmax
	dat_mean = mean

	cmin = disp_mean - std_show[0]*disp_stdev
	cmax = disp_mean + std_show[1]*disp_stdev

	#more diagnostic plots
	if band_sample:
		plt.figure()
		samples = 3
		row_inc = int(math.floor(float(dat.shape[0])/5.0))
		for i in xrange(0,samples):
			plt.plot(dat[row_inc*i - 1,:])

	if dist:
		plt.figure()
		hist, bedges = np.histogram(dat,bins=1000)
		plt.plot(bedges[0: len(bedges) - 1], hist)

	print cmin
	print cmax

	fig, ax = plt.subplots(figsize=(12,9))

	n_f = dat.shape[0]
	n_t = dat.shape[1]

	x = np.array(range(0,2400))

	#guides set to the approximate DM and displaced in time
	right_guide = n_f*(900.0 - np.sqrt(4.148808*dm*1000/(delta_t*x + 3.15099677037037)))/200.0
	left_guide = n_f*(900.0 - np.sqrt(4.148808*dm*1000/(delta_t*x + 3.19099677037037)))/200.0

	center = 158
	ax.plot((x + center + rad)/float(t_rebin),right_guide, color='1.0',linewidth=2)
	ax.plot((x + center - rad)/float(t_rebin),left_guide, color='1.0',linewidth=2)


	cax = ax.imshow(dat,interpolation='nearest',cmap=gs,aspect=aspct*float(dat.shape[1])/float(dat.shape[0]))

	print dat_min, dat_mean, dat_max

	cbar = fig.colorbar(cax, ticks=[np.min(dat), np.mean(dat), np.max(dat)])
	cbar.ax.set_yticklabels(['%.1f' % round(cmin,1), '0.0', '%.1f' % round(cmax,1)],fontsize=24)# vertically oriented colorbar
	cbar.set_label('Antenna Temperature (K)',fontsize=24)
	plt.ylabel('Frequency (MHz)', fontsize=24)
	plt.xlabel('Time (s)', fontsize=24)

	#since time is arbitrary,
	#choose offsets so rounding
	#is exact
	del_t = spec_tmax - spec_tmin
	dd_t = del_t%0.1
	del_t = del_t - dd_t

	delta_shape = dd_t/(t_rebin*delta_t)

	xlocs = np.linspace(0, dat.shape[1] - delta_shape ,5)
	xlabs = np.linspace(0, del_t, 5)
	xlabs = ['%.1f' % a for a in xlabs]
	ylocs = np.linspace(0, dat.shape[0], 3)
	ylabs = np.linspace(900, 700, 3)
	ylabs = [int(a) for a in ylabs]
	plt.xticks(xlocs,xlabs,fontsize=24,horizontalalignment='left')
	plt.yticks(ylocs,ylabs,fontsize=24)


	shape = dat.shape

	#for use with binned data
	orig_shape = f_len/(f_rebin/2), t_len/(t_rebin/2)

	f_len = f_len/math.pow(f_rebin,1.0)
	t_len = t_len/math.pow(t_rebin,1.0)

	#arbitrary scales chosen
	#to isolate a high-intensity
	#subregion
	x1,x2,y2,y1 = 0.465*t_len, 0.63*t_len, 0.87*f_len, 0.69*f_len

	lims = 0.465*orig_shape[1], 0.63*orig_shape[1], 0.87*orig_shape[0], 0.69*orig_shape[0]

	verts_l = [(x1,y1),
		(x2,y1)]
	verts_r = [(x1,y1),
		(x2,y1)]

	path = Path(verts_l)
	
	axins = zoomed_inset_axes(ax, zoom, loc=1, bbox_to_anchor=(x_anch,1.0), bbox_transform=ax.transAxes )

	plt.xticks(visible=False)
	plt.yticks(visible=False)

	axins.imshow(unbinned, interpolation='nearest', cmap=gs,aspect=aspct*float(dat.shape[1])/float(dat.shape[0]))
	axins.set_xlim(lims[0],lims[1])
	axins.set_ylim(lims[2],lims[3])

	print (x1,x2)
	print axins.get_xlim()

	#run the custom marking routine to
	#draw the box and outlines
	cust_mark_inset(ax, axins, x1,x2,y2,y1,loc11=3,loc12=3,loc21=4, loc22=4, fc="none", ec='1.0',linewidth=2)
	
	#save outputs
	if save:
		plt.savefig('ft_inset.png')
		plt.savefig('ft_inset.pdf')
		plt.savefig('ft_inset.eps')

	plt.show()

def integrated_power(dat,dm=dm,dt = 0.00102399999741,bins=10):
	inc = float(dat.shape[0] - 1)/bins
	ret = []
	dat,delay_ind = dedisperse(dat,dm,dt)
	for i in xrange(0,bins - 1):
		subdat = dat[i*inc:(i + 1)*inc,:]
		elem = np.zeros((1,subdat.shape[1]))
		for i in xrange(0,elem.shape[1]):
			elem[0,i] = np.sum(subdat[:,i])
		ret.append(elem)
	return ret

def mov_ave(dat, f_bins, t_bins):
	ret = np.zeros(dat.shape - np.array((f_bins - 1,t_bins - 1)))

	for i in xrange(0,ret.shape[0]):
		for j in xrange(0,ret.shape[1]):
			ret[i,j] = np.mean(dat[i:i + f_bins - 1,j:j + t_bins - 1])

		print 'moving average row {0} complete'.format(i)
	return ret

def dedisperse(dat,dm=dm,dt = 0.00102399999741):
	ret = np.zeros(dat.shape)
	df = -200.0/float(dat.shape[0])
	f0 = 900
	for i in xrange(0,dat.shape[0]):
		f = f0 + i*df
		delay_ind = int(round((disp_delay(f,dm) - disp_delay(f0,dm))/dt))
		for j in xrange(0,dat.shape[1] - delay_ind):
			ret[i,j] = dat[i,j + delay_ind]
	return ret,int((disp_delay(700,dm) - disp_delay(900,dm))/dt)

def disperse(dat,dm=dm):
	ret = np.zeros(dat.shape)
	dt = 0.00102399999741
	df = -200.0/float(dat.shape[0])
	f0 = 900
	for i in xrange(0,dat.shape[0]):
		f = f0 + i*df
		delay_ind = int(round((disp_delay(f,dm) - disp_delay(f0,dm))/dt))
		for j in xrange(0,dat.shape[1] - delay_ind - 1):
			ret[i,ret.shape[1] - 1 - j] = dat[i,dat.shape[1] - 1 -  j - delay_ind]
	return ret

def disp_delay(f,dm):
	"""Compute the dispersion delay (s) as a function of frequency (MHz) and DM"""
	return 4.148808*dm*(10.0**3)/(f**2)

def mean_sub(mat):
	for i in xrange(0,mat.shape[0]):
		mat[i,:] = mat[i,:] - np.mean(mat[i,:]) 

def cust_mark_inset(parent_axes, inset_axes, x1,x2,y1,y2, loc11,loc12,loc21,loc22,**kwargs):
	print parent_axes.transData
	print inset_axes.viewLim
	rect = TransformedBbox(Bbox.from_extents([x1,y1,x2,y2]), parent_axes.transData)
  
	pp = BboxPatch(rect, **kwargs)
	parent_axes.add_patch(pp)

	p1 = BboxConnector(inset_axes.bbox, rect, loc1=loc11,loc2=loc12, **kwargs)
	inset_axes.add_patch(p1)
	p1.set_clip_on(False)
	p2 = BboxConnector(inset_axes.bbox, rect, loc1=loc21,loc2=loc22, **kwargs)
	inset_axes.add_patch(p2)
	p2.set_clip_on(False)

def normalize(npmat):
	mn = npmat.min()
	mx = npmat.max()

	npmat -= mn
	npmat *= (1.0/(mx - mn))

	return mn, mx

def renorm_rebin(mat, unit_row, unit_col):
	ret = np.zeros((math.ceil(mat.shape[0]/unit_row), math.ceil(mat.shape[1]/unit_col)))
	for i in xrange(0,ret.shape[0] - 1):
		for j in xrange(0,ret.shape[1] - 1):
			submat = mat[unit_row*i:min(unit_row*(i+1), mat.shape[0] - 1),unit_col*j:min(unit_col*(j + 1), mat.shape[1] - 1)]
			ret[i,j] = (1 - above_mean(submat))*np.mean(submat)
	return ret

def above_mean(submat):
	mn = np.mean(submat)
	count = 0
	for r in submat:
		for elem in r:
			if elem >= mn:
				count += 1
	return round(float(count)/float(submat.shape[0]*submat.shape[1]))

def simple_rebin(dat,f_bins,t_bins):
	ret = np.zeros((dat.shape[0]/f_bins,dat.shape[1]/t_bins))
	for i in xrange(0,ret.shape[0]):
		for j in xrange(0,ret.shape[1]):
			ret[i,j] = dat[f_bins*i:f_bins*(i + 1),t_bins*j:t_bins*(j + 1)].mean()
	return ret

def rebin(mat,row=True,col=True):
	ret = np.zeros((mat.shape[0]/(1 + row),mat.shape[1]/(1 + col)))
	print mat.shape
	print ret.shape
	if row:
		for i in xrange(0,mat.shape[0]/2 - 1):
			ret[i,:] = 0.5*(mat[2*i,:] + mat[2*i + 1,:])
	if col:
		for i in xrange(0,mat.shape[1]/2 - 1):
			ret[:,i] = 0.5*(mat[:,2*i] + mat[:,2*i + 1])
	return ret