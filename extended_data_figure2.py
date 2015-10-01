import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import math
#import brewer2mpl
from scipy.special import erf
from scipy.optimize import leastsq
from scipy.interpolate import spline

delta_t = 0.00102399999741

def run(path,dm=623.297,bins=3,fit_color='black',colors=['blue','darkgreen','red'],linestyles=['--', '-.', ':','-'],labels=['high', 'middle','low'],fitlable='middle-fit',data_rep='default'):
    #raw = np.load(path)
    #raw = np.nan_to_num(raw)
    #Use the full intensity
    #dat =  raw[:,0,:]
    #del raw

    dat = np.fromfile(path)
    dat = np.reshape(dat,(3836,150))
    #dat.shape
    #dat = dat.T
    #plt.figure()
    #ax = plt.gca()
    #ax.imshow(dat,interpolation='nearest',aspect=float(dat.shape[1])/float(dat.shape[0]))
    print dat.shape
    power = integrated_power(dat,bins=bins,dm=dm)
    #print power

    f=plt.figure(1)

    offset = 5
    #colors = ('blue','red','darkgreen')
    #fit_color = 'black'
    #colors = cm.rainbow(np.linspace(0,1,len(power)))
    #bmap = brewer2mpl.get_map('Set1', 'qualitative', 4)
    #colors = bmap.mpl_colors

    q = power[1]
    q *= 1.0/float(max(q[0,:]))

    x = np.linspace(-10 + offset,20 + offset,1000)
    d_slice = q[0,35 + offset:66 + offset]
    fit_slice = q[0,45:60]
    x_old = xrange(0,len(d_slice))
    power_smooth = spline(x_old,d_slice,np.linspace(0,30,len(x)))

    params = optimize_exp(fit_slice)
    print params
    params = params[0]
    #plt.plot(p[0,:])


    plt.plot(np.linspace(0,30,len(x)),conv_exp(x,params[0],params[1],params[2],params[3]),color=fit_color,linewidth=2.0,linestyle=linestyles[3],label=fitlable)

    for i in xrange(0,3):
        p = power[i]
        p *= 1.0/float(max(p[0,:]))


        x = np.linspace(-10 + offset,20 + offset,1000)
        d_slice = p[0,35 + offset:66 + offset]
        fit_slice = p[0,45:60]
        x_old = xrange(0,len(d_slice))
        power_smooth = spline(x_old,d_slice,np.linspace(0,30,len(x)))

        params = optimize_exp(fit_slice)
        print params
        params = params[0]
        #plt.plot(p[0,:])


#        if i == 1:
#            #plt.gca().set_color_cycle(['black','red', 'green', 'blue'])
#            plt.plot(np.linspace(0,30,len(x)),conv_exp(x,params[0],params[1],params[2],params[3]),color=fit_color,linewidth=2.0,linestyle=linestyles[3],label=fitlable)


        if data_rep == 'scatter':
            plt.scatter(xrange(0,len(d_slice)),d_slice,marker='+',color=colors[i])
        elif data_rep == 'interp':
            plt.plot(np.linspace(0,50,len(x)),power_smooth,color=colors[i],linestyle=linestyles[i])
        elif data_rep == 'default':
            plt.plot(x_old,d_slice,color=colors[i],linewidth=2.0,linestyle=linestyles[i],label=labels[i])
            #plt.scatter(xrange(0,len(d_slice)),d_slice,marker='+')

#        params = optimize_exp(fit_slice)
#        print params 
#        params = params[0]
        #plt.plot(p[0,:])


#        if i == 1:
#            #plt.gca().set_color_cycle(['black','red', 'green', 'blue'])
#            plt.plot(np.linspace(0,50,len(x)),conv_exp(x,params[0],params[1],params[2],params[3]),color=fit_color,linewidth=2.0,linestyle=linestyles[3],label=fitlable)
        plt.ylabel('Normalized Intensity', fontsize=16)
        plt.xlabel('Time (ms)', fontsize=16)

        xlabs = ['%i' % int(j*delta_t*1000.0) for j in np.linspace(-10,20,7)]

        plt.xticks(np.linspace(0,30,7),xlabs, fontsize=16)
        plt.tick_params(axis='both', which='major', labelsize=16)

    plt.legend(loc='upper right', fontsize=16,frameon=False)
    plt.savefig('burst_profile_0929.eps', bbox_inches='tight')
    plt.close(f)
def integrated_power(dat,dm=623,dt = 0.00102399999741,bins=10):
    inc = float(dat.shape[0] - 1)/bins
    ret = []
    #dat,delay_ind = dedisperse(dat,dm,dt)
    for i in xrange(0,bins):
        subdat = dat[i*inc:(i + 1)*inc,:]
        elem = np.zeros((1,subdat.shape[1]))
        for i in xrange(0,elem.shape[1]):
            elem[0,i] = np.sum(subdat[:,i])
        ret.append(elem)
    return ret

def conv_exp(x, A, var, tau, mean):
    return A*np.exp(-(x - mean)/tau + var/(4.0*(tau**2)))*(1.0 + erf((2*(x - mean) - var/tau)/(2.0*np.sqrt(var))))/(2.0*np.sqrt(1.0/var))

def residuals(p, y, x):
    A, var, tau, mean = p
    return y - conv_exp(x, A, var, tau, mean)

def optimize_exp(dataset, p0 = [1.0,1.0,1.0,0.0]):
    x = np.linspace(0,len(dataset)-1,len(dataset))
    return leastsq(residuals, p0, args=(dataset, x))

def disp_delay(f,dm):
    """Compute the dispersion delay (s) as a function of frequency (MHz) and DM"""
    return 4.148808*dm*(10.0**3)/(f**2)

def dedisperse(dat,dm=623,dt = 0.00102399999741):
    ret = np.zeros(dat.shape)
    df = -200.0/float(dat.shape[0])
    f0 = 900
    for i in xrange(0,dat.shape[0]):
        f = f0 + i*df
        dm = 623.0
        delay_ind = int(round((disp_delay(f,dm) - disp_delay(f0,dm))/dt))
        for j in xrange(0,dat.shape[1] - delay_ind):
            ret[i,j] = dat[i,j + delay_ind]
    return ret,int((disp_delay(700,dm) - disp_delay(900,dm))/dt)
