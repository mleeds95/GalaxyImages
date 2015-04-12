import sys
import matplotlib
import numpy as np
import pickle
import matplotlib.pyplot as plt
from scipy.integrate import trapz
import pynbody

def main():
    path = '/Volumes/PROMISE PEGASUS/MUGS/'
    sim = 'g1536' # <-- To do

    #filelist,nfiles=getMUGsfilelistNdigits(path,sim,len(sim)-1L)
    filelist = ['00007', '00013', '00016', '00032', '00036', '00048', '00050', '00059', '00064', '00073', '00080', '00092', '00096', '00104', '00112', '00120', '00128', '00141', '00144', '00160', '00168', '00176', '00185', '00192', '00205', '00208', '00224', '00228', '00240', '00256', '00272', '00288', '00291', '00304', '00320', '00333', '00336', '00352', '00368', '00384', '00386', '00400', '00416', '00432', '00448', '00453', '00464', '00480', '00496', '00512', '00521', '00528', '00544', '00560', '00576', '00592', '00606', '00608', '00624', '00640', '00656', '00672', '00688', '00704', '00713', '00720', '00736', '00752', '00768', '00777', '00784', '00800', '00816', '00832', '00848', '00849', '00912', '00928', '00931', '00944', '00960', '00976', '00992', '01008', '01024']
    nfiles = 85

    snapids = np.zeros(nfiles, dtype=np.int8)
    snapidd = np.zeros(nfiles, dtype=np.int8)

    snapids[:]=-999
    snapidd[:]=-999


    haloidd = 1
    haloids = 1

    for i in range(nfiles-1):

        snap1 = pynbody.load(path+sim+'/'+sim+'.'+filelist[i])
        snap1.physical_units()
        idx = np.zeros(shape=(2,len(snap1.d)))
        halo1 = snap1.halos()
        halo1.make_grp()
        if len(halo1) == 0:
            print sim + '.' + filelist[i] + ' has no halos.'
            continue
        mycenter(halo1[1],1,mode='myhyb')
        pynbody.analysis.angmom.faceon(halo1[1],cen=[0,0,0])
            
        snap2 = pynbody.load(path+sim+'/'+sim+'.'+filelist[i+1])
        snap2.physical_units()
        halo2 = snap2.halos()
        halo2.make_grp()
        if len(halo2) == 0:
            print sim + '.' + filelist[i+1] + ' has no halos.'
            continue
        mycenter(halo2[1],1,mode='myhyb')
        pynbody.analysis.angmom.faceon(halo2[1],cen=[0,0,0])

        # Use dark matter to walk across the snapshots        
        ind = snap1.d['grp']==haloidd
        
        ids = snap2.d['grp'][ind]
        uids = np.unique(ids)
        
        if(np.size(uids)>0):
            muids = np.concatenate((uids,[uids.max()+1]))        
            hh,hids=np.histogram(ids,bins=muids)
            
            ind = hids[:-1]!=0
            if(hh[ind].sum()>0):
                imax = hh[ind].argmax()
                haloidd = hids[ind][imax]
                snapidd[i]=haloidd
                #print "==================="    
                #print "DD", snapidd
            
        # use the stars to walk across the snapshots
        ind = snap1.s['grp']==haloids 
       
        nstar = len(snap2.s)       
        ind2=ind[:nstar]
       
        ids = snap2.s['grp'][ind2]
        uids = np.unique(ids)
       
        if(np.size(uids)>0):
            muids = np.concatenate((uids,[uids.max()+1]))        
            hh,hids=np.histogram(ids,bins=muids)

            ind = hids[:-1]!=0
            if(hh[ind].sum()>0):
                imax = hh[ind].argmax()
                haloids = hids[ind][imax]
                snapids[i]=haloids
                #print "==================="    
                #print "SS", snapids
           
            
    np.savetxt('/Users/mwleeds/merger-trees/mainbranch_'+sim+'_star.txt',snapids,fmt='%d')
    np.savetxt('/Users/mwleeds/merger-trees/mainbranch_'+sim+'_dark.txt',snapidd,fmt='%d')
       
def mycenterpot(sim,haloid):
    i = sim["phi"][sim['grp']==haloid].argmin()
    return sim["pos"][sim['grp']==haloid][i].copy()
    
def my_hybrid_center(sim,haloid, r='3 kpc', **kwargs):
    """

    Determine the center of the halo by finding the shrink-sphere
    -center inside the specified distance of the potential minimum

    """

    try:
        cen_a = mycenterpot(sim,haloid)
    except KeyError:
        cen_a = pynbody.analysis.halo.center_of_mass(sim)
    return pynbody.analysis.halo.shrink_sphere_center(sim[pynbody.filt.Sphere(r, cen_a)], **kwargs)

def mycenter(sim,haloid=1, mode=None, retcen=False, vel=True, cen_size="1 kpc", move_all=True, wrap=False, **kwargs):
    """

    Determine the center of mass of the given particles using the
    specified mode, then recenter the particles (of the entire
    ancestor snapshot) accordingly

    Accepted values for *mode* are

      *mypot*: potential minimum including only host particles

      *myhyb*: for sane halos, returns the same as ssc, but works faster by
             starting iteration near potential minimum

    or a function returning the COM.

    **Other keywords:**

    *retcen*: if True only return the center without centering the
     snapshot (default = False)


    *vel*: if True, translate velocities so that the velocity of the
    central 1kpc (default) is zeroed. Other values can be passed with cen_size.

    *move_all*: if True (default), move the entire snapshot. Otherwise only move
    the particles in the halo passed in.

    *wrap*: if True, pre-centre and wrap the simulation so that halos on the edge
    of the box are handled correctly. Default False.
    """

    if mode is None:
        mode = 'myhyb'

    try:
        fn = {'mypot': mycenterpot,
              'myhyb': my_hybrid_center}[mode]
    except KeyError:
        fn = mode

    if move_all:
        target = sim.ancestor
    else:
        target = sim

    if wrap:
        # centre on something within the halo and wrap
        target = pynbody.transformation.inverse_translate(target, sim['pos'][0])
        target.sim.wrap()

    print fn(sim,haloid, **kwargs)
    if retcen:
        return fn(sim,haloid, **kwargs)
    else:
        cen = fn(sim,haloid, **kwargs)
        tx = pynbody.transformation.inverse_translate(target, cen)

    if vel:
        velc = pynbody.analysis.halo.vel_center(sim, cen_size=cen_size, retcen=True)
        tx = pynbody.transformation.inverse_v_translate(tx, velc)

    return tx


if __name__=='__main__':
    main()
