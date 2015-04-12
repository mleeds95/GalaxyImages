def getMUGsfilelistNdigits(path,simName,nDig):
    
    from subprocess import PIPE, Popen
    
    # Read the list of files from the simulation folder
    #cmd = "ls -r "+path+simName+'/'+simName+".0????z*_halos"
    cmd = "ls -r "+path+simName+'/'+simName+".0????.OxMassFrac"
    flist=Popen(cmd, stdout=PIPE, shell=True).stdout.read().split('\n')
    for i in xrange(np.size(flist)-1):
        flist[i]=flist[i][len(path+simName)+1:]
    
    
    str = np.empty(np.size(flist)-1).astype('a30')
    
    for i in xrange(np.size(flist)-1):
        if(nDig==5):
            str[i]=flist[i][7:12]
        elif(nDig==4):
            str[i]=flist[i][6:11]
        elif(nDig==3):
            str[i]=flist[i][5:10]
    return str,np.size(str)

def loadtipsy(path,simName,outNum,haloId,fRotCen,fhalos):
    
    snap=path+simName+'/'+simName+'.'+outNum
    f = pynbody.load(snap)
    f.physical_units()
    
    if(fhalos==1):
        h=f.halos()
        h.make_grp()
    
    if(fRotCen==1):
        mycenter(h[haloId],haloId,mode='myhyb')
        pynbody.analysis.angmom.faceon(h[haloId],cen=[0,0,0])
    elif(fRotCen==2):
        mycenter(h[haloId],haloId,mode='myhyb')
    
    if(fhalos==1):    
        return f,h
    else:
        return f
        
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
