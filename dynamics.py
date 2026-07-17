import numpy as np


# This code assumes that pbc_utils has already been used to
# rebuild the trajectories in time (removing PBC jumps)
def msd(rt, bx, nlogpoints=-1):
    nframes = np.shape(rt)[0]
    ns = np.shape(rt)[1]

    n_delts = nframes
    # initialize log space in here if needed

    ##

    msd = np.zeros((n_delts, 4))
    delt_list = np.zeros(n_delts)

    for dt_ind in range(1,n_delts):
        delt = dt_ind
        
        # switch to lgo space in here if needed

        ##


        delt_list[dt_ind] = delt

        for t in range(0, nframes - delt):
            
            dr = rt[t,:,:] - rt[t+delt,:,:]

            dr2 = dr**2


            msd[dt_ind,0:3] = np.sum(dr2,axis=0) / ns
            msd[dt_ind,3] = np.sum(msd[dt_ind,0:3])/3.0
            


    return msd, delt_list
