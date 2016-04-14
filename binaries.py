import numpy as np
import ctypes as C
import os
import copy
from os.path import join
import h5py


pc_in_AU = 2.06e5 
AU_to_parsec = 4.84e-6


try:
    lib_path=os.path.join( os.environ["BINARIESLIB"],"libbinaries.so")
except KeyError:
    lib_path="/usr/local/lib/libbinaries.so"
_lib = C.CDLL(lib_path)

double_pointer = C.POINTER(C.c_double)
int_pointer = C.POINTER(C.c_int)





class Binaries:
    """
    This object can be built from a StarFiddle HDF5 run:
      Run.WriteHDF5( "path/to/run" ) # writes path/to/run.hdf5
      BinStars = Binaries(10, "path/to/run.hdf5" )
    This writes a multi-step Binaries object, tracking binaries over time.

    Another way is to build a single step Binaries object directly from data
      Binaries.FromData(n1,n2,mt,m1,m2,cmr,a,p,e,density_ratio)
   
    Arrays :
      n1         : identities of all primary stars
      n2         : identities of all secondary stars
      a          : semi-major axis of binaries
      e          : eccentricities of binaries
      mt         : total masses
      m1         : primary masses
      m2         : secondary masses
      cmr        : distance of center of mass to origin
      p          : period in nbody time
    These are attributes of Binaries. For multi step objects, they are
    2-dimensionnal and contain zeros when a given binary is not defined at a given step.

"""
    def __init__(self, density_ratio, path=None, coordinates=None, Nnb=8,
                 G=1, RemoveTriples=False):
        self.attributes = ["n1","n2","mt","m1","m2","cmr","a","p","e"]
        if coordinates is not None and path is not None:
            raise Exception("Building Binaries object both from file and coordinates")
        if coordinates is not None:
            self.density_ratio = density_ratio
            self.nstep = 1
            self.Nbin, \
            self.n1,self.n2,self.a,self.e,self.mt,self.m1,self.m2,\
            self.cmr,self.p,self.Ebin,self.ratios = \
                    collect_binaries_from_data(*coordinates,
                                               density_ratio=density_ratio,
                                               Nnb=Nnb,G=G)
        if coordinates is None and density_ratio is not None:
            self.density_ratio = density_ratio
            self.Nbin, self.nstep, \
            self.n1,self.n2,self.a,self.e,self.mt,self.m1,self.m2,\
              self.cmr,self.p = follow_binaries(density_ratio,filename=path,G=G)
            self.snaprange = [0, self.nstep-1]
            if RemoveTriples:
                self.RemoveTriples()

    @classmethod
    def FromData(self, n1=None, n2=None, mt=None, m1=None, m2=None, 
                 cmr=None, a=None, p=None, e=None, density_ratio=None):
        Bin = Binaries(None,None)
        for name,arg in zip(["n1","n2","mt","m1","m2","cmr","a","p","e"],
                            [n1,n2,mt,m1,m2,cmr,a,p,e]):
            setattr(Bin,name,arg)
        Bin.Nbin = n1.shape[0]
        if np.size(n1) == Bin.Nbin:
            Bin.nstep = 1
        else: 
            Bin.nstep = Bin.n1.shape[1]
        Bin.density_ratio = density_ratio
        return Bin

    def copy(self):
        return copy.deepcopy(self)

    def LongerThan(self,n_occ, InternalTransformation = False):
        """
        Return a Binaries instance with binaries lasting for more
        timesteps than n_occ.
        """
        ind = self.n1 != 0
        validsteps = np.zeros(self.Nbin)
        for i,sub_ind in enumerate(ind):
            validsteps[i] = len(np.nonzero(sub_ind)[0])
        valid_index = np.nonzero(validsteps > n_occ)[0]
        newN = len(valid_index)
        n1,n2,mt,m1,m2,cmr,a,p,e = [ getattr(self,attr)[valid_index[:],:]
                                     for attr in self.attributes ]
        if InternalTransformation:
            self.n1 = n1; self.n2 = n2; self.mt = mt; 
            self.m1 = m1; self.m2 = m2; self.cmr = cmr; 
            self.a = a; self.p = p; self.e = e
            self.Nbin = newN
        else:
            return Binaries.FromData(n1,n2,mt,m1,m2,cmr,a,p,e,
                                     self.density_ratio)

    def RemoveZerosFromSingleStep(self):
        attributes = ["n1","n2","mt","m1","m2","cmr","a","p","e"]
        if self.nstep>1:
            raise Exception("This method cannot be applied to multiple steps objects")
        ind = self.n1!=0
        self.n1,self.n2,self.mt,self.m1,self.m2,self.cmr,self.a,self.p,self.e = \
                    [ getattr(self,attr)[ind] for attr in attributes ]
        self.Nbin = len(self.n1)

    def AtStep(self,step, TrackIndices=False):
        """
        Return a Binaries instance of specified timestep
        """
        if self.nstep <= 1:
            if TrackIndices:
                raise Exception("Asking to track indices on a single step object")
            return self

        n1s = self.n1[:,step]

        # Indexes of the binaries that exist at the requested step:
        ind =  n1s != 0  
        arr = [self.n1,self.n2,self.mt,self.m1,self.m2,
               self.cmr,self.a,self.p,self.e]
        n1,n2,mt,m1,m2,cmr,a,p,e = [x[ind,step] for x in arr]
        if TrackIndices:
            return np.arange(self.Nbin)[ind], \
                Binaries.FromData(n1,n2,mt,m1,m2,cmr,a,p,e, 
                                self.density_ratio)
        else:
            return Binaries.FromData(n1,n2,mt,m1,m2,cmr,a,p,e, 
                                     self.density_ratio)
            
            
    def Simplified(self,attr, scaling=None, log_axis=True):
        """
        Return a list of numpy arrays, one per time step, of 
        the requested attributes, cleaned of all zeros. This is for
        statistical studies as the identity of each pair is lost.
        For example, someone wants the periods of existing binaries
        at each steps:
            periods = Binaries.Simplied("p")
        Or the same for semi-major axis:
            axis = Binaries.Simplied("a")
        """
        if not hasattr(self,attr):
            raise Exception("Requested attribute was not found. "
                            "Please choose in the following list:\n"
                            " \"n1\" -> index of primary star"
                            " \"n2\" -> index of secondary star"
                            " \"mt\" -> total mass of the pair"
                            " \"m1\" -> primary mass "
                            " \"m2\" -> secondary mass"
                            " \"cmr\" -> center of mass distance to origin"
                            " \"a\" -> semi-major axis"
                            " \"p\" -> period"
                            " \"e\" -> eccentricity ")
        data=[]
        if attr is "a" and scaling is not None:
            factor = scaling / AU_to_parsec
        else:
            factor = 1
        func = np.log10 if log_axis else lambda x:x
        for step in range(self.nstep):
            n1s = self.n1[:,step]
            ind = n1s!=0
            data.append(func( factor * getattr(self,attr)[ind,step] ) ) 
        return data
        
    def ReturnMemberList(self,Unique=True):
        """
        Return a numpy array of all binary members
        """
        if Unique:
            return np.unique(np.array( list(self.n1) + list(self.n2) ))
        else:
            return np.array( list(self.n1) + list(self.n2) )

    def Permutation(self, ind):
        """
        Takes a numpy array of indices and changes the order of the
        arrays along it. For example, if you want to reorganize the
        instance to have everything sorted by semi-major axis:
             B = Binaries(.........).AtStep(0)
             ind = numpy.argsort(B.a)
             B.Permutation(ind)
        Only works for 1 step instances.
        """
        if self.nstep != 1:
            raise Exception("For now, only works for nstep=1")
        arr = [self.n1,self.n2,self.mt,self.m1, self.m2, 
               self.cmr, self.a,self.p,self.e]
        [self.n1,self.n2,self.mt,self.m1,self.m2,
              self.cmr,self.a,self.p,self.e] = [x[ind] for x in arr]
        self.Nbin = len(self.n1)



    def PrintIndexes(self,filename="Binaries_n1n2"):
        """
        Print all binaries members as one line per binary and one 
        column per step. At steps the binaries doesn't exist, blank 
        spaces are printed.
        """
        
        with open(filename,"w") as f:
            for i in range(self.nstep):
                for n in range(self.Nbin):
                    f.write(format_nstars(self.n1[n,i],self.n2[n,i]))
                f.write("\n")

    def RemoveTriples(self,ExtractTriples=True):
        """
        Remove from the object all binaries with members occuring in
        more than one pair. Effectively removing non-hierarchical
        triple, quadruple, etc.
        Stores them in self.Triples as a Binaries object
        """
        if hasattr(self,"OnlyTriples"):
            print "Triples were already extracted from the object."
            return
        print "Removing triple, quadruple, etc- systems"
        attributes = ["n1","n2","mt","m1","m2","cmr","a","p","e"]
        triple_container = []
        if ExtractTriples:
            self.OnlyTriples = self.copy()
            self.IncludeTriples = self.copy()
        # We go through each step and search for triples
        if self.nstep==1:
            n1,n2,mt,m1,m2,cmr,a,p,e = [ getattr(self,attr) 
                                         for attr in attributes ]
            n = np.hstack((n1,n2))
            name,count = np.unique(n, return_counts=True)
            triple_members = name [ np.nonzero(count-1)[0] ]
            triple_indices = np.nonzero( np.in1d(n,triple_members) )[0]
            triple_indices = np.unique( triple_indices % self.Nbin)
            if len(triple_indices) is not 0:
                for attr in attributes:
                    getattr(self,attr)[ triple_indices ] = 0
                if ExtractTriples:
                    non_triple = np.setdiff1d(np.arange(self.Nbin),triple_indices)
                    for attr in attributes:
                        getattr(self.OnlyTriples,attr)[ non_triple ] = 0
            self.RemoveZerosFromSingleStep()
            if ExtractTriples:
                self.OnlyTriples.RemoveZerosFromSingleStep()
            return
        for step in range(self.nstep):
            orig_ind,loc_bins = self.AtStep(step,TrackIndices=True)
            n1,n2,mt,m1,m2,cmr,a,p,e = [ getattr(loc_bins,attr) for attr in attributes ]
            n = np.hstack((n1,n2))
            name,count = np.unique(n, return_counts=True)
            triple_members = name [ np.nonzero(count-1)[0] ]
            triple_indices = np.nonzero( np.in1d(n,triple_members) )[0]
            triple_indices = np.unique( triple_indices % loc_bins.Nbin)

            # The indices of triples in the AtStep object are in
            # ind . To get their indices in the self object
            # I use the orig_ind array I created in the beginning 
            # which remembers the indices in self.
            # Now I set n1, n2, etc to 0 where the binary is a triple
            if len(triple_indices) is not 0:
                for attr in attributes:
                    if self.nstep>1:
                        getattr(self,attr)[ orig_ind[ triple_indices ], step] = 0
                    else:
                        getattr(self,attr)[ orig_ind[ triple_indices ]] = 0

            if ExtractTriples:
                all_ind = np.arange(len(orig_ind))
                non_triple = np.setdiff1d(all_ind,triple_indices)
                for attr in attributes:
                    if self.nstep>1:
                        getattr(self.OnlyTriples,attr)[orig_ind[non_triple], step] = 0    
                    else:
                        getattr(self.OnlyTriples,attr)[ orig_ind[non_triple] ] = 0    

        self.LongerThan(0,InternalTransformation=True)
        if ExtractTriples:
            self.OnlyTriples.LongerThan(0,InternalTransformation=True)
        
        

    def GetFullPopulation(self,scaling):
        if hasattr(self,"OnlyTriples"):
            return self.IncludeTriples.Simplified("a",scaling=scaling)
        else:
            return self.Simplified("a",scaling=scaling)


    def Remove(self,ind):
        if self.nstep==1:
            for attr in self.attributes:
                getattr(self,attr)[ ind ] = 0
            self.RemoveZerosFromSingleStep()
        else:
            for attr in self.attributes:
                for i in range(self.nstep):
                    getattr(self,attr)[ ind , i]  = 0
            self.LongerThan(0,InternalTransformation=True)


    def RemoveTriplesFromSubsample(self,Triples,ExtractTriples=True):
        if hasattr(self,"OnlyTriples"):
            print "Triples were already extracted from the object."
            return
        print "Removing triple, quadruple, etc- systems"
        triple_container = []
        if ExtractTriples:
            self.OnlyTriples = self.copy()
            self.IncludeTriples = self.copy()
        # We go through each step and search for triples
        for step in range(self.nstep):
            tB = Triples.AtStep(step) 
            orig_ind ,B = self.AtStep(step,TrackIndices=True)
            bool1 = np.in1d(B.n1, tB.n1)
            bool2 = np.in1d(B.n2, tB.n2)
            ind = np.nonzero(bool1*bool2)[0]

            if len(ind) is not 0:
                for attr in self.attributes:
                    getattr(self,attr)[ orig_ind[ ind ], step] = 0

            if ExtractTriples:
                all_ind = np.arange(len(orig_ind))
                non_triple = np.setdiff1d(all_ind,ind)
                for attr in self.attributes:
                    getattr(self.OnlyTriples,attr)[orig_ind[non_triple], step] = 0                

        self.LongerThan(0,InternalTransformation=True)
        if ExtractTriples:
            self.OnlyTriples.LongerThan(0,InternalTransformation=True)




#===============================================================================
# Various functions
#===============================================================================


def snapname(nsnap):
    """snapname(17) -> snap_0000017 """
    prefix = "snap_"
    suffix = "%07d" % nsnap
    name = prefix + suffix
    return name

def locate(array, target):
    """Find target in array, raise an exception if not found."""
    _lib.locate.argtypes = [C.c_int, int_pointer, C.c_int]
    _lib.locate.restype = C.c_int
    N = len(list(array))
    c_array = ( N*C.c_int )()
    for i in range(N): c_array[i] = int( array[i] )
    result =  _lib.locate(N,c_array,target)
    if result == -1:
        raise Exception("The target could not be found in the array.")
    else:
        return result


def get_a(r,v,mt,G=1):
    """Compute the semi major axis of a binary from r, v, mt and G"""
    _lib.get_a.argtypes = [C.c_double, C.c_double, C.c_double,
                                    C.c_double,C.c_double, C.c_double,
                                    C.c_double, C.c_double]
    _lib.get_a.restype = C.c_double
    return _lib.get_a(r[0],r[1],r[2],v[0],v[1],v[2],mt,G)



def inv_distance(xi,xj,yi,yj,zi,zj):
    """returns the inverse distance between two particles"""
    _lib.inv_distance.argtypes = [C.c_double, C.c_double, 
                                           C.c_double, C.c_double, 
                                           C.c_double, C.c_double]
    _lib.inv_distance.restype = C.c_double
    return _lib.inv_distance(xi,xj,yi,yj,zi,zj)
 

def ebin(x, y, z, vx, vy, vz, mt, G=1):
    """
    Returns the binding energy of two stars from relative distance, 
    velocity and total mass.
    """
    _lib.ebin.argtypes = [C.c_double, C.c_double, C.c_double,
                                   C.c_double, C.c_double, C.c_double, 
                                   C.c_double, C.c_double]
    _lib.ebin.restype = C.c_double
    return _lib.ebin(x, y, z, vx, vy, vz, mt, G)


#===========================================================================
# COLLECT_BINARIES_FROM_DATA
#===========================================================================
def collect_binaries_from_data(n,m,x,y,z,vx,vy,vz,
                               density_ratio=10,Nnb=5,G=1):
    """
    Nbin , n1, n2, a, e, mt, m1, m2, cmr, p, Ebin_array, ratio_array =
       collect_binaries_from_data(n, m, x, y, z, vx, vy, vz,
                                  density_ratio=10, Nnb=5, G=1)
    ---------------------------------------------------------------------
    INPUT
    ---------------------------------------------------------------------
      n            : Stars identities (integer>0)
      m            : Stars masses
      x,y,z        : Stars positions
      vx,vy,vz     : Stars velocities
      density_ratio: How much denser a binary has to be compared to its
                     environment to be registered as valid.
      Nnb          : Number of neighbours to compute local density
      G            : Gravitationnal constant. Nbody units: G=1
    ---------------------------------------------------------------------
    OUTPUT
    ---------------------------------------------------------------------
      Nbin       : Number of found binaries
      n1         : identities of all primary stars
      n2         : identities of all secondary stars
      a          : semi-major axis of binaries
      e          : eccentricities of binaries
      mt         : total masses
      m1         : primary masses
      m2         : secondary masses
      cmr        : distance of center of mass to origin
      p          : period in nbody time
      Ebin_array : binding energies
      ratio_array: ratio of binary density to local density

    Take a dynamical system, spot bound pairs of stars and compare the local
    neighbour density to the density defined by the two binaries components.
    If the ratio exceeds density_ratio, the binary is registered.

    """
#    
    _lib.collect_binaries_from_data.argtypes = [C.c_int,
                                                int_pointer,double_pointer,
                                                double_pointer, double_pointer, 
                                                double_pointer, double_pointer, 
                                                double_pointer, double_pointer, 
                                                int_pointer, int_pointer, 
                                                double_pointer, double_pointer, 
                                                double_pointer, double_pointer, 
                                                double_pointer, double_pointer, 
                                                double_pointer, double_pointer,
                                                double_pointer,
                                                C.c_double, C.c_int, C.c_double]
    _lib.collect_binaries.restype = C.c_int
    N = len(n)
    nc = (N*C.c_int)();  nc[:]=n; n=nc 
    mc = (N*C.c_double)();  mc[:]=m; m=mc 
    xc = (N*C.c_double)();  xc[:]=x; x=xc 
    yc = (N*C.c_double)();  yc[:]=y; y=yc 
    zc = (N*C.c_double)();  zc[:]=z; z=zc 
    vxc = (N*C.c_double)();  vxc[:]=vx; vx=vxc 
    vyc = (N*C.c_double)();  vyc[:]=vy; vy=vyc 
    vzc = (N*C.c_double)();  vzc[:]=vz; vz=vzc 
#
    size = 10*N
    n1 = (size*C.c_int)()
    n2 = (size*C.c_int)()
    Ebin_array = (size*C.c_double)()
    a = (size*C.c_double)()
#
    e = (size*C.c_double)()
    mt = (size*C.c_double)()
    m1 = (size*C.c_double)()
    m2 = (size*C.c_double)()
    cmr = (size*C.c_double)()
    p = (size*C.c_double)()
    ratio_array = (size*C.c_double)()
#    
    nbin = _lib.collect_binaries_from_data(N,n,m,x,y,z,vx,vy,vz,
                                           n1,n2,a,e,mt,m1,m2,cmr,p,Ebin_array, ratio_array,
                                           density_ratio, Nnb, G )
# 
    Ebin_array = np.asarray(Ebin_array[0:nbin] , dtype=np.double) 
    n1 = np.asarray(n1[0:nbin],dtype=np.int)
    n2 = np.asarray(n2[0:nbin],dtype=np.int)
    a = np.asarray(a[0:nbin] , dtype=np.double) 
    e = np.asarray(e[0:nbin] , dtype=np.double) 
    mt = np.asarray(mt[0:nbin] , dtype=np.double) 
    m1 = np.asarray(m1[0:nbin] , dtype=np.double) 
    m2 = np.asarray(m2[0:nbin] , dtype=np.double) 
    cmr = np.asarray(cmr[0:nbin] , dtype=np.double) 
    p   = np.asarray(  p[0:nbin] , dtype=np.double) 
    ratio_array = np.asarray(ratio_array[0:nbin] , dtype=np.double) 
#
    return int(nbin) , n1, n2, a, e, mt, m1, m2, cmr, p, Ebin_array, ratio_array



#===========================================================================
# FOLLOW_BINARIES
#===========================================================================
def follow_binaries(density_ratio,
                    filename="/home/dorval/heavy_data/nbody6/test/test/run64.hdf5",
                    Nnb=10,
                    G=1):
    """
    Nbin,nstep,n1,n2,a,e,mt,m1,m2,cmr,p = follow_binaries( density_ratio,
                                                           filename, Nnb, G)
    ---------------------------------------------------------------------
    INPUT
    ---------------------------------------------------------------------
      density_ratio: How much denser a binary has to be compared to its
                     environment to be registered as valid.
      filename     : Path to a HDF5 file of a StarFiddle Run
      G            : Gravitationnal constant. Nbody units: G=1
    ---------------------------------------------------------------------
    OUTPUT
    ---------------------------------------------------------------------
      Nbin       : Number of found binaries
      nstep      : Number of snapshots in the run
      [ The following are 2d arrays: (Nbin,nstep). If a binary does not exist
        at a given steps, the value is 0 ]
      n1         : identities of all primary stars
      n2         : identities of all secondary stars
      a          : semi-major axis of binaries
      e          : eccentricities of binaries
      mt         : total masses
      m1         : primary masses
      m2         : secondary masses
      cmr        : distance of center of mass to origin
      p          : period in nbody time

    Take the path to a HDF5 Run written by StarFiddle and go through all
    snapshots looking for binaries with the provided density ratio (see
    collect_binaries_from_data docstring).

    """
    _lib.follow_binaries.argtypes = [C.c_char_p, 
                                     C.c_int, C.c_int, C.c_double, 
                                     int_pointer, int_pointer, 
                                     double_pointer, double_pointer,
                                     double_pointer, double_pointer,
                                     double_pointer, double_pointer, 
                                     double_pointer, C.c_double]
    _lib.follow_binaries.restype = C.c_int
    # A binary should be at least denser than its environment.
    if density_ratio < 1: 
        raise Exception("The ratio must be superior to 1.")
#        
     # A filename was specified, but the file can't be found
    if not os.path.isdir(filename) and not os.path.isfile(filename):
            raise Exception(filename+" not found...")
    # Reading the HDF5 file to know how many steps there are in the run
    f =  h5py.File(filename, "r")
    nstep = f.attrs["nstep"];
    N = f["0"].attrs["N"];
    f.close()
    #raise(Exception("Stop here for now."))
#
    #Creating the arrays that will get filled by the c function
    nbinmax = 5*N
    size = nbinmax*nstep
    n1 = (size*C.c_int)()
    n2 = (size*C.c_int)()
    a = (size*C.c_double)()
    e = (size*C.c_double)()
    mt = (size*C.c_double)()
    m1 = (size*C.c_double)()
    m2 = (size*C.c_double)()
    cmr = (size*C.c_double)()
    p = (size*C.c_double)()
#
    nbintot = _lib.follow_binaries(filename,nbinmax,Nnb,
                                   density_ratio,
                                   n1,n2,a,e,mt,m1,m2,cmr,p,G)
    # We turn our ctypes objects into numpy arrays
    n1 =  np.asarray(n1[0:nbintot*nstep],dtype=np.int)
    n2 =  np.asarray(n2[0:nbintot*nstep],dtype=np.int)
    mt =  np.asarray(  mt[0:nbintot*nstep] , dtype=np.double) 
    m1 =  np.asarray(  m1[0:nbintot*nstep] , dtype=np.double) 
    m2 =  np.asarray(  m2[0:nbintot*nstep] , dtype=np.double) 
    cmr =  np.asarray(  cmr[0:nbintot*nstep] , dtype=np.double) 
    a  =  np.asarray(a[0:nbintot*nstep] , dtype=np.double) 
    p  =  np.asarray(  p[0:nbintot*nstep], dtype=np.double  )
    e  =  np.asarray(e[0:nbintot*nstep] , dtype=np.double) 
#
    if nstep != 1:
        # The c function fills 1d array, putting them back into 2d arrays
#        print "Reshaping...."
        n1 = np.reshape(n1,(nbintot,nstep))
        n2 = np.reshape(n2,(nbintot,nstep))
        mt = np.reshape(mt,(nbintot,nstep))
        m1 = np.reshape(m1,(nbintot,nstep))
        m2 = np.reshape(m2,(nbintot,nstep))
        cmr = np.reshape(cmr,(nbintot,nstep))
        a = np.reshape(a,(nbintot,nstep))
        p = np.reshape(p,(nbintot,nstep))
        e = np.reshape(e,(nbintot,nstep))
#   print "Done."
    return nbintot,nstep,n1,n2,a,e,mt,m1,m2,cmr,p



def format_nstars(n1,n2):
    """
    Returns \"  n1/n2\  " in a specific format. 
    Used by Binaries.PrintIndexes
    """
    if n1 == 0 and n2 == 0:
        return "      /       "
    elif n1 == 0 and n2 != 0:
        n2 = str(n2).ljust(7)
        return "       /"+n2
    elif n1!=0 and n2 == 0:
        return " %6d/       " %n1
    else:
        n2 = str(n2).ljust(7)
        string = "%6d/" % n1
        return string+n2







