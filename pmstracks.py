import os as _os
import numpy as _np
import astropy.units as _u
import astropy.constants as _c

year  = (1.0*_u.yr).cgs.value
R_sun = _c.R_sun.cgs.value
L_sun = _c.L_sun.cgs.value

class track:
    
    _zs   = ['01','02','03','04']
    _mask = [2,4,6] # column indices of L, Reff, and Teff in the data file
    _it   = -1      # column index of the age in the data file
    
    def __init__(self,m='0.5',z='01',track_dir=None):
        """
        Stellar evolutionary track based on Siess L., Dufour E., Forestini M. 2000, A&A, 358, 593.
        
        Data files need to reside in the folder `tracks/`, if not, the files will be downloaded from
        http://www.astro.ulb.ac.be/~siess/
        
        Keywords:
        ---------
        
        m : string
        :  stellar mass in solar masses as string
        
        z : string
        :  metallicity mass fraction, notation: \'01\' = 0.01. Possible values are"""+', '.join(self._zs)+\
        """
        
        track_dir : string
        :   path to the place where the tracks are located. Defaults to 'tracks/' in the module path.
        """
        
        
        from scipy.interpolate import interp1d
        import glob
        
        if track_dir is None: track_dir = _os.path.join(_os.path.dirname(__file__),'tracks')
        
        self._track_dir      = track_dir
        self._track_file     = _os.path.join(self._track_dir,'Z'+z,'m'+m+'z'+z+'.hrd')
        
        # download files if respective folder is missing
        
        if not _os.path.exists(_os.path.join(self._track_dir,'Z'+z)): self._download_files()
        
        # determine available stellar masses
        
        self._ms = [_os.path.basename(file_)[1:].split('z')[0] for file_ in glob.glob(_os.path.join(self._track_dir,'Z'+z,'*.hrd'))]
        
        if m not in self._ms:
            raise ValueError('Selected mass for z={} not in available masses: {}'.format(z,', '.join(self._ms)))
        
        
        # load data and create interpoation function
        
        self._track_data     = _np.loadtxt(self._track_file)
        self._track_function = interp1d(self._track_data[:,self._it],self._track_data[:,self._mask].T)
        
        
    def _download_files(self,delete=True):
        """
        Downloads all the grid files from Lionel Siess' website and extracts the evolutionary track data files
        
        Keywords:
        ---------
        
        delete : bool
        :   if True, delete the downloaded archives after extracting them
        """
        import urllib2, sys, tarfile
        
        if not _os.path.isdir(self._track_dir): _os.mkdir(self._track_dir)
        
        s = 'Downloading grid files'
        print(s+'\n'+'-'*len(s))
        
        downloaded_files = []
        
        for z in self._zs:
            
            url  = 'http://www.astro.ulb.ac.be/~siess/pmwiki/pmwiki.php/StellarModels/Z{0}?action=dirlistget&f=Grid_z{0}.tar.gz'.format(z)
            
            try:
                
                # open url, get file name and size
            
                remotefile = urllib2.urlopen(url)
                meta      = remotefile.info()
                file_name = meta.getheaders('Content-Disposition')[0].split('"')[1]
                file_size = int(meta.getheaders("Content-Length")[0])
                
                downloaded_files += [[file_name,z]]
            
            
                # download and write to file with progress bar
            
                with open(file_name, 'wb') as f:
                    file_size_dl = 0
                    block_sz = 8192
                    while True:
                        buf = remotefile.read(block_sz)
                        if not buf:
                            break
            
                        file_size_dl += len(buf)
                        f.write(buf)
            
                        sys.stdout.write("\r{}: {:2.2%}".format(file_name,float(file_size_dl)/file_size))
                        sys.stdout.flush()
                    print("")
            except:
                print('Failed to download grid file for z = '+z)
            finally:
                remotefile.close()
                
        # extracting into the specified directory
        
        s = 'Extracting tracks'
        print(s+'\n'+'-'*len(s))
        
        for file_name,z in downloaded_files:
            try:
                print(file_name)
                if file_name.endswith('tar.gz'):
                    mode = 'r:gz'
                elif file_name.endswith('tar'):
                    mode = 'r:'
                with tarfile.open(file_name,mode) as tar:
                    for file_ in tar:
                        if file_.isreg() and file_.name.endswith('.hrd'):
                            file_.name = _os.path.join(self._track_dir,'Z'+z,_os.path.basename(file_.name))
                            if _os.path.isfile(file_.name): _os.unlink(file_.name)
                            tar.extract(file_)
                            
            except Exception, _:  
                import traceback
                print('Could not extract '+file_name+'\n')  
                print('Traceback:')  
                print('----------')  
                for s in traceback.format_exc().split('\n'):  
                    print(4*' '+s)  
                print('----------')

            finally:
                if delete: _os.unlink(file_name)
                
                
    def get_stellar_params(self,t):
        """
        Returns stellar parameters L, Reff, Teff at the given time [s]
        
        Arguments:
        ----------
        
        t : float
        :   age of the star in seconds
        
        Output:
        -------
        L,Reff,Teff
        
        L : float
        :   stellar luminosity [erg/s]
        
        Reff : float
        :   effective radius [cm]
        
        Teff : float
        :   effective temperature [K]
        """
        
        # convert to years
        
        t = t/year
        
        if t<=self._track_data[0,self._it]:
            r = self._track_data[0,self._mask]
        elif t>= self._track_data[-1,self._it]:
            r = self._track_data[-1,self._mask]
        else:
            r = self._track_function(t)
            
        return r*[L_sun,R_sun,1.0]