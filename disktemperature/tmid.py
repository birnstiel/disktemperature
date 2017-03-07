##/usr/bin/env python
import numpy as np
import constants as _c

def t_function(Td,r,sig_g,sig_d,alpha,phi,M_star,T_star,L_star,kappa_r,kappa_p,T0):
    """
    Objective function for the temperature calculation. Solve equation
    t_function(Td)==0 to find Td.

    Arguments:
    ----------
    r : float
    : radial position [cm]

    sig_g,sig_d : float
    : gas and dust surface densities [g cm^-2]

    alpha : float
    : turbulence parameter

    phi :float
    : irradiation angle [rad]

    M_star,T_star,L_star : floats
    : stellar mass, temperature, luminosity [cgs]

    kappa_r,kappa_p : functions
    : scalar functions with temperature as only argument - returning the
      Rosseland and Planck mean opacities at that temperature, respectively
      
    T0 : float
    : minimum temperature

    Returns:
    --------
    objective function value, find the root of it to get the temperature
    """

    cs = np.sqrt(_c.k_b*Td/_c.mu/_c.m_p)
    om = np.sqrt(_c.Grav*M_star/r**3)

    kap_r_dust = kappa_r(Td)
    kap_p_dust = kappa_p(Td)
    kap_p_star = kappa_p(T_star)

    tau_r_dust = 0.5*sig_d*kap_r_dust
    tau_p_dust = 0.5*sig_d*kap_p_dust
    tau_p_star = 0.5*sig_d*kap_p_star

    return -Td**4 \
        + 9./(8.*_c.sig_sb)*(3./8.*tau_r_dust+0.5/tau_p_dust)*sig_g*alpha*cs**2*om \
        + L_star/(2.*_c.sig_sb*4*np.pi*r**2)* \
            ( kap_p_star/kap_p_dust*0.5*np.exp(-tau_p_star/np.sin(phi)) \
            + phi ) \
        + T0**4

class tmid:

        # these will be initialized in `__init__`

    kappa_abs = None
    kappa_sca = None
    rho_s     = None
    a         = None
    lam_mic   = None
    T_lookup  = None

    # all of those are common keywords, but need to have defaults

    T_star    = None
    R_star    = None
    L_star    = None

    phi       = 0.05 # irradiation angle
    Tmin      = 7.   # minimum temperature

    pseudo_gas_opacity = 0.001

        # those are the hidden attributes, that influence the initialization should
        # be updated to be "properties" at some point

    _t_max_dust = 1500.
    _delta_temp = 200.

    _amin     = 0.1
    _amax     = 1e5
    _na       = 100

    _Tmin     = 0.1
    _Tmax     = 10.
    _nT       = 110

    _lmin     = 0.1e-4
    _lmax     = 9.
    _nl       = 150

    sim       = None

    def __init__(self,**kwargs):
        """
        Initialize the opacities of the module.
        By default, the module is initialized with the default values from above.
        You can reinitialize by calling the function again with other parameters.

        Keywords:
        ---------

        All class variables can be given:

        _amin,_amax,_lmin,_lmax : float
        : minimum/maximum of the particle size (in cm)/wavelength (in micron)

        _na,_nl : int
        : number of grain sizes / wavelengths

        vol_fract : list
        : volume fractions which for silicate, carbon, water ice,
          and vacuum are by default [0.07,0.21,0.42,0.30] as in Ricci et al. 2010a.
        """
        import bhmie
        from uTILities import planck_dBnu_dT, planck_B_nu

        # update stellar parameters

        kwargs = self.update_stellar_parameters(doreturn=True,**kwargs)

        #
        # update global defaults and calculate size and wavelength arrays
        #
        for k,v in kwargs.items():
            if hasattr(self,k): setattr(self,k,v)
        self.a        = np.logspace(np.log10(self._amin),np.log10(self._amax),self._na)
        self.lam_mic  = np.logspace(np.log10(self._lmin)+4,np.log10(self._lmax)+4,self._nl)
        self.T_lookup = np.logspace(-1,np.log10(self._t_max_dust+self._delta_temp+10.),self._nT)

        #
        # calculate opacities using Lucas materials mix by default
        # default volume fractions are from Lucas thesis,
        # the fractions in Ricci+2010 are typos
        #
        vol_fract   = kwargs.pop('vol_fract',[0.07,0.21,0.42,0.30])
        c1          = bhmie.diel_luca('silicate',extrapol=True,lmax=self.lam_mic[-1]*1e-4,lmin=self.lam_mic[0]*1e-4)
        c2          = bhmie.diel_luca('carbon',  extrapol=True,lmax=self.lam_mic[-1]*1e-4,lmin=self.lam_mic[0]*1e-4)
        c3          = bhmie.diel_luca('ice',     extrapol=True,lmax=self.lam_mic[-1]*1e-4,lmin=self.lam_mic[0]*1e-4)
        c4          = bhmie.diel_vacuum()
        constants   = [c1,c2,c3,c4]
        densities   = [3.50,2.50,0.92,0.00]
        self.rho_s  = sum(densities*np.array(vol_fract))
        mix         = bhmie.diel_mixed(constants,vol_fract,rule='Bruggeman') # these are the mixed opacities
        m           = 4*np.pi/3.*self.rho_s*self.a**3
        q_abs,q_sca = bhmie.get_mie_coefficients(self.a,self.lam_mic*1e-4,mix)

        self.kappa_abs   = q_abs*np.tile(np.pi*self.a**2/m,[self._nl,1]).transpose()
        self.kappa_sca   = q_sca*np.tile(np.pi*self.a**2/m,[self._nl,1]).transpose()

        # convert to frequency sorting

        self.lam_mic   = self.lam_mic[::-1]
        self.freq      = _c.c_light/(self.lam_mic*1e-4)
        self.kappa_abs = self.kappa_abs[:,::-1]
        self.kappa_sca = self.kappa_sca[:,::-1]

        # calculate the opacity averaging arrays

        self.Bnu        = np.array([[planck_B_nu(f,T) for f in self.freq] for T in self.T_lookup])
        self.Bnu_int    = np.trapz(self.Bnu,x=self.freq)

        self.dBnudT     = np.array([[planck_dBnu_dT(f,T) for f in self.freq] for T in self.T_lookup])
        self.dBnudT_int = np.trapz(self.dBnudT,x=self.freq)

    def update_stellar_parameters(self,doreturn=False,**kwargs):
        """
        Updates the stellar properties by passing the keywords M_star,
        R_star, L_star, T_star, all in CGS units.

        M_star can be set individually. The other quantities need to updated
        by giving 2 out of the 3 quantities. The third one is updated
        according to

        sig_sb*T_star**4 = L_star/(4*np.pi*R_star**2)

        """
        M_star = kwargs.pop('M_star',None)
        T_star = kwargs.pop('T_star',None)
        R_star = kwargs.pop('R_star',None)
        L_star = kwargs.pop('L_star',None)

        if M_star is not None: self.M_star = M_star

        n_star_prop = sum([i is None for i in [T_star,R_star,L_star]])

        if n_star_prop > 0:
            if n_star_prop != 1:
                raise ValueError('Two stellar properties (T,R,L) need to be given, the third will be calculated')
            else:

                # calculate the missing quantity

                if T_star is None:
                    T_star = (L_star/(4*np.pi*R_star**2*_c.sig_sb))**0.25
                if L_star is None:
                    L_star = _c.sig_sb*T_star**4*4*np.pi*R_star**2
                if R_star is None:
                    R_star = (L_star/(4*np.pi*_c.sig_sb*T_star**4))**0.5

                # update all

                self.T_star = T_star
                self.R_star = R_star
                self.L_star = L_star

        if doreturn: return kwargs

    def get_t_mid(self,r,sig_g,sig_d,alpha,phi=None):
        """
         This subroutine updates the Temperature array

        r : array
        :   the grid
        
        sig_g , sig_d : array
        :   the total gas, dust surface density on grid r
        
        alpha : array
        :   turbulence alpha on grid r
        
        Keywords:
        ---------
        
        phi : float | array
        :   irradiation angle (constant or r-dependent), if None, the initialized value is taken 
        
        """
        from scipy.optimize import brentq
        
        if phi is None: phi = self.phi
        phi = phi*np.ones(len(r))

        # average the opacities of the grain sizes
        # this will be of shape (n_freq,n_r)

        sig_d_total = sig_d.sum(0)
        kappa_mean = (self.kappa_abs[:,:,np.newaxis]*sig_d[:,np.newaxis,:]).sum(0)/sig_d_total

        # update the Planck & Rosseland mean opacity lookup tables

        self.update_kappa_tables(kappa_mean)

        # solve the temperature equation at every radius

        T = np.zeros(len(r))
        for i in range(len(r)):
            k_r = lambda T: self.kappa_r(T,i)
            k_p = lambda T: self.kappa_p(T,i)
            T[i]  = brentq(t_function,0.1,1e4,args=(r[i],sig_g[i],sig_d_total[i],alpha[i],phi[i],self.M_star,self.T_star,self.L_star,k_r,k_p,self.Tmin))

        return T

    def _kappa_mean(self,temp,kappa):
        """
        This function interpolates amean opacity from the lookup tables
        and includes an arbitrary extrapolation to higher temperatures
        to mimick sublimation.

        Arguments:
        ----------

        temp : array
        : the temperature for which the opacity is searched

        kappa : array
        :    the rosseland or planck opacity array shape

        Output:
        -------

        k : float
        : the mean opacity
        """

        # interpolate the temperature in the temperature table

        k = np.interp(temp,self.T_lookup,kappa,left=kappa[0],right=kappa[1])

        # treat evaporation (interpolate down to the gas opacity)

        if(temp > self._t_max_dust):
            if(temp > self._t_max_dust+self._delta_temp):
                k = self.pseudo_gas_opacity
            else:
                k = (k+self.pseudo_gas_opacity)/2. * \
                        np.cos(np.pi/self._delta_temp*(temp-self._t_max_dust)) + \
                        (k+self.pseudo_gas_opacity)/2.
        return k

    def kappa_p(self,temp,i):
        """
        Returns the Planck mean opacity at the given temperature

        Arguments:
        ----------

        temp : float
        : Temperature [K]

        i : int
        : the radial index

        Returns:
        --------

        kappa_p : float
        : the Planck mean opacity at temperature `temp`
        """
        return self._kappa_mean(temp,self._kappa_p_lookup[:,i])

    def kappa_r(self,temp,i):
        """
        Returns the Rosseland mean opacity at the given temperature

        Arguments:
        ----------

        temp : float
        : Temperature [K]

        i : int
        : the radial index

        Returns:
        --------

        kappa_r : float
        : the Rosseland mean opacity at temperature `temp`
        """
        return self._kappa_mean(temp,self._kappa_r_lookup[:,i])

    def update_kappa_tables(self,kappa_mean):
        """
        Update the Planck and Rosseland mean opacitiy lookup tables
        for a given frequency dependent opacity array.

        Arguments:
        ----------

        kappa_mean : array
        :   absorption opacity at the given frequencies

        """

        self._kappa_p_lookup = np.trapz(kappa_mean[np.newaxis,:,:]*self.Bnu[:,:,np.newaxis],x=self.freq,axis=1)/self.Bnu_int[:,np.newaxis]
        self._kappa_r_lookup= np.trapz(1./kappa_mean[np.newaxis,:,:]*self.dBnudT[:,:,np.newaxis],x=self.freq,axis=1)/self.dBnudT_int[:,np.newaxis]



def temperature_iterator(r,T_of_phi,H_of_T,phi=0.05,n_i=30,n_poly=10,phi_min=0.01,do_plot=False,convergence=1e-3):
    """
    Iterates the temperature calculation and the irradiation angle calculation.
    
    Arguments:
    ----------
    
    r : array
    :   radial grid on which T is calculated [cm]
    
    T_of_phi : function
    :   given an irradiation angle(-array), return the calculated temperature
    
    H_of_T : function
    :   given a Temperature(-array), return the calculated scale height, the irradiation
        angle is then determined by

                dH      H
        phi =  ---- -  ---
                dr      r
                
    Keywords:
    ---------
    
    phi : float or array
    :   initial guess for the irradiation angle phi
    
    n_i : int
    :   maximum number of iterations to be done
    
    n_poly : int
    :   order of the polynomial that is used for fitting/smoothing phi
    
    phi_min : float
    :   minimum value for phi
    
    do_plot : bool
    :   whether or not to plot the iteration results
    
    convergence : float
    :   maximum relative change as exit criterion, default: 0.1%
    """
    
    # define the function that fits phi with a smooth curve to avoid oscillations
    
    from scipy.optimize import curve_fit
    import matplotlib.pyplot as plt
    def fitphi(r,phi,n):
        """
        Fit a n-th order poynomial to phi as function of log10(r).
        """

        def func(x, *p):
            return np.poly1d(p)(np.log10(x))

        popt, _ = curve_fit(func, r, phi, p0=np.ones(n))

        return func(r,*popt)
    

    phi    = np.minimum(phi,phi_min)*np.ones(len(r))
    T      = T_of_phi(phi)
    RES    = []

    if do_plot:
        _,axs = plt.subplots(2,sharex=True)
        cols  = plt.cm.get_cmap('Reds')(np.linspace(0,1,n_i+1))

    for i in range(n_i+1):
        phio = phi[:]
        To   = T[:]
        T    = T_of_phi(phi=phi)
        H    = H_of_T(T)
        dHdr = np.diff(H)/np.diff(r)
        dHdr = 0.5*(dHdr[1:]+dHdr[:-1])
        dHdr = np.hstack((dHdr[0],dHdr,dHdr[-1]))
        phi  = dHdr-H/r
        phi  = fitphi(r,phi,n_poly)
        phi  = np.maximum(phi,phi_min)
        phi  = (phi+phio)*0.5
        res  = np.abs(To/T-1).max()

        # plotting part
        
        if do_plot:
            RES += [res]
            if i in [0,n_i]:
                ls = '--'
            else:
                ls = '-'

            axs[0].loglog(r/_c.AU,T,    ls=ls,c=cols[i])
            axs[1].semilogx(r/_c.AU,phi,ls=ls,c=cols[i])
            
        # exit criterion
        
        if i>0 and res<convergence:
            break

    if do_plot:
        axs[0].set_xlabel('r [AU]')
        axs[0].set_ylabel('T [K]')
        axs[1].set_xlabel('r [AU]')
        axs[1].set_ylabel('$\phi$')

        _,ax = plt.subplots()
        ax.semilogy(np.array(RES)*100)
        ax.set_xlabel('iteration #')
        ax.set_ylabel('relative change in %');

    return T