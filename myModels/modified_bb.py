import astropy.units as u
import astropy.constants as ct
from astropy.modeling import Fittable1DModel, Parameter
#from astropy.modeling.blackbody import blackbody_nu

class ModifiedBB(Fittable1DModel):
    STD_UNITS = {'size_c':u.sr, 'size_w':u.sr, 'size_h':u.sr, 'Tc':u.K,
            'Tw':u.K, 'Th':u.K, 'nu0':u.GHz, 'Nw':1/u.cm**2, 'Nh':1/u.cm**2, 
            'x':u.GHz, 'y':u.Jy}
    inputs = ('nu',)
    outputs = ('F',)
    
    # Cold
    size_c = Parameter()
    Tc = Parameter()
    nu0 = Parameter()
    beta = Parameter()
    # Warm
    size_w = Parameter()
    Tw = Parameter()
    Nw = Parameter()
    # Hot
    size_h = Parameter()
    Th = Parameter()
    Nh = Parameter()

    def __init__(self, dust, **params):
        self.dust = dust
        size_c = params['size_c']
        Tc = params['Tc']
        nu0 = params['nu0']
        beta = params['beta']
        size_w = params['size_w']
        Tw = params['Tw']
        Nw = params['Nw']
        size_h = params['size_h']
        Th = params['Th']
        Nh = params['Nh']
        super(ModifiedBB, self).__init__(size_c, Tc, nu0, beta, size_w, Tw, Nw,
                size_h, Th, Nh)

    @staticmethod
    def blackbody_nu(nu, T):
        aux1 = 2.*ct.h*nu**3 / ct.c**2
        aux1 = aux / u.sr
        aux2 = np.exp(ct.h*nu/(ct.k_B*T))
        return aux1 / (aux2 - 1.)

    @staticmethod
    def cold_component(x, **params):
        """Modified black body cold component"""
        cold = params['size_c']*STD_UNITS['size_c'] * \
                blackbody_nu(x*STD_UNITS['x'], params['Tc']*STD_UNITS['Tc']) * \
                (1. - np.exp(-(x/params['nu0'])**params['beta']))
        return cold.to(STD_UNITS['y']).value

    @staticmethod
    def warm_component(x, rdg=1.E-2, mu=2.33, **params):
        """Modified black body warm component"""
        sigma = rdg * mu * ct.m_p * self.dust.interpolate(x*STD_UNITS['x'])
        sigma = sigma.to(1/STD_UNITS['Nw']).value
        warm = params['size_w']*STD_UNITS['size_w'] * \
                blackbody_nu(x*STD_UNITS['x'], params['Tw']*STD_UNITS['Tw']) * \
                (1. - np.exp(-params['Nw']*sigma)) * \
                np.exp(-0.5 * (x/params['nu0'])**params['beta'])
        return warm.to(STD_UNITS['y']).value

    @staticmethod
    def hot_component(x, rdg=1.e-2, mu=2.33, **params):
        """Modified black body hot component"""
        sigma = rdg * mu * ct.m_H * self.dust.interpolate(x*STD_UNITS['x'])
        sigma = sigma.to(1/STD_UNITS['Nh']).value
        hot = params['size_h']*STD_UNITS['size_h'] * \
                blackbody_nu(x*STD_UNITS['x'], params['Th']*STD_UNITS['Th']) * \
                (1. - np.exp(-params['Nh']*sigma)) * \
                np.exp(-0.5*((x/params['nu0'])**params['beta']+params['Nw']*sigma))
        return hot.to(STD_UNITS['y']).value

    @staticmethod
    def deriv_planck(x, T):
        aux1 = 2.*ct.h**2*(x*STD_UNITS['x'])**4 / \
                ((T*STD_UNITS['Tc'])**2*ct.c**2*ct.k_B)
        aux2 = np.exp(ct.h*(x*STD_UNITS['x'])/(ct.k_B*T*STD_UNITS['Tc']))
        deriv = aux1 * aux2/(aux2 - 1.)**2
        deriv = deriv.decompose()
        return deriv.to(STD_UNITS['y']/STD_UNITS['Tc']).value

    @staticmethod
    def evaluate(x, size_c, Tc, nu0, beta, size_w, Tw, Nw, size_h, Th, Nh):
        cold = cold_component(x, size_c=size_c, Tc=Tc, nu0=nu0, beta=beta) 
        warm = warm_component(x, size_w=size_w, Tw=Tw, Nw=Nw, 
                nu0=nu0, beta=beta)
        hot = hot_component(x, size_h=size_h, Th=Th, Nh=Nh, 
                nu0=nu0, beta=beta, Nw=Nw)
        return cold + warm + hot

    @staticmethod
    def fit_deriv(x, size_c, Tc, nu0, beta, size_w, Tw, Nw, size_h, Th, Nh):
        # Components
        cold = cold_component(x, size_c=size_c, Tc=Tc, nu0=nu0, beta=beta)
        warm = warm_component(x, size_w=size_w, Tw=Tw, Nw=Nw,
                nu0=nu0, beta=beta)
        hot = hot_component(x, size_h=size_h, Th=Th, Nh=Nh,
                nu0=nu0, beta=beta, Nw=Nw)
        sigma = 1.E-2 * 2.33 * ct.m_p * self.dust.interpolate(x*STD_UNITS['x'])
        sigma = sigma.to(1/STD_UNITS['Nc']).value

        # Blackbody functions
        bbcold = blackbody_nu(x*STD_UNITS['x'], Tc*STD_UNITS['Tc'])
        bbcold = bbcold.to(STD_UNITS['y']/STD_UNITS['size_c']).value
        bbwarm = blackbody_nu(x*STD_UNITS['x'], Tc*STD_UNITS['Tw'])
        bbwarm = bbcold.to(STD_UNITS['y']/STD_UNITS['size_w']).value
        bbhot = blackbody_nu(x*STD_UNITS['x'], Tc*STD_UNITS['Th'])
        bbhot = bbcold.to(STD_UNITS['y']/STD_UNITS['size_h']).value

        d_size_c =  cold / size_c
        d_size_w =  warm / size_w
        d_size_h =  hot / size_h
        d_Tc = cold * deriv_planck(x, Tc) / bbcold
        d_Tw = warm * deriv_planck(x, Tw) / bbwarm 
        d_Th = hot * deriv_planck(x, Th) / bbhot 
        d_beta = -size_c*bbcold * np.exp(-(x/nu0)**beta) * \
                (x/n0)**beta * np.log(x/nu0)
        d_beta += -0.5*(x/nu0)**beta * np.log(x/nu0) * (warm+hot)
        d_nu0 = -beta*(x/nu0)**beta * cold / (nu0*(np.exp((x/nu0)**beta)-1.))
        d_nu0 += 0.5*beta*(x/nu0)**beta / nu0 * (warm+hot)
        d_Nw = sigma * warm / (np.exp(sigma*Nw) - 1.)
        d_Nw += -0.5*sigma * hot 
        d_Nh = sigma * hot / (np.exp(sigma*Nh) - 1.)

        return [d_size_c, d_Tc, d_nu0, d_beta, d_size_w, d_Tw, d_Nw, d_size_h,
                d_Th, d_Nh]

