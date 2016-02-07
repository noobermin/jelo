import numpy as np;

e0 = 8.854e-12;
c  = 2.998e8;
e  = 1.602e-19;
m_e= 9.10938291e-31;


sq = lambda v: np.diagonal(np.dot(v.T,v));
az = lambda v: np.arctan2(v[1], v[0]);

a0 = lambda I, lm=800e-9: -e*(
    np.sqrt(I*1e4*2/e0/c) /
    (m_e*c*2*np.pi*c/lm)
);

gm0 = lambda a0, rho: np.sqrt(1+np.dot(rho,rho)*a0**2);
Rs  = lambda a0, rho: gm0(a0,rho) - abs(a0)*rho[0];
#garbage
sin_pm    = lambda h,al,sign: np.sin(h*(2*al+sign)/(2*al));
cos_pm    = lambda h,al,sign: np.cos(h*(2*al+sign)/(2*al));
sin_minus = lambda h,al: sin_pm(h,al,-1.0)
sin_plus  = lambda h,al: sin_pm(h,al, 1.0)
def f1_h0(h,al):
    return sin_minus(h,al)*al/(2*al-1) - sin_plus(h,al)*al/(2*al+1)
def f1(h,al):
    return f1_h0(h,al)-f1_h0(h[0],al);
cos = np.cos
def f1i(h,al):
    def f1i_term1(h,al):
        return 2*al**2*(cos_pm(h,al,1.0)/(1+2*al)**2-cos_pm(h,al,-1.0)/(2*al-1)**2);
    return -f1_h0(h[0],al)*(h-h[0]) + f1i_term1(h,al) - f1i_term1(h[0],al);

def f2i(h,al):
    def A_pm(h,al,sign):
        return al**2/(2.0*(2*al+sign)**2) *(
            (h-h[0])/2.0-(
                sin_pm(2.0*h,al,sign)-sin_pm(2.0*h[0],al,sign))*al/(2.0*(2*al+sign))
        );
    A = lambda h,al: A_pm(h,al,-1.0) + A_pm(h,al,1.0);
    def B(h,al):
        return al**2/(2.0*(4*al**2-1)) * (np.sin(2*h) - np.sin(2*h[0]) - al*(np.sin(h/al) - np.sin(h[0]/al)));
    C = lambda h,al: -f1_h0(h[0], al)* ( f1i(h,al) + f1_h0(h[0],al) * (h-h[0])/2.0);
    return A(h,al) + B(h,al) + C(h,al);

defal=30e-15/(8e-7/c);
def px(h,a0,rho,al=defal):
    return abs(a0)*rho[0] + (a0**2*f1(h,al)**2/2 +  abs(a0)*rho[1]*a0*f1(h,al))/Rs(a0,rho);
py = lambda h,a0,rho,al=defal: abs(a0)*rho[1] + a0*f1(h,al)
pz = lambda h,a0,rho,al=defal: abs(a0)*rho[2] + np.zeros(h.shape)
p  = lambda h,a0,rho,al=defal: np.array([
    px(h,a0,rho), py(h,a0,rho), pz(h,a0,rho)
]);

ke    = lambda h, a0, rho: np.sqrt(1+sq(p(h,a0,rho))) - 1;
angle = lambda h, a0, rho: az(p(h,a0,rho));

def kx(h, a0, rho, al=defal):
    return (h-h[0])*abs(a0)*rho[0]/Rs(a0,rho) + (a0**2 * f2i(h,al) + a0*f1i(h,al)*abs(a0)*rho[1])/Rs(a0,rho)**2;
ky = lambda h,a0,rho,al=defal: h*abs(a0)*rho[1]/Rs(a0,rho) + a0*f1i(h,al)/Rs(a0,rho);
kz = lambda h,a0,rho: h*abs(a0)*rho[2]/Rs(a0,rho);
wt = lambda h,a0,rho,al=defal: h + kx(h,a0,rho,al);
