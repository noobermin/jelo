import numpy as np;

#beta functions
beta_plus  = lambda u,bz0,phi_s,phi_t: (np.sin(phi_s+phi_t+u*bz0+u)-np.sin(phi_s+phi_t))/(bz0+1)
beta_minus = lambda u,bz0,phi_s,phi_t: (np.sin(phi_s-phi_t+u*bz0-u)-np.sin(phi_s-phi_t))/(bz0-1)
sinoid_beta = lambda u,bz0,phi_s,phi_t: beta_plus(u,bz0,phi_s,phi_t)+beta_minus(u,bz0,phi_s,phi_t);

#v functions
v_plus = lambda u,bz0,phi_s,phi_t: (np.cos(phi_s+phi_t+u*bz0+u)-np.cos(phi_s+phi_t))/(bz0+1)**2;
v_minus = lambda u,bz0,phi_s,phi_t: (np.cos(phi_s-phi_t+u*bz0-u)-np.cos(phi_s-phi_t))/(bz0-1)**2;
v_prime = lambda bz0, phi_s, phi_t: np.sin(phi_s+phi_t)/(bz0+1) + np.sin(phi_s-phi_t)/(bz0-1);

sinuoid_v = lambda u,bz0,phi_s,phi_t: v_plus(u,bz0,phi_s,phi_t)+v_minus(u,bz0,phi_s,phi_t)

xhat = np.array([1,0,0]);
yhat = np.array([0,1,0]);
zhat = np.array([0,0,1]);

def born(x0,b0, # position, speed
         a0,    # a_0
         Es,Et,Bs,Bt):
    b0_sq = np.dot(b0,b0);
    #shifting v to have it start at v=0.
    if not np.isclose(x0[2], 0.0):
        Es -= x0[2];
        Bs -= x0[2];
        x0[2] = 0.0;
    v = x0[2];
    pref = a0*(1-b0_sq);

    e = lambda u: sinoid_beta(u,b0[2],Es,Et)
    b = lambda u: sinoid_beta(u,b0[2],Bs,Bt)

    ebar   = lambda u: sinuoid_v(u,b0[2],Es,Et);
    eprime = v_prime(b0[2],Es,Et);
    bbar   = lambda u: sinuoid_v(u,b0[2],Bs,Bt);
    bprime = v_prime(b0[2],Bs,Bt);

    print(pref,eprime);
    def _beta(u):
        return b0+pref*( yhat*e(u) -b0[2]*e(u)*b0 + np.cross(b0,xhat)*b(u) );
    def _v(u):
        print(pref*(yhat)*(ebar(u)-u*eprime));
        return x0 + b0*u + pref*( (yhat-b0[2]*b0)*(ebar(u)-u*eprime) +np.cross(b0,xhat)*(bbar(u)-u*bprime));
    return _beta,_v
