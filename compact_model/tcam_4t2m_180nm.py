import numpy as np


tech_param={'v_t':0.4666,
            'k':2000e-6,
            'i_leak':8.4e-23,
            'lam':0.10,
            'm':0.55,
            'q':-0.15,
           'r_lrs':5e3,
           'r_hrs':5e6,
           'c_par':2e-15,
           'r_par':1.4,
           'sigma_g_hrs':0.05,
           'sigma_g_lrs':0.001}

circuit_param={'v_wl':1.5,
               'v_sl':0.75,
               'v_ml':0.6,
               'v_ns':0.1,
               'v_sense':0.3,
               't_clk':100e-9,
               'c_pc':50e-15,
               'c_sense':2e-12}


def i_ds(v_g,v_ds):
    v_t = tech_param['v_t']
    k = tech_param['k']
    lam = tech_param['lam']
    m = tech_param['m']
    q = tech_param['q']
    i_leak = tech_param['i_leak']
    if v_g>v_t:
        v_ds_sat_now = v_g*m+q
        if v_ds<v_ds_sat_now:
            # Linear region
            i_ds = k*((v_ds_sat_now)*v_ds-v_ds**2/2)
        else:
            i_ds = k/2*(v_ds_sat_now)**2*(1+lam*(v_ds-(v_ds_sat_now)))   
    else:
        i_ds=i_leak
    return i_ds

def row(store_vect,search_vect,returnPerf=False,verbose=False,r_low=0,r_high=0,t_clk=0,v_ml=0,v_sense=0):
    v_t = tech_param['v_t']
    k = tech_param['k']
    lam = tech_param['lam']
    m = tech_param['m']
    q = tech_param['q']
    if r_high==0:
        r_high = tech_param['r_hrs']
    if r_low==0:
        r_low = tech_param['r_lrs']
    c_par = tech_param['c_par']
    
    v_sl=circuit_param['v_sl']
    v_ns=circuit_param['v_ns']
    v_wl=circuit_param['v_wl']
    if t_clk==0:
        t_clk=circuit_param['t_clk']
    if v_ml ==0:
        v_ml=circuit_param['v_ml']
    c_pc=circuit_param['c_pc']
    if v_sense==0:
        v_sense=circuit_param['v_sense']
    c_sense=circuit_param['c_sense']

    
    
    width=len(store_vect)
    
    c_ml = c_pc + width*c_par + c_sense
    R_mos=k*(v_wl*m+q)

    v_match=[]
    for store,search in zip(store_vect,search_vect):
        if store==2 or search==2:
            v_match.append(v_sl*(R_mos+r_high)/(2*R_mos+r_high+r_high))
        elif store==search and store!=2:
            v_match.append(v_sl*(R_mos+r_low)/(2*R_mos+r_high+r_low))
        else:
            v_match.append(v_sl*(R_mos+r_high)/(2*R_mos+r_high+r_low))
    if verbose:
        print(v_match)
    i_ml=0
    for n in range(len(v_match)):
        i_ml=i_ml+i_ds(v_match[n]-v_ns,v_ml)
    
    # sense amplifier
    v_ml_final=v_ml-t_clk*i_ml/c_ml
    if v_ml_final<v_ns:
        v_ml_final=v_ns
    ml=1
    if v_ml_final<v_sense:
        ml=0
    
    if returnPerf==False:
        return i_ml,v_ml_final,ml    
    else:
        # dynamic power to discharge ML
        p_charge=0.5*(c_ml*(circuit_param["v_ml"]**2)/circuit_param["t_clk"])*width
        p_discharge=0.5*(c_ml*(circuit_param["v_ml"]-v_ml_final)**2/circuit_param["t_clk"])*width
        p_ml=p_charge+p_discharge
        
        # static power in the voltage divider
        p_static=width*(2*R_mos+r_high+r_low)*v_sl**2
        
        return i_ml,v_ml_final,ml,p_ml,p_static
    
def array(store_mat,search_vect,returnPerf=False,r_low=0,r_high=0,t_clk=0,v_ml=0,v_sense=0,v_wl=-1):
    width=store_mat.shape[1]
    height=store_mat.shape[0]
    c_ml = circuit_param['c_pc'] + width*tech_param['c_par'] + circuit_param['c_sense']
    c_wl=tech_param['c_par']
    if v_wl == -1:
        v_wl=circuit_param["v_wl"]
    if t_clk == 0:
        t_clk=circuit_param["t_clk"]

    result=[]
    for store_vect in store_mat:
        result.append(row(store_vect,search_vect,returnPerf=False,r_low=r_low,r_high=r_high,t_clk=t_clk,v_ml=v_ml,v_sense=v_sense))
    if returnPerf:
        p_static = np.sum(result[:][4])
        p_ml = np.sum(result[:][3])
        # charging and discharging WL
        p_wl=0.5*(c_ml*(v_wl**2)/t_clk)*width*height
        
        return result,p_ml,p_wl,p_static
    else:
        return result

    