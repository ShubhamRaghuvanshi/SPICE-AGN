import numpy as np
from yt.units import cm
from yt.utilities.physical_constants import mp

# ---- composition knob ----
Y_HE = 0.24   # helium mass fraction (adjust if needed)

def _get_XH(data):
    Z = data[('ramses','Metallicity')] if ('ramses','Metallicity') in data.ds.field_list else 0.0
    X_H = (1.0 - Y_HE - Z).clip(0.0, 1.0)
    return X_H

def _xHII(data):
    return data[('ramses','hydro_scalar_02')]

def _xHeII(data):
    return data[('ramses','hydro_scalar_03')] if ('ramses','hydro_scalar_03') in data.ds.field_list else data.ds.arr(0.0, '')

def _xHeIII(data):
    return data[('ramses','hydro_scalar_04')] if ('ramses','hydro_scalar_04') in data.ds.field_list else data.ds.arr(0.0, '')

# --- base element number densities ---
def _nH(field, data):
    rho = data[('gas','density')].to('g/cm**3')          # g/cm^3
    X_H = _get_XH(data)
    return (rho * X_H / mp).to('cm**-3')

def _nHe(field, data):
    rho = data[('gas','density')].to('g/cm**3') 
    return (rho * (Y_HE/4.0) / mp).to('cm**-3')

# --- species number densities ---
def _nHI(field, data):
    nH = data[('gas','nH')]
    xHI = (1.0 - _xHII(data)).clip(0.0, 1.0)
    return (xHI * nH).to('cm**-3')

def _nHII(field, data):
    nH = data[('gas','nH')]
    return (_xHII(data) * nH).to('cm**-3')

def _nHeII(field, data):
    nHe = data[('gas','nHe')]
    return (_xHeII(data) * nHe).to('cm**-3')

def _nHeIII(field, data):
    nHe = data[('gas','nHe')]
    return (_xHeIII(data) * nHe).to('cm**-3')

def _nHeI(field, data):
    nHe = data[('gas','nHe')]
    xHeI = (1.0 - (_xHeII(data) + _xHeIII(data))).clip(0.0, 1.0)
    return (xHeI * nHe).to('cm**-3')

# --- species mass DENSITIES (g/cm^3) ---
def _rhoH(field, data):
    rho = data[('gas','density')].to('g/cm**3') 
    X_H = _get_XH(data)
    return rho * X_H

def _rhoHe(field, data):
    rho = data[('gas','density')].to('g/cm**3') 
    return rho * Y_HE

def _rhoHI(field, data):
    rho = data[('gas','density')].to('g/cm**3') 
    X_H = _get_XH(data)
    xHI = (1.0 - _xHII(data)).clip(0.0, 1.0)
    return rho * (X_H * xHI)

def _rhoHII(field, data):
    rho = data[('gas','density')].to('g/cm**3') 
    X_H = _get_XH(data)
    return rho * (X_H * _xHII(data))

def _rhoHeI(field, data):
    rho = data[('gas','density')].to('g/cm**3') 
    xHeI = (1.0 - (_xHeII(data) + _xHeIII(data))).clip(0.0, 1.0)
    return rho * (Y_HE * xHeI)

def _rhoHeII(field, data):
    rho = data[('gas','density')].to('g/cm**3') 
    return rho * (Y_HE * _xHeII(data))

def _rhoHeIII(field, data):
    rho = data[('gas','density')].to('g/cm**3') 
    return rho * (Y_HE * _xHeIII(data))

# --- species cell MASSES (g) ---

def _m_from_rho(field, data, rho_field):
    rho  = data[rho_field].to('g/cm**3')
    V    = data[('gas','cell_volume')].to('cm**3')  # force proper cm^3
    mass = (rho * V).to('g')
    return mass.to('Msun')


def _mH(field, data):     return _m_from_rho(field, data, ('gas','rhoH'))
def _mHe(field, data):    return _m_from_rho(field, data, ('gas','rhoHe'))
def _mHI(field, data):    return _m_from_rho(field, data, ('gas','rhoHI'))
def _mHII(field, data):   return _m_from_rho(field, data, ('gas','rhoHII'))
def _mHeI(field, data):   return _m_from_rho(field, data, ('gas','rhoHeI'))
def _mHeII(field, data):  return _m_from_rho(field, data, ('gas','rhoHeII'))
def _mHeIII(field, data): return _m_from_rho(field, data, ('gas','rhoHeIII'))

# ---------- mass fraction fields (dimensionless) ----------

def _Z_massfrac(data):
    # Z may be missing; default to 0 with correct units/shape
    if ('ramses', 'Metallicity') in data.ds.field_list:
        Z = data[('ramses', 'Metallicity')]
    else:
        # make a dimensionless array shaped like density
        rho = data[('gas','density')]
        Z = (rho/rho) * 0.0
    return Z.in_units('')

def _X_H(field, data):
    return _get_XH(data).in_units('')  # hydrogen mass fraction

def _Y_He(field, data):
    # broadcast constant Y_HE to a field-shaped, dimensionless array
    rho = data[('gas','density')]
    return ((rho/rho) * Y_HE).in_units('')

def _Z(field, data):
    # metals mass fraction
    XH  = _get_XH(data)
    YHe = _Y_He(field, data)
    Z   = (1.0 - XH - YHe).clip(0.0, 1.0)
    return Z.in_units('')

# Ionization-state *mass* fractions (elemental mass fraction × ionization fraction)
def _X_HI(field, data):
    XH  = _get_XH(data)
    xHI = (1.0 - _xHII(data)).clip(0.0, 1.0)
    return (XH * xHI).in_units('')

def _X_HII(field, data):
    XH = _get_XH(data)
    return (XH * _xHII(data)).in_units('')

def _Y_HeI(field, data):
    rho = data[('gas','density')]
    YHe = (rho/rho) * Y_HE
    xHeI = (1.0 - (_xHeII(data) + _xHeIII(data))).clip(0.0, 1.0)
    return (YHe * xHeI).in_units('')

def _Y_HeII(field, data):
    rho = data[('gas','density')]
    YHe = (rho/rho) * Y_HE
    return (YHe * _xHeII(data)).in_units('')

def _Y_HeIII(field, data):
    rho = data[('gas','density')]
    YHe = (rho/rho) * Y_HE
    return (YHe * _xHeIII(data)).in_units('')

def _mass_fraction_sum(field, data):
    # Should be ≈ 1: (HI+HII) + (HeI+HeII+HeIII) + Z
    return (
        data[('gas','X_HI')] + data[('gas','X_HII')] +
        data[('gas','Y_HeI')] + data[('gas','Y_HeII')] + data[('gas','Y_HeIII')] +
        data[('gas','Z')]
    ).in_units('')

def register_number_density_and_mass_fields(ds):
    # ----- MASS FRACTIONS (dimensionless) -----
    ds.add_field(('gas','X_H'),      function=_X_H,      units='', sampling_type='cell', force_override=True)
    ds.add_field(('gas','Y_He'),     function=_Y_He,     units='', sampling_type='cell', force_override=True)
    ds.add_field(('gas','Z'),        function=_Z,        units='', sampling_type='cell', force_override=True)

    ds.add_field(('gas','X_HI'),     function=_X_HI,     units='', sampling_type='cell', force_override=True)
    ds.add_field(('gas','X_HII'),    function=_X_HII,    units='', sampling_type='cell', force_override=True)
    ds.add_field(('gas','Y_HeI'),    function=_Y_HeI,    units='', sampling_type='cell', force_override=True)
    ds.add_field(('gas','Y_HeII'),   function=_Y_HeII,   units='', sampling_type='cell', force_override=True)
    ds.add_field(('gas','Y_HeIII'),  function=_Y_HeIII,  units='', sampling_type='cell', force_override=True)

    # Optional sanity check: sum of all tracked mass fractions
    ds.add_field(('gas','mass_fraction_sum'), function=_mass_fraction_sum, units='', sampling_type='cell', force_override=True)

    # ----- NUMBER DENSITIES -----
    ds.add_field(('gas','nH'),     function=_nH,     units='cm**-3', sampling_type='cell', force_override=True)
    ds.add_field(('gas','nHe'),    function=_nHe,    units='cm**-3', sampling_type='cell', force_override=True)
    ds.add_field(('gas','nHI'),    function=_nHI,    units='cm**-3', sampling_type='cell', force_override=True)
    ds.add_field(('gas','nHII'),   function=_nHII,   units='cm**-3', sampling_type='cell', force_override=True)
    ds.add_field(('gas','nHeI'),   function=_nHeI,   units='cm**-3', sampling_type='cell', force_override=True)
    ds.add_field(('gas','nHeII'),  function=_nHeII,  units='cm**-3', sampling_type='cell', force_override=True)
    ds.add_field(('gas','nHeIII'), function=_nHeIII, units='cm**-3', sampling_type='cell', force_override=True)

    # ----- MASS DENSITIES -----
    ds.add_field(('gas','rhoH'),     function=_rhoH,     units='g/cm**3', sampling_type='cell', force_override=True)
    ds.add_field(('gas','rhoHe'),    function=_rhoHe,    units='g/cm**3', sampling_type='cell', force_override=True)
    ds.add_field(('gas','rhoHI'),    function=_rhoHI,    units='g/cm**3', sampling_type='cell', force_override=True)
    ds.add_field(('gas','rhoHII'),   function=_rhoHII,   units='g/cm**3', sampling_type='cell', force_override=True)
    ds.add_field(('gas','rhoHeI'),   function=_rhoHeI,   units='g/cm**3', sampling_type='cell', force_override=True)
    ds.add_field(('gas','rhoHeII'),  function=_rhoHeII,  units='g/cm**3', sampling_type='cell', force_override=True)
    ds.add_field(('gas','rhoHeIII'), function=_rhoHeIII, units='g/cm**3', sampling_type='cell', force_override=True)

    # ----- CELL MASSES -----
    ds.add_field(('gas','mH'),     function=_mH,     units='Msun', sampling_type='cell', force_override=True)
    ds.add_field(('gas','mHe'),    function=_mHe,    units='Msun', sampling_type='cell', force_override=True)
    ds.add_field(('gas','mHI'),    function=_mHI,    units='Msun', sampling_type='cell', force_override=True)
    ds.add_field(('gas','mHII'),   function=_mHII,   units='Msun', sampling_type='cell', force_override=True)
    ds.add_field(('gas','mHeI'),   function=_mHeI,   units='Msun', sampling_type='cell', force_override=True)
    ds.add_field(('gas','mHeII'),  function=_mHeII,  units='Msun', sampling_type='cell', force_override=True)
    ds.add_field(('gas','mHeIII'), function=_mHeIII, units='Msun', sampling_type='cell', force_override=True)
