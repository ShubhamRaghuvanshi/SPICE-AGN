import numpy as np
from yt.units import cm
from yt.utilities.physical_constants import mp

# --- alias AMR cell properties as particle fields on 'gas' ---

def _gas_particle_position_x(field, data):
    return data[('index','x')]  # cell centers

def _gas_particle_position_y(field, data):
    return data[('index','y')]

def _gas_particle_position_z(field, data):
    return data[('index','z')]

def _gas_particle_mass(field, data):
    # mass per cell = rho * cell_volume
    rho = data[('gas','density')].to('g/cm**3')
    V   = data[('gas','cell_volume')].to('cm**3')
    return (rho * V).to('g')

def _gas_particle_velocity_x(field, data):
    return data[('gas','velocity_x')]  # optional, handy for filtering

def _gas_particle_velocity_y(field, data):
    return data[('gas','velocity_y')]

def _gas_particle_velocity_z(field, data):
    return data[('gas','velocity_z')]

def register_gas_particle_aliases(ds):
    ds.add_field(('gas','particle_position_x'), function=_gas_particle_position_x,
                 units='cm', sampling_type='particle', force_override=True)
    ds.add_field(('gas','particle_position_y'), function=_gas_particle_position_y,
                 units='cm', sampling_type='particle', force_override=True)
    ds.add_field(('gas','particle_position_z'), function=_gas_particle_position_z,
                 units='cm', sampling_type='particle', force_override=True)

    ds.add_field(('gas','particle_mass'), function=_gas_particle_mass,
                 units='g', sampling_type='particle', force_override=True)

    # optional velocities
    ds.add_field(('gas','particle_velocity_x'), function=_gas_particle_velocity_x,
                 units='cm/s', sampling_type='particle', force_override=True)
    ds.add_field(('gas','particle_velocity_y'), function=_gas_particle_velocity_y,
                 units='cm/s', sampling_type='particle', force_override=True)
    ds.add_field(('gas','particle_velocity_z'), function=_gas_particle_velocity_z,
                 units='cm/s', sampling_type='particle', force_override=True)
