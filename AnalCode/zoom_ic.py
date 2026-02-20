import os
import numpy as np
import yt_funcs  # only used if we need to generate halodata.csv
import yt 
import math 

# ---- user paths (edit these) ----
simdir       = "/u/sraghu/sim/SHUB/b256/bigbox"
region_file = "/u/sraghu/sim/ic/b256/zoom_2/region_spheres.txt"

# ---- example call ----
outnum_1       = 1
outnum_6       = 9
halonum      = 2
rvir_mult   = 3.0
pad_lag =   1.0
boxlen_comov = 100.0     # cMpc/h
Ngridx_IC    = 256

outnum_1_char = str(outnum_1).zfill(5)
outnum_6_char = str(outnum_6).zfill(5)

outdir_1      = simdir+'/output_' + outnum_1_char
outdir_6      = simdir+'/output_' + outnum_6_char
halofile      = outdir_6 + '/info_' + outnum_6_char + "/halodata.csv"


ds_1          = yt.load(outdir_1)
ds_6          = yt.load(outdir_6)

ad_1          = ds_1.all_data()
ad_6          = ds_6.all_data()

# Snapshot (z~6) particle arrays (NPY, 1D)
snap_ids_npy = np.array ( ad_6[('DM', 'particle_identity')] )      
snap_x_npy   = np.array ( ad_6[('DM', 'particle_position_x')] )      
snap_y_npy   = np.array ( ad_6[('DM', 'particle_position_y')] )       
snap_z_npy   = np.array ( ad_6[('DM', 'particle_position_z')] )     

# IC (z_init) particle arrays (NPY, 1D) — same global IDs as snapshot
ic_ids_npy = np.array ( ad_1[('DM', 'particle_identity')] )      
ic_x_npy   = np.array ( ad_1[('DM', 'particle_position_x')] )      
ic_y_npy   = np.array ( ad_1[('DM', 'particle_position_y')] )      
ic_z_npy   = np.array ( ad_1[('DM', 'particle_position_z')] )     

print("Analyzing :", halofile)

#print("snap_x_npy :", min(snap_x_npy), max(snap_x_npy))
#print("snap_y_npy :", min(snap_y_npy), max(snap_y_npy))
#print("snap_z_npy :", min(snap_z_npy), max(snap_z_npy))

#print("ic_x_npy :", min(ic_x_npy), max(ic_x_npy))
#print("ic_y_npy :", min(ic_y_npy), max(ic_y_npy))
#print("ic_z_npy :", min(ic_z_npy), max(ic_z_npy))

if not os.path.exists(halofile):
  print("Halo catalog not found, creating one...")
  yt_funcs.write_halos(simdir, outnum_6)
else:
  print(halofile, "already exists")

# get the halonum-th line (1-based)
ihalo = 0
with open(halofile, 'r') as file:
  for line in file:
    if(ihalo < halonum-1) :
      ihalo = ihalo + 1
    else :
      break

halo_posx         = float(line.split()[1])
halo_posy         = float(line.split()[2])
halo_posz         = float(line.split()[3])  
halo_rvir_ckpc_h  = float(line.split()[4]) #kpccm/h
halo_mass         = float(line.split()[5])
halo_npart        = float(line.split()[6])
halo_z            = float(line.split()[7])
L = 1.0                                 # boxlength in code units 
halo_rvir = halo_rvir_ckpc_h/1000.0/boxlen_comov   # halo rvir in code units 
r_sphere = rvir_mult * halo_rvir        # selection radius at snapshot

#    Rvir_kpc    = halo_rvir * boxlen_comov * 1000.0 * 0.6774 /(1.0+halo_z)  # ckpc -> kpc


print("z, Rvir_physical[kpc], npart  :", halo_z,  halo_rvir_ckpc_h/0.6774/(1.0 + halo_z), halo_npart, halo_rvir_ckpc_h)
print("halo_center :", halo_posx, halo_posy, halo_posz)
print("halo_mass/1e12 [Msun]:", halo_mass/1e12)
print("zoom region radius :", halo_rvir_ckpc_h*rvir_mult/1000.0, r_sphere*boxlen_comov)

#create mask for particle in r_sphere 
dx = (snap_x_npy - halo_posx + 0.5*L) % L - 0.5*L
dy = (snap_y_npy - halo_posy + 0.5*L) % L - 0.5*L
dz = (snap_z_npy - halo_posz + 0.5*L) % L - 0.5*L
r2   = dx*dx + dy*dy + dz*dz
mask = r2 <= (r_sphere*r_sphere)
sphere_pid = snap_ids_npy[mask]
if(len(sphere_pid)==0):
   raise RuntimeError("No particles found within selection radius. Check units/paths.")

'''
#map particle positions to  IC 
xic = []
yic = []
zic = []
for i6 in range(0, len(sphere_pid)) :
  for i0 in range(0, len(ic_ids_npy)) :
    if(ic_ids_npy[i0] == sphere_pid[i6]) :
      xic.append(ic_x_npy[i0])
      yic.append(ic_y_npy[i0])
      zic.append(ic_z_npy[i0])
      break 
if(len(xic) != len(sphere_pid)) :
   raise ValueError("mapping not one to one")
'''   

# --- FAST mapping: selected IDs at z~6 -> indices in IC arrays ---
sphere_pid = np.asarray(sphere_pid, dtype=ic_ids_npy.dtype)

order = np.argsort(ic_ids_npy)          # sort IC IDs once (O(N log N))
ids0s = ic_ids_npy[order]

idx = np.searchsorted(ids0s, sphere_pid)     # vectorized binary search (O(M log N))
ok  = (idx < ids0s.size) & (ids0s[idx] == sphere_pid)

if not np.all(ok):
    missing = sphere_pid[~ok]
    raise ValueError(f"{missing.size} selected IDs missing in IC dump (first few: {missing[:5]})")

idx0 = order[idx]   # positions in the IC arrays

# Grab IC positions (code units)
Qx = ic_x_npy[idx0]
Qy = ic_y_npy[idx0]
Qz = ic_z_npy[idx0]

#print("len(Qx), len(sphere_pid) :", len(Qx), len(sphere_pid))
#print("Qxmin, Qxmax :", min(Qx), max(Qx))
#print("Qymin, Qymax :", min(Qy), max(Qy))
#print("Qzmin, Qzmax :", min(Qz), max(Qz))


# --- Robust COM across PBC (seed unwrap) ---
sx, sy, sz = Qx[0], Qy[0], Qz[0]
# unwrap everything relative to a seed particle
Qx_u0 = (Qx - sx + 0.5*L) % L - 0.5*L
Qy_u0 = (Qy - sy + 0.5*L) % L - 0.5*L
Qz_u0 = (Qz - sz + 0.5*L) % L - 0.5*L

# COM in *unwrapped* space, then wrap back to [0,1)
COM = (np.array([sx, sy, sz]) +
       np.array([Qx_u0.mean(), Qy_u0.mean(), Qz_u0.mean()])) % L

# unwrap once more around this COM to measure extent
Qx_u = (Qx - COM[0] + 0.5*L) % L - 0.5*L
Qy_u = (Qy - COM[1] + 0.5*L) % L - 0.5*L
Qz_u = (Qz - COM[2] + 0.5*L) % L - 0.5*L

# 4) Padded Lagrangian radius
R_raw = np.sqrt(Qx_u*Qx_u + Qy_u*Qy_u + Qz_u*Qz_u).max()  
R_pad = pad_lag * R_raw                                    

# 5) Wrap COM back to [0,L) and convert to box fractions for MUSIC
COM = COM%L
x_frac, y_frac, z_frac = COM
r_frac = R_pad / L

# 6) Write region_spheres.txt (one sphere; MUSIC takes union if you add more lines)
os.makedirs(os.path.dirname(region_file), exist_ok=True)
with open(region_file, "w") as f:
    f.write(f"{x_frac:.10f} {y_frac:.10f} {z_frac:.10f} {r_frac:.10f}\n")

#print(f"Wrote {region_file}")
#print(f"center(frac) = {x_frac:.6f} {y_frac:.6f} {z_frac:.6f}   radius(frac) = {r_frac:.6f}")
#print(f"Nsel = {len(Qx)}  R_raw = {R_raw:.3f}  pad = {pad_lag}")
#print("halo_npart : ", halo_npart)
#print(f"center int = {x_frac*Ngridx_IC:.6f} {y_frac*Ngridx_IC:.6f} {z_frac*Ngridx_IC:.6f}   radius int = {r_frac*Ngridx_IC:.6f}")

print("")
print("")
print("rvir [cMpc/h] :", halo_rvir*boxlen_comov)
print("rlagrange [cMpc/h] :", r_frac*boxlen_comov)
print("")
print("")
print("rvir [kpc at z=6] :", 1000.0*halo_rvir*boxlen_comov/0.6774/(7.0))
print("rlagrange [kpc at z=6] :", 1000.0*r_frac*boxlen_comov/0.6774/(7.0))
