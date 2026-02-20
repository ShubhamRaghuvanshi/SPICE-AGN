import yt
import matplotlib.pyplot as plt
from PIL import Image, ImageDraw, ImageFont
import os
import numpy as np 
import yt_funcs

Omega_m = 0.308899998664856
Omega_b = 0.04488
h = 0.667
rho_crit0 = 1.27e11 #Msun/Mpc^3
L_HR = 100.0

simdir_1  = "/u/sraghu/sim/SHUB/b256/unibox/ref_0"
outnum_1  = 6
center_1  = [0.5,0.5,0.5]
boxlen_1  = 6.67
extent_1  = [L_HR/boxlen_1, L_HR/boxlen_1, L_HR/boxlen_1]
simtyp_1  = "SmoothSN_A(SPICE run)"
pltclr_1  = "b"
marker_1  = "*" 

simdir_arr = [simdir_1]
simtyp_arr = [simtyp_1]
center_arr = [center_1]
boxlen_arr = [boxlen_1]
extent_arr = [extent_1]
pltclr_arr = [pltclr_1]
marker_arr = [marker_1]

width_val = -1#L_HR/boxlen_1
 
outnum_arr    = [6]
pltclr_arr    = ['b']
linestyle_arr = [ '-']
propval  = "hmf"

fig1, ax1 = plt.subplots()
iout =0 
for outnum in outnum_arr :
  if(propval=="dmdens") :
    #outnum       = outnum_arr[iout]
    outnum_char  = str(outnum).zfill(5)
    outdir       = simdir_1 + "/output_" + outnum_char
    ds           = yt.load(outdir)
    Z            = ds.current_redshift
    ad           = ds.all_data()

    center       = ds.arr(center_1, "code_length")
    width_arr    = [width_val, width_val, width_val]
    R            = ds.arr(width_arr, "code_length")          # half-width ("radius") of the cube
    left         = center - R/2.0
    right        = center + R/2.0
    ds_box       = ds.box(left, right)
    Mzoom        = ds_box[("DM", "particle_mass")].sum().in_units('Msun').value 
    Vzoom        = (L_HR * L_HR * L_HR)
  
    rho_mean     = (Mzoom / Vzoom)
    rho_mean_th  = ds.quan( Omega_m * rho_crit0, "Msun/Mpccm**3")
    dm_dens      = ds_box[("deposit", "DM_density")].in_units('Msun/Mpccm**3')
    dm_overdens  = np.array(dm_dens/rho_mean_th )
    dm_overdens  = np.log10(dm_overdens[dm_overdens > 0.001])

    delta   = np.array(dm_dens/rho_mean_th - 1.0)
    mean_delta = np.sum(delta)/float(len(delta)) 
    var_delta = np.var(delta)
    drat = Omega_m*rho_mean/rho_mean_th.value

    print("mean densities :", Omega_m*rho_mean/1e11, rho_mean_th/1e11)
    print("Mbox : ", rho_mean_th.v * (boxlen_1/h) * (boxlen_1/h) * (boxlen_1/h))

    label_str = rf"z: {Z:.2f}, $\langle \rho_{{\rm DM}} \rangle/\rho_{{\rm crit}}$: {drat:.2f}, $\langle \delta_{{\rm DM}} \rangle$: {mean_delta:.2f}"
    #mean_rho_dm = (Omega_m - Omega_b)*rho_crit0
    ax1.hist(dm_overdens, bins=40, color=pltclr_arr[iout], label=label_str  ,edgecolor=pltclr_arr[iout], histtype='step', density=True)
  elif(propval=='hmf') :
    # write halo catalog 
    overwrite=True 
    yt_funcs.create_halo_catalog(simdir_1, outnum, 0, overwrite)
  elif(propval=='amrres') :
    #outnum       = outnum_arr[iout]
    outnum_char  = str(outnum).zfill(5)
    outdir       = simdir_1 + "/output_" + outnum_char
    ds           = yt.load(outdir)
    Z            = ds.current_redshift
    ad           = ds.all_data()

    v, c = ds.find_max(('deposit', 'star_cic'))
    c_arr = c
    print("c_arr :", c_arr)
    #center       = ds.arr(center_1, "code_length")
    center = c_arr 
    width_val = 0.02/boxlen_1
    print("rwidth_kpc :", (boxlen_1*width_val*1000.0/0.6774)/(1+Z) )


    width_arr    = [width_val, width_val, width_val]
    R            = ds.arr(width_arr, "code_length")          # half-width ("radius") of the cube
    left         = center - R/2.0
    right        = center + R/2.0
    ds_box       = ds.box(left, right)


    amr_dx       = ad[('gas', 'dx')].in_units('Mpccm/h').value 
    amr_level    = np.log2(boxlen_1/amr_dx) 
    amr_level    = np.rint(amr_level).astype(int)

    amr_dx_box    = ds_box[('gas', 'dx')].in_units('Mpccm/h').value 
  
    dens_box      = ds_box[('gas', 'density')]
    vol_box       = ds_box[('gas', 'cell_volume')]
    mcell_box         = (dens_box * vol_box).to('Msun').v
    mjeans_box        = ds_box[('gas','jeans_mass')].to('Msun').v

    for icell in range(0, len(mjeans_box)) :
      if(amr_dx_box[icell]< boxlen_1/2**17.2 and dens_box[icell] > 1e-23) : 
        print("mjeans/mcell :", mjeans_box[icell]/mcell_box[icell], np.log2(boxlen_1/amr_dx_box[icell]))

    amr_level_box    = np.log2(boxlen_1/amr_dx_box) 
    amr_level_box    = np.rint(amr_level_box).astype(int)


    label_str     = f"z: {Z:.2f}"

    # --- Count cells per refinement level ---
    unique_levels, counts = np.unique(amr_level, return_counts=True)
    unique_levels_box, counts_box = np.unique(amr_level_box, return_counts=True)

    # Print them
    for lvl, cnt in zip(unique_levels_box, counts_box):
        print(f"Level {lvl:2d} : {cnt:10d} cells")

    #ax1.hist(amr_level, bins=32, color=pltclr_arr[iout], label=label_str  ,edgecolor=pltclr_arr[iout], histtype='step')
    ax1.plot(unique_levels_box, counts_box, marker="o", color=pltclr_arr[iout], label=label_str, linestyle="-", ms=3, zorder=3)
    # --- Annotate each point with dx_phys in kpc ---
    a = 1.0 / (1.0 + Z)                       # scale factor
    h = float(ds.hubble_constant)              # little h (dimensionless)

    # cell size at level L: dx_comoving = boxlen / 2^L  [in Mpc_comoving/h]
    dx_cMpc_h = (boxlen_1 / (2.0**unique_levels))      # array in Mpc_comoving/h (numbers)
    # convert to physical kpc: ckpc/h * a/h = kpc
    dx_phys_kpc = dx_cMpc_h * 1.0e3 * (a / h)          # -> kpc (proper)

    for L, cnt, dxk in zip(unique_levels, counts, dx_phys_kpc):
        ax1.annotate(f"{dxk:.2f} kpc",
                    xy=(L, cnt),
                    xytext=(6, 6), textcoords="offset points",
                    fontsize=8, ha='left', va='bottom',
                    bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.6))
    #ax1.plot(unique_levels_box, counts_box, marker="^", color=pltclr_arr[iout], label="", linestyle="--", ms=3, zorder=3)
    #y0 = np.zeros_like(counts)
    #ax1.vlines(unique_levels, y0, counts,
    #          color=pltclr_arr[iout], linewidth=1.2, alpha=0.6, zorder=2)    

  iout = iout + 1

if(propval=="dmdens") :
  ax1.set_yscale('log')
  #ax1.set_xscale('log')
  ax1.set_title('DM only Fullbox')
  ax1.set_xlabel(r"$\log10(1+\delta_{\rm DM}) = \log10(\rho_{\rm DM}/\rho_{\rm crit})$")
  ax1.set_ylabel(r"$P(\log10(1+\delta_{\rm DM}))$")
  ax1.set_xlim(-1,7)
  ax1.legend(loc='best')
  imgname_1 = 'hist_dmsens.png'
  print("Saving:", imgname_1)
  fig1.savefig(imgname_1, dpi=200, bbox_inches="tight")
  fig1.clf()

elif(propval=='hmf') :
  fitfunc = 3
  #fitting_function (int) Which fitting function to use. 1 = Press-Schechter, 2 = Jenkins, 3 = Sheth-Tormen, 4 = Warren, 5 = Tinker Default : 4.
  yt_funcs.plot_hmf(simdir_1, pltclr_arr, outnum_arr, fitfunc, boxlen_1)
elif(propval=='amrres') :
  ax1.set_yscale('log')
  #ax1.set_xscale('log')
  ax1.set_title('Gas AMR levels [ilevel = log2(Lbox/dx)]')
  ax1.set_xlabel("ilevel")
  ax1.set_ylabel("counts")
  #ax1.set_xlim(-1,7)
  ax1.legend(loc='best')

  fig1.savefig('amr_res.png')
  #fig.savefig(outfile_pdf)
  print("Saved: amr_res.png")
  