import yt 
from yt.extensions.astro_analysis.halo_analysis import HaloCatalog
import yt.visualization.eps_writer as eps
import os
import numpy as np
import matplotlib.pyplot as plt
from halo_mass_function import *
from sfr_spectrum import * 
from PIL import Image, ImageDraw, ImageFont
from matplotlib.ticker import FuncFormatter
from scipy.spatial import cKDTree

# HOP(HaloFinderOverDensity
finder_kwargs_hop = {
    "threshold":160,  
    "ptype":"DM"
}

# FOF
finder_kwargs_fof = {
    "padding": 0.02,       
    "link": 0.2, 
    "ptype": "DM"
}

# ROCKSTAR
finder_kwargs_rockstar = {
    "num_particles_min": 100,  
    "halo_finder_path": "/path/to/rockstar_executable"
}

FinderMethods = ['hop', 'fof', 'rockstar' ]
FinderKwargs  = [finder_kwargs_hop, finder_kwargs_fof, finder_kwargs_rockstar ] 

def create_halo_catalog(simdir, outnum, MethId, overwrite) :
  outnum_char=str(outnum).zfill(5)
  outdir = simdir+'/output_' + outnum_char
  halocatalog_file = outdir + '/info_'+ outnum_char + '/info_' + outnum_char +'.0.h5'

  if( not os.path.exists(halocatalog_file) or overwrite) :
    print("creating halo catalog: ", halocatalog_file)
    print("finder method : ", FinderMethods[MethId])
    print("finder kwargs : ", FinderKwargs[MethId])
    ds = yt.load(outdir)
    halo_catalog = HaloCatalog(data_ds=ds, finder_method=FinderMethods[MethId], finder_kwargs=FinderKwargs[MethId], output_dir = outdir)
    halo_catalog.create()
  else :
    print("halo catalog already exists :", halocatalog_file)

log_M_min = 7
log_M_max = 14
delta_log_M = 0.1
hmf_fitfunc_label_arr = ["PS", "JN", "ST", "WR", "TK" ]
cbar = ['m', 'r', 'g', 'b', 'y', 'm', 'c', 'k']


# compare hmf from with same outnum at same redshift
def compare_hmf(simdir_arr, simname_arr, cbar_arr, mbar_arr, outnum_arr, fitfunc_arr, Boxlength_comov): 

  plt.xscale('log')
  plt.yscale('log')
  plt.xlabel("Halo Mass [Msun]")
  plt.ylabel("Number of Halos")
  iout=0
  for simdir in simdir_arr :
    outnum           = outnum_arr[iout]
    outnum_char      = str(outnum).zfill(5)
    outdir           = simdir+'/output_' + outnum_char
    halocatalog_file = outdir + '/info_'+ outnum_char + '/info_' + outnum_char +'.0.h5'
    print("Analyzing : ", halocatalog_file)

    ds               = yt.load(halocatalog_file)
    Z                = ds.current_redshift
    ad               = ds.all_data()  
    NumPart          = np.array( ad['halos', 'particle_number'] )  
    HMass            = np.array( ad["halos", "particle_mass"].in_units("Msun"))
    halo_masses      = []

    
    print("z :", f"{Z:.2f}")
    print("Mhalo_0 :", HMass[0]/1e12 )
    print("Numpart :", NumPart[0])
    if (iout == 0 ):
      title = "DM Halo Mass function at" + f' z={Z:.2f}'
      plt.title(title)
      Zint  = round(Z)
      Zint2 = round(Z)

    if(Zint != round(Z)) :
      print("output reshift dont match :", outdir)
      raise ValueError("Compare Halo Mass Functions onlt at the same redshift")

    for ihalo in range (0,len(NumPart)):
      if( NumPart[ihalo] > 10  ):
        halo_masses.append(HMass[ihalo])   

    Bins        = 10**np.arange(log_M_min, log_M_max + delta_log_M, delta_log_M)
    hist, edges = np.histogram(halo_masses, bins=Bins)
    bin_centers = (edges[1:] + edges[:-1]) / 2
    bin_widths  = edges[1:] - edges[:-1]
    dn_dlogm    = hist
    err         = np.sqrt(hist)  
    labl        = simname_arr[iout]
    print("label :", labl)
    plt.errorbar(bin_centers, dn_dlogm, yerr=err, fmt=mbar_arr[iout], color=cbar_arr[iout],capsize=3, label=labl)
    iout = iout + 1
  
  plotanal = True 
  if(plotanal == True) :
    # now plot the analytical hmf
    ifunc=0
    for fit_func in fitfunc_arr : 
      hmf    = HaloMassFcn(halos_ds=ds, log_mass_min=math.log10(np.min(halo_masses)), log_mass_max=math.log10(np.max(halo_masses)), 
                        fitting_function=fit_func, num_sigma_bins= len(bin_widths) +1)
      n_anal = hmf.dndM_dM_analytic[:-1]/0.6774
      n_anal = n_anal*Boxlength_comov*Boxlength_comov*Boxlength_comov*np.log10(bin_widths)
      m_anal = hmf.masses_analytic[:-1]
      anal_label = label_arr[fit_func -1]
      plt.plot(m_anal, n_anal, color=cbar[iout-1],label=anal_label)
      ifunc = ifunc + 1
 
  plt.legend()
  plt.xlim(1e8, 1e12)
  plt.tight_layout()  
  imgname='hmf_cmp_z' + str(Zint2) + '.png' 
  print('Saving : ', imgname)
  plt.savefig(imgname)
  plt.clf()

# compare hmf from with same outnum at same redshift
def plot_hmf(simdir_arr, simtyp_arr, outnum_arr, cbar_arr, fitfunc, Boxlength_comov): 

  isim=0
  for simdir in simdir_arr :
    outnum           = outnum_arr[isim]
    outnum_char      = str(outnum).zfill(5)
    outdir           = simdir+'/output_' + outnum_char
    halocatalog_file = outdir + '/info_'+ outnum_char + '/info_' + outnum_char +'.0.h5'
    print("Analyzing : ", halocatalog_file)

    ds               = yt.load(halocatalog_file)
    Z                = ds.current_redshift
    ad               = ds.all_data()  
    NumPart          = np.array( ad['halos', 'particle_number'] )  
    HMass            = np.array( ad["halos", "particle_mass"].in_units("Msun"))
    halo_masses      = []

    
    print("z :", f"{Z:.2f}")
    print("Mhalo_0 :", HMass[0]/1e12 )
    print("Numpart :", NumPart[0])
   # print("mDM :", (0.315)*HMass[0]/NumPart[0]/1e9, HMass[1]/NumPart[1]/1e9, HMass[2]/NumPart[2]/1e9)
    if (isim == 0 ):
      title = "DM Halo Mass Function [fitting func : " + hmf_fitfunc_label_arr[fitfunc -1] + " ]"
      plt.title(title)
      Zint  = round(Z)
      Zint2 = round(Z)

   
    # build the log10 mass array (as you already do)
    halo_masses = []
    for ihalo in range(len(NumPart)):
      if NumPart[ihalo] > 10:
        halo_masses.append(math.log10(HMass[ihalo]))
    halo_masses = np.asarray(halo_masses)

    print("Num halos in catalog:", len(NumPart))
    print("Num halos with >10 particles:", np.sum(NumPart > 10))
    print("Mass range (Msun):", HMass.min(), HMass.max())


    if halo_masses.size:
      # same number of bins as before
      bins = np.linspace(halo_masses.min(), halo_masses.max(), 41)  # 40 bins
      counts, edges = np.histogram(halo_masses, bins=bins)
      centers = 0.5 * (edges[1:] + edges[:-1])

      # Poisson errors on counts
      y = counts
      yerr = np.sqrt(counts)

      # (optional) skip empty bins so you don't plot zero-height points
      m = counts > 0
      centers, y, yerr = centers[m], y[m], yerr[m]

      # convert x from log10(M) to M for plotting on a log-mass axis
      x_mass = 10**centers

      # points with symmetric Poisson error bars (x now in mass units)
      #labl = f"Z:{Z:.2f}"
      labl = simtyp_arr[isim]
      plt.errorbar(x_mass, y, yerr=yerr, fmt='o', ms=3, lw=1,
                    color=cbar_arr[isim], ecolor=cbar_arr[isim], label=labl)
        # (optional) if you want to connect the points:
        # plt.plot(x_mass, y, '-', lw=1, color=cbar_arr[isim])

      # --- vertical line at 10 * M_DM ---
      # NOTE: x_mass is in Msun, so make sure M_DM is also in Msun (if it's Msun/h, divide by h).
      if(isim==0) :
        M_DM = 0.315*HMass[0]/NumPart[0]
        print("M_DM :", M_DM/1e9)
        x_mark = 100.0 * M_DM
        ax = plt.gca()
        ax.axvline(x_mark, ls='--', lw=2, color='grey', alpha=0.7)
        ax.annotate(r'$100\,m_{\rm DM}$', xy=(x_mark, 4.0*ax.get_ylim()[1]),
                  xytext=(4, -4), textcoords='offset points',
                  rotation=90, va='top', ha='left', color='grey')

        x_mark = 1000.0 * M_DM
        ax.axvline(x_mark, ls='--', lw=2, color='grey', alpha=0.99)
        ax.annotate(r'$1000\,m_{\rm DM}$', xy=(x_mark, 4.0*ax.get_ylim()[1]),
                  xytext=(4, -4), textcoords='offset points',
                  rotation=90, va='top', ha='left', color='grey')


      # bin widths in dex (Δlog10 M)
      bin_widths_log10 = edges[1:] - edges[:-1]

      plotanal = True
      if (plotanal == True):
        # analytical HMF on the same mass range/bins
        hmf = HaloMassFcn(
            halos_ds=ds,
            log_mass_min=np.min(np.log10(x_mass/2.0)),                 # <-- use your log10(M) left edge
            log_mass_max=np.max(np.log10(x_mass)),                # <-- use your log10(M) right edge
            fitting_function=fitfunc,
            num_sigma_bins=len(x_mass)              # <-- edges count, not centers
        )

        # (dn/dM)*dM per bin, units: 1/Mpc^3 (already integrated over dM)
        n_anal = hmf.dndM_dM_analytic[:-1]

        # get h; fall back to Planck-ish if not present
        h = 0.6774
        # convert to expected plot units:
        # counts per bin in the box = number density * volume
        vol_Mpc3 = (Boxlength_comov / h)**3         # if Boxlength_comov is in Mpc/h
        n_anal_counts = n_anal * vol_Mpc3
        print("n_anal_counts :", n_anal_counts)
        print("vol_Mpc3 =", vol_Mpc3)
        print("n_anal[0] =", n_anal[0])
        print("Δlog10M[0] =", (edges[1]-edges[0]))
        edges_M = 10**edges
        print("ΔM[0] =", edges_M[1]-edges_M[0])
        print("N_bin_expected[0] =", n_anal[0]*vol_Mpc3)

        # mass bin midpoints for x in Msun (masses_analytic is in h^-1 Msun)
        m_anal = 0.5 * (hmf.masses_analytic[:-1] + hmf.masses_analytic[1:]) / h
        plt.plot(m_anal, n_anal_counts, color=cbar_arr[isim], label="")
    else:
      print("No halos passed the selection.")
    isim = isim + 1  


  plt.xscale('log')
  plt.yscale('log')
  plt.xlabel(r'$M\ \mathrm{[M_\odot]}$')   # or your preferred units/label
  plt.ylabel('Count')
  plt.legend()
 # plt.xlim(1e11, 1e13)
 # plt.ylim(0.1, 50)
  plt.tight_layout()  
  imgname='hmf_cmp_z' + str(Zint2) + '.png' 
  print('Saving : ', imgname)
  plt.savefig(imgname)
  plt.clf()


#star functions 
def write_starfile(simdir, outnum, c_arr, width_val, fullbox, overwrite):
  outnum_char      = str(outnum).zfill(5)
  outdir           = simdir+'/output_' + outnum_char
  sfr_file         = outdir + "/SFR_" +  outnum_char + ".out" 
  print("sfr_outdir : ", outdir)

  if( not os.path.exists(sfr_file) or overwrite ) :     
    ds         = yt.load(outdir)
    if fullbox:
      print("using fullbox")
      left  = ds.domain_left_edge
      right = ds.domain_right_edge
      ds_box = ds.box(left, right)
    else:
      extent = [width_val, width_val, width_val]
      Rbox   = ds.arr(extent, "code_length")  # extent is FULL width
      center = ds.arr(c_arr, "code_length")
      left   = center - Rbox/2.0
      right  = center + Rbox/2.0
      ds_box = ds.box(left, right)

    mass       = ds_box[("star", "particle_mass")].in_units('Msun')
    age        = ds_box[("star", "age")].in_units('Myr')
    ct         = ds_box[("star", "particle_birth_time")]
    Z          = ds.current_redshift
    L_code     = ds.domain_width.to("Mpccm/h")
    
    boxlen_comov = L_code[0]
    print("Box length, Lzoom  [Mpccm/h]:",boxlen_comov, width_val*boxlen_comov)

    threshold  = ds.quan(0.0, "Myr")
    mass_star  = mass[age > threshold]
    ct_star    = ct[age > threshold]

    print("len mass star :", len(mass_star), len(mass))

    sfr  = StarFormationRate(ds, star_mass=mass_star, star_creation_time=ct_star, volume=ds_box.volume())
    print("Writing : ", sfr_file)
    #print( "Number of stars particles :", len(star_mass))
    sfr.write_out(name=sfr_file)
  else :
    print(sfr_file, "already exists")  

def to_xy(pairs):
    if not pairs:
        return np.array([]), np.array([])
    a = np.array(pairs, dtype=float)      # shape (N, 2) → [z, logSFR]
    a = a[np.argsort(a[:, 0])]            # sort by z
    return a[:, 0], a[:, 1]

# compare sfr from with same outnum at same redshift
def compare_sfr(simdir_arr, simname_arr ,outnum_arr, cbar_arr, linestyle_arr ,plotanal, plotdata ): 

  fig, ax1 = plt.subplots()
  ax1.set_yscale('log')
  #ax1.set_xscale('log')
  ax1.set_xlabel("Redshift(z)")
  ax1.set_ylabel(r"SFR  [$M_\odot$/yr/Mpcc] ")
  ax1.set_title("SFRD vs redshift")

  iout=0
  for simdir in simdir_arr :
    outnum      = outnum_arr[iout]  
    outnum_char = str(outnum).zfill(5)
    outdir      = simdir+'/output_' + outnum_char
    sfr_file    = outdir + '/SFR_'+ outnum_char +'.out' 
    Redshift    = []
    SFR         = []
    with open(sfr_file, 'r') as file:
      numline     = 0
      for line in file:
        numline = numline + 1
        if(numline > 1):
          tage = float(line.split()[0]) 
          z    = float(line.split()[2])
          if(z > 6.0): 
            sfrcc = float(line.split()[4])
            Redshift.append(z)
            SFR.append(sfrcc)        
    ax1.plot(Redshift, SFR, linestyle=linestyle_arr[iout], color=cbar_arr[iout], label=simname_arr[iout])  
    iout = iout + 1    

  if(plotdata) :
    sfrd_datafile="sfrd_digitized_all.csv"
    with open(sfrd_datafile, 'r') as file:  
      numline     = 0
      sfrd_smoothSN = []
      sfrd_burstySN = []
      sfrd_hyperSN  = []
      sfrd_obs      = []
      for line in file:
        numline = numline + 1
        if(numline > 1):
          zs_str, logSFR_str, data_source = line.rstrip("\n").split(",", 2)
          zs      = float(zs_str.strip())
          logSFR  = float(logSFR_str.strip())
          source  = data_source.strip()      # <-- remove trailing spaces/newlines
          #print("data_source :", source)
          if source == "smooth-sn (blue)":
              sfrd_smoothSN.append([zs, logSFR])
          elif source == "hyper-sn (teal)":
              sfrd_hyperSN.append([zs, logSFR])
          elif source == "bursty-sn (red)":
              sfrd_burstySN.append([zs, logSFR])
          elif source == "obs (grey points)":
              sfrd_obs.append([zs, logSFR])

    print("len :", len(sfrd_smoothSN), len(sfrd_hyperSN), len(sfrd_burstySN), len(sfrd_obs))
    zs_arr, sfrd_arr = to_xy(sfrd_smoothSN)
    ax1.plot(zs_arr, 10**sfrd_arr,  linestyle="-"  ,color='grey', label="SPICE production run(smoothSN)")
    zs_arr_obs             = [5.91 , 7.00 , 7.80 , 7.90 , 8.01 , 9.00 , 10.00, 11.00, 12.01, 12.01, 13.75, 17.00 ]
    logsfrd_obs          = [-1.65, -1.92, -1.07, -2.07, -2.19, -2.59, -2.85, -2.70, -3.37, -3.30, -3.77, -3.62 ]  

    zs_arr_obs  = np.array(zs_arr_obs)
    logsfrd_obs = np.array(logsfrd_obs)

    #ax1.errorbar(Z_data, sfrd, yerr=[sfrd_m, sfrd_p] ,fmt='o' ,color='grey', label="Donnan(2022)+Harikane(2023)")
    #ax1.plot(Z_data, sfrd_ab_bursty, marker='*', linestyle="--"  ,color='r', label="SPICE, bursty")
    ax1.scatter(zs_arr_obs, 10**logsfrd_obs, marker='*' ,color='grey', label="observation")

    if(plotanal) :
      Z_anal   =   [8,10,12,14,16]
      SFR_anal =   []
      SFR_aniket_bursty=[-1.96, -2.33, -2.77,  ]
      for i in range(0, len(Z_anal)) :
        z   = Z_anal[i]
        sfr = 0.015 * (1 + z)**2.7 / (1 + ((1+z)/2.9)**5.6)
        SFR_anal.append(sfr)
      ax1.plot(Z_anal, SFR_anal, linestyle="-"  ,color='k', label="Madau, Dickinson")

  imgname = "sfr_cmp.png"
  plt.legend()
  print("Saving : ", imgname)
  plt.tight_layout()
  plt.savefig(imgname)
  plt.clf()

def compare_res(simdir_arr, simname_arr, outnum_arr, cbar_arr, linestyle_arr):
  fig, ax1 = plt.subplots()
  isim = 0

  for simdir in simdir_arr:
    outnum      = outnum_arr[isim]
    outnum_char = str(outnum).zfill(5)
    outdir      = simdir + '/output_' + outnum_char

    ds = yt.load(outdir)
    ad = ds.all_data()

    L_code       = ds.domain_width.to("Mpccm/h")
    boxlen_comov = float(L_code[0].value)   # <- make float in Mpccm/h

    # jeans mass and cell mass
    dm_particle_mass = np.min(ad[('DM','particle_mass')].to('Msun').v)
    dm_masses = ad[('DM','particle_mass')].to('Msun').v
    print("unique DM masses:", np.unique(dm_masses))

    dgas_cell   = ad[('gas', 'density')]
    vgas_cell   = ad[('gas', 'cell_volume')]
    mgas_cell   = (dgas_cell * vgas_cell).to('Msun').v

    dDM_cell   = ad[('deposit', 'DM_density')]
    mDM_cell   = (dDM_cell * vgas_cell).to('Msun').v

    Omega_m = 0.308899998664856
    Omega_b = 0.04488
    mcell_tot = mDM_cell + (Omega_m/Omega_b) * mgas_cell


    amr_dx       = ad[('gas', 'dx')].in_units('Mpccm/h').value 
    amr_level    = np.log2(boxlen_comov/amr_dx) 
    amr_level    = np.rint(amr_level).astype(int)
    print("levels:", amr_level.min(), amr_level.max(), np.unique(amr_level)[:10])

    # hydrogen number density
    #dgas_cell   = ad[('gas', 'density')]
    mp = 1.673e-24  # gram
    nH = (dgas_cell / (mp * yt.units.g)).to('1/cm**3').v


    mask = (amr_level < 15)
    mass_ratio = mcell_tot[mask]/(8.0*dm_particle_mass)
    nH_plot = nH[mask]

    bins = np.logspace(
        np.log10(mass_ratio.min()),
        np.log10(mass_ratio.max()),
        50
    )

    ax1.scatter(nH_plot, mass_ratio, c=amr_level[mask], cmap='viridis')
    ax1.set_xscale("log")
    ax1.set_yscale("log")
    ax1.set_xscale("log")
    imgname = "amr_res_" + simname_arr[isim] + ".png"
    ax1.set_xscale("log")
    plt.legend()
    print("Saving : ", imgname)
    plt.tight_layout()
    plt.savefig(imgname)
    plt.clf()
    isim = isim + 1

def compare_amr(simdir_arr, simname_arr, outnum_arr, cbar_arr, linestyle_arr):
  fig1, ax1 = plt.subplots()
  fig2, ax2 = plt.subplots()
  fig3, ax3 = plt.subplots()

  isim = 0
  nsim = len(simdir_arr)
  for simdir in simdir_arr:
    outnum      = outnum_arr[isim]
    outnum_char = str(outnum).zfill(5)
    outdir      = simdir + '/output_' + outnum_char

    ds = yt.load(outdir)
    ad = ds.all_data()

    L_code       = ds.domain_width.to("Mpccm/h")
    boxlen_comov = float(L_code[0].value)   # <- make float in Mpccm/h

    amr_dx       = ad[('gas', 'dx')].in_units('Mpccm/h').value 
    amr_level    = np.log2(boxlen_comov/amr_dx) 
    amr_level    = np.rint(amr_level).astype(int)

    # --- Count cells per refinement level ---
    unique_levels, counts = np.unique(amr_level, return_counts=True)
 
    avg_rhogas = []
    avg_rhoDM = []
    avg_nH  = []
    avg_gasnH = []
    # hydrogen number density
    dgas_cell   = ad[('gas', 'density')]
    dDM_cell   = ad[('deposit', 'DM_density')]

    for lvl in unique_levels:
      mask = (amr_level == lvl)
      avg_rhogas.append(np.mean(dgas_cell[mask]))
      avg_rhoDM.append(np.mean(dDM_cell[mask]))

    Omega_m = 0.308899998664856
    Omega_b = 0.04488
    mp = 1.673e-24 
    XH = 0.76
    dm_particle_mass = np.min(ad[('DM','particle_mass')].to('Msun').v)
 
    for il in range(0, len(avg_rhoDM)) :
      avg_nH.append(avg_rhoDM[il] + ( Omega_m/Omega_b )*avg_rhogas[il])
      avg_gasnH.append(XH*avg_rhogas[il]/mp)

    # Print them
    for lvl, cnt in zip(unique_levels, counts):
        print(f"Level {lvl:2d} : {cnt:10d} cells")

    ax1.plot(unique_levels, counts, marker="o", color=cbar_arr[isim], label=simname_arr[isim], linestyle="-", ms=3, zorder=3)
    ax2.plot(unique_levels, avg_nH, marker="o", color=cbar_arr[isim], label=simname_arr[isim], linestyle="-", ms=3, zorder=3)
    ax3.plot(unique_levels, avg_gasnH, marker="o", color=cbar_arr[isim], label=simname_arr[isim], linestyle="-", ms=3, zorder=3)
    ax3.axhline(y=10.0,linestyle='--',color='k',linewidth=1.5)
  
    if(isim == nsim-1):
      # --- Annotate each point with dx_phys in kpc ---
      Z  = ds.current_redshift
      a = 1.0 / (1.0 + Z)                       # scale factor
      h = float(ds.hubble_constant)              # little h (dimensionless)

      # cell size at level L: dx_comoving = boxlen / 2^L  [in Mpc_comoving/h]
      dx_cMpc_h = (boxlen_comov / (2.0**unique_levels))      # array in Mpc_comoving/h (numbers)
      # convert to physical kpc: ckpc/h * a/h = kpc
      dx_phys_kpc = dx_cMpc_h * 1.0e3 * (a / h)          # -> kpc (proper)

      for L, cnt, rhoeff, dxk in zip(unique_levels, counts, avg_nH, dx_phys_kpc):    
        ax1.annotate(
              f"{dxk:.2f} kpc",
              xy=(L, cnt),
              xytext=(6, 6),
              textcoords="offset points",
              fontsize=8,
              ha='left',
              va='bottom',
              bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.6)
          )
        ax2.annotate(
              f"{dxk:.2f} kpc",
              xy=(L, rhoeff),
              xytext=(6, 6),
              textcoords="offset points",
              fontsize=8,
              ha='left',
              va='bottom',
              bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.6)
          )
    isim = isim + 1
  ax1.set_yscale('log')
  #ax1.set_xscale('log')
  ax1.set_xlabel('Gas AMR levels [ilevel = log2(Lbox/dx_cell)]')
  #ax1.set_xlabel("ilevel")
  ax1.set_ylabel("counts")
  #ax1.set_xlim(-1,7)
  ax1.legend(loc='best')
  fig1.savefig('amr_count.png')
  #fig.savefig(outfile_pdf)
  print("Saved: amr_count.png")

  ax2.set_yscale('log')
  ax2.set_title(f"Average effective density per AMR level at Z = {Z:.2f}")
  ax2.set_xlabel('Gas AMR levels [ilevel = log2(Lbox/dx_cell)]')
  ax2.set_ylabel(r"$\langle \rho_{\rm eff} \rangle$ [g cm$^{-3}$]")
  ax2.legend(loc='best')
  fig2.savefig('amr_rhoeff.png')
  #fig.savefig(outfile_pdf)
  print("Saved: amr_rhoeff.png")

  ax3.set_yscale('log')
  ax3.set_title(f"Average nH per AMR level at Z = {Z:.2f}")
  ax3.set_xlabel('Gas AMR levels [ilevel = log2(Lbox/dx_cell)]')
  ax3.set_ylabel("nH [cm$^{-3}$]")
  ax3.legend(loc='best')
  fig3.savefig('amr_nH.png')
  #fig.savefig(outfile_pdf)
  print("Saved: amr_nH.png")


def plot_hist_dmddens(simdir_arr, outnum_arr, simtyp_arr, cbar_arr) :
  fig1, ax1 = plt.subplots()
  isim=0
  for simdir in simdir_arr :
    outnum      = outnum_arr[isim]  
    outnum_char = str(outnum).zfill(5)
    outdir      = simdir+'/output_' + outnum_char
    ds = yt.load(outdir)
    ad     = ds.all_data()
    dmdens = np.array ( ad[("deposit", "DM_density")].in_units('g/cm**3') )
    ax1.hist(dmdens, bins=40, color=cbar_arr[isim], label=simtyp_arr[isim]  ,edgecolor=cbar_arr[isim], histtype='step')
    isim = isim + 1

  ax1.set_yscale('log')
  ax1.set_xscale('log')
  ax1.set_xlabel("Redshift(z)")
  ax1.set_xlabel(r"DMdens [g/cm**3]")
  ax1.legend(loc='best')
  imgname_1='hist_dmsens.png'
  print("Saving : ", imgname_1)
  fig1.savefig(imgname_1)
  fig1.clf()


def plot_proj_dmdens(simdir, outnum, simtyp, radius, center, projdir, boxlen_comov): 
  outnum_char  = str(outnum).zfill(5)
  outdir       = simdir+'/output_' + outnum_char
  print("outdir :",outdir)
  imgname      = "temp_dmdens_proj_"+ simtyp + "_" + outnum_char  +".png"
  
  ds = yt.load(outdir)

  Z  = ds.current_redshift
  Zint2 = round(Z)
  imgname_z    =  "dmdens_proj_" + simtyp + "_z"+ str(Zint2) +  ".png"

  if(radius <= 0) :
    radius = boxlen_comov/2.0

  p  = yt.ProjectionPlot(ds, projdir, ("deposit", "DM_density"), weight_field=("deposit", "DM_density"),width=(2.0*radius,"Mpccm/h"))

  if(radius > 0) :
    if(projdir == 'z')   :
      p.set_center(    (center[0], center[1])  )
    elif(projdir == 'x') :
      p.set_center(    (center[1], center[2])  )
    elif(projdir == 'y') :
      p.set_center(    (center[2], center[0])  )
    else :
      raise ValueError("Invalid direction for projection")              
  p.set_cmap(("deposit", "DM_density"), "viridis")
  p.set_xlabel("cMpc/h")
  p.set_ylabel("cMpc/h")
  p.set_zlim( ("deposit", "DM_density"), 1e-26, 1e-21 )
  p.annotate_title(simtyp)
  p.save(imgname)
  
  image = Image.open(imgname)
  draw = ImageDraw.Draw(image)
  myfont = ImageFont.truetype('/usr/share/fonts/truetype/DejaVuSansMono-Bold.ttf', 42)
  text = f"Z:{Z:.2f}"
  text_color = (255, 0, 0)
  draw.text((730, 36), text, font=myfont, fill=(0, 0, 0))
  print("Saving plot: ", imgname_z)
  image.save(imgname_z)

def plot_proj_gasdens(simdir, outnum, simtyp, radius, center, projdir, boxlen_comov): 
  outnum_char  = str(outnum).zfill(5)
  outdir       = simdir+'/output_' + outnum_char
  print("outdir :",outdir)
  imgname      = "temp_gasdens_proj_"+ simtyp + "_" + outnum_char  +".png"
  
  ds = yt.load(outdir)
  Z  = ds.current_redshift
  Zint2 = round(Z)
  imgname_z    =  "gasdens_proj_" + simtyp + "_z"+ str(Zint2) +  ".png"

  if(radius <= 0) :
    radius = boxlen_comov/2.0

  p  = yt.ProjectionPlot(ds, projdir, ("gas", "density"), weight_field=("gas", "density"),width=(2.0*radius,"Mpccm/h"))

  if(radius > 0) :
    if(projdir == 'z')   :
      p.set_center(    (center[0], center[1])  )
    elif(projdir == 'x') :
      p.set_center(    (center[1], center[2])  )
    elif(projdir == 'y') :
      p.set_center(    (center[2], center[0])  )
    else :
      raise ValueError("Invalid direction for projection")              
  p.set_cmap(("gas", "density"), "viridis")
  p.set_xlabel("cMpc/h")
  p.set_ylabel("cMpc/h")
  p.set_zlim( ("gas", "density"), 1e-27, 1e-23 )
  p.annotate_title(simtyp)
  p.save(imgname)
  
  image = Image.open(imgname)
  draw = ImageDraw.Draw(image)
  myfont = ImageFont.truetype('/usr/share/fonts/truetype/DejaVuSansMono-Bold.ttf', 42)
  text = f"Z:{Z:.2f}"
  text_color = (255, 0, 0)
  draw.text((730, 36), text, font=myfont, fill=(0, 0, 0))
  print("Saving plot: ", imgname_z)
  image.save(imgname_z)

def plot_proj_gastemp(simdir, outnum, simtyp, radius, center, projdir, boxlen_comov): 
  outnum_char  = str(outnum).zfill(5)
  outdir       = simdir+'/output_' + outnum_char
  print("outdir :",outdir)
  imgname      = "temp_gastemp_proj_"+ simtyp + "_" + outnum_char  +".png"
  
  ds = yt.load(outdir)
  Z  = ds.current_redshift
  Zint2 = round(Z)
  imgname_z    =  "gastemp_proj_" + simtyp + "_z"+ str(Zint2) +  ".png"

  if(radius <= 0) :
    radius = boxlen_comov/2.0

  p  = yt.ProjectionPlot(ds, projdir, ("gas", "temperature"), weight_field=("gas", "density"),width=(2.0*radius,"Mpccm/h"))
  #p  = yt.ProjectionPlot(ds, projdir, ("gas", "velocity_magnitude"), weight_field=("gas", "density"),width=(2.0*radius,"Mpccm/h"))

  if(radius > 0) :
    if(projdir == 'z')   :
      p.set_center(    (center[0], center[1])  )
    elif(projdir == 'x') :
      p.set_center(    (center[1], center[2])  )
    elif(projdir == 'y') :
      p.set_center(    (center[2], center[0])  )
    else :
      raise ValueError("Invalid direction for projection")              
#  p.set_cmap(("gas", "velocity_magnitude"), "twilight_shifted")
  p.set_cmap(("gas", "temperature"), "twilight_shifted")
  p.set_xlabel("cMpc/h")
  p.set_ylabel("cMpc/h")
#  p.set_zlim( ("gas", "velocity_magnitude"), 1e6, 1e7 )
  p.set_zlim( ("gas", "temperature"), 1e0, 1e6 )
  p.annotate_title(simtyp)
  p.save(imgname)
  
  image = Image.open(imgname)
  draw = ImageDraw.Draw(image)
  myfont = ImageFont.truetype('/usr/share/fonts/truetype/DejaVuSansMono-Bold.ttf', 42)
  text = f"Z:{Z:.2f}"
  text_color = (255, 0, 0)
  draw.text((730, 36), text, font=myfont, fill=(0, 0, 0))

  print("Saving plot: ", imgname_z)
  image.save(imgname_z)

def tn_phaseplot(simdir, outnum, simtyp, boxlen_comov, overwrite): 
  outnum_char  = str(outnum).zfill(5)
  outdir       = simdir+'/output_' + outnum_char
  print("outdir :",outdir)
  imgname      = "temp_tn_"+ simtyp + "_" + outnum_char  +".png"
  
  ds = yt.load(outdir)
  ad = ds.all_data()
  Z  = ds.current_redshift
  Zint2 = round(Z,2)
  imgname_z    =  "tn_" + simtyp + "_z"+ str(Zint2) +  ".png"

  #family = np.array(ad['all', 'particle_family'])
  #alldens = np.array(ad['ramses', 'Density'])
  #alltemp = np.array(ad['ramses', 'Pressure']) 

  #print("len : ", len(family),  len(alldens), len(alltemp))
  #for ip in range(0, len(family)) :
  #  if(family[ip] == 2.0) :
  #    print("family :", family[ip], alldens[ip], alltemp[ip])

  if( not os.path.exists(imgname_z) or overwrite) :
    plot = yt.PhasePlot(ad, (('gas', 'number_density')), ("gas", "temperature"), ("gas", "mass"))
    plot.set_unit(("gas", "mass"), "Msun")
    plot.set_unit(("gas", "number_density"), "1/cm**3")
    plot.set_ylim( 1, 1e9 )
    plot.set_xlim( 1e-7,1e4 )
    plot.set_zlim( ("gas", "mass"), 1e0, 1e10 )
    plot.annotate_title(simtyp)
    plot.save(imgname)
  
    image = Image.open(imgname)
    draw = ImageDraw.Draw(image)
    myfont = ImageFont.truetype('/usr/share/fonts/truetype/DejaVuSansMono-Bold.ttf', 42)
    text = f"Z:{Z:.2f}"
    text_color = (255, 0, 0)
    draw.text((730, 36), text, font=myfont, fill=(0, 0, 0))
    print("Saving plot: ", imgname_z)
    image.save(imgname_z)
  else :
    print(imgname_z, "already exists")  

def location_phasepoints(simdir, outnum, overwrite, T1, T2, n1, n2):

  file_name='tn_xyz.txt'
  if(overwrite or not os.path.exists(file_name)) :
    outnum_char  = str(outnum).zfill(5)
    outdir       = simdir+'/output_' + outnum_char  
    ds           = yt.load(outdir)

    my_sphere    =   ds.sphere([0.17, 0.13, 0.44], (0.4, "Mpccm/h"))
    #my_sphere    = ds.all_data()
    Z            = ds.current_redshift

    temp = my_sphere["gas", "temperature"]
    nden = my_sphere["gas", "number_density"].in_units("1/cm**3")
    posx = np.array(my_sphere["gas", "x"].in_units("Mpccm/h"))
    posy = np.array(my_sphere["gas", "y"].in_units("Mpccm/h"))
    posz = np.array(my_sphere["gas", "z"].in_units("Mpccm/h"))
    mass = np.array(my_sphere['gas', 'cell_mass'].in_units("Msun")  )
    print("proceesing gas cells....")

    iline=0
    gasmass=  0.0
    with open(file_name, 'w') as f:
      for itemp in range(0, len(temp)) :
        if(temp[itemp] > T1 and temp[itemp] < T2 and nden[itemp] > n1 and nden[itemp] < n2 ) :
          print("temp, dens :", temp[itemp], nden[itemp], math.log10(temp[itemp]), math.log10(nden[itemp]) )
          x = posx[itemp]
          y = posy[itemp]
          z = posz[itemp]
          gasmass = gasmass + mass[itemp]
          f.write(f"{x}\t{y}\t{z}\n")
          if(iline < 5) :
            print(x,y,z)
            iline = iline + 1 
    print("done... total gas mass : ", gasmass)   
  else :
    print(file_name, "exists, do you want to overwrite ? ")   



def plot_proj_DMdens(simdir, outnum, simtyp): 
  outnum_char  = str(outnum).zfill(5)
  outdir       =  simdir + "/output_" + outnum_char
  imgname      =  simdir + "/temp_dmdens_proj_"+ simtyp + "_" + outnum_char  +".png"
  imgname_z    =  simdir + "/dmdens_proj_" + simtyp + "_"+outnum_char +  ".png"

  ds = yt.load(outdir)
  Z  = ds.current_redshift
  p = yt.ProjectionPlot(ds, "z", ('deposit', 'DM_density'), weight_field=("deposit", "DM_density") ,width=(10.0,"Mpccm/h")) 
  p.set_cmap(('deposit', 'DM_density'), "cividis")
  p.save(imgname)

  image = Image.open(imgname)
  draw = ImageDraw.Draw(image)
  myfont = ImageFont.truetype('/usr/share/fonts/truetype/DejaVuSansMono-Bold.ttf', 42)
  text = f"Z:{Z:.2f}"
  text_color = (255, 0, 0)
  draw.text((730, 36), text, font=myfont, fill=(0, 0, 0))
  print("Saving plot: ", imgname_z)
  image.save(imgname_z)

def plot_proj_halo_gasdens(simdir, outnum, simtyp, numplot): 
  outnum_char      = str(outnum).zfill(5)
  outdir           =  simdir + "/output_" + outnum_char
  halocatalog_file = outdir + '/info_'+ outnum_char + '/info_' + outnum_char +'.0.h5'

  ds        = yt.load(outdir)
  halo_ds   = yt.load(halocatalog_file)

  ad        = ds.all_data()
  halo_ad   = halo_ds.all_data()

  Z         = ds.current_redshift

  #NumPart   = np.array( halo_ad['halos', 'particle_number'] )
  halo_mass = np.array( halo_ad['halos', 'particle_mass'].in_units("Msun") )
  Rvir      = np.array( halo_ad['halos', 'virial_radius'].in_units("Mpccm/h"))
  halo_posx = np.array( halo_ad['halos', 'particle_position_x'])
  halo_posy = np.array( halo_ad['halos', 'particle_position_y'] )
  halo_posz = np.array( halo_ad['halos', 'particle_position_z'])

  sink_posx = np.array( ad['sink', 'particle_position_x'].in_units("Mpccm/h"))
  sink_posy = np.array( ad['sink', 'particle_position_y'].in_units("Mpccm/h") )
  sink_posz = np.array( ad['sink', 'particle_position_z'].in_units("Mpccm/h"))


  for isink in range(0, len(sink_posx)) :
    sep = 0 
    for idim in range(0,3):
      sep = sep + (halo_posx[0] - sink_posx[isink]) * (halo_posx[0] - sink_posx[isink])
    sep = math.sqrt(sep) 
    if(sep < 5.0*Rvir[0]): 
      print("********", isink+1, sep*1000.0, Rvir[0]*1000.0)


  for iplot in range(0, numplot) :

    halocenter     = [halo_posx[iplot], halo_posy[iplot], halo_posz[iplot]]
    #Rvir[iplot] = 0.001
    #left_edge  = [halo_posx[iplot] - Rvir[iplot], halo_posy[iplot] - Rvir[iplot], halo_posz[iplot] - Rvir[iplot] ]
    #right_edge = [halo_posx[iplot] + Rvir[iplot], halo_posy[iplot] + Rvir[iplot], halo_posz[iplot] + Rvir[iplot] ]
    #left_edge = [halo_posx[iplot] - 0.1, halo_posy[iplot] - 0.1, halo_posz[iplot] -0.1]
    #right_edge = [halo_posx[iplot] + 0.1, halo_posy[iplot] + 0.1, halo_posz[iplot] +0.1]
    #center = ds.domain_center.copy()
    #cds        = ds.region(halocenter, left_edge, right_edge)
    #print(center)

    print(halocenter)
    print(Rvir[iplot])
    print(len(sink_posx))
    #print(Rvir[iplot])
    #print(halo_ad['halos', 'particle_position_x'][iplot])
    #print(halo_ad['halos', 'virial_radius'][iplot])
    #print(left_edge)
    #print(right_edge)
    #cds = ds.region(ds.domain_center, ds.domain_left_edge, ds.domain_right_edge)

    L = [1, 1, 0]  # vector normal to cutting plane
    north_vector = [-1, 1, 0]

    print("******* plotting halo of mass", halo_mass[iplot])

    imgname      =  simdir + "/temp_halo_gasdens_proj_"+ simtyp + "_" + outnum_char+ "_"+ str(iplot) +".png"
    imgname_z    =  simdir + "/halo_gasdens_proj_" + simtyp + "_"+outnum_char+ "_" + str(iplot) +".png"

    #p  = yt.ProjectionPlot(ds, L, ('gas', 'density'), weight_field=('gas', 'density'), width=(5.0*Rvir[iplot],"kpccm/h"),north_vector=north_vector)
    p  = yt.ProjectionPlot(ds, "z", ('gas', 'density'), weight_field=('gas', 'density'), center=halocenter, width=(10.0*Rvir[iplot], "Mpccm/h"))
    #p  = yt.ParticlePlot(ds, ("sink", "particle_position_x"), ("sink", "particle_position_y"), ("sink", "particle_mass"), center=halocenter, width=(5.0*Rvir[iplot], "Mpccm/h"))
    p.annotate_particles(width=(10.0*Rvir[iplot], "Mpccm/h"),ptype='sink', p_size=5.0)
    #p.set_unit(("sink", "particle_mass"), "Msun")
    #eps_fig = eps.single_plot(p)
    #eps_fig.circle(radius=0.2, loc=(0.5, 0.5))
    #eps_fig.save_fig("zoom-circle", format="eps")

    #p.set_width((10, "Mpccm/h"))
    #p.set_xlim(halo_posz[iplot] - Rvir[iplot], halo_posz[iplot] + Rvir[iplot])
    #p.set_ylim(halo_posz[iplot] - Rvir[iplot], halo_posz[iplot] + Rvir[iplot])
    #p.set_zlim(ad, halo_posz[iplot] - Rvir[iplot], halo_posz[iplot] + Rvir[iplot])

    #p.set_center((halo_posx[iplot], halo_posy[iplot]))
    p.save(imgname)

    image = Image.open(imgname)
    draw = ImageDraw.Draw(image)
    myfont = ImageFont.truetype('/usr/share/fonts/truetype/DejaVuSansMono-Bold.ttf', 42)
    text = f"Z:{Z:.2f}"
    text_color = (255, 0, 0)
    draw.text((730, 36), text, font=myfont, fill=(0, 0, 0))
    #width, height = image.size
    #radius = (width /Rvir[iplot] ) * Rvir[iplot]
    #imgcenter = (int(halocenter[2] * width), int(halocenter[0] * height))
    #draw.ellipse([(imgcenter[0] - radius, imgcenter[1] - radius),
    #          (imgcenter[0] + radius, imgcenter[1] + radius)],
    #        outline="black", width=2)
    print("Saving plot: ", imgname_z)
    image.save(imgname_z)

import numpy as np
import yt

# Optional: if SciPy exists, KDTree is the fast way
try:
    from scipy.spatial import cKDTree
    _HAS_SCIPY = True
except Exception:
    _HAS_SCIPY = False


def _delta_pos_periodic(dx, boxlen_code):
    """
    Return minimum-image displacement for periodic box.
    dx and boxlen_code are in the same code-length units.
    """
    return dx - boxlen_code * np.rint(dx / boxlen_code)


def write_halos(simdir, outnum, boxlen_comov):
  outnum_char      = str(outnum).zfill(5)
  outdir           = simdir+'/output_' + outnum_char
  halocatalog_file = outdir + '/info_'+ outnum_char + '/info_' + outnum_char +'.0.h5'
  halo_ds          = yt.load(halocatalog_file)
  halo_ad          = halo_ds.all_data()
  Zcurr            = halo_ds.current_redshift
  L_code           = halo_ds.domain_width.to("Mpccm/h")

  boxlen_comov = round(boxlen_comov, 3)
  L_code[0]    = round(float(L_code[0]), 4)
  L_code[1]    = round(float(L_code[1]), 4)
  L_code[2]    = round(float(L_code[2]), 4)
  
  print("Box length [Mpccm/h]:", L_code, boxlen_comov)
  #if(L_code[0] != boxlen_comov or L_code[1] != boxlen_comov or L_code[2] != boxlen_comov) :
  #  raise ValueError("Incorrect boxlength in write_halos", boxlen_comov,L_code )

  halo_mass = np.array( halo_ad['halos', 'particle_mass'].in_units("Msun") )
  Rvir      = np.array( halo_ad['halos', 'virial_radius'].in_units("kpccm/h"))
  halo_posx = np.array( halo_ad['halos', 'particle_position_x'])
  halo_posy = np.array( halo_ad['halos', 'particle_position_y'] )
  halo_posz = np.array( halo_ad['halos', 'particle_position_z'])
  NumPart   = np.array( halo_ad['halos', 'particle_number'] )
  nhalo     = len(halo_mass)

  outfile_h   = outdir + '/info_'+ outnum_char + "/halodata.csv"
  outfile_g = outdir + '/info_'+ outnum_char + '/halogasinfo_' + outnum_char + '.txt'

  if  (os.path.exists(outfile_h) ) :
    command = "rm " + outfile_h
    print("deleting existing ", outfile_h)
    os.system(command)
  if  (os.path.exists(outfile_g) ) :
    command = "rm " + outfile_g
    print("deleting existing ", outfile_g)
    os.system(command)


  with open(outfile_h, 'w') as file:
    for ihalo in range(0, nhalo):
      line = f"{ihalo+1}\t{halo_posx[ihalo]:.6f}\t{halo_posy[ihalo]:.6f}\t{halo_posz[ihalo]:.6f}\t{Rvir[ihalo]:.6f}\t{halo_mass[ihalo]:.6f}\t{NumPart[ihalo]:.6f}\t{Zcurr:.4f}\n"
      file.write(line)

  write_haloprops = True  
  if(write_haloprops) :
    ds = yt.load(outdir)

    print("writing halo gas fractions")
    print("writing :", outfile_g)

    # Load gas once
    gad = ds.all_data()
    gas_m  = gad[("gas", "cell_mass")].to("Msun").d
    gas_x  = gad[("gas", "x")].to_value(ds.length_unit)
    gas_y  = gad[("gas", "y")].to_value(ds.length_unit)
    gas_z  = gad[("gas", "z")].to_value(ds.length_unit)

    # Your existing conversion (keep as-is, as requested)
    Rvir_code = (Rvir/1000.0)/boxlen_comov

    gas_pos = np.column_stack((gas_x, gas_y, gas_z))
    # Periodic boundaries: boxsize must be in SAME units as gas_pos and halo_pos
    # You are using boxlen_comov as the periodic box length (keep as-is).
    tree = cKDTree(gas_pos, boxsize=boxlen_comov)
  
    with open(outfile_g, "w") as f:
        #f.write("# id\tx\ty\tz\tRvir_kpccm_h\tMh_Msun\tNgas_0p5Rvir\tNgas_1p0Rvir\t"
        #        "Mgas_0p5Rvir_Msun\tMgas_1p0Rvir_Msun\tz\n")

        for ihalo in range(nhalo):
            print(f"\rProcessing halo: {ihalo+1}/{nhalo}", end="", flush=True)

            cx = float(halo_posx[ihalo])
            cy = float(halo_posy[ihalo])
            cz = float(halo_posz[ihalo])

            rvir = float(Rvir_code[ihalo])
            r05  = 0.5 * rvir
            r10  = 1.0 * rvir

            if tree is not None:
                # KDTree fast query
                idx05 = tree.query_ball_point([cx, cy, cz], r05)
                idx10 = tree.query_ball_point([cx, cy, cz], r10)

                Ngas05 = len(idx05)
                Ngas10 = len(idx10)

                Mgas05 = float(gas_m[idx05].sum()) if Ngas05 > 0 else 0.0
                Mgas10 = float(gas_m[idx10].sum()) if Ngas10 > 0 else 0.0

            else:
                # Slow fallback (your old method)
                dx = _delta_pos_periodic(gas_x - cx, boxlen_comov)
                dy = _delta_pos_periodic(gas_y - cy, boxlen_comov)
                dz = _delta_pos_periodic(gas_z - cz, boxlen_comov)
                r  = np.sqrt(dx*dx + dy*dy + dz*dz)

                m05_mask = (r <= r05)
                m10_mask = (r <= r10)

                Ngas05 = int(np.count_nonzero(m05_mask))
                Ngas10 = int(np.count_nonzero(m10_mask))

                Mgas05 = float(gas_m[m05_mask].sum()) if Ngas05 > 0 else 0.0
                Mgas10 = float(gas_m[m10_mask].sum()) if Ngas10 > 0 else 0.0

            line = (
                f"{ihalo+1}\t"
                f"{cx:.6f}\t{cy:.6f}\t{cz:.6f}\t"
                f"{float(Rvir[ihalo]):.6f}\t"
                f"{float(halo_mass[ihalo]):.6e}\t"
                f"{Ngas05}\t{Ngas10}\t"
                f"{Mgas05:.6e}\t{Mgas10:.6e}\t"
                f"{Zcurr:.4f}\n"
            )
            f.write(line)

    print("")
    print("Wrote:", outfile_g)



def volumefrac(simdir, outnum, lbox, t_frac):
  Mpc  = 3.0856e24
  h    = 0.67739997
  Mpch = Mpc/h  

  outnum_char = str(outnum).zfill(5)
  outdir      = simdir + "/output_" +  outnum_char
  if(os.path.exists(outdir)) :
    print("outdir : ", outnum_char)
    ds    = yt.load(outdir)
    ad    = ds.all_data()
    Z     = ds.current_redshift
    aexp  = 1.0/(1.0+Z)
    temp  = t_frac*np.array( ad[('gas', 'temperature')]  )
    cvol  = np.array( ad[('gas', 'cell_volume')] )
    ncell = len(temp)

    volfrac = 0.0 
    for icell in range(0, ncell) :
      celldx  = pow(cvol[icell], 1.0/3.0)/(Mpch*lbox*aexp)
      if( temp[icell] > 1e6 ) :
        volfrac = celldx*celldx*celldx + volfrac
      #print(- math.log2(celldx), volfrac)
    print(outnum_char, Z, volfrac)
    return[Z, volfrac]
  else :
    raise ValueError("My error: File" , outdir ,"does not exist")

def prop_diffsq(simdir_1, simdir_2, outnum):
  simdir_2 = simdir_1
  outnum_char      = str(outnum).zfill(5)
  outdir_1      = simdir_1+'/output_' + outnum_char
  outdir_2      = simdir_2+'/output_' + outnum_char

  ds_1 = yt.load(outdir_1)
  ds_2 = yt.load(outdir_2)
  ad_1 = ds_1.all_data()
  ad_2 = ds_2.all_data()

  Z_1 = ds_1.current_redshift
  Z_2 = ds_2.current_redshift

  prop_1 = np.array(ad_1['deposit', 'DM_density'])
  prop_2 = np.array(ad_2['deposit', 'DM_density'])

  print("Data lengths : ", len(prop_1), len(prop_2))
  if( len(prop_1) == len(prop_2) and Z_1 == Z_2) :
    sum_sq = 0  
    for i in range(0, len(prop_1)):
      diff = prop_1[i] - prop_2[i]
      sum_sq = sum_sq + diff*diff
    print("Sum of squares of difference : ", sum_sq)  
  else : 
    print( len(prop_1) ,len(prop_2), Z_1, Z_2 )
    raise ValueError("Apples and oranges")
