import yt
import yt_funcs
import matplotlib.pyplot as plt
from PIL import Image, ImageDraw, ImageFont
import os
import numpy as np 
import ion_fields
import gastracer
import matplotlib.colors as mcolors
from matplotlib.patches import Rectangle
from matplotlib.patches import Circle
from pylab import *
import math


# Build a FIRE/Illustris-like colormap
colors = [
    (0.0,  '#0a0a2a'),   # dark navy
    (0.25, '#213b7c'),   # blue
    (0.5,  '#5e4fa2'),   # indigo
    (0.7,  '#d94801'),   # orange
    (0.9,  '#fdae61'),   # light orange
    (1.0,  '#fefefe')    # white
]
cmap = mcolors.LinearSegmentedColormap.from_list("Hcol", colors)
norm = mcolors.LogNorm(vmin=1e20, vmax=1e23)


def rho_NFW(Mvir, Rvir_kpc, c, Npoints, rmin) :
  rmin = rmin * Rvir_kpc
  rs = Rvir_kpc/c
  f_c = np.log(1.0 + c) - c / (1.0 + c)
  #Mvir = 4.0*pi*rho_0*rs*rs*rs*f_c
  rho_0 = Mvir/( 4.0*np.pi*rs*rs*rs*f_c )
  r = np.linspace(rmin, Rvir_kpc, Npoints)
  x = r/rs
  rho_nfw = rho_0/(x*(1.0+x)*(1.0+x)) #Msun/kpc**3
  return [r/Rvir_kpc,rho_nfw]


def chi2_c(rho_fit, Mvir, Rvir_kpc, c, rmin ):
    """
    Simple chi^2 in log-space between sim profile and NFW for a given c.
    """
    Np = len(rho_fit)
    #if c <= 1.0 or c >= 30.0:  # sanity bounds for concentration
    #    return np.inf
    rho_nfw   = rho_NFW(Mvir, Rvir_kpc, c, Np, rmin)
    rho_model = rho_nfw[1]
    mask = (rho_model >0) & (rho_fit>0)
    rho_model = rho_model[mask]
    rho_fit   = rho_fit[mask] 
    # work in log-space to weight inner/outer radii comparably
    chi_sq = np.mean((np.log10(rho_fit) - np.log10(rho_model))**2)
    if(isinf(chi_sq) ) :
      raise ValueError("increase rmin")
    #print("chi^2 for c =:", c,  chi_sq)
    return math.sqrt(chi_sq)


Npart = 256
#simdir = "/ptmp/mpa/sraghu/smoothSN256L10/nosink"
simdir = "/ptmp/mpa/sraghu/smooth/sink_drag/star_sink"

hubble = 0.6774
c_arr = [0.5, 0.5, 0.5]
boxlen_comov = 6.774
projdir = "x"
width_val = 1.0#2.0/boxlen_comov

dsname = 'gas'
fieldname_arr = ["temperature"]
particleplot = False
phaseplot = False

#proval_arr = ["ndens", "temp", "dmdens"]

simtyp    = "star_S2"
outnum     = 11
outnum_char      = str(outnum).zfill(5)


outdir = simdir + "/output_" + outnum_char
ds     = yt.load(outdir)
ad     = ds.all_data()
Z      = ds.current_redshift
Z_int  = round(Z)

#v, c = ds.find_max(('deposit', 'DM_density'))
#c_arr = c.v


if(dsname == "gas"):
  ion_fields.register_number_density_and_mass_fields(ds)
#gastracer.register_gas_particle_aliases(ds)
print("Fields registered")


# --- define the cube (code units) ---
#v, cmax = ds.find_max(("gas", "density"))

#v, c = ds.find_max(("gas", "density"))
##p.set_center((c[0], c[1]))
#c_arr = c

#halonum =   int(sys.argv[1])
#halonum =   2
halonum = -1
  

if(halonum>0) : 
  halofile    = outdir + '/info_'+ outnum_char + "/halodata.csv"
  #read halos
  if  (not os.path.exists(halofile) ) : 
    print('outdir :', outdir)
    print('creating halodata.csv, this needs halo catalog')
    yt_funcs.create_halo_catalog(simdir, outnum, 0, 0)
    yt_funcs.write_halos(simdir, outnum, boxlen_comov, c_arr , [boxlen_comov, boxlen_comov, boxlen_comov])
  else :
    #yt_funcs.write_halos(simdir, outnum, boxlen_comov, c_arr , [boxlen_comov, boxlen_comov, boxlen_comov])
    print(halofile, "already exists")  

  if  (not os.path.exists(halofile) ) :
    raise ValueError( halofile, "not found")
  ihalo=1
  with open(halofile, 'r') as file:
    for line in file:
      ihalo = ihalo + 1
      if(ihalo > halonum) :
        halo_posx = float(line.split()[1])
        halo_posy = float(line.split()[2])
        halo_posz = float(line.split()[3])  
        halo_rvir = float(line.split()[4]) / 1000.0 / boxlen_comov  #code
        #halo_rvir = halo_rvir*0.6774
        halo_mass = float(line.split()[5])
        numpart   = float(line.split()[6])
        halo_z    = float(line.split()[7])
        break

  c_arr = [halo_posx, halo_posy, halo_posz]
  width_val = 2.0*halo_rvir
  #c_arr = [0.0, 0.076, 0.971]
  #width_val =8.0/boxlen_comov

  print("float(line.split()[4]) :", float(line.split()[4]))
  print("halo center       :", c_arr)
  print("halo rvir[cMpc/h] :", halo_rvir*boxlen_comov)
  print("rvir[kpc] :"        , halo_rvir * boxlen_comov * 1000.0 / 0.6774/(1.0+halo_z))
  print("Lbox[cMpc/h]      :", width_val*boxlen_comov)
  print("halo rvir[ckpc]   :", halo_rvir*boxlen_comov*1000.0*0.6774)
  print("Mhalo_12[Msun]    :", halo_mass/1e12)
  print("mDM_7             :", halo_mass/numpart/1e7)
  print("Numpart_halo      :", numpart)
  rmin  = 0.1
  print("rmin [kpc]      :", rmin*halo_rvir*boxlen_comov*1000.0*0.6774)

  plothaloprops=False
  if(plothaloprops) :
    #halo_rvir   = (1+halo_z)*(63.441 / 1000.0/boxlen_comov)/0.6774
    my_sphere   = ds.sphere(c_arr, (halo_rvir*boxlen_comov, "Mpccm/h"))
    M_arr       = my_sphere[("DM", "particle_mass")].in_units("Msun").v
    x_arr       = my_sphere[("DM", "particle_position_x")].v - c_arr[0]
    y_arr       = my_sphere[("DM", "particle_position_y")].v - c_arr[1]
    z_arr       = my_sphere[("DM", "particle_position_z")].v - c_arr[2]

    Lh  = 1.0
    Lhb = 1.0/2.0
    for ip in range(0, len(x_arr)) :
      if(x_arr[ip] > Lhb) :
        x_arr[ip] = x_arr[ip] - Lh
      if(x_arr[ip] < -Lhb) :
        x_arr[ip] = x_arr[ip] + Lh 
      if(y_arr[ip] > Lhb) :
        y_arr[ip] = y_arr[ip] - Lh
      if(y_arr[ip] < -Lhb) :
        y_arr[ip] = y_arr[ip] + Lh 
      if(z_arr[ip] > Lhb) :
        z_arr[ip] = z_arr[ip] - Lh
      if(z_arr[ip] < -Lhb) :
        z_arr[ip] = z_arr[ip] + Lh 


    r_arr       = np.sqrt(x_arr*x_arr + y_arr*y_arr + z_arr*z_arr)
    r_over_Rvir = r_arr/halo_rvir
    print( " min, max :", np.min(r_arr), np.max(r_arr)  )

    m_sphere = my_sphere[("DM", "particle_mass")].in_units("Msun").sum().v
    print("halo_mass (catalog) [Msun]:", halo_mass)
    print("enclosed mass in sphere [Msun]:", m_sphere)
    print("ratio (catalog / enclosed):", halo_mass / m_sphere)
    Mh = halo_mass
    halo_mass = m_sphere


    Nbins = 32
    mass_shell = np.zeros(Nbins)

    #edges = np.linspace(rmin, max(r_over_Rvir), Nbins+1)       # bin edges in r/Rvir
    #rmid_arr = 0.5 * (edges[:-1] + edges[1:])         # midpoints
  
    edges = np.logspace(np.log10(rmin), np.log10(max(r_over_Rvir)), Nbins + 1)
    rmid_arr = 0.5 * (edges[:-1] + edges[1:])
  
  
    for i in range(Nbins):
        in_shell = (r_over_Rvir >= edges[i]) & (r_over_Rvir < edges[i+1])
        mass_shell[i] = M_arr[in_shell].sum()
        #print("ibin :", i, mass_shell[i], r_over_Rvir[i], edges[i], edges[i+1])
    # ----- convert to physical densities -----
  
    # Rvir_kpc is physical virial radius
    Rvir_kpc    = halo_rvir * boxlen_comov * 1000.0 / 0.6774 # ckpc
    Rvir_kpc    = Rvir_kpc/(1.0+halo_z)  # ckpc -> kpc
    r_inner_kpc = edges[:-1] * Rvir_kpc
    r_outer_kpc = edges[1:]  * Rvir_kpc
    r_mid_kpc   = rmid_arr   * Rvir_kpc
    vol_shell   = (4.0/3.0) * np.pi * (r_outer_kpc**3 - r_inner_kpc**3)  # kpc^3
    rho_shell   = mass_shell / vol_shell                                  # Msun/kpc^3

      # crude grid search for best-fit c
    c_grid = np.linspace(1.0, 30.0, 4096)    
    chi_grid = np.array([chi2_c(rho_shell, halo_mass, Rvir_kpc, c, rmin ) for c in c_grid])
    c_best = c_grid[np.argmin(chi_grid)]

    print("chi_grid :", chi_grid)
    print(f"Best-fit NFW concentration: c_NFW ≈ {c_best:.2f}")
    c_nfw = c_best 
 
    rho_nfw = rho_NFW(halo_mass, Rvir_kpc, c_nfw, Nbins+1, rmin)
    #print("rho_nfw :", rho_nfw[1])
    figdr, axdr = plt.subplots()
    label = (
        f"\ncenter = ({c_arr[0]:.3f}, {c_arr[1]:.3f}, {c_arr[2]:.3f})"
        f"\nRvir [kpc] = {Rvir_kpc:.3f}"
        f"\nMvir [1e12Msun] = {Mh/1e12:.3f}"
    )
    axdr.plot(rmid_arr , rho_shell, drawstyle="steps-mid", label=label)
    label = (
        fr"NFW (c = {c_nfw:.1f})"
        f"\nrho_rms = {np.min(chi_grid):.3f}"
    )
    axdr.plot(rho_nfw[0], rho_nfw[1], "--", label=label)
  
    #c=1.6
    rho_nfw = rho_NFW(halo_mass, Rvir_kpc, 1.6, Nbins+1, rmin)
    label = (
        fr"NFW (c = {1.6:.1f})"
    )
    axdr.plot(rho_nfw[0], rho_nfw[1],"--", label=label, color='grey')


    axdr.set_xscale("log")
    axdr.set_yscale("log")
    axdr.set_xlabel(r"$r / R_{\rm vir}$")
    axdr.set_ylabel(r"$\rho_{\rm DM}\ [{\rm M_\odot\,kpc^{-3}}]$")
    axdr.set_xlim(rmin, 1.0)  # tweak as you like
    axdr.legend()
    figdr.tight_layout()

    imgname = "halo_rho_r_"+ str(Npart) + "_" + str(halonum)  + ".png"
    print("saving :", imgname)
    figdr.savefig(imgname, dpi=200)    
    #raise ValueError("feeef")

    #plot halo radial density profile here as function fo r/rvir


Ncell = len( ad[('gas', 'density')].v)

v, c = ds.find_max(("gas", "density"))
#p.set_center((c[0], c[1], c[2]))
c_arr = c.to_value()
print("c_arr :", c_arr)
center = ds.arr(c_arr, "code_length")

width_arr =  [width_val, width_val, width_val]
R      = ds.arr(width_arr, "code_length")          # half-width ("radius") of the cube
left   = center - R/2.0
right  = center + R/2.0
box    = ds.box(left, right)


projwidth_arr =  [width_val, width_val, width_val]
Rproj      = ds.arr(projwidth_arr, "code_length")          # half-width ("radius") of the cube
projleft   = center - Rproj/2.0
projright  = center + Rproj/2.0
projbox    = ds.box(projleft, projright)

#temp = projbox[('gas', 'temperature')]   # YTArray
#Ncell = int(temp.size)
#Ncell = 10002440

# 2) With numeric formatting
#tmin, tmax = box.quantities.extrema(("gas", "temperature"))
#tmin_K = tmin.to_value("K")   # strip units -> float
#tmax_K = tmax.to_value("K")
#print(f"T_min = {tmin_K:.3g} K, T_max = {tmax_K:.3g} K")



for fieldname in fieldname_arr :
  print("fieldname :", fieldname )
  ad = ds.all_data()

  if(fieldname == "metallicity") :
    Zmetal = ad[(dsname,fieldname)]
    print("Z min/mean/max:", Zmetal.min().to_value(''), Zmetal.mean().to_value(''), Zmetal.max().to_value(''))
  #if(fieldname=="X_HII") :
   # rho     = ad[('ramses','Density')].to('g/cm**3')
   # mask = (rho > 1e-26)
   # Zf      = ad[('ramses','Metallicity')].to_value()[mask]  # mass fraction
   # x_HII   = ad[('ramses','hydro_scalar_02')].to_value()[mask]
   # x_HeII  = ad[('ramses','hydro_scalar_03')].to_value()[mask]
   # x_HeIII = ad[('ramses','hydro_scalar_04')].to_value()[mask]
   # T       = ad[('gas','temperature')].to_value('K')[mask]
   # rho     = rho[mask]  
    # n = nHI + nHII + nHeI + nHeII + nHeIII + nZ
   # print("corr( log10 T , x_HII ) =", np.corrcoef(np.log10(T+1), x_HII)[0,1])
   # mw_xHII = np.sum(rho * x_HII) / np.sum(rho)
   # print("<x_HII>_mass =", mw_xHII)
   # print("x_HII ~", np.min(x_HII), np.max(x_HII))
   # print("x_HeII ~", np.min(x_HeII), np.max(x_HeII))
   # print("x_HeIII ~", np.min(x_HeIII), np.max(x_HeIII))
   # print("x_HeII + x_HeIII ~", np.min(x_HeII+x_HeIII), np.max(x_HeII+x_HeIII))

  if(phaseplot) :
    if(dsname != "gas") :
      raise ValueError("dsname must be gas for phaseplot")
    p = yt.PhasePlot(box, ('gas','number_density'), ("gas", "temperature"), (dsname, fieldname))
   # p.set_unit((dsname, fieldname), "Msun")
    p.set_unit(("gas", "number_density"), "1/cm**3")
    p.annotate_title(simtyp)
    p.set_ylim( 1e0, 1e9 )
    p.set_xlim( 1e-5,1e4 )
    p.set_zlim( (dsname, fieldname), 1e1, 1e5 )
    imgname = "tn_" + fieldname + "_" + simtyp + ".png" 
    p.save(imgname)
    # annotate redshift
    image  = Image.open(imgname)
    draw   = ImageDraw.Draw(image)
    myfont = ImageFont.truetype('/usr/share/fonts/truetype/DejaVuSansMono-Bold.ttf', 50)
    draw.text((730, 35), f"z:{Z:.2f}", font=myfont, fill=(0, 0, 0))
    print("Saving plot:", imgname)
    image.save(imgname)

  else : 
    if(particleplot) :
      p = yt.ParticleProjectionPlot(
        ds, projdir, (dsname, fieldname),
        center=center,
        width=(width_val*boxlen_comov, "Mpccm/h"),
        data_source=projbox, density=False
      )

      #p.annotate_particles(
      #  (width_val*boxlen_comov, "Mpccm/h"),  # thickness along LOS
      #  p_size=2.0,                          # increase this for bigger points
      #  col="white",
      #  marker="o",
      #  alpha=0.5,
      #  data_source=projbox,                 # same region you used above
      #)

      p.set_unit((dsname, fieldname), "Msun")
      #p.set_zlim((dsname, fieldname), 1e9, 1e10)
      p.set_cmap((dsname, fieldname), "viridis")
      p.set_background_color((dsname, fieldname), "black")
      #p.set_plot_args((dsname, fieldname), dict(marker='o', s=6, linewidths=0, alpha=0.8))
      p.set_log((dsname,fieldname), True)
      imgname = "pp_" + dsname[:2] + "_" + fieldname +"_" + simtyp +  "_" + str(Npart) + "_" + str(halonum) + ".png"
      #imgname = "pp_" + dsname[:2] + "_" + fieldname +"_" + simtyp  + ".png"
    else :
      if(fieldname == "DM_density") : 
        p = yt.ProjectionPlot(
            ds, projdir, (dsname, fieldname),
            weight_field=(dsname, 'DM_mass'),
            center=center,
            width=(width_val*boxlen_comov, "Mpccm/h"),
            #data_source=projbox
        )
      else :  
        p = yt.ProjectionPlot(
            ds, projdir, (dsname, fieldname),
            weight_field=(dsname, 'density'),
            center=center,
            width=(width_val*boxlen_comov, "Mpccm/h")
            #,data_source=projbox
        )
      
      nfields  = ["number_density","nH","nHI","nHII","nHe","nHeI","nHeII","nHeIII"]
      xfields  = ["X_H","X_HI","X_HII","Y_He","Y_HeI","Y_HeII","Y_HeIII","Z","mass_fraction_sum"]
      rhoflds  = ["DM_density","density","rhoH","rhoHI","rhoHII","rhoHe","rhoHeI","rhoHeII","rhoHeIII"]
      mfields  = ["mH","mHI","mHII","mHe","mHeI","mHeII","mHeIII"]

      if fieldname in nfields:
        p.set_cmap((dsname, fieldname), cmap)
        p.set_unit((dsname, fieldname), "1/cm**3")
        p.set_zlim((dsname, fieldname), 1e-4,1e0)
      if fieldname in xfields:
        p.set_cmap((dsname, fieldname), "cividis")
        #p.set_unit((dsname, fieldname), "1/cm**3")
        p.set_zlim((dsname, fieldname), 1e-5, 1)
      if fieldname in rhoflds:
        p.set_cmap((dsname, fieldname), "cividis")
        #p.set_unit((dsname, fieldname), "1/cm**3")
        #p.set_zlim((dsname, fieldname), 1e-30, 1e-23)
      if fieldname in mfields:
        p.set_cmap((dsname, fieldname), "ocean")
        #p.set_unit((dsname, fieldname), "1/cm**3")
        #p.set_zlim((dsname, fieldname), 1e-30, 1e-23)
      if ( fieldname == "temperature" ) :
        p.set_cmap((dsname, fieldname), "gist_heat") # gist_heat, berlin,  bwr, coolwarm
        #p.set_unit((dsname, fieldname), "1/cm**3")
        p.set_zlim((dsname, fieldname), 10.0, 2e6)
      imgname = "proj_" + fieldname + "_" + simtyp + ".png"
    
    #common params
    p.annotate_scale()
    p.annotate_title(simtyp)
    p.set_xlabel("cMpc/h")
    p.set_ylabel("cMpc/h")
    #p.set_buff_size((2048, 2048))
    #p.hide_colorbar()
    #p.annotate_grids(alpha=0.0025, edgecolors='black')

    p._setup_plots()
    ax  = p.plots[(dsname, fieldname)].axes
    fig = p.plots[(dsname, fieldname)].figure
    # store current limits
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    ax.autoscale(False)
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)

    #overplot stuff 
    plotcircle = False
    if plotcircle :
      Rcirc   = halo_rvir * boxlen_comov     # circle radius in cMpc/h
      xc      = 0.0                           # circle center x (your zoom region is centered)
      yc      = 0.0                           # circle center y
      # draw circle without changing limits
      circ = Circle((xc, yc), Rcirc, fill=False, linewidth=2, color='white')
      ax.add_patch(circ)
      # prevent autoscaling caused by patch

    markhalos=False
    nmark=3000
    if(markhalos) :
      with open(halofile, 'r') as file:
        nhalo = 0 
        for line in file:
          nhalo = nhalo + 1
        print("nhalo, ", nhalo, nmark)
   
      with open(halofile, 'r') as file:
        ihalo = 0   
        for line in file:
          ihalo = ihalo + 1
          if(ihalo <= min(nmark, nhalo) ) :
            xh = float(line.split()[1])
            yh = float(line.split()[2])
            zh = float(line.split()[3])  
            rh = float(line.split()[4]) / 1000.0 / boxlen_comov 
            Nh = float(line.split()[6])
            halo_mass = float(line.split()[5])
            if(Nh > 400  ) :
              #print("Nh/numpart :", ihalo, Nh/numpart, Nh)
              xcirc      = (xh - c_arr[0]) * boxlen_comov                          
              ycirc      = (yh - c_arr[1]) * boxlen_comov                          
              zcirc      = (zh - c_arr[2]) * boxlen_comov  
              Lh  = boxlen_comov
              Lhb = boxlen_comov/2.0
              if(xcirc > Lhb) :
                xcirc = xcirc - Lh
              if(xcirc < -Lhb) :
                xcirc = xcirc + Lh 
              if(ycirc > Lhb) :
                ycirc = ycirc - Lh
              if(ycirc < -Lhb) :
                ycirc = ycirc + Lh 
              if(zcirc > Lhb) :
                zcirc = zcirc - Lh 
              if(zcirc < -Lhb) :
                zcirc = zcirc + Lh       
              rcirc      = rh * boxlen_comov    
              # draw circle without changing limits
              color_arr = [ "red", "yellow", "green", "blue", "purple",  "orange", "orange" , "m"]
              #print("marking halo :", ihalo, xh, yh, zh, rcirc)
              if(ihalo == 1) :
                circ = Circle((xcirc, ycirc), 5.0*rcirc, fill=False, linewidth=2, color="red")
              else :  
                circ = Circle((xcirc, ycirc), 5.0*rcirc, fill=False, linewidth=2, color="white")
                
              ax.add_patch(circ)
              ax.text(xcirc, ycirc, f"{ihalo}", color="white", fontsize=8,
                      ha='center', va='center', weight='bold')

            # prevent autoscaling caused by patch

    plotsquare = False
    if(plotsquare ) :
      Lzoom_comov = 2.9
      xleft   =   - Lzoom_comov/2.0
      yleft   =   - Lzoom_comov/2.0
      # keep current limits
      xlim = ax.get_xlim(); ylim = ax.get_ylim()
      # axis units (your colorbar/labels say "Mpc")
      rx = xleft 
      ry = yleft 
      # draw without changing limits; save with matplotlib (not p_big.save)
      ax.add_patch(Rectangle((rx, ry), Lzoom_comov, Lzoom_comov, fill=False, linewidth=2))

    plotsinks = False 
    overwrite_sink = True 
    if(plotsinks) :
      xsink = []
      ysink = []
      zsink = []
      msink = []
      nsink=0
      sinkfile = outdir + "/sinkpos.csv"
      if ( not os.path.exists(sinkfile)  or overwrite_sink == True ) :  
        if( os.path.exists(sinkfile)) :
          command = "rm  " + sinkfile  
          print("command :", command)
          os.system(command)
        command = "./read_sink -inp " + outdir + " -out " + sinkfile  
        print("writing : ", sinkfile)
        os.system(command)
      else :
        print(sinkfile, "already exists")

      with open(sinkfile, 'r') as file:
        for line in file:
          mass = float(line.split(',')[1]) * 1.3220385697284638e+17
          posx = float(line.split(',')[2]) 
          posy = float(line.split(',')[3]) 
          posz = float(line.split(',')[4]) 
          if( posx>projleft[0] and posx<projright[0] and posy>projleft[1] and posy<projright[1] and posz>projleft[2] and posz<projright[2] ) :
            sink_posx = (posx - c_arr[0])*boxlen_comov 
            sink_posy = (posy - c_arr[1])*boxlen_comov 
            sink_posz = (posz - c_arr[2])*boxlen_comov

            xsink.append( sink_posx )
            ysink.append( sink_posy ) 
            zsink.append( sink_posz )
            msink.append(mass)
            print("SINK id in this projection :", nsink, float(line.split(',')[0]), sink_posx, sink_posy, sink_posz, mass/1e5)
            nsink=nsink+1

      if(nsink > 0) :
        figh, axh = plt.subplots()
        axh.hist(msink, bins=len(msink), color='b', label=""  ,edgecolor='b', histtype='step')
        imgname_1 = "hist_msink.png"
        axh.set_xlabel("Msink[Msun]")  
        axh.set_ylabel("count")
        #axh.set_title("redshift of sink spawn")
        #ax1.set_xscale('log')
        #ax1.set_yscale('log')
        axh.legend(loc='best')
        figh.gca().invert_xaxis()
        figh.tight_layout()
        print("Saving : ", imgname_1)
        figh.savefig(imgname_1)

        xlim = ax.get_xlim(); ylim = ax.get_ylim()
        #print(f"ax xlim: [{xlim[0]:.6g}, {xlim[1]:.6g}]  (units: cMpc/h)")
        #print(f"ax ylim: [{ylim[0]:.6g}, {ylim[1]:.6g}]  (units: cMpc/h)")

        if projdir == 'z':
            ax.scatter(xsink, ysink, s=40, facecolors='crimson',
                      edgecolors='black', linewidths=1.0, zorder=99)
        elif projdir == 'x':
            ax.scatter(ysink, zsink, s=40, facecolors='black',
                      edgecolors='black', linewidths=1.0, zorder=99)
        elif projdir == 'y':
            ax.scatter(zsink, xsink, s=40, facecolors='black',
                      edgecolors='black', linewidths=1.0, zorder=99)
        else :
          print("No sink in the given region")
        #ax.autoscale(False)
        #ax.set_xlim(*xlim); ax.set_ylim(*ylim)
        #fig.tight_layout() 
    print("saving :", imgname)     
    fig.savefig(imgname, dpi=300)
    # annotate redshift
    image  = Image.open(imgname)
    draw   = ImageDraw.Draw(image)
    Z = 7.0
    myfont = ImageFont.truetype('/usr/share/fonts/truetype/DejaVuSansMono-Bold.ttf', 80)
    draw.text((1900, 100),  f" \n z = {Z:.2f}\n Ngrid = {Ncell} \n chi^2 = 0",font=myfont, fill=(255, 255, 255))
    image.save(imgname)

