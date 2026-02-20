import yt_funcs
import cfuncs

boxlen_comov = 6.774 
sidelen = -1
center = [0.5,0.5,0.5]

#simdir_1  = "/u/sraghu/sim/SHUB/b256/unibox/ref_0"
#simtyp_1  = "SmoothSN_A"

#simdir_2  = "/u/sraghu/sim/SHUB/b256/unibox/ref_0_lmax17"
#simtyp_2  = "SmoothSN_S2_lmax17"

#simdir_3  = "/u/sraghu/sim/SHUB/b256/unibox/ref_0"
#simtyp_3  = "SmoothSN_S2_lmax18"

#simdir_4  = "/u/sraghu/sim/SHUB/b256/unibox/ref_0_lmax19"
#simtyp_4  = "SmoothSN_S2_lmax19"

#simdir_5  = "/u/sraghu/sim/SHUB/b256/unibox/ref_0_lmax20"
#simtyp_5  = "SmoothSN_S2_lmax20"


#simdir_1  = "/ptmp/mpa/sraghu/smoothSN256L10/nosink" #10cMpc/h

simdir_1  = "/u/sraghu/sim/SHUB/b256/unibox/ref_0_lmax17"
simdir_2  = "/u/sraghu/sim/SHUB/b256/unibox/ref_0"
simdir_3  = "/u/sraghu/sim/SHUB/b256/unibox/ref_0_lmax19"
simdir_4  = "/u/sraghu/sim/SHUB/b256/unibox/ref_0_lmax20"
simdir_5  = "/ptmp/mpa/sraghu/unibox/star_0"

#simdir_2 = "/ptmp/mpa/sraghu/smooth/sink_drag/star_sink"
#simdir_2  = "/ptmp/mpa/sraghu/smoothSN256L10/nosink" 
#simdir_1  = "/ptmp/mpa/sraghu/unibox/ref_0_lmax17"


simtyp_1  = "SmoothSN_S2_lmax17"
simtyp_2  = "SmoothSN_S2_lmax18"
simtyp_3  = "SmoothSN_S2_lmax19"
simtyp_4  = "SmoothSN_S2_lmax20"

outnum_arr    = [8,8, 9, 9]
#outnum_arr    = [10,9]
#outnum_arr    = [8,8]
pltclr_arr    = ['orange','g','purple','r']
linestyle_arr = ['-', '-', '-', '-', '-']
marker_arr    = ['*', '*', '*', '*', '*']
#simdir_arr = [simdir_1, simdir_2, simdir_3, simdir_4, simdir_5] 
#simtyp_arr = [simtyp_1, simtyp_2, simtyp_3, simtyp_4, simtyp_5]

simdir_arr = [simdir_1, simdir_2, simdir_3, simdir_4] #simdir_2, simdir_3, simdir_4, simdir_5] 
simtyp_arr = [simtyp_1, simtyp_2, simtyp_3, simtyp_4] #simtyp_2, simtyp_3, simtyp_4, simtyp_5]


prop_arr = ['haloprops']

for propval in prop_arr:
  if(propval == "sfrd") :
    #compare SFR
    plotanal=False
    plotdata=True
    overwrite=True
    fullbox = True
    #write star files
    isim = 0 
    for simdir in simdir_arr:
      yt_funcs.write_starfile(simdir, outnum_arr[isim], center, sidelen , fullbox, overwrite)
      isim = isim + 1
    # plot sfr  
    yt_funcs.compare_sfr(simdir_arr, simtyp_arr ,outnum_arr, pltclr_arr, linestyle_arr ,plotanal, plotdata )
  elif(propval =='amrres') :
    #yt_funcs.compare_amr(simdir_arr, simname_arr, outnum_arr, cbar_arr, linestyle_arr)
    yt_funcs.compare_amr(simdir_arr, simtyp_arr ,outnum_arr, pltclr_arr, linestyle_arr )   
  elif(propval =='haloprops') :
    overwrite_sink = False
    overwrite_star = False
    plotsinks      = False
    cfuncs.haloprops(simdir_arr, simtyp_arr, outnum_arr, pltclr_arr, marker_arr, linestyle_arr, boxlen_comov,  40000, overwrite_star, overwrite_sink, plotsinks)
  elif(propval =='hmf') :
    # write halo catalog 
    overwrite=False
    fitfunc = 3
    yt_funcs.plot_hmf(simdir_arr, simtyp_arr, outnum_arr, pltclr_arr, fitfunc, boxlen_comov)