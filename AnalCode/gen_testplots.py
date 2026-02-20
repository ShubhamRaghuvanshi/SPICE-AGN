import yt_funcs
import cfuncs


simdir_1 = "/ptmp/mpa/sraghu/smoothSN256L10/nosink"

simtyp_arr     = [  "snagn" ]
simdir_arr     = [  simdir_1]

cbar_arr       = [ 'g'  ]
linestyle_arr  = ['-'  ]
mbar_arr       = ['o' ]


sfr_outnum_arr = [9,9]
halonum_arr    = [-1,-1]
outnum_max     = [-1,-1]
outnum_5       = [9]
outnum_6       = [9]
outnum_arr     = [outnum_5, outnum_6]
Nout           = len(outnum_6)


'''
simdir_5 = "/ptmp/mpa/sraghu/zoom/bigbox"

simtyp_arr     = [  "big_box" ]
simdir_arr     = [  simdir_5, ]

cbar_arr       = [ 'm'  ]
linestyle_arr  = ['-'  ]
mbar_arr       = ['o' ]


sfr_outnum_arr = [17,17]
halonum_arr    = [-1,-1]
outnum_max     = [-1,-1]
outnum_5       = [6]
outnum_arr     = [outnum_5]
Nout           = len(outnum_5)
'''


'''
#simdir_1 = "/ptmp/mpa/sraghu/smooth/ref_opt/star_sn_sink"  #16
#simdir_2 = "/ptmp/mpa/sraghu/smooth/ref_opt/star_sn_sink_avgacc" #19


simdir_1 = "/ptmp/mpa/sraghu/smooth/sink_drag/star_sink" #10
simdir_2  = "/u/sraghu/sim/SHUB/b256/star_sink_avgacc" #19

#simdir_3 = "/u/sraghu/sim/SHUB/b256/star_sinkAGN" #19
# simdir_4 = "/u/sraghu/sim/SHUB/b256/star_sinkAGN_avgacc" #19

simdir_3 = "/ptmp/mpa/sraghu/smooth/ref_opt/star_sn_sink" #16
simdir_4 = "/ptmp/mpa/sraghu/smooth/ref_opt/star_sn_sink_avgacc" #19

simdir_5 = "/ptmp/mpa/sraghu/smooth/ref_opt/star_sn_sink_rt"
simdir_6 = simdir_4

simtyp_arr     = [ "star-sink", "star-sink2(acc_mode=2)", "smoothSN_sink", "smoothSN_sink2", "smoothSN_sink_RT" , "smoothSN_sink2_RT"]
simdir_arr     = [ simdir_1, simdir_2, simdir_3, simdir_4, simdir_5, simdir_6]

cbar_arr       = ['y','y', 'g', 'g', 'b','b','b','b'    ]
linestyle_arr  = ['-', '--','-', '--','-', '--'  ]
mbar_arr       = ['*'  ,'*' ,'*', '*','*', '*'   ]
sfr_outnum_arr = [11,15,13,13]
halonum_arr    = [-1,-1,-1,-1]
outnum_max     = [-1,-1,-1,-1]
outnum_1       = [11]
outnum_2       = [13]
outnum_arr     = [outnum_1, outnum_2]
Nout           = len(outnum_1)
'''

'''
simdir_1       = "/u/sraghu/sim/SHUB/b256/L100/dmonly/AGNtest"
#simdir_1       = "/ptmp/mpa/sraghu/ref"
simtyp_arr     = ["ref"]
simdir_arr     = [ simdir_1]
cbar_arr       = [ 'b'  , 'r' ]
linestyle_arr  = ['-'   , '-'       ]
mbar_arr       = ['*'   , '*'     ]
sfr_outnum_arr = [3]
halonum_arr    = [-1]
outnum_max     = [-1,-1]
outnum_1       = [4]
outnum_arr     = [outnum_1]
Nout           = len(outnum_1)
'''

#prop_arr      = ['sfr', 'hmf' , 'dmdens_proj' , 'dens_proj', 'temp_proj' , 'phaseplot']#,  'globalsink_vs_z',  'haloprops', 'simtime']
prop_arr = ['globalsink_vs_z']



for prop in prop_arr :
  print('######################## prop :', prop, "##########################")
  if(prop=='sfr') :
    #compare SFR
    plotanal=False
    plotdata=True
    overwrite=False
    physical = True
    fullbox = True
    #write star files
    iout=0
    for simdir in simdir_arr:
      outnum = sfr_outnum_arr[iout]
      if(halonum_arr[iout] > 0) :
        halo_attrb = yt_funcs.give_haloprops(simdir, outnum, halonum_arr[iout], physical, False)
        center = halo_attrb[0]
        radius = halo_attrb[1]
        mass   = halo_attrb[2] 
        fullbox = False     
      else :
        center = [-1,-1,-1]
        radius = -1  
        mass = -1e9
      print("halo center, radius, mass : ", center, radius, mass/1e9)  
      yt_funcs.write_starfile(simdir, outnum, overwrite, radius, center)
      iout = iout+1
    # plot sfr  
    yt_funcs.compare_sfr(simdir_arr, simtyp_arr ,sfr_outnum_arr, cbar_arr, linestyle_arr ,plotanal, plotdata, fullbox )
  if(prop=='hmf') :
    # write halo catalog 
    overwrite=False
    Boxlength_comov=100.0
    for iout in range(0,Nout) :
      isim=0
      hmf_outnum_arr = []
      for simdir in simdir_arr:
        outnum = outnum_arr[isim][iout]  
        hmf_outnum_arr.append(outnum)
        yt_funcs.create_halo_catalog(simdir, outnum, 0, overwrite)
        isim = isim + 1
      fitfunc_arr = [2]
      #fitting_function (int) Which fitting function to use. 1 = Press-Schechter, 2 = Jenkins, 3 = Sheth-Tormen, 4 = Warren, 5 = Tinker Default : 4.
      print('hmf_outnum_arr : ', hmf_outnum_arr)
      yt_funcs.compare_hmf(simdir_arr, simtyp_arr, cbar_arr, mbar_arr, hmf_outnum_arr, fitfunc_arr, Boxlength_comov)
  if(prop=='dmdens_proj') :
    boxlen_comov=6.77
    overwrite=False
    physical = False
    fullbox = True
    projdir = 'y'
    
    for iout in range(0,Nout) :
      isim=0
      hist_outnum_arr = []
      for simdir in simdir_arr:
        outnum = outnum_arr[isim][iout]
        simtyp = simtyp_arr[isim]
        hist_outnum_arr.append(outnum)
        if(fullbox or halonum_arr[iout] < 0) :
          radius = -1
          center = [0.5, 0.5, 0.5]
        else :
          halo_attrb = yt_funcs.give_haloprops(simdir, outnum, halonum_arr[iout], physical, False)
          center = halo_attrb[0]
          radius = halo_attrb[1]
          aexp   = halo_attrb[3]
          #if(physical) :
          #  radius = radius / aexp  
        #yt_funcs.plot_proj_dmdens(simdir, outnum, simtyp, radius, center, projdir, boxlen_comov)
        isim = isim + 1
      yt_funcs.plot_hist_dmddens(simdir_arr, hist_outnum_arr, simtyp_arr, cbar_arr)  
  if(prop=='dens_proj') :
    boxlen_comov=6.77
    overwrite=False
    physical = False
    fullbox = False
    projdir = 'y'
    
    for iout in range(0,Nout) :
      isim=0
      for simdir in simdir_arr:
        outnum = outnum_arr[isim][iout]
        simtyp = simtyp_arr[isim]
        if(fullbox or halonum_arr[iout] < 0) :
          radius = -1
          center = [0.5, 0.5, 0.5]
        else :
          halo_attrb = yt_funcs.give_haloprops(simdir, outnum, halonum_arr[iout], physical, False)
          center = halo_attrb[0]
          radius = halo_attrb[1]
          aexp   = halo_attrb[3]
          #if(physical) :
          #  radius = radius / aexp  
        yt_funcs.plot_proj_gasdens(simdir, outnum, simtyp, radius, center, projdir, boxlen_comov)
        isim = isim + 1
  if(prop=='temp_proj') :
    boxlen_comov=6.77
    overwrite=True
    physical = False
    fullbox = True
    projdir = 'z'
    
    for iout in range(0,Nout) :
      isim=0
      for simdir in simdir_arr:
        outnum = outnum_arr[isim][iout]
        simtyp = simtyp_arr[isim]
        if(fullbox or halonum_arr[iout] < 0) :
          radius = -1
          center = [0.5, 0.5, 0.5]
        else :
          halo_attrb = yt_funcs.give_haloprops(simdir, outnum, halonum_arr[iout], physical, False)
          center = halo_attrb[0]
          radius = halo_attrb[1]
          aexp   = halo_attrb[3]
          #if(physical) :
          #  radius = radius / aexp  
        yt_funcs.plot_proj_gastemp(simdir, outnum, simtyp, radius, center, projdir, boxlen_comov)
        isim = isim + 1
  if(prop=='metalicity_proj') :
    boxlen_comov=6.77
    overwrite=False
    for iout in range(0,Nout) :
      isim=0
      for simdir in simdir_arr:
        outnum = outnum_arr[isim][iout]
        simtyp = simtyp_arr[isim]
        yt_funcs.plot_proj_gastemp(simdir, outnum, simtyp, boxlen_comov, overwrite)
        isim = isim + 1
  if(prop=='phaseplot') :
    boxlen_comov=6.77
    overwrite=False
    for iout in range(0,Nout) :
      isim=0
      for simdir in simdir_arr:
        outnum = outnum_arr[isim][iout]
        simtyp = simtyp_arr[isim]
        yt_funcs.tn_phaseplot(simdir, outnum, simtyp, boxlen_comov, overwrite)
        isim = isim + 1
  if(prop=='globalsink_vs_z')  :
    boxlen_comov = 6.774
    #cfuncs.compare_sinkmasstot(simdir_arr, simtyp_arr, cbar_arr, linestyle_arr, boxlen_comov)
    #cfuncs.plot_sink_pos_vs_z(simdir_arr, cbar_arr, simtyp_arr, 4, boxlen_comov, True, 100)
    #cfuncs.plot_sink_pos_vs_z(simdir_arr, cbar_arr, simtyp_arr, boxlen_comov, 100)
    #cfuncs.plot_sink_relpos_vs_z(simdir_arr, cbar_arr, simtyp_arr, 9, boxlen_comov, 100, outnum_max)
    #cfuncs.hist_z_sinkspawn(simdir_arr, simtyp_arr, cbar_arr)

    #isink_arr = cfuncs.give_massive_sinkid(simdir_arr)
    isink_arr=  [1]
    #print("isink_arr :", isink_arr)
    cfuncs.sink_m_dmdt(simdir_arr, isink_arr, simtyp_arr, cbar_arr, linestyle_arr,  boxlen_comov,1)
      #cfuncs.plot_sinkmassonly(simdir_arr, isink, simtyp_arr, cbar_arr, boxlen_comov)

  if(prop=='haloprops') :
    overwrite_sink = False
    overwrite_star = False
    plotsinks      = False
    boxlen_comov = 6.774
   
    for iout in range(0,Nout) :
      isim = 0
      haloprop_outnum_arr = []
      for simdir in simdir_arr:
        outnum = outnum_arr[isim][iout] 
        haloprop_outnum_arr.append(outnum)
        isim = isim + 1
      print('haloprop_outnum_arr :', haloprop_outnum_arr) 
      cfuncs.haloprops(simdir_arr, simtyp_arr, haloprop_outnum_arr, cbar_arr, mbar_arr, linestyle_arr, 40000, boxlen_comov, overwrite_star, overwrite_sink, plotsinks)
  if(prop=='simtime') :
    #cfuncs.read_runtime(simdir_arr, simtyp_arr, cbar_arr, linestyle_arr)
    cfuncs.read_finestep(simdir_arr, simtyp_arr, cbar_arr)  



   








