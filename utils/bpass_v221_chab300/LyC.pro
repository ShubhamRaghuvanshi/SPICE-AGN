

	savefile = 'phr_BPASSv2_100.sav'
 
	; user-supplied ionisation bins
	ev_HI = 13.60
	ev_HeI = 24.59
	ev_HeII = 54.42
	ev_LW = 11.2

	ev_ion = [ev_HI, ev_HeI, ev_HeII, ev_LW]
	ev_ion_sta = ev_ion
	ev_ion_end = [ev_HeI, ev_HeII, 1d6, ev_HI]

	nion = n_elements(ev_ion)

	ytit = strarr(nion)
	ytit[0] = 'N_{ph,HI} [#/s/Msun]'
	ytit[1] = 'N_{ph,HeI} [#/s/Msun]'
	ytit[2] = 'N_{ph,HeII} [#/s/Msun]'
	ytit[3] = 'N_{ph,LW} [#/s/Msun]'

	ytit2 = strarr(nion)
	ytit2[0] = 'tot N_{ph,HI} [#/Msun]'
	ytit2[1] = 'tot N_{ph,HeI} [#/Msun]'
	ytit2[2] = 'tot N_{ph,HeII} [#/Msun]'
	ytit2[3] = 'tot N_{ph,LW} [#/Msun]'
	output_file1 = 'LyC.eps'
	output_file2 = 'LyC_cum.eps'
	output_file3 = 'LyC_cum2.eps'


	;=================================================
	; no need to touch from here
	;=================================================
	Lsun = 3.8256d33
	c_cms = 2.9980000d+10 ; cm / s
	h_pl = 6.6260755d-27 ; erg s
	ev2erg = 1.6021772d-12 ; erg
	yr2s = 365d0*24*3600

	openr,lun,'age_bins.dat',/get_lun
	nage=0L
	readf,lun,nage
	age = dblarr(nage)
	readf,lun,age
	close,lun
	free_lun,lun


	openr,lun,'metallicity_bins.dat',/get_lun
	nmet = 0L
	readf,lun,nmet
	metal = dblarr(nmet)
	readf,lun,metal
	close,lun
	free_lun,lun


	openr,lun,'all_seds.dat',/get_lun, /f77_unfor
	nwave=0L
	readu,lun,nwave
	wave=dblarr(nwave)
	readu,lun,wave
	sed = dblarr(nwave,nage,nmet)
	tmp = dblarr(nwave)
	for iz=0,nmet-1 do begin
	for it=0,nage-1 do begin
		readu,lun,tmp
		sed[*,it,iz]=tmp*Lsun ; [erg/s/A]
		; because tmp is normalised to Lsolar
	endfor
	endfor

	erg_ion_sta = dblarr(nion)
	erg_ion_end = dblarr(nion)
	wave_ion_sta = dblarr(nion)
	wave_ion_end = dblarr(nion)
	iwave_min = dblarr(nion)
	iwave_max = dblarr(nion)
	for i=0,nion-1 do begin
		erg_ion_sta [i] = ev_ion_sta[i]*ev2erg
		erg_ion_end [i] = ev_ion_end[i]*ev2erg
		wave_ion_sta[i] = h_pl * c_cms / erg_ion_end[i] * 1d8 ; Ang ; needed to reverse 
		wave_ion_end[i] = h_pl * c_cms / erg_ion_sta[i] * 1d8 ; Ang ; needed to reverse

		dum = min( abs(wave-wave_ion_end[i]), imin)
		iwave_max[i] = imin
		dum = min( abs(wave-wave_ion_sta[i]), imin)
		iwave_min[i] = imin

		print, iwave_min[i], iwave_max[i]
	endfor

	nph = dblarr(nage,nmet,nion)
	cnph = dblarr(nage,nmet,nion)
	
	lam = wave * 1d-8 ; wave is in A	
	nu = c_cms/lam

	for iz=0,nmet-1 do begin
		tprev = 0d0
		nprev = dblarr(nion)
		for it=0,nage-1 do begin
			L_lam = sed[*,it,iz] ; [erg/s/A]
			L_lam = L_lam * 1d8  ; [erg/s/cm]
			
			; L_lam dlam = L_nu dnu
			; L_nu = L_lam c /nu^2 = L_lam lam /nu
			L_nu = L_lam*lam/nu
			y    = L_nu/(h_pl*nu) ; number of photons of this spectrum
   
			; for cumulative quantities
			dt = (age[it]-tprev)*yr2s ; [s]		
	
			; count number of photons between iwave_min and iwave_max
			; Np = Fnu / (h*nu)
			for j=0, nion-1L do begin
				ista = iwave_min[j]	
				iend = iwave_max[j]
   
				dnu  =  nu[ista:iend-1] - nu[ista+1:iend]
				yavg = ( y[ista:iend-1] +  y[ista+1:iend])/2d0
				nph[it,iz,j] = total(yavg*dnu,/double)
				cnph[it,iz,j] = nprev[j] + nph[it,iz,j]*dt	

				;if iz eq 0 and it eq 62 then stop

				if nph[it,iz,j] eq 0 then nph[it,iz,j] += 1d-20		
				nprev[j] = cnph[it,iz,j]	
				tprev = age[it]				
			endfor
   
		endfor
	endfor

	tgrid = age
	tgrid[0] = 1
	zgrid = metal
	phrdata = nph
	save, tgrid,zgrid,phrdata,nph,cnph,  file=savefile

	; plot
   SET_PLOT,'PS'
   !P.FONT=0
   !X.THICK=3
   !Y.THICK=3
   !P.THICK=3
   !P.CHARTHICK=3
   !P.CHARSIZE=2
	!P.MULTI=[0,1,4]

	yr = [1d40,1d47]
	xr = [1d4,max(age)]
	device,file=output_file1,/encap,/color,xs=14,ys=18
	loadct,33
	dum = getcolor('black',0)
	dum = getcolor('white',255)

	x0 = 0.14
	x1 = 0.96
	y0 = 0.07
	y1 = 0.96
	ypad = 0.05
	dy = (y1-y0-ypad*(nion-1))/float(nion)
	
	ytit = textoidl(ytit)

	for iion=0,nion-1 do begin
		yl = y0 + (nion-iion-1)*(dy+ypad)
		yh = yl + dy
		pos = [x0, yl, x1, yh]
		if iion eq nion-1 then xtit='time [yr]'
		plot, xr, yr, /nodata, xr=xr,yr=yr,/xlog,/ylog,/xs,/ys,pos=pos,ytit=ytit[iion],xtit=xtit
		for iz=0,nmet-1 do begin
			col = (iz+1)/float(nmet)*254
			oplot, age, nph[*,iz,iion], col=col
		endfor
	endfor

	cols = (findgen(nmet)+1)/nmet*254
	items = 'Z='+string(metal,format='(f6.4)')
	pos = [1d7,1d46]
	legend, items, col=cols,/short,spacing=1.6,charsize=1.2,linestyle=0,pos=pos,thick=4,box=0

	device,/close

	; part 2
	yr = [1d58,1d62]
	device,file=output_file2,/encap,/color,xs=14,ys=18
	loadct,33
	dum = getcolor('black',0)
	dum = getcolor('white',255)

	x0 = 0.14
	x1 = 0.96
	y0 = 0.07
	y1 = 0.96
	ypad = 0.05
	dy = (y1-y0-ypad*(nion-1))/float(nion)
	
	ytit = textoidl(ytit2)
	xtit = ''
	for iion=0,nion-1 do begin
		yl = y0 + (nion-iion-1)*(dy+ypad)
		yh = yl + dy
		pos = [x0, yl, x1, yh]
		if iion eq nion-1 then xtit='time [yr]'
		plot, xr, yr, /nodata, xr=xr,yr=yr,/xlog,/ylog,/xs,/ys,pos=pos,ytit=ytit[iion],xtit=xtit
		for iz=0,nmet-1 do begin
			col = (iz+1)/float(nmet)*254
			oplot, age, cnph[*,iz,iion], col=col
		endfor
	endfor

	device,/close

	; part 3
	yr = [0,12]
	xr = [1d5,1d10]
	device,file=output_file3,/encap,/color,xs=14,ys=14
	loadct,33
	dum = getcolor('black',0)
	dum = getcolor('white',255)

	x0 = 0.14
	x1 = 0.96
	y0 = 0.07
	y1 = 0.96
	ypad = 0.05
	dy = (y1-y0-ypad*(nion-1))/float(nion)
	
	ytit = 'cum N_{LyC} [10^{60} #/Msun]'
	ytit = textoidl(ytit)
	xtit='time [yr]'
	pos = [x0, y0, x1, y1]
		
	plot, xr, yr, /nodata, xr=xr,yr=yr,/xlog,/xs,/ys,pos=pos,ytit=ytit,xtit=xtit
	for iz=0,nmet-1 do begin
		col = (iz+1)/float(nmet)*254
		yy  = cnph[*,iz,0] + cnph[*,iz,1] + cnph[*,iz,2]
		oplot, age, yy/1d60, col=col
	endfor
	legend, items, col=cols,/short,spacing=1.6,charsize=0.8,linestyle=0,thick=5,box=0

	device,/close




end
