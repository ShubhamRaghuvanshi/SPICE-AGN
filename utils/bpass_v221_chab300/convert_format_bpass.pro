pro convert_format_bpass, correct=correct, show=show, ps=ps
	
	correct = keyword_set(correct)

	show = keyword_set(show)
	do_ps = keyword_set(ps)



	; note that BPASS_v2.2.1 has too many wavelength info after > 5000 A
	; if we simply convert the format for RAMSES-RT, the information alone is ~500MB
	; therefore we decide to sample the small number of wavelength bins, as in BC03
	reference_direc = '~/soft/lib/SED/bc03_Chabrier/'
	reference_sed   = reference_direc+'all_seds.dat'
	if correct and ~file_test(reference_sed) then begin
		print, '>> ERR: no such file '+reference_file
		retall
	endif
	;read the wavelength bin for bc03_Chabrier
	openr,2,reference_sed,/f77
	nWave_bc03 = 0L
	readu,2,nWave_bc03
	wave_bc03 = dblarr(nWave_bc03)
	readu,2,Wave_bc03
	close,2



	; Each file has 52 columns and 100000 rows. The first column is the wavelength of the SED and 
	; each subsequent n-th column has the flux for the population at an age of 10^(6+0.1*(n-2)) years.
	; Each row than the flux at the wavelength given in the first column, 
	; the units of flux are Solar Luminosities per Angstrom for either a cluster of 1e6Msun or 1Msun/yr of star-formation.

	path = './'

	output_direc = './'	
	output_met = output_direc+'metallicity_bins.dat'
	output_age = output_direc+'age_bins.dat'
	output_sed = output_direc+'all_seds.dat'

	if ~file_test(output_direc) then file_mkdir, output_direc

	cd, path, curr=curr_dir

	; for binary stars
	filelist = file_search('spectra-bin-*.dat')
	; for single stars
	;filelist = file_search('spectra.z???.dat')

	; This is from bpass_v2.2.1
	; https://drive.google.com/file/d/1IYCYf5Bxt1WmqPuFTYLQ7kpN-hKY2SAp/view	
	Zsun = 0.02
	Lsun = 3.848d26 * 1d7 ; erg/s
	M_cluster = 1d6      ; for BPASS_v2

	nWave = 100000L

	nAge  = 51 ; 41 -> 51 (bpsss_v2 -> bpass_v2.2.1)
	age_yr = 10d0^(6+0.1*(lindgen(nAge)))

	cd, curr_dir

	nZ = n_elements(filelist)
	sZgrid = strarr(nZ)
	Zgrid = dblarr(nZ)

	for iz=0,nZ-1 do begin
		filename = filelist[iz]
		dum = strsplit(filename,'.',/extract)
		sZgrid[iz] = dum[1]
		dum2 = strsplit(dum[1],'z',/extract)

		if strmatch(dum2,'em?') then begin
			dum3 = strsplit(dum2,'em',/extract)
			Zgrid[iz] = 10d0^('-'+dum3[0])
		endif else begin
			sz = '0.'+dum2[0]
			Zgrid[iz] =  double(sz)
		endelse
	endfor


	;sorting is necessary
	iso = sort(Zgrid)
	if n_elements(iso) ne n_elements(filelist) then begin
		print, 'cannot be sorted bc n_elements(zgrid) != n_element(filelist)'
		stop
	endif
	filelist = filelist[iso]
	Zgrid = Zgrid[iso] 	

	;------------------------------------------------------
	; save metallicity grid	
	;------------------------------------------------------
	openw,1,output_met
	printf,1,fix(nZ)
	for iz=0,nZ-1 do begin
		printf,1,string(Zgrid[iz],format='(E14.6)')
	endfor
	close,1


	wave = dblarr(nWave)
	sed = dblarr(nWave,nAge,nZ)

	; read the data 
	for iz=0,nZ-1L do begin
		filename = path+filelist[iz]
		print, 'processing...'+filename
		openr,1,filename
		text = ''
		for iw = 0L, nWave-1L do begin
			readf,1,text
			text_split = strsplit(text,/extract)
			wave[iw] = double(text_split[0])
			sed [iw,*,iz]= double(text_split[1:nAge])
		endfor
		close,1
	endfor

	sed = sed  / M_cluster ; Lsun / A / Msun
	


	; 1) Change the wavelength bins
	ok_short = where(wave lt 5000, nwave_short)
	ok_long  = where(wave_bc03 gt 5000, nwave_long)

	nWave_reduced = nwave_short + nwave_long
	wave_reduced = dblarr(nWave_reduced)
	wave_reduced[0:nwave_short-1L] = wave[ok_short]
	wave_reduced[  nwave_short:nWave_reduced-1L] = wave_bc03[ok_long]
	
	sed_reduced = dblarr(nWave_reduced,nAge,nZ)
	for it=0L,nAge-1L do begin
		for iz=0L,nZ-1L do begin
			if min(sed[*,it,iz]) le 0 then stop
			reduced_spectrum = 10d0^interpol( alog10(sed[*,it,iz]),alog10(wave),alog10(wave_reduced) )
			sed_reduced[*,it,iz] = reduced_spectrum
		endfor
	endfor


	; 2) Add the first set of SED if age_yr[0] ne 0
	if age_yr[0] gt 0 then begin

		sed_new = dblarr(nWave_reduced,nAge+1,nZ)
		for it=0,nAge-1L do begin
			for iz=0,nZ-1L do begin
				sed_new [*,it+1,iz] = sed_reduced[*,it,iz]
			endfor
		endfor
		; make t=0 bin equivalent to the first age bin from BPASS
		for iz=0,nZ-1L do begin
			sed_new[*,0,iz] = sed_reduced[*,1,iz]
		endfor
		nAge_new  = nAge + 1 ; 51 -> 52
		age_yr_new = [0.0, age_yr]
	endif


	if show then plot,wave_corr, sed_new[*,nAge-2,0], xr=[100,1d5],/xs,yr=[1d-5,1d-2],/ys, /xlog,/ylog, linestyle=2
	
	if show then begin
		for kk=4,1,-1 do begin
			oplot, wave,sed_reduced[*,nAge-kk,0],color=getcolor('orange',2)
			wait,0.1
		endfor
	endif

	;------------------------------------------------------
	; save age grid	
	;------------------------------------------------------
	openw,10,output_age
	printf,10,fix(nAge_new)
	for it=0,nAge_new-1 do begin
		printf,10,string(age_yr_new[it],format='(E14.6)')
	endfor
	close,10


	;------------------------------------------------------
	; save sed	
	;------------------------------------------------------
	openw,10,output_sed,/f77_unformatted
	writeu,10,long(nWave_reduced),long(nAge_new)
	writeu,10,double(wave_reduced)
	for iZ=0,nZ-1 do begin
		for it=0,nAge_new-1L do begin
			writeu,10, double(sed_new[*,it,iZ])
		endfor
	endfor
	close,10

	if do_ps then begin
		loadct,33
		dum = getcolor('black',0)
		dum = getcolor('black',255)
		set_plot,'ps'
		!p.charsize=1.2
		!p.thick=3	
		device,file='SEDs.ps',/color,xs=14,ys=14
		for iz=0,nZ-1 do begin
			xr = [10,1d5]
			yr = [1d-6,2]
			plot, xr, yr, xr=xr,/xs,yr=yr,/ys, /xlog,/ylog, /nodata
			for it=0,40,10 do begin
				icol = (it+1)/float(nAge_corr)*254
				oplot,wave,sed_new[*,it,iZ],color=icol
			endfor
			legend,/right,sZgrid[iz],/bottom
		endfor
		device,/close
	endif


end	
