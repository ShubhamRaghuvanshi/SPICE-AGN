
;README: Spectral Energy Distributions
;These files contain the primary output of BPASS, which is the stellar spectral energy
;distribution (SED). Flux values are given every 1 Angstrom in the range 1 - 100,000 A.
;Most users will wish to resample to lower resolution, depending on their use case. We
;caution that some of the stellar atmospheres we use also have slightly lower spectral
;resolution.
;Each file has 52 columns and 106 rows. The first column lists a wavelength in angstroms,
;and each remaining column n (n>1) holds the model flux for the population at an age of
;10^(6+0.1*(n-2)) years at that wavelength.
;The units of flux are Solar Luminosities per Angstrom, normalised for a cluster of 1e6
;Msun formed in a single instantaneous burst. The total luminosity of the SED can be
;simply calculated by summing all the rows together. Or the flux over the wavelength
;range from, for example, 2000 to 3000 Angstroms can be calculated by summing the
;2000th to 3000th rows.
;Filenames: spectra-<opt>-<imf>.<z>.dat



	path = './'

	cd, path, curr=curr_dir

	; for single stars
	filelist = file_search('*.gz')

	for i=0L,n_elements(filelist)-1L do begin
		file = filelist[i]
		command = 'gunzip '+file
		print, command
		spawn, command
	endfor



end	
