from PIL import Image, ImageEnhance

W  = 1920
H  = 640

new = Image.new("RGBA", (int(W),int(H)))

img1 = Image.open("sfr_vs_mstar_z10.png")
img2 = Image.open("mstar_vs_mhalo_z10.png")
img3 = Image.open("mstar_to_mhalo_z10.png")

#img1 = Image.open("proj_DM_density_star_S2.png")
#img2 = Image.open("proj_nH_star_S2.png")
#img3 = Image.open("proj_temperature_star_S2.png")

#img1 = Image.open("proj_DM_density_smoothSN_S2.png")
#img2 = Image.open("proj_nH_smoothSN_S2.png")
#img3 = Image.open("proj_temperature_smoothSN_S2.png")


img1 = img1.resize((int(W/3),int(H)))
img2 = img2.resize((int(W/3),int(H)))
img3 = img3.resize((int(W/3),int(H)))

new.paste(img1, (0,0))
new.paste(img2, (int(W/3.0),0))
new.paste(img3, (int(2.0*W/3.0),0))

figname = "haloprops_z10.png"
print("saving :", figname) 
new.save(figname)


'''
new = Image.new("RGBA", (int(W),int(H)))
#img1 = Image.open("sink_mass_1.png")
#img2 = Image.open("sink_mass_2.png")
#img3 = Image.open("sink_mass_5.png")

#img1 = Image.open("sinkAGN_mass_1.png")
#img2 = Image.open("sinkAGN_mass_2.png")
#img3 = Image.open("sinkAGN_mass_5.png")

img1 = Image.open("ndens_AGN_9.png")
img2 = Image.open("temp_AGN_9.png")
img3 = Image.open("sink_mass.png")


img1 = img1.resize((int(W/3),int(H)))
img2 = img2.resize((int(W/3),int(H)))
img3 = img3.resize((int(W/3),int(H)))


new.paste(img1, (0,0))
new.paste(img2, (int(W/3),0))
new.paste(img3, (int(2*W/3),0))

figname = "priminary_zoom_AGN.png"
print("saving :", figname) 
new.save(figname)
'''

