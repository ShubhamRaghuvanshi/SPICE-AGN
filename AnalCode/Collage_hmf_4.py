from PIL import Image, ImageEnhance

W  = 1920
H  = 1920

new = Image.new("RGBA", (int(W),int(H)))

img1 = Image.open("pp_DM_particle_mass_dmonly_256_h1.png")
img2 = Image.open("pp_DM_particle_mass_dmonly_256_h2.png")
img3 = Image.open("halo_rho_r_256_h1.png")
img4 = Image.open("halo_rho_r_256_h2.png")


img1 = img1.resize((int(W/2),int(H/2)))
img2 = img2.resize((int(W/2),int(H/2)))
img3 = img3.resize((int(W/2),int(H/2)))
img4 = img4.resize((int(W/2),int(H/2)))

new.paste(img1, (0,0))
new.paste(img2, (int(W/2.0),0))
new.paste(img3, (0,int(H/2)))
new.paste(img4, (int(W/2.0),int(H/2)))

figname = "h1h2.png"
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

