from PIL import Image, ImageEnhance

W  = 1080
H  = 1920

new = Image.new("RGBA", (int(W),int(H)))

img1 = Image.open("pp_DM_particle_mass_dmonly_256_1.png")
img2 = Image.open("pp_DM_particle_mass_dmonly_256_2.png")
img3 = Image.open("pp_DM_particle_mass_dmonly_256_3.png")
img4 = Image.open("pp_DM_particle_mass_dmonly_256_4.png")
img5 = Image.open("pp_DM_particle_mass_dmonly_256_5.png")
img6 = Image.open("pp_DM_particle_mass_dmonly_256_6.png")
img7 = Image.open("pp_DM_particle_mass_dmonly_256_7.png")
img8 = Image.open("pp_DM_particle_mass_dmonly_256_8.png")
img9 = Image.open("pp_DM_particle_mass_dmonly_256_9.png")

#img1 = Image.open("pp_DM_particle_mass_dmonly_512_10.png")
#img2 = Image.open("pp_DM_particle_mass_dmonly_512_11.png")
#img3 = Image.open("pp_DM_particle_mass_dmonly_512_12.png")
#img4 = Image.open("pp_DM_particle_mass_dmonly_512_13.png")
#img5 = Image.open("pp_DM_particle_mass_dmonly_512_13.png")
#img6 = Image.open("pp_DM_particle_mass_dmonly_512_15.png")
#img7 = Image.open("pp_DM_particle_mass_dmonly_512_16.png")
#img8 = Image.open("pp_DM_particle_mass_dmonly_512_17.png")
#img9 = Image.open("pp_DM_particle_mass_dmonly_512_18.png")

#img1 = Image.open("pp_DM_particle_mass_dmonly_512_19.png")
#img2 = Image.open("pp_DM_particle_mass_dmonly_512_19.png")
#img3 = Image.open("pp_DM_particle_mass_dmonly_512_21.png")
#img4 = Image.open("pp_DM_particle_mass_dmonly_512_22.png")
#img5 = Image.open("pp_DM_particle_mass_dmonly_512_23.png")
#img6 = Image.open("pp_DM_particle_mass_dmonly_512_24.png")
#img7 = Image.open("pp_DM_particle_mass_dmonly_512_24.png")
#img8 = Image.open("pp_DM_particle_mass_dmonly_512_26.png")
#img9 = Image.open("pp_DM_particle_mass_dmonly_512_27.png")


img1 = img1.resize((int(W/3),int(H/3)))
img2 = img2.resize((int(W/3),int(H/3)))
img3 = img3.resize((int(W/3),int(H/3)))

img4 = img4.resize((int(W/3),int(H/3)))
img5 = img5.resize((int(W/3),int(H/3)))
img6 = img6.resize((int(W/3),int(H/3)))

img7 = img7.resize((int(W/3),int(H/3)))
img8 = img8.resize((int(W/3),int(H/3)))
img9 = img9.resize((int(W/3),int(H/3)))


new.paste(img1, (0,0))
new.paste(img2, (int(W/3.0),0))
new.paste(img3, (int(2.0*W/3.0),0))

new.paste(img4, (0,int(H/3)))
new.paste(img5, (int(W/3.0),int(H/3)))
new.paste(img6, (int(2.0*W/3.0),int(H/3)))

new.paste(img7, (0,int(2.0*H/3)))
new.paste(img8, (int(W/3.0),int(2.0*H/3)))
new.paste(img9, (int(2.0*W/3.0),int(2.0*H/3)))

figname = "proj_dmdens_512.png"
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

