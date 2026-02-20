from PIL import Image, ImageEnhance

W  = 1280
H  = 1920

new = Image.new("RGBA", (int(W),int(H)))

img1 = Image.open("pp_DM_particle_mass_snsink_256_1.png")
img2 = Image.open("pp_DM_particle_mass_snagn_256_1.png")
img3 = Image.open("proj_nHI_snsink.png")
img4 = Image.open("proj_nHI_snagn.png")
img5 = Image.open("pp_st_particle_mass_snsink_256_1.png")
img6 = Image.open("pp_st_particle_mass_snagn_256_1.png")
img7 = Image.open("tn_mH_snsink.png")
img8 = Image.open("tn_mH_snagn.png")

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


img1 = img1.resize((int(W/2),int(H/4)))
img2 = img2.resize((int(W/2),int(H/4)))

img3 = img3.resize((int(W/2),int(H/4)))
img4 = img4.resize((int(W/2),int(H/4)))

img5 = img5.resize((int(W/2),int(H/4)))
img6 = img6.resize((int(W/2),int(H/4)))

img7 = img7.resize((int(W/2),int(H/4)))
img8 = img8.resize((int(W/2),int(H/4)))


new.paste(img1, (0,0))
new.paste(img2, (int(W/2.0),0))

new.paste(img3, (0,int(H/4)))
new.paste(img4, (int(W/2.0),int(H/4)))

new.paste(img5, (0,int(2.0*H/4)))
new.paste(img6, (int(W/2.0),int(2.0*H/4)))

new.paste(img7, (0,int(3.0*H/4)))
new.paste(img8, (int(W/2.0),int(3.0*H/4)))


figname = "tesplots_2.png"
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

