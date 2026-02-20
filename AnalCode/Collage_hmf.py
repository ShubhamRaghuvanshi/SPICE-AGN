from PIL import Image, ImageEnhance

W  = 1920
H  = 920

new = Image.new("RGBA", (int(W),int(H)))
#img1 = Image.open("sink_mass_1.png")
#img2 = Image.open("sink_mass_2.png")
#img3 = Image.open("sink_mass_5.png")

#img1 = Image.open("sinkAGN_mass_1.png")
#img2 = Image.open("sinkAGN_mass_2.png")
#img3 = Image.open("sinkAGN_mass_5.png")

img1 = Image.open("hmf_cmp_z6.png")
img2 = Image.open("hmf_cmp_z6_512.png")

#img1 = Image.open("proj_nH_snsink_bigbox_z100.png")
#img2 = Image.open("proj_nH_snsink_bigbox_z9.png")

img1 = img1.resize((int(W/2),int(H)))
img2 = img2.resize((int(W/2),int(H)))

new.paste(img1, (0,0))
new.paste(img2, (int(W/2),0))

figname = "proj_nH_snsink_vs_z.png"
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