#Colour palettes script 

#colour brewer palettes
require(RColorBrewer)
display.brewer.all()

#Baseline colours pastel
wild = "#C7E9C0"
wild_like = "#C6DBEF"
lab_like = "#DECBE4"
CD1 = "#FBB4AE"
MF1 = "#FED9A6"
#Baseline colours prominent
wild = "#74C476"
wild_like = "#6BAED6"
lab_like = "#BC80BD"
CD1 = "#FB8072"
MF1 = "#FDB462"

Reds<-brewer.pal(9, c("Reds"))
#"#FFF5F0" "#FEE0D2" "#FCBBA1" "#FC9272" "#FB6A4A" "#EF3B2C"
#"#CB181D" "#A50F15" "#67000D"

OrgRd <-brewer.pal(9,c("OrRd"))
#"#FFF7EC" "#FEE8C8" "#FDD49E" "#FDBB84" "#FC8D59" "#EF6548"
#"#D7301F" "#B30000" "#7F0000"

Oranges <-brewer.pal(9, c("Oranges"))
#"#FFF5EB" "#FEE6CE" "#FDD0A2" "#FDAE6B" "#FD8D3C" "#F16913"
#"#D94801" "#A63603" "#7F2704"

YelOrgRb<-brewer.pal(9, c("YlOrRd"))
# "#FFFFCC" "#FFEDA0" "#FED976" "#FEB24C" "#FD8D3C" "#FC4E2A"
# "#E31A1C" "#BD0026" "#800026"

YelOrgBrn <-brewer.pal(9, c("YlOrBr"))
# "#FFFFE5" "#FFF7BC" "#FEE391" "#FEC44F" "#FE9929" "#EC7014"
# "#CC4C02" "#993404" "#662506"

YelGrn<-brewer.pal(9, c("YlGn"))
# "#FFFFE5" "#F7FCB9" "#D9F0A3" "#ADDD8E" "#78C679" "#41AB5D"
# "#238443" "#006837" "#004529"

Greens <- brewer.pal(9, c("Greens"))
# "#F7FCF5" "#E5F5E0" "#C7E9C0" "#A1D99B" "#74C476" "#41AB5D"
# "#238B45" "#006D2C" "#00441B"

GrnBlu <- brewer.pal(9, c("BuGn"))
# "#F7FCFD" "#E5F5F9" "#CCECE6" "#99D8C9" "#66C2A4" "#41AE76"
# "#238B45" "#006D2C" "#00441B"

YelGrnBlu<-brewer.pal(9, c("YlGnBu"))
# "#FFFFD9" "#EDF8B1" "#C7E9B4" "#7FCDBB" "#41B6C4" "#1D91C0"
# "#225EA8" "#253494" "#081D58"

BluGrn <- brewer.pal(9, c("GnBu"))
# "#F7FCF0" "#E0F3DB" "#CCEBC5" "#A8DDB5" "#7BCCC4" "#4EB3D3"
# "#2B8CBE" "#0868AC" "#084081"

PurBluGrn <-brewer.pal(9, c("PuBuGn"))
# "#FFF7FB" "#ECE2F0" "#D0D1E6" "#A6BDDB" "#67A9CF" "#3690C0"
# "#02818A" "#016C59" "#014636"

Blues <- brewer.pal(9, c("Blues"))
# "#F7FBFF" "#DEEBF7" "#C6DBEF" "#9ECAE1" "#6BAED6" "#4292C6"
# "#2171B5" "#08519C" "#08306B"

PurBlu <-brewer.pal(9, c("PuBu"))
# "#FFF7FB" "#ECE7F2" "#D0D1E6" "#A6BDDB" "#74A9CF" "#3690C0"
# "#0570B0" "#045A8D" "#023858"

BluPur <-brewer.pal(9, c("BuPu"))
# "#F7FCFD" "#E0ECF4" "#BFD3E6" "#9EBCDA" "#8C96C6" "#8C6BB1"
# "#88419D" "#810F7C" "#4D004B"

Purples <-brewer.pal(9, c("Purples"))
# "#FCFBFD" "#EFEDF5" "#DADAEB" "#BCBDDC" "#9E9AC8" "#807DBA"
# "#6A51A3" "#54278F" "#3F007D"

PurRd <-brewer.pal(9, c("PuRd"))
# "#F7F4F9" "#E7E1EF" "#D4B9DA" "#C994C7" "#DF65B0" "#E7298A"
# "#CE1256" "#980043" "#67001F"

RdPur<-brewer.pal(9, c("RdPu"))
# "#FFF7F3" "#FDE0DD" "#FCC5C0" "#FA9FB5" "#F768A1" "#DD3497"
# "#AE017E" "#7A0177" "#49006A"


Greys <- brewer.pal(9, c("Greys"))
# "#FFFFFF" "#F0F0F0" "#D9D9D9" "#BDBDBD" "#969696" "#737373"
# "#525252" "#252525" "#000000"

#68 colours for Genus - currently 50 -add 18 more 
gen_cols<-c("#FEE0D2","#FCBBA1","#FB8072","#EF3B2C","#CB181D","#A50F15","#67000D","#B30000","#FC4E2A","#FC8D59",
            "#FDAE6B","#FDD0A2","#FDB462","#FD8D3C","#F16913","#D94801","#A63603","#FE9929","#FEC44F","#FED976",
            "#FFED6F","#FFEDA0","#F7FCB9","#D9F0A3","#ADDD8E","#B3DE69","#74C476","#41AB5D","#41AE76","#006D2C",
            "#00441B","#014636","#016C59","#66C2A4","#8DD3C7","#7BCCC4","#4EB3D3","#4292C6","#045A8D","#2171B5",
            "#08519C","#253494","#08306B","#0570B0","#3690C0","#6BAED6","#80B1D3","#9ECAE1","#A6BDDB","#D0D1E6",
            "#BEBADA","#8C6BB1","#88419D","#54278F","#3F007D","#4D004B","#BC80BD","#C994C7","#D4B9DA","#FDE0DD",
            "#FCC5C0","#FA9FB5","#F768A1","#DD3497","#AE017E","#CE1256","#980043","#800026","#67001F")
 

#40 colours for Family 
fam_cols<-c("#FCBBA1","#FB8072","#EF3B2C","#A50F15","#67000D","#FC9272","#FDD49E","#FDB462",
            "#FD8D3C","#F16913","#D94801","#FEC44F","#FED976","#FFED6F","#D9F0A3","#B3DE69",
            "#74C476","#41AE76","#006D2C","#00441B","#014636","#016C59","#66C2A4","#8DD3C7",
            "#4EB3D3","#2171B5","#08306B","#80B1D3","#BEBADA","#8C6BB1","#54278F","#3F007D",
            "#BC80BD","#D4B9DA","#FCC5C0","#FA9FB5","#DD3497","#AE017E","#980043","#67001F")

#15 cols for hpoly inf expt 
gen_cols <-c("#FCBBA1","#FB8072","#FDB462","#FDD49E","#FFED6F","#D9F0A3","#B3DE69","#41AE76", 
             "#8DD3C7","#3690C0","#80B1D3","#BEBADA","#BC80BD","#FA9FB5","#F768A1")


#SETS 
pastel1 <- brewer.pal(9, c("Pastel1"))
display.brewer.pal(9, "Pastel1")
pastel1
# "#FBB4AE" "#B3CDE3" "#CCEBC5" "#DECBE4" "#FED9A6" "#FFFFCC" "#E5D8BD" "#FDDAEC" "#F2F2F2"

pastel2 <- brewer.pal(8, "Pastel2")
display.brewer.pal(8, "Pastel2")
pastel2
# "#B3E2CD" "#FDCDAC" "#CBD5E8" "#F4CAE4" "#E6F5C9" "#FFF2AE" "#F1E2CC" "#CCCCCC"

set3 <- brewer.pal(12, "Set3")
display.brewer.pal(12, "Set3")
set3
# "#8DD3C7" "#FFFFB3" "#BEBADA" "#FB8072" "#80B1D3" "#FDB462" "#B3DE69" "#FCCDE5" "#D9D9D9" "#BC80BD" "#CCEBC5" "#FFED6F"

paired <-brewer.pal(12, "Paired")
display.brewer.pal(12, "Paired")
paired
# "#A6CEE3" "#1F78B4" "#B2DF8A" "#33A02C" "#FB9A99" "#E31A1C" "#FDBF6F" "#FF7F00" "#CAB2D6"
# "#6A3D9A" "#FFFF99" "#B15928"

purples <-brewer.pal(9, "Purples")
display.brewer.pal(9, "Purples")
purples
# "#FCFBFD" "#EFEDF5" "#DADAEB" "#BCBDDC" "#9E9AC8" "#807DBA" "#6A51A3" "#54278F" "#3F007D"

greens <-brewer.pal(9, "Greens")
display.brewer.pal(9, "Greens")
greens
# "#F7FCF5" "#E5F5E0" "#C7E9C0" "#A1D99B" "#74C476" "#41AB5D" "#238B45" "#006D2C" "#00441B"


lw_colours <- c(set3, purples, greens)
lw_colours
# "#8DD3C7" "#FFFFB3" "#BEBADA" "#FB8072" "#80B1D3" "#FDB462" "#B3DE69" "#FCCDE5" "#D9D9D9"
# "#BC80BD" "#CCEBC5" "#FFED6F" "#FCFBFD" "#EFEDF5" "#DADAEB" "#BCBDDC" "#9E9AC8" "#807DBA"
# "#6A51A3" "#54278F" "#3F007D" "#F7FCF5" "#E5F5E0" "#C7E9C0" "#A1D99B" "#74C476" "#41AB5D"
# "#238B45" "#006D2C" "#00441B"
list(lw_colours)

spectral <- RColorBrewer::brewer.pal(8, "Spectral")
spectral
# "#D53E4F" "#F46D43" "#FDAE61" "#FEE08B" "#E6F598" "#ABDDA4" "#66C2A5" "#3288BD"

ColoursRamp <- colorRampPalette(spectral)
ColoursExpand <- ColoursRamp(n=12)
ColoursExpand
# "#D53E4F" "#E85B47" "#F67E4B" "#FCA85E" "#FDC977" "#F9E38D" 
# "#EAF195" "#CBEA9D" "#A4DAA4" "#78C9A4" "#53ACAD" "#3288BD"
