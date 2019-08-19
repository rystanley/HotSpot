## Function to do a Hotspot analysis  

HotspotFun <- function(xyz,resolution=20,nearneighdist=50,poly=NULL,region=NULL,
                       rasterFun=function(x,...){length(unique(na.omit(x)))},
                       ap=0.1,facet=NULL,bathy=NULL,otherMPAs=FALSE,variable="",trans=NULL,
                       forceextent = NULL,buffer=100,extended=F,
                       cropregion=F,facetorder=NULL,ncol=2,linewidth=1.05){
  
  # xyzi - dataframe with columns ordred by Longitude, Latitude, Value, SampleSite ID
  #resolution - is minimum grid cell resolution in km (default is 20 or 20 km2)
  #rasterFun - function for the rasterize function -- see ?raster::rasterize for details 
  
  #https://stackoverflow.com/questions/48650274/spatial-efficient-way-of-finding-all-points-within-x-meters-of-a-point
  #https://www.rdocumentation.org/packages/spdep/versions/0.7-7/topics/localG
  #https://rpubs.com/adam_dennett/126356
  
  #load required libraries
  require(sp)
  require(rgdal)
  require(raster)
  require(shape)
  require(rgeos)
  require(raster)
  require(ggplot2)
  require(mapdata)
  require(marmap)
  require(maptools)
  require(spdep)
  require(dplyr)
  
  #Map data for ggplot -- regions can be used to trim the data to a descrete planning region. We can make the 'null' 
  #option either a no trim or an object from the workshace. In this example it will just try to force a download. 
  if(is.null(region)){MaritimeRegion <- readOGR("data/Maritimes Planning/MaritimesPlanningArea.shp")}
  if(!is.null(region)){MaritimeRegion <- readOGR(region)}

  ## get extent of all data
  tempextent <- xyz
  colnames(tempextent) <- c("Long","Lat","Z")
  coordinates(tempextent) <- c("Long", "Lat")
  proj4string(tempextent) <- CRS("+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0") 
  if(is.null(forceextent)){latlong_extent <- extent(tempextent)} else {latlong_extent <- forceextent}
  tempextent <- spTransform(tempextent,CRS("+proj=utm +zone=20 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"))
  
  ## create raster grid according to the resolution parameter
  if(is.null(forceextent)){grid <- raster(extent(tempextent))}
  
  if(!is.null(forceextent)){
    temppoly <- as(forceextent, "SpatialPolygons")
    proj4string(temppoly) <- CRS("+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0") 
    temppoly <- spTransform(temppoly,CRS("+proj=utm +zone=20 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"))
    newextent <- extent(temppoly)
    
    if(!is.null(buffer)){newextent[1] <- newextent[1]-buffer*1000
    newextent[2] <- newextent[2]+buffer*1000
    newextent[3] <- newextent[3]-buffer*1000
    newextent[4] <- newextent[4]+buffer*1000}
    
    grid <- raster(newextent)
  }
  
  res(grid) <- resolution*1000 
  proj4string(grid)<-proj4string(tempextent)
  
  if(cropregion){
    tempregion <- spTransform(MaritimeRegion,CRS("+proj=utm +zone=20 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"))
    regionextent <- extent(tempregion)
    if(!is.null(buffer)){
      regionextent[1] <- regionextent[1]-buffer*1000
      regionextent[2] <- regionextent[2]+buffer*1000
      regionextent[3] <- regionextent[3]-buffer*1000
      regionextent[4] <- regionextent[4]+buffer*1000}
    grid <- crop(grid,regionextent)
    xyz2 <- xyz
    colnames(xyz2) <- c("Long","Lat","Z")
    coordinates(xyz2) <- c("Long", "Lat")
    proj4string(xyz2) <- CRS("+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0") # must be entered as long and lat in decimal degrees
    xyz2 <- spTransform(xyz2,CRS("+proj=utm +zone=20 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"))
    #make a box for the regional extent
    p <- as(regionextent,'SpatialPolygons')
    proj4string(p) <- proj4string(tempregion)
    xyzindex <- complete.cases(over(xyz2,p))
    xyz <- xyz[xyzindex,]
    if(!is.null(facet)){facet <- facet[xyzindex]}
  }
  
  
  #rasterize - no facet
  if(is.null(facet)){
    #Set up data -- Put into UTM coordinate system
    colnames(xyz) <- c("Long","Lat","Z")
    coordinates(xyz) <- c("Long", "Lat")
    proj4string(xyz) <- CRS("+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0") # must be entered as long and lat in decimal degrees
    xyz <- spTransform(xyz,CRS("+proj=utm +zone=20 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"))
    ras <- raster::rasterize(coordinates(xyz), grid, field = xyz$Z,fun=rasterFun)
    ras <- projectRaster(ras,crs="+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0")
    md <- as.data.frame(rasterToPoints(ras))
    colnames(md) <- c("Longitude", "Latitude", "MAP")
    
    ## Hotspot analysis -
    writeLines("Working on the Neighbourhood contiguity by distance matrix. This can take a some time depending on resolution and the nearneighdist paramter.")
    
    xycoords <- md%>%dplyr::select(Longitude,Latitude)%>%dplyr::rename(x=Longitude,y=Latitude)
    coordinates(xycoords) <- ~x+y
    proj4string(xycoords) <- CRS("+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0")
    xycoords <- spTransform(xycoords,CRS("+proj=utm +zone=20 +datum=NAD83 +units=km +no_defs +ellps=GRS80 +towgs84=0,0,0"))
    nbdata <- rasterToPolygons(ras)%>%poly2nb(.,queen=T)
    #nbdata <- coordinates(xycoords)%>%dnearneigh(.,0,nearneighdist)
    locG <- localG(md$MAP, nb2listw(nbdata, style="B",zero.policy=TRUE))
    md$val <- locG%>%as.matrix()
    
    writeLines("Finished hot spot analysis")
    
  }
  
  
  
  #rasterize - with facet grouping variable
  if(!is.null(facet)){
    carfacet <- as.character(facet)
    md <- NULL
    for (i in unique(carfacet)){
      txyz <- xyz[carfacet==i,]
      colnames(txyz) <- c("Long","Lat","Z")
      coordinates(txyz) <- c("Long", "Lat")
      proj4string(txyz) <- CRS("+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0") # must be entered as long and lat in decimal degrees
      txyz <- spTransform(txyz,CRS("+proj=utm +zone=20 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"))
      ras <- raster::rasterize(coordinates(txyz), grid, field = txyz$Z,fun=rasterFun)
      ras <- projectRaster(ras,crs="+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0")
      tmd <- as.data.frame(rasterToPoints(ras))
      colnames(tmd) <- c("Longitude", "Latitude", "MAP")
      tmd$yeargroup <- i
      
      writeLines(paste0("Working on facet ",i))
      
      nbdata <- rasterToPolygons(ras)%>%poly2nb(.,queen=T)
      #nbdata <- coordinates(xycoords)%>%dnearneigh(.,0,nearneighdist)
      locG <- localG(tmd$MAP, nb2listw(nbdata, style="B",zero.policy = TRUE))
      #globG <- globalG.test(round(md$MAP),nb2listw(nbdata, style="B"))
      tmd$val <- locG%>%as.matrix()
      
      md <- rbind(md,tmd)
    }
    
  }
  
  if(!is.null(facetorder) & !is.null(facet)){
    md$yeargroup <- factor(md$yeargroup,levels=facetorder)
  }
  
  #Set plotting limits
  rextent <- latlong_extent
  Long.lim  <-  c(rextent[1], rextent[2])
  Lat.lim <-  c(rextent[3], rextent[4])
  
  #change projections of the shape files to lat long 
  MaritimeRegion <- spTransform(MaritimeRegion,CRS("+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0"))
  
  states <- map_data("state")
  if(!extended){usa <- subset(states,region == "maine")} #just need maine
  
  if(extended){usa <- subset(states,region %in% c("maine","new hampshire",
                                                  "massachusetts","connecticut",
                                                  "rhode island","vermont"))}
  canada <- map_data("worldHires", "Canada")
  FortMar <- fortify(MaritimeRegion)
  
  #read in bathymetry data
  
  if(is.null(bathy)){
    if(!dir.exists("data")){dir.create("data/")} #create a data directory if it doesn't exist
    curdir <- getwd()
    setwd("data/")
    bathy <- getNOAA.bathy(Long.lim[1],Long.lim[2],Lat.lim[1],Lat.lim[2],res=1,keep=T)
    setwd(curdir)
  }
  
  #format bathy for ggplot
  bathy <- fortify(bathy)
  
  #Break MD into levels
  #breaklevels <- c(-2.58,-1.96,1.96,2.58)
  #breaklabs <- c("Very low","Low","High","Very high")
  
  md$breaks <- "Very low"
  md[md$val>-2.58 & md$val<=(-1.96) & !is.na(md$val),"breaks"] <- "Low"
  md[md$val>-1.96 & md$val<1.96 & !is.na(md$val),"breaks"] <- "Insignificant"
  md[md$val>=1.96 & md$val<2.58& !is.na(md$val),"breaks"] <- "High"
  md[md$val>=2.58& !is.na(md$val),"breaks"] <- "Very high"
  
  md$breaks <- factor(md$breaks,levels=c("Very low","Low","Insignificant","High","Very high"))
  
  
  #Make the plot
  p1 <- ggplot() +
    geom_tile(data=md,aes(x=Longitude,y=Latitude,fill=breaks))+
    geom_polygon(data = usa, 
                 aes(x=long, y = lat, group = group), 
                 fill = "white", 
                 color="black") +
    geom_polygon(data = canada, aes(x=long, y = lat, group = group), 
                 fill = "white", color="black") + 
    geom_polygon(data=FortMar,
                 aes(long, lat, group = group),fill=NA,col="black")+
    geom_contour(data=bathy,aes(x=x,y=y,z=z),breaks=c(-200),lwd=0.05,colour="grey20")+
    coord_fixed(xlim = Long.lim,  ylim = Lat.lim, ratio = 1.2)+
    theme_bw()+
    theme(legend.position = "bottom",
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = "white", colour = "black"),
          plot.background = element_rect(colour = "white"),
          strip.background = element_rect(colour = "black", fill = "white"))+
    labs(x=expression(paste("Longitude ",degree,"W",sep="")),
         y=expression(paste("Latitude ",degree,"N",sep="")))+
    scale_fill_manual(values=colorRampPalette(rev(c("red","yellow","springgreen","royalblue")))(5))
  
  
  if(!is.null(facet)){
    # to save room we will crop canada before we facet
    canadabox <- latlong_extent
    canadabox[1] <- canadabox[1]-1
    canadabox[4] <- canadabox[4]+4
    
    boundbox <- as(canadabox, 'SpatialPolygons')
    proj4string(boundbox) <- proj4string(Fundian)
    
    coordinates(canada) <- ~long + lat
    proj4string(canada) <- proj4string(Fundian)
    canada <- as(canada,'SpatialPointsDataFrame')
    ind_canada <- over(canada,boundbox)
    
    #reload and index
    canada <- map_data("worldHires", "Canada")
    canada <- canada[!is.na(ind_canada),]
    
    p1=p1+facet_wrap(~yeargroup,ncol=ncol)} # end if facet
  
  
  #return the plot back
  return(p1)
  
}