#initialize the libraries
#TO ANAND: This should be the same format as the Gcross code set up for snakemake. For now, I have some dummy code for testing.

library(spatstat)
library(sp)
library(reshape2)
library(ggplot2)
xmin = 0
xmax = 1392
ymin = 0
ymax = 1040

R = 40 #radius in pixels-convert from micron before applying to code)
##oad the image and section out the different phenotypes of interest

fn1='2_cell_seg_data.txt' #sample data
fn1=read.table(fn1,header=TRUE,sep="\t")

# Identify cell types
fn1$CellOfInterest = ""
fn1[which(fn1$Epithelial == "pos"), "CellOfInterest"]= "Tumor"
fn1[which(fn1$CTL == "pos"), "CellOfInterest"]= "CD8"
fn1[which(fn1$CD4 == "pos"), "CellOfInterest"]= "CD4"

dataPoints=data.frame(Cell.X.Position=(fn1$Cell.X.Position),
                      Cell.Y.Position=(fn1$Cell.Y.Position), 
                      CellOfInterest= fn1$CellOfInterest)

dataPoints <- dataPoints[dataPoints$CellOfInterest!="", ]
#due to the way the code is setup, we have to interrogate each pair of cells separately. It would mean having some additional steps at the end
#to concatenate all the results for cell neighbours for an individual instance of the reference cell.

Cell1 = "Tumor"
List_Cell2 = c("CD4","CD8")

#the convex hull of the current set of coordinates is determined
cw <- convexhull.xy(dataPoints$Cell.X.Position,dataPoints$Cell.Y.Position)
ow <- cw$bdry[[1]]
#an observation window is created 
ww <- owin(poly = ow)

#convert the x,y coordinates into the form of spatial point patterns
pp <- as.ppp(cbind(dataPoints$Cell.X.Position,dataPoints$Cell.Y.Position), W = owin(c(xmin,xmax),c(ymin,ymax))) # As used in NatComm paper

#creates a list with the phenotype and the spatial points, along with their frequencies which can
#be used by the spatstat package
pp <- pp %mark% factor(dataPoints$CellOfInterest)
#the bounding box for polygon window is added as a factor to the list
pp$window <- ww

#Split each cell into a list
Y <- split(pp)

CrossPairCoordData=list()
#Run crosspairs on each combo
for (i in 1:length(List_Cell2)){
  CrossPairCoordData[[paste0(Cell1,"vs",List_Cell2[i])]] <- crosspairs(Y$Tumor, Y[[paste0(List_Cell2[i])]], R)
CrossPairCount <- crosspaircounts(Y$Tumor, Y[[paste0(List_Cell2[i])]], R)
}

#The next step would be to consolidate all the results from the crosspairs in a reference cell-by-cell manner, 
#in such a way that for each ref cell in an image, we have a count of how many cells of interest are present at the given distance
#from it. We can make a table with Cell X,Y coords and the other columns containing counts of neighbours of each type.


