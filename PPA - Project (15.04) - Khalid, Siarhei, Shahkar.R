###############
# Call spatstat
###############

library(spatstat)
library(maptools)

#################
# Inspecting data
#################

X  # information on the point pattern. Point patterns in spatstat are objects of class "ppp"

# Note : you can create your own point pattern from a text file with XY coordinates with the command ppp, or can be converted from GIS formats with readShapePoints

# To extract or manipulate the data in a point pattern object:
x <- Pop
plot(x)
x
# The population data is a Marked planar point pattern with 296 points
# the Average intensity is 4.775791e-07 points per square unit
popdata <- Pop

is.multitype(Pop)
summary(popdata)
plot(popdata)
popdata
class(popdata)
npoints(popdata)
chin <- cut(popdata$marks, breaks = 4)
chin <- as.factor(chin)
marks(popdata) <- chin
marks(popdata) <- factor(chin)
is.multitype(popdata)
plot(popdata)

M <- owin(c(514117.1, 536080.4), c(1827049,1857177))
popdata1 <- popdata[M]
plot(popdata1)
# Intensity
plot(density(popdata1,1000))   # local spatial variations in density
contour(density(popdata1,1000), axes = F)   # to have a contour plot

###########################
# Exploratory data analysis
###########################

# Quadrat counting

Qu <- quadratcount(popdata1, nx = 3, ny = 3)
Qu
plot(popdata1)
plot(Qu, add=TRUE, cex = 2)
# since the number of points in each quadrant is not the same, we dont have 
#homogeneous intensity, hence, inhomogeneous intensity. 

# Ripley function

K <- Kest(popdata1)
K
plot(K)

# envelopes of the K function
E <- envelope(popdata1,Kest,nsim=200)
plot(E)
# Here we compute some confidence interval unsing boostrapping. 


########
# Models
########
#THis is a poisson model. Process fitting 
fit <- ppm(popdata1,~1, Strauss(0.1))   # Fit a Strauss process to the data
fit
plot(simulate(fit))  # Plot the simulated data with the process
#Interaction distance:	0.1
#Fitted interaction parameter gamma:	 1.12e-05
plot(envelope(fit,Kest,nsim=200)) # Make a goodness of fit test for this fitted model



##########################
# Multitype point patterns (categorical marks)
##########################

levels(chin) # to see the types of the categorical mark
# There are 4 levels of population based on the popualtion value. 
table(chin) # Contingency table stating how many population we have in each categpry

plot((popdata1))

plot(split(popdata1))   # Split the point pattern depending on type

################### INTENSITY #################

# To have an estimate of lambda in the case of homogenous intensity just do a summary of the pp
lamb <- summary(popdata1)$intensity   # Extract the intensity value
lamb
#the intensity is 0.0000004322133

# Quadrat counting for inhomogenous intensity

Q <- quadratcount(popdata1, nx=3,ny=3)
Q
plot(Q,add=TRUE, cex = 1.5, col="red")
# There is inhomogeneity becuase the number of points per quadrant is not equal
#There are more concentration of points in some palces tahn the other. 

#Kernel Smoothing
den <- density(popdata1,kernel = "epanechnikov",sigma=300)  # The resulting object is a pixel image. This class has methods for print, summary, plot, contour and persp
plot(den)
plot(popdata,add=TRUE)
#more global and visual represenatation. but there is not test involved. 
#it is a fancy way of doing quadrat counting, a nonparametric way. 
#There is more concentration of points at the middle than other place. Hence, not homo
summary(den)
#the density is a real-valued pixel image of size 128 x 128 pixel. 
contour(den)
# 3D plot 
persp(den)

mean(den$v, na.rm = TRUE)
# the average intensity is 4.759229e-07
median(den$v, na.rm = TRUE)
#the median value is 1.257512e-22
box <- !is.na(den$v)
boxplot(box) # First convert the v matrix into a vector

H <- quadrat.test(popdata1, nx = 3, ny = 3 )
#This is a basic chi-sqauree test. We compare the observed realisation with the expected value (H0). 
#H0 = complete random distribution. We reject the H0 of complete random distribution since the p-value is less.
#there could be some other information here. It could actually mean inhomogeneity or interdependence. or more. 
# The chi-sqaure test here is 218.5
# However, the test has be criticsed for so many reasons. 
plot(popdata1)
plot(H, add=TRUE, cex = 0.8, col="red")


# Particular case : Distance maps : intensity depending on distance

Z <- distmap(popdata1)
plot(popdata1,lwd=2,main="")
contour(Z,add=TRUE)
plot(rhohat(popdata1,Z),xlab="Distance to nearest fault") # Estimate of the intensity as a function of distance to nearest fault

# Fitting models
exp(-16.04064)
popdata3 <- rescale(popdata1, 1000)
popdata4 <- rescale(popdata3, 1000)
ppm(popdata3,~1)               # Homogeneous Poisson model. 1 means we have constant intensity. The Intensity here is 0.0000001080534
ppm(popdata3,~x+y)             # Inhomogeneous Poisson model : intensity is log-linear in the cartesian coordinates
#The result shows that the log intensity varies linear based on x and y as the result is significant. 
ppm(popdata4,~polynom(x,y,2))  # Inhomogeneous Poisson model : intensity is log-quadratic in the cartesian coordinates
#The result here are also significant. 
fit <- ppm(popdata4,~polynom(x,y,2))
plot(effectfun(fit,"x", y=0.5)) #This shows how intensity vary based on Latitude. It remains in the same level throughout.
plot(effectfun(fit,"y", x=0.5)) # This shows how intensity varies based on longitude. It was first low, 
#It went high afterwards to above 0.8 and it went down again to around 0. 

Z <- distmap(popdata)
fit <- ppm(popdata,~Z)               # Inhomogeneous Poisson model : intensity depends on a covariate
fit
plot(effectfun(fit)) # Plots fitted curve or lambda against Z

#grad <- popdata$marks
#plot(grad)
#ppm(popdata, ~slope, covariates = list(slope = grad))  # Fits a Poisson model with intensity that is loglinear function of slope l(u)=exp(b0+b1.Z(u))
# Here the estimated intensity is about exp(-5.391) = 0.004559 trees per square metre or 45.59 trees per hectare and increases by a factor of exp(5.026) = 152.4 if the slopes increases by 1. 

# Predict fitted values

fit <- ppm(popdata1,~x+y)  
lam <- predict(fit)
plot(lam)
#With a fitted model where the intensity depends on both coordinates, it was 
#realised that the intensity depends on x (latitude) and y (longitude) as it changes 
# as we move from south to north and west to east. 

# It is possible to retrieve elements from ppm objects

fit <- ppm(popdata3, ~x + y)
fit
coef(fit)
vcov(fit)
sqrt(diag(vcov(fit)))
round(vcov(fit, what = "corr"), 2)

SE <- predict(fit, se = TRUE) # standard error of the fitted intensity at each location
plot(SE, main = "standard error of fitted intensity")
# This plot shows us the part that we make most error in our fitting. 
# The blue part shows little or no variation in what is observed and what was fitted. 
# However, the yellow part show bigger variation (errors). 

# Likelihood ratio test - analysis of deviance / always use nested models
# H0 = Homogeneous poisson process. 
fit0 <- ppm(popdata3,~1) #-H0
fit1 <- ppm(popdata3,~polynom(x,y,2)) #-H1
#Analysis of Variance test which is formally the same as likelihood ratio test
anova(fit0,fit1,test="Chi")
#While the model 1 is a homogeneous poisson model, The model two is a polynomial
#form of degree 2 of both coordinate x and y. 
# The likelihood ratio test result here is 189.27 and the p-value is 0.00000000000000022 
# less than 5% hence the second model is significant, the inhomogeneous result is significant. 
# The Poisson process is inhomogeneous. 

step(fit)

# Goodness of fit
# We also try the non parametric quarat test.(This test may not be necessary, since 
#the last result is inhomogeneous)
fit <- ppm(Pop,~x)
M <- quadrat.test(fit,nx=3,ny=3)
M  # the gof test rejects the fitted model: what is the departure of the model ? Inspect residual counts
plot(bei,pch=".")
plot(M,add= TRUE, cex=1.5,col="red")

# validation using residuals

fit <- ppm(popdata,~x+y)  
plot(predict(fit))
plot(bei,add=TRUE,pch="+")          # Fitted versus observed

fitx <- ppm(popdata,~x)  
diagnose.ppm(fitx,which="smooth")   # smoothed residuals

grad <- bei.extra$grad
lurking(fitx,grad,type="raw")       # Lurking variable plot to see if fitted intensity missing a dependence to another covariate

miplot(popdata1) #This is a morisita plot. This shows cluster since they are high values. subjective of spatial clustering
#This could be an indication of interdependence. 
fryplot(popdata1, axes=TRUE) #This fryplot is completely unreadable. it suggests nothing. But could mean clustering since
#there is no distance ayt the middle. 

# Distance methods
#we do analysis of the location of populations categories here.
emp <- distmap(popdata1)
plot(emp,main="Empty space distance") #This plot shows the minimum distance around each points. 
#for most point the distance is quite small. Also indicating clustering. 
plot(popdata1, add= TRUE)

d <- nndist(popdata1)
plot(popdata1 %mark% nndist(popdata1), markscale=1) #This also gives the minimum distance for 
#each point to another. There are more of small distances. #

d <- pairdist(popdata1)
d3 <- pairdist(popdata1, periodic=TRUE)
d2 <- pairdist(popdata1, squared=TRUE)
plot(d3) 
plot(d2)

# K and L functions; envelopes and bootstrat for CI

Gc <- Kest(popdata1)
Gc
par(pty = "s", cex=0.6)
plot(Gc)
# The blue line is h0, and the black line is estimated value, red and green are evelops
# This signifies that it is not an homogeneous poision process and some form for clustering
# since the e,perica K function is above the h0 line. 

#Lest is the transformatiom of K. The square of K divided by pie. 

L <- Lest(popdata1)
plot(L, main = "L function") # This give us better visualization of what happened in Kest. 

E <- envelope(popdata1, Kest, nsim = 200, rank = 1)
E
plot(E, main = "pointwise envelope")  
# Since all the emperical k function are not completely inside the grey area, it is not 
# comepletely random. There are some clustering under the homogenous poission assumptions. 

E <- envelope(popdata1, Lest, nsim = 200, rank = 1, global = TRUE)  # This stabilizes the variance
E
plot(E, main = "global envelopes of L(r")
# The same as above is observed here. 

v <- varblock(popdata1, Kest, nx=3, ny=3) # This is for the K function The confidence interval 
# for the estimated k function. 
plot(v) # 95% confidence interval (shaded) for teh true value of K function
# obtained using boostrap


################### Marked pp #################

# Example of a multitype pp

plot(split(popdata1))
plot(density(split(popdata1))) #This shows that the 4 different category behave differently. 
# There is indeed segregation based on the marks. 
table(marks(popdata1))

M <- marktable(popdata1,R = 0.1)  # Contingency table of the marks of all points within a given radio of each data point
plot(Gcross(popdata1)) # Distribution function of the distance from a point ot type i to the nearest point of type j


# Segregration

lansP <- relrisk(popdata1) #we compute the relative risk ratio here. 
plot(lansP)
#we compared what we have to what we'll have if we dont have segregation. H0 is the absense of segregation. 
#The relative risk map also suggest that we have segregation since we have differences in 
#the colors i.e not monocolor. There is differences bewtween the observed intensity and 
#the intensity under H0, hence we don't have a stationary multitype poission process. 

#Interaction between the different categories
plot(Gcross(popdata1), cex = 0.2)

#Modelling Multitype inhomogenous dependent point pattern model
fit0 <- ppm(popdata1,~marks + polynom(x,y,3))
fit1 <- ppm(popdata1,~marks * polynom(x,y,3))
anova(fit0,fit1,test="Chi")
plot(predict(fit0))
plot(predict(fit1))

segregation.test(popdata1,nsim=19)


# Inhomogeneous K multitype K function

fit1 <- ppm(popdata1,~marks * polynom(x,y,3))
lamb <- predict(fit1)
plot(lamb)
plot(Kcross.inhom(popdata1), cex = 0.2)
 
