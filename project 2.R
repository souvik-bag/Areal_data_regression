### Load the required libraries
suppressMessages(library(sf))
suppressMessages(library(tidyverse)) 
suppressMessages(library(magrittr)) 
suppressMessages(library(dplyr)) 
suppressMessages(library(sf)) 
suppressMessages(library(spdep))
suppressMessages(library(ggplot2))
suppressMessages(library(spatialreg))

# load the data set
MO_county_sf <- st_read("./MO_2022_County_Boundaries/MO_2022_County_Boundaries_shp.shp",
                        stringsAsFactors = FALSE,
                        as_tibble = TRUE)
# Load the aditional  data downloaded from internet
Mo_income <- read.csv("MO_income.csv",header = F)
Mo_age_senior <- read.csv("Age_senior.csv",header = T)
Mo_unempl <- read.csv("C:\\Users\\Admin\\Downloads\\Missouri_Unemployment_Rates_For_Various_Regions.csv", header = T)

# Pre processing 
Mo_unempl %<>%
  filter(Year %in% c(2013), Area.Type %in% c(4), Period %in% c(2))
Mo_unempl %<>% select(Area.Name,Unemployment.Rate,RegionGeometry) %<>% rename(geometry = RegionGeometry) 
Mo_unempl =  Mo_unempl [order(Mo_unempl $Area.Name),]
Mo_sf_order <- Mo_sf[order(Mo_sf$COUNTYN),]
Mo_sf_order$Income = (Mo_income$V1)
Mo_sf_order$senior = Mo_age_senior$Age
Mo_sf_order$Unempl = Mo_unempl$Unemployment.Rate


# Visualize Vaccine rate
ggplot(Mo_sf_order) + geom_sf(aes(fill=(vac_rate))) +
  theme_bw() + labs(fill="vaccination rate \n (in '000s)") +
  scale_fill_continuous(low="khaki", high = "firebrick") + labs( x = ("Longitude"), y = ('Latitude'))

# Visualize test rate
ggplot(Mo_sf_order) + geom_sf(aes(fill=(test_rate))) +
  theme_bw() + labs(fill="Testing rate \n (in '000s)") +
  scale_fill_continuous(low="khaki", high = "firebrick") + labs( x = ("Longitude"), y = ('Latitude'))

# Visualize positive rate
ggplot(Mo_sf_order) + geom_sf(aes(fill=(positive_rate))) +
  theme_bw() + labs(fill="positive test rate \n (in '000s)") +
  scale_fill_continuous(low="khaki", high = "firebrick") + labs( x = ("Longitude"), y = ('Latitude'))




# Convert to data frame

Mo_df <- as.data.frame(MO_county_sf)

Mo_df %<>%  mutate(test_rate = (PCR_Tst + Antgn_T)/Popultn * 1000, 
                  positive_rate = (PCR_Pos + Antgn_P)/Popultn * 1000,
                  vac_rate = (FemalVc + MaleVac + UnknwSV)/Popultn * 1000)
Mo_sf <- st_as_sf(Mo_df)
## 
# Calculate the proximity matrix # Obtain the coordinates of the areal centroids 
centroids <- st_coordinates(st_centroid(Mo_sf$geometry))
# Label each neighborhood with their indices (POLYID) 
ggplot(Mo_sf) + geom_sf() + geom_label(aes(x=centroids[,1], y=centroids[,2], label = COUNTYG ), size=2) + 
  theme_bw() + xlab("Longitude") + ylab("Latitude")

# Neighbours using distance
Mo_dist1 <- dnearneigh(centroids, d1 = 0, d2 = 70000, row.names = Mo_sf$COUNTYG)
# Calculating the list of neighbours 
# nb_q <- poly2nb(Mo_sf, queen = TRUE) 
# First 2 nearest neighbors
# col.knn <- knearneigh(centroids, k=2)
# col.knn.nb <- knn2nb(col.knn, row.names= Mo_sf$COUNTYG)

# Calculating the proximity matrix 
W = nb2listw(Mo_dist1, style="W", zero.policy = TRUE)
W2 = nb2listw(Mo_dist1, style="B", zero.policy = TRUE)
W3 = nb2listw(col.knn.nb, style="W", zero.policy = TRUE)

# Create W matrix for CAR
W.sym=similar.listw(W)
W_mat=as.matrix(as_dgRMatrix_listw(W))
W.sym_mat=as.matrix(as_dgRMatrix_listw(W.sym))




plot(st_geometry(Mo_sf), border = "grey60", reset = FALSE)
plot(Mo_dist1, coords = centroids, add=T, col = "red")


# Moran's I statistic for test rate
moran.test(Mo_sf$test_rate, listw = W, randomisation=FALSE, alternative="greater", zero.policy = TRUE)
# Moran's I statistic for positive rate
moran.test(Mo_sf$positive_rate, listw = W3, randomisation=FALSE, alternative="greater", zero.policy = TRUE)
# Moran's I statistic for vaccine rate
m1 = moran.test(Mo_sf$vac_rate, listw = W, randomisation=FALSE, alternative="greater", zero.policy = TRUE)
moran.test(res, listw = W, randomisation=FALSE, alternative="greater", zero.policy = TRUE)
lm.morantest(model,listw = W3, alternative = "greater", zero.policy = TRUE)

## Fit a linear model
model <- lm(positive_rate ~ hoval  + test_rate + vac_rate + Income + Popultn + senior, data = Mo_sf_order)
summary(model)
# print the summary of linear model in latex
print(xtable(summary(model), type = "latex"))

# Plot the residual
library(viridis)
library(ggplot2)
colpalette <- viridis_pal(option = "magma")
colpal_rev <- rev(colpalette(20)) # sequential scheme: light -> low values, dark -> high values
# calculate the residual
res=resid(model)
# Standardize them
Mo_sf_order$res=(res-min(res))/diff(range(res))
# plot
ggplot(Mo_sf_order) + geom_sf(aes(fill= res)) +
  theme_bw() + ggtitle("Standardized residuals") +
  scale_fill_gradientn(name = "Residuals", colours = colpal_rev, na.value = "black", limits = c(0,1))


# Geary's C statistic
geary.test(residuals(model), listw = W, randomisation=FALSE,
alternative="greater", zero.policy = TRUE)

# print the result of linear model in latex
print(xtable(summary(model), type = "latex"))

# Residual analysis 
par(mfrow=c(1,2))
plot(fitted(model), res, ylab = "Residuals", xlab = "Y predicted")
abline(0,0, col ='red')
qqPlot( resid(model), ylab = 'Residuals')
par(mfrow=c(1,1))


# Fit a CAR model


# See the image of weight matrix
image(W.sym_mat)

# Fit CAR
model.car=spautolm(positive_rate ~ hoval  + test_rate + vac_rate + Income + Popultn + senior, data = Mo_sf_order, W.sym, family="CAR")
# Summary of the model
summary(model.car)


# extract and format results for LaTeX
summary_table <- summary(model.car)
coef_table <- xtable(summary_table$Coef)


## Q3 Plot the Fitted data

model.fitted <-(model.car$fit)$fitted.values

Mo_sf_order$fitted <- model.fitted


# Plot


# Visualize one variable
ggplot(Mo_sf_order) + geom_sf(aes(fill=( model.fitted))) +
  theme_bw() + labs(fill="Fitted positive test rate \n (in '000s)") +
  scale_fill_continuous(low="khaki", high = "firebrick") + labs( x = ("Longitude"), y = ('Latitude'))






