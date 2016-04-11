#check this thread to install Boom http://stackoverflow.com/questions/36034316/r-package-boom-fails-to-install-on-ubuntu-linux
install.packages((Boom) #dependency for BoomSpikeSlab
install.packages((BoomSpikeSlab) #dependency for bsts
install.packages((bsts) #package written by Scott 

library(bsts)
data(iclaims)
plot(initial.claims)

bsts_obj <- bsts()
#generate draws from posterior predictive distributions of a bsts object
prediction <- predict.bsts(bsts_obj, 
             horizon = 1, 
             )
#Plot the posterior predictive distribution from a bsts prediction object
plot.bsts.prediction(prediction, 
                     )
