#check this thread to install Boom http://stackoverflow.com/questions/36034316/r-package-boom-fails-to-install-on-ubuntu-linux
install.packages((Boom) #dependency for BoomSpikeSlab
install.packages((BoomSpikeSlab) #dependency for bsts
install.packages((bsts) #package written by Scott 
#for error "`GLIBCXX_3.4.20' not found", check http://askubuntu.com/questions/575505/glibcxx-3-4-20-not-found-how-to-fix-this-error

library(BoomSpikeSlab)
library(bsts) #package written by Scott 
data(iclaims)
plot(initial.claims)

#Create a spike and slab prior for use with lm.spike
prior <- SpikeSlabPrior()

#spike and slab regression
lm.spike()

bsts_obj <- bsts()
#generate draws from posterior predictive distributions of a bsts object
prediction <- predict.bsts(bsts_obj, 
             horizon = 1, 
             )
#Plot the posterior predictive distribution from a bsts prediction object
plot.bsts.prediction(prediction, 
                     )




#REFERENCES:
#Scott, S., Varian, H. (2012) International Journal of Mathematical Modelling and Numerical Optimisation, vol. 5 (2014), pp. 4-23.
#Giannone, D., Reichlin, L., Small, D.. Journal of Monetary Economics 55 (2008) pp. 665-676
#Nakata, T., Tonetti, C. (2010) Kalman Filtere and Kalman Smoother
