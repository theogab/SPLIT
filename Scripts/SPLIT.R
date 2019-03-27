## GSM Marine Tropical
require(raster)
require(sp)
require(GISTools)
require(Rcpp)

# Set directories
datadir <- 'Data'
resdir <- 'Results'
setwd("cluster/weekly/tgaboria")

Rcpp::sourceCpp("Scripts/kernels.cpp")

# Get arguments
args <- as.numeric(commandArgs(trailingOnly = TRUE))
cat("Launched SPLIT with the following arguments :\n")
ds <- args[1]
cat("ds = ", ds) 
d <- args[2]
cat("; d = ", d)
ps <- args[3]
cat("; ps = ", ps)
K <- args[4]
cat("; K = ", K)
rk <- args[5]
cat("; rk = ", rk)
pa <- args[6]
cat("; pa = ", pa)
rs <- args[7]
cat("; rs = ", rs)
xmin <- args[8]
cat("; xmin = ", xmin)
xmax <- args[9]
cat("; xmax = ", xmax)
ymin <- args[10]
cat("; ymin = ", ymin)
ymax <- args[11]
cat("; ymax = ", ymax)
TimeStart <- args[12]
cat("; TimeStart = ", TimeStart)
simnb <- args[13]
cat("; simnb = ", simnb, "\n")

#Load steps matrix
load(file.path(datadir, "Matrix_steps.rda"))
start_sp <- Matrix_steps[, c(1, 2, 141-TimeStart)]

#The rasters contains cells with different values:
#1 = Land
#2 = Deep water
#3 = Shallow water
#10 = Land tropical
#20 = Deep water tropical
#30 = Shallow water tropical

#The species can only be in 30. They will be able to cross deep water tropical 20, 
start_sp <- start_sp[ start_sp[ ,3 ] %in% c(3,30), ]
start_sp[,3] <- 1

#Create the distribution of the ancestral species
x <- start_sp[,1]
y <- start_sp[,2]
start_sp[x > xmax | x < xmin | y > ymax | y < ymin, 3] <- 0
colnames(start_sp) <- c("X","Y","1")

# When to save the species data.
saveSteps <- TimeStart:1

# Number of steps
steps <- TimeStart:1

dir <- paste0("d-", d,"_ds-", ds, "_ps-", ps, "_K-", K, "_pa-", pa, "_rk-", rk, "_rs-", rs)

if (!file.exists(file.path(resdir,dir))){
  dir.create(file.path(resdir, dir))
}

# Prepare the different simulations
sdir <- paste("Sim_nb", simnb, sep="")
dir.create(file.path(resdir, dir, sdir))

#################################
##### START THE SIMULATION ######
################################

# Saving the initial occurence matrix
PA <- start_sp
save(PA, file=file.path(resdir, dir, sdir, sprintf("PA_%s.rda",TimeStart+1)))

# Setting initial conditions
Nsp <- 1
Extinction <- c()
Hab <- Matrix_steps[, 141 - TimeStart]
Adaptation <- TestRange(as.matrix(PA[,3]), Hab[Hab %in% c(3,30)])

for(i in steps){
  #In case the speciation is unrealistic
  cat('i =', i,"; ")

    ## Get the transition matrices Internal for speciation and external for dispersion with cell changes
    Hab <- Hab[Hab %in% c(3,30)]
    NewHab <- Matrix_steps[, 142 - i]
    NewHab <- NewHab[NewHab %in% c(3,30)]
    load(paste(datadir, "/Distance_matrices_internal/Distance_", i, ".rda", sep=""))
    load(paste(datadir, "/Distance_matrices_external/Distance_", i, "_", i - 1, ".rda", sep=""))
    
    if(Nsp == 1){
      Extinction[sum(PA[, 3]) > 0] <- i
      Adaptation <- TestRange(as.matrix(PA[,3]), Hab)
    } else {
      Extinction[colSums(PA[, 3:(Nsp + 2)]) > 0] <- i
      Adaptation <- TestRange(as.matrix(PA[,-c(1,2)]), Hab)
    }
    
    cat('Ntrop =', sum(Adaptation == 1), '; Ntemp =', sum(Adaptation == 2), '; Nboth =', sum(Adaptation == 3), '; Ext =', sum(Adaptation == 0), '\n')
    
    # Simulating speciation and dispersion for all extant species at time i
    sp_disp <- SpDiv(ds, d, ps, as.matrix(PA[, -c(1, 2)]), dist1, dist2, K, rk, rs, pa, Hab, NewHab)
    
    # Updating Phylogeny with new speciation events
    descendants <- grep('d', colnames(sp_disp))
    nb_des <- sapply(1:(ncol(PA) - 2), function(x) length(grep(paste0('a', x, '_'), colnames(sp_disp))))
    ancestors <- rep(1:(ncol(PA) - 2), nb_des)
    age <- rep(i, length(ancestors))
    spec <- rep('allopatry', length(ancestors))
    spec[grep('sym', colnames(sp_disp)[descendants])] <- 'sympatry'
    
    if(!exists("phylogeny")){
      phylogeny <- data.frame(Ancestral = ancestors, Derived = descendants, Age = age, Speciation = spec)
    } else {
      
      phylogeny <- rbind(phylogeny, cbind(Ancestral = ancestors, Derived = descendants, Age = age, Speciation= spec))
    }
    
    # Preparing Matrix for next steps
    Nsp <- sum(colSums(sp_disp) > 0)
    Hab <- Matrix_steps[, 142 - i]
    PA <- cbind(Matrix_steps[Hab %in% c(3,30), 1:2], sp_disp)
    colnames(PA) <- c("X", "Y", as.character(1:ncol(sp_disp)))
    
    if (i %in% saveSteps){ 
      save(PA, file = file.path(resdir, dir, sdir, sprintf("PA_%s.rda",i-1)))
    }
}

# Preparing Matrix for next steps
Nsp <- sum(colSums(sp_disp) > 0)
Hab <- Matrix_steps[, 141]
load(file.path(datadir, "Distance_matrices_external/Distance_0_0.rda"))
sp_disp <- SpDiv(0, d, 0, sp_disp, dist1, dist2, K, rk, rs, pa, NewHab, NewHab)
PA <- cbind(Matrix_steps[Hab %in% c(3,30), 1:2], sp_disp)
colnames(PA) <- c("X", "Y", as.character(1:ncol(sp_disp)))

save(PA, file = file.path(resdir, dir, sdir, sprintf("PA_%s.rda",0)))

#Save the phylogeny as a text file
write.table(phylogeny, file.path(resdir, dir, sdir, "Phylogeny_raw.txt"), sep="\t")

# The last step Extinction and Adaptation record have not been saved yet
Extinction[((ncol(PA) - 2) - length(Extinction)):(ncol(PA) - 2)] <-  min(steps)
Hab <- Hab[Hab %in% c(3,30)]
Adaptation[colSums(PA[, -c(1,2)]) > 0] <- TestRange(as.matrix(PA[,-c(1,2)]), Hab)[colSums(PA[, -c(1,2)]) > 0]

# The species that went extinctint right after their appearence are NA
if(!all(!is.na(Extinction))){
  Extinction[is.na(Extinction)] <- phylogeny[which(is.na(Extinction)) - 1, 'Age']
}

#Save the Extinction record as a text file
Extinction <- data.frame(cbind(1:(ncol(PA) - 2), Extinction, Adaptation))
colnames(Extinction) <- c("Species","Survival","Habitat")
write.table(Extinction, file=file.path(resdir, dir, sdir, "Extinction.txt"), sep="\t")
