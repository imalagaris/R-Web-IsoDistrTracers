isotope <- readRDS(file = "./data/isotope.RDS")
dat <- readRDS(file = "./data/dat.RDS")
dat_wt <- dat[dat$atoms != "", c("atoms", "weight")]
