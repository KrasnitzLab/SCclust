

data_dir <- Sys.getenv("SGAINS_DATA")
sim_filename <- file.path(data_dir, "nyu003/results/nyu003.benign.1simP.txt")
expect_true(file.exists(sim_filename))

sim_mat <- as.matrix(scan(sim_filename))
a = sim_mat[1:10,1]


pinmat_filename <- file.path(data_dir, "nyu003/results/nyu003.benign.1smear1bpPinMat.txt")
expect_true(file.exists(pinmat_filename))
pinmat_df <- load_table(pinmat_filename)

pins_filename <- file.path(data_dir, "nyu003/results/nyu003.benign.1smear1bpPins.txt")
expect_true(file.exists(pins_filename))
pins_df <- load_table(pins_filename)

# res <- sim_fisher_wrapper(pinmat_df, pins_df, njobs=30, nsim=500, nsweep=200)
res <- sim_fisher_wrapper(pinmat_df, pins_df, njobs=30, nsim=2, nsweep=200)

sim_res <- res$sim
dim(sim_res)
d <- nrow(sim_res)*ncol(sim_res)
d
a <- sim_mat[1:d,1]
dim(a) <- dim(sim_res)
expect_equal(a, sim_res, tolerance=0.5)
