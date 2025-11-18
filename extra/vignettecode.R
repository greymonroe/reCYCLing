
source("R/functions.R")
# test --------------------------------------------------------------------

startseq<-paste0(sample(c("A", "C", "G", "T"), size = 178, replace = TRUE), collapse = "")
ps_resultsmu3 <- run_sim_ps(init_sequence_type = "given",
                         ancestor_seq=startseq,
                         init_l=178,
                         n = 5,
                         verbose = F,
                         p_del_chunk=0.000,
                         init_k0 = 10,
                         max_t = 1000000,
                         mu_total = 0.00003)

ps_resultsmu2 <- run_sim_ps(init_sequence_type = "given",
                            ancestor_seq=startseq,
                            init_l=178,
                            n = 5,
                            verbose = F,
                            p_del_chunk=0.000,
                            init_k0 = 10,
                            max_t = 1000000,
                            mu_total = 0.00002)


ps_resultsmu1 <- run_sim_ps(init_sequence_type = "given",
                            ancestor_seq=startseq,
                            init_l=178,
                            n = 5,
                            verbose = F,
                            p_del_chunk=0.000,
                            init_k0 = 10,
                            max_t = 1000000,
                            mu_total = 0.00001)

ps_resultsmu9 <- run_sim_ps(init_sequence_type = "given",
                            ancestor_seq=startseq,
                            init_l=178,
                            n = 5,
                            verbose = F,
                            p_del_chunk=0.000,
                            init_k0 = 10,
                            max_t = 1000000,
                            mu_total = 0.00009)

ps_resultsmu3nolong <- run_sim_ps(init_sequence_type = "given",
                            ancestor_seq=startseq,
                            init_l=178,
                            n = 5,
                            verbose = F,
                            p_del_chunk=0.000,
                            p_distal_dup = 0,
                            init_k0 = 10,
                            max_t = 1000000,
                            mu_total = 0.00003)

plot_ps_summary(ps_resultsmu3nolong, i = 1,
                     title = "mu 3e-5, long=0",
                     rel_widths_bottom = c(1.2, 1, 1, 1),
                     rel_heights = c(3, 1))

plot_ps_summary(ps_resultsmu3, i = 1,
                title = "mu 3e-5",
                rel_widths_bottom = c(1.2, 1, 1, 1),
                rel_heights = c(3, 1))

plot_ps_summary(ps_resultsmu2, i = 1,
                title = "mu 2e-5",
                rel_widths_bottom = c(1.2, 1, 1, 1),
                rel_heights = c(3, 1))

plot_ps_summary(ps_resultsmu1, i = 2,
                title = "mu 1e-5",
                rel_widths_bottom = c(1.2, 1, 1, 1),
                rel_heights = c(3, 1))

plot_ps_summary(ps_resultsmu9, i = 1,
                title = "mu 9e-5",
                rel_widths_bottom = c(1.2, 1, 1, 1),
                rel_heights = c(3, 1))


