############################################## uncorrected regressions
#regressions.pdl=regressions.pdl

##pdl 

regressions.pdl <- Qe1.2_pdl %>%
  gather(v, value, -CID,) %>%
  group_by(CID) %>% 
  nest() %>%
  # mutate(data = map(data, ~ .x %>%  mutate(pdl_n = c(rep(.16,4),rep(.26,4),rep(.43,4),rep(.56,4),rep(.86,4),rep(1,2))))) %>% # add numeric value for PDL or tert TP to percent of timecourse between PDL20 to PDL 50
  mutate(data = map(data, ~ .x %>%  mutate(pdl_n = c(rep(25,4),rep(28,4),rep(32,4),rep(37,4),rep(46,4),rep(50,2))))) %>% # add numeric value for PDL or tert TP to percent of timecourse between PDL20 to PDL 50
  
  mutate(
    lin.fit = map(data, ~lm(log(value) ~ pdl_n, data = .)),
    lin.tidied = map(lin.fit, tidy),
    lin.glanced = map(lin.fit, glance))  %>%
  mutate(
    log.fit = map(data, ~lm(log(value) ~ pdl_n + (pdl_n)^2, data = .)),
    log.tidied = map(log.fit, tidy),
    log.glanced = map(log.fit, glance)) #%>%
# unnest(lin.glanced,names_sep = ".") %>%
#unnest(log.glanced,names_sep = ".")


#regressions.pdl.lin=regressions %>%
# unnest(lin.glanced,names_sep = ".")  

##### object with slopes and pvalues for each metabolite
regressions.pdl.log= regressions.pdl %>%
  unnest(log.tidied,names_sep = ".")  %>%
  filter(log.tidied.term == "pdl_n")


regressions.pdl.lin= regressions.pdl %>%
  unnest(lin.tidied,names_sep = ".")  %>%
  filter(lin.tidied.term == "pdl_n")



#### get fdr values

pvalues=regressions.pdl.log$log.tidied.p.value
qobj.pdl <- qvalue(p = pvalues)
regressions.pdl.log$fdr=qobj.pdl$lfdr
regressions.pdl.log$qval=qobj.pdl$qvalues


regressions.pdl.log=regressions.pdl.log %>% dplyr::select(c("CID","log.tidied.estimate","log.tidied.std.error","log.tidied.statistic","log.tidied.p.value","fdr","qval"))

sig.metab.pdl.up=regressions.pdl.log %>% filter(qval < 0.001 & log.tidied.estimate > 0)  %>% pull(CID)
#ig.metab.pdl.down=regressions.pdl.log %>% filter(fdr < 0.05 & log.tidied.estimate < 0)  %>% pull(CID)

sig.metab.pdl.down=regressions.pdl.log %>% filter(qval < 0.001 & log.tidied.estimate < 0)  %>% pull(CID)




#####################
pvalues=regressions.pdl.lin$lin.tidied.p.value
qobj.pdl <- qvalue(p = pvalues)
regressions.pdl.lin$fdr=qobj.pdl$lfdr
regressions.pdl.lin$qval=qobj.pdl$qvalues


regressions.pdl.lin=regressions.pdl.lin %>% dplyr::select(c("CID","lin.tidied.estimate","lin.tidied.std.error","lin.tidied.statistic","lin.tidied.p.value","fdr","qval"))

sig.metab.pdl.up=regressions.pdl.lin %>% filter(qval < 0.001 & lin.tidied.estimate > 0)  %>% pull(CID)
#ig.metab.pdl.down=regressions.pdl.lin %>% filter(fdr < 0.05 & lin.tidied.estimate < 0)  %>% pull(CID)

sig.metab.pdl.down=regressions.pdl.lin %>% filter(qval < 0.001 & lin.tidied.estimate < 0)  %>% pull(CID)



View(regressions.pdl.log)


####################Now TERT


regressions.tert <- Qe1.2_tert %>%
  gather(v, value, -CID,) %>%
  group_by(CID) %>% 
  nest() %>%
  mutate(data = map(data, ~ .x %>% mutate(pdl_n = c(rep(.16,4),rep(.43,4),rep(.56,4),rep(.86,4),rep(1,3))))) %>% # add numeric value for PDL or tert TP to percent of timecourse between PDL20 to PDL 50
  mutate(
    lin.fit = map(data, ~lm(value ~ pdl_n, data = .)),
    lin.tidied = map(lin.fit, tidy),
    lin.glanced = map(lin.fit, glance))  %>%
  mutate(
    log.fit = map(data, ~lm(log(value) ~ pdl_n, data = .)),
    log.tidied = map(log.fit, tidy),
    log.glanced = map(log.fit, glance)) #%>%
# unnest(lin.glanced,names_sep = ".") %>%
#unnest(log.glanced,names_sep = ".")


#regressions.pdl.lin=regressions %>%
# unnest(lin.glanced,names_sep = ".")  

##### object with slopes and pvalues for each metabolite
regressions.tert.lin= regressions.tert %>%
  unnest(lin.tidied,names_sep = ".")  %>%
  filter(lin.tidied.term == "pdl_n")

#### get fdr values

pvalues=regressions.tert.lin$lin.tidied.p.value
qobj.tert <- qvalue(p = pvalues)
regressions.tert.lin$fdr=qobj.tert$lfdr
regressions.tert.lin$qval=qobj.tert$qvalues

regressions.tert.lin=regressions.tert.lin %>% dplyr::select(c("CID","lin.tidied.estimate","lin.tidied.std.error","lin.tidied.statistic","lin.tidied.p.value","fdr","qval"))



sig.metab.pdl.up=regressions.pdl.lin %>% filter(qval < 0.001 & lin.tidied.estimate > 0)  %>% pull(CID)
sig.metab.pdl.down=regressions.pdl.lin %>% filter(qval < 0.001 & lin.tidied.estimate < 0)  %>% pull(CID)


sig.metab.tert.up=regressions.tert.lin %>% filter(qval < 0.001 & lin.tidied.estimate > 0)  %>% pull(CID)
sig.metab.tert.down=regressions.tert.lin %>% filter(qval < 0.001 & lin.tidied.estimate < 0)  %>% pull(CID)

trd.up=setdiff(sig.metab.pdl.up,sig.metab.tert.up)

trd.down =setdiff(sig.metab.pdl.down,sig.metab.tert.down)


################## FC for plotting



