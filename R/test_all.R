setwd(dirname(parent.frame(2)$ofile))
source("test_commons.R")
source("test_simulationstudy1.R")

test.commons.runAll()
test.study1.runAll()