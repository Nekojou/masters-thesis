this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)
source("test_commons.R")
source("test_simulationstudy1.R")

test.commons.runAll()
test.study1.runAll()