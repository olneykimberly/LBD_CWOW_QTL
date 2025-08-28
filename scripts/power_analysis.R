if (!requireNamespace("genpwr", quietly = TRUE)) {
  install.packages("genpwr")
}
library(genpwr)

meta <- read.delim("metadata/filtered_metadata.txt")
table(meta$TYPE)


# Define parameters
sample_size <- 589
maf <- 0.1
alpha <- .05
or <- 1.5

help(genpwr.calc)
# Calculate power
genpwr.calc(calc = "power", 
            model = "logistic", 
            N = sample_size, 
            MAF = 0.1, 
            Alpha = alpha, 
            k = .4, # number of controls per cases
            OR = 2)
