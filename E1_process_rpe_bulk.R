
## Load libraries ####

library(data.table)

##

test <- fread(
    "../../../../../lcyril_data/sraDir/star/ReadsPerGene.out.tab",
    skip = 4
)

p <- ggplot(
    test,
    aes(
        x = log(V2),
        y = log(V4)
    )
) + geom_point()

sum(test$V2)

ggplot(
    test,
    aes(x = log10(V2 + 1))
) + geom_histogram(bins = 50)
