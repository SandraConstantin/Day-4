library(TraMineR)
data(biofam)
biofam$cohort <- cut(biofam$birthyr, c(1900,1930,1940,1950,1960),labels=c("1900-1929", "1930-1939", "1940-1949", "1950-1959"), right=FALSE)
biofam.states <- c("Parent", "Left", "Married", "Left/Married", "Child", "Left/Child","Left/Married/Child", "Divorced")
biofam.shortlab <- c("P","L","M","LM","C","LC", "LMC", "D")
biofam.seq <- seqdef(biofam[,10:25], states=biofam.shortlab, labels=biofam.states,weights=biofam$wp00tbgs)

biofam_table <- data.frame(seqlength(biofam.seq), seqtransn(biofam.seq), seqsubsn(biofam.seq), seqient(biofam.seq), seqST(biofam.seq), seqici(biofam.seq))
head(biofam_table)

summary(biofam_table)

par (mfrow = c(2,3))
hist(seqtransn(biofam.seq), col="red")
hist(seqsubsn(biofam.seq), col="green")
hist(seqient(biofam.seq), col="purple")
hist(seqST(biofam.seq), col="LightBlue")
hist(seqici(biofam.seq), col="yellow")

tail(seqdss(biofam.seq))
tail(seqdur(biofam.seq))

summary(apply(seqdur(biofam.seq),1,mean,na.rm=TRUE))
summary(apply(seqdur(biofam.seq),1,var,na.rm=TRUE))

names(biofam_table) <- c("length","transitions","subsequences","entropy","turbulence","complexity")
plot(biofam_table[,c("entropy","turbulence","complexity")])

boxplot(biofam_table$complexity ~ biofam$cohort, col="blue")

lm_complexity <- lm(biofam_table$complexity ~ cohort + sex + plingu02, data = biofam)
summary(lm_complexity)

