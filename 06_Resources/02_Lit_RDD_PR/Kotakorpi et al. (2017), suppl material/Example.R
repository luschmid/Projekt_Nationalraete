

# Example for party BA in district 5 in year 2007
# Kotakorpi, K., Poutvaara, P. and Tervi??, M. (2017) ???Returns to Office in National and Local Politics: A Bootstrap Method and Evidence from Finland???, The Journal of Law, Economics, and Organization, 33(3), pp. 413???442. 

M=10000
m=100

votes1=4024
votes2=4388

votes_all=c(rep(0,votes1),rep(1,votes2))

p1=votes1/(votes1+votes2)
(p2=votes2/(votes1+votes2))

elected <- matrix(NA,nrow = M, ncol = 1)

for (i in 1:M){
votes <- sample(votes_all,m,replace=TRUE)
elected[i] <- round(mean(votes))
}

mean(elected)