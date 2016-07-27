# setwd("~/Desarrollo/RProjects/Wth Simulation/1_setup/")

# PP2 = read.table('data/Salado_prcp.dat',header=T)
# MX2 = read.table('data/Salado_tmax.dat',header=T)
# MN2 = read.table('data/Salado_tmin.dat',header=T)
PP2 = read.table('data/Pampas_prcp.csv',header=T, sep = ',')
MX2 = read.table('data/Pampas_tmax.csv',header=T, sep = ',')
MN2 = read.table('data/Pampas_tmin.csv',header=T, sep = ',')

colnames(PP2) <- stringr::str_replace(colnames(PP2), 'pr\\.', '')
colnames(MX2) <- stringr::str_replace(colnames(MX2), 'tx\\.', '')
colnames(MN2) <- stringr::str_replace(colnames(MN2), 'tn\\.', '')

prcp <- reshape2::melt(PP2, id.vars='date', value.name='prcp')
tx <- reshape2::melt(MX2, id.vars='date', value.name='tx')
tn <- reshape2::melt(MN2, id.vars='date', value.name='tn')

colnames(prcp) <- c('date', 'station', 'prcp')
colnames(tx) <- c('date', 'station', 'tx')
colnames(tn) <- c('date', 'station', 'tn')

climate <- prcp %>% left_join(tx) %>% left_join(tn)

write.table(climate, file='data/Pampas.dat')

