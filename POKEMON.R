#importing the dataset
library(purrr)
library(bnlearn)
library(igraph)
library(corrplot)
library(gRain)
library(graph)
library(catnet)
library(gRbase)
library(Rgraphviz)

pokemon=read.csv("C:\\Users\\Paola\\Desktop\\PROBABILISTIC MODELLING\\PROBABILISTIC MODELLING\\pokemon.csv", 
                   header = TRUE,sep = ',')

head(pokemon)

summary(pokemon)

length(unique(pokemon$classfication))
length(unique(pokemon$name))

#PREPROCESSING____________________________________________________________________________
sum(is.na(pokemon))
sum(is.na(pokemon$abilities))

sum(is.na(pokemon$name))
#replaced 40 values with respective mean
pokemon$weight_kg[is.na(pokemon$weight_kg)] = mean(pokemon$weight_kg,na.rm=TRUE)
pokemon$height_m[is.na(pokemon$height_m)] = mean(pokemon$height_m,na.rm=TRUE)

#98 values 

pokemon$is_genderless = with(
  pokemon, ifelse(is.na(pokemon$percentage_male), 1, 0))
pokemon$percentage_male[is.na(pokemon$percentage_male)]= 0

#not null but maybe useful for our analysis
pokemon$is_type2 = with(
  pokemon, ifelse(pokemon$type2 == "", 0, 1))
#chnge the capture rate

length(pokemon$weight_kg[which(pokemon$weight_kg > 300)])    #25 pokemon over 300kg
pokemon[pokemon$weight_kg[which(pokemon$weight_kg > 700)]]
pokemon$name[pokemon$weight_kg > 700,] 
pokemon$name[which(pokemon$weight_kg > 700)]  #[1] "Groudon"    "Giratina"   "Mudsdale"   "Cosmoem"    "Celesteela" "Guzzlord"  #confirmed by literature
pokemon$capture_rate <- sapply(pokemon$capture_rate,as.numeric)
pokemon$capture_rate[is.na(pokemon$capture_rate)] <- 30

#create number of abilities
pokemon$num_abilities = lengths(strsplit(pokemon$abilities, ","))


pokemon$attack <- sapply(pokemon$attack,as.numeric)
pokemon$base_egg_steps <- sapply(pokemon$base_egg_steps,as.numeric)
pokemon$base_happiness <- sapply(pokemon$base_happiness,as.numeric)
pokemon$base_total<- sapply(pokemon$base_total,as.numeric)
pokemon$sp_defense<- sapply(pokemon$sp_defense,as.numeric)
pokemon$experience_growth<- sapply(pokemon$experience_growth,as.numeric)
pokemon$hp<- sapply(pokemon$hp,as.numeric)
pokemon$sp_attack <- sapply(pokemon$sp_attack,as.numeric)
pokemon$defense <- sapply(pokemon$defense,as.numeric)
pokemon$speed <- sapply(pokemon$speed,as.numeric)
pokemon$weight_kg <- sapply(pokemon$weight_kg,as.numeric)
pokemon$generation <- sapply(pokemon$generation,as.numeric)
pokemon$is_legendary <- sapply(pokemon$is_legendary,as.numeric)
pokemon$is_type2 <- sapply(pokemon$is_type2,as.numeric)
pokemon$is_genderless <- sapply(pokemon$is_genderless,as.numeric)
pokemon$num_abilities <- sapply(pokemon$num_abilities,as.numeric)


corrplot(cor(pokemon[c(20:24,26:29,32,34:36,39:44)]), method = 'number', order = 'alphabet',type = "lower",tl.col = "black")
boxplot(pokemon[c(20,26,29,34:36)])

#Discretize the data_______________________
dpokemon2 = bnlearn::discretize(pokemon[c(20:24,26:29,34:36,39,40,44)], method = "interval",
                                #dpokemon = bnlearn::discretize(npokemn, method = "interval",
                                breaks= 4) #


dpokemon2 = bnlearn::discretize(pokemon[c(20:24,26:29,34:36,39,40,44)], method = "interval",breaks= 4) #
for (i in names(dpokemon2))
  levels(dpokemon2[, i]) = c( "very low","low", "medium","high")


dpokemon1 = bnlearn::discretize(pokemon[c(32,41:43)], method = "interval",breaks= 2) #
for (i in names(dpokemon1))
  levels(dpokemon1[, i]) = c( "low","high")

dpokemon=cbind(dpokemon2, dpokemon1)

dpokemon$type1 <-as.factor(pokemon$type1)
levels(dpokemon$percentage_male)
levels(dpokemon$type1)

#LEARNING PHASE____________________________________________________________________________________

whitelist=matrix(c(
  "attack","base_total",
  "defense","base_total",
  "sp_attack","base_total",
  "sp_defense","base_total",
  "speed","base_total",
   "hp","base_total",
   "percentage_male","is_genderless",
   "is_legendary","base_egg_steps"
   ),,2,byrow = TRUE)

colnames(whitelist)=c("from","to")
whitelist

#Bootstrap in hc model_______________________________________________________________________
set.seed(12)
str.diff_hc = boot.strength(dpokemon, R = 400, algorithm = "hc",  cpdag = TRUE,
                            algorithm.args = list(score="aic", whitelist = whitelist))

avg.diff_hc = averaged.network(str.diff_hc, threshold = 0.35)
strength.plot(avg.diff_hc, str.diff_hc, shape = "ellipse", highlight = list(arcs = whitelist, lwd = 0.35))

head(str.diff_hc)
hc_mat = amat(avg.diff_hc)
hc_mat

#Bootstrap in IAMB model_______________________________________________________________________
set.seed(12)
str.diff_iamb = boot.strength(dpokemon, R = 400, algorithm = "iamb", cpdag = TRUE,
                            algorithm.args = list(whitelist = whitelist))

avg.diff_iamb = averaged.network(str.diff_iamb, threshold = 0.35)
strength.plot(avg.diff_iamb, str.diff_iamb, shape = "ellipse", highlight = list(arcs = whitelist, lwd = 0.70))

head(str.diff_iamb)
iamb_mat = amat(avg.diff_iamb)
iamb_mat


#Bootstrap in MMHC modelL_______________________________________________________
set.seed(12)
str.diff_rsmax2 = boot.strength(dpokemon, R = 400, algorithm = "rsmax2",  
                                algorithm.args = list(whitelist = whitelist))

avg.diff_rsmax2 = averaged.network(str.diff_rsmax2, threshold = 0.35)
strength.plot(avg.diff_rsmax2, str.diff_rsmax2, shape = "ellipse", highlight = list(arcs = whitelist, lwd = 0.35))
rsmax2_mat = amat(avg.diff_rsmax2)
rsmax2_mat

#Union of the matrices
mat_sum = hc_mat + iamb_mat   + rsmax2_mat
mat_sum


bn_from_mat = function(adj_mat, tshl){
  adj = adj_mat
  adj[which(mat_sum < tshl)] = 0
  adj[which(mat_sum >= tshl)] = 1
  
  model = empty.graph(colnames(adj))
  amat(model) = adj
  
  return(model)
}

var = colnames(dpokemon)
e = empty.graph(var)
arcs(e) = whitelist


#model_1: threshold = 1
model_1 = bn_from_mat(mat_sum, 1)


Dmodel_1 = cextend(model_1)
#some unidirected arcs were created so try to test them e choose a direction

undirected.arcs(model_1)

# from              to               
# [1,] "attack"          "hp"             
# [2,] "capture_rate"    "percentage_male"
# [3,] "defense"         "sp_defense"     
# [4,] "defense"         "speed"          
# [5,] "hp"              "attack"         
# [6,] "percentage_male" "capture_rate"   
# [7,] "sp_defense"      "defense"        
# [8,] "speed"           "defense"        
# [9,] "generation"      "num_abilities"  
# [10,] "num_abilities"   "generation" 
#[11,] "generation"        "num_abilities"    
#[12,] "num_abilities"     "generation" 

choose.direction(model_1, dpokemon, arc = c("generation", "num_abilities"), debug = TRUE)
choose.direction(model_1, dpokemon, criterion = "aic",arc = c("generation", "num_abilities"), debug = TRUE)
choose.direction(model_1, dpokemon, criterion = "bic",arc = c("percentage_male", "capture_rate"), debug = TRUE)

UPDATE = set.arc(model_1,"generation","num_abilities")
UPDATE1 = set.arc(UPDATE,"percentage_male","capture_rate")
UPDATE2 = set.arc(UPDATE1,"percentage_male","type1")
UPDATE3 = drop.arc(UPDATE2,"attack","hp")
UPDATE4 = drop.arc(UPDATE3,"defense","speed")
UPDATE5=drop.arc(UPDATE4,"defense","sp_defense")

undirected.arcs(UPDATE5)

#"sp_defense"        "defense
choose.direction(model_1, dpokemon, arc = c("base_happiness", "experience_growth"), debug = TRUE)
choose.direction(model_1, dpokemon, criterion = "bde",arc = c("base_happiness", "experience_growth"), debug = TRUE)

Dmodel_1 = cextend(UPDATE5)

choose.direction(UPDATE2, dpokemon,criterion = "aic", arc = c("defense", "sp_defense"), debug = TRUE)
#model_2: threshold = 2
model_2 = bn_from_mat(mat_sum, 2)
model_2 = cextend(model_2)

#model_3: threshold = 3

model_3 = bn_from_mat(mat_sum, 3)
model_3 = cextend(model_3)

par(mfrow = c(1,2))
G_compare <-graphviz.compare(e, model_3, model_2, Dmodel_1, shape = 'ellipse',
                 main = c("knowledge based", "knowledge based + th = 3", "knowledge based + th = 2", "knowledge based + th = 1"))



cv.bic_base = bn.cv(dpokemon, bn = e, runs = 10, 
                    algorithm.args = list(score = "bic"))

cv.bic_1 = bn.cv(dpokemon, bn = Dmodel_1, runs = 10, 
                 algorithm.args = list(score = "bic"))

cv.bic_2 = bn.cv(dpokemon, bn = model_2, runs = 10, 
                 algorithm.args = list(score = "bic"))

cv.bic_3 = bn.cv(dpokemon, bn = model_3, runs = 10, 
                 algorithm.args = list(score = "bic"))


plot(cv.bic_base, cv.bic_1, cv.bic_2, cv.bic_3, xlab = c("BIC KB", "BIC t=1", "BIC t = 2", "BIC t = 3"))

losses1=c(mean(loss(cv.bic_base)),mean(loss(cv.bic_1)),mean(loss(cv.bic_2)),mean(loss(cv.bic_3)))
losses=data.frame(losses1)
rownames(losses)=c('base','1','2','3')
losses

#losses1
# losses1
# base 16.27250
# 1    14.34642
# 2    15.54045
# 3    15.87061

fit = bn.fit(Dmodel_1, dpokemon)
fit

par(mfrow = c(1,1))
gR = graphviz.plot(fit, layout = "dot", shape = "ellipse", 
                   highlight = list(nodes = c("type1","base_total"), col = c("green"), fill = "green"))

node.attrs = nodeRenderInfo(gR)
#evidence 
#The other important limitation of gRain, compared to bnlearn, is that it does not allow for conditional probabilities to be NaN.
type1_values = unique(dpokemon$type1)
experience_values = unique(dpokemon$experience_growth)
type2_values = unique(dpokemon$is_type2)
capturerate_values = unique(dpokemon$capture_rate)
generation_values = unique(dpokemon$generation)

query <- cpquery(fit, (percentage_male == 'low'),(base_total == 'high'),method = 'ls', debug = TRUE)
query


query2 <- cpquery(fit, (is_type2 == "high"),(type1 == "normal" & experience_growth == "very low"),method = 'ls', debug = TRUE)
query2

query3- cpquery(fit, (is_type2 == "high"),(type1 == "normal"))
#query about type                           
                                             
type2_h_G_type1 = c()
for (i in 1:length(type1_values)){
  value = toString(i)
  prob = cpquery(fit, (is_type2 == "high") , type1 == type1_values[i])
  type2_h_G_type1[i] = prob
}


type2_l_G_type1 = c()
for (i in 1:length(type1_values)){
  value = toString(i)
  prob = cpquery(fit, (is_type2 == "low") , type1 == type1_values[i])
  type2_l_G_type1[i] = prob
}
# > type1_values
# [1] grass    fire     water    bug      normal   poison   electric ground   fairy    fighting psychic  rock    
# [13] ghost    ice      dragon   dark     steel    flying  
# 18 Levels: bug dark dragon electric fairy fighting fire flying ghost grass ground ice normal poison ... water
# > type2_values
# [1] high low 
# Levels: low high
#check 
querycheck = cpquery(fit, (is_type2 == "low") , type1 == "bug")
querycheck



df <-cbind.data.frame(type2_h_G_type1,type2_l_G_type1, row.names = as.vector(type1_values))
df
g_range = range(0.0, 1.0)
barplot(t(as.matrix(df)), beside = TRUE, col = c("green", "red"),ylim = g_range,legend.text = c("Yes","No"))


##query 2

#queying about capture_rate

queryCgB=cpquery(fit, (capture_rate == "high") , (base_total== "high" ))

cpquery(fit, (capture_rate == "high") , (base_total== "very low" ))
#it is confirmed that the strongest a pokemon is the more difficult to hatch is

CP_L_G_%male_and_bt = c()
for (i in 1:length(type1_values)){
  value = toString(i)
  prob = cpquery(fit, (capture_rate == "very low") , (percentage_male == "low" ))
  type2_h_G_type1[i] = prob
}

 cpquery(fit, (capture_rate == "very low") , (percentage_male == "high" ))
# 0.8063837 seems that if a pokemon is male it has a low capture rate
 cpquery(fit, (capture_rate == "very low") , (percentage_male == "high" & base_total == "very low"))
 #[1] 0.6046512 drop quite
cpquery(fit, (capture_rate == "high") , (percentage_male == "high" ))
 cpquery(fit, (capture_rate == "very low") , (percentage_male == "low" ))
#[1] 0.4055761  genderless or more female easily to be captured
cpquery(fit, (capture_rate == "very low") , (percentage_male == "low" & base_total == "very low"))
#0.07285145

 cpquery(fit, (capture_rate == "very low") , (percentage_male == "low" & base_total == "high"))
# [1] 1
 
cpquery(fit, (is_legendary == "high") , (is_genderless == "high"))
 #[1] 0.6673718 confirmed 


cpquery(fit, (base_total == "high") , (hp == "high" & sp_attack == "low" & sp_defense== "low"))
#[1] 0.06451613

##query about generation 
generationG_bt = c()
for (i in 1:length(generation_values)){
  value = toString(i)
  prob = cpquery(fit, (generation == generation_values[i]) , (base_total== "high" ))
  generationG_bt[i] = prob
}
generationG_bt 
#firsts generation seems to have an higher bt [1] 0.2960000 0.3228070 0.1887550 0.1742424
#Generation VII added the most Pokémon species whose gender is unknown, with a total of 29.

typeG_bt = c()
for (i in 1:length(type1_values)){
  value = toString(i)
  prob = cpquery(fit, (type1 == type1_values[i]) , (base_total== "high" ))
  typeG_bt[i] = prob
}
typeG_bt
#grass water  normal

typeG_bt2 = c()
for (i in 1:length(type1_values)){
  value = toString(i)
  prob = cpquery(fit, (type1 == type1_values[i]) , (base_total== "medium" ))
  typeG_bt2[i] = prob
}
typeG_bt2 
#grass bug normal

typeG_bt3 = c()
for (i in 1:length(type1_values)){
  value = toString(i)
  prob = cpquery(fit, (type1 == type1_values[i]) , (base_total== "low" ))
  typeG_bt2[i] = prob
}
typeG_bt3 




query4 = cpquery(fit, (is_type2 == "high") , (base_total== "high" ))
query4
#[1] 0.4882812
#overall queries
typeG_ALL = c()
for (i in 1:length(type1_values)){
 value = toString(i)
prob = cpquery(fit, (type1 == type1_values[i]) , (experience_growth == 'medium' &
                      base_total== "high" & is_type2 == 'high' ))
typeG_ALL[i] = prob
}
typeG_ALL
#water   ghst rock  steel

typeG_ALL = c()
for (i in 1:length(type1_values)){
     value = toString(i)
    prob = cpquery(fit, (type1 == type1_values[i]) , (
       base_total== "medium"& base_egg_steps== 'medium'))
     typeG_ALL[i] = prob
   }
typeG_ALL

#water normal electric ground dragon ice
typeG_ALL = c()
for (i in 1:length(type1_values)){
    value = toString(i)
    prob = cpquery(fit, (type1 == type1_values[i]) , (experience_growth == 'medium' &
                         base_total== "high"& base_egg_steps== 'high'))
    typeG_ALL[i] = prob
  }
typeG_ALL

#grass steel dragon rock


typeG_ALL = c()
 for (i in 1:length(type1_values)){
    value = toString(i)
     prob = cpquery(fit, (type1 == type1_values[i]) , (experience_growth == 'medium' &
            base_total== "high" & is_legendary == 'high' ))
    typeG_ALL[i] = prob
 }
typeG_ALL

#grass water normal physic steel e dragon 
