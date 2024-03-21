library(randomForestSRC)
library(raster)
library(dismo)
library(spdep)
library(viridis)
library(gridExtra)
library(stats)
library(colorspace)
library(dplyr)
library(ggplot2)
library(sf)
library(rnaturalearth)
library(binsmooth)
library(RColorBrewer)
library(akima)
library(fields)


dir.create('FIGURES')
dir.create('RESULTS')
set.seed(19820408)
###
stdize<-function(x){
  ul<-max(x,na.rm=T)
  ll<-min(x,na.rm=T)
  d<-ul-ll
  return((x-ll)/d)
  }

###

radians <- function(deg) {(deg * pi) / (180)}
deg_to_km<-function(lat){
  res<-1
  R<-6371.0072
  return ((sin(radians(lat+(res/2.0)))-sin(radians(lat-(res/2.0))))*radians(res/2.0)*(R**2))
}



###CONVERT GDP TO GNI
gni_to_gdp<-read.csv('./INPUT/GDP/gni_vs_gdp.csv',header=T) ###data from https://ourworldindata.org/grapher/gni-per-capita-vs-gdp-per-capita?tab=table
plot(gni_to_gdp$gni,gni_to_gdp$gdp)

gni_mod<-lm(gni_to_gdp$gdp~gni_to_gdp$gni)
summary(gni_mod)
gdp_a<-gni_mod$coefficients[1]
gdp_b<-gni_mod$coefficients[2]
####

aa<-read.csv('./OUTPUT/matched_climate_crisis.csv',header=T)
aa<-aa[is.finite(rowSums(aa[,6:ncol(aa)])),]
aa<-aa[aa$fam_lev!=-99999,]
aa<-aggregate(.~year+row+col,data=aa,FUN='max')

gni_tre<-727.48

gni_safe<-1*(((gdp_b*aa$gdp/aa$pop)+gdp_a)>gni_tre)

sum(gni_safe*as.numeric(as.character(aa$fam_lev)),na.rm=T)/sum(as.numeric(as.character(aa$fam_lev)))



###link coordinates to countries

lat_lon<-data.frame('latitude'=aa$row,'longitude'=aa$col,id=1:nrow(aa))
dim(unique(cbind(aa$row,aa$col)))
length(unique(aa$year))
sum(as.numeric(as.character(aa$fam_lev)))
sf::sf_use_s2(FALSE)
sfRegion <- st_as_sf(lat_lon, coords = c("longitude","latitude"),crs=4326)
sfCountry <- ne_countries(returnclass='sf')
sfJoined <- st_join(sfRegion,sfCountry) ###MODIFIED, REVERSED JOINT
matched_countries<-data.frame('id'=sfJoined$id,'continent'=sfJoined$continent,'country'=sfJoined$name_long)


cont_area<-data.frame('continent'=world$continent,'area'=world$area_km2)
cont_area<-aggregate(area~continent,data=cont_area,FUN='sum')

countries<-unique(matched_countries$country)

cell_areas<-deg_to_km(lat_lon$latitude)
all_res_continent<-data.frame('year'=aa$year[matched_countries$id],'country'=matched_countries$country,
                              'continent'=matched_countries$continent,'new_crisis'=aa$fam_lev[matched_countries$id],
                              'continued_crisis'=aa$crisis_continued_val[matched_countries$id],
                              'crisis'=aa$crisis_famine_val[matched_countries$id],
                              'area'=cell_areas[matched_countries$id],'pop'=aa$pop_density[matched_countries$id])


all_res_continent<-cbind(all_res_continent,'cont_area'=cont_area[match(all_res_continent$continent,cont_area$continent),2])
all_res_continent$new_crisis<-(as.numeric(as.character(all_res_continent$new_crisis)))*(all_res_continent$pop/1000000)
all_res_continent$continued_crisis<-(as.numeric(as.character(all_res_continent$continued_crisis)))*(all_res_continent$pop/1000000)
all_res_continent$crisis<-(as.numeric(as.character(all_res_continent$crisis)))*(all_res_continent$pop/1000000)
all_res_continent<-all_res_continent[is.finite(all_res_continent$year),]

cont_sums<-aggregate(cbind(all_res_continent$new_crisis,all_res_continent$continued_crisis,all_res_continent$crisis)~all_res_continent$continent+all_res_continent$year,FUN='sum')
colnames(cont_sums)<-c('continent','year','exposed_new','exposed_continued','exposed_crisis')



p1<-ggplot(cont_sums,                                      # Grouped barplot using ggplot2
           aes(x = factor(year),
               y = exposed_new,
               fill = continent)) +
  ggtitle("new crisis")+
  ylim(0,140)+
  scale_color_manual(values = viridis(3), guide = "none") +
  scale_fill_manual(values = viridis(3)) + xlab("year") + ylab("exposed population (million people)")+
  geom_bar(stat = "identity")+           #,position = "dodge")
  theme(legend.position="top")


p2<-ggplot(cont_sums,                                      # Grouped barplot using ggplot2
           aes(x = factor(year),
               y = exposed_crisis,
               fill = continent)) +
  ggtitle("ongoing crisis")+
  ylim(0,140)+
  scale_color_manual(values = viridis(3), guide = "none") +
  scale_fill_manual(values = viridis(3)) + xlab("year") + ylab("exposed population (million people)")+
  geom_bar(stat = "identity")+ theme(legend.position = "none")+theme(legend.position="top")




pdf('./FIGURES/exposed_people_by_continent_from_dataset.pdf',width = 12,height=4.5)
grid.arrange(p1, p2, nrow = 1)
dev.off()





###make cross spatial validation by splitting the dataset by countries and year

t_cols<-match(paste0('t',0:23),colnames(aa))
p_cols<-match(paste0('p',0:23),colnames(aa))
sum_t<-rowSums(aa[,t_cols])
sum_p<-rowSums(aa[,p_cols])
sd_t<-apply(aa[,t_cols],1,sd)
sd_p<-apply(aa[,p_cols],1,sd)
max_t<-apply(aa[,t_cols],1,max)
max_p<-apply(aa[,p_cols],1,max)
min_t<-apply(aa[,t_cols],1,min)
min_p<-apply(aa[,p_cols],1,min)
range_t<-max_t-min_t
range_p<-max_p-min_p
aa<-cbind(aa,sum_t,sum_p,sd_t,sd_p,max_t,max_p,min_t,min_p,range_t,range_p)
add_var_ids<-match(c("sum_t","sum_p","sd_t","sd_p","max_t","max_p","min_t","min_p","range_t","range_p"),colnames(aa))

dep_var_name<-'fam_lev'#'crisis_famine_val'#
dep_var<-which(colnames(aa)==dep_var_name)
aa<-aa[which(aa[,dep_var]!=-99999),]
aa$id<-1:nrow(aa)

sf::sf_use_s2(FALSE)
sfRegion <- st_as_sf(aa, coords = c("col","row"),crs=4326)
sfCountry <- ne_countries(returnclass='sf')
sfJoined <- st_join(sfRegion,sfCountry)
sfJoined<-sfJoined[which(!is.na(sfJoined$id)),]

aa<-aa[sfJoined$id,]
aa$locs<-sfJoined$name_long
aa$continent<-sfJoined$continent
locs<-unique(aa$locs)

aa[,dep_var]<-as.factor(aa[,dep_var])
aa$continent<-as.factor(aa$continent)


crisis<-which(aa[,dep_var]==1)
not_crisis<-which(aa[,dep_var]==0)
years<-unique(aa$year)
null_vals<-which(aa[,dep_var]==-99999)



y_loc<-paste0(round(aa$row/1),'_',round(aa$col/1))
y_loc_u<-unique(y_loc)
length(y_loc_u)

DT <- data.frame(longitude=aa$col,
                  latitude=aa$row+runif(nrow(aa))/10000)



pts = st_as_sf(DT, coords = c("longitude", "latitude"), 
                 crs = 4326, agr = "constant")



models<-list()
for (replicate in 1:100){
  for (loc_frac in seq(0.1,0.9,0.1)){
      for (smote in seq(0.5,3,0.5)){
       for (sample_frac in seq(0.5,0.75,0.1)){
          rand_locs<-sample(y_loc_u,length(y_loc_u)*loc_frac)
          train_ids<-which(!(y_loc%in%rand_locs))
          test_ids<-setdiff(1:nrow(aa),c(train_ids,null_vals))
          if (length(test_ids)>15000){
           test_ids<-sample(test_ids,length(test_ids))[1:15000]}
          crisis_<-intersect(train_ids,crisis)
          if (length(crisis_)>0){
            not_crisis_<-intersect(train_ids,not_crisis)
            samp_size<-length(crisis_)*sample_frac
            train_ids<-c(sample(crisis_,samp_size),sample(not_crisis_,samp_size*smote)) #SMOTE: Synthetic Minority Over-sampling Technique
            train_data<-aa[train_ids,]
            test_data<-aa[test_ids,]
            rf_train<-train_data[,c(dep_var,t_cols,p_cols,add_var_ids)]
            rf_test<-test_data[,c(dep_var,t_cols,p_cols,add_var_ids)]
            colnames(rf_train)[1]<-'fam_lev'
            colnames(rf_test)[1]<-'fam_lev'
            rf <-rfsrc(as.factor(fam_lev)~.,data=rf_train,ntree=250)
            pred_data<-predict(rf,rf_test)$predicted[,2]
            err<-c()
            for (tre in seq(0.3,0.8,0.1)){
              pred_test<-(1*(pred_data>tre))-as.numeric(as.character(rf_test$fam_lev))
              fn<-sum(pred_test==-1)/sum(rf_test$fam_lev==1)
              fp<-sum(pred_test==1)/sum(rf_test$fam_lev==0)
              tss<-(1-fp)+(1-fn)-1
              err<-rbind(err,c(tre,fp,fn,tss))
            }
            best_tre<-err[which(err[,4]==max(err[,4])),]
            models[[length(models)+1]]<-list('train_data'=train_ids,
                                             'tre'=best_tre[1],
                                             'false_pos'=best_tre[2],
                                             'false_neg'=best_tre[3],
                                             'tss'=best_tre[4],
                                             'loc_f' = loc_frac,
                                             'smote' = smote,
                                             'sample_f' = sample_frac
                                             )
              tss<-as.numeric(lapply(models,function(x) x$tss))
              print (paste(replicate,length(models),paste('fp=',round(best_tre[2],2),
                                                'fn=',round(best_tre[3],2),
                                                'tss=',round(best_tre[4],2),
                                                'max_tss=',round(max(tss),2)
                                                
                                                )))

            } 
            
            }
          }
        }
      }



#identify the best parametrization, i.e. the one leading to the higher mean TSS
tss<-as.numeric(lapply(models,function(x) x$tss))
smote<-as.numeric(lapply(models,function(x) x$smote))
loc_f<-as.numeric(lapply(models,function(x) x$loc_f))
sample_f<-as.numeric(lapply(models,function(x) x$sample_f))
eval_mod<-aggregate(tss~smote+loc_f+sample_f,FUN='mean')
best_pars<-eval_mod[which(eval_mod$tss==max(eval_mod$tss)),]
smote<-as.numeric(best_pars[1])
loc_frac<-as.numeric(best_pars[2])
sample_frac<-as.numeric(best_pars[3])

###now obtain 100 models with the best parametrization and no spatial autocorrelation
best_models<-list()
n_spatial_bias<-0
n_tried_models<-0
while (length(best_models)<100){
    n_tried_models<-n_tried_models+1
    rand_locs<-sample(y_loc_u,length(y_loc_u)*loc_frac)
    train_ids<-which(!(y_loc%in%rand_locs)) 
    test_ids<-setdiff(1:nrow(aa),c(train_ids,null_vals))
    if (length(test_ids)>15000){
      test_ids<-sample(test_ids,length(test_ids))[1:15000]}
    crisis_<-intersect(train_ids,crisis)
    if (length(crisis_)>0){
      not_crisis_<-intersect(train_ids,not_crisis)
      samp_size<-length(crisis_)*sample_frac
      train_ids<-c(sample(crisis_,samp_size),sample(not_crisis_,samp_size*smote)) #SMOTE: Synthetic Minority Over-sampling Technique
      train_data<-aa[train_ids,]
      test_data<-aa[test_ids,]
      rf_train<-train_data[,c(dep_var,t_cols,p_cols,add_var_ids)]
      rf_test<-test_data[,c(dep_var,t_cols,p_cols,add_var_ids)]
      colnames(rf_train)[1]<-'fam_lev'
      colnames(rf_test)[1]<-'fam_lev'
      rf <-rfsrc(as.factor(fam_lev)~.,data=rf_train,ntree=500)
      pred_data<-predict(rf,rf_test)$predicted[,2]
      err<-c()
      for (tre in seq(0.3,0.8,0.1)){
        pred_test<-(1*(pred_data>tre))-as.numeric(as.character(rf_test$fam_lev))
        fn<-sum(pred_test==-1)/sum(rf_test$fam_lev==1)
        fp<-sum(pred_test==1)/sum(rf_test$fam_lev==0)
        acc<-sum(pred_test==0)/length(pred_test)
        tss<-(1-fp)+(1-fn)-1
        err<-rbind(err,c(tre,fp,fn,tss,acc))
      }
      best_tre<-err[which(err[,4]==max(err[,4])),]
      #test spatial autocorrelation
      pred_data<-1*(predict(rf,rf_train)$predicted[,2]>best_tre[1])
      res<-pred_data-as.numeric(as.character(rf_train$fam_lev))
      if (sum(abs(res))>0){
         res<-as.factor(res)
         levels(res)<-c("-1","0","1")
         nn5 = knn2nb(knearneigh(pts[train_ids,],5))
         w = nb2listw(nn5, style="B")
         p<-joincount.test(res, w)[[1]]$p.value
         } else {p<-1}
      if (p>0.05){
        best_models[[length(best_models)+1]]<-list(
          'train_data'=train_ids,
          'p'=p,
          'tre'=best_tre[1],
          'false_pos'=best_tre[2],
          'false_neg'=best_tre[3],
          'tss'=best_tre[4],
          'acc' = best_tre[5],
          'loc_f' = loc_frac,
          'smote' = smote,
          'sample_f' = sample_frac
          )
        tss<-as.numeric(lapply(best_models,function(x) x$tss))
        print (paste
               (length(best_models),paste('fp=',round(best_tre[2],2),
                                          'fn=',round(best_tre[3],2),
                                          'tss=',round(best_tre[4],2),
                                          'acc=',round(best_tre[5],2)
                                          )
                 )
               )
        }
      else{
        n_spatial_bias<-n_spatial_bias+1
        print(paste('spatial autocorrelation',n_spatial_bias/n_tried_models))
        }   
    } 
  }



train_data_n<-as.numeric(lapply(best_models,function(x) length(x$train_data)))
p<-as.numeric(lapply(best_models,function(x) x$p))
tss<-as.numeric(lapply(best_models,function(x) x$tss))
fp<-as.numeric(lapply(best_models,function(x) x$false_pos))
fn<-as.numeric(lapply(best_models,function(x) x$false_neg))
tre<-as.numeric(lapply(best_models,function(x) x$tre))
acc<-as.numeric(lapply(best_models,function(x) x$acc))


sink('./RESULTS/cross_validation_stats.txt')
print('accuracy')
summary(acc)
sd(acc)
print('tss')
summary(tss)
sd(tss)
print('false positives')
summary(fp)
sd(fp)
print('false positives')
summary(fn)
sd(fn)
print('threshold')
summary(tre)
sd(tre)
sink()



train_ids<-1:nrow(aa) 
crisis_<-intersect(train_ids,crisis)
not_crisis_<-intersect(train_ids,not_crisis)
samp_size<-length(crisis_)
train_ids<-c(sample(crisis_,samp_size),sample(not_crisis_,samp_size*smote)) 
train_data<-aa[train_ids,]
rf_train<-train_data[,c(dep_var,t_cols,p_cols,add_var_ids)]
colnames(rf_train)[1]<-'fam_lev'
final_rf <-rfsrc(as.factor(fam_lev)~.,data=rf_train,ntree=1000)


####future projections
scenarios<-c('cur','ssp1_rcp26','ssp2_rcp45','ssp3_rcp7','ssp4_rcp6')
preds<-list()
for (sce in scenarios){
      pred<-c()
      pop_y<-c()
      gdp_y<-c()
      temp<-read.csv(paste0('./OUTPUT/future_data_for_projections/future_temp_',sce,'.csv'),header=F)
      prec<-read.csv(paste0('./OUTPUT/future_data_for_projections/future_prec_',sce,'.csv'),header=F)
      pop_dens<-read.csv(paste0('./OUTPUT/future_data_for_projections/future_pop_',sce,'.csv'),header=F)
      gdp<-read.csv(paste0('./OUTPUT/future_data_for_projections/future_gdp_',sce,'.csv'),header=F)
      for (i in seq(25,nrow(gdp),12)){
        te<-temp[(i-24):(i-1),]
        pe<-prec[(i-24):(i-1),]
        sum_t<-colSums(te)
        sum_p<-colSums(pe)
        sd_t<-apply(te,2,sd)
        sd_p<-apply(pe,2,sd)
        max_t<-apply(te,2,max)
        max_p<-apply(pe,2,max)
        min_t<-apply(te,2,min)
        min_p<-apply(pe,2,min)
        range_t<-max_t-min_t
        range_p<-max_p-min_p
        test_data<-cbind(t(te),t(pe),sum_t,sum_p,sd_t,sd_p,max_t,max_p,min_t,min_p,range_t,range_p)
        colnames(test_data)<-colnames(aa)[c(t_cols,p_cols,add_var_ids)]
        test_data <- rbind(final_rf$xvar[1, ] , test_data) #add and remove one line from the original data to ensure correct format
        test_data <- test_data[-1,]
        pred<-rbind(pred,predict(final_rf,test_data)$predicted[,2])
        pop_y<-rbind(pop_y,pop_dens[i,])
        gdp_y<-rbind(gdp_y,gdp[i,])
        print (c(sce,i))  
      }
      preds[[length(preds)+1]]<-list('pred'=pred,'pop'=pop_y,'gdp'=gdp_y)
    }




tre<-median(as.numeric(lapply(best_models,function(x) x$tre)))
gni_tre<-727.48####converting 1135 (2023 threshold) to 2005 value https://blogs.worldbank.org/opendata/new-world-bank-group-country-classifications-income-level-fy24
###inflation conversion https://www.officialdata.org/us/inflation/2023?endYear=2005&amount=1135

pop_cur<-preds[[1]]$pop
pred_cur<-preds[[1]]$pred
pred_cur[pop_cur<=0]<-NA
gdp_cur<-preds[[1]]$gdp
h_inc<-which(((gdp_b*gdp_cur/pop_cur)+gdp_a)>gni_tre)
pred_cur_no_gdp<-pred_cur
pred_cur[h_inc]<-0



pop_26<-preds[[2]]$pop
pred_26<-preds[[2]]$pred
pred_26[pop_26<=0]<-NA
gdp_26<-preds[[2]]$gdp
h_inc<-which(((gdp_b*gdp_26/pop_26)+gdp_a)>gni_tre)
pred_26_no_gdp<-pred_26
pred_26[h_inc]<-0



pop_45<-preds[[3]]$pop
pred_45<-preds[[3]]$pred
pred_45[pop_45<=0]<-NA
gdp_45<-preds[[3]]$gdp
h_inc<-which(((gdp_b*gdp_45/pop_45)+gdp_a)>gni_tre)
pred_45_no_gdp<-pred_45
pred_45[h_inc]<-0

pop_7<-preds[[4]]$pop
pred_7<-preds[[4]]$pred
pred_7[pop_26<=0]<-NA
gdp_7<-preds[[4]]$gdp
h_inc<-which(((gdp_b*gdp_7/pop_7)+gdp_a)>gni_tre)
pred_7_no_gdp<-pred_7
pred_7[h_inc]<-0

pop_6<-preds[[5]]$pop
pred_6<-preds[[5]]$pred
pred_6[pop_6<=0]<-NA
gdp_6<-preds[[5]]$gdp
h_inc<-which(((gdp_b*gdp_6/pop_6)+gdp_a)>gni_tre)
pred_6_no_gdp<-pred_6
pred_6[h_inc]<-0


tre_ll<-tre+tre*0.10
tre_ul<-tre-tre*0.10
pred_cur_ll<-1*(pred_cur>tre_ll)
pred_26_ll<-1*(pred_26>tre_ll)
pred_45_ll<-1*(pred_45>tre_ll)
pred_7_ll<-1*(pred_7>tre_ll)
pred_6_ll<-1*(pred_6>tre_ll)

pred_cur_ul<-1*(pred_cur>tre_ul)
pred_26_ul<-1*(pred_26>tre_ul)
pred_45_ul<-1*(pred_45>tre_ul)
pred_7_ul<-1*(pred_7>tre_ul)
pred_6_ul<-1*(pred_6>tre_ul)


pred_cur<-1*(pred_cur>tre)
pred_26<-1*(pred_26>tre)
pred_45<-1*(pred_45>tre)
pred_7<-1*(pred_7>tre)
pred_6<-1*(pred_6>tre)


pred_cur_no_gdp_ll<-1*(pred_cur_no_gdp>tre_ll)
pred_26_no_gdp_ll<-1*(pred_26_no_gdp>tre_ll)
pred_45_no_gdp_ll<-1*(pred_45_no_gdp>tre_ll)
pred_7_no_gdp_ll<-1*(pred_7_no_gdp>tre_ll)
pred_6_no_gdp_ll<-1*(pred_6_no_gdp>tre_ll)

pred_cur_no_gdp_ul<-1*(pred_cur_no_gdp>tre_ul)
pred_26_no_gdp_ul<-1*(pred_26_no_gdp>tre_ul)
pred_45_no_gdp_ul<-1*(pred_45_no_gdp>tre_ul)
pred_7_no_gdp_ul<-1*(pred_7_no_gdp>tre_ul)
pred_6_no_gdp_ul<-1*(pred_6_no_gdp>tre_ul)

pred_cur_no_gdp<-1*(pred_cur_no_gdp>tre)
pred_26_no_gdp<-1*(pred_26_no_gdp>tre)
pred_45_no_gdp<-1*(pred_45_no_gdp>tre)
pred_7_no_gdp<-1*(pred_7_no_gdp>tre)
pred_6_no_gdp<-1*(pred_6_no_gdp>tre)




#export for maps
data_for_maps<-rbind(
  cbind(0,2002:(2002+nrow(pred_cur)-1),as.matrix(pred_cur*pop_cur)),
  cbind(26,2017:(2017+nrow(pred_26)-1),as.matrix(pred_26*pop_26)),
  cbind(45,2017:(2017+nrow(pred_45)-1),as.matrix(pred_45*pop_45)),
  cbind(7,2017:(2017+nrow(pred_7)-1),as.matrix(pred_7*pop_7)),
  cbind(6,2017:(2017+nrow(pred_6)-1),as.matrix(pred_6*pop_6))
)


write.table(data_for_maps,'./RESULTS/data_4_maps.csv',row.names=F,col.names = F,
            sep=',',quote = F)





my_pal<-plasma(5)[1:4]
my_pal_a<-plasma(5,alpha=0.3)[1:4]
fac<-1

lat_lon<-read.csv('./OUTPUT/lat_lon.csv',header=T)
area<-lat_lon$area
tot_area<-sum(area)

pdf('./FIGURES/mean_and_total_fut_risk_low_income.pdf',width=9,height=4)
par(mfrow=c(1,3))
all_data<-rowSums(t(t(rbind(
  pred_cur_no_gdp,
  pred_26_no_gdp,
  pred_45_no_gdp,
  pred_7_no_gdp,
  pred_6_no_gdp,
  pred_cur_no_gdp_ll,
  pred_26_no_gdp_ll,
  pred_45_no_gdp_ll,
  pred_7_no_gdp_ll,
  pred_6_no_gdp_ll,
  pred_cur_no_gdp_ul,
  pred_26_no_gdp_ul,
  pred_45_no_gdp_ul,
  pred_7_no_gdp_ul,
  pred_6_no_gdp_ul))*100*area/tot_area),na.rm = T)


ll<-min(all_data)
ul<-max(all_data)


plot(0,0,type='n',
     xlab='year',
     ylab='potential climate-driven\n famine crisis risk (% of land area)',
     xlim=c(2000,2100),
     ylim=c(ll,ul),
     cex.lab=1.3,cex.axis=1.2,las=1)


col_n<-1
count_val<-function(x){sum(is.finite(x))}
for (pred_data in list(
  list(pred_26_no_gdp,pred_26_no_gdp_ll,pred_26_no_gdp_ul),
  list(pred_45_no_gdp,pred_45_no_gdp_ll,pred_45_no_gdp_ul),
  list(pred_7_no_gdp,pred_7_no_gdp_ll,pred_7_no_gdp_ul),
  list(pred_6_no_gdp,pred_6_no_gdp_ll,pred_6_no_gdp_ul))){
  
  pred_data_y<-pred_data[[1]][1:nrow(pred_data[[1]]),]
  pred_data_y<-t(t(pred_data_y)*area)
  start_y<-2017
  x<-start_y:(start_y+nrow(pred_data_y)-1)
  y<-100*rowSums(pred_data_y,na.rm = T)/tot_area
  
  pred_data_ll<-pred_data[[2]][1:nrow(pred_data[[2]]),]
  pred_data_ll<-t(t(pred_data_ll)*area)
  y_ll<-100*rowSums(pred_data_ll,na.rm = T)/tot_area
  
  pred_data_ul<-pred_data[[3]][1:nrow(pred_data[[3]]),]
  pred_data_ul<-t(t(pred_data_ul)*area)
  y_ul<-100*rowSums(pred_data_ul,na.rm = T)/tot_area
  
  polygon(c(rev(x), x), c(rev(y_ll),(y_ul)),
          col = my_pal_a[col_n], border = NA)
  
  lines(x,y,lwd=2,col=my_pal[col_n])
  col_n<-col_n+1
  }


pred_data<-list(pred_cur_no_gdp,pred_cur_no_gdp_ll,pred_cur_no_gdp_ul)

start_y<-2002
pred_data_y<-as.matrix(pred_data[[1]])
pred_data_y<-t(t(pred_data_y)*area)
x<-start_y:(start_y+nrow(pred_data_y)-1)
y<-100*rowSums(pred_data_y,na.rm = T)/tot_area


pred_data_ll<-pred_data[[2]][1:nrow(pred_data[[2]]),]
pred_data_ll<-t(t(pred_data_ll)*area)
y_ll<-100*rowSums(pred_data_ll,na.rm = T)/tot_area

pred_data_ul<-pred_data[[3]][1:nrow(pred_data[[3]]),]
pred_data_ul<-t(t(pred_data_ul)*area)
y_ul<-100*rowSums(pred_data_ul,na.rm = T)/tot_area

polygon(c(rev(x), x), c(rev(y_ll),(y_ul)),
        col = '#00000020', border = NA)
lines(x,y,col='black',lwd=1.5)

legend(2020,ul,c('historical','SSP1-26','SSP2-45','SSP3-7','SSP4-6'),col=c('black',my_pal),lwd=4,cex=1)



###exposed populations
all_data<-rowSums(rbind(
  pred_cur_ll*pop_cur,
  pred_26_ll*pop_26,
  pred_45_ll*pop_45,
  pred_7_ll*pop_7,
  pred_6_ll*pop_6,
  pred_cur_ul*pop_cur,
  pred_26_ul*pop_26,
  pred_45_ul*pop_45,
  pred_7_ul*pop_7,
  pred_6_ul*pop_6
  ),na.rm = T)


ll<-min(all_data)/10**6
ul<-max(all_data)/10**6


plot(0,0,type='n',
     xlab='year',
     ylab='potential climate-driven\n famine crisis risk (% of land area)',
     xlim=c(2000,2100),
     ylim=c(ll,ul),
     cex.lab=1.3,cex.axis=1.2,las=1)


col_n<-1
count_val<-function(x){sum(is.finite(x))}
for (pred_data in list(
  list(pred_26*pop_26,pred_26_ll*pop_26,pred_26_ul*pop_26),
  list(pred_45*pop_45,pred_45_ll*pop_45,pred_45_ul*pop_45),
  list(pred_7*pop_7,pred_7_ll*pop_7,pred_7_ul*pop_7),
  list(pred_6*pop_6,pred_6_ll*pop_6,pred_6_ul*pop_6))){
  
  start_y<-2017
  pred_data_y<-pred_data[[1]][1:nrow(pred_data[[1]]),]
  x<-start_y:(start_y+nrow(pred_data_y)-1)
  y<-rowSums(pred_data_y,na.rm = T)/10**6
  
  pred_data_ll<-pred_data[[2]][1:nrow(pred_data[[2]]),]
  y_ll<-rowSums(pred_data_ll,na.rm = T)/10**6
  
  pred_data_ul<-pred_data[[3]][1:nrow(pred_data[[3]]),]
  y_ul<-rowSums(pred_data_ul,na.rm = T)/10**6
  
  polygon(c(rev(x), x), c(rev(y_ll),(y_ul)),
          col = my_pal_a[col_n], border = NA)
  
  lines(x,y,lwd=2,col=my_pal[col_n])
  col_n<-col_n+1
}


pred_data<-list(pred_cur*pop_cur,pred_cur_ll*pop_cur,pred_cur_ul*pop_cur)

start_y<-2002
pred_data_y<-as.matrix(pred_data[[1]])
x<-start_y:(start_y+nrow(pred_data_y)-1)
y<-rowSums(pred_data_y,na.rm = T)/10**6


pred_data_ll<-pred_data[[2]][1:nrow(pred_data[[2]]),]
y_ll<-rowSums(pred_data_ll,na.rm = T)/10**6

pred_data_ul<-pred_data[[3]][1:nrow(pred_data[[3]]),]
y_ul<-rowSums(pred_data_ul,na.rm = T)/10**6

polygon(c(rev(x), x), c(rev(y_ll),(y_ul)),
        col = '#00000020', border = NA)
lines(x,y,col='black',lwd=1.5)

legend(2020,ul,c('historical','SSP1-26','SSP2-45','SSP3-7','SSP4-6'),col=c('black',my_pal),lwd=4,cex=1)



###exposed population %

all_data<-rowSums(rbind(
  pred_cur_ll*pop_cur,
  pred_26_ll*pop_26,
  pred_45_ll*pop_45,
  pred_7_ll*pop_7,
  pred_6_ll*pop_6,
  pred_cur_ul*pop_cur,
  pred_26_ul*pop_26,
  pred_45_ul*pop_45,
  pred_7_ul*pop_7,
  pred_6_ul*pop_6
),na.rm = T)

tot_pop<-rowSums(rbind(pop_cur,pop_26,pop_45,pop_7,pop_6),na.rm = T)
all_data<-100*(all_data/tot_pop)


ll<-min(all_data)
ul<-max(all_data)


plot(0,0,type='n',
     xlab='year',
     ylab='potential climate-driven\n famine crisis risk (% of land area)',
     xlim=c(2000,2100),
     ylim=c(ll,ul),
     cex.lab=1.3,cex.axis=1.2,las=1)


col_n<-1
count_val<-function(x){sum(is.finite(x))}
for (pred_data in list(
  list(pred_26*pop_26,pred_26_ll*pop_26,pred_26_ul*pop_26,pop_26),
  list(pred_45*pop_45,pred_45_ll*pop_45,pred_45_ul*pop_45,pop_45),
  list(pred_7*pop_7,pred_7_ll*pop_7,pred_7_ul*pop_7,pop_7),
  list(pred_6*pop_6,pred_6_ll*pop_6,pred_6_ul*pop_6,pop_6))){
  
  start_y<-2017
  pred_data_y<-pred_data[[1]][1:nrow(pred_data[[1]]),]
  x<-start_y:(start_y+nrow(pred_data_y)-1)
  y<-100*rowSums(pred_data_y,na.rm = T)/rowSums(pred_data[[4]],na.rm = T)
  
  pred_data_ll<-pred_data[[2]][1:nrow(pred_data[[2]]),]
  y_ll<-100*rowSums(pred_data_ll,na.rm = T)/rowSums(pred_data[[4]],na.rm = T)
  
  pred_data_ul<-pred_data[[3]][1:nrow(pred_data[[3]]),]
  y_ul<-100*rowSums(pred_data_ul,na.rm = T)/rowSums(pred_data[[4]],na.rm = T)
  
  polygon(c(rev(x), x), c(rev(y_ll),(y_ul)),
          col = my_pal_a[col_n], border = NA)
  
  lines(x,y,lwd=2,col=my_pal[col_n])
  col_n<-col_n+1
}


pred_data<-list(pred_cur*pop_cur,pred_cur_ll*pop_cur,pred_cur_ul*pop_cur,pop_cur)

start_y<-2002
pred_data_y<-as.matrix(pred_data[[1]])
x<-start_y:(start_y+nrow(pred_data_y)-1)
y<-100*rowSums(pred_data_y,na.rm = T)/rowSums(pred_data[[4]],na.rm = T)


pred_data_ll<-pred_data[[2]][1:nrow(pred_data[[2]]),]
y_ll<-100*rowSums(pred_data_ll,na.rm = T)/rowSums(pred_data[[4]],na.rm = T)

pred_data_ul<-pred_data[[3]][1:nrow(pred_data[[3]]),]
y_ul<-100*rowSums(pred_data_ul,na.rm = T)/rowSums(pred_data[[4]],na.rm = T)

polygon(c(rev(x), x), c(rev(y_ll),(y_ul)),
        col = '#00000020', border = NA)
lines(x,y,col='black',lwd=1.5)

legend(2020,ul,c('historical','SSP1-26','SSP2-45','SSP3-7','SSP4-6'),col=c('black',my_pal),lwd=4,cex=1)



dev.off()


####store summary data


y<-2017:2099

table_1<-rbind(
  cbind(26,y,rowSums(pred_26*pop_26,na.rm=T)/10**6,rowSums(pred_26_ll*pop_26,na.rm=T)/10**6,rowSums(pred_26_ul*pop_26,na.rm=T)/10**6),
  cbind(45,y,rowSums(pred_45*pop_45,na.rm=T)/10**6,rowSums(pred_45_ll*pop_45,na.rm=T)/10**6,rowSums(pred_45_ul*pop_45,na.rm=T)/10**6),
  cbind(7,y,rowSums(pred_7*pop_7,na.rm=T)/10**6,rowSums(pred_7_ll*pop_7,na.rm=T)/10**6,rowSums(pred_7_ul*pop_7,na.rm=T)/10**6),
  cbind(6,y,rowSums(pred_6*pop_6,na.rm=T)/10**6,rowSums(pred_6_ll*pop_6,na.rm=T)/10**6,rowSums(pred_6_ul*pop_6,na.rm=T)/10**6),
  cbind(0,2002:2022,rowSums(pred_cur*pop_cur,na.rm=T)/10**6,rowSums(pred_cur_ll*pop_cur,na.rm=T)/10**6,rowSums(pred_cur_ul*pop_cur,na.rm=T)/10**6))


colnames(table_1)<-c('ssp','year','exposed','exposed_ll','exposed_ul')
write.table(table_1,"./RESULTS/Table1.csv",sep=',',quote=F,row.names = F,
            col.names=T)




####CONTINENT
lat_lon<-read.csv('./OUTPUT/lat_lon.csv',header=T)
lat_lon<-cbind(lat_lon,'id'=1:nrow(lat_lon))
sf::sf_use_s2(FALSE)
sfRegion <- st_as_sf(lat_lon, coords = c("longitude","latitude"),crs=4326)
sfCountry <- ne_countries(returnclass='sf')
sfCountry<-sfCountry[sfCountry$iso_a3!='-99',]
sfJoined <- st_join(sfCountry, sfRegion)


matched_countries<-data.frame('id'=sfJoined$id,'continent'=sfJoined$continent,'country'=sfJoined$name_long,'iso3'=sfJoined$iso_a3)
matched_countries<-matched_countries[which(!is.na(matched_countries$id)),]

countries<-unique(matched_countries$country)
continents<-unique(matched_countries$continent)
continents<-c("Africa","Asia","Europe","North America","Oceania","South America")

matched_countries<-matched_countries[which(matched_countries$continent %in% continents),]

matched_countries$continent<-as.factor(matched_countries$continent)
area<-lat_lon$area


cont_fut_means<-c()
sce_n<-1

for (pred_data_n in 2:5){
  pred_data<-1*(preds[[pred_data_n]]$pred>tre)
  pred_data_ll<-1*(preds[[pred_data_n]]$pred>tre_ll)
  pred_data_ul<-1*(preds[[pred_data_n]]$pred>tre_ul)
  pop_data<-as.matrix(preds[[pred_data_n]]$pop)
  gdp<-preds[[pred_data_n]]$gdp
  h_inc<-which(((gdp_b*gdp/pop_data)+gdp_a)>gni_tre)
  pred_data[h_inc]<-NA
  pred_data<-as.matrix(pred_data)*as.matrix(pop_data)
  
  pred_data_ll[h_inc]<-NA
  pred_data_ll<-as.matrix(pred_data_ll)*as.matrix(pop_data)
  pred_data_ul[h_inc]<-NA
  pred_data_ul<-as.matrix(pred_data_ul)*as.matrix(pop_data)
  
  for (year in 1:nrow(pred_data)){
    tot_cont_pred<-aggregate(pred_data[year,matched_countries$id],list(matched_countries$continent),FUN='sum',na.rm=T)
    tot_cont_pred_ll<-aggregate(pred_data_ll[year,matched_countries$id],list(matched_countries$continent),FUN='sum',na.rm=T)[,2]
    tot_cont_pred_ul<-aggregate(pred_data_ul[year,matched_countries$id],list(matched_countries$continent),FUN='sum',na.rm=T)[,2]
    tot_cont_pop<-aggregate(pop_data[year,matched_countries$id],list(matched_countries$continent),FUN='sum',na.rm=T)[,2]
    cont_fut_means<-rbind(cont_fut_means,cbind(tot_cont_pred,tot_cont_pred_ll,tot_cont_pred_ul,tot_cont_pop,year,sce_n)
      )
    } 
    sce_n<-sce_n+1
  }



colnames(cont_fut_means)<-c('continent','tot_exposure','tot_exposure_ll','tot_exposure_ul','tot_population','year','scenario')
cont_fut_means<-data.frame(cont_fut_means)
head(cont_fut_means)
####


my_pal<-brewer.pal(n = length(continents), name = "Set2")
rcp_scenarios<-c('SSP1-26','SSP2-45','SSP3-7','SSP4-6')

options(scipen=500)
pdf('./FIGURES/continent_tot_exposure.pdf',width=8,height=8)
par(mfrow=c(2,2))

###first the actual data, in M people exposed
for (rcp_scen in 1:4){
  data<-cont_fut_means[cont_fut_means$scenario==rcp_scen,]

  data_for_l<-c()
  for (cont in continents){
    cont_data<-data[data$continent==cont,]
    y_ll<-cont_data$tot_exposure_ll
    y_ul<-cont_data$tot_exposure_ul
    data_for_l<-c(data_for_l,y_ll,y_ul)
    }
  
  ll<-min(data_for_l,na.rm = T)/1000000
  ul<-max(data_for_l,na.rm = T)/1000000
  plot(1,1,type='n',
       xlab='year',
       ylab='exposed individuals (M)',
       xlim=c(2010,2100),
       ylim=c(ll,ul),
       cex.lab=1.3,cex.axis=1.2,las=1,cex.main=1.5,
       main=rcp_scenarios[rcp_scen])
  
  col_n<-1
  for (cont in continents){
    cont_data<-data[data$continent==cont,]
    x<-seq(2017,2099)
    y<-cont_data$tot_exposure/1000000
    y_ll<-cont_data$tot_exposure_ll/10**6
    y_ul<-cont_data$tot_exposure_ul/10**6
    
    polygon(c(rev(x), x), c(rev(y_ll),(y_ul)),
            col = paste0(my_pal[col_n],70), border = NA)
    lines(x,y,col=my_pal[col_n],lwd=1.5)
    col_n<-col_n+1
  }
  
  if (rcp_scen==1){
  legend(2010,ul,continents,col=my_pal,lwd=4,cex=1)
   } 
}


dev.off()




#export data to supplementary table
to_export<-data.frame('ssp'=cont_fut_means$scenario,
                      'year'=cont_fut_means$year,
                      'continent'=cont_fut_means$continent,
                      'exposed'=cont_fut_means$tot_exposure,
                      'exposed_ll'=cont_fut_means$tot_exposure_ll,
                      'exposed_ul'=cont_fut_means$tot_exposure_ul
                      )


to_export$year<-to_export$year+2016
to_export$ssp<-rcp_scenarios[to_export$ssp]
to_export$exposed<-round(to_export$exposed)
head(to_export)

write.table(to_export,"./RESULTS/TableS1.csv",sep=',',quote=F,row.names = F,
            col.names=T)


tot_cont_ssp<-aggregate(to_export$exposed~to_export$ssp+to_export$continent,FUN='sum')
cont_perc_exp<-c()
for (sce in unique(tot_cont_ssp$`to_export$ssp`)){
  data<-tot_cont_ssp[which(tot_cont_ssp$`to_export$ssp`==sce),]
  data$`to_export$exposed`<-round(100*data$`to_export$exposed`/sum(data$`to_export$exposed`),1) 
  cont_perc_exp<-rbind(cont_perc_exp,data)
}

colnames(cont_perc_exp)<-c('ssp','continent','exposed_perc')

write.table(cont_perc_exp,'./RESULTS/exposure_by_continent.csv',col.names=T,row.names=F,quote=F,sep=',')



####population estimates https://population.un.org/wpp/Download/Standard/CSV/
###https://population.un.org/wpp/DefinitionOfProjectionScenarios/

est_pop<-read.csv('./INPUT/pop_density/WPP2022_Demographic_Indicators_OtherVariants.csv',header=T)
est_pop<-est_pop[est_pop$Location=='World',]
variants<-unique(est_pop$Variant)
pops_v<-c()
for (v in variants){
  pop_y<-est_pop[est_pop$Variant==v,]
  pop_y<-pop_y[pop_y$Time>2016,]
  pop_y<-pop_y[pop_y$Time<2100,]$TPopulation1Jan 
  pops_v<-cbind(pops_v,pop_y)
}



###add the medium scenario

pop_med<-read.csv('./INPUT/pop_density/WPP2022_Demographic_Indicators_Medium.csv',header=T)


pop_med<-pop_med[pop_med$Location=='World',]
pop_med<-pop_med[pop_med$Time>2021,]
pop_med<-pop_med[pop_med$Time<2100,]$TPopulation1Jan 
pops_v<-cbind(pops_v,pop_med)

variants<-c(variants,'Medium')

best_pop_mod<-c()
for (pop_ssp in list(pop_26,pop_45,pop_7,pop_6)){
  y<-rowSums(pop_ssp,na.rm=T)[6:nrow(pop_26)]
  pop_cor<-c()
  for (i in 1:ncol(pops_v)){
    pop_cor<-c(pop_cor,cor(pops_v[,i],y/1000))
    }
  best_pop_mod<-c(best_pop_mod,variants[which(pop_cor==max(pop_cor))])
}

best_pop_mod[2]<-'Medium' ###identical to zero migration, but makes more sense



###now take birth and mortality rate
######OTHER APPROACH, COMPUTE AGE DISTRIBUTION

all_countries<-unique(matched_countries$iso3)

est_pop_age<-read.csv('./INPUT/pop_density/WPP2022_PopulationByAge5GroupSex_OtherVariants.csv',header=T)
medium_age<-read.csv('./INPUT/pop_density/WPP2022_PopulationByAge5GroupSex_Medium.csv',header=T)
est_pop_age<-rbind(est_pop_age,medium_age)
est_pop<-est_pop[est_pop$ISO3_code %in% all_countries,]
est_pop_age<-est_pop_age[est_pop_age$Time>2021,]
est_pop_age<-est_pop_age[est_pop_age$Time<2100,]

#citation("binsmooth")
age_pop<-c()
all_y<-unique(est_pop_age$Time)
for (mod in best_pop_mod){
  for (loc in all_countries){
    for (y in all_y){
      sel_rows<-intersect(intersect(which(est_pop_age$Variant==mod),which(est_pop_age$ISO3_code==loc)),which(est_pop_age$Time==y))
      sel_data<-est_pop_age[sel_rows,]      
      n_ac<-c(rep(sel_data$PopTotal[1],4)/5,sel_data$PopTotal[2:length(sel_data$PopTotal)],0)
      n_right<-c(1:4,seq(9,104,5),NA)
      mean_age<-weighted.mean(n_right, n_ac,na.rm=T)
      s<-splinebins(n_right,n_ac,mean_age)
    all_age_classes<-round(s$splinePDF(0:100)*sum(n_ac),3)
    all_age_classes<-(all_age_classes/sum(all_age_classes))*sum(n_ac)
    age_pop<-rbind(age_pop,c(mod,loc,y,all_age_classes))
    print(c(mod,loc,y))
    }
  }
}

####

age_pop<-data.frame(age_pop)

colnames(age_pop)<-c('Variant','ISO3_code','Time',0:100)
age_pop_c<-age_pop

for (i in  4:ncol(age_pop)){
  age_pop[,i]<-as.numeric(age_pop[,i])
}

tot_rows<-rowSums(age_pop[,4:ncol(age_pop)],na.rm=T)[which(rowSums(age_pop[,4:ncol(age_pop)],na.rm=T)!=1)]
age_pop[,4:ncol(age_pop)]<-age_pop[,4:ncol(age_pop)]/tot_rows


best_pop_mod[4]<-'High'
all_countries<-unique(matched_countries$iso3)

pred_datas<-list(pred_26,pred_45,pred_7,pred_6)
pred_datas_ll<-list(pred_26_ll,pred_45_ll,pred_7_ll,pred_6_ll)
pred_datas_ul<-list(pred_26_ul,pred_45_ul,pred_7_ul,pred_6_ul)
pop_datas<-list(pop_26,pop_45,pop_7,pop_6)

pred_datas_ci<-list(pred_datas,pred_datas_ll,pred_datas_ul)

#####first for all ages

max_age<-100
accum_exp_ci<-list()
for (pred_data_ci in pred_datas_ci){
  accum_exp<-list()
  for (v in 1:4){
    pred_data<-as.matrix(pred_data_ci[[v]])
    pred_data<-pred_data[6:nrow(pred_data),matched_countries$id]
    with_crisis<-which(colSums(pred_data,na.rm=T)>0)
    pop_data<-as.matrix(pop_datas[[v]])
    pop_data<-pop_data[6:nrow(pop_data),matched_countries$id]
    pop_y<-age_pop[age_pop$Variant==best_pop_mod[v],] #[1] "Low"         "Medium"      "High"        "Upper 95 PI"
    b_table<-matrix(0,nrow=78,ncol=nrow(matched_countries))    
    for (country in all_countries){
      target_cols<-intersect(which(matched_countries$iso3==country),with_crisis)
      target_rows<-which(pop_y$ISO3_code==country)
      pop_age_loc<-pop_y[pop_y$ISO3_code==country,4:ncol(pop_y)]
    ###normalize pop_age_loc
      for (col in target_cols){
        cris_series<-pred_data[,col]
        n<-max_age
        exp_pop<-c()
        for (row in 1:nrow(pred_data)){
          n<-n+1
          if (cris_series[row]==0){
            exp_pop<-c(exp_pop,0)
            } 
          else {
                if (n>max_age){
                   n<-max_age
                  }
              exp_pop<-c(exp_pop,sum(as.numeric(pop_age_loc[row,1:(n+1)]*pop_data[row,col]),na.rm=T))
              n<-0
            }        
          }
        b_table[,col]<-exp_pop
      }
    }  
    y<-cumsum(rowSums(b_table,na.rm=T))/1000000000
    accum_exp[[length(accum_exp)+1]]<-y
    print (v)
  }
  accum_exp_ci[[length(accum_exp_ci)+1]]<-accum_exp
}


all_data<-unlist(accum_exp_ci)
y_min<-0#min(all_data)
y_max<-max(all_data)
x<-2022:2099
id_2050<-which(x==2050)

pdf('./FIGURES/cumulative_exposed.pdf',width=4.5,height=4.5)
plot(1,1,type='n',ylim=c(y_min,y_max),xlim=c(2020,2100),
     ylab='cumulative exposed people (B)',
     xlab='year',las=1,cex.axis=1.2,cex.lab=1.2,
     #     log='y'
)

for (i in 1:4){
  lines(x,accum_exp_ci[[1]][[i]],lwd=2,col=plasma(5)[i])
  polygon(c(rev(x), x), c(rev(accum_exp_ci[[2]][[i]]),(accum_exp_ci[[3]][[i]])),
          col=plasma(5,alpha=0.2)[i], border = NA)
  segments(2000,accum_exp_ci[[1]][[i]][id_2050],2050,accum_exp_ci[[1]][[i]][id_2050],lwd=1,col=plasma(5)[i],lty=1)
  }

lines(x,accum_exp_ci[[1]][[4]]-accum_exp_ci[[1]][[1]],lwd=2,lty=2,col='black')

segments(2050,0,2050,accum_exp_ci[[1]][[4]][id_2050])
text(2050,0,'2050',cex=1)

legend(2022,1.3,c('SSP1-2.6','SSP2-4.5','SSP3-7.0','SSP4-6.0'),col=plasma(5),lwd=4,cex=0.8,bty='n')
legend(2022,1,c('excess SSP4-SSP1'),lwd=2,lty=2,cex=0.7,bty='n')

dev.off()



##########
#####then for children<5

max_age<-4

accum_exp_ci_ch<-list()
for (pred_data_ci in pred_datas_ci){
  accum_exp<-list()
  for (v in 1:4){
    pred_data<-as.matrix(pred_data_ci[[v]])
    pred_data<-pred_data[6:nrow(pred_data),matched_countries$id]
    with_crisis<-which(colSums(pred_data,na.rm=T)>0)
    pop_data<-as.matrix(pop_datas[[v]])
    pop_data<-pop_data[6:nrow(pop_data),matched_countries$id]
    pop_y<-age_pop[age_pop$Variant==best_pop_mod[v],] #[1] "Low"         "Medium"      "High"        "Upper 95 PI"
    b_table<-matrix(0,nrow=78,ncol=nrow(matched_countries))    
    for (country in all_countries){
      target_cols<-intersect(which(matched_countries$iso3==country),with_crisis)
      target_rows<-which(pop_y$ISO3_code==country)
      pop_age_loc<-pop_y[pop_y$ISO3_code==country,4:ncol(pop_y)]
      ###normalize pop_age_loc
      for (col in target_cols){
        cris_series<-pred_data[,col]
        n<-max_age
        exp_pop<-c()
        for (row in 1:nrow(pred_data)){
          n<-n+1
          if (cris_series[row]==0){
            exp_pop<-c(exp_pop,0)
          } 
          else {
            if (n>max_age){
              n<-max_age
            }
            exp_pop<-c(exp_pop,sum(as.numeric(pop_age_loc[row,1:(n+1)]*pop_data[row,col]),na.rm=T))
            n<-0
          }        
        }
        b_table[,col]<-exp_pop
      }
    }  
    y<-cumsum(rowSums(b_table,na.rm=T))/1000000000
    accum_exp[[length(accum_exp)+1]]<-y
    print (v)
  }
  accum_exp_ci_ch[[length(accum_exp_ci_ch)+1]]<-accum_exp
}


all_data<-unlist(accum_exp_ci_ch)
y_min<-0#min(all_data)
y_max<-max(all_data)
x<-2022:2099
id_2050<-which(x==2050)

pdf('./FIGURES/cumulative_exposed_children.pdf',width=4.5,height=4.5)
plot(1,1,type='n',ylim=c(y_min,y_max),xlim=c(2020,2100),
     ylab='cumulative exposed children (B)',
     xlab='year',las=1,cex.axis=1.2,cex.lab=1.2,
)


for (i in 1:4){
  lines(x,accum_exp_ci_ch[[1]][[i]],lwd=2,col=plasma(5)[i])
  polygon(c(rev(x), x), c(rev(accum_exp_ci_ch[[2]][[i]]),(accum_exp_ci_ch[[3]][[i]])),
          col=plasma(5,alpha=0.2)[i], border = NA)
  segments(2000,accum_exp_ci_ch[[1]][[i]][id_2050],2050,accum_exp_ci_ch[[1]][[i]][id_2050],lwd=1,col=plasma(5)[i],lty=1)
}

lines(x,accum_exp_ci_ch[[1]][[4]]-accum_exp_ci_ch[[1]][[1]],lwd=2,lty=2,col='black')

segments(2050,0,2050,accum_exp_ci_ch[[1]][[4]][id_2050])
text(2050,0,'2050',cex=1)

legend(2022,1.3,c('SSP1-2.6','SSP2-4.5','SSP3-7.0','SSP4-6.0'),col=plasma(5),lwd=4,cex=0.8,bty='n')
legend(2022,1,c('excess SSP4-SSP1'),lwd=2,lty=2,cex=0.7,bty='n')

dev.off()



#####then for newborns

max_age<-0

accum_exp_ci_nb<-list()
for (pred_data_ci in pred_datas_ci){
  accum_exp<-list()
  for (v in 1:4){
    pred_data<-as.matrix(pred_data_ci[[v]])
    pred_data<-pred_data[6:nrow(pred_data),matched_countries$id]
    with_crisis<-which(colSums(pred_data,na.rm=T)>0)
    pop_data<-as.matrix(pop_datas[[v]])
    pop_data<-pop_data[6:nrow(pop_data),matched_countries$id]
    pop_y<-age_pop[age_pop$Variant==best_pop_mod[v],] #[1] "Low"         "Medium"      "High"        "Upper 95 PI"
    b_table<-matrix(0,nrow=78,ncol=nrow(matched_countries))    
    for (country in all_countries){
      target_cols<-intersect(which(matched_countries$iso3==country),with_crisis)
      target_rows<-which(pop_y$ISO3_code==country)
      pop_age_loc<-pop_y[pop_y$ISO3_code==country,4:ncol(pop_y)]
      ###normalize pop_age_loc
      for (col in target_cols){
        cris_series<-pred_data[,col]
        n<-max_age
        exp_pop<-c()
        for (row in 1:nrow(pred_data)){
          n<-n+1
          if (cris_series[row]==0){
            exp_pop<-c(exp_pop,0)
          } 
          else {
            if (n>max_age){
              n<-max_age
            }
            exp_pop<-c(exp_pop,sum(as.numeric(pop_age_loc[row,1:(n+1)]*pop_data[row,col]),na.rm=T))
            n<-0
          }        
        }
        b_table[,col]<-exp_pop
      }
    }  
    y<-cumsum(rowSums(b_table,na.rm=T))/1000000000
    accum_exp[[length(accum_exp)+1]]<-y
    print (v)
  }
  accum_exp_ci_nb[[length(accum_exp_ci_nb)+1]]<-accum_exp
}


all_data<-unlist(accum_exp_ci_nb)
y_min<-0
y_max<-max(all_data)
x<-2022:2099
id_2050<-which(x==2050)


pdf('./FIGURES/cumulative_exposed_newborns.pdf',width=4.5,height=4.5)
plot(1,1,type='n',ylim=c(y_min,y_max),xlim=c(2020,2100),
     ylab='cumulative exposed newborn (B)',
     xlab='year',las=1,cex.axis=1.2,cex.lab=1.2,
)


for (i in 1:4){
  lines(x,accum_exp_ci_nb[[1]][[i]],lwd=2,col=plasma(5)[i])
  polygon(c(rev(x), x), c(rev(accum_exp_ci_nb[[2]][[i]]),(accum_exp_ci_nb[[3]][[i]])),
          col=plasma(5,alpha=0.2)[i], border = NA)
  segments(2000,accum_exp_ci_nb[[1]][[i]][id_2050],2050,accum_exp_ci_nb[[1]][[i]][id_2050],lwd=1,col=plasma(5)[i],lty=1)
}

lines(x,accum_exp_ci_nb[[1]][[4]]-accum_exp_ci_nb[[1]][[1]],lwd=2,lty=2,col='black')

segments(2050,0,2050,accum_exp_ci_nb[[1]][[4]][id_2050])
text(2050,0,'2050',cex=1)

legend(2022,1.3,c('SSP1-2.6','SSP2-4.5','SSP3-7.0','SSP4-6.0'),col=plasma(5),lwd=4,cex=0.8,bty='n')
legend(2022,1,c('excess SSP4-SSP1'),lwd=2,lty=2,cex=0.7,bty='n')

dev.off()







exp_age_data<-data.frame(cbind(2022:2099,accum_exp_ci[[1]][[1]],accum_exp_ci[[1]][[2]],accum_exp_ci[[1]][[3]],accum_exp_ci[[1]][[4]],
                               accum_exp_ci_ch[[1]][[1]],accum_exp_ci_ch[[1]][[2]],accum_exp_ci_ch[[1]][[3]],accum_exp_ci_ch[[1]][[4]],
                               accum_exp_ci_nb[[1]][[1]],accum_exp_ci_nb[[1]][[2]],accum_exp_ci_nb[[1]][[3]],accum_exp_ci_nb[[1]][[4]]))

colnames(exp_age_data)<-c('year','SSP1-2.6','SSP2-4.5','SSP3-7.0','SSP4-6.0',
                          paste('CH_',c('SSP1-2.6','SSP2-4.5','SSP3-7.0','SSP4-6.0')),
                          paste('NB_',c('SSP1-2.6','SSP2-4.5','SSP3-7.0','SSP4-6.0')))


write.table(exp_age_data,'./RESULTS/exp_age_data_summary.csv',
            row.names=F,col.names=T,sep=',',quote=F)



round(exp_age_data[exp_age_data$year==2099,],4)


##################SENSITIVITY GNI
start_y<-2017

sce_names<-c("cur","SSP1-2.6","SSP2-4.5","SSP3-7.0","SSP4-6.0") 
gni_sens<-c()
for (n in 2:5){
  for (gni_tre in seq(50,5000,50)){
    pop_data<-preds[[n]]$pop
    pred_data<-preds[[n]]$pred
    pred_data[pop_data<=0]<-NA
    gdp_data<-preds[[n]]$gdp
    h_inc<-which(((gdp_b*gdp_data/pop_data)+gdp_a)>gni_tre)
    pred_data[h_inc]<-0
    pred_data<-1*(pred_data>tre)*(pop_data)
    x<-start_y:(start_y+nrow(pred_data)-1)
    y<-rowSums(pred_data,na.rm = T)/10**6
    gni_sens<-rbind(gni_sens,cbind(n,gni_tre,x,y))
    print(c(n,gni_tre))
  }
}

head(gni_sens)
pdf('./FIGURES/sensitivity_gni.pdf',height=8,width=8)
par(mfrow=c(2,2))
for (n in 2:5){
  gni_sens_<-gni_sens[gni_sens[,1]==n,]
image.plot(interp(gni_sens_[,3],gni_sens_[,4]/1000,gni_sens_[,2],duplicate='mean'),col=viridis(30),
           xlab='year',ylab='exposed population (B)',main=sce_names[n],
           las=1,cex.axis=1.5,cex.lab=1.5,cex.main=2,legend.lab='GNI treshold (2005 $)',horizontal=T)
}
dev.off()




######SENSITIVITY ANALYSIS MONTHS BEFORE CRISIS
###make cross spatial validation by splitting the dataset by countries and year
aa<-read.csv('./OUTPUT/matched_climate_crisis.csv',header=T)
aa<-aa[is.finite(rowSums(aa[,6:ncol(aa)])),]
aa<-aa[aa$fam_lev!=-99999,]
aa<-aggregate(.~year+row+col,data=aa,FUN='max')
dep_var_name<-'fam_lev'#'crisis_famine_val'#
dep_var<-which(colnames(aa)==dep_var_name)
aa<-aa[which(aa[,dep_var]!=-99999),]
aa$id<-1:nrow(aa)
sf::sf_use_s2(FALSE)
sfRegion <- st_as_sf(aa, coords = c("col","row"),crs=4326)
sfCountry <- ne_countries(returnclass='sf')
sfJoined <- st_join(sfCountry, sfRegion)
sfJoined<-sfJoined[which(!is.na(sfJoined$id)),]
aa<-aa[sfJoined$id,]
aa$locs<-sfJoined$name_long
aa$continent<-sfJoined$continent
locs<-unique(aa$locs)
aa[,dep_var]<-as.factor(aa[,dep_var])
aa$continent<-as.factor(aa$continent)
crisis<-which(aa[,dep_var]==1)
not_crisis<-which(aa[,dep_var]==0)
years<-unique(aa$year)
null_vals<-which(aa[,dep_var]==-99999)
y_loc<-paste0(round(aa$row/1),'_',round(aa$col/1))#paste0(aa$year,aa$locs)
y_loc_u<-unique(y_loc)
DT <- data.frame(longitude=aa$col,
                 latitude=aa$row+runif(nrow(aa))/10000)
pts = st_as_sf(DT, coords = c("longitude", "latitude"), 
               crs = 4326, agr = "constant")


mod_sensitivity<-list()
for (month_n in c(11,35)){
  t_cols<-match(paste0('t',0:month_n),colnames(aa))
  p_cols<-match(paste0('p',0:month_n),colnames(aa))
  sum_t<-rowSums(aa[,t_cols])
  sum_p<-rowSums(aa[,p_cols])
  sd_t<-apply(aa[,t_cols],1,sd)
  sd_p<-apply(aa[,p_cols],1,sd)
  max_t<-apply(aa[,t_cols],1,max)
  max_p<-apply(aa[,p_cols],1,max)
  min_t<-apply(aa[,t_cols],1,min)
  min_p<-apply(aa[,p_cols],1,min)
  range_t<-max_t-min_t
  range_p<-max_p-min_p
  aaa<-cbind(aa,sum_t,sum_p,sd_t,sd_p,max_t,max_p,min_t,min_p,range_t,range_p)
  add_var_ids<-match(c("sum_t","sum_p","sd_t","sd_p","max_t","max_p","min_t","min_p","range_t","range_p"),colnames(aaa))
  models<-list()
  for (replicate in 1:10){
    for (loc_frac in seq(0.1,0.9,0.1)){
      for (smote in seq(0.5,3,0.5)){
        for (sample_frac in seq(0.5,0.75,0.1)){
          rand_locs<-sample(y_loc_u,length(y_loc_u)*loc_frac)
          train_ids<-which(!(y_loc%in%rand_locs))
          test_ids<-setdiff(1:nrow(aa),c(train_ids,null_vals))
          if (length(test_ids)>15000){
            test_ids<-sample(test_ids,length(test_ids))[1:15000]}
          crisis_<-intersect(train_ids,crisis)
          if (length(crisis_)>0){
            not_crisis_<-intersect(train_ids,not_crisis)
            samp_size<-length(crisis_)*sample_frac
            train_ids<-c(sample(crisis_,samp_size),sample(not_crisis_,samp_size*smote)) #SMOTE: Synthetic Minority Over-sampling Technique
            train_data<-aaa[train_ids,]
            test_data<-aaa[test_ids,]
            rf_train<-train_data[,c(dep_var,t_cols,p_cols,add_var_ids)]
            rf_test<-test_data[,c(dep_var,t_cols,p_cols,add_var_ids)]
            colnames(rf_train)[1]<-'fam_lev'
            colnames(rf_test)[1]<-'fam_lev'
            rf <-rfsrc(as.factor(fam_lev)~.,data=rf_train,ntree=250)
            pred_data<-predict(rf,rf_test)$predicted[,2]
            err<-c()
            for (tre in seq(0.3,0.8,0.1)){
              pred_test<-(1*(pred_data>tre))-as.numeric(as.character(rf_test$fam_lev))
              fn<-sum(pred_test==-1)/sum(rf_test$fam_lev==1)
              fp<-sum(pred_test==1)/sum(rf_test$fam_lev==0)
              tss<-(1-fp)+(1-fn)-1
              err<-rbind(err,c(tre,fp,fn,tss))
            }
            best_tre<-err[which(err[,4]==max(err[,4])),]
            models[[length(models)+1]]<-list('train_data'=train_ids,
                                             'tre'=best_tre[1],
                                             'false_pos'=best_tre[2],
                                             'false_neg'=best_tre[3],
                                             'tss'=best_tre[4],
                                             'loc_f' = loc_frac,
                                             'smote' = smote,
                                             'sample_f' = sample_frac
            )
            tss<-as.numeric(lapply(models,function(x) x$tss))
            print (paste(replicate,length(models),paste('fp=',round(best_tre[2],2),
                                                        'fn=',round(best_tre[3],2),
                                                        'tss=',round(best_tre[4],2),
                                                        'max_tss=',round(max(tss),2)
                                                        
            )))
            
          } 
          
        }
      }
    }
  }
  #identify the best parametrization, i.e. the one leading to the higher mean TSS
  tss<-as.numeric(lapply(models,function(x) x$tss))
  smote<-as.numeric(lapply(models,function(x) x$smote))
  loc_f<-as.numeric(lapply(models,function(x) x$loc_f))
  sample_f<-as.numeric(lapply(models,function(x) x$sample_f))
  eval_mod<-aggregate(tss~smote+loc_f+sample_f,FUN='mean')
  best_pars<-eval_mod[which(eval_mod$tss==max(eval_mod$tss)),]
  smote<-as.numeric(best_pars[1])
  loc_frac<-as.numeric(best_pars[2])
  sample_frac<-as.numeric(best_pars[3])
  ###now obtain 100 models with the best parametrization and no spatial autocorrelation
  best_models<-list()
  while (length(best_models)<100){
    rand_locs<-sample(y_loc_u,length(y_loc_u)*loc_frac)
    train_ids<-which(!(y_loc%in%rand_locs)) 
    test_ids<-setdiff(1:nrow(aa),c(train_ids,null_vals))
    if (length(test_ids)>15000){
      test_ids<-sample(test_ids,length(test_ids))[1:15000]}
    crisis_<-intersect(train_ids,crisis)
    if (length(crisis_)>0){
      not_crisis_<-intersect(train_ids,not_crisis)
      samp_size<-length(crisis_)*sample_frac
      train_ids<-c(sample(crisis_,samp_size),sample(not_crisis_,samp_size*smote)) #SMOTE: Synthetic Minority Over-sampling Technique
      train_data<-aaa[train_ids,]
      test_data<-aaa[test_ids,]
      rf_train<-train_data[,c(dep_var,t_cols,p_cols,add_var_ids)]
      rf_test<-test_data[,c(dep_var,t_cols,p_cols,add_var_ids)]
      colnames(rf_train)[1]<-'fam_lev'
      colnames(rf_test)[1]<-'fam_lev'
      rf <-rfsrc(as.factor(fam_lev)~.,data=rf_train,ntree=500)
      pred_data<-predict(rf,rf_test)$predicted[,2]
      err<-c()
      for (tre in seq(0.3,0.8,0.1)){
        pred_test<-(1*(pred_data>tre))-as.numeric(as.character(rf_test$fam_lev))
        fn<-sum(pred_test==-1)/sum(rf_test$fam_lev==1)
        fp<-sum(pred_test==1)/sum(rf_test$fam_lev==0)
        acc<-sum(pred_test==0)/length(pred_test)
        tss<-(1-fp)+(1-fn)-1
        err<-rbind(err,c(tre,fp,fn,tss,acc))
      }
      best_tre<-err[which(err[,4]==max(err[,4])),]
      #test spatial autocorrelation
      pred_data<-1*(predict(rf,rf_train)$predicted[,2]>best_tre[1])
      res<-pred_data-as.numeric(as.character(rf_train$fam_lev))
      if (sum(abs(res))>0){
        res<-as.factor(res)
        levels(res)<-c("-1","0","1")
        nn5 = knn2nb(knearneigh(pts[train_ids,],5))
        w = nb2listw(nn5, style="B")
        p<-joincount.test(res, w)[[1]]$p.value
      } else {p<-1}
      if (p>0.05){
        best_models[[length(best_models)+1]]<-list(
          'train_data'=train_ids,
          'p'=p,
          'tre'=best_tre[1],
          'false_pos'=best_tre[2],
          'false_neg'=best_tre[3],
          'tss'=best_tre[4],
          'acc' = best_tre[5],
          'loc_f' = loc_frac,
          'smote' = smote,
          'sample_f' = sample_frac
        )
        tss<-as.numeric(lapply(best_models,function(x) x$tss))
        print (paste
               (length(best_models),paste('fp=',round(best_tre[2],2),
                                          'fn=',round(best_tre[3],2),
                                          'tss=',round(best_tre[4],2),
                                          'acc=',round(best_tre[5],2)
               )
               )
        )
      }else{print('spatial autocorrelation')}   
    } 
  }
  
  
  train_data_n<-as.numeric(lapply(best_models,function(x) length(x$train_data)))
  p<-as.numeric(lapply(best_models,function(x) x$p))
  tss<-as.numeric(lapply(best_models,function(x) x$tss))
  fp<-as.numeric(lapply(best_models,function(x) x$false_pos))
  fn<-as.numeric(lapply(best_models,function(x) x$false_neg))
  tre<-as.numeric(lapply(best_models,function(x) x$tre))
  acc<-as.numeric(lapply(best_models,function(x) x$acc))
  mod_sensitivity[[as.character(month_n)]]<-list('tss'=tss,'acc'=acc,'fp'=fp,'fn'=fn)
}



sink("./RESULTS/sensitivity_analysis_months_before_crisis.txt")
print('with 12 months')
print('tss')
print(summary(mod_sensitivity[['11']]$tss))
print(sd(mod_sensitivity[['11']]$tss))

print('acc')
print(summary(mod_sensitivity[['11']]$acc))
print(sd(mod_sensitivity[['11']]$acc))

print('false pos')
print(summary(mod_sensitivity[['11']]$fp))
print(sd(mod_sensitivity[['11']]$fp))

print('false neg')
print(summary(mod_sensitivity[['11']]$fn))
print(sd(mod_sensitivity[['11']]$fn))


print('with 36 months')
print('tss')
print(summary(mod_sensitivity[['35']]$tss))
print(sd(mod_sensitivity[['35']]$tss))

print('acc')
print(summary(mod_sensitivity[['35']]$acc))
print(sd(mod_sensitivity[['35']]$acc))

print('false pos')
print(summary(mod_sensitivity[['35']]$fp))
print(sd(mod_sensitivity[['35']]$fp))

print('false neg')
print(summary(mod_sensitivity[['35']]$fn))
print(sd(mod_sensitivity[['35']]$fn))


sink()


#########SENSITIVITY ANALYSIS ADMIN AREAS
aa<-read.csv('./OUTPUT/matched_climate_crisis.csv',header=T)
aa<-aa[is.finite(rowSums(aa[,6:ncol(aa)])),]
aa<-aa[aa$fam_lev!=-99999,]
aa<-aggregate(.~year+row+col,data=aa,FUN='max')
dep_var_name<-'fam_lev'#'crisis_famine_val'#
dep_var<-which(colnames(aa)==dep_var_name)
aa<-aa[which(aa[,dep_var]!=-99999),]
aa$id<-1:nrow(aa)
sf::sf_use_s2(FALSE)

fews_admin_bounds<-read_sf(dsn = "./FEWS_NET_Admin_Boundaries", layer = "FEWS_Admin_LZ_v3")
fews_admin_bounds<-st_transform(fews_admin_bounds, 4326)
sfRegion <- st_as_sf(aa, coords = c("col","row"),crs=4326)

#join raster data to admin boundaries used by FEWS
admin_bound_join <- st_join(sfRegion,fews_admin_bounds)
sfCountry <- ne_countries(returnclass='sf')
sfJoined <- st_join(sfRegion,sfCountry)

aa<-aa[sfJoined$id,]
aa$locs<-sfJoined$name_long
aa$continent<-sfJoined$continent
aa<-aa[admin_bound_join$id,]
aa$admin_codes<-admin_bound_join$admin_code
locs<-unique(aa$locs)
admin_codes<-unique(aa$admin_codes)

aa[,dep_var]<-as.factor(aa[,dep_var])
aa$continent<-as.factor(aa$continent)
crisis<-which(aa[,dep_var]==1)
not_crisis<-which(aa[,dep_var]==0)
years<-unique(aa$year)
null_vals<-which(aa[,dep_var]==-99999)
y_loc<-paste0(round(aa$row/1),'_',round(aa$col/1))
y_loc_u<-unique(y_loc)
DT <- data.frame(longitude=aa$col,
                 latitude=aa$row+runif(nrow(aa))/10000)

pts = st_as_sf(DT, coords = c("longitude", "latitude"), 
               crs = 4326, agr = "constant")

month_n<-24
t_cols<-match(paste0('t',0:month_n),colnames(aa))
p_cols<-match(paste0('p',0:month_n),colnames(aa))
sum_t<-rowSums(aa[,t_cols])
sum_p<-rowSums(aa[,p_cols])
sd_t<-apply(aa[,t_cols],1,sd)
sd_p<-apply(aa[,p_cols],1,sd)
max_t<-apply(aa[,t_cols],1,max)
max_p<-apply(aa[,p_cols],1,max)
min_t<-apply(aa[,t_cols],1,min)
min_p<-apply(aa[,p_cols],1,min)
range_t<-max_t-min_t
range_p<-max_p-min_p
aaa<-cbind(aa,sum_t,sum_p,sd_t,sd_p,max_t,max_p,min_t,min_p,range_t,range_p)
add_var_ids<-match(c("sum_t","sum_p","sd_t","sd_p","max_t","max_p","min_t","min_p","range_t","range_p"),colnames(aaa))

sample_frac<-0.7
loc_frac<-0.1
smote<-2.5
best_models<-list()
n_spatial_bias<-0
n_tried_models<-0
while (length(best_models)<100){
  n_tried_models<-n_tried_models+1
  rand_locs<-sample(admin_codes,length(admin_codes)*loc_frac)
  train_ids<-which(!(aa$admin_codes%in%rand_locs))
  test_ids<-setdiff(1:nrow(aa),c(train_ids,null_vals))
  if (length(test_ids)>15000){
    test_ids<-sample(test_ids,length(test_ids))[1:15000]}
  crisis_<-intersect(train_ids,crisis)
  if (length(crisis_)>0){
    not_crisis_<-intersect(train_ids,not_crisis)
    samp_size<-length(crisis_)*sample_frac
    train_ids<-c(sample(crisis_,samp_size),sample(not_crisis_,samp_size*smote)) #SMOTE: Synthetic Minority Over-sampling Technique
    train_data<-aaa[train_ids,]
    test_data<-aaa[test_ids,]
    rf_train<-train_data[,c(dep_var,t_cols,p_cols,add_var_ids)]
    rf_test<-test_data[,c(dep_var,t_cols,p_cols,add_var_ids)]
    colnames(rf_train)[1]<-'fam_lev'
    colnames(rf_test)[1]<-'fam_lev'
    rf <-rfsrc(as.factor(fam_lev)~.,data=rf_train,ntree=500)
    pred_data<-predict(rf,rf_test)$predicted[,2]
    err<-c()
    for (tre in seq(0.3,0.8,0.1)){
      pred_test<-(1*(pred_data>tre))-as.numeric(as.character(rf_test$fam_lev))
      fn<-sum(pred_test==-1)/sum(rf_test$fam_lev==1)
      fp<-sum(pred_test==1)/sum(rf_test$fam_lev==0)
      acc<-sum(pred_test==0)/length(pred_test)
      tss<-(1-fp)+(1-fn)-1
      err<-rbind(err,c(tre,fp,fn,tss,acc))
    }
    best_tre<-err[which(err[,4]==max(err[,4])),]
    #test spatial autocorrelation
    pred_data<-1*(predict(rf,rf_train)$predicted[,2]>best_tre[1])
    res<-pred_data-as.numeric(as.character(rf_train$fam_lev))
    if (sum(abs(res))>0){
      res<-as.factor(res)
      levels(res)<-c("-1","0","1")
      nn5 = knn2nb(knearneigh(pts[train_ids,],5))
      w = nb2listw(nn5, style="B")
      p<-joincount.test(res, w)[[1]]$p.value
    } else {p<-1}
    if (p>0.05){
      best_models[[length(best_models)+1]]<-list(
        'train_data'=train_ids,
        'p'=p,
        'tre'=best_tre[1],
        'false_pos'=best_tre[2],
        'false_neg'=best_tre[3],
        'tss'=best_tre[4],
        'acc' = best_tre[5],
        'loc_f' = loc_frac,
        'smote' = smote,
        'sample_f' = sample_frac
      )
      tss<-as.numeric(lapply(best_models,function(x) x$tss))
      print (paste
             (length(best_models),paste('fp=',round(best_tre[2],2),
                                        'fn=',round(best_tre[3],2),
                                        'tss=',round(best_tre[4],2),
                                        'acc=',round(best_tre[5],2)
             )
             )
      )
    }else{
      
      n_spatial_bias<-n_spatial_bias+1
      print(paste('spatial autocorrelation',n_spatial_bias/n_tried_models))
      
    }   
  } 
}


train_data_n<-as.numeric(lapply(best_models,function(x) length(x$train_data)))
p<-as.numeric(lapply(best_models,function(x) x$p))
tss<-as.numeric(lapply(best_models,function(x) x$tss))
fp<-as.numeric(lapply(best_models,function(x) x$false_pos))
fn<-as.numeric(lapply(best_models,function(x) x$false_neg))
tre<-as.numeric(lapply(best_models,function(x) x$tre))
acc<-as.numeric(lapply(best_models,function(x) x$acc))


sink('/RESULTS/cross_validation_stats_admin_areas.txt')
print('accuracy')
summary(acc)
sd(acc)
print('tss')
summary(tss)
sd(tss)
print('false positives')
summary(fp)
sd(fp)
print('false positives')
summary(fn)
sd(fn)
print('threshold')
summary(tre)
sd(tre)
sink()




