# Load the libraries
library(tidyverse)
library(dplyr)
library(rEDM)
library(readr)
library(castr)


# Load the datasets 
zoo <- read_tsv("data/zoo_abundances.tsv.gz")
env <- read_tsv ("data/std_all_depths.tsv") 


##### Step 1 : Plankton

# Calculate biovolumes of each organism
plankton_biovol<-zoo%>%
  # convert the pixel variables in metric form
  mutate(
    area = area*(0.0106)^2,
    major = axis_major_length*0.0106,
    minor = axis_minor_length*0.0106
  ) |>
  # compute biovolume volume for each individual
  mutate(
    radius=sqrt(area/pi),
    spherical_vol=(4/3)*pi*(radius)^3,
    biovol=spherical_vol*concentration,
    biovol_ellipsoidal=((4/3) * pi * (major / 2) * ((minor / 2)^2))*concentration)


######Groups
# What groups do we have
levels(as.factor(zoo$label))

# Do some groupings
# On the feeding mode
gelatinous_filter_feeder <- c("Fritillariidae","Oikopleuridae","Salpida","Doliolida","tail<Appendicularia","nucleus<Salpidae", "trunk<Appendicularia" )
cruise_feeders_cops <- c("Centropagidae","Oncaeidae","Harpacticoida","Metridinidae")
ambush_feeders_cops <- c("Candaciidae","Oithonidae","Corycaeidae")
filter_feeder_cops <- c("Calanoida","Calanidae","Temoridae")
mix_feeder_cops <- c("Acartiidae")

gelatinous_carnivorous <- c("Aglaura", "Cnidaria", "Diphyidae", "Abylidae","Hydrozoa",
                         "Rhopalonema velatum", "Physonectae","part<Siphonophorae")
chaetognaths <- c("Chaetognatha","tail<Chaetognatha","head<Chaetognatha")

# On the source of feeding

carnivorous_cops_large <- c("Candaciidae")
carnivorous_cops_small <- c("Corycaeidae","Oncaeidae")
omnivorous_cops <- c("Centropagidae","Metridinidae")
herbivorous_cops <- c("Calanoida","Acartiidae","Oithonidae","Temoridae","Calanidae","Harpacticoida")



# Recode in the table
biovol_group <- plankton_biovol %>%
  mutate(
    feeding_mode = case_when(
      label %in% gelatinous_filter_feeder ~ "gelatinous_filter_feeder",
      label %in% cruise_feeders_cops ~ "cruise_feeders_cops",
      label %in% ambush_feeders_cops ~ "ambush_feeders_cops",
      label %in% mix_feeder_cops ~ "mix_feeder_cops",
      label %in% gelatinous_carnivorous ~ "gelatinous_carnivorous",
      label %in% filter_feeder_cops ~ "filter_feeder_cops",
      label %in% chaetognaths ~ "chaetognaths",
      TRUE ~ label 
    ),
    # feeding_source
    feeding_source = case_when(
      label %in% carnivorous_cops_large ~ "carnivorous_cops_large",
      label %in% carnivorous_cops_small ~ "carnivorous_cops_small",
      label %in% omnivorous_cops ~ "omnivorous_cops",
      label %in% herbivorous_cops ~ "herbivorous_cops",
      TRUE ~ feeding_mode  
    )
  )

# Did it work
levels(as.factor(biovol_group$feeding_mode))
levels(as.factor(biovol_group$feeding_source))
biovol_group%>%filter(feeding_mode=="cruise_feeders_cops")%>%select(label)%>%distinct()

##### TO DO : Travailler en biovolumes mais checker que concentration assez élevée !!
#Calculate the BIOVOLUME for each taxonomic group
taxo_FM <- biovol_group %>%
  group_by(feeding_mode, date.x)%>%
  dplyr::summarise(biovol=sum(biovol),
                   biovol_ellipsoidal=sum(biovol_ellipsoidal)) %>%
  ungroup()

taxo_FS <-  biovol_group %>%
  group_by(feeding_source, date.x)%>%
  dplyr::summarise(biovol=sum(biovol),
                   biovol_ellipsoidal=sum(biovol_ellipsoidal)) %>%
  ungroup()


#Replace plankton NA per 0
#crossing of all the possibilities
#uu<- bind_cols(crossing(date=taxo_data$object_date, taxo=taxo_data$label))
#uu$object_date<- uu$date
#uu$label<- uu$taxo

#merge with my data
#test<- merge(uu, taxo_data, by=c("object_date","label"), all.x=TRUE)
#dim(test)

#pheno_dataa<-test
#colnames(pheno_dataa)

#Remove duplicated columns
#pheno_dataa<- pheno_dataa%>% select(-taxo, -date)

#Replace NA by 0
#pheno_dataa<-pheno_dataa %>%
#  mutate(across(everything(), ~replace_na(.x, 0)))

#Works this way
#pheno_dataa <- pheno_dataa %>%
#  mutate(conc = replace_na(conc, 0))

#Represent the abundances of each taxonomic group with time
#taxo_data%>%ggplot(aes(x=object_date,y=conc))+geom_line(col="blue")+facet_wrap(~label, scales="free_y")+theme_bw()

#pheno_dataa%>%ggplot(aes(x=object_date,y=conc))+geom_line(col="blue")+facet_wrap(~label, scales="free_y")+theme_bw()

#For now, keep only a few taxa
categories <- c("gelatinous_filter_feeder","cruise_feeders_cops","ambush_feeders_cops","filter_feeder_cops","mix_feeder_cops","gelatinous_carnivorous","chaetognaths",
                "carnivorous_cops_large","carnivorous_cops_small","omnivorous_cops","herbivorous_cops") #maybe add gelatinous as well

p_abundances <- taxo_FS %>%
  filter(feeding_source %in% categories)

#Represent the abundances of each taxonomic group with time
p_abundances%>%ggplot(aes(x=date.x,y=biovol))+geom_line(col="blue")+facet_wrap(~feeding_source, scales="free_y")+theme_bw()+labs(xlab="date", ylab="biovolume")

#####################################################################################################################
#####Regularize the data #####
p_abundances <- p_abundances %>% mutate(object_date=date.x)
unique(p_abundances$object_date)

ref <- tibble(
  target_date=seq(from=as.Date("1967-01-05"), to=as.Date("2022-12-31"), by=14),
  year=year(target_date)
)
# identify years in which the number of obs is larger than usual
pbs <- ref |>
  count(year) |>
  filter(n>26)
# ->this is often an extra date in very late decembre => just remove it
ref <- filter(ref, !(year %in% pbs$year & yday(target_date) > 360))
ref %>%
  count(year) %>%
  filter(n>26)
# -> all OK

# match data based on these reference dates
avail_dates_p <- unique(p_abundances$object_date)
avail_dates_p<-as.data.frame(avail_dates_p)
avail_dates_p<- avail_dates_p%>%mutate(year=year(avail_dates_p))
ref_p<-ref%>%
  mutate(
    closest_date = castr::closest(ref$target_date, avail_dates_p$avail_dates_p),
    date_diff = abs(closest_date - target_date) |> as.numeric()
  )

ref_p%>%group_by(year)%>%summarise(count=n())%>%filter(count!=26)

plankton_match <- left_join(ref_p,p_abundances, by=c("closest_date"="object_date"), relationship="many-to-many")

# erase data for matches that are too far from the target
plankton_match <- plankton_match %>%
  mutate(conc = if_else(date_diff > 14, NA, biovol_ellipsoidal))

ggplot(plankton_match) +
  geom_point(aes(x=target_date, y=biovol_ellipsoidal), size=0.2) + theme_bw()+labs(x="date", y="valeur")+facet_wrap(~feeding_source)+scale_y_log10()

#Reput in the right format
plankton_ccm<- plankton_match%>% pivot_wider(names_from="feeding_source", values_from="biovol_ellipsoidal")


#####Linear interpolation if no value available
#Linear interpolation in order not to have NA

interpolate_column <- function(column) {
  approx(x = plankton_ccm$target_date, y = column, xout = plankton_ccm$target_date, method = "linear", rule=1)$y
}

plankton<- plankton_ccm%>% select(-target_date,-year,-closest_date,-date_diff)

# For all columns
interpolated_val <- lapply(plankton, interpolate_column)

# Convert to a dataframe
interpolated_valuess_p <- as.data.frame(interpolated_val)

# Rename them
colnames(interpolated_valuess_p) <- colnames(plankton)

plankton_for_ccm <- interpolated_valuess_p%>% mutate(date= plankton_ccm$target_date,
                                                     year=plankton_ccm$year)%>% select(-date.x, -biovol,-obejct_date,-conc)%>%distinct()


plankton_for_ccm<- plankton_for_ccm%>%mutate(across(
    .cols = -c(date, year),       
    .fns = ~log10(.x + 1),
    .names = "{.col}_log10"))



##############################################################################################################
##### Step 2: Environment ###################################################################################

head(env)

#Pivot it here
std_pivot <- env %>%pivot_wider(names_from = "name", values_from="value" )
std_pivot<- as.data.frame(std_pivot)

#Plot relationships between variables
#env_mean_pivot%>%dplyr::select(-date)%>%ggpairs()

#First step: Regularization

# Define a sequence of dates 
ref <- tibble(
  target_date=seq(from=as.Date("1967-01-05"), to=as.Date("2022-12-31"), by=14),
  year=year(target_date)
)

#Start in 1992 if you want to be consistent with CCM analyses but you can also remove this line if you want to have it on the whole time series
#ref<- ref%>% filter(target_date>"1992-01-16")

# identify years in which the number of obs is larger than usual
pbs <- ref %>%
  count(year) %>%
  filter(n>26)

# ->this is often an extra date in very late decembre => just remove it
ref <- filter(ref, !(year %in% pbs$year & yday(target_date) > 360))
ref %>%
  count(year) %>%
  filter(n>26)
# -> all OK

# Match data based on these reference dates
avail_dates <- unique(std_pivot$date)
ref <- ref %>%
  mutate(
    closest_date = castr::closest(ref$target_date, avail_dates),
    date_diff = abs(closest_date - target_date) %>% as.numeric()
  )

# Insert the data based on the matched dates
table_reg <- left_join(ref, std_pivot, by=c("closest_date"="date"), relationship="many-to-many")

# erase data for matches that are too far from the target
table_reg<- table_reg %>%
  mutate(value = if_else(date_diff > 6, NA, value))

#env_ccm <- table_reg%>%filter(depth=="10")%>%select("T","S","CHLA","NO3","O","target_date","year")%>%filter(target_date>"1992-01-01")

env_ccm <- table_reg %>%
  filter(target_date > as.Date("1992-01-01")) %>%
  group_by(target_date, year) %>%
  summarise(
    T = mean(T, na.rm = TRUE),
    S = mean(S, na.rm = TRUE),
    CHLA = mean(CHLA, na.rm = TRUE),
    logchla=log10(CHLA),
    NO3 = mean(NO3, na.rm = TRUE),
    O = mean(O, na.rm = TRUE),
    .groups = "drop"
  )


# Linear interpolation

env_ccm_filled <- env_ccm %>%
  pivot_longer(
    cols = c(T, S, CHLA, NO3, O),       
    names_to = "name",
    values_to = "value"
  ) %>%
  group_by(name)%>%
  arrange(target_date) %>%
  mutate(
    values_filled = interpolate(x = target_date, y = value, xout = target_date)
  ) %>% select(-value)%>%
  ungroup() %>%
  pivot_wider(
    names_from = name,
    values_from = values_filled
  )

env_ccm_filled <- env_ccm_filled%>% rename("date"="target_date")


####Step 3: Merge the two#####
#####Merge it with environment based on date #####
ccm_plankton_env<- merge(env_ccm_filled, plankton_for_ccm, by=c("date","year"))

#Plot it
ccm_plankton_env%>%ggplot()+geom_point(aes(x=date, y=T))

tail(ccm_plankton_env)
length(is.na(ccm_plankton_env))
length(ccm_plankton_env$date)

colSums(is.na(ccm_plankton_env))#no NA en fait

#Save it
ccm_plankton_env <- ccm_plankton_env %>% write_tsv("data/ccm_plankton_env")
