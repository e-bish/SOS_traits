###########################################################################################################################################################
#### R-script producing data release for TOFF users. The data release includes three CSV files: one with trait measures, one with environmental measures and one with both measure types.  
#### It requires a login and password that can be obtained from the TOFF editorial board (see http://toff-project.univ-lorraine.fr).
#### R-script version 8 - January 2021
########################################################################################################################################################### 

# Libraries loading
if (!require(RPostgreSQL)) install.packages('RPostgreSQL')
library(RPostgreSQL) # to be installed in R (in R --> tab "packages" --> Install Package(s) --> select a mirror --> select the package RPostgreSQL))
if (!require(sf)) install.packages('sf')
library(sf) # to be installed in R (see before) 
if (!require(reshape2)) install.packages('reshape2')
library(reshape2) # to be installed in R (see before)) 
if (!require(writexl)) install.packages('writexl')
library(writexl) # to be installed in R (see before))
if (!require(easycsv)) install.packages('easycsv')
library(easycsv) # to be installed in R (see before))

###########################################################################################################################################################
#### 1. TOFF connection and setup
########################################################################################################################################################### 
 folder_path = easycsv::choose_dir() # cannot be modified. Here, you select the folder where you want to download your TOFF extraction
 hote="pggeodb.nancy.inra.fr" # cannot be modified
 port="5432" # cannot be modified
 database="db_toff" # cannot be modified
 schema_defaut="public" # cannot be modified
 user_name="embishop" # specify, here, your user name
 pswd="569752mdfesp4169fS"# specify, here, your password
 connection_db <- dbConnect(dbDriver("PostgreSQL"),dbname=database, user=user_name,password=pswd,host=hote)
 dbGetInfo(connection_db)
 setwd(folder_path)
 
###########################################################################################################################################################
#### 2. Functional trait dataset
###########################################################################################################################################################
 #SQL query for functional traits 
 requete_traits = "SELECT 
  public.t_measure_mea.mea_id,
  public.tr_family_trait_fam.fam_name,
  public.tr_trait_tra.tra_name,
  public.tr_genus_gen.gen_name,
  public.tr_species_spe.spe_name,
  public.t_measure_mea.mea_value,
  public.tr_sex_sex.sex_wording,
  public.tr_type_measure_tym.tym_wording,
  public.tr_stage_development_dev.dev_name,
  public.t_measure_mea.mea_date,
  public.t_measure_mea.mea_precision_date,
  public.t_measure_mea.mea_age,
  public.t_measure_mea.mea_nsize,
  public.t_measure_mea.mea_memo,
  public.t_measure_mea.mea_original_species,
  public.t_measure_mea.mea_degree_polyploidy,
  public.t_measure_mea.mea_encoding_date,
  public.t_measure_mea.mea_valide,
  public.t_publi_pub.pub_name,
  public.t_set_metadata_measure_met.met_id
  
FROM
  public.t_measure_mea
  INNER JOIN public.tr_trait_tra ON (public.t_measure_mea.mea_tra_id = public.tr_trait_tra.tra_id)
  INNER JOIN public.tr_species_spe ON (public.t_measure_mea.mea_spe_id = public.tr_species_spe.spe_id)
  INNER JOIN public.tr_family_trait_fam ON (public.tr_trait_tra.tra_fam_id = public.tr_family_trait_fam.fam_id)
  INNER JOIN public.tr_genus_gen ON (public.tr_species_spe.spe_gen_id = public.tr_genus_gen.gen_id)
  INNER JOIN public.tr_sex_sex ON (public.t_measure_mea.mea_sex_id = public.tr_sex_sex.sex_id)
  INNER JOIN public.tr_type_measure_tym ON (public.t_measure_mea.mea_tym_id = public.tr_type_measure_tym.tym_id)
  INNER JOIN public.tr_stage_development_dev ON (public.t_measure_mea.mea_dev_id = public.tr_stage_development_dev.dev_id)
  INNER JOIN public.t_set_metadata_measure_met ON (public.t_measure_mea.mea_met_id = public.t_set_metadata_measure_met.met_id)
  INNER JOIN public.t_publi_pub ON (public.t_set_metadata_measure_met.met_pub_id = public.t_publi_pub.pub_id)

ORDER BY
  fam_name,
  tra_name,
  gen_name,
  spe_name"
 
# Extracting and exporting functional trait dataframe
df_traits <- dbGetQuery (connection_db, requete_traits)
df_traits$FullName <- paste (df_traits$gen_name, df_traits$spe_name)
head(df_traits)
dim (df_traits)
summary(as.factor(df_traits$FullName), maxsum = 10000) # Displays a list of species for which functional trait datasets are already available in TOFF
unique (as.factor(df_traits$FullName))  # Displays the number of species for which functional trait datasets are already available in TOFF
unique (as.factor(df_traits$pub_name)) # Displays the number of publications used to obtain currently available functional trait datasets in TOFF
write.table(df_traits, file = paste(Sys.Date(), user_name, "df_traits.csv", sep = "."))# Export dataframe as a csv file 
write_xlsx(df_traits, path = paste(Sys.Date(), user_name, "df_traits.xlsx",sep = "."))# Export dataframe as an xlsx file


###########################################################################################################################################################
#### 3. Environmental features 
###########################################################################################################################################################
# All functional trait values recorded in TOFF are associated with information about measurements of their environmental context. The match between df_traits and df_metadata is possible through the mea_id column and/or met_id.

# SQL query for environmental features 
query_metadata = "SELECT 
 public.t_set_metadata_measure_met.met_id,
 public.t_measure_mea.mea_id,
 public.tr_characteristic_char.char_name,
 public.tj_measure_characteristic_mch.mch_value

FROM
  public.t_measure_mea
  INNER JOIN public.t_set_metadata_measure_met ON (public.t_measure_mea.mea_met_id = public.t_set_metadata_measure_met.met_id)
  INNER JOIN public.tj_measure_characteristic_mch ON (public.t_set_metadata_measure_met.met_id = public.tj_measure_characteristic_mch.mch_met_id)
  INNER JOIN public.tr_characteristic_char ON (public.tr_characteristic_char.char_id = public.tj_measure_characteristic_mch.mch_char_id)"

# Extracting and exporting environmental feature dataframe
df_metadata <- dbGetQuery (connection_db, query_metadata)
head(df_metadata)
dt_metadata=dcast(df_metadata,formula=mea_id + met_id ~ char_name, value.var='mch_value')
head(dt_metadata)
write.table(dt_metadata, file = paste(Sys.Date(), user_name, "dt_metadata.csv", sep = "."))# Export dataframe as a csv file  
write_xlsx(dt_metadata, path = paste(Sys.Date(), user_name, "df_metadata.xlsx",sep = "."))# Export dataframe as an xlsx file

###########################################################################################################################################################
#### 4. Merging Functional trait dataset Environmental features into one dataframe (Facultative)
############################################################################################################################################################ The two previous parts of the script return functional trait dataset and environmental dataset as two distinct files. Having trait values and their associated environmental measures in the same file can be useful. This fourth part of the script allows merging the two datasets into one file.

df_traitsandmetadata <- merge(dt_metadata,df_traits, all.x = TRUE, all.y = TRUE)# Merge dataframes into one
write.table(df_traitsandmetadata, file = paste(Sys.Date(), user_name, "df_merged_metadata_traits.csv", sep = "."))# Export dataframe 
write_xlsx(df_traitsandmetadata, path = paste(Sys.Date(), user_name, "df_merged_metadata_traits.xlsx",sep = "."))# Export dataframe in xlsx

###########################################################################################################################################################
#### 5. TOFF disconnection
########################################################################################################################################################### 
# Please, do not forget to logout when your dataset extraction is done.
dbDisconnect(connection_db)