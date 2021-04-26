
##### integration of the bayclone_C code
rm(list=ls(all=TRUE))

### command to run this R script
### R CMD BATCH --no-save --no-restore '--args 9c949c63-aa63-4a84-ac5a-b5260775cee6 
### ./RUN/9c949c63-aa63-4a84-ac5a-b5260775cee6/inputData/9c949c63-aa63-4a84-ac5a-b5260775cee6.dat 
### ./bayclone_C_code ./RUN/9c949c63-aa63-4a84-ac5a-b5260775cee6/rda_out ./RUN/9c949c63-aa63-4a84-ac5a-b5260775cee6/fig_out ./RUN/9c949c63-aa63-4a84-ac5a-b5260775cee6/txt_out' 
### ./bayclone_C_code/bayclone_C_main.R ./RUN/9c949c63-aa63-4a84-ac5a-b5260775cee6/logs/9c949c63-aa63-4a84-ac5a-b5260775cee6.log

find_mult_for_each_SNV <-function(data,pur,consensus_ploidy,dir_path,config_param){
  
  CHR_NAME_INDX = config_param$CHR_NAME_INDX
  CHR_POS_INDX = config_param$CHR_POS_INDX
  
  TOTAL_RD_INDX = config_param$TOTAL_RD_INDX       
  VARIANT_RD_INDX = config_param$VARIANT_RD_INDX
  SAMPLE_T_CN_INDX = config_param$SAMPLE_T_CN_INDX
  
  MAX_CN_CUTOFF = config_param$MAX_CN_CUTOFF
  PLOIDY_MAX_CN_CUTOFF = config_param$PLOIDY_MAX_CN_CUTOFF
  
  
  f1 = sprintf("./%s/%s_mult_ccf.dat",dir_path,sample_id)
  
  VAF = data[,VARIANT_RD_INDX]/data[,TOTAL_RD_INDX]
  
  vaf_adj = (pur*data[,SAMPLE_T_CN_INDX] + (1-pur)*2)/pur
  
  #### mutation CN
  mutation_CN_for_SNV = VAF*vaf_adj
  
  nSNV = length(data[,VARIANT_RD_INDX])          
  mult = rep(0,nSNV)
  
  ##### this part is not needed for PCAWG @Jan21 
  ##### calculates ploidy from sample_avg_tumor_CN #####
  id = which(data[,SAMPLE_T_CN_INDX] < PLOIDY_MAX_CN_CUTOFF)
  S = length(id)
  ploidy = sum(data[id,SAMPLE_T_CN_INDX] )/S
  
  ##############
  
  mult = mutation_CN_for_SNV
  
  data1 = cbind(data,mutation_CN_for_SNV,mult)
  names(data1) = c("chr","pos","N","n","M_B","nMaj1_A","nMin1_A","frac1_A","nMaj2_A","nMin2_A","frac2_A","mutation_CN","multiplicity")
  write.table(data1,row.names=F,sep="\t",quote=F,file=f1)
  
  ###############
  
  f2 = sprintf("./%s/%s_multiplicity.txt",dir_path,sample_id)
  chr_name = as.character(data[,CHR_NAME_INDX])
  chr_loc = data[,CHR_POS_INDX]
  tumor_CN = data[,SAMPLE_T_CN_INDX]
  
  col4=rep(NA,length(mult))
  col5=rep(NA,length(mult))
  T1 = cbind(chr_name,chr_loc,tumor_CN,mult,col4,col5)
  write.table(T1,col.names=c("chr","pos","tumour_copynumber","multiplicity","multiplicity_options","probabilities"),row.names=F,sep="\t",quote=F,file=f2)
  
  
  f3 = sprintf("./%s/%s_purity_ploidy.txt",dir_path,sample_id)
  T1 = cbind(sample_id,pur,ploidy)#consensus_ploidy)  ### ploidy is our estimate
  write.table(T1,col.names=F,row.names=F,sep="\t",quote=F,file=f3)
  
  res = NULL
  res$mult = mult
  res$ploidy = ploidy
  
  return(res)
}

find_nclust_mult_of_SNV_no_merge <- function(mult,config_param){
  
  MAX_CLUSTER = config_param$MAX_CLUSTER  ### 7 for pcawg
  MIN_CLUST_DIST = config_param$MIN_CLUST_DIST
  nClust = -1
  MULT_Info = NULL
  
  MULT_Info$mult = mult
  
  MC = Mclust(mult,MAX_CLUSTER)
  
  MULT_Info$model = MC
  
  #### compute the prob(mult_s > 1) and prob(mult_s < 1)
  prob_gt_1_vec = MC$z%*%(1-pnorm(1,MC$parameters$mean,sqrt(MC$parameters$variance$sigmasq)))
  prob_lt_1_vec = MC$z%*%(pnorm(1,MC$parameters$mean,sqrt(MC$parameters$variance$sigmasq)))
  
  #cat("all means ",MC$parameters$mean,"\n")
  
  uniq_class = sort(unique(MC$classification))
  uniq_class_len =  length(uniq_class)
  #cat("uniq_class ", uniq_class, "\n")
  
  MULT_Info$old_label = 1:length(MC$parameters$mean)
  MULT_Info$mean_mult = MC$parameters$mean
  
  MULT_Info$mut_assign_vec = MC$classification
  MULT_Info$old_size = as.vector(table(MC$classification))
  
  
  MULT_Info$prob_gt_1_vec = prob_gt_1_vec
  MULT_Info$prob_lt_1_vec = prob_lt_1_vec
  
  return(MULT_Info)
}

# Edward(Dehua)_190313: I added the find_thre_vec functioni to calculate the probability greather than
#                       a threshold or less than a threshold value. It is based on the saved GMM model,
#                       the MC variable I added in the find_nclust_mult_of_SNV_no_merge function. It
#                       takes one extra input variable, the threshold_value. We can put q^tile[k] there,
#                       or we can also put p[1]q^tile1[k] + p[2]q^tile2[k] there.

find_thre_vec <- function(MULT_Info, threshold_value){
  
  Threshold_val_Info = NULL
  MC = MULT_Info$model
  prob_gt_thre_vec = MC$z%*%(1-pnorm(threshold_value,MC$parameters$mean,sqrt(MC$parameters$variance$sigmasq)))
  prob_lt_thre_vec = MC$z%*%(pnorm(threshold_value,MC$parameters$mean,sqrt(MC$parameters$variance$sigmasq)))
  Threshold_val_Info$prob_gt_thre_vec = prob_gt_thre_vec
  Threshold_val_Info$prob_lt_thre_vec = prob_lt_thre_vec
  
  return(Threshold_val_Info)
}

# Edward(Dehua)_190327: function to get mode

getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

generate_result_based_on_mclust <- function(data,MULT_Info,pur,dir_path,input_dir,config_param,box1_thres,box2_thres){
  
  lda_flag1 = FALSE
  lda_flag2 = FALSE
  lda_flag3 = FALSE
  
  CHR_NAME_INDX = config_param$CHR_NAME_INDX
  CHR_POS_INDX = config_param$CHR_POS_INDX
  
  TOTAL_RD_INDX = config_param$TOTAL_RD_INDX
  VARIANT_RD_INDX = config_param$VARIANT_RD_INDX
  
  SAMPLE_T_CN_INDX = config_param$SAMPLE_T_CN_INDX
  
  MAJOR1_CN_INDX = config_param$MAJOR1_CN_INDX
  MINOR1_CN_INDX = config_param$MINOR1_CN_INDX
  FRAC1_INDX = config_param$FRAC1_INDX
  
  MAJOR2_CN_INDX = config_param$MAJOR2_CN_INDX
  MINOR2_CN_INDX = config_param$MINOR2_CN_INDX
  FRAC2_INDX = config_param$FRAC2_INDX
  MEAN_CENTER_CUTOFF = config_param$MEAN_CENTER_CUTOFF
  
  new_number_of_cluster = length(MULT_Info$mean_mult)-1
  
  mean_mult = MULT_Info$mean_mult
  old_label = MULT_Info$old_label
  
  mut_assign_vec = MULT_Info$mut_assign_vec
  old_size = MULT_Info$old_size
  
  prob_gt_1_vec = MULT_Info$prob_gt_1_vec
  prob_lt_1_vec = MULT_Info$prob_lt_1_vec
  
  
  #cat("mean_mult from mclust without merge = ",mean_mult,"\n")
  
  f3_1 = sprintf("./%s/%s_clustered_mult.txt",dir_path,sample_id)
  chr_name = as.character(data[,CHR_NAME_INDX])
  chr_loc = data[,CHR_POS_INDX]
  T3_1 = data.frame(chr_name,chr_loc,mut_assign_vec,mean_mult[mut_assign_vec])
  write.table(T3_1,col.names=c("chr","pos","cluster","est_mean_mult"),row.names=F,quote=F,sep="\t",file=f3_1)
  
  ##########################################################################
  ##### new logic to assign clonal(1) or subclonal(0) or unknown(-1) status
  ##########################################################################
  ### this is the probability that multiplicity is less than 1
  
  # Edward190724: Gridsearch idea: threshold value depends on the result of the grid search
  #               We will have two thresholds, one for m_s < 1, one for other two conditions
  if(config_param$REAL_DATA_FLAG == 0){
    PROB_THRESHOLD_1 = box1_thres
  }
  # Edward190724: Not sure about real data, have NOT change the logic for real data yet!
  #               Needs further investigation.
  #               TODO: Find a way for real data.
  if(config_param$REAL_DATA_FLAG == 1){
    # Edward: This parts generate the s* list.
    prob_less_than_1 = MULT_Info$prob_lt_1_vec
    prob_less_than_1_sorted_des = prob_less_than_1[order(prob_less_than_1,decreasing=T),]
    L = length(prob_less_than_1_sorted_des)
    cum_vec = cumsum(1-prob_less_than_1_sorted_des)/(1:L)
    id_less_threshold = which(cum_vec < config_param$FDR_THRESHOLD)
    L_thr = length(id_less_threshold)
    if(L_thr > 0){
      PROB_THRESHOLD_1 = prob_less_than_1_sorted_des[L_thr+1]
      cat("L_thr+1 = ",L_thr+1, "val = ", prob_less_than_1_sorted_des[L_thr+1],"\n")
      cat("L_thr  = ",L_thr, "val = ", prob_less_than_1_sorted_des[L_thr],"\n")
      cat("cum_vec[L_thr+1] ",cum_vec[L_thr+1], "cum_vec[L_thr] = ", cum_vec[L_thr],"\n")
    }else{
      PROB_THRESHOLD_1 = 1.0 ### no SNV is subclonal	
    }
    
  }
  
  #cat("PROB_THRESHOLD_1 = ",PROB_THRESHOLD_1,"\n") 
  ### based on 10 samples this is the best threshold
  #PROB_THRESHOLD_1 = 0.7942314
  #PROB_THRESHOLD_1 = 0.735
  
  # Edward: This part runs through the logic of section 3.3
  if(!is.na(PROB_THRESHOLD_1)){
    # Edward190724: Gridsearch idea: threshold value depends on the result of the grid search
    #               We will have two thresholds, one for m_s < 1, one for other two conditions
    if(config_param$REAL_DATA_FLAG == 0){
      PROB_THRESHOLD_qs1 = box2_thres
      PROB_THRESHOLD_qs2 = box2_thres
    }
    # Edward190724: Not sure about real data, have NOT change the logic for real data yet!
    #               Needs further investigation.
    #               TODO: Find a way for real data.
    if(config_param$REAL_DATA_FLAG == 1){
      PROB_THRESHOLD_qs1 = box2_thres
      PROB_THRESHOLD_qs2 = box2_thres
    }
    
    mult_vec = MULT_Info$mult
    snv_status = rep(-1,length(mut_assign_vec))
    qs_list = rep(0, length(mut_assign_vec))   
    prob_qs_list = rep(0, length(mut_assign_vec))  
    
    for (i in 1:length(mut_assign_vec)){
      
      if(MULT_Info$prob_lt_1_vec[i] > PROB_THRESHOLD_1){            
        # Edward: save qs for LDA
        qs_list[i] = max(data[i,MAJOR1_CN_INDX],data[i,MINOR1_CN_INDX])
        lda_flag1 = TRUE
        snv_status[i] = 0  ### subclonal
      }else{    
        ### we need to see CN file
        if(data[i,FRAC1_INDX] == 1 ){ ## CN segment is clonal         
          q_2s = max(data[i,MAJOR1_CN_INDX],data[i,MINOR1_CN_INDX])   
          qs_list[i] = q_2s
          
          Threshold_val_Info = find_thre_vec(MULT_Info, q_2s)         
          prob_qs_list[i] = Threshold_val_Info$prob_gt_thre_vec[i]
          if(Threshold_val_Info$prob_gt_thre_vec[i] > PROB_THRESHOLD_qs1){
            lda_flag2 = TRUE
            snv_status[i] = 1 ### clonal
          }else{
            lda_flag3 = TRUE
            snv_status[i] = -1                                        
          }
        }else{ ### CN segment is subclonal                            
          f1_times_q1 = data[i,FRAC1_INDX]*max(data[i,MAJOR1_CN_INDX],data[i,MINOR1_CN_INDX])
          f2_times_q2 = data[i,FRAC2_INDX]*max(data[i,MAJOR2_CN_INDX],data[i,MINOR2_CN_INDX])
          f_times_q = f1_times_q1 + f2_times_q2
          qs_list[i] = f_times_q
          Threshold_val_Info = find_thre_vec(MULT_Info, f_times_q)
          prob_qs_list[i] = Threshold_val_Info$prob_gt_thre_vec[i]
          if(Threshold_val_Info$prob_gt_thre_vec[i] > PROB_THRESHOLD_qs2){
            lda_flag2 = TRUE
            snv_status[i] = 1   ### clonal
          }else{
            lda_flag3 = TRUE
            snv_status[i] = -1                                       
          }
          
        }       
      }
    }
  }
  
  if(config_param$REAL_DATA_FLAG == 0 & !is.na(PROB_THRESHOLD_1)){
    
    f3_3 = sprintf("./%s/%s.est_mutclass.dat",dir_path,sample_id)
    chr_name = as.character(data[,CHR_NAME_INDX])
    chr_loc = data[,CHR_POS_INDX]
    T3_3 = cbind(chr_name,chr_loc,snv_status)
    write.table(T3_3,col.names=c("chr","pos","status"),row.names=F,quote=F,sep="\t",file=f3_3)
    
    # Edward: construct the data for LDA
    f4_3 = sprintf("./%s/%s.all_data_est_mutclass_realdata0.dat",dir_path,sample_id)
    char_N = data[,TOTAL_RD_INDX]
    char_n = data[,VARIANT_RD_INDX]
    char_M = data[,SAMPLE_T_CN_INDX]
    char_maj1 = data[,MAJOR1_CN_INDX]
    char_min1 = data[,MINOR1_CN_INDX]
    char_farc1 = data[,FRAC1_INDX]
    char_maj2 = data[,MAJOR2_CN_INDX]
    char_min2 = data[,MINOR2_CN_INDX]
    char_frac2 = data[,FRAC2_INDX]
    T4_3 = cbind(chr_name,chr_loc,char_N,char_n,char_M,mean_mult[mut_assign_vec],qs_list,prob_qs_list,char_maj1,char_min1,char_farc1,char_maj2,char_min2,char_frac2,snv_status)
    write.table(T4_3,col.names=c("chr","pos","N", "n", "M", "mult", "Qs", "Prob_Qs", "maj1", "min1", "frac1", "maj2", "min2", "frac2", "status"),row.names=F,quote=F,sep="\t",file=f4_3)
  }
  
  if(!is.na(PROB_THRESHOLD_1)){	
    f3_4 = sprintf("./%s/%s.est_mutclass_with_prob.dat",dir_path,sample_id)
    chr_name = as.character(data[,CHR_NAME_INDX])
    chr_loc = data[,CHR_POS_INDX]
    T3_4 = cbind(chr_name,chr_loc,MULT_Info$prob_lt_1_vec,snv_status)
    write.table(T3_4,col.names=c("chr","pos","Pr(mult<1)","status"),row.names=F,quote=F,sep="\t",file=f3_4)
    
    # Edward: construct the data for LDA
    f4_4 = sprintf("./%s/%s.all_data_est_mutclass_realdata1.dat",dir_path,sample_id)
    char_N = data[,TOTAL_RD_INDX]
    char_n = data[,VARIANT_RD_INDX]
    char_M = data[,SAMPLE_T_CN_INDX]
    char_maj1 = data[,MAJOR1_CN_INDX]
    char_min1 = data[,MINOR1_CN_INDX]
    char_farc1 = data[,FRAC1_INDX]
    char_maj2 = data[,MAJOR2_CN_INDX]
    char_min2 = data[,MINOR2_CN_INDX]
    char_frac2 = data[,FRAC2_INDX]
    T4_4 = cbind(chr_name,chr_loc,char_N,char_n,char_M,mean_mult[mut_assign_vec],qs_list,prob_qs_list,char_maj1,char_min1,char_farc1,char_maj2,char_min2,char_frac2,snv_status)
    write.table(T4_4,col.names=c("chr","pos","N", "n", "M", "mult", "Qs", "Prob_Qs", "maj1", "min1", "frac1", "maj2", "min2", "frac2", "status"),row.names=F,quote=F,sep="\t",file=f4_4)
    if(lda_flag1 == FALSE || lda_flag2 == FALSE || lda_flag3 == FALSE){
      f4_4 = ""
    }
    #cat("f4_4 = ", f4_4,"\n")
    return(f4_4) # Edward: save the final output table for LDA
  }
}

# Edward(Dehua)_190319: The lda_classify_unknown function implements the LDA algorithm. It is 
#                       currently training on N, n and M values. Looking to include other features.

lda_classify_unknown<-function(data_file, dir_path, sample_id){
  data = read.table(data_file,header=T)
  training_data_idx <- which(data$status != '-1')
  training_data <- as.data.frame(data[training_data_idx,c(2:7,15)])
  predict_data_idx <- which(data$status == '-1')
  predict_data <- as.data.frame(data[predict_data_idx,c(2:7)])
  result_lda <- lda(formula = status ~ ., data = training_data)
  # predictions
  predict_lda <- predict(result_lda, newdata=predict_data)$class
  final_data <- cbind(data,data$status)
  final_data[predict_data_idx,16] <- as.numeric(levels(predict_lda))[predict_lda]
  f5 = sprintf("./%s/%s.all_data_lda_est_mutclass.dat",dir_path,sample_id)
  write.table(final_data, col.names=c("chr","pos","N", "n", "M", "mult", "Qs","Prob_Qs", "maj1", "min1", "frac1", "maj2", "min2", "frac2", "status", 'lda_status'), row.names=F,quote=F,sep="\t",file=f5)
  data_result_2 <- cbind(data[,1], data[,2], final_data[,16])
  f6 = sprintf("./%s/%s.lda_est_mutclass.dat",dir_path,sample_id)
  write.table(data_result_2, col.names=c("chr","pos", 'lda_status'), row.names=F,quote=F,sep="\t",file=f6)
  # check CV accuracy, when fit, if we set CV=True, we cannot use the predict function.
  lda_fit_crossvalid <- lda(formula = status ~ ., data = training_data, CV=T)
  ct <- table(training_data$status, lda_fit_crossvalid$class)
  lda_acc = sum(diag(prop.table(ct)))
  #cat("LDA accuracy = ",lda_acc,"\n")
}


###############################

require(mclust)
require(fpc)
require(pROC)
require(MASS)    # Edward_190318: Required library for LDA
args <- commandArgs(TRUE)
set.seed(17345)

###############################

config_param = NULL

config_param$CHR_NAME_INDX = 1
config_param$CHR_POS_INDX = 2

config_param$TOTAL_RD_INDX = 3
config_param$VARIANT_RD_INDX = 4

config_param$SAMPLE_T_CN_INDX = 5

config_param$MAJOR1_CN_INDX = 6
config_param$MINOR1_CN_INDX = 7
config_param$FRAC1_INDX = 8
config_param$MAJOR2_CN_INDX = 9
config_param$MINOR2_CN_INDX = 10
config_param$FRAC2_INDX = 11

config_param$MAX_CLUSTER = 6 ### 7 for PCAWG
config_param$MIN_CLUST_DIST = 0.1

config_param$MAX_CN_CUTOFF = 6.0  ### for multiplicty calculation
config_param$PLOIDY_MAX_CN_CUTOFF = 15.0  ### for ploidy calculation previously 20
config_param$MEAN_CENTER_CUTOFF = 0.9 #0.8 0.95

config_param$SAMPLE_NAME_INDX = 1
config_param$PURITY_INDX = 2
config_param$PLOIDY_INDX = 3
config_param$REAL_DATA_FLAG = 0
config_param$FDR_THRESHOLD = 0.06

##############################
sample_id = args[1]
input.file = args[2]
R_code_src_dir = args[3]
rda_output_dir = args[4]
out_fig_folder = args[4]#args[5]
out_txt_folder = args[4]#args[6]
##############################

if(length(grep("A15E",sample_id)) > 0 ){
  config_param$REAL_DATA_FLAG = 1
  cat("This is a real sample !!\n")
}

SAMPLE_NAME_INDX = config_param$SAMPLE_NAME_INDX
PURITY_INDX = config_param$PURITY_INDX
PLOIDY_INDX = config_param$PLOIDY_INDX


##############################
data = data.frame(read.table(input.file))

#### Jan 28 update to read from individual purity file
dir_input_file = dirname(input.file)
pur_ploidy_file_name = sprintf("%s/%s.pur.ploidy.dat",dir_input_file,sample_id)
pur_ploidy_val_str = read.table(pur_ploidy_file_name,header=T)
pur_val_from_file = as.numeric(unlist(pur_ploidy_val_str)[PURITY_INDX])
ploidy_val_from_file = as.numeric(unlist(pur_ploidy_val_str)[PLOIDY_INDX])
##########

cat("purity and ploidy = ",pur_val_from_file,ploidy_val_from_file,"\n")

Sys.time()
time1 = Sys.time()

res_mult = find_mult_for_each_SNV(data,pur_val_from_file,ploidy_val_from_file,out_txt_folder,config_param) 

MULT_Info = NULL
MULT_Info = find_nclust_mult_of_SNV_no_merge(res_mult$mult,config_param)

input_dir = dirname(input.file)

best_box1 = 0.86 #0.845 #0.845 #0.84
best_box2 = 0.16 #0.33 #0.53 #0.33

output_file_dir = generate_result_based_on_mclust(data,MULT_Info,pur_val_from_file,out_txt_folder,input_dir,config_param,best_box1,best_box2)
if(output_file_dir != ""){
  cat("outfile = ",output_file_dir,"\n")
  lda_classify_unknown(output_file_dir, out_txt_folder, sample_id)
}


######################################
orig_mutclass_file = sprintf("%s/%s.mutclass.dat",input_dir,sample_id)
est_mutclass_file = sprintf("%s/%s.est_mutclass.dat",out_txt_folder,sample_id)
LDA_est_mutclass_file = sprintf("%s/%s.lda_est_mutclass.dat",out_txt_folder,sample_id)
  
if(config_param$REAL_DATA_FLAG == 0){
    
    tryCatch({
      A_0 = read.table(orig_mutclass_file,header=T)
      # compare with original
      A_1 = read.table(est_mutclass_file,header=T)
      
      #cat(A_1[1:20,3],"\n")
      
      id_not_unknown = which(A_1[,3] != -1)
      match_id = which(A_0[,3] == A_1[,3])  #### which are identical class
      true_pos_id = which((A_0[,3] == A_1[,3])&(A_1[,3] == 1))
      true_neg_id = which((A_0[,3] == A_1[,3])&(A_1[,3] == 0))
      false_pos_id = which((A_0[,3] != A_1[,3])&(A_1[,3] == 1))
      false_neg_id = which((A_0[,3] != A_1[,3])&(A_1[,3] == 0))
      
      not_unknown_cnt = length(id_not_unknown)
      mat_cnt = length(match_id)
      total_cnt = length(A_0[,1])
      tp_cnt = length(true_pos_id)
      tn_cnt = length(true_neg_id)
      fp_cnt = length(false_pos_id)
      fn_cnt = length(false_neg_id)
      
      prop_in_not_unknown = mat_cnt/not_unknown_cnt
      prop_in_total = mat_cnt/total_cnt
      sensitivity = tp_cnt/(tp_cnt + fn_cnt)
      specificity = tn_cnt/(tn_cnt+fp_cnt)
      precision = tp_cnt/(tp_cnt+fp_cnt)
      
      f_out = sprintf("./%s/%s.logic_result.dat",out_txt_folder,sample_id)
      T3_2 = cbind(sample_id,prop_in_not_unknown,prop_in_total,sensitivity,specificity,precision)
      write.table(T3_2,col.names=c("sample_id","prop_in_not_unknown","prop_in_total","sens","spec","prec"),row.names=F,quote=F,sep="\t",file=f_out)
      
      result_summary = NULL
      result_summary$tot_SNV_orig = length(A_0[,3])
      result_summary$orig_n_clonal_SNV = length(which(A_0[,3]==1))
      result_summary$orig_n_subclonal_SNV = length(which(A_0[,3]==0))
      
      result_summary$tot_SNV_called = length(A_1[,3])
      result_summary$called_n_clonal_SNV = length(which(A_1[,3]==1))
      result_summary$called_n_subclonal_SNV = length(which(A_1[,3]==0))
      result_summary$called_n_unknown_SNV = length(which(A_1[,3]==-1))
      
      write.table(result_summary,col.names=c("total_orig_SNVs","orig_clonal_cnt","orig_subclonal_cnt","total_called_SNVs","called_clonal_cnt","called_subclonal_cnt","called_n_unknown_SNV"),row.names=F,quote=F,sep="\t",file=f_out,append=TRUE)
      
      # compare with lda
      A_2 = read.table(LDA_est_mutclass_file,header=T)
      
      #cat(A_2[1:20,3],"\n")
      
      id_not_unknown = which(A_2[,3] != -1)
      match_id = which(A_0[,3] == A_2[,3])  #### which are identical class
      
      not_unknown_cnt = length(id_not_unknown)
      mat_cnt = length(match_id)
      total_cnt = length(A_0[,1])
      true_pos_id = which((A_0[,3] == A_1[,3])&(A_1[,3] == 1))
      true_neg_id = which((A_0[,3] == A_1[,3])&(A_1[,3] == 0))
      false_pos_id = which((A_0[,3] != A_1[,3])&(A_1[,3] == 1))
      false_neg_id = which((A_0[,3] != A_1[,3])&(A_1[,3] == 0))
      
      prop_in_not_unknown = mat_cnt/not_unknown_cnt
      prop_in_total = mat_cnt/total_cnt
      sensitivity = tp_cnt/(tp_cnt + fn_cnt)
      specificity = tn_cnt/(tn_cnt+fp_cnt)
      precision = tp_cnt/(tp_cnt+fp_cnt)
      
      f_out_1 = sprintf("./%s/%s.lda_result.dat",out_txt_folder,sample_id)
      T3_3 = cbind(sample_id,prop_in_not_unknown,prop_in_total,sensitivity,specificity,precision)
      write.table(T3_3,col.names=c("sample_id","prop_in_not_unknown","prop_in_total","sens","spec","prec"),row.names=F,quote=F,sep="\t",file=f_out_1)
      
      cat("accuracy after LDA: ",prop_in_total, "\n")
      
      result_summary = NULL
      result_summary$tot_SNV_orig = length(A_0[,3])
      result_summary$orig_n_clonal_SNV = length(which(A_0[,3]==1))
      result_summary$orig_n_subclonal_SNV = length(which(A_0[,3]==0))
      
      result_summary$tot_SNV_called = length(A_2[,3])
      result_summary$called_n_clonal_SNV = length(which(A_2[,3]==1))
      result_summary$called_n_subclonal_SNV = length(which(A_2[,3]==0))
      
      
      write.table(result_summary,col.names=c("total_orig_SNVs","orig_clonal_cnt","orig_subclonal_cnt","total_called_SNVs","called_clonal_cnt","called_subclonal_cnt"),row.names=F,quote=F,sep="\t",file=f_out,append=TRUE)
    }, error = function(e_r){
      cat("NO FILE Exists !!\n")
    })	
}

##############################
cat("\n analysis finished and generated files for sample ",sample_id,"!!!\n")
##############################


####################################


