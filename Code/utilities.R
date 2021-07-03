
#################################
##Performance Evaluation Functions
#################################

is_positive <- function(estimated_cp, possible_positives, threshold)
{
  #Input:
  #estimated_cp -  1 estimate of a changepoint 
  #possible positives - actual changepoints with those already counted as true positives removed
  #threshold - acceptable tolerance for estimates to still be counted as a true positive
  
  #Output: 
  #Returns 0 if false positive or the index of the closest true positive it corresponds to    
  
  output = 0
  
  diff_vec = possible_positives - estimated_cp
  diff_vec = abs(diff_vec)
  
  logical = diff_vec <= threshold
  logical = any(logical)
  
  if(is.na(logical)){return(output)}
  
  else if(logical)
  {
    min_diff = min(diff_vec)
    index = match(min_diff, diff_vec)
    output = index
  }
  
  return(output)
  
}

update_positives <- function(possible_positives, index)
{
  #simple function to stop double counting of true positives
  m = length(possible_positives)
  l = index - 1
  r = index + 1
  
  left = possible_positives[1:l]
  right = possible_positives[r:m]
  if(m == 1 & index == 1){remaining = c(NA)}
  
  else if(index == m){remaining = left}
  
  else if(index == 1){remaining = right}
  
  else{remaining = c(left, right)}
  
  return(remaining)
}


missed_positives <- function(estimated_cps, actual_cps, threshold)
{
  #function which returns the number of changepoints not correctly estimated by the procedure
  #used to count total true / false positives
  
  remaining_positives = actual_cps
  
  for(estimate in estimated_cps)
  {
    int = is_positive(estimate, remaining_positives, threshold)
    
    if(int > 0)
    {
      remaining_positives = update_positives(remaining_positives, int)
    }
    
  }
  
  m = length(remaining_positives)
  
  if(all(is.na(remaining_positives))){return(0)}
  else{return(m)}
  
}


confusion_matrix <- function(data_length, estimated_cps, actual_cps, threshold)
{
  #this func returns the confusion matrix in vector form
  false_negatives = missed_positives(estimated_cps, actual_cps, threshold)
  
  true_positives = length(actual_cps) - false_negatives
  
  false_positives = length(estimated_cps) - true_positives
  
  true_negatives = data_length - false_negatives - true_positives - false_positives
  
  confusion_vec = c(true_positives, false_positives, true_negatives, false_negatives)
  
  return(confusion_vec)
  
}

precision <- function(confusion_matrix)
{
 
  if(confusion_matrix[1] == 0 & confusion_matrix[2] == 0){pr = 1}
  else{pr = confusion_matrix[1] / (confusion_matrix[1] + confusion_matrix[2])}
  
  return(pr)
  
}


recall <- function (confusion_matrix)
{

  if(confusion_matrix[1] == 0 & confusion_matrix[4] == 0){r = 0}
  else(r = confusion_matrix[1] / (confusion_matrix[1] + confusion_matrix[4]))
  
  return(r)
  
}


f1_score <- function(confusion_matrix)
{
  p = precision(confusion_matrix)
  r = recall(confusion_matrix)
 
  f1 = 2 * (p * r) / (p + r)
  return(f1)
}

# mean_square_error <- function(estimated_cps, actual_cps)
# {
#   M=length(estimated_cps)
#   sum=0
#   
#   for(est_cp in estimated_cps)
#   {
#     diff = (actual_cps - est_cp)^2
#     min_diff = min(diff)
#     sum = sum + min_diff
#     
#   }
#   
#   return(sum/M)
#   
# }

evaluation_matrix <- function(experiment_changepoints, evaluation_method, evaluation_params)
{
  #Saving the matrix of true changepoints
  changepoints = evaluation_params$actual_cps
  
  #setting dimension 
  cp_dim = dim(experiment_changepoints)
  num_simulations = cp_dim[1]
  num_penalties = cp_dim[3]
  
  for (i in 1 : num_penalties)
  {
    for(j in 1 : num_simulations)
    {
      
      #converting from sparse matrix back to a vector filled with the locations
      
      estimated_instance = experiment_changepoints[j,,i]
      estimated_instance = which(estimated_instance>0)
      
      actual_instance = changepoints[j,]
      actual_instance = which(actual_instance>0)
      
      evaluation_params$estimated_cps = estimated_instance
      evaluation_params$actual_cps = actual_instance
      
      #initialising the matrix on the first loop as we do not know the dimension of the evaluation matrix
      if(i == 1 & j == 1)
      {
        first_eval = do.call(evaluation_method, evaluation_params)
        eval_dim = length(first_eval)
        init_dim = c(num_simulations, eval_dim, num_penalties)
        eval_matrix = array(0, dim = init_dim)
        
        eval_matrix[j,,i] = first_eval
      }
      
      else
      {
        eval_matrix[j,,i] = do.call(evaluation_method, evaluation_params)
      }
      
    }
  }
  
  return(eval_matrix)     
  
}

evaluation_averages <- function(evaluation_matrix, eval_transform)
{
  #input should be the result of calling the evaluation_matrix with the confusion matrix method
  #the result is a matrix consisting of columns of the monte carlo average precision / recall for 
  #each set of experimental parameters 
  
  mat_dim = dim(evaluation_matrix)
  
  num_simulations = mat_dim[1]
  num_penalties = mat_dim[3]
  num_cps = mat_dim[4]
  num_means = mat_dim[5]
  
  dim_transform = length(eval_transform(evaluation_matrix[1,,1,1,1]))
  
  if(dim_transform == 1){avg_eval_matrix = array(0, dim = c(num_penalties, num_cps, num_means))}
  
  else{avg_eval_matrix = array(0, dim = c(num_penalties, dim_transform, num_cps, num_means))}
  
  for(i in 1 : num_penalties)
  {
    for(k in 1 : num_cps)
    {
      for(z in 1 : num_means)
      {
        m = array(0, dim = c(num_simulations, dim_transform))  
        
        for(j in 1 : num_simulations)
        {
          eval_instance = evaluation_matrix[j,, i, k, z]
          
          instance_transform = eval_transform(eval_instance)
          
          m[j,] = instance_transform  
        }
        
        avg_vec = colMeans(m)
        
        if(dim_transform == 1){avg_eval_matrix[i, k, z] = avg_vec}
        
        else{avg_eval_matrix[i,, k, z] = avg_vec}
      }
    }
  }
  
  return(avg_eval_matrix)
  
}

######################################
#####Data generating functions
######################################

gauss_data <- function(data_length, changepoints, means, variances)
{
  #A function to generate normal data with given parameter values and changepoints
  #change_points is the indexes at which the changepoints occur
  
  data = rnorm(changepoints[1]-1, means[1], variances[1])
  
  n = length(changepoints)
  
  if(n > 1)
  {
    for(i in 1 : length(changepoints) - 1)
    {
      new_data = rnorm(changepoints[i+1] - changepoints[i], means[i+1], variances[i+1])
    
      data = c(data, new_data)
    } 
  }
  
    final_data = rnorm(data_length + 1 - changepoints[length(changepoints)], means[length(means)], variances[length(variances)])
  
    data=c(data, final_data)
  
  
  return(data)
}

t_data <- function(data_length, changepoints, means, variances)
{
  #A function to generate normal data with given parameter values and changepoints
  #change_points is the indexes at which the changepoints occur
  
  data = rt(changepoints[1]-1, df = variances[1]) + means[1]
  
  n = length(changepoints)
  
  if(n > 1)
  {
    for(i in 1 : length(changepoints) - 1)
    {
      new_data = rt(changepoints[i+1] - changepoints[i], df = variances[i+1]) + means[i+1]
      
      data = c(data, new_data)
    } 
  }
  
  final_data = rt(data_length + 1 - changepoints[length(changepoints)], df = variances[length(variances)]) + means[length(means)]
  
  data=c(data, final_data)
  
  return(data)
}



cp_generator <- function(dist, params)
{
  segment_length = do.call(dist, params)
  
  num_cps = params$n
  
  cp_locations = rep(0 , num_cps)
  
  cp_locations[1] = segment_length[1]
  
  if(num_cps > 1){
  
    for(i in 1 : (num_cps - 1))
    {
      
      cp_locations[i+1] = segment_length[i+1] + cp_locations[i]
      
    }
  }
  
  return(cp_locations)
  
}

rate <- function(length, cps){return(floor(length / (cps + 1)))}


#Draws with 1 or -1 with equal probability, draws from unif[0.5, 1]
mean_sample <- function(diff, num_cps)
{
  up_down = sample(x = c(1, -1), size = num_cps + 1, replace = TRUE)
  
  mean_diff = diff * up_down
  
  final_means = rep(0, num_cps + 1)
  
  final_means[1] = mean_diff[1]
  
  for(i in 1 : num_cps)
  {
    final_means[i+1] = final_means[i] + mean_diff[i+1]
  }
  
  return(final_means)
}


tensor_index <- function(cp_instance, sim, param, pen)
{
  index_list = list()
  
  for(j in 1 : length(cp_instance))
  {
    cp = cp_instance[j]
    
    if(pen == -1){index = c(sim, cp, param)}
    else{index = c(sim, cp, param, pen)}
    
    index_list[[j]] = index
  }
  
  return(index_list)    
}

poisson_generator <- function(dist, params, max_length)
{
  
  segment_length = do.call(dist, params)
  
  num_cps = params$n
  
  cp_locations = segment_length[1]
  
  if(num_cps > 1){
    
    for(i in 1 : (num_cps - 1))
    {
      l = segment_length[i+1] + cp_locations[length(cp_locations)]
      
      if(l <= max_length){cp_locations = c(cp_locations, l)}
      
    }
  }
  
  return(cp_locations)
  
}

###################################
#Functions to call different algorithms
###################################

cp_storage <- function(cp_finder, data, penalties)
{
  #cp finder takes input of a dataset / penalty and returns the locations of changepoints
  #as identified by a given algoriithm
  #these locations are then stored as 1s in a sparse tensor
  
  tensor_dim = dim(data)
  tensor_dim = c(tensor_dim, length(penalties))
  
  num_pens = length(penalties)
  num_simulations = tensor_dim[1]
  
  #INITIALISE array
  stored_cps = array(0, dim = tensor_dim)
  
  
  for(j in 1 : num_pens)
  {
      for(i in 1 : num_simulations)
      {
        data_instance = data[i,]
        pen_instance = penalties[j]
        cp_instance = cp_finder(data_instance, pen_instance)
        
        #cp_finder to return -1 if no changepoints were found
        if(all(cp_instance != -1))
        {
          for(cp in cp_instance){stored_cps[i, cp, j] = 1}
        } 
      }
    }
  
  return(stored_cps)
}

library(changepoint.np)

NP_PELT_finder <- function(data_instance, penalty)
{
  cps = cpt.np(data = data_instance,
               method = 'PELT',
               test.stat = 'empirical_distribution',
               class = FALSE,
               penalty = 'Manual',
               nquantiles = 4*log(10000),
               pen.value = penalty)

  if(all(cps == length(data_instance))){cps = -1}

  else(cps = cps[cps != length(data_instance)])

  return(cps)

}

library(ecp)

KS_cp3o_finder <- function(data_instance, penalty)
{
  data_instance = matrix(data_instance)

  cps = ks.cp3o_delta(data_instance,
                      minsize = 1,
                      verbose = FALSE)$estimates

  if(all(cps == length(data_instance))){cps = -1}

  else(cps = cps[cps != length(data_instance)])

  return(cps)
}

library(wbs)

WBS_finder_250 <- function(data_instance, penalty)
{
  wbs_object = wbs(data_instance, M = 250)

  cps = try(changepoints(wbs_object, th = penalty), silent = TRUE)

  if(inherits(cps, "try-error")){cps = -1}

  else
  {
    cps = cps$cpt.th

    cps = cps[[1]]

    if(is.na(cps)){cps = - 1}
  }

  return(cps)
}


############################################
#Setting Experimental Parameters 
############################################

##Number of repetitions for the monte carlo simulations
num_simulations = 100

##Set data length
data_length = 10000

##Set the distribution of the data 
data_generator = t_data

init_dim = c(num_simulations, data_length)

segment_dist = rpois

##Set the 'correctness' threshold
threshold = 5

##Set the number of cps to have in the dataset on average
cps_range = c() 

##Set the absolute value of the difference in mean over each segment
mean_differences = c() 

##Set range of penalties to perform grid search over
experiment_penalties = c()

############################################
##Experiment - adjust the change detection algorithm used
############################################

for(diff in mean_differences)
{
        
  for(num_cps in cps_range)
  {
  
    l = rate(length = data_length, cps = num_cps)
  
    cp_params = list("n" = num_cps, "lambda" = l)
    
    experiment_data = array(0, dim = init_dim)
    
    experiment_actual_cps = array(0, dim = init_dim)
    

    for(i in 1 : num_simulations)
    {
      set.seed(51621 * i )
    
      cp_instance = poisson_generator(segment_dist, params = cp_params, max_length = data_length)
    
      for(cp in cp_instance){experiment_actual_cps[i, cp] = 1}
      
      
      ##Var only changes
      # means_instance = rep(0, num_cps + 1)
      # variances_instance = variance_sample(var_differences, num_cps + 1)
      
      ##Simultaneous mean and var changes
      # mean_var_instance = offsetting_mean_variance_sample(mean_differences[j], var_differences[j], num_cps + 1)
      # means_instance = mean_var_instance$means
      # variances_instance = mean_var_instance$variances
      
      ##Mean only, var only, or simultaneous mean and var
      # mean_var_instance = mean_variance_sample(mean_differences[j], var_differences[j], num_cps + 1)
      # means_instance = mean_var_instance$means
      # variances_instance = mean_var_instance$variances
      
      ##mean only changes
      means_instance = mean_sample(diff, num_cps)
      variances_instance = rep(4, num_cps + 1)
      
      data_params = list("data_length" = data_length, "changepoints" = cp_instance, "means" = means_instance, 
                       "variances" = variances_instance)
    
      data_instance = do.call(data_generator, data_params)
    
      experiment_data[i,] = data_instance
    
    }
    
    for(pen in experiment_penalties)
    {
      
    eval_params = list("data_length" = data_length, "threshold" = threshold, "actual_cps" = experiment_actual_cps)
  
    t1 = Sys.time()

    NPL_BS_est_cps = cp_storage(cp_finder = NPL_BS,
                                penalties = pen,
                                data = experiment_data)
    
    
    NPL_BS_confusion_vec = evaluation_matrix(experiment_changepoints = NPL_BS_est_cps,
                                             evaluation_method = confusion_matrix,
                                             evaluation_params = eval_params)
    
    NPL_BS_MSE_vec = evaluation_matrix(experiment_changepoints = NPL_BS_est_cps,
                                       evaluation_method = mean_square_error,
                                       evaluation_params = list("actual_cps" = experiment_actual_cps))
    
    # WBS_est_cps = cp_storage(cp_finder = WBS_finder,
    #                                penalties = pen,
    #                                data = experiment_data)
    #  
    # WBS_confusion_vec = evaluation_matrix(experiment_changepoints = WBS_est_cps,
    #                                             evaluation_method = confusion_matrix,
    #                                             evaluation_params = eval_params)
    #  
 
    NPL_BS_confusion_vec_file = paste('NPL - BS - mean diff', diff, 'num_cps', num_cps, 'penalty', pen, 'confusion vec.RData')
    save(NPL_BS_confusion_vec, file = NPL_BS_confusion_vec_file)  

    NPL_BS_mse_vec_file = paste('NPL - BS - mean diff', diff, 'num_cps', num_cps, 'penalty', pen, 'mse vec.RData')
    save(NPL_BS_MSE_vec, file = NPL_BS_mse_vec_file)
    
    t2 = Sys.time()
    print(t2-t1)
    gc()
    
    }  
  }
}

############################################
###function to retreive saved data files from a given directory
############################################

# loadRData <- function(fileName){
#   #loads an RData file, and returns it
#   load(fileName)
#   get(ls()[ls() != "fileName"])
# }
# 
# 
# #directory = "\\ \\"
# 
# setwd(directory)
# 
# #range of parameters which the experiment was run over 
# cps_range = c(0) 
# 
# mean_differences = (0) 
# 
# experiment_penalties = c(0)
# 
# 
# #initialising evaluation matrices
# 
# num_simulations = 100
# 
# dim_confusion = 4
# 
# confusion_matrix_dim = c(num_simulations, dim_confusion, length(experiment_penalties), length(cps_range), length(mean_differences))
# 
# WBS_5000_t2_confusion = array(0, dim = confusion_matrix_dim)
# 
# #loading the data files and converting to the desired format
# 
# for(i in 1 : length(mean_differences))
# {
#   diff = mean_differences[i]
#   
#   for(j in 1 : length(cps_range))
#   {
#     num_cps = cps_range[j]
#     
#     for(k in 1 : length(experiment_penalties))
#     {
#       pen = experiment_penalties[k]
#       
#       #setting up file name
#       confusion_vec_file = paste('WBS 5000 - mean diff', diff, 'num_cps', num_cps, 'penalty', pen, 'confusion vec.RData')
#       
#       #loading data 
#       confusion_vec_instance = loadRData(paste(directory, confusion_vec_file, sep=''))
#       
#       #storing in the final matrix
#       WBS_5000_t2_confusion[,,k,j,i] = confusion_vec_instance
#       
#     }
#   }
# }


###################################
#Functions to call different algorithms 2
###################################

#library(changepoint)
# PELT_mean_finder <- function(data_instance, penalty)
# {
#   cps = cpt.mean(data = data_instance,
#                  method = 'PELT',
#                  test.stat = 'Normal',
#                  class = FALSE,
#                  penalty = 'Manual',
#                  pen.value = penalty)
#   
#   if(all(cps == length(data_instance))){cps = -1}
#   
#   else(cps = cps[cps != length(data_instance)])
#   
#   return(cps)
#   
# }
# 
# PELT_var_finder <- function(data_instance, penalty)
# {
#   cps = cpt.var(data = data_instance,
#                 method = 'PELT',
#                 test.stat = 'Normal',
#                 class = FALSE,
#                 penalty = 'Manual',
#                 pen.value = penalty)
#   
#   if(all(cps == length(data_instance))){cps = -1}
#   
#   else(cps = cps[cps != length(data_instance)])
#   
#   return(cps)
#   
# }
# 
# 
# 
# PELT_mean_var_finder <- function(data_instance, penalty)
# {
#   cps = cpt.meanvar(data = data_instance,
#                     method = 'PELT',
#                     test.stat = 'Normal',
#                     class = FALSE,
#                     penalty = 'Manual',
#                     pen.value = penalty)
#   
#   if(all(cps == length(data_instance))){cps = -1}
#   
#   else(cps = cps[cps != length(data_instance)])
#   
#   return(cps)
#   
# }
# 
# BS_var_finder <- function(data_instance, penalty)
# {
#   cps = cpt.var(data = data_instance,
#                 method = 'BinSeg',
#                 test.stat = 'Normal',
#                 class = FALSE,
#                 penalty = 'Manual',
#                 pen.value = penalty,
#                 Q = 10000)
#   
#   if(all(cps == length(data_instance))){cps = -1}
#   
#   else(cps = cps[cps != length(data_instance)])
#   
#   return(cps)
#   
# }
# 
# BS_mean_var_finder <- function(data_instance, penalty)
# {
#   cps = cpt.meanvar(data = data_instance,
#                     method = 'BinSeg',
#                     test.stat = 'Normal',
#                     class = FALSE,
#                     penalty = 'Manual',
#                     pen.value = penalty,
#                     Q = 10000)
#   
#   if(all(cps == length(data_instance))){cps = -1}
#   
#   else(cps = cps[cps != length(data_instance)])
#   
#   return(cps)
#   
# }
# 
# 
# BS_norm_mean_finder <- function(data_instance, penalty)
# {
#   cps = cpt.mean(data = data_instance,
#                  method = 'BinSeg',
#                  test.stat = 'Normal',
#                  class = FALSE,
#                  penalty = 'Manual',
#                  pen.value = penalty,
#                  Q = 10000)
#   
#   if(all(cps == length(data_instance))){cps = -1}
#   
#   else(cps = cps[cps != length(data_instance)])
#   
#   return(cps)
#   
# }
# 
# BS_CUSUM_finder <- function(data_instance, penalty)
# {
#   wbs_object = sbs(data_instance)
#   
#   cps = changepoints(wbs_object, th = penalty)
#   
#   cps = cps$cpt.th
#   
#   cps = cps[[1]]
#   
#   if(is.na(cps)){cps = - 1}
#   #cps = as.vector(cps)
#   
#   #cp_vec = as.vector(cps)
#   
#   return(cps)
# }
