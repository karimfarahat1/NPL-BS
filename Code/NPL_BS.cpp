#include <Rcpp.h>
using namespace Rcpp;


std::set<int> changepoint_scan(const std::vector<double> &time_sorted_segment,
                               std::map<double, int> &rank_sorted_segment,
                               const double* log_lookup,
                               const double* log_diff_lookup,
                               const double* log_sum_lookup,
                               const double &penalty,
                               const int changepoint_offset,
                               std::set<int> &changepoint_list)
{
  ///Initialisation
  int current_cp_index = 1;
  int max_cp_index = current_cp_index;
  
  std::vector<double>::const_iterator iter = time_sorted_segment.begin()+1;
  int unseg_data_length = time_sorted_segment.size();
  int alternate_seg_length = unseg_data_length - 1;
  
  //Initial Rank
  int right_rank_iter = rank_sorted_segment[time_sorted_segment[0]];
  int left_rank_iter = 1;
  int alternate_seg_rank = right_rank_iter - 1;
  
  //Initial Changepoint Cost
  double current_cp_cost = log_lookup[alternate_seg_rank - 1]  +
    log_lookup[alternate_seg_length - alternate_seg_rank - 1] +
    2 * log_sum_lookup[alternate_seg_length - 2] -
    (unseg_data_length - 1) *  (log_lookup[alternate_seg_length - 1]);
  
  double max_cp_cost = current_cp_cost;
  
  //Initialising binary search trees & iterators
  std::map<double, int> current_left_seg;
  current_left_seg.insert(std::make_pair(time_sorted_segment[0], right_rank_iter));
  std::map<double, int> max_cp_left_seg = current_left_seg;
  std::map<double, int> current_right_seg = rank_sorted_segment;
  current_right_seg.erase(time_sorted_segment[0]);
  std::map<double, int> max_cp_right_seg = current_right_seg;
  std::map<double, int>::iterator right_map_iter;
  std::map<double, int>::iterator left_map_iter;
  std::map<double, int>::iterator map_iter;
  
  ///Main loop to find the changepoint candidate
  for(; iter != time_sorted_segment.end()-2; ++iter)
  {
    ++current_cp_index;
    --alternate_seg_length;
    
    //Updating segment with new point & preparing iters
    current_left_seg.insert(std::make_pair(*iter, rank_sorted_segment[*iter]));
    left_map_iter = current_left_seg.find(*iter);
    left_rank_iter = 0;
    
    map_iter = current_left_seg.begin();
    
///Computing cost associated with the left tree
    
  //Ranks of points unchanged - likelihood update required due to segment length
    for(; map_iter != left_map_iter; ++map_iter)
    {
      //Rank of point within own segment
      ++left_rank_iter;
      
      //Rank of point within alternate segment
      alternate_seg_rank = (*map_iter).second - left_rank_iter;
      
      //Partial cost associated with point
      current_cp_cost -= log_diff_lookup[alternate_seg_length - alternate_seg_rank];
    }
    
    //Cost update of the point just switched
    ++left_rank_iter;
    alternate_seg_rank = (*left_map_iter).second - left_rank_iter;
    
    //As a special case - difference cost of this point directly
    current_cp_cost += log_lookup[alternate_seg_rank - 1] +
      log_lookup[alternate_seg_length - alternate_seg_rank - 1] -
      log_lookup[left_rank_iter - 2] -
      log_lookup[current_cp_index - left_rank_iter - 1];
    
    ++left_map_iter;
    map_iter = current_left_seg.end();
    
    ++left_rank_iter;
    alternate_seg_rank = (*left_map_iter).second - left_rank_iter;
    
  //Ranks of points updated by 1 - likelihood update required due to this 
    //While loop ensures that ranks do not exceed the number of points within the alternate segment
    while((alternate_seg_length - alternate_seg_rank > 0) & (left_map_iter != map_iter))
    {
      //Partial cost associated with point
      current_cp_cost -= log_diff_lookup[alternate_seg_rank];
      
      //Rank of point within own segment
      ++left_rank_iter;
      ++left_map_iter;
      
      //Rank of point within alternate segment
      alternate_seg_rank = (*left_map_iter).second - left_rank_iter;
    }
    
    right_rank_iter = 0;
    right_map_iter = current_right_seg.upper_bound(*iter);
    current_right_seg.erase(*iter);
    
    map_iter = current_right_seg.begin();
    
///Computing cost associated with the right tree
  
  //Ranks of points unchanged - likelihood update required due to segment length
    for(; map_iter != right_map_iter; ++map_iter)
    {
      //Rank of point within own segment
      ++right_rank_iter;
      
      //Rank of point within alternate segment
      alternate_seg_rank = (*map_iter).second - right_rank_iter;
      
      //Partial cost associated with point
      current_cp_cost += log_diff_lookup[current_cp_index - alternate_seg_rank - 1];
    }
    
    //Iterator triggers change of update
    map_iter = current_right_seg.end();
    ++right_rank_iter;
    alternate_seg_rank = (*right_map_iter).second - right_rank_iter;
    
  //Ranks of points updated by 1 - likelihood update required due to this 
    //While loop ensures that ranks do not exceed the number of points within the alternate segment
    while((current_cp_index - alternate_seg_rank > 0) & (right_map_iter != map_iter))
    {
      //Partial cost associated with point
      current_cp_cost += log_diff_lookup[alternate_seg_rank - 1];
      
      //Rank within own segment
      ++right_rank_iter;
      ++right_map_iter;
      //Rank of point within alternate segment
      alternate_seg_rank = (*right_map_iter).second - right_rank_iter;
    }
    
    //Total changepoint cost:
      //including data independent terms
      //points with rank greater than the maximum of the alternate segment   
    current_cp_cost += 2 * (log_lookup[current_cp_index - 2] -
      log_lookup[alternate_seg_length - 1]) -
      (unseg_data_length - 1) * (log_diff_lookup[current_cp_index - 1] -
      log_diff_lookup[alternate_seg_length])
      - (current_cp_index - left_rank_iter + 1) * log_diff_lookup[alternate_seg_length - 1]
    + (alternate_seg_length - right_rank_iter + 1) * log_diff_lookup[current_cp_index - 1];
    
    //Tracking maximum changepoint candidate
    if(current_cp_cost - max_cp_cost > 0)
    {
      max_cp_index = current_cp_index;
      max_cp_cost = current_cp_cost;
      max_cp_left_seg = current_left_seg;
      max_cp_right_seg = current_right_seg;
    }
  }
  
  //Including baseline and weighting
  max_cp_cost -= 2 * log_sum_lookup[unseg_data_length - 2] - ((unseg_data_length - 1))*log_lookup[unseg_data_length - 1];
  max_cp_cost *= (1 / double(unseg_data_length));
  
///Recursion stage
  
  //Checking significance of maximum
  if(max_cp_cost - penalty > 0)
  {
    //Add changepoint
    changepoint_list.insert(max_cp_index + changepoint_offset);
    
    //Check minimum segment size of 2
    if(unseg_data_length - max_cp_index > 2)
    {
      //Initialise tree
      std::map<double, int> max_cp_right_ranks;
      right_rank_iter = 0;
      
      //Update global tree for use in next iteration
      for(right_map_iter = max_cp_right_seg.begin(); right_map_iter != max_cp_right_seg.end(); ++right_map_iter)
      {
        ++right_rank_iter;
        max_cp_right_ranks.insert(std::make_pair((*right_map_iter).first, right_rank_iter));
      }
      
      //Implementing recursion
      std::set<int> right_changepoints = changepoint_scan(std::vector<double>(time_sorted_segment.begin() + max_cp_index, time_sorted_segment.end()),
                                                          max_cp_right_ranks,
                                                          log_lookup,
                                                          log_diff_lookup,
                                                          log_sum_lookup,
                                                          penalty,
                                                          changepoint_offset + max_cp_index,
                                                          changepoint_list);
    }
    
    //Check minimum segment size of 2
    if(max_cp_index > 2)
    {
      //Initialise tree
      std::map<double, int> max_cp_left_ranks;
      left_rank_iter = 0;
      
      //Update global tree for use in next iteration
      for(left_map_iter = max_cp_left_seg.begin(); left_map_iter != max_cp_left_seg.end(); ++left_map_iter)
      {
        ++left_rank_iter;
        max_cp_left_ranks.insert(std::make_pair((*left_map_iter).first, left_rank_iter));
      }
      
      //Implementing recursion
      std::set<int> left_changepoints = changepoint_scan(std::vector<double>(time_sorted_segment.begin(), time_sorted_segment.begin() + max_cp_index),
                                                         max_cp_left_ranks,
                                                         log_lookup,
                                                         log_diff_lookup,
                                                         log_sum_lookup,
                                                         penalty,
                                                         changepoint_offset,
                                                         changepoint_list);
    }
    
    return changepoint_list;
  }
  else
  {
    return changepoint_list;
  }
}



// [[Rcpp::export]]

NumericVector NPL_BS(const NumericVector data, const double penalty)
{
  std::set<int> changepoints;
  std::vector<double> time_ordered_data = as< std::vector<double> >(data);
  std::map<double, int> rank_ordered_data;
  std::vector<double>::iterator iter;
  int rank_iter = 0;
  
  double log_sum_lookup[time_ordered_data.size()] = {};
  double log_lookup[time_ordered_data.size()] = {};
  double log_diff_lookup[time_ordered_data.size()] = {};
  
  //Creating first order stat tree
  std::vector<double> sorted_data = time_ordered_data;
  std::sort(sorted_data.begin(), sorted_data.end());
  
  for(iter = sorted_data.begin(); iter != sorted_data.end(); ++iter)
  {
    //Build tree
    rank_iter +=1;
    rank_ordered_data.insert(std::make_pair(*iter, rank_iter));
    
    //Storing log terms
    log_lookup[rank_iter - 1] = rank_iter * logf(rank_iter);
    log_diff_lookup[rank_iter - 1] = log_lookup[rank_iter - 1 ] - log_lookup[rank_iter - 2];
    log_sum_lookup[rank_iter - 1] = log_sum_lookup[rank_iter - 2] + log_lookup[rank_iter - 1];
    
  }
  
  
  int offset = 0;
  
  //Calling c++ function
  changepoints = changepoint_scan(time_ordered_data,
                                  rank_ordered_data,
                                  log_lookup,
                                  log_diff_lookup,
                                  log_sum_lookup,
                                  penalty,
                                  offset,
                                  changepoints);
  
  NumericVector final_cps = wrap(changepoints);
  
  return final_cps;
}



