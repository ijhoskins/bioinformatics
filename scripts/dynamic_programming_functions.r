# Dynamic programming functions

# Note: the paradigm in the LCS function can be used to solve other dp problems
# the only difference should be the scoring functions used and the determination of the max predecessors

# Longest common subsequence between two strings
# Subsequences need not be consecutive matches in the two strings (like substrings must be)
lcs<- function(v,w){
  
  # Plus 1 for the dpt because we need an additional row and column
  len_v<- nchar(v) + 1
  len_w<- nchar(w) + 1
  
  # Init a dynamic programming table and a backtracking matrix
  dpt<- matrix(data=0, nrow=len_v, ncol=len_w)
  backtrack_matrix<- matrix(nrow=len_v, ncol=len_w)
  
  # Fill in the dpt and store path in the backtracking matrix
  for(i in 2:len_v){
    for(j in 2:len_w){
      
      # Determine the score of the predecessors to the vertex in question
      top_pred<- dpt[i - 1, j]
      left_pred<- dpt[i, j - 1]
      diag_pred<- dpt[i - 1, j - 1]
      
      # If we have a match, add 1 to the diag
      if(substr(v, i - 1, i - 1) == substr(w, j - 1, j - 1)){
        diag_pred<- diag_pred + 1
      }
      
      # Determine the highest scoring predecessor, and store in the table
      all_preds<- c(top_pred, left_pred, diag_pred)
      best_pred_index<- which.max(all_preds)
      max_pred<- all_preds[best_pred_index]
      dpt[i,j]<- max_pred
      
      # Also store some backtracking pointers to recreate the best path
      if(max_pred == top_pred){
        backtrack_matrix[i,j]<- "top"
      } else if (max_pred == left_pred) {
        backtrack_matrix[i,j]<- "left"
      } else {
        backtrack_matrix[i,j]<- "diag"
      }
    }
  }
  
  # Annotate the backtrack with the sequence strings
  v_list<- unlist(strsplit(v, ""))
  c_list<- unlist(strsplit(w, ""))
  backtrack_matrix<- backtrack_matrix[-1,-1]  # remove the unnecessary row and column
  rownames(backtrack_matrix)<- v_list
  colnames(backtrack_matrix)<- c_list
  
  # Finally compute some results
  lcs<- dpt[len_v, len_w]  # the LCS is just the score of the corner of the dpt
  results<- list(LCS=lcs, Backtrack=backtrack_matrix)
  return(results)
}

