make_edge <- function(st1, n1, st2, n2, tr_prob) {
  res = c(paste(st1, toString(n1), sep=""), 
            paste(st2, toString(n2), sep=""), tr_prob)
  
  res
}

# Constructs the edge between the states
append_edge <- function(st1,n1,st2,n2, edges) {
  
  if(is.null(edges)) {
    edges = make_edge(st1,n1,st2,n2)
  }
  else {
    edges = cbind(edges, make_edge(st1,n1,st2,n2))
  }
  
  edges
}

# Computes the transition probability from 
# insertion to match state
comp_ins_to_mch <-function (ins_mat, vec) {
  
  br = cbind(vec,ins_mat)
  no_letters = 0
  
  for (i in 1:length(vec)) {
      if(br[i,1] != "-") {
        for (j in 2:length(br[1,])) {
          if(br[i,j] != "-") {
            no_letters = no_letters + 1
            break
          }
        }
      }
  }
  no_letters
}

comp_mch_ins_mch <- function(vec1, ins_mat, vec2) {
  
  br = cbind(vec1,ins_mat, vec2)
  no_letters = 0
  
  for (i in 1:length(br[,1])) {
    no_mch_fl = 0
    for (j in 2:(length(br[1,])-1)) {
      if(br[i,j] != "-") {
        no_mch_fl = 1
        break;
      }
    }
    
    if(no_mch_fl == 0 && vec1[i] != "-" && vec2[i] != "-") {
      no_letters = no_letters + 1
    }
  }
  
  no_letters
}


# Computes the transition probability from 
# insertion to insertion state
comp_ins_to_ins <-function (ins_mat) {
  
  no_letters = 0
  
  for (i in 1:length(ins_mat[,1])) {
      if(ins_mat[i,1] != "-") {
        for (j in 2:length(ins_mat[1,])) {
          if(ins_mat[i,j] != "-") {
            no_letters = no_letters + 1
            found_fl = 1
            break
          }
        }
      }
  }
  no_letters
}

# Computes the transition probability from 
# match to insertion state
comp_mch_to_ins <- function(vec, insert_mat ) {
  
  br = cbind(vec,insert_mat)
  
  found_fl = 0
  no_letters = 0
  
  for (i in 1:length(vec)) {
      if(br[i,1] != "-") {
        
        for (j in 2:length(br[1,])) {
          
          if(br[i,j] != "-") {
            
            no_letters = no_letters + 1
            
            found_fl = 1
            break
          }
        }
      }
  }
  no_letters
}



profile_hmm2 <- function(matrixData, dashes) {
  

  edges = NULL
  lin = 0 # last insertion number
  lmn = 0 # last match number
  ldn = 0 # last delete number
  
  ins_flag = 0
  insert_mat = NULL
  
  for (i in 1:(length(matrixData[1,])-1)) {

    tr_prob = (length(matrixData[,1]) - dashes[i])/length(matrixData[,1])
        
    if(dashes[i+1] >= 3) { # if current is insertion state
    
        ii = i
        while(dashes[ii+1] >= 3) {
          
            insert_mat = cbind(insert_mat,matrixData[ ,ii+1])
          
          ii= ii+1
        }
        
        if(dashes[i] >= 3 && i > 1) {} 
        else if(dashes[i] < 3 && dashes[i] > 0) {
          
            
          
            edges = cbind(edges, make_edge("m", lmn, "i", lin+1, tr_prob))
            edges = cbind(edges, make_edge("d", lmn, "i", lin+1, tr_prob))
            lin = lin + 1
          
        }
        else if(dashes[i] == 0 ) {
          
          # if(is.null(edges)) {
          #   edges = make_edge("m",lmn,"i", lin+1)
          #   lin = lin + 1
          # } else {
            edges = cbind(edges, make_edge("m", lmn, "i", lin+1, tr_prob))
            lin = lin + 1
          # }
        }
    } 
    else if(dashes[i+1] < 3 && dashes[i+1] > 0) {
    
        if(dashes[i] >= 3 && i > 1) { #if previous state is insertion
          
          edges = cbind(edges, make_edge("i", lin, "m", lmn+1, tr_prob))
          edges = cbind(edges, make_edge("i", lin, "d", ldn+1, tr_prob))
        }
        else if(dashes[i] < 3 && dashes[i] > 0) {
          # if(is.null(edges)) {
          #   edges = make_edge("m",lmn,"m", lmn+1)
          # } else if (!is.null(edges)) {
            edges = cbind(edges, make_edge("m", lmn, "m", lmn+1, tr_prob))
          # }
          edges = cbind(edges, make_edge("d", ldn,"m", lmn+1, tr_prob))
          edges = cbind(edges, make_edge("d", ldn,"d", ldn+1, tr_prob))
        } 
        else if(dashes[i] == 0) {
          
          # if(is.null(edges)) {
          #   edges = make_edge("m",lmn,"m", lmn+1)
          # } else {
            edges = cbind(edges, make_edge("m", lmn, "m", lmn+1, tr_prob))
          # }
        }
        
        edges = cbind(edges, make_edge("m", lmn, "d", ldn+1, tr_prob))
        lmn = lmn + 1
        ldn = ldn + 1
    } 
    else if(dashes[i+1] == 0) {
        if(dashes[i] >= 3 && i > 1) { #if previous state is insertion
            edges = cbind(edges, make_edge("i", lin, "m", lmn+1, tr_prob))
            lmn = lmn + 1
        }
        else if(dashes[i] < 3 && dashes[i] > 0) {
            # if(is.null(edges)) {
            #   edges = make_edge("m",lmn,"m", lmn+1)
            # } 
            # else {
              edges = cbind(edges, make_edge("m", lmn, "m", lmn+1, tr_prob))
            # }
            edges = cbind(edges, make_edge("d", ldn,"m", lmn+1, tr_prob))
            lmn = lmn + 1
        }
        else if(dashes[i] == 0) {
          # if(is.null(edges)) {
          #   edges = make_edge("m",lmn,"m", lmn+1)
          # } 
          # else {
            edges = cbind(edges, make_edge("m", lmn, "m", lmn+1, tr_prob)) 
          # }
          lmn = lmn + 1
        }
    } 
  }
  edges
}