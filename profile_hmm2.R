make_edge <- function(st1, n1, st2, n2) {
  res = c(paste(st1, toString(n1), sep=""), 
            paste(st2, toString(n2), sep=""))
  
  res
}

append_edge <- function(st1,n1,st2,n2, edges) {
  
  if(is.null(edges)) {
    edges = make_edge(st1,n1,st2,n2)
  }
  else {
    edges = cbind(edges, make_edge(st1,n1,st2,n2))
  }
  
  edges
}

profile_hmm2 <- function(matrixData, dashes) {
  

  edges = NULL
  lin = 0 # last insertion number
  lmn = 0 # last match number
  ldn = 0 # last delete number
  
  ins_flag = 0
  insert_mat = NULL
  
  for (i in 1:(length(matrixData[1,])-1)) {
    
    if(dashes[i+1] >= 3) { # if current is insertion state
    
        ii = i
        while(dashes[ii+1] >= 3) {
          if(is.null(insert_mat)) {
            insert_mat = matrixData[ ,ii+1]
          } else {
            insert_mat = cbind(insert_mat,matrixData[ ,ii+1])
          }
          ii= ii+1
        }
        
        if(dashes[i] >= 3 && i > 1) {} 
        else if(dashes[i] < 3 && dashes[i] > 0) {
          if(is.null(edges)) {
            edges = make_edge("m",lmn,"i", lin+1)
            lin = lin + 1
          } else {
            edges = cbind(edges, make_edge("m", lmn, "i", lin+1))
            edges = cbind(edges, make_edge("d", lmn, "i", lin+1))
            lin = lin + 1
          }
        }
        else if(dashes[i] == 0 ) {
          
          if(is.null(edges)) {
            edges = make_edge("m",lmn,"i", lin+1)
            lin = lin + 1
          } else {
            edges = cbind(edges, make_edge("m", lmn, "i", lin+1))
            lin = lin + 1
          }
        }
    } 
    else if(dashes[i+1] < 3 && dashes[i+1] > 0) {
    
        if(dashes[i] >= 3 && i > 1) { #if previous state is insertion
          edges = cbind(edges, make_edge("i", lin, "m", lmn+1))
          edges = cbind(edges, make_edge("i", lin, "d", ldn+1))
        }
        else if(dashes[i] < 3 && dashes[i] > 0) {
          if(is.null(edges)) {
            edges = make_edge("m",lmn,"m", lmn+1)
          } else if (!is.null(edges)) {
            edges = cbind(edges, make_edge("m", lmn, "m", lmn+1))
          }
          edges = cbind(edges, make_edge("d", ldn,"m", lmn+1))
          edges = cbind(edges, make_edge("d", ldn,"d", ldn+1))
        } 
        else if(dashes[i] == 0) {
          
          if(is.null(edges)) {
            edges = make_edge("m",lmn,"m", lmn+1)
          } else {
            edges = cbind(edges, make_edge("m", lmn, "m", lmn+1))
          }
        }
        
        edges = cbind(edges, make_edge("m", lmn, "d", ldn+1))
        lmn = lmn + 1
        ldn = ldn + 1
    } 
    else if(dashes[i+1] == 0) {
        
        if(dashes[i] >= 3 && i > 1) { #if previous state is insertion
            edges = cbind(edges, make_edge("i", lin, "m", lmn+1))
            lmn = lmn + 1
        }
        else if(dashes[i] < 3 && dashes[i] > 0) {
            if(is.null(edges)) {
              edges = make_edge("m",lmn,"m", lmn+1)
            } 
            else {
              edges = cbind(edges, make_edge("m", lmn, "m", lmn+1))
            }
            edges = cbind(edges, make_edge("d", ldn,"m", lmn+1))
            lmn = lmn + 1
        }
        else if(dashes[i] == 0) {
          if(is.null(edges)) {
            edges = make_edge("m",lmn,"m", lmn+1)
          } 
          else {
            edges = cbind(edges, make_edge("m", lmn, "m", lmn+1)) 
          }
          lmn = lmn + 1
        }
    } 
    
    
  }
  edges
}