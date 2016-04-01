vertices_m = NULL
vertices_i = NULL
vertices_d = NULL

comp_emis_prob <- function(vec, letter) {
  
  sum = 0
  no_let = length(vec)
  
  
    for (v in vec) {
      if (v == letter) {
        sum = sum + 1
      }
    }
  
  sum/no_let
}



emis_prob <-function(vec) {
  c(comp_emis_prob(vec,"A"), comp_emis_prob(vec,"T"), comp_emis_prob(vec,"C"),comp_emis_prob(vec,"G"))
}

add_vertex <-function(vertices, v, num, vec_data) {

  st = paste(v, toString(num), sep="")
  
  if(v == "d" && !is.element(st, vertices)) {
    vertices = c(vertices, st)
  } else if(v == "i" && !is.element(st, vertices)) {
    vertices = cbind(vertices,c( st, comp_emis_ins_wrapper(vec_data)))
  }
  else if(!is.element(st, vertices)) {
    em_vec1 = emis_prob(vec_data)
    vertices = cbind(vertices, c(st, em_vec1))
  }

  vertices
}

make_edge <- function(st1, n1, st2, n2, tr_prob) {

  res = c(paste(st1, toString(n1), sep=""), 
            paste(st2, toString(n2), sep=""), tr_prob)
  
  res
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
# delete to insertion state
comp_del_to_ins <- function(vec, insert_mat ) {
  
  br = cbind(vec,insert_mat)
  
  found_fl = 0
  no_letters = 0
  
  for (i in 1:length(vec)) {
      if(br[i,1] == "-") {
        
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
  
  vertices_m = NULL
  vertices_d = NULL
  vertices_i = NULL
  
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
            
            edges = cbind(edges, make_edge("m", lmn, "i", lmn, comp_mch_to_ins(matrixData[,i], insert_mat)))
            edges = cbind(edges, make_edge("d", lmn, "i", ldn, comp_del_to_ins(matrixData[,i], insert_mat)))
            lin = i+1
        }
        else if(dashes[i] == 0 ) {  
            edges = cbind(edges, make_edge("m", lmn, "i", lmn, comp_mch_to_ins(matrixData[,i], insert_mat)))
            lin = lmn
        }
        
        vertices_i = add_vertex(vertices_i, "i", lmn, insert_mat)
    } 
    else if(dashes[i+1] < 3 && dashes[i+1] > 0) {
    
        if(dashes[i] >= 3 && i > 1) { #if previous state is insertion

          edges = cbind(edges, make_edge("i", lin, "m", i+1, comp_ins_to_mch(insert_mat, matrixData[,i+1])))
          edges = cbind(edges, make_edge("i", lin, "d", i+1, comp_del_to_ins(matrixData[,i+1],insert_mat)))
          edges = cbind(edges, make_edge("m", lmn, "m", i+1, comp_mch_ins_mch(matrixData[,lmn],insert_mat, matrixData[,i+1])))
          
          vertices_m = add_vertex(vertices_m, "m", i+1, matrixData[,i+1])
          vertices_d = add_vertex(vertices_d, "d", i+1, matrixData[,i+1])
          
          lmn = i+1
          ldn = i+1
          next;
          
          ins_mat=NULL
        }
        else if(dashes[i] < 3 && dashes[i] > 0) {
            
            edges = cbind(edges, make_edge("m", lmn, "m", i+1, tr_prob))  
            edges = cbind(edges, make_edge("d", ldn,"m", i+1, tr_prob))
            edges = cbind(edges, make_edge("d", ldn,"d", i+1, tr_prob))
            
           # lmn = i + 1
            #ldn = lmn
        } 
        else if(dashes[i] == 0) {
            edges = cbind(edges, make_edge("m", lmn, "m", i+1, tr_prob))
            
        }
        edges = cbind(edges, make_edge("m", lmn, "d", i+1, tr_prob))
        
        vertices_m = add_vertex(vertices_m, "m", i+1, matrixData[,i+1])
        vertices_d = add_vertex(vertices_d, "d", i+1, matrixData[,i+1])
        
        lmn = i+1
        ldn = i+1
    } 
    else if(dashes[i+1] == 0) {
        if(dashes[i] >= 3 && i > 1) { #if previous state is insertion
            edges = cbind(edges, make_edge("i", lin, "m", i+1, comp_ins_to_mch(insert_mat, matrixData[,i+1])))
            lmn = i+1
        }
        else if(dashes[i] < 3 && dashes[i] > 0) {
            edges = cbind(edges, make_edge("m", lmn, "m", i+1, tr_prob))
            edges = cbind(edges, make_edge("d", ldn,"m", i+1, tr_prob))
            lmn = i+1
        }
        else if(dashes[i] == 0) {
            edges = cbind(edges, make_edge("m", lmn, "m", i+1, tr_prob)) 
            lmn = i+1
        }
      
        vertices_m = add_vertex(vertices_m, "m", i+1, matrixData[,i+1])
    } 
  }
  list(edges, vertices_i, vertices_d, vertices_m)
}


comp_emis_ins <- function(vec, letter) {
  
  sum = 0
  no_let = length(vec)
  
  for (i in 1:length(vec[,1])) {
      for (j in 1:length(vec[1,])) {
        if(vec[i,j] == letter) {
          sum = sum + 1
          no_let = no_let + 1
        } else if(vec[i,j] != "-") {
          no_let = no_let + 1
        }
      }
    }
  sum/no_let
}

comp_emis_ins_wrapper <- function(vec) {
  
  c(comp_emis_ins(vec,"A"), comp_emis_ins(vec,"T"), comp_emis_ins(vec,"C"),comp_emis_ins(vec,"G"))
  
}


comp_dashes <-function(mat_data) {
  
  count_dashes = 1:length(mat_data[1,]);
  
  count_dash = 0;
  for (i in 1:length(mat_data[1,])) {
    for (el in mat_data[,i]) {
      
      if(el == '-') {
        count_dash = count_dash + 1; 
      }
      
    }
    count_dashes[i] = count_dash
    count_dash = 0
  }
  
  count_dashes
}


