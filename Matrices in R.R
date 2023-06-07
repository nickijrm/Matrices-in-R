
# This project aims to create some simple matrix-related functions for R


#Returns spectral radius of a given matrix
spectral_radius = function(M){
  largest_eigenvalue = max(eigen(M)$values)
  smallest_eigenvalue = min(eigen(M)$values)
  max_eigenvalue = max(abs(largest_eigenvalue),abs(smallest_eigenvalue))
  return(max_eigenvalue)
}

#Returns spectral abscissa of a given matrix (largest real-part of all eigenvalues)
spectral_abscissa = function(M){
  if(class(eigen(M)$values)==complex){
    X= as.numeric(eigen(M)$values)
    largest_eigenvalue = max(X)
    smallest_eigenvalue = min(X)
    max_eigenvalue = max(abs(largest_eigenvalue),abs(smallest_eigenvalue))
    return(max_eigenvalue)
  }
  else{
    return(spectral_radius(M))
  }
}

#Based on the theorem that a matrix is irreducible iff it corresponds to a strongly connected graph
#Returns TRUE, FALSE, or "ERROR_NOTSQUARE"
#Will try to do this by going through a loop of four steps-
#Step 1 - convert main diagonal to zeroes
#Step 2 - in first row, find first non-zero entry column number, then go to row corresponding to column number and find first non-zero entry in that row.
# Keep repeating this process until either a row has only zeroes or the row/column number has already been visited.
#Step three - if the row has all zeroes -> check matrix size. If matrix size is 1, return true. If not, return false.
# if the row has been repeated -> for all row/column numbers visited between the repeated rows, merge them together, adding row entries and column entries, then return to step one.
# (ex. 1->4->5->6->4 would merge rows 4,5 & 6) 
is_irreducible_matrix = function(M){
  
  if(ncol(M)!=nrow(M))
  {
    return("ERROR_NOTSQUARE")
    #Check to see it's not square
  }
  else{
    #begins by formatting matrix for search algorithm
    n = ncol(M)
   
    for(i in 1:n){
      for(j in 1:n){
        M[i,i] = 0
        M[i,j]= abs(M[i,j])
        
      }
    }
    rowloops = c(1)
    
    return(M)
  }
  
}


improved_eigen = function(M){
  EG= eigen(M)
  values= unique(EG$values)
  
  vectors= EG$vectors/colSums(EG$vectors)[col(EG$vectors)]
  
  alg_mult=  table(EG$values)
  # this gives a table with the first row being the eigenvalues, and the next being the corresponding algebraic multiplicity
  
  geo_mult = table(dim(M-diag(ncol(M))*EG$values))
  
  properties = list("values"= values,"vectors"=vectors,"alg_mult"=alg_mult,"geo_mult"= geo_mult)
  
  return(properties)
  
}

is_defective = function(M){
  if(ncol(M)==nrow(M)){
    #checks to see if the matrix is square. Otherwise, returns FALSE
    if(ncol(eigen(M)$vectors)==ncol(M)){
      #checks to see if there are n eigenvectors for an n by n matrix, otherwise, returns true
      if((det(eigen(M)$vectors)>-0.000000001)&(det(eigen(M)$vectors)<0.000000001)){
        #if there are n eigenvectors, checks whether they form a singular matrix, if so, returns true. (i.e checks whether the eigenvectors are linearly independent)
        return(TRUE)
      }
      else{
        return(FALSE)
      }
    }
    else{
      return(TRUE)
    }
  }
  else{
    return(FALSE)
  }
}



is_diagonalizable = function(M){
  #according to Horn and Johnson (pg.77) all non-defective, square matrices are diagonalizable
  if((is_defective(M)==FALSE)&(ncol(M)==nrow(M))){
    return(TRUE)
  }
  else{
    return(FALSE)
  }
}

diagonalize = function(M){
  if(is_diagonalizable(M)==TRUE){
    S= improved_eigen(M)$vectors
    D= S%*%M%*%(solve(S))
    answer = list('S'=S,"D"=D)
    return(answer)
  }
  else{
    return(FALSE)
  }
}
diagonalize(C)

is_block_diagonal = function(M){
  if(ncol(M)== nrow(M)){
    if((M[1,ncol(M)]==M[nrow(M),1])&(M[1,ncol(M)]==0)){
      #checks to see if the matrix is square and has a zero in the top right and bottom left corner, otherwise, returns false
      x=1
      
      while(x<ncol(M)){
        i= nrow(M)
        while(i>1){
          if(M[i,x]!=0){
            break
          }
          i= i-1
        }
        
        j= ncol(M)
        while(j>1){
          if(M[x,j]!=0){
            break
          }
          j= j-1
        }
        #i and j are the fist row number and column number going from bottom to top and right to left that are nonzero entries in the matrix
        
        y = max(i,j)+1
        
        if(y>=ncol(M)){
          return(FALSE)
        }
        # if y>=ncol(M), this implies that the nth bottom or right entry is not zero, despite the fact that no block has been constructed yet.
        
        top_submatrix = M[1:(y-1),y:ncol(M)]
        bottom_submatrix = M[y:nrow(M),1:(y-1)]
        # bottom and top submatrix give the top right and bottom left matrices to the centre point (y,y)
        # both bottom and top submatrices should have all zero entries if the matrix is block diagonal
        if((sum(top_submatrix)!=sum(bottom_submatrix))|(sum(top_submatrix)!=0)){
          x=x+1
        }
        else{
          return(TRUE)
        }
      }
      #end of the first while loop. it continues to loop through each row, going up from the bottom and right to left 
    }
    else{
      return(FALSE)
    }
  }
  else{
    return(FALSE)
  }
}


is_block_upper_triangular = function(M){
  if(ncol(M)== nrow(M)){
    x=1
    
    while(x<ncol(M)){
      i= nrow(M)
      while(i>1){
        if(M[i,x]!=0){
          break
        }
        i= i-1
      }
      #i is the fist row number going from bottom to top that is a nonzero entry in the matrix
      
      y = i+1
      
      if(y>=ncol(M)){
        return(FALSE)
      }
      # if y>=ncol(M), this implies that the nth bottom or right entry is not zero, despite the fact that no block has been constructed yet.
      
      bottom_submatrix = M[y:nrow(M),1:(y-1)]
      # bottom_submatrix gives the  bottom left matrix to the centre point (y,y)
      # should have all zero entries if the matrix is block upper triangular
      if(sum(bottom_submatrix)!=0){
        x=x+1
      }
      else{
        return(TRUE)
      }
    }
    #end of the first while loop. it continues to loop through each row, going up from the bottom
  }
  else{
    return(FALSE)
  }
}



is_block_lower_triangular = function(M){
  if(ncol(M)== nrow(M)){
    x=1
    
    while(x<ncol(M)){
      i= nrow(M)
      while(i>1){
        if(M[x,i]!=0){
          break
        }
        i= i-1
      }
      #i is the fist column number going from right to left that is a nonzero entry in the matrix
      
      y = i+1
      
      if(y>=ncol(M)){
        return(FALSE)
      }
      # if y>=ncol(M), this implies that the nth right entry is not zero, despite the fact that no block has been constructed yet.
      
      top_submatrix = M[1:(y-1),y:ncol(M)]
      # top_submatrix gives the  top right matrix to the centre point (y,y)
      # should have all zero entries if the matrix is block lower triangular
      if(sum(top_submatrix)!=0){
        x=x+1
      }
      else{
        return(TRUE)
      }
    }
    #end of the first while loop. it continues to loop through each row, going up from the right to left
  }
  else{
    return(FALSE)
  }
}

is_block_triangular = function(M){
  if((is_block_upper_triangular(M)==TRUE)|(is_block_lower_triangular(M)==TRUE)){
    return(TRUE)
  }
  else{
    return(FALSE)
  }
  
}



