appendGenes <- function(mouseGenesVector, GSEMatrix)
{
  rownamesGSEMatrix <- rownames(GSEMatrix) #Get rownames from GSEMatrix (new GSE file)
  
  rowCountHumanGenes <- nrow(mouseGenesVector) #Calculate number of rows from list of full human genes
  rowCountNewGSEFile <- nrow(GSEMatrix) #Calculate number of rows of GSE matrix
  
  missing_rows <- setdiff(mouseGenesVector, rownamesGSEMatrix) #Use setdiff function to figure out rows which are different/missing from GSE matrix
  missing_rows #Display missing rows
  
  zeroExpressionMatrix <- matrix(0, nrow = length(missing_rows), ncol = ncol(GSEMatrix)) #Create a placeholder matrix with zeroes and missing_rows length as row length
  zeroExpressionMatrix[1:3, 1:3] #Check first three entries of zeroExpressionMatrix
  dim(zeroExpressionMatrix) #Check dimensions of matrix
  
  rownames(zeroExpressionMatrix) <- missing_rows #Assign row names
  colnames(zeroExpressionMatrix) <- colnames(GSEMatrix) #Assign column names
  zeroExpressionMatrix[1:3, 1:3] #Check col and row names assigned correctly
  
  fullMatrix <- rbind(GSEMatrix, zeroExpressionMatrix) #Bind GSEMatrix and zeroExpressionMatrix together
  dim(fullMatrix) #Check dimensions of fullMatrix
  
  #Reorder matrix
  fullMatrix <- fullMatrix[mouseGenesVector, ] #Reorder fullMatrix to preserve gene order
  fullMatrix[1:10, 1:3] #Check reordering
  return(fullMatrix) #Return fullMatrix
}

checkRawCounts <- function(GSEMatrix, max_log_value = 50)
{
  if (is.integer(GSEMatrix))
  {
    return("raw counts")
  }
  else if (is.double(GSEMatrix))
  {
    if(max(GSEMatrix) > max_log_value)
    {
      return("normalized")
    } 
    else if (min(GSEMatrix) < 0) 
    {
      stop("negative values detected, likely scaled data")
    } 
    else 
    {
      return("log-normalized")
    }
  }
  else 
  {
    stop("unknown matrix format: ", typeof(GSEMatrix))
  }
}