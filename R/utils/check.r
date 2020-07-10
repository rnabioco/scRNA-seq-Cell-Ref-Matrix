appendGenes <- function(geneVector, GSEMatrix)
{
  rownamesGSEMatrix <- rownames(GSEMatrix) #Get rownames from GSEMatrix (new GSE file)

  rowCountHumanGenes <- nrow(geneVector) #Calculate number of rows from list of full human genes
  rowCountNewGSEFile <- nrow(GSEMatrix) #Calculate number of rows of GSE matrix

  missing_rows <- setdiff(geneVector, rownamesGSEMatrix) #Use setdiff function to figure out rows which are different/missing from GSE matrix

  zeroExpressionMatrix <- matrix(0, nrow = length(missing_rows), ncol = ncol(GSEMatrix)) #Create a placeholder matrix with zeroes and missing_rows length as row length

  rownames(zeroExpressionMatrix) <- missing_rows #Assign row names
  colnames(zeroExpressionMatrix) <- colnames(GSEMatrix) #Assign column names

  fullMatrix <- rbind(GSEMatrix, zeroExpressionMatrix) #Bind GSEMatrix and zeroExpressionMatrix together

  #Reorder matrix
  fullMatrix <- fullMatrix[geneVector, ] #Reorder fullMatrix to preserve gene order
  return(fullMatrix) #Return fullMatrix
}

checkRawCounts <- function(GSEMatrix, max_log_value = 50)
{
  if(!is.matrix(GSEMatrix))
  {
    GSEMatrix <- as.matrix(GSEMatrix)
  }
  if (is.integer(GSEMatrix))
  {
    return("raw counts")
  }
  else if (is.double(GSEMatrix))
  {
    if (all(GSEMatrix == floor(GSEMatrix)))
    {
      return("raw counts")
    }
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