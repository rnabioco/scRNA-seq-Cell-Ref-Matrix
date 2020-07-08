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