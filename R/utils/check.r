checkRawCounts <- function(GSEMatrix, max_log_value = 50)
{
  if (is.integer(GSEMatrix))
  {
    return("raw counts")
  }
  else if (is.double(GSEMatrix))
  {
    if(max(x) > max_log_val)
    {
      return("normalized")
    } 
    else if (min(x) < 0) 
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
    stop("unknown matrix format: ", typeof(x))
  }
}