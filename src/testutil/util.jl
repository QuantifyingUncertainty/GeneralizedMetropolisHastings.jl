###Function to display the size of the memory used by an object
function objectsizetostr(a::Any)
  bytes = Base.summarysize(a)
  result = ""
  if bytes < 10*1024
    result = @sprintf("%6d bytes  ",bytes)
  elseif bytes < 1024*1024
    result = @sprintf("%5.3f kB",bytes/(1024))
  else
    result = @sprintf("%5.3f MB",bytes/1024/1024)
  end
  result
end
