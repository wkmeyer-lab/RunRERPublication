trimGmtFilenameToSetname = function(string){
  gmtSetName = sub("\\.gmt$", "", sub(".*/", "", string))
  gmtSetName
}
