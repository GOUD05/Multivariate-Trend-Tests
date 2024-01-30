# =============================================================================
#                     FONCTION : SET UP du DOCUMENT EXCEL
# =============================================================================

excel_setup = function(Dir){
  
  setwd(Dir)
  
  # loader des packages necessaires:
  RequiredPackages <- c("copula","Kendall","VGAM","zoo","gtools","openxlsx",
                        "resample","openxlsx","foreach","doParallel",
                        "parallel","MASS")
  
  for (i in RequiredPackages) { 
    if (!require(i, character.only = TRUE)) require(i)
  }
  
  row_names = t(c("","MDT","MOT","CIT","CET"))
  col_name = paste("% p.val < alpha",sep="")
  
  compare = createWorkbook()
  addWorksheet(compare,sheetName="test_power")
  
  writeData(compare,sheet="test_power",x=row_names,colNames = FALSE,rowNames = FALSE)
  writeData(compare,sheet="test_power",x=col_name,colNames = FALSE,rowNames = FALSE,startRow = 2)
  
  top_row_style = createStyle(border="bottom",borderStyle="thick",textDecoration = "bold")
  addStyle(compare,"test_power",style=top_row_style,rows=1,cols=1:4)
  addStyle(compare,"test_power",style=createStyle(textDecoration = "bold",),rows=2:50,cols=1)
  
  saveWorkbook(compare,file="COMPARE.xlsx",overwrite = TRUE)

}