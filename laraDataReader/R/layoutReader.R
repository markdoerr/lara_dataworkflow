#_____________________________________________________________________________
#
# PROJECT: LARA
# CLASS: 
# FILENAME: layoutReader.R
#
# CATEGORY:
#
# AUTHOR: mark doerr
# EMAIL: mark@ismeralda.org
#
# VERSION: 0.0.2
#
# CREATION_DATE: 2015/04/27
# LASTMODIFICATION_DATE: 2015/06/01
#
# BRIEF_DESCRIPTION: Library for reading microtiter plate layouts 
# DETAILED_DESCRIPTION: 
# HISTORY: 
#
# ____________________________________________________________________________
#
#   Copyright:
#   This program is free software; you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation; either version 2 of the License, or
#   (at your option) any later version.
#
#   This file is provided "AS IS" with NO WARRANTY OF ANY KIND,
#   INCLUDING THE WARRANTIES OF DESIGN, MERCHANTABILITY AND FITNESS FOR
#   A PARTICULAR PURPOSE.
#
#   For further Information see COPYING file that comes with this distribution.
#_______________________________________________________________________________

printDebug <- function(...) cat(sprintf(...), '\n' ,sep='', file=stdout())

#' loadPlateLayout
#'
#' @title Reading Plate Layout Information from Plate Layout csv File for Microtiter Plates.
#' @description Reading plate layout information from plate layout csv file for any microtiter well plates.
#'              Plate geometry is defined in layout file. Default geometry 8x12.
#'              The function is looking for a file "0*"barcode"_plate_layout.*csv" in the given directory.
#'              Default substrate, cofactor and enzyme concentrations overwrite individual concentrations.
#' @param barcode="0000" - this is used for finding the layout file BC_plate_layout_*.csv
#' @param dir="./" (string) - directory of layout file with leading / 
#' @param as_list=TRUE (boolean) this returns not only the plate layout, but also the barcodes, description and plate geometry (= layout meta information)
#' @return data.frame containing Wells, Type and Description, Substrates, Cofactors and Enzymes as columns
#' @keywords plate readers layout microtiter plate
#' @export    

loadPlateLayout <- function(as_list=TRUE, barcode="0000",  dir="./")
{
  padding <- 4 # overall number of digits, leading digits are filled with zeros
  
  file_pattern = paste("0*", barcode, "_plate_layout.*csv", sep="")
  printDebug("loadPlateLayout pattern: %s ", file_pattern)
  layout_filename <- list.files(path = dir, pattern = file_pattern  , all.files = FALSE,
                                full.names = FALSE, recursive = FALSE,
                                ignore.case = FALSE, include.dirs = FALSE)[1]
  
  #reading file only once
  if ( is.na(layout_filename) ) {print("ERROR (loadPlateLayout): no layout file found !"); return(FALSE) } else {
    data_file <- readLines(layout_filename,encoding= "UTF-8") }
  
  # as_list returns a list of layout meta information, like barcode, matrix size, and layout
  if (isTRUE(as_list)) {
    # retrieving barcode or barcode interval
    barcodes_position <- grep('Barcodes:',data_file, value=FALSE)
    
    layout_barcodes <- grep('Barcodes:',data_file, value=TRUE)
    
    if (length(layout_barcodes) > 0L){
      bc_str <- strsplit(layout_barcodes, split=":")[[1]][2]
      bc_str <- sub("\".*","", bc_str)
      bc_str_vec <- strsplit(bc_str, split=";")[[1]]
      
      # it is possible to denote an barcode interval, e.g. 4510-4555 
      check_interval <- function(bc_str){ 
        bc_interval <- as.numeric(strsplit(bc_str, split="-")[[1]])
        #printDebug("bardode interval %s",bc_interval)
        if ( length(bc_interval) == 1) return(sprintf(paste("%0",padding,"d", sep=""), bc_interval) ) else {
          # !!! works only with integer barcodes !!!
          bc_padded <- sapply(bc_interval[1]:bc_interval[2], 
                              function(x) sprintf(paste("%0",padding,"d", sep=""), x) )
          return(bc_padded)
        }
      }
      all_barcodes_vec <- unlist(lapply(bc_str_vec, check_interval))
    }
    # retrieving description of layout
    layout_description <- grep('Description:',data_file, value=TRUE)
    descr_str <- strsplit(layout_description, split=":")[[1]][2]
    description <- sub("\".*","",descr_str) # strsplit(descr_str, split="\"")[[1]][1]
  }
  
  # plate size info, if present, else using default 8x12 layout
  plate_size_info <- grep('colums=', data_file, value=TRUE)
  
  if( length(plate_size_info) > 0 )  {
    # \\2 means 2nd match of parthesized pattern
    psi <- sub("([\"# *rows=]+)(\\d+)(\\s*;\\s*)(colums=\\s*)(\\d+)(\\s*)(\",*)", "\\2;\\5", plate_size_info)
    nrows <- as.numeric(strsplit(psi, ";")[[1]][1])
    ncols <- as.numeric(strsplit(psi, ";")[[1]][2])
    if(is.na(ncols)) {
      # default values for 96 well plate
      nrows <- 8
      ncols <- 12
    }
  }

  descr_tab_position <- grep('HEAD_DESCR',data_file, value=FALSE) 
  raw_well_info_df <- read.csv(textConnection(data_file), skip=descr_tab_position-1, nrows=nrows, 
                         header=TRUE, row.names=1, 
                         stringsAsFactors=TRUE, 
                         sep=",", encoding="UTF-8", comment.ch="#")
  wells_df <- NULL
  well_info <- function(raw_well_info)
  {
    wi <- strsplit(raw_well_info, ":")[[1]] 
    temp_df <- data.frame('Type'=as.factor(wi[1]), 'Description'=wi[2] )
    wells_df <<- rbind(wells_df, temp_df )
  }
  invisible( apply(raw_well_info_df, c(2,1), well_info ))
  
  wells_df$Well=genWellNumbers(nrows=nrows, ncols=ncols)
  
  # processing subtrate 1 info 
  substr1_pos <- grep('Substrates 1:', data_file, value=FALSE)
  
  if( length(substr1_pos) > 0 )  {
    substr1_raw <- grep('Substrates 1:', data_file, value=TRUE)
    #substr1_info <- sub("(\"#\\s*Substrates\\s*1:\\s*)(\\s*[[:alpha:]\\(\\)\\+-]+)(\\s*=\\s*)(\\d*.\\d*)(\\s*)(\",*)", "\\2;\\4", substr1_raw)
    substr1_info <- sub("(\"#\\s*Substrates\\s*1:\\s*)([^=]+)(\\s*=\\s*)(\\d*.\\d*)(\\s*)(\",*)", "\\2;\\4", substr1_raw)
    substr1_name <- strsplit(substr1_info, ";")[[1]][1]
    substr1_conc <- as.numeric(strsplit(substr1_info, ";")[[1]][2])
     
    if ( (! is.na(substr1_conc))  & (substr1_conc > 0)) {
      wells_df$Substrate1 = as.factor(substr1_name)
      wells_df$Substr1Conc = substr1_conc
    } else {
      raw_substr1_info_df <- read.csv(textConnection(data_file), skip=substr1_pos, nrows=nrows, 
                                 header=TRUE, row.names=1, 
                                 stringsAsFactors=TRUE, 
                                 sep=",", encoding="UTF-8", comment.ch="#")
      temp_df <- NULL
      substr1info <- function(raw_substr1_info)
      {
        s1info <- strsplit(raw_substr1_info, ":")[[1]] 
        temp_df <<- rbind(temp_df, data.frame('Substrate1'=as.factor(s1info[1]), 'Substr1Conc'=s1info[2] ))
      }
      invisible( apply(raw_substr1_info_df, c(2,1), substr1info ))
    
      wells_df <- cbind(wells_df,temp_df)
    }
  }
  
  # processing subtrate 2 info 
  substr2_pos <- grep('Substrates 2:', data_file, value=FALSE)
  
  if( length(substr2_pos) > 0 )  {
    substr2_raw <- grep('Substrates 2:', data_file, value=TRUE)
    substr2_info <- sub("(\"#\\s*Substrates\\s*2:\\s*)([^=]+)(\\s*=\\s*)(\\d*.\\d*)(\\s*)(\",*)", "\\2;\\4", 
                        substr2_raw)
    substr2_name <- strsplit(substr2_info, ";")[[1]][1]
    substr2_conc <- as.numeric(strsplit(substr2_info, ";")[[1]][2])
    
    if ( (! is.na(substr2_conc))  & (substr2_conc > 0)) {
      wells_df$Substrate2 = as.factor(substr2_name)
      wells_df$Substr2Conc = substr2_conc
    } else {
      raw_substr2_info_df <- read.csv(textConnection(data_file), skip=substr2_pos, nrows=nrows, 
                                      header=TRUE, row.names=1, 
                                      stringsAsFactors=TRUE, 
                                      sep=",", encoding="UTF-8", comment.ch="#")
      temp_df <- NULL
      substr2info <- function(raw_substr2_info)
      {
        s2info <- strsplit(raw_substr2_info, ":")[[1]] 
        temp_df <<- rbind(temp_df, data.frame('Substrate2'=as.factor(s2info[1]), 'Substr2Conc'=s2info[2] ))
      }
      invisible( apply(raw_substr2_info_df, c(2,1), substr2info ))
      
      wells_df <- cbind(wells_df,temp_df)
    }
  }
  
  # processing Cofactors info 
  cofact_pos <- grep('Cofactors:', data_file, value=FALSE)

  if( length(cofact_pos) > 0 )  {
    cofact_raw <- grep('Cofactors:', data_file, value=TRUE)
    cofact_info <- sub("(\"#\\s*Cofactors\\s*:\\s*)([^=]+)(\\s*=\\s*)(\\d*.\\d*)(\\s*)(\",*)", "\\2;\\4", 
                        cofact_raw)
    cofact_name <- strsplit(cofact_info, ";")[[1]][1]
    cofact_conc <- as.numeric(strsplit(cofact_info, ";")[[1]][2])
    
    if ( (! is.na(cofact_conc))  & (cofact_conc > 0)) {
      wells_df$Cofactor = as.factor(cofact_name)
      wells_df$CofactConc = cofact_conc
    } else {
      raw_cofact_info_df <- read.csv(textConnection(data_file), skip=cofact_pos, nrows=nrows, 
                                      header=TRUE, row.names=1, 
                                      stringsAsFactors=TRUE, 
                                      sep=",", encoding="UTF-8", comment.ch="#")
      temp_df <- NULL
      cofactInfo <- function(raw_cofact_info)
      {
        cf_info <- strsplit(raw_cofact_info, ":")[[1]] 
        temp_df <<- rbind(temp_df, data.frame('Cofactor'=as.factor(cf_info[1]), 'CofactConc'=cf_info[2] ))
      }
      invisible( apply(raw_cofact_info_df, c(2,1), cofactInfo ))
      
      wells_df <- cbind(wells_df,temp_df)
    }
  }
  
  # processing Enzyme Concentration info 
  enzyme_pos <- grep('Enzyme Concentrations:', data_file, value=FALSE)
  
  if( length(enzyme_pos) > 0 )  {
    enzyme_raw <- grep('Enzyme Concentrations:', data_file, value=TRUE)
    enzyme_info <- sub("(\"#\\s*Enzyme\\sConcentrations\\s*:\\s*)([^=]+)(\\s*=\\s*)(\\d*.\\d*)(\\s*)(\",*)", "\\2;\\4", 
                       enzyme_raw)
    enzyme_name <- strsplit(enzyme_info, ";")[[1]][1]
    enzyme_conc <- as.numeric(strsplit(enzyme_info, ";")[[1]][2])
    
    if ( (! is.na(enzyme_conc))  & (enzyme_conc > 0)) {
      wells_df$Enzyme = as.factor(enzyme_name)
      wells_df$EnzymeConc = enzyme_conc
    } else {
      raw_enzyme_info_df <- read.csv(textConnection(data_file), skip=enzyme_pos, nrows=nrows, 
                                     header=TRUE, row.names=1, 
                                     stringsAsFactors=TRUE, 
                                     sep=",", encoding="UTF-8", comment.ch="#")
      temp_df <- NULL
      enzymeInfo <- function(raw_enzyme_info)
      {
        enz_info <- strsplit(raw_enzyme_info, ":")[[1]] 
        temp_df <<- rbind(temp_df, data.frame('Enzyme'=as.factor(enz_info[1]), 'EnzymeConc'=enz_info[2] ))
      }
      invisible( apply(raw_enzyme_info_df, c(2,1), enzymeInfo ))
      
      wells_df <- cbind(wells_df,temp_df)
    }
  }
  
  if (isTRUE(as_list)) {
    return(list('Barcodes'=all_barcodes_vec, 'LayoutDescription'=description, 
                 'Rows'=nrows, 'Columns'=ncols, 'Layout'=wells_df))
  }   else return(wells_df)
}

#' addPlateLayout
#'
#' @title Adding Plate Layout Information to reader data frame 
#' @description Convenient function to load plate layout from a layout file or database 
#'              and merging it with the reader data frame.
#' @param reader_df=NULL (data.frame) - plate reader data frame from import reader file
#' @param barcode="0000" (string) 
#' @param set_Value_NA=TRUE (boolean) - set values to NA for empty plates of type 0
#' @param set_Slope_NA=FALSE (boolean) - set slope and intercept values to NA for empty plates of type 0
#' @keywords plate readers, plate layout
#' @export 
#'   
#' @note todo : merging bug with multiple measurements per data frame (e.g. groth data)

addPlateLayout <- function(reader_df=NULL, barcode="0000", 
                           use_db_layout=FALSE, 
                           set_Value_NA=FALSE, set_Slope_NA=FALSE)
{  
  padding <- 4  # used for adding leading zeros to barcode
  # auto choose barcode from reader_df
  if (barcode == "0000") barcode <- sprintf(paste("%0", padding,"d", sep=""), as.numeric(levels(reader_df$Barcode)[1]))

  if(isTRUE(use_db_layout)) { 
    require("laraDB") 
    printDebug("addPlate layout: getting layout from DB")
    wells_df <- getPlateLayoutDB(barcode=barcode) 
  }
  else { layout_list <- loadPlateLayout(barcode=barcode)
         wells_df <- layout_list$Layout
  }
  
  # remove old Type and Description
  reader_df$Type = NULL
  reader_df$Description = NULL
  
  # merging new layout info into original data frame
  reader_df <- merge(reader_df, wells_df )
 
  #print(head(reader_df))
  # set values to NA for empty plates
  if(isTRUE(set_Value_NA)) reader_df[reader_df$Type == '0',]$Value = NA
  if(set_Slope_NA) {
    reader_df[reader_df$Type == '0',]$Slope = NA
    reader_df[reader_df$Type == '0',]$Intercept = NA
  } 
  return(reader_df)
}

