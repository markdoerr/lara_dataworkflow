#_____________________________________________________________________________
#
# PROJECT: LARA
# CLASS: 
# FILENAME: plateReader.R
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
# BRIEF_DESCRIPTION: Library for reading microtiter plate reader files 
# DETAILED_DESCRIPTION: currently supported devices BMG omega and Thermo Varioskan
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

#' readerDF2Array
#'
#' @title Converts a data frame to an easy plottable 3D array
#' @description This convenient function can be used to transform a plate reader dataframe into an 3D array
#'              in the format 12 x 8 x num.
#'              It is used to quickly iterate through kinetic data for each well. 
#' @param reader_data_df (data.frame) - source data frame
#' @param num=0 (integer) - number of measurements
#' @param wavelength=0
#' @param add_zero=TRUE
#' @param use_slope=FALSE
#' @return array (12 x 8 x num)
#' @keywords plate readers
#' @export 
#' @note todo:
#' 

readerDF2Array <- function(reader_data_df, num=0, wavelength=0, add_zero=TRUE, use_slope=FALSE)
{
  if(use_slope) { reader_data.value <- reader_data_df$Slope
                  if (num == 0) num <-1 }
  else {
    if( wavelength > 0) reader_data.value <- reader_data_df[reader_data_df$Wavelength == wavelength,]$Value
    else reader_data.value <- reader_data_df$Value
  }
  
  if (num == 0) {
    print("auto determine size")
    
    date_time_ref = min(reader_data_df$DateTime)
    reader_data_df$DiffTime <-  as.factor(difftime(reader_data_df$DateTime, date_time_ref , units="mins"))
    
    reader_data_df <- reader_data_df[order(reader_data_df$DiffTime, reader_data_df$Well),]
    
    #print(head(reader_data_df, n=20))
    #print(tail(reader_data_df, n=20))
    
    num = length(levels(reader_data_df$DiffTime))
    reader_data_arr <- array(reader_data.value, c(12,8,num)) 
    
    #transform
    reader_data_lst <- lapply(c(1:num),function(n) t(reader_data_arr[,,n]))
  }
  else {
    reader_data_arr <- array(reader_data.value, c(12,8,num)) 
    #transform
    reader_data_lst <- lapply(c(1:num),function(n) t(reader_data_arr[,,n]))
  }
  
  # adding zeros if required
  if (add_zero) {
    reader_data_lst <- append( reader_data_lst, list(array(c(rep(0,96)),dim=c(8,12))), after=0)
  }
  
  #  return(array(unlist(full_data_lst),c(8,12,num_meas_cycles)))
  return(array(unlist(reader_data_lst),c(8,12,length(reader_data_lst))))
  #return(reader_data_lst)
}

### ------------- individual plate reader import routines -----------------------

#' LA_ImportData.varioskan
#'
#' @title Importing Data from Microtiter Plate Readers
#' @description This function is the main dispatcher function for reading several file formats into dataframes.
#'              Currently supported devices:
#'                BMG Readers (e.g. Omega ) - device='omega'
#'                Thermo Scientific (e.g. Varioskan) - device='varioskan'
#'              Currently supported formats (method):
#'                SPabs : Single point absorption measurement (Omega, Varioskan),
#'                SPabsMatr: Single point absorption measurement in Matrix Layout format (Omega),
#'                SPfl:   Single point fluorescense measurement (Varioskan),
#'                KINabs: Kinetic absorption measurement (Varioskan) 
#' @param    filename/pattern (string)
#' @param    path (string)
#' @param    method (string) - one of the above described methods 
#' @param    device (string) - one of the above mentioned devices
#' @param    layout=TRUE (boolean) - automatically read plate layout file
#' @param    layout_dir="./" (string) - plate layout file directory with leading / 
#' @param    use_db_layout=FALSE (boolean) - reading the layout from a database
#' @param    PLC=FALSE  (boolean) - automatically calculate pathlength correction 
#'             under the assumptions of a 96 round well SBS layout (r=034 cm)
#'             and 200ul liquid volume
#' @param    well_radius=0.34 (double) well radius for pathlength calculation in cm
#' @param    liquid_volume=0.2 (double) liquid volume for pathlength calculation in uL
#' @return   data frame with reader data
#' @keywords plate readers
#' @examples
#'           kin_df <- importReaderData("_KINabs_245_265_", method="KINabs", device='varioskan', layout=FALSE, PLC=TRUE)
#' 
#' @note     todo: debugging,  error handling, filename patterns 
#' @export 

#' LA_ImportData.varioskan
#'
#' @title Thermo Varioskan: Kinetic absorption measurement with multiple wavelengths - list file format
#' @description This function requires a table/list Varioskan  output file starting 
#'              from the keyword "Photometric1" with the following columns:
#'              'Barcode', 'Well', 'Type', 'Description', 'SampleNo', 'Value', 'Time', 'Wavelength', 'Read'
#'              
#' @param filename (string) - single varioskan outputfile of kinetic data
#' @param barcode="0000" (string) - for defining a new barcode
#' @return data.frame
#' @keywords plate readers, Thermo Varioskan
#' @export 
#' 
#' @note todo - adding number of measurements, Time as factors ??

LA_ImportData.varioskan <- function(filename, coords=FALSE, ...)
{
  if(isTRUE(coords)) printDebug("reading Varioskan SP abs file with coordinates: %s", filename_list[1]) else
    printDebug("reading Varioskan KIN abs file: %s \nWarning: reading possibly multiple files", filename_list[1])
 
  table_offset = 2  # offset between the keyword 'Photometric1' and the begin of the data table
 
  switch(current_method,
         SPabs={  table_start_string <- 'Photometric1\t\t'
                  # column names for output data frame
                  if ( isTRUE(coords)) col_names <- c('Barcode', 'Well', 'Type', 'Description', 'SampleNo', 
                                                      'Value', 'Time', 'Wavelength', 'C1', 'C2', 'C3') else
                                       col_names <- c('Barcode', 'Well', 'Type', 'Description', 'SampleNo', 
                                                               'Value', 'Time', 'Wavelength', 'Read')
         },
         KINabs={ table_start_string <- 'Photometric1\t\t'
                  # column names for output data frame
                  if (isTRUE(coords)) col_names <- c('Barcode', 'Well', 'Type', 'Description', 'SampleNo', 
                                                    'Value', 'Time', 'Wavelength', 'C1', 'C2', 'C3') else
                                      col_names <- c('Barcode', 'Well', 'Type', 'Description', 'SampleNo', 
                                                                'Value', 'Time', 'Wavelength', 'Read')},
         SPfl={table_start_string <- 'Fluorometric1\t\t'
               col_names <- c('Barcode', 'Well', 'Type', 'Description', 'SampleNo', 'Value', 'Duration', 'ExWL', 'EmWL') },
         { cat(current_method, " - method not found Varioskan device\n"); return(help_info); }
  )
  
  read_all_data <- function(filename)
  {  
    # reading input file only once
    data_file <- readLines(filename,encoding= "UTF-8")
    raw_date <- gsub("\t", "", grep('Run started',data_file, value=TRUE), fixed=T)
    print(raw_date)
    # %I hours in 1-12, %p AM/PM
    date_time <- strptime(raw_date, format = "Run started%m/%d/%Y %H:%M:%S+02:00", tz = "CET")
    print(date_time)
    # detecting number of wavelengths
    num_wavelengths <<- length(grep("Wavelength \\[nm\\]", data_file))
    # detecting number of readings
    
    num_readings <<- as.numeric(strsplit(grep('Readings',data_file, value=TRUE)[1], "\t")[[1]][6])
    if (is.na(num_readings) ) num_readings <<- 1 
    
    # finding start of data
    table_start_row <- grep(table_start_string, data_file)[2] + table_offset 
    
    # printDebug(num_readings, num_wavelengths,  current_sample_num)
    raw_tab <- read.table(filename, skip = table_start_row, header=FALSE, col.names=col_names,
                          nrows=num_readings * num_wavelengths * current_sample_num)
    
    # converting Barcode colum to factos
    if( current_barcode == "0000" ) raw_tab$Barcode <- as.factor(raw_tab$Barcode)
    else raw_tab$Barcode <- as.factor(current_barcode)
    # adding num for number of measurements - for compatiblity reasons - todo: fill with real number
    raw_tab$NumReads <- as.factor(file_num) 
    raw_tab$DateTime <- date_time 
    
    full_abs_df <<- rbind(full_abs_df,raw_tab)
    file_num <<- file_num + 1
  }
  
  file_num <- 1
  full_abs_df <- NULL
  sapply(filename_list, read_all_data )
  
  output_df <- full_abs_df[order(full_abs_df[,2]),]
  
  output_df$Rown = as.factor(rep(LETTERS[1:8], each=12*num_readings*num_wavelengths))
  output_df$Coln = as.factor(rep(rep(1:12), each=num_readings*num_wavelengths))
  
  if (isTRUE(coords)) {
    # removing coord columns
    output_df$C1 <- NULL
    output_df$C2 <- NULL
    output_df$C3 <- NULL
  }
  
  if ( isTRUE(load_layout) ) {
    reader_data <- addPlateLayout(reader_df=output_df, barcode=current_barcode, 
                                  use_db_layout=load_layout_from_db )
  }
  
  # **** optical pathlength correction
  if (isTRUE(pathlength_correction)) {
    # liquid volume in ml, radius in cm
    pathlen_fac <- pathlengthCorrectionFactor( liquid_volume=liquid_volume,
                                               well_radius=well_radius)
    reader_data$Value <- reader_data$Value * pathlen_fac
  }
  
  # returning data frame ordered by well names
  return(output_df)
}

LA_ImportData.omega <- function(filename, coords=FALSE, ...)
{
  printDebug("plateReader - omegaSPabs: reading Omega SPabs - table format")
  
  read_all_data <- function(filename)
  {  
    table_offset = 2  # offset between the keyword 'Chromatic' and the begin of the data table
    #reading file only once  ?collapse=" ",
    data_file <- readLines(filename, encoding= "UTF-8")
    # reading time
    raw_date <- grep('Date:',data_file, value=TRUE)
    # %I hours in 1-12, %p AM/PM
    date_time <- strptime(raw_date, format = "Date: %Y/%m/%d Time: %I:%M:%S %p", tz = "CET")
    # reading barcode of plate
    barcode <- strsplit(grep('ID1:.*ID2:', data_file, value=TRUE), " ")[[1]][2]
    # finding wavelengths
    wavelength_lines <- grep('\\d{3}nm',data_file, value=TRUE)
    # removing "nm" and gain values  from wavelengths
    wavelengths <- sub("nm.*", "", sub(".*: ", "",  wavelength_lines))
    # finding all data tables
    abs_table_start_rows <- grep('Chromatic', data_file)
    
    readMultipleWavelengthData <- function(wavelength)
    {
      # reading absorption data
      table_start_row <- abs_table_start_rows[i] + table_offset
      col_names <- c('Well', 'Value')
      raw_tab <- read.table(textConnection(data_file), skip = table_start_row, header=TRUE,
                            col.names=col_names, nrows=current_sample_num)  
      
      raw_tab$Barcode <- as.factor(barcode)
      raw_tab$Num <- as.factor(file_num)
      raw_tab$DateTime <- date_time
      raw_tab$Type <- as.factor("sample")
      raw_tab$Wavelength <- as.factor(wavelength) 
      
      i <<- i+1 
      full_abs_df <<- rbind(full_abs_df,raw_tab)
    }
    # iterating reading through all wavelength tables
    i <- 1  # "static" variable for iterating through table starts
    abs_data <- sapply(wavelengths, readMultipleWavelengthData )
    file_num <<- file_num + 1
  } 
  
  # iterating through all files in the file list
  file_num <- 1
  full_abs_df <- NULL
  sapply(filename_list, read_all_data )
  
  output_df <- full_abs_df[order(full_abs_df[,2]),]
  
  # pathlen corr and layout needed !
  
  return(output_df)
}
