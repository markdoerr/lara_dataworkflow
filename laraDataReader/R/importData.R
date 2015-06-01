#_____________________________________________________________________________
#
# PROJECT: LARA
# CLASS: 
# FILENAME: importData.R
#
# CATEGORY:
#
# AUTHOR: mark doerr
# EMAIL: mark@ismeralda.org
#
# VERSION: 0.0.2
#
# CREATION_DATE: 2015/06/01
# LASTMODIFICATION_DATE: 2015/06/01
#
# BRIEF_DESCRIPTION: Library for reading
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


#' LA_ImportData
#'
#' @title Importing Data from Microtiter Plate Readers - generic method 
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

LA_ImportData <- function ( filename="default_file_name", data_path="./",
                           method='SPabs',  
                           barcode='0000',  wavelength=0, numSamples = 96,
                           layout=FALSE, useDBlayout=FALSE, layoutDir="./", 
                           PLC=FALSE, wellRadius=0.34, liquidVolume=0.2, ...)
{  
  current_data_path <- data_path
  current_filename <- filename
  current_method <- method
  current_barcode <- barcode
  current_wavelength <- wavelength
  current_sample_num <- numSamples
  load_layout <- layout 
  current_layout_dir <- layoutDir
  load_layout_from_db <- useDBlayout
  
  pathlength_correction <- PLC
  current_well_radius <- wellRadius 
  current_liquid_volume <- liquidVolume 
  
  if(current_filename == "" ) {print("ERROR (LA_ReaderData): no filename specified !"); return()}
  help_info <- "please type: ?importReaderData"
  
  filename_list <- list.files( pattern = current_filename, all.files = FALSE,
                           full.names = FALSE, recursive = FALSE,
                           ignore.case = FALSE, include.dirs = FALSE)
  
  if (length(filename_list) == 0L) {print("ERROR (LA_ReaderData): No files found or wrong filname pattern specified"); return(FALSE)}

  if (isTRUE(load_layout)) reader_data <- addPlateLayout(reader_df=reader_data, barcode=barcode, 
                                                    use_db_layout=use_db_layout )
  
  UseMethod("LA_ImportData")
}

LA_ImportData.default <- function (filename, ...) {
  
  printDebug("No import class specified. Please specify import reading class %s fn: %s", current_method, current_filename)
  #invisible(x)
}

#' pathlengthCorrectionFactor
#'
#' @title A simple pathlenght correction factor for microtiterplaters
#' @description A simple pathlenght correction factor is generated, based on volume and radius of a cylinder
#'              More accurate pathlength correction is achieved by measureing the absorption of a solution 
#'              with know concentration with different per well
#'              volumes and recording of a calibration curve.
#' @param liquid_volume=0.2
#' @param well_radius=0.34 (standard radius of well in 96 well MTPs)
#' @keywords plate readers, pathlength correction
#' @export 
#' @examples
#'     pl_factor <- pathlengthCorrectionFactor(liquid_volume=0.2, well_radius=0.34)
#' 

pathlengthCorrectionFactor<-function(liquid_volume=0.2, well_radius=0.34)
{
  optical_pathlength = liquid_volume / (pi * well_radius^2)
  pathlength_correction_factor = 1 / optical_pathlength  #signif((1/optical_pathlength),digits=3)
  
  return(pathlength_correction_factor)
}

#' genWellNumbers
#'
#' @title Well Number Generator
#' @description Generator for microtiter plate well coordinates
#' @param padding=2 (overall number of digits, leading digits are filled with zeros)  
#' @param nrows=8 (number of rows on plate)
#' @param ncols=12 (number of rows on plate)
#' @return as.factor(well_numbers) (factors of well numbers)
#' @keywords plate readers, microtiter plates
#' @export 
#' @examples  genWellNumbers(padding=1)
#'              >  A1, A2, ...., H12
#'            genWellNumbers(padding=3)
#'              > A001, A002, ..., H012

genWellNumbers <- function(padding=2, nrows=8, ncols=12)
{
  well_numbers <- sapply(LETTERS[1:nrows], 
                         function(x) sprintf(paste("%s%0",padding,"d", sep=""), x, c(1:ncols)) )
  return(as.factor(well_numbers)) 
}
