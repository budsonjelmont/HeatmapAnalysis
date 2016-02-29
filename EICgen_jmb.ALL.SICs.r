#!/usr/bin/Rscript
#Similar to Andy's EICgen program, except that it plots all SICs for an experiment without looking for matches in another dataset (i.e. no xcms vs xdk comparison)
# TO USE
# Rscript path_to_EICgen.r path_to_files data_type num_rep num_col time_point peptide_counter
# Put mzXML's into "mzXML" folder nested in path_to_files
# Put autoFill xml input files into "xml" folder nested in path_to_files
# Put complete output files into "txt" folder nested in path_to_files

suppressPackageStartupMessages({
    library(XML)
    library(xcms)
    # library(gdata)
    library(data.table)
    })


parseFMXML <- function(path) {
    fieldCounter <- 0
    rowCounter <- 0
    colCounter <- 0
    datCounter <- 0
    colBool <- FALSE
    datBool <- FALSE
    fieldSizes <- list()
    parsedXML <- list()
    invisible(xmlEventParse(
        path, 
        handlers = list(
            startElement = function(name, attr) {
                if (name == "FIELD") {
                    fieldCounter <<- fieldCounter + 1
                    fieldSizes[[fieldCounter]] <<- as.numeric(attr[[2]])
                }
                if (name == "RESULTSET") {
                    rowCounter <<- 0
                    parsedXML <<- lapply(1:attr, function(x) lapply(fieldSizes, function(y) rep(NA, y)))
                }
                if (name == "ROW") {
                    rowCounter <<- rowCounter + 1
                    colCounter <<- 0
                    datCounter <<- 0
                }
                if (name == "COL") {
                    colBool <<- TRUE
                    colCounter <<- colCounter + 1
                }
                if (name == "DATA") {
                    datBool <<- TRUE
                    datCounter <<- datCounter + 1
                }
            }, text = function(text) {
                if (datBool) {
                    parsedXML[[rowCounter]][[colCounter]][datCounter] <<- text
                }
            }, endElement = function(name, uri) {
                if (name == "COL") {
                    colBool <<- FALSE
                    datCounter <<- 0
                }
                if (name == "DATA") {
                    datBool <<- FALSE
                }
            }
        )))
    return(do.call(rbind, parsedXML))
    }

parseFile1 <- function(file1Path, mzXMLPath) {
    # print(paste("Parsing ", file1Path, sep=""))
    if (!file.exists(file1Path)) waitForInput(paste(file1Path, " does not exist. Please check input folder.", sep=""))

    parsedXML <- parseFMXML(file1Path)
    parsedF1 <- cbind(parsedXML, NA)
    parsedF1[,1:2] <- as.numeric(parsedF1[,1:2])
    colnames(parsedF1) <- c("column", "fraction", "path", "xcmsRaw")
    # rownames(parsedF1) <- c(unlist(parsedF1[,"path"]))

    columns <- unlist(parsedF1[,"column"])
    fractions <- unlist(parsedF1[,"fraction"])

    minCol <- min(columns)
    maxCol <- max(columns)
    minFrac <- min(fractions)
    maxFrac <- max(fractions)

    dat <- matrix(list(), nrow=maxCol, ncol=maxFrac)

    for (i in 1:nrow(parsedF1)) {
        # print(paste(mzXMLPath, unlist(fileNameNoExt(basename(parsedF1[i,"path"][[1]]))), ".mzXML", sep=""))
        dat[parsedF1[i,"column"][[1]], parsedF1[i,"fraction"][[1]]][1] <- xcmsRaw(paste(mzXMLPath, unlist(fileNameNoExt(basename(parsedF1[i,"path"][[1]]))), ".mzXML", sep=""), profstep=0)
    }

    # apply(parsedF1, 1, function(x) dat[x[[1]], x[[2]]] <<- xcmsRaw(paste(mzXMLPath, unlist(fileNameNoExt(basename(x[[3]]))), ".mzXML", sep=""), profstep=0))

    return(dat)
}

parseFile2 <- function(file2Path, SILAC) {
    # print(paste("Parsing ", file2Path, sep=""))
    if (!file.exists(file2Path)) waitForInput(paste(file2Path, " does not exist. Please check input folder.", sep=""))

    parsedXML <- parseFMXML(file2Path)
    counterIndex <- 1
    rtIndex <- 2
    pathIndex <- ifelse(SILAC, 3, 5)
    massIndex <- 4
    if (SILAC) {
        parsedXML[, massIndex] <- lapply(1:nrow(parsedXML), function(i) as.numeric(unlist(strsplit(unlist(parsedXML[i, massIndex]), ","))))
        parsedXML[, counterIndex] <- as.numeric(parsedXML[, counterIndex])
        parsedXML[, rtIndex] <- lapply(1:nrow(parsedXML), function(j) as.numeric(unlist(parsedXML[j, rtIndex])))
    } else {
        parsedXML[, massIndex] <- as.numeric(parsedXML[, massIndex])
        parsedXML[, counterIndex] <- as.numeric(parsedXML[, counterIndex])
        parsedXML[, rtIndex] <- lapply(1:nrow(parsedXML), function(i) as.numeric(unlist(parsedXML[i, rtIndex])))
        parsedXML <- parsedXML[,c(counterIndex, rtIndex, pathIndex, massIndex)]
    }

    colnames(parsedXML) <- c("counter", "rt", "path", "mass")

    return(parsedXML)
}

parseTxtFile <- function(txtfilePath) {
    # print(paste("Parsing ", txtfilePath, sep=""))
    if (!file.exists(txtfilePath)) waitForInput(paste(txtfilePath, " does not exist. Please check input folder.", sep=""))

    alignFile <- read.table(txtfilePath, sep="\t")
    return(alignFile)
}

fileNameNoExt <- function(fileName) {
    sub("^([^.]*).*", "\\1", fileName)
}

waitForInput <- function(error="") {
    cat(paste("Error occured: ", error, "\n", sep=""))
    cat("To continue, type [y]. To cancel, type [n]. Then press [enter].\n")
    key <- scan("stdin", n=1, what="")
    if (key == "y") {
        cat("Continuing job.\n")
    } else if (key == "n") {
        cat("Canceling job.\n")
        quit()
    } else {
        cat("Invalid key input, try again\n")
        waitForInput()
    }
}

pmMz <- function (x, ppm=10) {
    c(x - x * ppm * 0.000001, x + x * ppm * 0.000001)
}

pmRt <- function (x, rtwidth=10) {
    c(x - rtwidth, x + rtwidth)
}

rtRangeToScanRange <- function (xcmsRaw, rtrange) {
    rtrange <- range(rtrange)
    scanidx <- (xcmsRaw@scantime >= rtrange[1]) & (xcmsRaw@scantime <= rtrange[2])
    c(match(TRUE, scanidx), length(scanidx) - match(TRUE, rev(scanidx)))
}

RTWINDOW <<- 30
txtdir <<- "\\TXT\\"
mzXMLdir <<- "\\mzXML\\"
xmldir <<- "\\XML\\"

args <- commandArgs(trailingOnly = TRUE)

path <- args[1]

xcms <- read.table(args[2],sep="\t",header=TRUE,fill=TRUE,comment.char="")
dataType <- args[3]
repCount <- as.numeric(args[4])
colCount <- as.numeric(args[5])
timeP <- as.numeric(args[6])

if (dataType != "LABELFREE") {
    quit("Wrong type of data")
}

setwd(path)

#xcms$id <- paste(xcms[,'protein.name.manual'],xcms[,'phosphosite.annotated'],xcms[,'charge.state.peptide'])
xcms$id <- paste(xcms[,'protein.name.manual'],xcms[,'charge.state.peptide'])

peptides <- as.numeric(xcms$counter)

colCountHalf <- ceiling(as.numeric(colCount)/2)
adjtimeP <- timeP + 1

txtPath <- paste(path, txtdir, sep="")
mzXMLPath <- paste(path, mzXMLdir, sep="")
xmlPath <- paste(path, xmldir, sep="")

for (i in 1:repCount) {
    txtFile <- paste(txtPath, "completeLabelFree", i, ".txt", sep="")
    xml1File <- paste(xmlPath, "replicate", i, "_file1_labelfree.xml", sep="")
    xml2File <- paste(xmlPath, "replicate", i, "_file2_labelfree.xml", sep="")
    if (file.exists(txtFile) && file.exists(xml1File) && file.exists(xml2File)) {
        print(paste("Replicate", i))
        parsedXML1 <- parseFile1(xml1File, mzXMLPath)
        parsedXML2 <- parseFile2(xml2File, FALSE)
        parsedtxt <- parseTxtFile(txtFile)
        for (pep in peptides) {
		    if(!is.na(pep)){
				print(pep)
				matchIndexXml <- which(as.numeric(parsedXML2[,"counter"]) == pep)
				mass <- parsedXML2[matchIndexXml,"mass"][[1]]
				matchIndexTxt <- which(as.numeric(parsedtxt[,"V1"]) == pep)
				row <- parsedtxt[matchIndexTxt,]
				col1 <- adjtimeP
				frac <- 1

				tmp <- as.numeric(unlist(strsplit(as.character(unlist(row[col1])), ",")))

				alRT <- tmp[1] * 60
				#ignore columns with missing data (alRT should = 0 if there's no data)
				if(alRT != 0){
					RTl <- tmp[2]
					detectedRT <- tmp[3]
					RTr <- tmp[4]
					intensity <- tmp[5]
					rtw <- tmp[6]
					
					if(rtw == 0){
						rtw <- 10
					}

					xr1 <- parsedXML1[col1-1, frac][[1]]
					
					plot1Name <- paste("pep ", pep, " rep ", i, " col ", timeP, " alRT ", alRT, " RTl ",RTl, " detectedRT ", detectedRT, " RTr ", RTr," mass ",mass, " int ", intensity, " 10ppm.png", sep="")
					SIC1Name <- paste("pep ", pep, " rep ", i, " col ", timeP, " alRT ", alRT, " RTl ",RTl, " detectedRT ", detectedRT, " RTr ", RTr, " mass ",mass," int ", intensity, " 10ppm.txt", sep="")

					png(plot1Name)
					tryCatch(plotEIC(xr1, mzrange=pmMz(mass, 10), rtrange=pmRt(alRT, RTWINDOW)), error=function(ex) NA)
					dev.off()
					if(RTl - rtw >= 0){
						RTrange <- c(RTl-(RTr-RTl),RTr+(RTr-RTl))
					} else {
						RTrange <- c(0,RTr+(RTr-RTl))
					}
					output <- rawEIC(xr1,mzrange=pmMz(as.numeric(mass), 10),rtrange=RTrange)
					output <- cbind(xr1@scantime[output$scan],output$intensity)
					write.table(output, file=SIC1Name, sep="\t", quote=FALSE, col.names = FALSE, row.names = FALSE)
				}
	        }
        }
    }
}