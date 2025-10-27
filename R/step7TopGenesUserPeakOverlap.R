## Internal functions used by step7TopGenesUserPeakOverlap

## Internal function which determines whether each methylation site overlaps
## with each peak in the given dataset
.peakDatasetOverlapFunction <- function(
  peakDataGRanges,
  methSiteGRanges
) {
    ## If the peak dataset has no names, generate them from the coordinates
    if (is.null(names(peakDataGRanges))) {
        names(peakDataGRanges) <- paste(
            GenomicRanges::seqnames(peakDataGRanges),
            GenomicRanges::start(peakDataGRanges),
            GenomicRanges::end(peakDataGRanges),
            sep = "_"
        )
    }

    ## Make sure the names are unique, since they will be used as column names
    names(peakDataGRanges) <- make.unique(names(peakDataGRanges))

    ## For each peak in the peak dataset, determine whether it overlaps with
    ## each methylation site
    overlaps <- GenomicRanges::findOverlaps(peakDataGRanges, methSiteGRanges)
    overlappedSites <- S4Vectors::subjectHits(overlaps)
    overlappedPeaks <- S4Vectors::queryHits(overlaps)

    ## Create a Boolean matrix with methylation sites in the rows and peaks in
    ## the columns, whose values are all initially FALSE
    overlapMatrix <- matrix(
        nrow = length(methSiteGRanges),
        ncol = length(peakDataGRanges),
        data = FALSE,
        dimnames = list(names(methSiteGRanges), names(peakDataGRanges))
    )

    ## Set the values at the intersection of the overlapped sites and peaks to
    ## TRUE
    overlapMatrix[cbind(overlappedSites, overlappedPeaks)] <- TRUE

    ## Convert the matrix to a data frame and return it
    return(as.data.frame(overlapMatrix))
}

## For the top genes/TFs, get the RE DNA methylation sites linked to them,
## overlap them with each of the datasets, then create data frames noting which
## peaks the RE DNA methylation sites overlapped with
.userPeakOverlapInternalDatasetOutputFunction <- function(
  TENETMultiAssayExperiment,
  hyperHypo,
  geneOrTF,
  geneIDdf,
  methSiteIDdf,
  peakData,
  topGeneNumber,
  distanceFromREDNAMethylationSites,
  coreCount
) {
    ## Generate the quadrant result name to grab data for
    quadrantResultsName <- paste0(hyperHypo, "methGplusResults")

    ## Ensure the quadrant's results are present in step 5
    .ensureStepPresent(
        TENETMultiAssayExperiment,
        stepName = "step5OptimizeLinks",
        substepName = quadrantResultsName
    )

    ## Get the IDs of the top genes/TFs. If there are fewer genes/TFs than the
    ## topGeneNumber specified by the user, get all the genes/TFs available.
    topQuadrantGeneOrTFIDs <- .getQuadrantTopGenesOrTFs(
        TENETMultiAssayExperiment = TENETMultiAssayExperiment,
        geneOrTF = geneOrTF,
        hyperHypo = hyperHypo,
        topGeneNumber = topGeneNumber
    )$geneID
    if (.isSingleNA(topQuadrantGeneOrTFIDs)) {
        return(NA)
    }

    ## Convert the gene IDs to gene names
    topQuadrantGeneName <- geneIDdf[
        topQuadrantGeneOrTFIDs,
        "geneName"
    ]

    ## Get all unique RE DNA methylation sites linked to at least one of the top
    ## genes selected
    quadrantMethSitesLinkedToSignificantGenes <- unique(
        TENETMultiAssayExperiment@metadata$step5OptimizeLinks[[
            quadrantResultsName
        ]][
            TENETMultiAssayExperiment@metadata$step5OptimizeLinks[[
                quadrantResultsName
            ]]$geneID %in% topQuadrantGeneOrTFIDs,
            "DNAMethylationSiteID"
        ]
    )

    quadrantMethSitesLinkedToSignificantGenes <- sort(
        quadrantMethSitesLinkedToSignificantGenes
    )

    ## Initialize a data frame with the methylation site IDs, info about
    ## the RE DNA methylation sites, and the search start and end coordinates
    methSitesLinkedToGenesDF <- data.frame(
        "DNAMethylationSiteID" = quadrantMethSitesLinkedToSignificantGenes,
        "chromosome" = methSiteIDdf[
            quadrantMethSitesLinkedToSignificantGenes,
            "chromosome"
        ],
        "start" = methSiteIDdf[
            quadrantMethSitesLinkedToSignificantGenes,
            "start"
        ],
        "end" = methSiteIDdf[
            quadrantMethSitesLinkedToSignificantGenes,
            "end"
        ],
        "searchStart" = methSiteIDdf[
            quadrantMethSitesLinkedToSignificantGenes,
            "start"
        ] - distanceFromREDNAMethylationSites,
        "searchEnd" = methSiteIDdf[
            quadrantMethSitesLinkedToSignificantGenes,
            "end"
        ] + distanceFromREDNAMethylationSites,
        stringsAsFactors = FALSE
    )

    ## Add columns to that data frame indicating
    ## which of the RE DNA methylation sites is linked to each of the top genes
    for (i in seq_along(topQuadrantGeneOrTFIDs)) {
        ## Identify if the quadrantMethSitesLinkedToSignificantGenes are among
        ## RE DNA methylation sites linked to the specific gene of interest
        TFVector <- quadrantMethSitesLinkedToSignificantGenes %in%
            TENETMultiAssayExperiment@metadata$step5OptimizeLinks[[
                quadrantResultsName
            ]][
                TENETMultiAssayExperiment@metadata$step5OptimizeLinks[[
                    quadrantResultsName
                ]]$geneID %in% topQuadrantGeneOrTFIDs[i],
                "DNAMethylationSiteID"
            ]

        methSitesLinkedToGenesDF[i + 6] <- TFVector
    }

    ## Reset the colnames and rownames of the DF
    colnames(methSitesLinkedToGenesDF) <- c(
        c(
            "DNAMethylationSiteID",
            "chromosome",
            "start",
            "end",
            "searchStart",
            "searchEnd"
        ),
        paste(
            topQuadrantGeneName,
            topQuadrantGeneOrTFIDs,
            "linked",
            sep = "_"
        )
    )
    rownames(methSitesLinkedToGenesDF) <- methSitesLinkedToGenesDF$
        DNAMethylationSiteID

    ## Create a GRanges copy of this dataset. This maintains the methylation
    ## site IDs in the names() of the GRanges object.
    methSitesLinkedToGenesDFGRanges <- GenomicRanges::makeGRangesFromDataFrame(
        df = methSitesLinkedToGenesDF[
            ,
            c("chromosome", "searchStart", "searchEnd")
        ],
        keep.extra.columns = FALSE,
        starts.in.df.are.0based = FALSE
    )

    ## For each of the peak datasets, identify if each RE DNA methylation site
    ## is found in the vicinity of each peak in each dataset
    peakOverlapInfoList <- parallel::mclapply(
        peakData,
        .peakDatasetOverlapFunction,
        methSiteGRanges = methSitesLinkedToGenesDFGRanges,
        mc.cores = coreCount
    )

    ## Create a nested list containing the peak overlap information and linked
    ## DNA methylation site information
    returnList <- list(
        "peakDatasetOverlapInfo" = peakOverlapInfoList,
        "linkedDNAMethylationSiteInfo" = methSitesLinkedToGenesDF
    )

    ## Return the list
    return(returnList)
}

## Main step7TopGenesUserPeakOverlap function

#' Identify if RE DNA methylation sites linked to top genes and transcription
#' factors are located within a specific distance of specified genomic regions
#'
#' This function takes the top genes and transcription factors (TFs) by number
#' of linked RE DNA methylation sites identified by the
#' `step6DNAMethylationSitesPerGeneTabulation` function, up to the number
#' specified by the user, and identifies if the RE DNA methylation sites linked
#' to those genes/TFs from the hyper- and/or hypomethylated G+ analysis
#' quadrants are found in the vicinity of genomic regions (peaks) of interest,
#' supplied by the user in the form of .bed, .narrowPeak, .broadPeak, and/or
#' gappedPeak files, directories containing these files, data frames, and/or
#' GRanges objects.
#'
#' @param TENETMultiAssayExperiment Specify a MultiAssayExperiment object
#' containing expression and methylation SummarizedExperiment objects, such as
#' one created by the TCGADownloader function. The object's metadata must
#' contain the results from the `step5OptimizeLinks` and
#' `step6DNAMethylationSitesPerGeneTabulation` functions.
#' @param peakData Specify a data frame, matrix, or GRanges object with
#' genomic regions (peaks) of interest, organized in a BED-like manner (see
#' <https://genome.ucsc.edu/FAQ/FAQformat.html#format1>), a path to a .bed,
#' .narrowPeak, .broadPeak, and/or .gappedPeak file with peaks of interest,
#' a path to a directory containing one or more of these file types,
#' or a named list of any of these types of input. Peak names are
#' taken from the fourth column of the input if it exists, or, if the input is a
#' GRanges object, the names of the ranges. Additional columns can be included,
#' but are not used by this function. If no names are present, they are
#' generated from peak coordinates and take the form
#' `<chromosome>\_<start>\_<end>[.<optionalDuplicateNumber>]`. Input files may
#' optionally be compressed (.gz/.bz2/.xz).
#' @param geneAnnotationDataset Specify a gene annotation dataset which is
#' used to identify names for genes by their Ensembl IDs. The argument must be
#' either a GRanges object (such as one imported via `rtracklayer::import`) or a
#' path to a GFF3 or GTF file. Both GENCODE and Ensembl annotations are
#' supported. Other annotation datasets may work, but have not been tested.
#' See the "Input data" section of the vignette for information on the required
#' dataset format.
#' Specify NA to use the gene names listed in the "geneName" column of the
#' elementMetadata of the rowRanges of the "expression" SummarizedExperiment
#' object within the TENETMultiAssayExperiment object. Defaults to NA.
#' @param DNAMethylationArray Specify the name of a DNA methylation probe
#' array supported by the sesameData package (see
#' `?sesameData::sesameData_getManifestGRanges`). If an array is specified,
#' RE DNA methylation sites and their locations in that array's manifest are
#' cross-referenced with RE DNA methylation site IDs included in the rownames
#' of the methylation dataset provided in the "methylation"
#' SummarizedExperiment object within the TENETMultiAssayExperiment object, and
#' only those overlapping will be considered for analysis. If set to NA, all RE
#' DNA methylation sites with locations listed in the rowRanges of the
#' "methylation" SummarizedExperiment object are used. Defaults to NA.
#' @param hypermethGplusAnalysis Set to TRUE to create data frames with the peak
#' overlap information for the unique hypermethylated RE DNA methylation sites
#' linked to the top genes and TFs by most hypermethylated RE DNA methylation
#' sites with G+ links. Defaults to TRUE.
#' @param hypomethGplusAnalysis Set to TRUE to create data frames with the peak
#' overlap information for the unique hypomethylated RE DNA methylation sites
#' linked to the top genes and TFs by most hypomethylated RE DNA methylation
#' sites with G+ links. Defaults to TRUE.
#' @param topGeneNumber Specify the number of top genes and TFs, based on the
#' most linked RE DNA methylation sites of a given analysis type, for which to
#' generate data showing overlap with the specified peak datasets for the RE DNA
#' methylation sites linked to those genes. Defaults to 10.
#' @param distanceFromREDNAMethylationSites Specify the distance from the linked
#' RE DNA methylation sites within which an RE DNA methylation site will be
#' considered to overlap a peak. Must be a nonnegative integer. Defaults to 100.
#' @param coreCount Argument passed as the mc.cores argument to mclapply. See
#' `?parallel::mclapply` for more details. Defaults to 1.
#' @return Returns the MultiAssayExperiment object given as the
#' TENETMultiAssayExperiment argument with an additional list named
#' 'step7TopGenesUserPeakOverlap' in its metadata containing the output of this
#' function. This list contains `hypermethGplus` and/or `hypomethGplus` lists,
#' as selected by the user, which contain lists for the top overall genes and
#' top TF genes. Each of these lists contains two elements. The first,
#' `peakDatasetOverlapInfo`, is a list of data frames named
#' after the peak datasets (without file extensions). If a single R object was
#' provided as input, the list will contain a single element named 'peakData'.
#' Each data frame contains peak names in the column names and RE DNA
#' methylation site IDs in the row names. The Boolean values indicate whether
#' each RE DNA methylation site overlaps with each peak. The second,
#' `linkedDNAMethylationSiteInfo`, is a data frame containing a row for each of
#' the unique RE DNA methylation sites linked to the top genes/TFs for the
#' specified analysis types. The columns note the location of the RE DNA
#' methylation site, the specified search window for the site, and whether the
#' site is linked to each of the top genes/TFs.
#' @export
#'
#' @examplesIf interactive()
#' ## This example uses the example MultiAssayExperiment provided in the
#' ## TENET.ExperimentHub package to overlap example peaks with all unique RE
#' ## DNA methylation sites linked to the top 10 genes by number of linked
#' ## hyper- and hypomethylated RE DNA methylation sites, using a GRanges object
#' ## containing the genomic coordinates of peaks of interest. Gene names and
#' ## the locations of RE DNA methylation sites will be retrieved from the
#' ## rowRanges of the 'expression' and 'methylation' SummarizedExperiment
#' ## objects in the example MultiAssayExperiment. A window of 100 base pairs
#' ## will be used to identify if the RE DNA methylation sites lie within the
#' ## vicinity of peaks. The analysis will be performed using one CPU core.
#'
#' ## Load the example TENET MultiAssayExperiment object
#' ## from the TENET.ExperimentHub package
#' exampleTENETMultiAssayExperiment <-
#'     TENET.ExperimentHub::exampleTENETMultiAssayExperiment()
#'
#' ## Load the example peak GRanges object from the TENET.ExperimentHub package
#' exampleTENETPeakRegions <- TENET.ExperimentHub::exampleTENETPeakRegions()
#'
#' ## Use the example datasets to perform the peak overlapping
#' returnValue <- step7TopGenesUserPeakOverlap(
#'     TENETMultiAssayExperiment = exampleTENETMultiAssayExperiment,
#'     peakData = exampleTENETPeakRegions
#' )
#'
#' ## This example uses the example MultiAssayExperiment provided in the
#' ## TENET.ExperimentHub package to overlap specified peaks with all unique RE
#' ## DNA methylation sites linked to only the top 5 genes by number of linked
#' ## hypomethylated RE DNA methylation sites. The genomic coordinates of peaks
#' ## of interest will be loaded from BED-like files located in the user's R
#' ## working directory. Gene names will be retrieved from the rowRanges of the
#' ## 'expression' SummarizedExperiment object in the example
#' ## MultiAssayExperiment, and RE DNA methylation sites and their locations
#' ## will be retrieved from the HM450 array via the sesameData package. A
#' ## window of 500 base pairs will be used to identify if the RE DNA
#' ## methylation sites lie within the vicinity of peaks. The analysis will be
#' ## performed using 8 CPU cores.
#'
#' ## Load the example TENET MultiAssayExperiment object from the
#' ## TENET.ExperimentHub package
#' exampleTENETMultiAssayExperiment <-
#'     TENET.ExperimentHub::exampleTENETMultiAssayExperiment()
#'
#' ## Use the example datasets to perform the peak overlapping
#' returnValue <- step7TopGenesUserPeakOverlap(
#'     TENETMultiAssayExperiment = exampleTENETMultiAssayExperiment,
#'     peakData = ".",
#'     DNAMethylationArray = "HM450",
#'     hypermethGplusAnalysis = FALSE,
#'     topGeneNumber = 5,
#'     distanceFromREDNAMethylationSites = 500,
#'     coreCount = 8
#' )
step7TopGenesUserPeakOverlap <- function(
  TENETMultiAssayExperiment,
  peakData,
  geneAnnotationDataset = NA,
  DNAMethylationArray = NA,
  hypermethGplusAnalysis = TRUE,
  hypomethGplusAnalysis = TRUE,
  topGeneNumber = 10,
  distanceFromREDNAMethylationSites = 100,
  coreCount = 1
) {
    ## Validate the analysis types and get a vector of the ones selected
    analysisTypes <- .validateAnalysisTypes(
        hypermethGplusAnalysis, hypomethGplusAnalysis
    )

    ## Return an error message if the input MultiAssayExperiment is invalid
    .validateMultiAssayExperiment(
        TENETMultiAssayExperiment,
        needGeneNames = is.na(geneAnnotationDataset)
    )

    ## peakData can be a character string (path to a directory), a single
    ## GRanges object, matrix, or data frame, or a named list of any of these.

    ## If one object was given, convert it into a one-element list so special
    ## cases are not needed in later code
    if (!is.list(peakData)) {
        peakData <- list("peakData" = peakData)
    } else {
        ## Ensure the list has names; later code won't work if it doesn't
        if (is.null(names(peakData))) {
            .stopNoCall(
                "The list given as the peakData argument must have names."
            )
        }
    }

    ## We don't use for(in) because we need the index to get the name
    for (i in seq_along(peakData)) {
        peakDataset <- peakData[[i]]

        ## Check whether peakDataset is a single peak dataset or file or
        ## directory path
        if (!is.character(peakDataset)) {
            ## Ensure the file is a matrix or data frame
            if (!inherits(peakDataset, "GRanges")) {
                if (!is.data.frame(peakDataset)) {
                    if (!is.matrix(peakDataset)) {
                        .stopNoCall(
                            "The input \"", names(peakData)[[i]], "\" is not ",
                            "a data frame, matrix, GRanges object, BED-like ",
                            "file, or path to a directory containing BED-like ",
                            "files."
                        )
                    } else {
                        ## It's a matrix; convert it to a data frame
                        peakData <- as.data.frame(peakData)
                    }
                }

                ## Since the dataset is not a GRanges object, convert it to one

                ## Change the column names for the first four columns in
                ## the object
                colnames(peakDataset)[seq_len(3)] <- c("chr", "start", "end")

                ## Set the row names of the data frame so they are copied to
                ## the GRanges object. peakData[, 4] will be NA if there are
                ## only three columns, correctly resulting in unset row names.
                rownames(peakDataset) <- peakData[, 4]

                ## Create a GRanges object. Assume starts are 0-based
                peakDataset <- GenomicRanges::makeGRangesFromDataFrame(
                    df = peakDataset,
                    keep.extra.columns = FALSE,
                    starts.in.df.are.0based = TRUE
                )

                ## Replace the entry in the list with its GRanges equivalent
                peakData[[i]] <- peakDataset
            }
        } else {
            ## It is a character string, so it must be a file or directory path.
            ## Ensure that the supplied path exists. If it does, load any
            ## .bed, .narrowPeak, .broadPeak, and/or .gappedPeak files found
            ## there.
            peakFileList <- .listExtBedFiles(
                extPaths = peakData,
                paramName = names(peakData)[[i]],
                paramDescription = paste(
                    "peaks for factors of interest to overlap with",
                    "linked RE DNA methylation sites"
                )
            )

            ## Create a list to store the loaded peak files as GRanges objects
            peakDataGRangesList <- list()

            for (i in peakFileList) {
                ## Load the first four columns of the file as a GRanges object
                peaksGRanges <- rtracklayer::import.bed(
                    i,
                    colnames = c("chrom", "start", "end", "name")
                )

                ## Add that GRanges object to the list
                peakDataGRangesList <- c(peakDataGRangesList, peaksGRanges)
            }

            ## Set the names of the peaks GRanges list to the names of the
            ## files, with the file extensions removed
            names(peakDataGRangesList) <- make.unique(sub(
                "\\.[^.]*(\\.gz)?$",
                "",
                basename(unlist(c(peakFileList))),
                ignore.case = TRUE
            ))

            ## Add the new peak datasets to the overall list
            peakData <- c(peakData, peakDataGRangesList)
        }
    }

    ## Get gene IDs and names from the MAE, or gene annotation dataset if
    ## provided
    geneIDdf <- .getGeneIDsAndNames(
        TENETMultiAssayExperiment, geneAnnotationDataset
    )

    ## Get methylation site IDs and names from the MAE, or methylation array if
    ## provided
    methSiteIDdf <- .getMethSiteIDsAndLocations(
        TENETMultiAssayExperiment, DNAMethylationArray
    )

    ## Create an empty list to hold the results from this step 7 function
    resultsList <- list()

    ## Generate results for the selected analysis types
    for (hyperHypo in analysisTypes) {
        ## Return results for all genes then TFs for each analysis type
        for (geneOrTF in c("Gene", "TF")) {
            ## Return DFs for each analysis type and top genes/TFs
            resultsList[[
                paste0(hyperHypo, "methGplusResults")
            ]][[
                paste0("top", geneOrTF, "s")
            ]] <- .userPeakOverlapInternalDatasetOutputFunction(
                TENETMultiAssayExperiment = TENETMultiAssayExperiment,
                hyperHypo = hyperHypo,
                geneOrTF = geneOrTF,
                geneIDdf = geneIDdf,
                methSiteIDdf = methSiteIDdf,
                peakData = peakData,
                topGeneNumber = topGeneNumber,
                distanceFromREDNAMethylationSites =
                    distanceFromREDNAMethylationSites,
                coreCount = coreCount
            )
        }
    }

    ## Add the results list to the MultiAssayExperiment
    TENETMultiAssayExperiment@metadata$
        step7TopGenesUserPeakOverlap <- resultsList

    return(TENETMultiAssayExperiment)
}
