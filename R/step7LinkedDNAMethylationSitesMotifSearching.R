## Internal functions used by step7LinkedDNAMethylationSitesMotifSearching

## Internal function to validate a PWM provided by the user in the TFMotifList
## Checks if the PWM is a matrix, and if it has 4 rows
.internalPWMvalidator <- function(PWMMmatrix) {
    return(is.matrix(PWMMmatrix) && nrow(PWMMmatrix) == 4)
}

## Internal function to pull gene names and IDs from booleans
.internalGeneNameIDFromBool <- function(IDBool,
                                        nameBool,
                                        TFMotifListName,
                                        geneIDDF) {
    ## First identify if the IDBool is is TRUE
    if (IDBool) {
        ## If IDBool is TRUE, that means the TFMotifListName is the gene ID
        geneIDValue <- TFMotifListName

        ## We can get the gene name by getting the row where the TFMotifList
        ## name matches its geneName, then getting the geneName for that gene
        geneNameValue <- geneIDDF[
            which(TFMotifListName == geneIDDF$geneID),
            "geneName"
        ]
    } else if (nameBool) {
        ## If nameBool is TRUE, that means the TFMotifListName is the gene name
        geneNameValue <- TFMotifListName

        ## We can get the gene ID by getting the row where the TFMotifList
        ## name matches its geneName, then getting the geneID for that gene
        geneIDValue <- geneIDDF[
            which(TFMotifListName == geneIDDF$geneName),
            "geneID"
        ]
    } else {
        ## If both are FALSE, then the name/ID can't be found in the
        ## TFMotifListName at all
        geneIDValue <- NA
        geneNameValue <- NA
    }

    ## Return a vector with the gene name then ID
    return(c(geneIDValue, geneNameValue))
}

## Internal function to get the DNA string of the motif matches
.internalMotifGrabber <- function(start, end, DNAString) {
    return(as.character(DNAString[start:end]))
}

## Internal function to find the motif occurrences in the specified vicinity of
## each linked RE DNA methylation site. Note: All instances of + 1 and -1 are
## necessary and correct.
.findMotifSurroundingMethSite <- function(index,
                                          GRangesObject,
                                          genome,
                                          motifPWM,
                                          matchPWMMinScore) {
    ## Get site information from the GRangesObject
    methSiteChr <- as.character(GenomicRanges::seqnames(GRangesObject))[index]
    methSiteStart <- as.numeric(GenomicRanges::start(GRangesObject))[index]
    methSiteEnd <- as.numeric(GenomicRanges::end(GRangesObject))[index]
    DNAMethylationSiteID <- names(GRangesObject)[index]

    ## Get the DNA string sequence for that segment from the genome
    DNAString <- genome[[methSiteChr]][methSiteStart:methSiteEnd]

    ## Do the motif search
    search <- Biostrings::matchPWM(motifPWM, DNAString, matchPWMMinScore)

    ## Return a data frame with the information about the found motifs

    ## Get the start and end locations where the motif was found
    starts <- IRanges::start(search@ranges)
    ends <- IRanges::end(search@ranges)

    ## Get the motif matches
    motifMatches <- mapply(
        .internalMotifGrabber,
        starts,
        ends,
        MoreArgs = list("DNAString" = DNAString)
    )

    ## Create a data frame with results and return it. The start and end
    ## position calculations look strange but are correct, since matchPWM
    ## is being run only on the DNA string of the region surrounding the RE
    ## DNA methylation site, so we need to add the start of the region
    ## surrounding the RE DNA methylation site to the position of the motif
    ## in the string to get the location of the motif on the whole
    ## chromosome.
    if (length(as.character(search)) == 0) {
        ## No results are found, so try creating an empty data frame
        ## with the columns of interest, but 0 rows
        return(data.frame(matrix(nrow = 0, ncol = 5)))
    } else {
        returnDF <- data.frame(
            rep(DNAMethylationSiteID, length(starts)),
            rep(methSiteChr, length(starts)),
            (starts + methSiteStart - 1),
            (ends + methSiteStart - 1),
            motifMatches
        )
        colnames(returnDF) <- paste0("X", seq_len(ncol(returnDF)))
        return(returnDF)
    }
}

## Main step7LinkedDNAMethylationSitesMotifSearching function

#' Perform motif searching for transcription factor motifs in the vicinity of
#' DNA methylation sites or custom regions defined by the user.
#'
#' This function takes a user-specified named list of transcription factors
#' (TFs) and their binding motifs combined with search terms passed to MotifDb's
#' `query()` function to identify additional TF binding motifs. For each
#' of the TFs, this function identifies if the specified motif, in the from of
#' a position weight matrix (PWM), is found within a user-specified distance to
#' RE DNA methylation sites from the hyper- or hypomethylated G+ analysis
#' quadrants and/or specified RE DNA methylation sites, and/or within custom
#' genomic regions, as selected by the user.
#'
#' **Note**: When running this function, it is recommended to either select a
#' small number of TFs and larger number of DNA methylation sites, or a
#' small number of sites and larger number of TFs, as motif analysis can take
#' a significant amount of time to run.
#'
#' @param TENETMultiAssayExperiment Specify a MultiAssayExperiment object
#' containing expression and methylation SummarizedExperiment objects, such as
#' one created by the TCGADownloader function. This MultiAssayExperiment object
#' should also contain the results from the `step5OptimizeLinks` function
#' in its metadata if `hypermethGplusAnalysis` or `hypomethGplusAnalysis` are
#' TRUE.
#' @param hypermethGplusAnalysis Set to TRUE to do motif searching in the
#' vicinity of hypermethylated RE DNA methylation sites with at least one linked
#' TF. **Note**: if `useOnlyDNAMethylationSitesLinkedToTFs` is also TRUE, only
#' RE DNA methylation sites linked specifically to TFs specified via the
#' `TFMotifList` argument will be used. Defaults to TRUE.
#' @param hypomethGplusAnalysis Set to TRUE to do motif searching in the
#' vicinity of hypomethylated RE DNA methylation sites with at least one linked
#' TF. **Note**: if `useOnlyDNAMethylationSitesLinkedToTFs` is also TRUE, only
#' RE DNA methylation sites linked specifically to TFs specified via the
#' `TFMotifList` argument will be used. Defaults to TRUE.
#' @param DNAMethylationSites Supply a vector of DNA methylation site IDs
#' to perform motif searching in the vicinity of. These sites will be combined
#' with any RE DNA methylation sites selected by the `hypermethGplusAnalysis`
#' and `hypomethGplusAnalysis` arguments. If set to NA, no additional DNA
#' methylation sites in the TENETMultiAssayExperiment will be included in the
#' motif search. Defaults to NA.
#' @param distanceFromREDNAMethylationSites Specify a positive integer in base
#' pairs to be the distance from the DNA methylation sites identified by the
#' `hypermethGplusAnalysis`, `hypomethGplusAnalysis`, and `DNAMethylationSites`
#' arguments within which motif searching will be performed. Defaults to 100.
#' @param GRangesToSearch Specify a GRanges object which contains user-specified
#' genomic coordinates of regions to perform motif searching on. The coordinates
#' should correspond to the human hg38 genome. Any regions included in this
#' GRanges object will be combined with regions defined by the
#' `hypermethGplusAnalysis`, `hypomethGplusAnalysis`, `DNAMethylationSites`, and
#' `distanceFromREDNAMethylationSites` arguments. If set to NA, no additional
#' user specified ranges will be included in the motif search. Defaults to NA.
#' @param andStrings Specify a vector of values which will be provided to the
#' `andStrings` argument of the `query()` function in the MotifDb package. Good
#' values to provide to this vector include species and transcription factor
#' database names to refine the search. Set to NULL to include no terms in this
#' search (**Note:**: if both `andStrings` and `orStrings` are set to NULL, only
#' the PWMs specified by the `TFMotifList` argument will be used). Defaults to
#' NULL.
#' @param orStrings Specify a vector of values which will be provided to the
#' `orStrings` argument of the `query()` function in the MotifDb package. Good
#' values to provide to this vector include names of specific TFs to limit the
#' search to. The user may also specify "humanTranscriptionFactors" to use all
#' TFs identified in 'The Human Transcription Factors' by Lambert et al. 2018.
#' Set to NULL to include no terms in this search (note, if both `andStrings`
#' and `orStrings` are set to NULL, only the PWMs specified by the `TFMotifList`
#' argument will be used). TFs with valid PWMs found in combination with the
#' `andStrings` object will be included along with PWMs specified by the user in
#' a list given to the `TFMotifList` argument in the final search assessing TF
#' motifs in the vicinity of the regions specified previously. Defaults to NULL.
#' @param notStrings Specify a vector of values which will be provided to the
#' `notStrings` argument of the `query()` function in the MotifDb package. Set
#' to NULL to exclude no terms from this search. Defaults to NULL.
#' @param TFMotifList Specify a named list mapping transcription factor gene
#' names and/or IDs to their respective motif position weight matrix (PWM). The
#' PWMs should be in the form of a 4xN matrix. PWMs specified in this list are
#' combined with any TF PWMs identified by the MotifDb package using the
#' `andStrings` and `orStrings` arguments. Set to NA to not use any user
#' specified TF motif PWMs.
#' @param useOnlyDNAMethylationSitesLinkedToTFs If set to TRUE, only
#' hypomethylated or hypermethylated RE DNA methylation sites, as selected by
#' the `hypermethGplusAnalysis` and `hypomethGplusAnalysis` arguments, which are
#' found to be linked to the TFs in the given `TFMotifList` by TENET will be
#' analyzed. To use this functionality, at least one of `hypermethGplusAnalysis`
#' or `hypomethGplusAnalysis` must be set to TRUE, `DNAMethylationSites`,
#' `andStrings`, and `orStrings` must be NA, and the name of each PWM in the
#' list given to `TFMotifList` must match the gene name or Ensembl ID of a gene
#' in the TENETMultiAssayExperiment with RE DNA methylation sites linked to it
#' for the specified analysis types. Defaults to TRUE.
#' @param geneAnnotationDataset Specify a gene annotation dataset which is
#' used to identify names for genes by their Ensembl IDs. The argument must be
#' either a GRanges object (such as one imported via `rtracklayer::import`) or a
#' path to a GFF3 or GTF file. Both GENCODE and Ensembl annotations are
#' supported. Other annotation datasets may work, but have not been tested.
#' See the "Input data" section of the vignette for information on the required
#' dataset format.
#' Specify NA to use the names for genes listed in the "geneName" column of the
#' elementMetadata of the rowRanges of the "expression" SummarizedExperiment
#' object within the TENETMultiAssayExperiment object. Defaults to NA.
#' @param DNAMethylationArray Specify the name of a DNA methylation probe array
#' supported by the sesameData package (see
#' `?sesameData::sesameData_getManifestGRanges`). If an array is specified, RE
#' DNA methylation sites and their locations in that array's manifest are
#' cross-referenced with RE DNA methylation site IDs included in the rownames
#' of the methylation dataset provided in the "methylation"
#' SummarizedExperiment object within the TENETMultiAssayExperiment object, and
#' only those overlapping will be considered for analysis. If set to NA, all RE
#' DNA methylation sites with locations listed in the rowRanges of the
#' "methylation" SummarizedExperiment object are used. Defaults to NA.
#' @param matchPWMMinScore Specify the `min.score` argument passed to the
#' matchPWM function for motif searching. See `?Biostrings::matchPWM` for more
#' details. Defaults to "75%".
#' @param coreCount Argument passed as the mc.cores argument to mclapply. See
#' `?parallel::mclapply` for more details. Defaults to 1.
#' @return Returns the MultiAssayExperiment object given as the
#' TENETMultiAssayExperiment argument with an additional list of information
#' named 'step7LinkedDNAMethylationSitesMotifSearching' in its metadata with
#' the output of this function. This list includes the GRanges object
#' "DNAMethylationSitesGRanges" with the regions motif searching was
#' performed on, "TFMotifPWMList" which includes the TF PWMs used in the
#' searching, "TFMotifSeqLogoList" which includes the visual seqLogo
#' representations of the PWMs included in the "TFMotifPWMList", the
#' "DNAMethylationSitesMotifOccurrences" data frame, which notes all
#' predicted motifs with their location, the DNA methylation region they were
#' found within, and the TF PWM that was found, as well as a second
#' "totalMotifOccurrencesPerDNAMethylationSite" data frame, which notes how
#' many times each TF PWM listed in the "TFMotifPWMList" was found in each
#' region in the "DNAMethylationSitesGRanges". If
#' useOnlyDNAMethylationSitesLinkedToTFs was set to TRUE, a final data frame
#' "linkedUniqueDNAMethylationSitesTFOverlap" is included, which notes which of
#' the TFs in the "TFMotifPWMList" the identified hyper- or hypomethylated RE
#' DNA methylation sites used in the analysis were linked to (otherwise this
#' will be NA).
#' @export
#'
#' @examplesIf interactive()
#' ## Show available motifs for example TF FOXA1
#' names(MotifDb::query(MotifDb::MotifDb, "FOXA1"))
#'
#' ## The seqLogos for all input motifs will also be included in the output
#' ## of this function. Alternatively, individual motifs can be visualized
#' ## with the seqLogo function from the seqLogo package.
#' seqLogo::seqLogo(MotifDb::query(MotifDb::MotifDb, "FOXA1")[[3]])
#'
#' ## Once PWMs have been selected for use, a list containing them must be
#' ## created
#' exampleTFMotifList <- list(
#'     "FOXA1" = MotifDb::query(MotifDb::MotifDb, "FOXA1")[[3]],
#'     "ESR1" = MotifDb::query(MotifDb::MotifDb, "ESR1")[[4]]
#' )
#'
#' ## This example uses the example MultiAssayExperiment provided in the
#' ## TENET.ExperimentHub package to perform motif overlapping for all
#' ## hyper- and hypomethylated RE DNA methylation sites linked to the
#' ## FOXA1 and ESR1 genes using the motifs for each gene specified in the
#' ## exampleTFMotifList. Gene names and locations, and the locations of RE
#' ## DNA methylation sites, will be retrieved from the rowRanges of the
#' ## 'expression' and 'methylation' SummarizedExperiment objects in the
#' ## example MultiAssayExperiment. Regions within 100 bp of linked RE DNA
#' ## methylation sites will be checked for motifs, and a similarity
#' ## threshold of 75% will be used to identify motifs. The analysis will
#' ## be performed using one CPU core.
#'
#' ## Load the example TENET MultiAssayExperiment object
#' ## from the TENET.ExperimentHub package
#' exampleTENETMultiAssayExperiment <-
#'     TENET.ExperimentHub::exampleTENETMultiAssayExperiment()
#'
#' ## Use the example dataset to do the motif searching
#' returnValue <- step7LinkedDNAMethylationSitesMotifSearching(
#'     TENETMultiAssayExperiment = exampleTENETMultiAssayExperiment,
#'     TFMotifList = exampleTFMotifList
#' )
#'
#' ## This example uses the example MultiAssayExperiment provided in the
#' ## TENET.ExperimentHub package to perform motif overlapping for only
#' ## hypomethylated RE DNA methylation sites linked to the FOXA1 and ESR1
#' ## genes using the motifs for each gene specified in the
#' ## exampleTFMotifList. Gene names and locations, and the locations of RE
#' ## DNA methylation sites, will be retrieved from the rowRanges of the
#' ## 'expression' and 'methylation' SummarizedExperiment objects in the
#' ## example MultiAssayExperiment. Regions within 50 bp of linked RE DNA
#' ## methylation sites will be checked for motifs, and a similarity
#' ## threshold of 80% will be used to identify motifs. The analysis will
#' ## be performed using 8 CPU cores.
#'
#' ## Load the example TENET MultiAssayExperiment object
#' ## from the TENET.ExperimentHub package
#' exampleTENETMultiAssayExperiment <-
#'     TENET.ExperimentHub::exampleTENETMultiAssayExperiment()
#'
#' ## Use the example dataset to do the motif searching
#' returnValue <- step7LinkedDNAMethylationSitesMotifSearching(
#'     TENETMultiAssayExperiment = exampleTENETMultiAssayExperiment,
#'     TFMotifList = exampleTFMotifList,
#'     hypermethGplusAnalysis = FALSE,
#'     distanceFromREDNAMethylationSites = 50,
#'     matchPWMMinScore = "80%",
#'     coreCount = 8
#' )
#'
#' ## This final example uses the example MultiAssayExperiment provided in
#' ## the TENET.ExperimentHub package to perform motif overlapping for a
#' ## just a pair of specific DNA methylation sites, using the FOXA1 and
#' ## MYBL2 motifs, as well as motifs for all human transcription factors
#' ## found in the SwissRegulon database accessed by the MotifDb::query()
#' ## function. As before, Gene names and locations, and the locations of RE
#' ## DNA methylation sites, will be retrieved from the rowRanges of the
#' ## 'expression' and 'methylation' SummarizedExperiment objects in the
#' ## example MultiAssayExperiment. Regions within 100 bp of linked RE DNA
#' ## methylation sites will be checked for motifs, and a similarity
#' ## threshold of 75% will be used to identify motifs. The analysis will
#' ## be performed using one CPU core.
#'
#' ## Create a new list of example PWMs
#' exampleTFMotifList2 <- list(
#'     "FOXA1" = MotifDb::query(MotifDb::MotifDb, "FOXA1")[[3]],
#'     "MYBL2" = MotifDb::query(MotifDb::MotifDb, "MYBL2")[[5]]
#' )
#'
#' ## Load the example TENET MultiAssayExperiment object
#' ## from the TENET.ExperimentHub package
#' exampleTENETMultiAssayExperiment <-
#'     TENET.ExperimentHub::exampleTENETMultiAssayExperiment()
#'
#' ## Use the example dataset to do the motif searching
#' returnValue <- step7LinkedDNAMethylationSitesMotifSearching(
#'     TENETMultiAssayExperiment = exampleTENETMultiAssayExperiment,
#'     hypermethGplusAnalysis = FALSE,
#'     hypomethGplusAnalysis = FALSE,
#'     DNAMethylationSites = c("cg04134755", "cg10216151"),
#'     andStrings = c("Hsapiens", "SwissRegulon"),
#'     orStrings = "humanTranscriptionFactors",
#'     TFMotifList = exampleTFMotifList,
#'     useOnlyDNAMethylationSitesLinkedToTFs = FALSE
#' )
step7LinkedDNAMethylationSitesMotifSearching <- function(
    TENETMultiAssayExperiment,
    hypermethGplusAnalysis = TRUE,
    hypomethGplusAnalysis = TRUE,
    DNAMethylationSites = NA,
    distanceFromREDNAMethylationSites = 100,
    GRangesToSearch = NA,
    andStrings = NULL,
    orStrings = NULL,
    notStrings = NULL,
    TFMotifList,
    useOnlyDNAMethylationSitesLinkedToTFs = TRUE,
    geneAnnotationDataset = NA,
    DNAMethylationArray = NA,
    matchPWMMinScore = "75%",
    coreCount = 1) {
    ## Return an error message if the input MultiAssayExperiment is invalid
    .validateMultiAssayExperiment(
        TENETMultiAssayExperiment,
        needGeneName = is.na(geneAnnotationDataset)
    )

    ## Validate that the supplied GRangesToSearch is a GRanges object, if
    ## specified
    if (!.isSingleNA(GRangesToSearch)) {
        if (!inherits(GRangesToSearch, "GRanges")) {
            .stopNoCall(
                "The object given as the GRangesToSearch argument is not a ",
                "a GRanges object."
            )
        }
    }

    ## Get methylation site IDs and names from the MAE, or methylation array if
    ## provided, and validate that the DNAMethylationSiteID column is present if
    ## analysis types where it is needed are specified
    methSiteIDDF <- .getMethSiteIDsAndLocations(
        TENETMultiAssayExperiment, DNAMethylationArray
    )

    if (!.isSingleNA(DNAMethylationSites) ||
        hypomethGplusAnalysis ||
        hypermethGplusAnalysis
    ) {
        if (!("DNAMethylationSiteID" %in% colnames(methSiteIDDF))) {
            .warningNoCall(
                "The column named 'DNAMethylationSiteID' was not found in the ",
                "specified DNAMethylationArray. Thus, DNA methylation sites ",
                "linked to the TFs specified in the names of the TFMotifList ",
                "cannot be identified as the locations of the DNA methylation ",
                "sites cannot be assessed. Regions provided in the ",
                "GRangesToSearch object can still be assessed however."
            )
        }
    }

    ## Check that the correct combination of arguments have been specified if
    ## useOnlyDNAMethylationSitesLinkedToTFs is set to TRUE.
    if (useOnlyDNAMethylationSitesLinkedToTFs) {
        ## First check that at least one of hypermethGplusAnalysis or
        ## hypomethGplusAnalysis is TRUE
        if (!any(hypermethGplusAnalysis, hypomethGplusAnalysis)) {
            .stopNoCall(
                "If useOnlyDNAMethylationSitesLinkedToTFs is set to TRUE, at ",
                "least one of the hypermethGplusAnalysis or ",
                "hypomethGplusAnalysis arguments must also be TRUE in order ",
                "to identify DNA methylation sites which are linked to TFs by ",
                "TENET for motif searching."
            )
        }

        ## Next check that the DNAMethylationSites, andStrings, and orStrings
        ## arguments are all NA (if they aren't set, they will default to NA)
        if (!all(is.na(c(DNAMethylationSites, andStrings, orStrings)))) {
            .stopNoCall(
                "If useOnlyDNAMethylationSitesLinkedToTFs is set to TRUE, the ",
                "DNAMethylationSites, andStrings, and orStrings arguments ",
                "must be NA, as no additional DNA methylation sites or TFs ",
                "can be analyzed beyond the TFs specified by the user in the ",
                "TFMotifList argument and the hyper- or hypomethyalted DNA ",
                "methylation sites linked to those TFs by TENET analyses."
            )
        }

        ## Validate the analysis types and get a vector of the ones selected,
        ## since the user will need to select at least one of these
        analysisTypes <- .validateAnalysisTypes(
            hypermethGplusAnalysis, hypomethGplusAnalysis
        )

        ## Then, check that a TFMotifList has been provided
        if (missing(TFMotifList)) {
            .stopNoCall(
                "The TFMotifList parameter must be specified if ",
                "useOnlyDNAMethylationSitesLinkedToTFs is set to TRUE."
            )
        }

        ## Also check that the TFMotifList has names for each PWM, as these will
        ## be used to find linked probes to each TF
        if (is.null(names(TFMotifList))) {
            .stopNoCall(
                "The names of the TFs for each motif PWM must be given as the ",
                "names of the TFMotifList so DNA methylation sites linked to ",
                "the TFs by TENET can be identified to search."
            )
        }

        ## Get gene IDs and names from the MAE, or gene annotation dataset if
        ## provided and validate if the geneID and geneName columns are present
        geneIDDF <- .getGeneIDsAndNames(
            TENETMultiAssayExperiment, geneAnnotationDataset
        )

        if (!(("geneID" %in% colnames(geneIDDF)) &&
            ("geneName" %in% colnames(geneIDDF))
        )) {
            .stopNoCall(
                "Columns named 'geneID' and 'geneName' were not found in the ",
                "specified geneAnnotationDataset. Thus, DNA methylation sites ",
                "linked to the TFs specified in the names of the TFMotifList ",
                "cannot be identified as the names/IDs of the TFs cannot be ",
                "assessed."
            )
        }
    } else {
        ## Validate the analysis types and get a vector of the ones selected,
        ## allowing the user to select no analysis types
        analysisTypes <- .validateAnalysisTypes(
            hypermethGplusAnalysis, hypomethGplusAnalysis,
            allowNone = TRUE
        )
    }

    ## Let's start by assembling a list of PWMs for TFs we want to search

    ## First, check if the user has defined a TFMotifList and ensure it is
    ## properly formatted. If they haven't specified an object create an empty
    ## list to append to

    ## First ensure the PWMs given are properly formatted
    ## Should be matrices and have 4 rows
    if (!.isSingleNA(TFMotifList)) {
        if (!all(unlist(
            parallel::mclapply(
                TFMotifList,
                FUN = .internalPWMvalidator,
                mc.cores = coreCount
            )
        ))
        ) {
            .stopNoCall(
                "At least one PWM given in the TFMotifList is improperly ",
                "formatted. PWMs must be matrices with 4 rows, representing ",
                "the position weights for A, C, G, and T respectively for ",
                "each base in the TF motif (in the columns)."
            )
        }
    } else {
        ## Create an empty TFMotifList to append to later
        TFMotifList <- list()
    }

    ## Now if the user has specified either andStrings or orStrings, use those
    ## to find additional PWMs using the MotifDb::query function
    if (any(c(!is.null(andStrings), !is.null(orStrings)))) {
        ## First, if the user has included "humanTranscriptionFactors" in the
        ## orStrings argument, replace it with the names of all the identified
        ## human TF names identified previously
        if ("humanTranscriptionFactors" %in% orStrings) {
            ## Remove the "humanTranscriptionFactors" element since MotifDb does
            ## not understand it
            orStrings <- orStrings[orStrings != "humanTranscriptionFactors"]

            ## Create an environment to store data from the TENET package
            TENETDataEnv <- new.env(parent = emptyenv())

            ## Get the IDs and names of the human transcription factors

            ## Load the humanTranscriptionFactorDb object from data
            utils::data(
                "humanTranscriptionFactorDb",
                package = "TENET",
                envir = TENETDataEnv
            )
            humanTFDb <- TENETDataEnv$humanTranscriptionFactorDb

            ## Get the confirmed human TFs
            humanTFDb <- humanTFDb[humanTFDb$Is.TF. == "Yes", ]

            ## Add both the ensembl ID and names for the human TFs
            orStrings <- c(orStrings, humanTFDb$HGNC.symbol)
        }

        ## Now use the andStrings and orStrings arguments to find TF PWMs
        motifDbQueryPWMList <- as.list(
            MotifDb::query(
                MotifDb::MotifDb,
                andStrings = andStrings,
                orStrings = orStrings,
                notStrings = notStrings
            )
        )

        ## Add the motifDbQueryPWMList to the TFMotifList and ensure the names
        ## are unique
        TFMotifList <- c(TFMotifList, motifDbQueryPWMList)
        names(TFMotifList) <- make.unique(names(TFMotifList))
    }

    ## Do a quick check to ensure there are TFs to search on
    if (length(TFMotifList) == 0) {
        .stopNoCall(
            "No valid TF PWM motifs were found. Please check the list given ",
            "as TFMotifList and check settings of other function arguments ",
            "and try again."
        )
    }

    ## Now we need to create a GRanges object with regions we want to perform TF
    ## motif searching on

    ## First, since they'll be used in the specific
    ## useOnlyDNAMethylationSitesLinkedToTFs analysis, let's create a data frame
    ## with step 5 results for the selected analysisTypes.

    ## Next, if there are analysis types selected, get the step 5 results from
    ## each analysis type, and combine them with an additional column noting if
    ## the results are hyper or hypo
    if (!is.null(analysisTypes)) {
        ## Create an empty object
        hyperHypoResultsList <- list()

        for (i in seq_along(analysisTypes)) {
            quadrantResultsName <- paste0(analysisTypes[i], "methGplusResults")

            ## Ensure the quadrant's results are present in step 5
            .ensureStepPresent(
                TENETMultiAssayExperiment,
                stepName = "step5OptimizeLinks",
                substepName = quadrantResultsName
            )

            ## Load the quadrant's significant links from step 5
            quadrantSigLinkZScores <- TENETMultiAssayExperiment@metadata$
                step5OptimizeLinks[[quadrantResultsName]]

            ## Add a column noting the analysis types of this result
            quadrantSigLinkZScores$quadrant <- quadrantResultsName

            ## Add dummy column names to allow the data frames from hyper- and
            ## hypomethylated results to rbind
            colnames(quadrantSigLinkZScores) <- LETTERS[
                seq_along(colnames(quadrantSigLinkZScores))
            ]

            ## Add them to the list
            hyperHypoResultsList[[i]] <- quadrantSigLinkZScores
        }

        ## Bind the list into a single data frame
        hyperHypoResultsDF <- do.call(rbind, hyperHypoResultsList)
    } else {
        ## Create an empty data frame which mimics the ones that would hold
        ## results so we can combine it as normal with later sites specified by
        ## the user. Note: if the step 5 results are ever adjusted, this will
        ## need to be updated too, particularly the ncol value
        hyperHypoResultsDF <- data.frame(matrix(nrow = 0, ncol = 14))
        colnames(hyperHypoResultsDF) <- LETTERS[seq_len(14)]
    }

    ## Start with the case where useOnlyDNAMethylationSitesLinkedToTFs is TRUE.
    if (useOnlyDNAMethylationSitesLinkedToTFs) {
        ## First, for each TF name/ID in the names of the TFMotifList, get the
        ## corresponding ID/name
        geneNameCheckDF <- data.frame(
            "geneIDsBool" = (names(TFMotifList) %in% geneIDDF$geneID),
            "geneNamesBool" = (names(TFMotifList) %in% geneIDDF$geneName),
            stringsAsFactors = FALSE
        )
        rownames(geneNameCheckDF) <- names(TFMotifList)

        ## Mapply the .internalGeneNameIDFromBool function to this df to
        ## get both gene names and IDs no matter which notation was given in
        ## TFMotifList
        tfMotifListNamesDF <- as.data.frame(
            t(parallel::mcmapply(
                .internalGeneNameIDFromBool,
                IDBool = geneNameCheckDF$geneIDsBool,
                nameBool = geneNameCheckDF$geneNamesBool,
                TFMotifListName = rownames(geneNameCheckDF),
                MoreArgs = list("geneIDDF" = geneIDDF),
                mc.cores = coreCount
            ))
        )
        colnames(tfMotifListNamesDF) <- c("geneIDs", "geneNames")

        ## Get elements which are NA (meaning there is a missing TFMotifList
        ## name or the gene's name/ID can't be identified) and return a warning
        ## message about them
        elementNumberNotPresent <- which(is.na(tfMotifListNamesDF$geneIDs))

        .warningNoCall(
            "The gene names or IDs of TFs at position(s) ",
            paste(elementNumberNotPresent, sep = ", "),
            " were not identified in the specified geneAnnotationDataset. ",
            "Motif searching will still be performed on these TF PWMs, but ",
            "regions of DNA methylation sites linked to that TF cannot be ",
            "found and will not be included in the motif searching."
        )

        ## Remove any NA genes from the list
        tfMotifListNamesDF <- tfMotifListNamesDF[
            !is.na(tfMotifListNamesDF$geneIDs),
        ]

        ## Now, subset the hyperHypoResultsDF to just the DNA methylation sites
        ## linked to at least one of the named TFs
        hyperHypoResultsDF <- hyperHypoResultsDF[
            hyperHypoResultsDF$A %in% tfMotifListNamesDF$geneIDs,
        ]

        ## Next, check that the remaining DNA methylation sites are found in the
        ## methSiteIDdf so that we can get coordinates for those later
        MethSitesLinkedToGenesNotFound <- hyperHypoResultsDF[
            !(hyperHypoResultsDF$B %in% methSiteIDDF$DNAMethylationSiteID),
            "B"
        ]

        if (length(MethSitesLinkedToGenesNotFound) > 0) {
            .warningNoCall(
                "Genomic locations for the following RE DNA methylation sites ",
                "linked to TF genes specified in the TFMotifList were not ",
                "identified: ",
                paste(
                    MethSitesLinkedToGenesNotFound,
                    collapse = ", "
                ),
                ". Please check the specified DNAMethylationArray argument. ",
                "These methylation sites have been excluded from this analysis."
            )
        }

        ## Remove those sites from the analysis now
        hyperHypoResultsDF <- hyperHypoResultsDF[
            !(hyperHypoResultsDF$B %in% MethSitesLinkedToGenesNotFound),
        ]

        ## Issue an error if no remaining sites are identified
        if (nrow(hyperHypoResultsDF) == 0) {
            .stopNoCall(
                "No valid DNA methylation sites linked to TFs in the ",
                "TFMotifList were identified. Please ensure there are DNA ",
                "methylation sites linked to the TFs specified in the ",
                "TFMotifList, and these probes have valid coordinates."
            )
        }

        ## Now for each TF and analysis type selected, let's create a vector
        ## going over each DNA methylation site and analysis type and see if it
        ## is linked to the TF in that type
        linkedTFTypes <- paste(
            rep(
                paste(
                    tfMotifListNamesDF$geneIDs,
                    tfMotifListNamesDF$geneNames,
                    sep = "_"
                ),
                each = length(analysisTypes)
            ),
            paste0(analysisTypes, "methGplusResults"),
            sep = "_"
        )

        linkedTFTypesDNAMethylationBoolList <- list()

        for (j in seq_along(linkedTFTypes)) {
            ## Split the name of the object into the TF ID and the analysis type
            ## again
            individualTFID <- strsplit(linkedTFTypes[j], "_")[[1]][1]
            individualAnalysisType <- strsplit(linkedTFTypes[j], "_")[[1]][3]

            ## For each unique CpG, isolate the rows for that CpG and identify
            ## if the unique TF ID and analysis type combo is found for that TF.
            ## This is done because we want a new data frame with just results
            ## for each unique TF
            siteSpecificVector <- NULL

            for (k in seq_along(unique(hyperHypoResultsDF$B))) {
                ## Isolate the results for the given DNA methylation site
                hyperHypoResultsDFSiteSpecific <- hyperHypoResultsDF[
                    hyperHypoResultsDF$B == unique(hyperHypoResultsDF$B)[k],
                ]

                ## For each DNA methylation site in hyperHypoResultsDF, create a
                ## boolean noting if it matches the TF and analysis type and add
                ## it to the list
                siteSpecificVector <- c(
                    siteSpecificVector,
                    ifelse(
                        individualTFID %in% hyperHypoResultsDFSiteSpecific$A,
                        ifelse(
                            individualAnalysisType %in%
                                hyperHypoResultsDFSiteSpecific$N,
                            TRUE,
                            FALSE
                        ),
                        FALSE
                    )
                )
            }

            ## Add the vector to the list
            linkedTFTypesDNAMethylationBoolList[[j]] <- siteSpecificVector
        }

        ## cbind the linkedTFTypesDNAMethylationBoolList into a single data
        ## frame noting whether each TF/analysis type is linked to each DNA
        ## methylation site in the motif analysis
        linkedUniqueDNAMethylationSitesTFOverlap <- as.data.frame(
            do.call(
                cbind,
                linkedTFTypesDNAMethylationBoolList
            )
        )

        ## Set the rownames to be the names of the methylation sites, then add
        ## the gene IDs/names and analysis types as the columns (adding 'Linked'
        ## to the end)
        rownames(linkedUniqueDNAMethylationSitesTFOverlap) <- unique(
            hyperHypoResultsDF$B
        )
        colnames(linkedUniqueDNAMethylationSitesTFOverlap) <- paste0(
            linkedTFTypes,
            "Linked"
        )

        ## Lastly, create a GRanges object with the range around each DNA
        ## methylation site for the actual TF motif searching
        ## GRangesToSearch should be NA (empty) so we can overwrite it without
        ## problem

        ## Create a subset of the methSiteIDDF with just the DNA methylation
        ## sites of interest and sort it to match the order of sites in
        ## hyperHypoResultsDF
        methSiteIDDFHyperHypoResults <- methSiteIDDF[
            methSiteIDDF$DNAMethylationSiteID %in% hyperHypoResultsDF$B,
        ]

        methSiteIDDFHyperHypoResults <- methSiteIDDFHyperHypoResults[
            order(
                match(
                    methSiteIDDFHyperHypoResults$DNAMethylationSiteID,
                    hyperHypoResultsDF$B
                )
            ),
        ]

        ## Then adjust the size of these regions according to the
        ## distanceFromREDNAMethylationSites value
        ## NOTE: This makes the assumption the start and end postions of DNA
        ## methylation sites are in the 3rd and 4th columns of the methSiteIDDF
        methSiteIDDFHyperHypoResults[, 3] <- c(
            methSiteIDDFHyperHypoResults[, 3] -
                distanceFromREDNAMethylationSites
        )

        methSiteIDDFHyperHypoResults[, 4] <- c(
            methSiteIDDFHyperHypoResults[, 4] +
                distanceFromREDNAMethylationSites
        )

        ## Now create the GRanges object
        GRangesToSearch <- GenomicRanges::makeGRangesFromDataFrame(
            df = methSiteIDDFHyperHypoResults,
            keep.extra.columns = FALSE,
            ignore.strand = TRUE
        )
    } else {
        ## Create a dummy linkedUniqueDNAMethylationSitesTFOverlap variable to
        ## avoid if statements in later code
        linkedUniqueDNAMethylationSitesTFOverlap <- NA

        ## When useOnlyDNAMethylationSitesLinkedToTFs is set to FALSE, we don't
        ## have to give as much concern to checking the names of the genes in
        ## the TFMotifList, but we do have to collate regions of interest from
        ## several sources now

        ## First, check that all the DNA methylation sites from the hyper- and
        ## hypometh analyses are found in the methSiteIDdf so that we can get
        ## coordinates for those later
        MethSitesLinkedToGenesNotFound <- hyperHypoResultsDF[
            !(hyperHypoResultsDF$B %in% methSiteIDDF$DNAMethylationSiteID),
            "B"
        ]

        if (length(MethSitesLinkedToGenesNotFound) > 0) {
            .warningNoCall(
                "Genomic locations for the following RE DNA methylation sites ",
                "were not identified: ",
                paste(
                    MethSitesLinkedToGenesNotFound,
                    collapse = ", "
                ),
                ". Please check the specified DNAMethylationArray argument. ",
                "These methylation sites have been excluded from this analysis."
            )
        }

        ## Remove those sites from the analysis now
        hyperHypoResultsDF <- hyperHypoResultsDF[
            !(hyperHypoResultsDF$B %in% MethSitesLinkedToGenesNotFound),
        ]

        ## To get sites, create a subset of the methSiteIDDF with just the DNA
        ## methylation sites of interest and sort it to match the order of sites
        ## in hyperHypoResultsDF
        methSiteIDDFHyperHypoResults <- methSiteIDDF[
            methSiteIDDF$DNAMethylationSiteID %in% hyperHypoResultsDF$B,
        ]

        methSiteIDDFHyperHypoResults <- methSiteIDDFHyperHypoResults[
            order(
                match(
                    methSiteIDDFHyperHypoResults$DNAMethylationSiteID,
                    hyperHypoResultsDF$B
                )
            ),
        ]

        ## Then adjust the size of these regions according to the
        ## distanceFromREDNAMethylationSites value
        ## NOTE: This makes the assumption the start and end postions of DNA
        ## methylation sites are in the 3rd and 4th columns of the methSiteIDDF
        methSiteIDDFHyperHypoResults[, 3] <- c(
            methSiteIDDFHyperHypoResults[, 3] -
                distanceFromREDNAMethylationSites
        )

        methSiteIDDFHyperHypoResults[, 4] <- c(
            methSiteIDDFHyperHypoResults[, 4] +
                distanceFromREDNAMethylationSites
        )

        ## Create the GRanges object if it is NA (hasn't been given by the user)
        ## or add these ranges to the existing ranges specified by the user
        if (.isSingleNA(GRangesToSearch)) {
            GRangesToSearch <- GenomicRanges::makeGRangesFromDataFrame(
                df = methSiteIDDFHyperHypoResults,
                keep.extra.columns = FALSE,
                ignore.strand = TRUE
            )
        } else {
            GRangesToSearch <- c(
                GenomicRanges::makeGRangesFromDataFrame(
                    df = methSiteIDDFHyperHypoResults,
                    keep.extra.columns = FALSE,
                    ignore.strand = TRUE
                ),
                GRangesToSearch
            )
        }

        ## Next, let's add any custom regions the user has specified with the
        ## DNAMethylationSites argument, if any
        if (!.isSingleNA(DNAMethylationSites)) {
            ## First, check that all the DNA methylation sites the user has
            ## selected are found in the methSiteIDdf so that we can get
            ## coordinates for those later
            userDNAMethylationSitesNotFound <- DNAMethylationSites[
                !(DNAMethylationSites %in% methSiteIDDF$DNAMethylationSiteID)
            ]

            if (length(userDNAMethylationSitesNotFound) > 0) {
                .warningNoCall(
                    "Genomic locations for the following sites specified in ",
                    "userDNAMethylationSitesNotFound were not identified: ",
                    paste(
                        userDNAMethylationSitesNotFound,
                        collapse = ", "
                    ),
                    ". Please check the values specified in ",
                    "DNAMethylationSites. ",
                    "For now, these methylation sites have been excluded ",
                    "from this analysis."
                )
            }

            ## Remove those sites from the analysis now
            DNAMethylationSites <- DNAMethylationSites[
                !(DNAMethylationSites %in% userDNAMethylationSitesNotFound)
            ]

            ## To get sites, create a subset of the methSiteIDDF with just the
            ## DNA methylation sites of interest and sort it to match the order
            ## of sites in the user's DNAMethylationSites
            DNAMethylationSitesResults <- methSiteIDDF[
                methSiteIDDF$DNAMethylationSiteID %in% DNAMethylationSites,
            ]

            DNAMethylationSitesResults <- DNAMethylationSitesResults[
                order(
                    match(
                        DNAMethylationSitesResults$DNAMethylationSiteID,
                        DNAMethylationSites
                    )
                ),
            ]

            ## Then adjust the size of these regions according to the
            ## distanceFromREDNAMethylationSites value.
            ## NOTE: This makes the assumption the start and end postions of DNA
            ## methylation sites are in the 3rd and 4th columns of the
            ## methSiteIDDF.
            DNAMethylationSitesResults[, 3] <- c(
                DNAMethylationSitesResults[, 3] -
                    distanceFromREDNAMethylationSites
            )

            DNAMethylationSitesResults[, 4] <- c(
                DNAMethylationSitesResults[, 4] +
                    distanceFromREDNAMethylationSites
            )

            ## Now create the GRanges object if it is NA (hasn't been given by
            ## the user). Otherwise, add these ranges to the existing ranges
            ## specified by the user (after existing ones)
            if (.isSingleNA(GRangesToSearch)) {
                GRangesToSearch <- GenomicRanges::makeGRangesFromDataFrame(
                    df = DNAMethylationSitesResults,
                    keep.extra.columns = FALSE,
                    ignore.strand = TRUE
                )
            } else {
                GRangesToSearch <- c(
                    GRangesToSearch,
                    GenomicRanges::makeGRangesFromDataFrame(
                        df = DNAMethylationSitesResults,
                        keep.extra.columns = FALSE,
                        ignore.strand = TRUE
                    )
                )
            }
        }
    }

    ## Do a check to make sure there are GRangesToSearch remaining and issue an
    ## error if not
    if (length(GRangesToSearch) == 0) {
        .stopNoCall(
            "No ranges to perform motif searching on were identified. Please ",
            "check the `hypermethGplusAnalysis`, `hypomethGplusAnalysis`, ",
            "`DNAMethylationSites`, and `GRangesToSearch` arguments, as well ",
            "as the settings for the `useOnlyDNAMethylationSitesLinkedToTFs` ",
            "and `DNAMethylationArray` arguments to ensure regions aren't ",
            "being excluded."
        )
    }

    ## Check if there are any duplicated rownames, and set them to be unique if
    ## there are
    if (any(duplicated(names(GRangesToSearch)))) {
        .warningNoCall(
            "Duplicated DNA methylation sites were identified. All regions ",
            "will be kept, but the names of the duplicated regions will be ",
            "made unique."
        )

        names(GRangesToSearch) <- make.unique(names(GRangesToSearch))
    }

    ## Now, go across each of the TFs in the TFMotifList and do motif binding
    ## predictions on each of the sites in GRangesToSearch.

    ## First create a list which will hold seqLogos for each TF
    TFSeqLogoList <- list()

    ## Also create a list for storing the motif occurences data frames as well
    ## as the counts
    listOfMethSiteInfo <- list()
    listOfMethSiteCountInfo <- list()

    for (l in seq_along(TFMotifList)) {
        ## Start by creating a seqLogo of the motif PWM
        .newInvisibleRecordablePlot(height = 3.5, width = 6)
        seqLogo::seqLogo(TFMotifList[[l]])
        thisSeqLogo <- .recordTENETSavedSizePlot()

        ## Close the plot
        grDevices::dev.off()

        ## Add it to the relevant list
        TFSeqLogoList[[l]] <- thisSeqLogo

        ## Do the analysis of predicted motifs for the given TF across the
        ## DNA methylation sites
        TFSpecificListOfMethSiteInfo <- unname(parallel::mclapply(
            X = seq_along(GRangesToSearch),
            FUN = .findMotifSurroundingMethSite,
            genome = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38,
            GRangesObject = GRangesToSearch,
            motifPWM = TFMotifList[[l]],
            matchPWMMinScore = matchPWMMinScore,
            mc.cores = coreCount
        ))

        ## Transform the results into a data frame
        TFResultsDF <- as.data.frame(
            do.call(rbind, TFSpecificListOfMethSiteInfo)
        )

        ## Add a final column noting the name of the TF
        TFResultsDF$X6 <- rep(
            names(TFMotifList)[l],
            nrow(TFResultsDF)
        )

        ## Add that data frame to the listOfMethSiteInfo list
        listOfMethSiteInfo[[l]] <- TFResultsDF

        ## Now get a vector noting the number of motifs for the given
        ## TF found near each site and add it to the listOfMethSiteCountInfo
        ## list
        givenTFPerSiteCount <- NULL

        for (m in seq_along(GRangesToSearch)) {
            ## If the DNA methylation site had a TF linked in it,
            ## Get that section of the dataset and tabulate how many motifs of
            ## each TF were found
            if (names(GRangesToSearch)[m] %in% TFResultsDF$X1) {
                ## Get that section of the data frame
                siteTFResultsDF <- TFResultsDF[
                    TFResultsDF$X1 == names(GRangesToSearch)[m],
                ]

                ## Get a count of the number of times the given TF is found
                ## It's just the length of this section of the dataset, since
                ## it's done individually per TF
                givenTFPerSiteCount <- c(
                    givenTFPerSiteCount,
                    nrow(siteTFResultsDF)
                )
            } else {
                givenTFPerSiteCount <- c(
                    givenTFPerSiteCount,
                    0
                )
            }
        }

        ## Add the count to its list
        listOfMethSiteCountInfo[[l]] <- givenTFPerSiteCount
    }

    ## Create a final data frame of the TF binding predictions and counts from
    ## the listOfMethSiteInfo andlistOfMethSiteCountInfo objects
    TFBindingPredictionResultsDF <- do.call(rbind, listOfMethSiteInfo)
    colnames(TFBindingPredictionResultsDF) <- c(
        "DNAMethylationSiteID",
        "Chromosome",
        "MotifStart",
        "MotifEnd",
        "MotifSequence",
        "TFName"
    )

    TFBindingPredictionCountsDF <- as.data.frame(
        do.call(rbind, listOfMethSiteCountInfo),
        stringsAsFactors = FALSE
    )
    colnames(TFBindingPredictionCountsDF) <- paste0(
        names(GRangesToSearch),
        "PredictedMotifCount"
    )
    rownames(TFBindingPredictionCountsDF) <- names(TFMotifList)

    ## Assemble a list of all the objects we want to return to the user
    resultsList <- list(
        "DNAMethylationSitesGRanges" = GRangesToSearch,
        "TFMotifPWMList" = TFMotifList,
        "TFMotifSeqLogoList" = TFSeqLogoList,
        "DNAMethylationSitesMotifOccurrences" = TFBindingPredictionResultsDF,
        "totalMotifOccurrencesPerDNAMethylationSite" =
            TFBindingPredictionCountsDF,
        "linkedUniqueDNAMethylationSitesTFOverlap" =
            linkedUniqueDNAMethylationSitesTFOverlap
    )

    ## Add the results list to the MultiAssayExperiment
    TENETMultiAssayExperiment@metadata$
        step7LinkedDNAMethylationSitesMotifSearching <- resultsList

    return(TENETMultiAssayExperiment)
}
