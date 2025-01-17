#' Run the step 1 through step 6 functions with default arguments
#'
#' This function combines the step 1 through step 6 functions
#' (`step1MakeExternalDatasets`, `step2GetDifferentiallyMethylatedSites`,
#' `step3GetAnalysisZScores`,
#' `step4SelectMostSignificantLinksPerDNAMethylationSite`, `step5OptimizeLinks`,
#' and `step6DNAMethylationSitesPerGeneTabulation`) into a single function.
#' Arguments for this function generally reflect the arguments of the component
#' functions without clearly defined defaults, with the exception of the
#' `step1MakeExternalDatasets` function where all arguments have been included
#' to support all of the options available to the user to define regions with
#' relevant regulatory elements. All remaining arguments of the component
#' functions are set to their default values.
#'
#' @param TENETMultiAssayExperiment Specify a MultiAssayExperiment object
#' containing expression and methylation SummarizedExperiment objects,
#' such as one created by the TCGADownloader function. Coordinates for
#' genes and DNA methylation sites should be included in the rowRanges
#' of their respective SummarizedExperiment objects and should be annotated
#' to the human hg38 genome.
#' An argument of all functions except `step1MakeExternalDatasets`.
#' @param extHM To use custom histone modification datasets, specify a path
#' to a directory containing .bed, .narrowPeak, .broadPeak, and/or
#' .gappedPeak files with these datasets. The files may optionally be compressed
#' (.gz/.bz2/.xz). Otherwise, specify NA or do not specify this argument.
#' An argument of the `step1MakeExternalDatasets` function.
#' @param extNDR To use custom open chromatin or NDR datasets, specify a
#' path to a directory containing .bed, .narrowPeak, .broadPeak, and/or
#' .gappedPeak files with these datasets. The files may optionally be compressed
#' (.gz/.bz2/.xz). Otherwise, specify NA or do not specify this argument.
#' An argument of the `step1MakeExternalDatasets` function.
#' @param consensusEnhancer Set to TRUE to use the consensus enhancer data
#' included in TENET.AnnotationHub. Defaults to TRUE.
#' An argument of the `step1MakeExternalDatasets` function.
#' @param consensusPromoter Set to TRUE to use the consensus promoter data
#' included in TENET.AnnotationHub. Defaults to FALSE.
#' An argument of the `step1MakeExternalDatasets` function.
#' @param consensusNDR Set to TRUE to use the consensus open chromatin
#' data included in TENET.AnnotationHub. Defaults to TRUE.
#' An argument of the `step1MakeExternalDatasets` function.
#' @param publicEnhancer Set to TRUE to use the preprocessed publicly available
#' enhancer (H3K27ac) datasets included in TENET.AnnotationHub. If set to TRUE,
#' `cancerType` must be specified. Defaults to FALSE.
#' An argument of the `step1MakeExternalDatasets` function.
#' @param publicPromoter Set to TRUE to use the preprocessed publicly available
#' promoter (H3K4me3) datasets included in TENET.AnnotationHub. If set to TRUE,
#' `cancerType` must be specified. Defaults to FALSE.
#' An argument of the `step1MakeExternalDatasets` function.
#' @param publicNDR Set to TRUE to use the preprocessed publicly available
#' open chromatin (ATAC-seq, DNase-seq) datasets included in
#' TENET.AnnotationHub. If set to TRUE, `cancerType` must be specified.
#' Defaults to FALSE. An argument of the `step1MakeExternalDatasets` function.
#' @param cancerType If `publicEnhancer`, `publicPromoter`, and/or `publicNDR`
#' is TRUE, specify a vector of cancer types from 'BLCA', 'BRCA', 'COAD',
#' 'ESCA', 'HNSC', 'KIRP', 'LIHC', 'LUAD', 'LUSC', and 'THCA' to include the
#' public data relevant to those cancer types. Defaults to NA.
#' An argument of the `step1MakeExternalDatasets` function.
#' @param ENCODEPLS Set to TRUE to use the ENCODE promoter-like elements
#' dataset included in TENET.AnnotationHub. Defaults to FALSE.
#' An argument of the `step1MakeExternalDatasets` function.
#' @param ENCODEpELS Set to TRUE to use the ENCODE proximal enhancer-like
#' elements dataset included in TENET.AnnotationHub. Defaults to FALSE.
#' An argument of the `step1MakeExternalDatasets` function.
#' @param ENCODEdELS Set to TRUE to use the ENCODE distal enhancer-like
#' elements dataset included in TENET.AnnotationHub. Defaults to FALSE.
#' An argument of the `step1MakeExternalDatasets` function.
#' @param assessPromoter Set to TRUE to identify DNA methylation sites
#' that mark promoter regions or FALSE to identify distal enhancer regions.
#' Defaults to FALSE.
#' An argument of the `step2GetDifferentiallyMethylatedSites` function.
#' @param TSSDist Specify a positive integer distance in base pairs to any
#' transcription start site within which DNA
#' methylation sites are considered promoter DNA methylation sites. DNA
#' methylation sites outside of the TSSDist from any transcription start site
#' will be considered enhancer methylation sites. Defaults to 1500.
#' An argument of the `step2GetDifferentiallyMethylatedSites` function.
#' @param minCaseCount Specify a positive integer to be the
#' minimum number of case samples to be considered for the
#' hyper- or hypomethylated groups. Should be less than the total number
#' of case samples.
#' An argument of the `step2GetDifferentiallyMethylatedSites` function.
#' @param coreCount Argument passed as the mc.cores argument to mclapply. See
#' `?parallel::mclapply` for more details. Defaults to 1.
#' Used by the `step3GetAnalysisZScores`,
#' `step4SelectMostSignificantLinksPerDNAMethylationSite`, and
#' `step5OptimizeLinks` functions.
#' @return Returns the created MultiAssayExperiment object containing data from
#' all step 1 through 6 functions.
#' @export
#'
#' @examplesIf interactive()
#' ## This example creates a dataset of putative enhancer regulatory elements
#' ## from consensus datasets and breast invasive carcinoma-relevant sources
#' ## collected in the TENET.AnnotationHub package, then runs the step 2 through
#' ## step 6 TENET functions analyzing RE DNA methylation sites in potential
#' ## enhancer elements located over 1500 bp from transcription start sites
#' ## listed for genes and transcripts in the GENCODE v36 human genome
#' ## annotations, using a minimum case sample count of 5 and one CPU core
#' ## to perform the analysis.
#'
#' ## Load the example TENET MultiAssayExperiment object from the
#' ## TENET.ExperimentHub package
#' exampleTENETMultiAssayExperiment <-
#'     TENET.ExperimentHub::exampleTENETMultiAssayExperiment()
#'
#' ## Use the example TENET MultiAssayExperiment to run the step 1 through
#' ## step 6 TENET functions
#' returnValue <- easyTENET(
#'     TENETMultiAssayExperiment = exampleTENETMultiAssayExperiment,
#'     extHM = NA,
#'     extNDR = NA,
#'     publicEnhancer = TRUE,
#'     publicNDR = TRUE,
#'     cancerType = "BRCA",
#'     ENCODEdELS = TRUE,
#'     minCaseCount = 5
#' )
#'
#' ## This example creates a dataset of putative promoter regulatory elements
#' ## using bed-like files contained in the user's working directory, consensus
#' ## NDR and promoter regions, and regions with promoter-like signatures from
#' ## the ENCODE SCREEN project, but excluding cancer type-specific public
#' ## datasets. This dataset is then used to analyze DNA methylation sites in
#' ## promoter elements within 2000 bp of all transcription start sites
#' ## provided in the MultiAssayExperiment only, identifying alterations found
#' ## in at least 10 samples, and using 8 CPU cores to perform the analysis.
#'
#' ## Load the example TENET MultiAssayExperiment object from the
#' ## TENET.ExperimentHub package
#' exampleTENETMultiAssayExperiment <-
#'     TENET.ExperimentHub::exampleTENETMultiAssayExperiment()
#'
#' ## Use the example TENET MultiAssayExperiment to run the step 1 through
#' ## step 6 TENET functions
#' returnValue <- easyTENET(
#'     TENETMultiAssayExperiment = exampleTENETMultiAssayExperiment,
#'     extHM = ".",
#'     extNDR = ".",
#'     consensusEnhancer = FALSE,
#'     consensusPromoter = TRUE,
#'     ENCODEPLS = TRUE,
#'     assessPromoter = TRUE,
#'     TSSDist = 2000,
#'     minCaseCount = 10,
#'     coreCount = 8
#' )
easyTENET <- function(
    TENETMultiAssayExperiment,
    extHM = NA,
    extNDR = NA,
    consensusEnhancer = TRUE,
    consensusPromoter = FALSE,
    consensusNDR = TRUE,
    publicEnhancer = FALSE,
    publicPromoter = FALSE,
    publicNDR = FALSE,
    cancerType = NA,
    ENCODEPLS = FALSE,
    ENCODEpELS = FALSE,
    ENCODEdELS = FALSE,
    assessPromoter = FALSE,
    TSSDist = 1500,
    minCaseCount,
    coreCount = 1) {
    ## Return an error message if the input MultiAssayExperiment is invalid
    .validateMultiAssayExperiment(TENETMultiAssayExperiment)

    ## Bail early if minCaseCount is missing
    if (missing(minCaseCount)) {
        .stopNoCall("The minCaseCount parameter must be specified.")
    }

    ## Run the step1MakeExternalDatasets function to create a
    ## GRanges object with relevant regulatory elements, passing on the
    ## arguments relevant for that function to it
    REGRangesObject <- TENET::step1MakeExternalDatasets(
        extHM = extHM,
        extNDR = extNDR,
        consensusEnhancer = consensusEnhancer,
        consensusPromoter = consensusPromoter,
        consensusNDR = consensusNDR,
        publicEnhancer = publicEnhancer,
        publicPromoter = publicPromoter,
        publicNDR = publicNDR,
        cancerType = cancerType,
        ENCODEPLS = ENCODEPLS,
        ENCODEpELS = ENCODEpELS,
        ENCODEdELS = ENCODEdELS
    )

    ## Run the step2GetDifferentiallyMethylatedSites function on the
    ## specified TENETMultiAssayExperiment object to identify differentially
    ## methylated RE sites, using default values for arguments not specified in
    ## this wrapper function
    TENETMultiAssayExperiment <- TENET::step2GetDifferentiallyMethylatedSites(
        TENETMultiAssayExperiment = TENETMultiAssayExperiment,
        regulatoryElementGRanges = REGRangesObject,
        assessPromoter = assessPromoter,
        TSSDist = TSSDist,
        minCaseCount = minCaseCount
    )

    ## Run the step3GetAnalysisZScores function on the TENETMultiAssayExperiment
    ## to calculate Z-scores comparing the mean expression of each gene in the
    ## case samples that are hyper- or hypomethylated for each RE DNA
    ## methylation site, using default arguments other than the specified
    ## coreCount value
    TENETMultiAssayExperiment <- TENET::step3GetAnalysisZScores(
        TENETMultiAssayExperiment = TENETMultiAssayExperiment,
        coreCount = coreCount
    )

    ## Run the step4SelectMostSignificantLinksPerDNAMethylationSite function on
    ## the TENETMultiAssayExperiment to select the most significant RE DNA
    ## methylation site-gene links to each RE DNA methylation site, using
    ## default arguments other than the specified coreCount value
    TENETMultiAssayExperiment <-
        TENET::step4SelectMostSignificantLinksPerDNAMethylationSite(
            TENETMultiAssayExperiment = TENETMultiAssayExperiment,
            coreCount = coreCount
        )

    ## Run the step5OptimizeLinks function on the TENETMultiAssayExperiment to
    ## find final RE DNA methylation site-gene links using various optimization
    ## metrics and using default arguments other than the specified coreCount
    ## value
    TENETMultiAssayExperiment <- TENET::step5OptimizeLinks(
        TENETMultiAssayExperiment = TENETMultiAssayExperiment,
        coreCount = coreCount
    )

    ## Run the step6DNAMethylationSitesPerGeneTabulation function on the
    ## TENETMultiAssayExperiment to tabulate the total number of RE DNA
    ## methylation sites linked to each of the genes, using default arguments
    TENETMultiAssayExperiment <-
        TENET::step6DNAMethylationSitesPerGeneTabulation(
            TENETMultiAssayExperiment = TENETMultiAssayExperiment
        )

    ## Return the final TENETMultiAssayExperiment with the results from the
    ## step 1 through step 6 functions
    return(TENETMultiAssayExperiment)
}
