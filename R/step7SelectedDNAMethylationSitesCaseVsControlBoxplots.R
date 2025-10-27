#' Generate boxplots or violin plots comparing the methylation level of the
#' specified RE DNA methylation sites in case and control samples
#'
#' This function takes a vector of RE DNA methylation sites specified by the
#' user and generates boxplots or violin plots displaying the methylation level
#' of each of these DNA methylation sites in the case compared to control
#' samples, along with the results of a Student's t-test comparing the
#' methylation level between these two groups.
#'
#' @param TENETMultiAssayExperiment Specify a MultiAssayExperiment object
#' containing a methylation SummarizedExperiment object, such as one created by
#' the TCGADownloader function.
#' @param DNAMethylationSites Supply a vector of RE DNA methylation site IDs
#' for which to create boxplots or violin plots with the methylation of those RE
#' DNA methylation sites.
#' @param violinPlots Set to TRUE to generate violin plots instead of boxplots.
#' Defaults to FALSE.
#' @param coreCount Argument passed as the mc.cores argument to mcmapply. See
#' `?parallel::mcmapply` for more details. Defaults to 1.
#' @return Returns the MultiAssayExperiment object given as the
#' TENETMultiAssayExperiment argument with an additional list named
#' 'step7SelectedDNAMethylationSitesCaseVsControlBoxplots' in its metadata,
#' which contains boxplots or violin plots comparing the methylation of the RE
#' DNA methylation sites of interest in the case and control samples. The titles
#' of the plots contain the ID of the RE DNA methylation site and the Student's
#' t-test p-value.
#' @export
#'
#' @examplesIf interactive()
#' ## This example uses the example MultiAssayExperiment provided in the
#' ## TENET.ExperimentHub package to generate boxplots for several selected
#' ## RE DNA methylation sites, using one CPU core.
#'
#' ## Load the example TENET MultiAssayExperiment object
#' ## from the TENET.ExperimentHub package
#' exampleTENETMultiAssayExperiment <-
#'     TENET.ExperimentHub::exampleTENETMultiAssayExperiment()
#'
#' ## Use the example dataset to create the boxplots
#' returnValue <- step7SelectedDNAMethylationSitesCaseVsControlBoxplots(
#'     TENETMultiAssayExperiment = exampleTENETMultiAssayExperiment,
#'     DNAMethylationSites = c("cg03095778", "cg24011501", "cg12989041"),
#'     coreCount = 1
#' )
step7SelectedDNAMethylationSitesCaseVsControlBoxplots <- function(
  TENETMultiAssayExperiment,
  DNAMethylationSites,
  violinPlots = FALSE,
  coreCount = 1
) {
    ## Return an error message if the input MultiAssayExperiment is invalid
    .validateMultiAssayExperiment(TENETMultiAssayExperiment)

    if (missing(DNAMethylationSites)) {
        .stopNoCall("The DNAMethylationSites parameter must be specified.")
    }

    ## Get the methylation data
    methylationData <- .getExpOrMetSamplesOfType(
        TENETMultiAssayExperiment,
        "methylation"
    )

    ## Exclude RE DNA methylation sites that are not present in the methylation
    ## data
    methSiteList <- .excludeMissingMethylationSites(
        DNAMethylationSites, methylationData
    )

    ## Get the names and types of the methylation samples in the data
    sampleNames <- .getExpOrMetSamplesOfType(
        TENETMultiAssayExperiment,
        "methylation",
        namesOnly = TRUE
    )

    sampleTypes <- TENETMultiAssayExperiment@sampleMap[
        match(sampleNames, TENETMultiAssayExperiment@sampleMap$colname),
        "sampleType"
    ]

    ## Create a data frame containing group info for the case and control
    ## samples
    groupInfo <- data.frame(
        group = sampleNames,
        cluster = sampleTypes,
        stringsAsFactors = FALSE
    )

    ## Generate the plots for all RE DNA methylation sites of interest and add
    ## them to the MultiAssayExperiment
    TENETMultiAssayExperiment@metadata$
        step7SelectedDNAMethylationSitesCaseVsControlBoxplots <-
        parallel::mcmapply(
            FUN = .quadrantBoxplotFunction,
            geneOrMethSiteID = methSiteList,
            MoreArgs = list(
                expOrMet = "methylation",
                expOrMetData = methylationData,
                groupInfo = groupInfo,
                violinPlot = violinPlots
            ),
            mc.cores = coreCount,
            USE.NAMES = FALSE
        )

    return(TENETMultiAssayExperiment)
}
