## Internal functions used by step7ExpressionVsDNAMethylationScatterplots

## Internal function to generate a scatterplot for an RE DNA methylation
## site-gene combination
.internalScatterplotFunction <- function(
  geneOfInterest,
  methSiteOfInterest,
  expressionData,
  methylationData,
  sampleInfo,
  simpleOrComplex,
  geneIDNameDF,
  CNVDataset,
  SMDataset,
  purityValues,
  hyperHypo
) {
    ## Get gene expression values
    geneExpression <- c(unlist(
        expressionData[geneOfInterest, sampleInfo$expNames]
    ))

    ## Get DNA methylation values
    DNAMethylationSiteMethylation <- c(unlist(
        methylationData[methSiteOfInterest, sampleInfo$metNames]
    ))

    ## Convert the gene ID into the gene name
    geneName <- geneIDNameDF[geneOfInterest, "geneName"]

    ## Manually coloring samples
    groupColors <- c("Control" = "dodgerblue3", "Case" = "red3")

    ## Manually shaping samples
    mutationCopyNumberShape <- c(
        "Copy gain 1" = 0, # open square
        "Copy gain 2+" = 15, # colored square
        "Copy loss 1" = 2, # open upward triangle
        "Copy loss 2+" = 17, # colored upward triangle
        "No alteration" = 16, # colored circle
        "Mutation" = 8 # 8 pointed asterisk
    )

    ## Assemble a data frame of information for the boxplot
    boxplotDF <- data.frame(
        "geneExpression" = geneExpression,
        "DNAMethylationSiteMethylation" = DNAMethylationSiteMethylation,
        "primaryNames" = sampleInfo$primaryNames,
        "sampleType" = sampleInfo$type,
        "combinedSMCNVValue" = "No alteration",
        "purityLabels" = 1,
        stringsAsFactors = FALSE
    )

    ## Set the type of plot that will be done. If complex is selected but CNV
    ## or SM data aren't available for the gene, a simple plot will be done
    ## instead
    plotType <- "simple"

    if (simpleOrComplex == "complex") {
        ## Since data are available, do a complex plot later
        plotType <- "complex"

        ## Get the expected names for the genes in the CNV and SM data
        geneCNVColumnName <- paste0(geneOfInterest, "_CNV")
        geneSMColumnName <- paste0(geneOfInterest, "_SM")

        ## If data are present for the given gene in both the CNV and SM data,
        ## Get those values from their respective datasets. Otherwise
        ## use default values presented.
        if (
            (geneCNVColumnName %in% colnames(CNVDataset)) &
                (geneSMColumnName %in% colnames(SMDataset))
        ) {
            ## Add CNV status to the boxplotDF
            boxplotDF$rawCNVValues <- CNVDataset[
                boxplotDF$primaryNames,
                paste0(geneOfInterest, "_CNV")
            ]

            ## Transform CNV values for control samples to 0
            ## Again, we are supressing warnings since we want to
            ## use as.numeric to coerce non-numeric values to NAs
            boxplotDF$CNVValuesControl0 <- suppressWarnings(
                as.numeric(
                    ifelse(
                        boxplotDF$sampleType == "Control",
                        0,
                        boxplotDF$rawCNVValues
                    )
                )
            )

            ## Create grouped labels for CNV
            boxplotDF$CNVLabels <- ifelse(
                boxplotDF$CNVValuesControl0 < -1,
                "Copy loss 2+",
                ifelse(
                    boxplotDF$CNVValuesControl0 < 0,
                    "Copy loss 1",
                    ifelse(
                        boxplotDF$CNVValuesControl0 == 0,
                        "No alteration",
                        ifelse(
                            boxplotDF$CNVValuesControl0 > 1,
                            "Copy gain 2+",
                            ifelse(
                                boxplotDF$CNVValuesControl0 > 0,
                                "Copy gain 1",
                                "No alteration"
                            )
                        )
                    )
                )
            )

            ## Add SM status to the boxplotDF
            boxplotDF$rawSMValues <- SMDataset[
                boxplotDF$primaryNames,
                paste0(geneOfInterest, "_SM")
            ]

            ## Transform SM values for control samples to 0
            boxplotDF$SMValuesControl0 <- ifelse(
                boxplotDF$sampleType == "Control",
                0,
                boxplotDF$rawSMValues
            )

            ## Create grouped labels for SM

            ## Create a small conversion table
            SMConversionDF <- data.frame(
                "values" = c("0", "no mutation", "1", "mutation"),
                "return" = c(
                    "No alteration", "No alteration", "Mutation", "Mutation"
                ),
                stringsAsFactors = FALSE
            )
            rownames(SMConversionDF) <- SMConversionDF$values

            ## Use the conversion table
            boxplotDF$SMLabels <- ifelse(
                tolower(as.character(boxplotDF$SMValuesControl0)) %in%
                    SMConversionDF$values,
                SMConversionDF[
                    tolower(as.character(boxplotDF$SMValuesControl0)),
                    "return"
                ],
                "No alteration"
            )

            ## Combine the SM and CNV labels
            ## Give precedence to mutations over CNV
            boxplotDF$combinedSMCNVValue <- ifelse(
                boxplotDF$SMLabels == "Mutation",
                boxplotDF$SMLabels,
                boxplotDF$CNVLabels
            )
        }

        ## Add the raw purity values and convert out of bounds values to NA
        boxplotDF$rawPurity <- ifelse(
            purityValues[boxplotDF$primaryNames, 1] > 1,
            NA,
            ifelse(
                purityValues[boxplotDF$primaryNames, 1] < 0,
                NA,
                purityValues[boxplotDF$primaryNames, 1]
            )
        )

        ## Add values for proper sizing later
        boxplotDF$purityScaled <- (1 + (2 * boxplotDF$rawPurity))

        ## Convert NA values to 1s
        boxplotDF$purityScaledNA1 <- ifelse(
            is.na(boxplotDF$purityScaled),
            1,
            boxplotDF$purityScaled
        )

        ## Convert control samples to 3s even if they are NA
        boxplotDF$purityLabels <- ifelse(
            boxplotDF$sampleType == "Control",
            3,
            boxplotDF$purityScaledNA1
        )
    }

    if (plotType == "simple") {
        ## Create a simple scatterplot which lacks changes to point size
        ## and shape
        scatterplot <- ggplot2::ggplot(
            boxplotDF,
            ggplot2::aes(
                x = get("geneExpression"),
                y = get("DNAMethylationSiteMethylation"),
                color = get("sampleType")
            )
        ) +
            ggplot2::geom_point()
    } else if (plotType == "complex") {
        ## Create a complex scatterplot which will adjust shape based
        ## on the combined CNV + SM status and size based on sample purity
        scatterplot <- ggplot2::ggplot(
            boxplotDF,
            ggplot2::aes(
                x = get("geneExpression"),
                y = get("DNAMethylationSiteMethylation"),
                color = get("sampleType"),
                shape = get("combinedSMCNVValue")
            )
        ) +
            ggplot2::geom_point(ggplot2::aes(size = get("purityLabels"))) +
            ggplot2::scale_shape_manual(
                values = mutationCopyNumberShape,
                name = "Copy number variation or mutation",
                labels = sort(unique(boxplotDF$combinedSMCNVValue))
            ) +
            ggplot2::scale_size_continuous(
                name = "Sample purity",
                breaks = c(
                    min(boxplotDF$purityLabels),
                    2,
                    max(boxplotDF$purityLabels)
                ),
                labels = c("Low", "Medium", "High (or control sample)")
            ) +
            ggplot2::guides(
                color = ggplot2::guide_legend(order = 1),
                shape = ggplot2::guide_legend(order = 2),
                size = ggplot2::guide_legend(order = 3)
            )
    }

    ## Add more information to the scatterplot, regardless of its type
    scatterplot <- scatterplot + ggplot2::ylab(
        paste(methSiteOfInterest, "DNA methylation")
    ) + ggplot2::xlab(
        paste(geneName, "-", geneOfInterest, "expression")
    ) + ggplot2::theme_bw() +
        ggplot2::scale_color_manual(
            values = groupColors,
            name = "Sample type"
        ) + ggplot2::theme(
            plot.title = ggplot2::element_text(hjust = 0.5, size = 18),
            legend.title = ggplot2::element_text(hjust = 0.5, size = 12),
            legend.text = ggplot2::element_text(size = 10, color = "black"),
            panel.border = ggplot2::element_rect(
                color = "black", fill = NA, linewidth = 1
            ),
            plot.background = ggplot2::element_rect(fill = "white"),
            axis.title.x = ggplot2::element_text(size = 16),
            axis.title.y = ggplot2::element_text(size = 16),
            axis.text.x = ggplot2::element_text(size = 14, color = "black"),
            axis.text.y = ggplot2::element_text(size = 14, color = "black"),
            panel.grid.major = ggplot2::element_blank(),
            panel.grid.minor = ggplot2::element_blank()
        )

    ## To display this plot correctly, the plot area should be at least 10
    ## inches wide by 7 inches tall
    return(TENETSavedSizePlot(scatterplot, width = 10, height = 7))
}

## Internal function that, when given a list of RE DNA methylation sites linked
## to a gene, or genes linked to an RE DNA methylation site, will plot
## scatterplots showing the methylation of the RE DNA methylation site on the Y
## axis and expression of that gene on the X axis across all the case and
## control samples. listValues must be either "genes" or "methSites".
.quadrantScatterplotFunction <- function(
  geneOrMethSiteListOfInterest,
  geneOrMethSiteVectorOfInterest,
  listValues,
  geneIDNameDF,
  expressionData,
  methylationData,
  sampleInfo,
  metToExpSampleConversion,
  simpleOrComplex,
  CNVDataset,
  SMDataset,
  purityValues,
  hyperHypo
) {
    ## Unlist the list values
    unlistedValues <- c(unlist(geneOrMethSiteListOfInterest))

    ## Use the converter to get the names of the expression samples that
    ## correspond with the identified methylation samples
    matchedMetNames <- colnames(methylationData)
    matchedExpNames <- unname(metToExpSampleConversion[matchedMetNames])

    ## Depending on whether the listValues is "methSites" or "genes", execute
    ## the plotting function with different layering. I.e., if doing the normal
    ## analysis, listValues will be "methSites" each time the function is run,
    ## a single gene will be provided, with a list of RE DNA methylation sites
    ## linked to it. Otherwise, if listValues is "genes", when dealing with the
    ## RE DNA methylation sites given to the 'DNAMethylationSites' argument of
    ## the top level function, instead a single RE DNA methylation site will be
    ## provided, with a list of genes linked to it.
    if (listValues == "methSites") {
        result <- lapply(
            unlistedValues, # These will be DNA methylation sites
            .internalScatterplotFunction,
            geneOfInterest = geneOrMethSiteVectorOfInterest,
            expressionData = expressionData,
            methylationData = methylationData,
            sampleInfo = sampleInfo,
            simpleOrComplex = simpleOrComplex,
            geneIDNameDF = geneIDNameDF,
            CNVDataset = CNVDataset,
            SMDataset = SMDataset,
            purityValues = purityValues,
            hyperHypo = hyperHypo
        )
    } else {
        result <- lapply(
            unlistedValues, # These will be gene IDs
            .internalScatterplotFunction,
            methSiteOfInterest = geneOrMethSiteVectorOfInterest,
            expressionData = expressionData,
            methylationData = methylationData,
            sampleInfo = sampleInfo,
            simpleOrComplex = simpleOrComplex,
            geneIDNameDF = geneIDNameDF,
            CNVDataset = CNVDataset,
            SMDataset = SMDataset,
            purityValues = purityValues,
            hyperHypo = hyperHypo
        )
    }

    ## Name each returned scatterplot after the listed values, then
    ## return the list of created scatterplots
    names(result) <- unname(unlistedValues)

    ## If the results are empty, instead display a message noting that no linked
    ## RE DNA methylation sites/genes were identified
    if (length(result) == 0) {
        thisEndOfLink <- ifelse(
            listValues != "methSites", "gene", "RE DNA methylation site"
        )
        otherEndOfLink <- ifelse(
            listValues == "methSites", "gene", "RE DNA methylation site"
        )
        return(paste0(
            "No ", thisEndOfLink, "s linked to the given ", otherEndOfLink,
            " were identified, so no scatterplots have been created for this ",
            otherEndOfLink, "."
        ))
    } else {
        return(result)
    }
}

## Main step7ExpressionVsDNAMethylationScatterplots function

#' Create scatterplots displaying the expression of the top genes and the
#' methylation levels of each of their linked RE DNA methylation sites,
#' optionally incorporating copy number variation, somatic mutation, and purity
#' data
#'
#' This function takes the top genes and transcription factors by number of
#' linked RE DNA methylation sites identified by the
#' `step6DNAMethylationSitesPerGeneTabulation` function up to a number
#' specified by the user, or all genes linked to selected RE DNA methylation
#' sites specified by the user, and generates scatterplots displaying
#' the expression level of each of these genes in the X-axis and the
#' methylation level of each RE DNA methylation site linked to them in the
#' Y-axis for the hyper- and/or hypomethylated G+ analysis quadrants.
#' The scatterplots may optionally incorporate provided copy number variation
#' (CNV), somatic mutation (SM), and purity information for each sample.
#'
#' @param TENETMultiAssayExperiment Specify a MultiAssayExperiment object
#' containing expression and methylation SummarizedExperiment objects, such as
#' one created by the TCGADownloader function. The object's metadata
#' must contain the results from the `step5OptimizeLinks` and
#' `step6DNAMethylationSitesPerGeneTabulation` functions.
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
#' @param hypermethGplusAnalysis Set to TRUE to create scatterplots for genes
#' with hypermethylated RE DNA methylation sites with G+ links and each of their
#' linked RE DNA methylation sites. Defaults to TRUE.
#' @param hypomethGplusAnalysis Set to TRUE to create scatterplots for genes
#' with hypomethylated RE DNA methylation sites with G+ links and each of their
#' linked RE DNA methylation sites. Defaults to TRUE.
#' @param topGeneNumber Specify the number of top genes and TFs, based on the
#' most linked RE DNA methylation sites of a given analysis type, for which to
#' create scatterplots. Defaults to 10.
#' @param DNAMethylationSites Supply a vector of RE DNA methylation site IDs for
#' which scatterplots will be generated, if these sites have any linked
#' genes/TFs with expression in each specified analysis type.
#' @param simpleOrComplex Set to 'complex' to incorporate copy number variation,
#' somatic mutation, and purity data into the scatterplots. Otherwise, set to
#' 'simple'. If set to 'complex', copy number variation, somatic mutation, and
#' purity data must be provided via the `CNVData`, `SMData`, and `purityData`
#' arguments respectively. **Note:** At this time, either all or none of these
#' optional data types must be provided. Defaults to 'simple'.
#' @param CNVData Specify a dataset containing CNV status for each of the top
#' genes, as selected by the analysis type and 'topGeneNumber' arguments, in
#' each sample in the TENETMultiAssayExperiment. CNV status must be an
#' integer representing the change in copy number for each gene, with negative
#' numbers representing a loss and positive numbers representing a gain.
#' **Note:** Copy number changes of 2 or more will be grouped together. The
#' dataset may be given as a data frame, matrix, or TSV file path. If it is a
#' data frame or matrix, its rownames must contain sample names. If a TSV file
#' is provided, the first column must contain sample names, and the first row
#' must contain column headers. Sample names must match those in the colData of
#' the TENETMultiAssayExperiment object. Column names must contain gene IDs
#' followed by "_CNV". If set to NA, the data will be loaded from the colData of
#' the TENETMultiAssayExperiment object. **Note:** If data are missing for a
#' given gene, the plot will be generated without considering its CNV status.
#' Defaults to NA, and is only considered if `simpleOrComplex` is set to
#' "complex".
#' @param SMData Specify a dataset containing the somatic mutation status for
#' each of the top genes in each sample in the TENETMultiAssayExperiment. This
#' argument behaves the same way as the `CNVData` argument, except that the
#' names of the columns containing SM status must end with "_SM", and the status
#' must be an integer 0 or 1 or a string "no mutation" or "mutation". Defaults
#' to NA.
#' @param purityData Specify the cellularity/purity data for each sample in the
#' TENETMultiAssayExperiment. Purity values must range from 0 to 1. The dataset
#' may be given as a vector, data frame, matrix, or TSV file path. If a vector
#' is given, the names of the vector elements must correspond to the names of
#' the samples in the rownames of the colData of the TENETMultiAssayExperiment
#' object. If no names are provided for the vector, then the number of elements
#' in the vector must equal the number of samples in the colData, and it is
#' assumed to align with the samples as they are ordered in the colData. If a
#' data frame, matrix, or TSV file is given, it must be in the same format as
#' for the `CNVData` argument, except that the first column of data (excluding
#' the rownames) must contain the purity data. If this argument is set to NA,
#' purity data will be loaded from the "purity" column of the colData of the
#' TENETMultiAssayExperiment object. Defaults to NA, and is only considered if
#' 'simpleOrComplex' is set to "complex".
#' @param coreCount Argument passed as the mc.cores argument to mcmapply. See
#' `?parallel::mcmapply` for more details. Defaults to 1.
#' @return Returns the MultiAssayExperiment object given as the
#' TENETMultiAssayExperiment argument with an additional list
#' named 'step7ExpressionVsDNAMethylationScatterplots' in its metadata with the
#' output of this function. This list is subdivided into hypermethGplus or
#' hypomethGplus results as selected by the user, which are further subdivided
#' into lists with data for the top overall genes, and for top TF genes only.
#' Each of these lists contains a final list for each of the top genes/TFs
#' containing scatterplots for each RE DNA methylation site linked to the gene.
#' If the user has specified RE DNA methylation sites of interest, an additional
#' list named 'selectedDNAMethylationSites' is generated for each quadrant
#' containing scatterplots for each gene linked to each specified RE DNA
#' methylation site. In each scatterplot, the expression of the gene is plotted
#' on the X-axis, and the methylation of the linked RE DNA methylation site is
#' plotted on the Y-axis. If complex plots are being created, the CNV and SM
#' status of each sample, if present, will be represented by each point's shape
#' (with SM status taking precedence over CNV), and the purity of each sample
#' will be reflected in each point's size.
#'
#' @export
#'
#' @examplesIf interactive()
#' ## This example uses the example MultiAssayExperiment provided in the
#' ## TENET.ExperimentHub package to create scatterplots for the top 10
#' ## genes and TFs by number of linked hyper- and hypomethylated RE DNA
#' ## methylation sites, showing expression of these genes and the DNA
#' ## methylation level of their linked RE DNA methylation sites. Gene names
#' ## will be retrieved from the rowRanges of the 'expression'
#' ## SummarizedExperiment object in the example MultiAssayExperiment. No CNV,
#' ## SM, or purity data will be incorporated, and the analysis will be
#' ## performed using one CPU core.
#'
#' ## Load the example TENET MultiAssayExperiment object
#' ## from the TENET.ExperimentHub package
#' exampleTENETMultiAssayExperiment <-
#'     TENET.ExperimentHub::exampleTENETMultiAssayExperiment()
#'
#' ## Use the example dataset to create the scatterplots
#' returnValue <- step7ExpressionVsDNAMethylationScatterplots(
#'     TENETMultiAssayExperiment = exampleTENETMultiAssayExperiment
#' )
#'
#' ## This example demonstrates many of the analysis options, creating
#' ## scatterplots for the top 5 genes/TFs as well as some example RE DNA
#' ## methylation sites of interest. As before, gene names will be retrieved
#' ## from the rowRanges of the 'expression' SummarizedExperiment object.
#' ## Complex scatterplots are created which display each sample's CNV and SM
#' ## status for each gene, as well as purity data, where available. The CNV,
#' ## SM, and purity data will be taken from specific columns of the
#' ## exampleTENETClinicalDataFrame object. The analysis will be performed using
#' ## 8 CPU cores.
#'
#' ## Load the example TENET MultiAssayExperiment object
#' ## from the TENET.ExperimentHub package
#' exampleTENETMultiAssayExperiment <-
#'     TENET.ExperimentHub::exampleTENETMultiAssayExperiment()
#'
#' ## Load the data frame with example clinical data for patients in the TENET
#' ## MultiAssayExperiment object from the TENET.ExperimentHub package
#' exampleTENETClinicalDataFrame <-
#'     TENET.ExperimentHub::exampleTENETClinicalDataFrame()
#'
#' ## Use the example datasets to create the scatterplots
#' returnValue <- step7ExpressionVsDNAMethylationScatterplots(
#'     TENETMultiAssayExperiment = exampleTENETMultiAssayExperiment,
#'     hypermethGplusAnalysis = FALSE,
#'     topGeneNumber = 5,
#'     DNAMethylationSites = c("cg03095778", "cg24011501", "cg12989041"),
#'     simpleOrComplex = "complex",
#'     CNVData = exampleTENETClinicalDataFrame[seq(4, 42, by = 2)],
#'     SMData = exampleTENETClinicalDataFrame[seq(5, 43, by = 2)],
#'     purityData = exampleTENETClinicalDataFrame[3],
#'     coreCount = 8
#' )
step7ExpressionVsDNAMethylationScatterplots <- function(
  TENETMultiAssayExperiment,
  geneAnnotationDataset = NA,
  hypermethGplusAnalysis = TRUE,
  hypomethGplusAnalysis = TRUE,
  topGeneNumber = 10,
  DNAMethylationSites = NA,
  simpleOrComplex = "simple",
  CNVData = NA,
  SMData = NA,
  purityData = NA,
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

    ## Get gene IDs and names from the MAE, or gene annotation dataset if
    ## provided
    geneIDNameDF <- .getGeneIDsAndNames(
        TENETMultiAssayExperiment, geneAnnotationDataset
    )

    ## Get the expression and methylation data for all samples
    expressionData <- .getExpOrMetSamplesOfType(
        TENETMultiAssayExperiment,
        "expression"
    )

    methylationData <- .getExpOrMetSamplesOfType(
        TENETMultiAssayExperiment,
        "methylation"
    )

    ## Get the names of the samples in the expression and methylation data
    expressionSampleNames <- .getExpOrMetSamplesOfType(
        TENETMultiAssayExperiment,
        "expression",
        namesOnly = TRUE
    )

    methylationSampleNames <- .getExpOrMetSamplesOfType(
        TENETMultiAssayExperiment,
        "methylation",
        namesOnly = TRUE
    )

    if (!.isSingleNA(DNAMethylationSites)) {
        ## Exclude RE DNA methylation sites that are not present in the
        ## methylation data
        DNAMethylationSites <- .excludeMissingMethylationSites(
            DNAMethylationSites, methylationData
        )

        ## If there are no more eligible DNA methylation sites after they have
        ## been excluded, return an error message
        if (length(DNAMethylationSites) == 0) {
            .stopNoCall(
                "All specified RE DNA methylation sites are missing from the ",
                "methylation data in the specified TENETMultiAssayExperiment. ",
                "Please check the specified RE DNA methylation sites and ",
                "rerun this function."
            )
        }
    }

    ## Get the methylation values that match with expression values
    ## using the mapping data. This assumes the methylation and expression
    ## values share a clinical data match within the mapping.
    metToExpSampleConversion <- .createMetToExpSampleConversionVector(
        TENETMultiAssayExperiment
    )

    ## Create a data frame containing collapsed info from the sample mapping.
    ## This is used to pair expression and methylation data, and obtain sample
    ## purity, CNV, and SM data, if wanted
    sampleInfo <- data.frame(
        primaryNames = TENETMultiAssayExperiment@sampleMap[
            match(
                methylationSampleNames,
                TENETMultiAssayExperiment@sampleMap$colname
            ),
            "primary"
        ],
        expNames = unname(metToExpSampleConversion[methylationSampleNames]),
        metNames = methylationSampleNames,
        type = TENETMultiAssayExperiment@sampleMap[
            match(
                methylationSampleNames,
                TENETMultiAssayExperiment@sampleMap$colname
            ),
            "sampleType"
        ],
        stringsAsFactors = FALSE
    )

    ## Create an empty list to hold the results from this step 7 function
    resultsList <- list()

    ## Generate results for the selected analysis types
    for (hyperHypo in analysisTypes) {
        quadrantResultsName <- paste0(hyperHypo, "methGplusResults")

        resultsList[[quadrantResultsName]] <- list()

        ## Ensure the quadrant's results are present in step 5
        .ensureStepPresent(
            TENETMultiAssayExperiment,
            stepName = "step5OptimizeLinks",
            substepName = quadrantResultsName
        )

        ## Load the quadrant's significant links from step 5
        quadrantSigLinkZScores <- TENETMultiAssayExperiment@metadata$
            step5OptimizeLinks[[quadrantResultsName]]

        for (geneOrTF in c("Gene", "TF")) {
            ## Get the IDs of the top genes/TFs in this quadrant
            quadrantGenesOfInterest <- .getQuadrantTopGenesOrTFs(
                TENETMultiAssayExperiment, geneOrTF, hyperHypo,
                topGeneNumber
            )$geneID

            resultsMainName <- "step7ExpressionVsDNAMethylationScatterplots"
            resultsSubName <- paste0("top", geneOrTF, "s")

            if (.isSingleNA(quadrantGenesOfInterest)) {
                resultsList[[quadrantResultsName]][[resultsSubName]] <- NA
                next
            }

            if (simpleOrComplex == "complex") {
                ## Process the copy number variation data of the samples if
                ## CNVData is not NA. Otherwise, look into the supplied colData
                ## of the TENETMultiAssayExperiment to see if CNV data are
                ## available for each of the top genes/TFs.
                CNVDataResults <- .importClinicalData(
                    userInput = CNVData,
                    argumentName = "CNVData",
                    clinicalDataColumn = paste0(
                        quadrantGenesOfInterest, "_CNV"
                    ),
                    returnType = "multiple",
                    TENETMultiAssayExperiment = TENETMultiAssayExperiment
                )

                ## Process the somatic mutation data of the samples if SMData
                ## is not NA. Otherwise, look into the supplied colData of the
                ## TENETMultiAssayExperiment to see if SM data are available
                ## for each of the top genes/TFs.
                SMDataResults <- .importClinicalData(
                    userInput = SMData,
                    argumentName = "SMData",
                    clinicalDataColumn = paste0(
                        quadrantGenesOfInterest, "_SM"
                    ),
                    returnType = "multiple",
                    TENETMultiAssayExperiment = TENETMultiAssayExperiment
                )

                ## Process the purity data of the samples if
                ## purityData is not NA. Otherwise, look into the supplied
                ## colData of the TENETMultiAssayExperiment to see if purity
                ## data are available for each of the top genes/TFs. Unlike CNV
                ## and SM, purity data should be a single column.
                purityDataResults <- .importClinicalData(
                    userInput = purityData,
                    argumentName = "purityData",
                    clinicalDataColumn = "purity",
                    returnType = "single",
                    TENETMultiAssayExperiment = TENETMultiAssayExperiment
                )
            } else {
                CNVDataResults <- NA
                SMDataResults <- NA
                purityDataResults <- NA
            }

            ## Create an empty list for the RE DNA methylation sites linked to
            ## the genes of interest
            quadrantMethSitesOfInterest <- list()

            ## For each gene of interest, get the list of RE DNA methylation
            ## sites associated with it
            for (i in seq_along(quadrantGenesOfInterest)) {
                ## Get the RE DNA methylation sites linked to the gene and add
                ## the RE DNA methylation sites to the list
                quadrantMethSitesOfInterest[[i]] <- unique(
                    quadrantSigLinkZScores[
                        quadrantSigLinkZScores$geneID ==
                            quadrantGenesOfInterest[i],
                        "DNAMethylationSiteID"
                    ]
                )
            }

            ## Add the names of the genes to the lists
            names(quadrantMethSitesOfInterest) <- quadrantGenesOfInterest

            ## Generate the plots for the genes of interest and add them to the
            ## results list
            resultsList[[quadrantResultsName]][[resultsSubName]] <-
                parallel::mcmapply(
                    FUN = .quadrantScatterplotFunction,
                    geneOrMethSiteListOfInterest = quadrantMethSitesOfInterest,
                    geneOrMethSiteVectorOfInterest = names(
                        quadrantMethSitesOfInterest
                    ),
                    MoreArgs = list(
                        "listValues" = "methSites",
                        "geneIDNameDF" = geneIDNameDF,
                        "expressionData" = expressionData,
                        "methylationData" = methylationData,
                        "sampleInfo" = sampleInfo,
                        "metToExpSampleConversion" =
                            metToExpSampleConversion,
                        "simpleOrComplex" = simpleOrComplex,
                        "CNVDataset" = CNVDataResults,
                        "SMDataset" = SMDataResults,
                        "purityValues" = purityDataResults,
                        "hyperHypo" = hyperHypo
                    ),
                    mc.cores = coreCount
                )
        }

        ## If selected RE DNA methylation sites were used, add a separate list
        ## to the quadrant results, which will include the specified RE DNA
        ## methylation sites instead of genes, with individual plots to each
        ## gene under them with the names of any genes linked to the RE DNA
        ## methylation sites
        if (!.isSingleNA(DNAMethylationSites)) {
            ## Create an empty list for the genes linked to the RE DNA
            ## methylation sites of interest
            quadrantGenesOfInterest <- list()

            ## Get the genes linked to each RE DNA methylation site and add them
            ## to the list
            for (i in seq_along(DNAMethylationSites)) {
                quadrantGenesOfInterest[[i]] <- unique(
                    quadrantSigLinkZScores[
                        quadrantSigLinkZScores$DNAMethylationSiteID ==
                            DNAMethylationSites[i],
                        "geneID"
                    ]
                )
            }

            ## Add the names of the genes to the lists
            names(quadrantGenesOfInterest) <- DNAMethylationSites
            resultsSubName <- "selectedDNAMethylationSites"

            ## Generate the plots for the RE DNA methylation sites of interest
            ## and add them to the results list
            resultsList[[quadrantResultsName]][[resultsSubName]] <-
                parallel::mcmapply(
                    FUN = .quadrantScatterplotFunction,
                    geneOrMethSiteListOfInterest = quadrantGenesOfInterest,
                    geneOrMethSiteVectorOfInterest = names(
                        quadrantGenesOfInterest
                    ),
                    MoreArgs = list(
                        "listValues" = "genes",
                        "geneIDNameDF" = geneIDNameDF,
                        "expressionData" = expressionData,
                        "methylationData" = methylationData,
                        "sampleInfo" = sampleInfo,
                        "metToExpSampleConversion" =
                            metToExpSampleConversion,
                        "simpleOrComplex" = simpleOrComplex,
                        "CNVDataset" = CNVDataResults,
                        "SMDataset" = SMDataResults,
                        "purityValues" = purityDataResults,
                        "hyperHypo" = hyperHypo
                    ),
                    mc.cores = coreCount
                )
        }
    }

    ## Add the results list to the MultiAssayExperiment
    TENETMultiAssayExperiment@metadata[[resultsMainName]] <- resultsList

    return(TENETMultiAssayExperiment)
}
