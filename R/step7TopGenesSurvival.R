## Internal functions used by step7TopGenesSurvival

## Internal function to generate KM survival statistics
## for gene expression of genes or DNA methylation of sites
## and create plots if specified
.survivalFunction <- function(
    inputID, ## The gene or RE DNA methylation site's ID
    expressionOrMethylation, ## "Expression" or "Methylation"
    geneIDdf = NULL, ## Specify when expressionOrMethylation is "Expression"
    clinicalObject,
    TENETMultiAssayExperiment,
    survivalGroupingCutoffs,
    useJenksBreaks,
    jenksBreaksGroupCount,
    createPlot ## TRUE or FALSE - affects plots for KM only
    ) {
    ## If a gene was specified, get the gene name corresponding to the gene ID
    if (expressionOrMethylation == "Expression") {
        inputName <- geneIDdf[
            inputID,
            "geneName"
        ]
    } else if (expressionOrMethylation == "Methylation") {
        ## The RE DNA methylation site's name is its ID
        inputName <- inputID
    }

    ## Get the expression/methylation of the gene/RE DNA methylation site and
    ## add it to the clinical object
    if (expressionOrMethylation == "Expression") {
        clinicalObject$inputValue <- TENETMultiAssayExperiment@
        ExperimentList$expression@assays@data$expression[
            inputID,
            clinicalObject$expressionSampleNames
        ]
    } else if (expressionOrMethylation == "Methylation") {
        clinicalObject$inputValue <- TENETMultiAssayExperiment@
        ExperimentList$methylation@assays@data$methylation[
            inputID,
            clinicalObject$methylationSampleNames
        ]
    }

    ## Calculate some basic data
    controlSampleN <- nrow(
        clinicalObject[clinicalObject$sampleType == "Control", ]
    )

    caseSampleN <- nrow(
        clinicalObject[clinicalObject$sampleType == "Case", ]
    )

    ## Count the number samples with NA expression/methylation
    NACountControlInputValue <- sum(
        is.na(
            clinicalObject[
                clinicalObject$sampleType == "Control",
                "inputValue"
            ]
        )
    )

    NACountCaseInputValue <- sum(
        is.na(
            clinicalObject[
                clinicalObject$sampleType == "Case",
                "inputValue"
            ]
        )
    )

    ## Calculate mean expression/methylation for all samples
    ## that have expression/methylation
    controlMeanInputValue <- mean(
        clinicalObject[
            clinicalObject$sampleType == "Control",
            "inputValue"
        ],
        na.rm = TRUE
    )

    caseMeanInputValue <- mean(
        clinicalObject[
            clinicalObject$sampleType == "Case",
            "inputValue"
        ],
        na.rm = TRUE
    )

    ## Calculate the number of case samples remaining with missing
    ## clinical data
    NACountCaseClinical <- sum(
        !stats::complete.cases(
            clinicalObject[
                clinicalObject$sampleType == "Case",
                c("vitalStatus", "time")
            ]
        )
    )

    ## Calculate the number of case/control samples considered
    ## Case samples considered have complete expression/methylation + survival
    ## clinical data.
    ## Control samples only need expression/methylation (they aren't included
    ## in the survival analyses, only for expression/methylation statistics)
    controlPresentSampleN <- sum(
        !is.na(
            clinicalObject[
                clinicalObject$sampleType == "Control",
                "inputValue"
            ]
        )
    )

    casePresentSampleN <- sum(
        stats::complete.cases(
            clinicalObject[
                clinicalObject$sampleType == "Case",
                c("vitalStatus", "time", "inputValue")
            ]
        )
    )

    ## Create a subsetted dataset with only the case samples
    ## that have complete information for the vitalStatus, time, and
    ## expression/methylation variables
    completeCasesClinicalObject <- clinicalObject[
        clinicalObject$sampleType == "Case",
    ]

    completeCasesClinicalObject <- completeCasesClinicalObject[
        stats::complete.cases(
            completeCasesClinicalObject[
                , c("vitalStatus", "time", "inputValue")
            ]
        ),
    ]

    ### Now determine the groups of interest:

    ## First, if Jenks breaks are TRUE, we need to calculate where they are:
    if (useJenksBreaks) {
        ## First, add one to the jenksBreaksGroupCount. This
        ## is because when breaks are calculated, it actually generates that
        ## many breaks, including lowest and highest bounds, so it actually
        ## creates one fewer groups:
        jenksBreaksGroupCountInt <- (jenksBreaksGroupCount + 1)

        ## Calculate the Jenks breaks for the number of groups specified by the
        ## user
        breaksValues <- BAMMtools::getJenksBreaks(
            completeCasesClinicalObject[, "inputValue"],
            k = jenksBreaksGroupCountInt
        )

        ## Sometimes breaks aren't unique - if this is the case, make the breaks
        ## unique by adding a small offset value to the breaks that aren't
        ## unique:
        while (any(duplicated(breaksValues))) {
            breaksValues[
                duplicated(breaksValues)
            ] <- breaksValues[
                duplicated(breaksValues)
            ] + (max(breaksValues) * 0.01)
        }

        ## Create a vector with each of the groups' cutoff values:
        cutoffVector <- c(
            rbind(breaksValues[-length(breaksValues)], breaksValues[-1])
        )

        ## Create names for the Jenks Groups:
        if (jenksBreaksGroupCount == 2) {
            GroupNames <- c("jenksGroupA_(Lowest)", "jenksGroupB_(Highest)")
        } else {
            ## If there are more than two groups, then there are groups in the
            ## middle that don't get a special label:
            GroupNames <- c(
                "jenksGroupA_(Lowest)",
                paste0(
                    "jenksGroup",
                    LETTERS[seq_len(jenksBreaksGroupCount)[
                        c(-1, -jenksBreaksGroupCount)
                    ]]
                ),
                paste0(
                    "jenksGroup",
                    LETTERS[jenksBreaksGroupCount],
                    "_(Highest)"
                )
            )
        }

        ## Use the cut() function to assign groups based on the Jenks breaks,
        ## then adjust the names of the groups to be the names we set above.
        ## This will ensure the first group includes the lowest sample, up to
        ## and including samples at the first break point, and every group after
        ## will exclude samples at the lower bound, but include those at the
        ## next breakpoint.
        jenksGroupsValuesFactor <- factor(
            cut(
                completeCasesClinicalObject[, "inputValue"],
                breaks = breaksValues,
                include.lowest = TRUE,
                right = TRUE
            )
        )
        levels(jenksGroupsValuesFactor) <- GroupNames

        ## Then assign these values back to the completeCasesClinicalObject
        completeCasesClinicalObject$grouping <- jenksGroupsValuesFactor
    } else {
        ## Use the row names of the survivalGroupingCutoffs as the group names
        GroupNames <- rownames(survivalGroupingCutoffs)

        ## For each proportion value listed in the columns of
        ## survivalGroupingCutoffs, get an actual value correlating to it
        survivalGroupingCutoffs[, 3] <- unname(
            stats::quantile(
                completeCasesClinicalObject[, "inputValue"],
                survivalGroupingCutoffs[, 1],
                na.rm = TRUE
            )
        )

        survivalGroupingCutoffs[, 4] <- unname(
            stats::quantile(
                completeCasesClinicalObject[, "inputValue"],
                survivalGroupingCutoffs[, 2],
                na.rm = TRUE
            )
        )

        ## Create a vector with each of the groups' cutoff values:
        cutoffVector <- c(
            rbind(survivalGroupingCutoffs[, 3], survivalGroupingCutoffs[, 4])
        )

        ## For each of the inputValues, see what survival grouping it falls in
        setGroupsValuesFactor <- NULL

        for (i in seq_len(nrow(completeCasesClinicalObject))) {
            ## Get the value for i:
            iValue <- completeCasesClinicalObject[i, "inputValue"]

            ## Find the row, if any, the value is between:
            rowOverlapBool <- NULL

            for (j in seq_len(nrow(survivalGroupingCutoffs))) {
                if (j == 1) {
                    rowOverlapBool <- c(
                        rowOverlapBool,
                        (
                            iValue >= survivalGroupingCutoffs[j, 3] &
                                iValue <= survivalGroupingCutoffs[j, 4]
                        )
                    )
                } else {
                    rowOverlapBool <- c(
                        rowOverlapBool,
                        (
                            iValue > survivalGroupingCutoffs[j, 3] &
                                iValue <= survivalGroupingCutoffs[j, 4]
                        )
                    )
                }
            }

            ## Use the row overlap Bool to identify which group the sample is
            ## in
            if (!any(rowOverlapBool)) {
                setGroupsValuesFactor <- c(
                    setGroupsValuesFactor,
                    NA
                )
            } else {
                setGroupsValuesFactor <- c(
                    setGroupsValuesFactor,
                    which(rowOverlapBool)
                )
            }
        }

        ## Now set the levels of the factor to be the GroupNames
        setGroupsValuesFactor <- factor(setGroupsValuesFactor)
        levels(setGroupsValuesFactor) <- GroupNames

        ## Then assign these values back to the completeCasesClinicalObject
        completeCasesClinicalObject$grouping <- setGroupsValuesFactor
    }

    ## Now get the counts for each group in the dataset:

    ## First get a count of all the NA values, if there are any
    NACountGroup <- sum(is.na(completeCasesClinicalObject$grouping))

    ## Next get the counts from the "Freq" column when converting a table of the
    ## grouping column to a data frame (this will automatically ignore all NA
    ## values which is why we needed to get an NA count first)
    groupCounts <- as.data.frame(table(completeCasesClinicalObject$grouping))$
        Freq

    ## Next calculate a mean expression/methylation value for each group in the
    ## dataset:

    ## First do it for samples that are NA:
    NAGroupMean <- mean(
        completeCasesClinicalObject[
            is.na(completeCasesClinicalObject$grouping),
            "inputValue"
        ]
    )

    ## Then get the means for the groups in order
    groupMeans <- stats::aggregate(
        . ~ grouping,
        completeCasesClinicalObject[, c("inputValue", "grouping")],
        mean
    )$inputValue

    ## Calculate the proportion of samples which have reached the event
    ## rather than not/censored in the groupings
    ## This can be done largely with aggregate means on vitalStatus, subtracting
    ## 1 since events are listed as 2, while not/censored as 1.
    NAProportionEvent <- mean(
        completeCasesClinicalObject[
            is.na(completeCasesClinicalObject$grouping),
            "vitalStatus"
        ]
    ) - 1

    groupProportionEvents <- stats::aggregate(
        . ~ grouping,
        completeCasesClinicalObject[, c("vitalStatus", "grouping")],
        mean
    )$vitalStatus - 1

    ## Now check whether the group with the lowest or highest
    ## expression/methylation had higher event proportion occuring:
    if (expressionOrMethylation == "Expression") {
        highestEventProportionGroup <- ifelse(
            groupProportionEvents[length(groupProportionEvents)] >
                groupProportionEvents[1],
            "highestExpressionLowSurvival",
            ifelse(
                groupProportionEvents[length(groupProportionEvents)] <
                    groupProportionEvents[1],
                "lowestExpressionLowSurvival",
                "unclear"
            )
        )
    } else if (expressionOrMethylation == "Methylation") {
        highestEventProportionGroup <- ifelse(
            groupProportionEvents[length(groupProportionEvents)] >
                groupProportionEvents[1],
            "highestMethylationLowSurvival",
            ifelse(
                groupProportionEvents[length(groupProportionEvents)] <
                    groupProportionEvents[1],
                "lowestMethylationLowSurvival",
                "unclear"
            )
        )
    }

    ## Add a numeric value for the groups to the completeCasesClinicalObject.
    ## Samples with the lowest grouping will be 1, incrementing by 1 for the
    ## next highest expression/methylation group.
    completeCasesClinicalObject$groupingNumerical <- as.numeric(
        completeCasesClinicalObject$grouping
    )

    ## Now remove the NA values for the purposes of the survival analysis
    completeCasesClinicalObjectNoNA <- completeCasesClinicalObject[
        !is.na(completeCasesClinicalObject$grouping),
    ]

    ## Start with KM statistics, comparing only the highest to lowest groups

    ## Get the completeCasesClinicalObjectNoNA dataset with just the highest and
    ## lowest groups
    completeCasesClinicalObjectNoNAKM <- completeCasesClinicalObjectNoNA[
        completeCasesClinicalObjectNoNA$groupingNumerical %in% c(
            1,
            max(completeCasesClinicalObjectNoNA$groupingNumerical)
        ),
    ]

    ## Create a survival object for KM analyses
    KMSurvivalObject <- survival::Surv(
        completeCasesClinicalObjectNoNAKM$time,
        completeCasesClinicalObjectNoNAKM$vitalStatus
    )
    rownames(KMSurvivalObject) <- rownames(completeCasesClinicalObjectNoNAKM)

    ## Perform the survival analysis.
    ## This uses the expression grouping as the x variable
    ## and the KM survival object as the y variable to create
    ## a table with information about the test
    ## including chi-squared p-value
    KMSurvivalTable <- survival::survdiff(
        KMSurvivalObject ~ completeCasesClinicalObjectNoNAKM$grouping
    )

    ## Get the chi-squared test statistic from the analysis above
    KMChiSquared <- KMSurvivalTable$chisq

    ## Calculate a p-value based on the test statistic to get
    ## a precise p-value for KM analysis
    KMSurvivalPvalue <- as.numeric(
        1 - stats::pchisq(abs(KMChiSquared), df = 1)
    )

    ## Next do Cox Regression analyses considering the groups both as a
    ## categorical and numerical variable

    ## Create a survival object for Cox analyses
    coxSurvivalObject <- survival::Surv(
        completeCasesClinicalObjectNoNA$time,
        completeCasesClinicalObjectNoNA$vitalStatus
    )
    rownames(coxSurvivalObject) <- rownames(completeCasesClinicalObjectNoNA)

    ## Set up separate Cox regression objects, one with the groups as a
    ## categorical variable and one with them as a continuous numerical variable
    ## set up earlier. suppressWarnings is here as the Loglik may not converge
    ## for some analyses.
    coxSurvivalTableCategorical <- suppressWarnings(survival::coxph(
        coxSurvivalObject ~ completeCasesClinicalObjectNoNA$grouping
    ))

    coxSurvivalTableContinuousNumerical <- suppressWarnings(survival::coxph(
        coxSurvivalObject ~ completeCasesClinicalObjectNoNA$groupingNumerical
    ))

    ## Get the coefficient, Hazard Ratio (exp(coeff)) and individual p-values
    ## for each of the groups in the Categorical analysis.
    ## Hazard ratios are found in the second column of the
    ## summary()$coefficients object, while individual p values are found in the
    ## 5th
    coxSurvivalCategoricalGroupCoefficients <- unname(
        summary(
            coxSurvivalTableCategorical
        )$coefficients[, 1]
    )

    coxSurvivalCategoricalGroupHazardRatios <- unname(
        summary(
            coxSurvivalTableCategorical
        )$coefficients[, 2]
    )

    coxSurvivalCategoricalGroupPvalues <- unname(
        summary(
            coxSurvivalTableCategorical
        )$coefficients[, 5]
    )

    ## Get the Coefficient and Hazard Ratio (exp(coeff)) for the continuous
    ## numerical treatment of the groups
    ## Singular in variable name here since only one value is present
    coxSurvivalContinuousNumericalGroupCoefficient <- unname(
        summary(
            coxSurvivalTableContinuousNumerical
        )$coefficients[, 1]
    )

    coxSurvivalContinuousNumericalGroupHazardRatio <- unname(
        summary(
            coxSurvivalTableContinuousNumerical
        )$coefficients[, 2]
    )

    coxSurvivalContinuousNumericalGroupPvalue <- unname(
        summary(
            coxSurvivalTableContinuousNumerical
        )$coefficients[, 5]
    )

    ## Get the p-values for the overall models
    ## p-values are found in the 3rd element of each of the summary()$"testname"
    ## objects
    coxSurvivalCategoricalLikelihoodRatioPvalue <- unname(summary(
        coxSurvivalTableCategorical
    )$logtest[3])

    coxSurvivalCategoricalScoreLogRankPvalue <- unname(summary(
        coxSurvivalTableCategorical
    )$sctest[3])

    coxSurvivalCategoricalWaldPvalue <- unname(summary(
        coxSurvivalTableCategorical
    )$waldtest[3])

    coxSurvivalContinuousNumericalLikelihoodRatioPvalue <- unname(summary(
        coxSurvivalTableContinuousNumerical
    )$logtest[3])

    coxSurvivalContinuousNumericalScoreLogRankPvalue <- unname(summary(
        coxSurvivalTableContinuousNumerical
    )$sctest[3])

    coxSurvivalContinuousNumericalWaldPvalue <- unname(summary(
        coxSurvivalTableContinuousNumerical
    )$waldtest[3])

    ## Create and return KM plots if createPlot is TRUE
    ## Otherwise, return vector of statistics for both KM and Cox
    ## to later combine into a data frame
    if (createPlot) {
        ## Create a survfit formatted survival object for the Cox Regression
        ## analysis with the grouping considered as a categorical variable:
        survfitObject <- survival::survfit(
            survival::Surv(
                time,
                vitalStatus
            ) ~ grouping,
            data = completeCasesClinicalObjectNoNA
        )

        ## Create a vector of legend labels which note the group names, and the
        ## number of samples present in each group:
        legendLabels <- paste0(
            levels(completeCasesClinicalObjectNoNA$grouping),
            "_(n=",
            groupCounts,
            ")"
        )

        ## Format two of the p-values for display in the title
        ## Round the p-value displayed the graph to 3 digits
        coxSurvivalCategoricalScoreLogRankPvalueFormatted <- formatC(
            coxSurvivalCategoricalScoreLogRankPvalue,
            format = "e",
            digits = 3
        )

        KMSurvivalPvalueFormatted <- formatC(
            KMSurvivalPvalue,
            format = "e",
            digits = 3
        )

        ## Create a title to display:
        ## Create the plot title
        ## with gene name and p-value included
        ## If expression is specified, include the gene name and ESNG.
        ## If it's methylation, include only the methylation site ID
        if (expressionOrMethylation == "Expression") {
            SurvivalTitle <- paste0(
                inputName,
                " - ",
                inputID,
                "\nCox Regression Log-Rank p = ",
                coxSurvivalCategoricalScoreLogRankPvalueFormatted,
                "\nKaplan-Meier Lowest vs. Highest p = ",
                KMSurvivalPvalueFormatted
            )
        } else if (expressionOrMethylation == "Methylation") {
            SurvivalTitle <- paste0(
                inputName,
                "\nCox Regression Log-Rank p = ",
                coxSurvivalCategoricalScoreLogRankPvalueFormatted,
                "\nKaplan-Meier Lowest vs. Highest p = ",
                KMSurvivalPvalueFormatted
            )
        }

        ## Create the plot:
        .newInvisibleRecordablePlot()

        basePlot <- survminer::ggsurvplot(
            survfitObject,
            censor.shape = "",
            size = 1,
            xlab = "Time",
            ylab = "Survival Proportion",
            legend.labs = legendLabels,
            title = SurvivalTitle,
            data = completeCasesClinicalObjectNoNA
        )

        basePlot$plot + ggplot2::theme(
            plot.title = ggplot2::element_text(
                hjust = 0.5, face = "bold", size = 8
            ),
            legend.text = ggplot2::element_text(size = 8)
        )

        ## Save the plot to an object
        SurvivalPlotObject <- .recordTENETSavedSizePlot()

        ## Close the plot
        grDevices::dev.off()

        ## Return the plot
        return(SurvivalPlotObject)
    } else {
        ## Convert NaN values for samples lacking a group to NA
        ## (Which can happen as there may be no samples without a group)
        if (is.nan(NAGroupMean)) {
            NAGroupMean <- NA
            NAProportionEvent <- NA
        }

        ## Assemble a vector of results with information relevant to
        ## the given gene/RE DNA methylation site
        survivalReturnVector <- c(
            controlSampleN,
            caseSampleN,
            NACountControlInputValue,
            NACountCaseInputValue,
            controlMeanInputValue,
            caseMeanInputValue,
            NACountCaseClinical,
            controlPresentSampleN,
            casePresentSampleN,
            cutoffVector,
            NACountGroup,
            groupCounts,
            NAGroupMean,
            groupMeans,
            NAProportionEvent,
            groupProportionEvents,
            highestEventProportionGroup,
            KMSurvivalPvalue,
            coxSurvivalCategoricalGroupCoefficients,
            coxSurvivalCategoricalGroupHazardRatios,
            coxSurvivalCategoricalGroupPvalues,
            coxSurvivalContinuousNumericalGroupCoefficient,
            coxSurvivalContinuousNumericalGroupHazardRatio,
            coxSurvivalContinuousNumericalGroupPvalue,
            coxSurvivalCategoricalLikelihoodRatioPvalue,
            coxSurvivalCategoricalScoreLogRankPvalue,
            coxSurvivalCategoricalWaldPvalue,
            coxSurvivalContinuousNumericalLikelihoodRatioPvalue,
            coxSurvivalContinuousNumericalScoreLogRankPvalue,
            coxSurvivalContinuousNumericalWaldPvalue
        )

        namesTemplate <- c(
            "controlSampleCount",
            "caseSampleCount",
            "controlSampleCount@TYPE@Missing",
            "caseSampleCount@TYPE@Missing",
            "controlMean@TYPE@",
            "caseMean@TYPE@",
            "caseSampleCountClinicalMissing",
            "controlSampleCountWithData",
            "caseSampleCountWithData",
            c(
                rbind(
                    paste0(
                        levels(completeCasesClinicalObject$grouping),
                        "_min@TYPE@CutoffValue"
                    ),
                    paste0(
                        levels(completeCasesClinicalObject$grouping),
                        "_max@TYPE@CutoffValue"
                    )
                )
            ),
            "caseSampleCountWithout@TYPE@Group",
            paste0(
                "caseSampleCount_",
                levels(completeCasesClinicalObject$grouping),
                "_@TYPE@Group"
            ),
            "caseMean@TYPE@Without@TYPE@Group",
            paste0(
                "caseMean@TYPE@_",
                levels(completeCasesClinicalObject$grouping),
                "_@TYPE@Group"
            ),
            "caseProportionEventWithout@TYPE@Group",
            paste0(
                "caseProportionEvent_",
                levels(completeCasesClinicalObject$grouping),
                "_@TYPE@Group"
            ),
            "KMSurvivalDirectionOfEffect",
            "KMSurvivalPValue",
            paste0(
                "coxRegressionGroupCategoricalAnalysis_",
                levels(completeCasesClinicalObject$grouping)[
                    2:length(levels(completeCasesClinicalObject$grouping))
                ],
                "_Coefficient"
            ),
            paste0(
                "coxRegressionGroupCategoricalAnalysis_",
                levels(completeCasesClinicalObject$grouping)[
                    2:length(levels(completeCasesClinicalObject$grouping))
                ],
                "_HazardRatio"
            ),
            paste0(
                "coxRegressionGroupCategoricalAnalysis_",
                levels(completeCasesClinicalObject$grouping)[
                    2:length(levels(completeCasesClinicalObject$grouping))
                ],
                "_GroupPvalue"
            ),
            "coxRegressionGroupContinuousNumericalAnalysisCoefficient",
            "coxRegressionGroupContinuousNumericalAnalysisHazardRatio",
            "coxRegressionGroupContinuousNumericalAnalysisGroupPvalue",
            "coxRegressionGroupCategoricalAnalysisLikelihoodRatioPvalue",
            "coxRegressionGroupCategoricalAnalysisScore(LogRank)Pvalue",
            "coxRegressionGroupCategoricalAnalysisWaldPvalue",
            "coxRegressionGroupContinuousNumericalAnalysisLikelihoodRatioPvalue",
            "coxRegressionGroupContinuousNumericalAnalysisScore(LogRank)Pvalue",
            "coxRegressionGroupContinuousNumericalAnalysisWaldPvalue"
        )

        if (expressionOrMethylation == "Expression") {
            ## For expression, add gene ID and gene name columns
            survivalReturnVector <- c(
                inputID,
                inputName,
                survivalReturnVector
            )

            namesTemplate <- c(
                "geneID",
                "geneName",
                namesTemplate
            )
        }

        ## Add names to the return values
        names(survivalReturnVector) <- gsub(
            "@TYPE@", expressionOrMethylation, namesTemplate
        )

        ## Return the vector:
        return(survivalReturnVector)
    }
}

## Internal function to return survival statistics or graphs for a given
## quadrant
.returnSurvivalStatisticsOrGraphs <- function(
    hyperHypo,
    geneIDdf,
    clinicalObject,
    TENETMultiAssayExperiment,
    topGeneNumber,
    geneOrTF, ## Return info for top genes ("Gene") or TFs ("TF")
    ## Return results for genes ("Genes") or RE DNA methylation sites linked to
    ## genes ("DNAMethylationSites")
    genesOrMethSites,
    statsOrPlots, ## Return stats ("Stats") or plots ("Plots")
    survivalGroupingCutoffs,
    useJenksBreaks,
    jenksBreaksGroupCount,
    coreCount) {
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
    topQuadrantGeneName <- geneIDdf[topQuadrantGeneOrTFIDs, "geneName"]

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

    ## Perform and return analyses for genes
    if (genesOrMethSites == "Genes") {
        if (statsOrPlots == "Stats") {
            ## Return survival statistics for genes
            returnValue <- as.data.frame(
                do.call(
                    rbind,
                    parallel::mclapply(
                        X = topQuadrantGeneOrTFIDs,
                        FUN = .survivalFunction,
                        expressionOrMethylation = "Expression",
                        geneIDdf = geneIDdf,
                        clinicalObject = clinicalObject,
                        TENETMultiAssayExperiment = TENETMultiAssayExperiment,
                        survivalGroupingCutoffs = survivalGroupingCutoffs,
                        useJenksBreaks = useJenksBreaks,
                        jenksBreaksGroupCount = jenksBreaksGroupCount,
                        createPlot = FALSE,
                        mc.cores = coreCount
                    )
                )
            )

            rownames(returnValue) <- topQuadrantGeneOrTFIDs

            return(returnValue)
        } else {
            ## Return survival plots for genes
            returnValue <- parallel::mclapply(
                X = topQuadrantGeneOrTFIDs,
                FUN = .survivalFunction,
                expressionOrMethylation = "Expression",
                geneIDdf = geneIDdf,
                clinicalObject = clinicalObject,
                TENETMultiAssayExperiment = TENETMultiAssayExperiment,
                survivalGroupingCutoffs = survivalGroupingCutoffs,
                useJenksBreaks = useJenksBreaks,
                jenksBreaksGroupCount = jenksBreaksGroupCount,
                createPlot = TRUE,
                mc.cores = coreCount
            )

            names(returnValue) <- topQuadrantGeneOrTFIDs

            return(returnValue)
        }
    } else {
        if (statsOrPlots == "Stats") {
            ## For RE DNA methylation sites, initialize a data frame with the
            ## methylation site IDs
            returnValue <- data.frame(
                "DNAMethylationSiteID" =
                    quadrantMethSitesLinkedToSignificantGenes,
                stringsAsFactors = FALSE
            )

            ## Add columns to that data frame indicating which of the RE DNA
            ## methylation sites is linked to each of the top genes
            for (i in seq_along(topQuadrantGeneOrTFIDs)) {
                ## Identify if the quadrantMethSitesLinkedToSignificantGenes
                ## are among RE DNA methylation sites linked to the specific
                ## gene of interest
                TFVector <- quadrantMethSitesLinkedToSignificantGenes %in%
                    TENETMultiAssayExperiment@metadata$step5OptimizeLinks[[
                        quadrantResultsName
                    ]][
                        TENETMultiAssayExperiment@
                        metadata$step5OptimizeLinks[[
                            quadrantResultsName
                        ]]$geneID %in% topQuadrantGeneOrTFIDs[i],
                        "DNAMethylationSiteID"
                    ]

                returnValue[i + 1] <- TFVector
            }

            ## Reset the colnames of the DF as the methylation site IDs,
            ## then the combined gene names and IDs
            colnames(returnValue) <- c(
                "DNAMethylationSiteID",
                paste(
                    topQuadrantGeneName,
                    topQuadrantGeneOrTFIDs,
                    "linked",
                    sep = "_"
                )
            )

            ## Return survival statistics for RE DNA methylation sites linked to
            ## top genes
            returnValue2 <- as.data.frame(
                do.call(
                    rbind,
                    parallel::mclapply(
                        X = quadrantMethSitesLinkedToSignificantGenes,
                        FUN = .survivalFunction,
                        expressionOrMethylation = "Methylation",
                        clinicalObject = clinicalObject,
                        TENETMultiAssayExperiment = TENETMultiAssayExperiment,
                        survivalGroupingCutoffs = survivalGroupingCutoffs,
                        useJenksBreaks = useJenksBreaks,
                        jenksBreaksGroupCount = jenksBreaksGroupCount,
                        createPlot = FALSE,
                        mc.cores = coreCount
                    )
                )
            )

            ## Combine the two return data frames before doing the final return
            returnValue <- cbind(returnValue, returnValue2)

            rownames(returnValue) <-
                quadrantMethSitesLinkedToSignificantGenes

            return(returnValue)
        } else {
            ## Return survival plots for RE DNA methylation sites linked to top
            ## genes
            returnValue <- parallel::mclapply(
                X = quadrantMethSitesLinkedToSignificantGenes,
                FUN = .survivalFunction,
                expressionOrMethylation = "Methylation",
                clinicalObject = clinicalObject,
                TENETMultiAssayExperiment = TENETMultiAssayExperiment,
                survivalGroupingCutoffs = survivalGroupingCutoffs,
                useJenksBreaks = useJenksBreaks,
                jenksBreaksGroupCount = jenksBreaksGroupCount,
                createPlot = TRUE,
                mc.cores = coreCount
            )

            names(returnValue) <- quadrantMethSitesLinkedToSignificantGenes

            return(returnValue)
        }
    }
}

#' Perform Kaplan-Meier and Cox regression analyses to assess the association
#' of top gene expression and linked RE DNA methylation site methylation with
#' patient survival
#'
#' This function takes the top genes/transcription factors (TFs) by number of
#' linked RE DNA methylation sites identified by the
#' `step6DNAMethylationSitesPerGeneTabulation` function up to the number
#' specified by the user and generates survival plots and tables with statistics
#' from survival analyses assessing the survival association of the expression
#' level of each gene as well as the methylation level of each RE DNA
#' methylation site linked to them, using groupings based on either percentile
#' cutoffs or Jenks natural breaks as specified by the user, for Kaplan-Meier
#' analyses.
#'
#' @param TENETMultiAssayExperiment Specify a MultiAssayExperiment object
#' containing expression and methylation SummarizedExperiment objects, such as
#' one created by the TCGADownloader function. This MultiAssayExperiment object
#' should also contain the results from the
#' `step2GetDifferentiallyMethylatedSites`, `step5OptimizeLinks`, and
#' `step6DNAMethylationSitesPerGeneTabulation` functions in
#' its metadata. Additionally, the colData of the MultiAssay object must contain
#' a 'vital_status' and 'time' column, containing data on the patients' survival
#' status and time to event/censorship, respectively.
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
#' @param hypermethGplusAnalysis Set to TRUE to perform survival analyses on
#' top genes/TFs by most hypermethylated RE DNA methylation sites with G+ links,
#' as well as their linked RE DNA methylation sites.
#' @param hypomethGplusAnalysis Set to TRUE to perform survival analyses on
#' top genes/TFs by most hypomethylated RE DNA methylation sites with G+ links,
#' as well as their linked RE DNA methylation sites.
#' @param topGeneNumber Specify the number of top genes/TFs, based on the most
#' linked RE DNA methylation sites of a given analysis type, for which to
#' perform survival analyses. Defaults to 10.
#' @param vitalStatusData Specify the vital status data for samples in the
#' TENETMultiAssayExperiment. Vital status should be given in the form of
#' either "alive" or "dead" (case-insensitive), or 1 or 2, indicating that the
#' sample was collected from a patient who was alive/censored or dead/reached
#' the outcome of interest, respectively. These
#' data can be given in a variety of forms, including a vector, data frame/
#' matrix, or a path to a file that contains the vital status data. If a vector
#' is given, the names of the vector elements must correspond to the names of
#' the samples in the rownames of the colData of the TENETMultiAssayExperiment
#' object. If no names are provided for the vector, then the number of elements
#' in the vector must equal the number of samples in the colData, and are
#' assumed to align with the samples as they are ordered in the colData. If a
#' data frame or matrix is given, its rownames must include the sample names as
#' they appear in the colData of the TENETMultiAssayExperiment object, and its
#' first column must include the vital status data. If a single string
#' is provided, then it is assumed to be a path to a tab-delimited file
#' containing vital status data in the second column, and the names of the
#' samples, again corresponding with the sample names in the colData, in the
#' first column, which will be loaded as the row names. The first row of the
#' file must contain column names. If this variable is set to NA, then
#' the vital status data will be assumed to already be contained in the colData
#' of the TENETMultiAssayExperiment under a column titled "vital_status".
#' Defaults to NA.
#' @param survivalTimeData Specify the survival time data for samples in the
#' TENETMultiAssayExperiment. Survival time should be given in the form of a
#' numeric variable. These data can be given in a variety of forms, including a
#' vector, data frame/matrix, or a path to a file that contains the survival
#' time data. See the documentation for the vitalStatusData argument for more
#' information. If this variable is set to NA, then the survival time
#' data will be assumed to already be contained in the colData of the
#' TENETMultiAssayExperiment under a column titled "time". Defaults to NA.
#' @param highProportion Set a number ranging from 0 to 1, indicating the
#' proportion of all samples to include in the high expression/methylation
#' group for Kaplan-Meier survival analyses. If values are specified for this
#' and `lowProportion`, splitting the samples in this manner will supersede
#' any arguments which are given for `survivalGroupingCutoffs`, `useJenksBreaks`
#' and `jenksBreaksGroupCount`. Defaults to 0.5.
#' @param lowProportion Set a number ranging from 0 to 1, indicating the
#' proportion of all samples to include in the low expression/methylation
#' group for Kaplan-Meier survival analyses. The total value of the
#' highProportion and lowProportion arguments should not exceed 1. **Note:**
#' If both this value and the highProportion value are set to 0.5, samples at
#' exactly the 50th percentile will be assigned to the "Low" group. If values
#' are specified for this and `highProportion`, splitting the samples in this
#' manner will supersede any arguments which are given for
#' `survivalGroupingCutoffs`, `useJenksBreaks` and `jenksBreaksGroupCount`.
#' Defaults to 0.5.
#' @param survivalGroupingCutoffs Specify a data frame or matrix object with two
#' columns and n rows, where n represents the number of groups the expression/
#' methylation samples should be broken into. Values in the object should range
#' from 0 to 1, reflecting the proportion of samples to include in each given
#' group. Values in the first column should reflect the minimum proportion to
#' include in each group, while values in the second column should reflect the
#' max proportion (up to, but not including) for samples in the group. Row names
#' can be given to the object to specify the names the user wishes to be used
#' for the groups. If a valid data frame or matrix is given, it will supersede
#' any arguments given for `useJenksBreaks` and `jenksBreaksGroupCount`.
#' Defaults to NA (to not use custom grouping).
#' @param useJenksBreaks Set to TRUE to automatically set cutoffs for the a
#' number of groups as specified by the `jenksBreaksGroupCount` using Jenks
#' natural breaks optimization. If this is TRUE, a value for the number of
#' groups to be assessed must also be specified for the `jenksBreaksGroupCount`
#' argument. Additionally, the `highProportion`, `lowProportion`, and
#' `survivalGroupingCutoffs` arguments must be NA, as they supersede this
#' analysis type. Defaults to FALSE to not use Jenks breaks.
#' @param jenksBreaksGroupCount Set to a positive integer to specify the number
#' of groups the survival data will be broken into, with cutoffs set between
#' groups using Jenks natural breaks optimization. To use, `useJenksBreaks` must
#' be TRUE. Defaults to NA to not use Jenks breaks.
#' @param generatePlots Set to TRUE to create and save plots displaying the
#' Kaplan-Meier survival results for the genes/TFs of interest, as well as the
#' RE DNA methylation sites linked to them. Defaults to TRUE.
#' @param coreCount Argument passed as the mc.cores argument to mclapply. See
#' `?parallel::mclapply` for more details. Defaults to 1.
#' @return Returns the MultiAssayExperiment object given as the
#' TENETMultiAssayExperiment argument with an additional list of information
#' named 'step7TopGenesSurvival' in its metadata with the output of this
#' function. This list is subdivided into hypermethGplus or hypomethGplus
#' results as selected by the user, which are further subdivided into lists with
#' data for the top overall genes, and for top TF genes only. Each contains a
#' list with data frames containing survival statistics for the top genes/TFs as
#' well as their linked RE DNA methylation sites, from both Kaplan-Meier and Cox
#' regression analyses. Additionally, they will also include a list of
#' Kaplan-Meier plots if generatePlots is TRUE.
#' @export
#'
#' @examplesIf interactive()
#' ## This example uses the example MultiAssayExperiment provided in the
#' ## TENET.ExperimentHub package to perform Kaplan-Meier and Cox regression
#' ## survival analyses for the top 10 genes/TFs, by number of linked hyper- or
#' ## hypomethylated RE DNA methylation sites, as well as for all unique RE DNA
#' ## methylation sites linked to those
#' ## 10 genes/TFs. The vital status and survival time of patients will be
#' ## taken from the "vital_status" and "time" columns present in the colData of
#' ## the example MultiAssayExperiment. Gene names will be retrieved from the
#' ## rowRanges of the 'expression' SummarizedExperiment object in the example
#' ## MultiAssayExperiment. For Kaplan-Meier analyses, the patient samples with
#' ## complete clinical information in the highest half of
#' ## expression/methylation will be compared to the patient samples with
#' ## complete clinical information in the lowest half. Kaplan-Meier plots will
#' ## be saved for the genes and RE DNA methylation sites, and the analysis will
#' ## be performed using one CPU core.
#'
#' ## Load the example TENET MultiAssayExperiment object
#' ## from the TENET.ExperimentHub package
#' exampleTENETMultiAssayExperiment <-
#'     TENET.ExperimentHub::exampleTENETMultiAssayExperiment()
#'
#' ## Use the example dataset to do the survival analysis
#' returnValue <- step7TopGenesSurvival(
#'     TENETMultiAssayExperiment = exampleTENETMultiAssayExperiment
#' )
#'
#' ## This example uses the example MultiAssayExperiment provided in the
#' ## TENET.ExperimentHub package to perform Kaplan-Meier and Cox regression
#' ## survival analyses for only the top 5 genes/TFs, by number of linked
#' ## hypomethylated RE DNA methylation sites, as well as for all unique
#' ## RE DNA methylation sites linked to those 5 genes/TFs only. The vital
#' ## status and survival time of patients will be
#' ## taken from specific columns in a separate data frame with example patient
#' ## data from the TENET.ExperimentHub package. Gene names will be retrieved
#' ## from the rowRanges of the 'expression' SummarizedExperiment object in the
#' ## example MultiAssayExperiment. For Kaplan-Meier analyses, the patient
#' ## samples with complete clinical information in the highest quartile of
#' ## expression/methylation will be compared to the patient samples with
#' ## complete clinical information in the lowest quartile. Kaplan-Meier plots
#' ## will not be saved for the genes and RE DNA methylation sites, and the
#' ## analysis will be performed using 8 CPU cores.
#'
#' ## Load the example TENET MultiAssayExperiment object
#' ## from the TENET.ExperimentHub package
#' exampleTENETMultiAssayExperiment <-
#'     TENET.ExperimentHub::exampleTENETMultiAssayExperiment()
#'
#' ## Again, this loads the data frame with example clinical data for patients
#' ## in the TENET MultiAssayExperiment object from the TENET.ExperimentHub
#' ## package
#' exampleTENETClinicalDataFrame <-
#'     TENET.ExperimentHub::exampleTENETClinicalDataFrame()
#'
#' ## Use the example datasets to do the survival analysis
#' returnValue <- step7TopGenesSurvival(
#'     TENETMultiAssayExperiment = exampleTENETMultiAssayExperiment,
#'     hypermethGplusAnalysis = FALSE,
#'     topGeneNumber = 5,
#'     vitalStatusData = exampleTENETClinicalDataFrame["vital_status"],
#'     survivalTimeData = exampleTENETClinicalDataFrame["time"],
#'     highProportion = 0.25,
#'     lowProportion = 0.25,
#'     generatePlots = FALSE,
#'     coreCount = 8
#' )
#'
#' ## This example uses the example MultiAssayExperiment provided in the
#' ## TENET.ExperimentHub package to perform Kaplan-Meier and Cox regression
#' ## survival analyses for the top 10 genes/TFs, by number of linked hyper- or
#' ## hypomethylated RE DNA methylation sites, as well as for all unique RE DNA
#' ## methylation sites linked to those 10 genes/TFs. The vital status and
#' ## survival time of patients will be taken from the "vital_status" and "time"
#' ## columns present in the colData of the example MultiAssayExperiment. Gene
#' ## names will be retrieved from the rowRanges of the 'expression'
#' ## SummarizedExperiment object in the example MultiAssayExperiment. For
#' ## survival analyses, custom cutoffs representing quartiles will be used, and
#' ## Kaplan-Meier plots will be saved for the genes and RE DNA methylation
#' ## sites, and the analysis will be performed using one CPU core.
#'
#' ## Load the example TENET MultiAssayExperiment object
#' ## from the TENET.ExperimentHub package
#' exampleTENETMultiAssayExperiment <-
#'     TENET.ExperimentHub::exampleTENETMultiAssayExperiment()
#'
#' ## Create a custom cutoffsMatrix which will split the samples into quartiles
#' ## for the purposes of the survival analyses and will also define custom
#' ## names for these groups:
#' cutoffMatrix <- data.frame(
#'     "Low" = c(0, (1 / 4), (1 / 2), (3 / 4)),
#'     "High" = c((1 / 4), (1 / 2), (3 / 4), 1)
#' )
#' rownames(cutoffMatrix) <- c(
#'     "GroupOne",
#'     "GroupTwo",
#'     "GroupThree",
#'     "GroupFour"
#' )
#'
#' ## Use the example dataset and cutoffMatrix to do the survival analysis
#' returnValue <- step7TopGenesSurvival(
#'     TENETMultiAssayExperiment = exampleTENETMultiAssayExperiment,
#'     highProportion = NA,
#'     lowProportion = NA,
#'     survivalGroupingCutoffs = cutoffMatrix
#' )
#'
#' ## This final example uses the example MultiAssayExperiment provided in the
#' ## TENET.ExperimentHub package to perform Kaplan-Meier and Cox regression
#' ## survival analyses for the top 10 genes/TFs, by number of linked hyper- or
#' ## hypomethylated RE DNA methylation sites, as well as for all unique RE DNA
#' ## methylation sites linked to those 10 genes/TFs. The vital status and
#' ## survival time of patients will be taken from the "vital_status" and "time"
#' ## columns present in the colData of the example MultiAssayExperiment. Gene
#' ## names will be retrieved from the rowRanges of the 'expression'
#' ## SummarizedExperiment object in the example MultiAssayExperiment. For
#' ## survival analyses, Jenks natural breaks optimization will be used to
#' ## determine cutoffs for 3 groups in the survival analysis.
#' ## Kaplan-Meier plots will be saved for the genes and RE DNA methylation
#' ## sites, and the analysis will be performed using one CPU core.
#'
#' ## Load the example TENET MultiAssayExperiment object
#' ## from the TENET.ExperimentHub package
#' exampleTENETMultiAssayExperiment <-
#'     TENET.ExperimentHub::exampleTENETMultiAssayExperiment()
#'
#' ## Use the example dataset to do the survival analysis
#' returnValue <- step7TopGenesSurvival(
#'     TENETMultiAssayExperiment = exampleTENETMultiAssayExperiment,
#'     highProportion = NA,
#'     lowProportion = NA,
#'     useJenksBreaks = TRUE,
#'     jenksBreaksGroupCount = 3
#' )
step7TopGenesSurvival <- function(
    TENETMultiAssayExperiment,
    geneAnnotationDataset = NA,
    hypermethGplusAnalysis = TRUE,
    hypomethGplusAnalysis = TRUE,
    topGeneNumber = 10,
    vitalStatusData = NA,
    survivalTimeData = NA,
    highProportion = 0.5,
    lowProportion = 0.5,
    survivalGroupingCutoffs = NA,
    useJenksBreaks = FALSE,
    jenksBreaksGroupCount = NA,
    generatePlots = TRUE,
    coreCount = 1) {
    ## Validate the analysis types and get a vector of the ones selected
    analysisTypes <- .validateAnalysisTypes(
        hypermethGplusAnalysis, hypomethGplusAnalysis
    )

    ## Return an error message if the input MultiAssayExperiment is invalid
    .validateMultiAssayExperiment(
        TENETMultiAssayExperiment,
        needGeneName = is.na(geneAnnotationDataset)
    )

    ## Validate settings of the highProportion, lowProportion,
    ## survivalGroupingCutoffs, useJenksBreaks, and jenksBreaksGroupCount
    ## values.
    ## We want to use the following hierarchy: If the user has specified valid
    ## highProportion and lowProportion values, use them regardless of the other
    ## settings. If they aren't set (are NA), first check to see if the user
    ## has provided a valid survivalGroupingCutoffs dataset. If that hasn't been
    ## provided, then check if the user has specified valid useJenksBreaks
    ## and jenksBreaksGroupCount values.

    ## First let's check if any of the high/lowProportion values are specified
    if (!is.na(highProportion) || !is.na(lowProportion)) {
        ## If only one is NA, return an error
        if (is.na(highProportion) || is.na(lowProportion)) {
            .stopNoCall(
                "One of the highProportion or lowProportion arguments has ",
                "been set, but the other has not been set. Please ensure ",
                "that both arguments are set as a number ranging from 0 to 1, ",
                "with their total not exceeding 1, or that both are set to NA ",
                "while valid settings are provided for the ",
                "survivalGroupingCutoffs, useJenksBreaks, and ",
                "jenksBreaksGroupCount arguments."
            )
        }

        ## Check for nonsensical high/lowProportion values which would cause
        ## invalid results
        if (!is.numeric(highProportion) || !is.numeric(lowProportion)) {
            .stopNoCall(
                "Invalid highProportion and/or lowProportion specified. ",
                "Both must be numeric values between 0 and 1."
            )
        }

        if ((highProportion + lowProportion) > 1 ||
            highProportion <= 0 || lowProportion <= 0
        ) {
            .stopNoCall(
                "Invalid highProportion and/or lowProportion specified. Both ",
                "must be positive, and their sum may not be greater than 1."
            )
        }

        ## If the values look valid, format a survivalGroupingCutoffs matrix for
        ## the cutoffs specified, since they work as if the user has specified
        ## two groups with the specified cutoffs in a survivalGroupingCutoffs
        ## matrix/data frame
        survivalGroupingCutoffs <- data.frame(
            "min" = c(0, highProportion),
            "max" = c(lowProportion, 1)
        )
        rownames(survivalGroupingCutoffs) <- c("low", "high")
    }

    ## Validate the survivalGroupingCutoffs:
    if (!.isSingleNA(survivalGroupingCutoffs)) {
        ## First ensure that if supplied, survivalGroupingCutoffs is either a
        ## matrix or data frame:
        if (
            !inherits(survivalGroupingCutoffs, "matrix") &
                !inherits(survivalGroupingCutoffs, "data.frame")
        ) {
            .stopNoCall(
                "The object given as the survivalGroupingCutoffs argument ",
                "must be either a matrix or data frame."
            )
        }

        ## Next, ensure the matrix/data frame is properly formatted with two
        ## columns:
        if (!ncol(survivalGroupingCutoffs) == 2) {
            .stopNoCall(
                "The survivalGroupingCutoffs object must have two columns, ",
                "the first with the minimum proportion cutoff for each ",
                "group in the rows, and the second with the maximum ",
                "proportion cutoff."
            )
        }

        ## Also ensure there are at least two rows, representing two groups, in
        ## the object
        if (nrow(survivalGroupingCutoffs) < 2) {
            .stopNoCall(
                "The survivalGroupingCutoffs object must have at least two ",
                "rows, representing at least two groups to compare ",
                "in the survival analyses."
            )
        }

        ## Make sure that the values in the specified survivalGroupingCutoffs
        ## are between 0 and 1
        if (min(survivalGroupingCutoffs) < 0 ||
            max(survivalGroupingCutoffs) > 1
        ) {
            .stopNoCall(
                "All values within the survivalGroupingCutoffs object must be ",
                "between 0 and 1, representing the proportion cutoffs for the ",
                "groups in the rows."
            )
        }

        ## Check to make sure the values in the second column are larger than
        ## those in the first column
        if (!all(
            (survivalGroupingCutoffs[, 2] - survivalGroupingCutoffs[, 1]) > 0
        )) {
            .stopNoCall(
                "Values in the second column of the survivalGroupingCutoffs ",
                "object are not all larger than the respective value per row ",
                "in the first column. Since values in the second column ",
                "represent the maximum proportion cutoff for each group (in ",
                "the rows), they should be larger than the values in the ",
                "first column."
            )
        }

        ## Next, ensure that the groups in the survivalGroupingCutoffs are
        ## ordered by increasing value of the minimum proportional cutoffs in
        ## the first column
        survivalGroupingCutoffs <- survivalGroupingCutoffs[
            order(survivalGroupingCutoffs[, 1], decreasing = FALSE),
        ]

        ## Then check that the minimum value of every row past the first in the
        ## dataset is equal or larger than the max value of the previous column.
        ## If the value is not equal to or larger for every row past the first,
        ## issue an error since it implies there are overlaps in the groups
        ## specified. If they are not equal, prepare a warning noting that there
        ## might be gaps between the groups, causing a potential loss of samples
        rowMinEqualOrLargerThanPrevRowMax <- NULL
        rowMinEqualToPrevRowMax <- NULL

        for (i in seq_len(nrow(survivalGroupingCutoffs))) {
            ## If it's the first row, return TRUE, since there is no
            ## previous row to compare it to:
            if (i == 1) {
                rowMinLargerThanPrevRowMax <- TRUE
                rowMinEqualToPrevRowMax <- TRUE
            } else {
                ## Check that the min value in column 1 of the given row is
                ## equal or larger than the max value in column 2 from the
                ## previous row
                rowMinLargerThanPrevRowMax <- c(
                    rowMinLargerThanPrevRowMax,
                    (survivalGroupingCutoffs[i, 1] >=
                        survivalGroupingCutoffs[(i - 1), 2])
                )

                ## Also check that the min value in column 1 of the given row is
                ## equal to the max value in column 2 from the previous row
                rowMinEqualToPrevRowMax <- c(
                    rowMinEqualToPrevRowMax,
                    (survivalGroupingCutoffs[i, 1] ==
                        survivalGroupingCutoffs[(i - 1), 2])
                )
            }
        }

        ## Now, check that all the values in rowMinEqualOrLargerThanPrevRowMax.
        ## If they are not, it implies there is overlap in the group.
        if (!all(rowMinLargerThanPrevRowMax)) {
            .stopNoCall(
                "The proportional values that define each group appear to ",
                "overlap. Please check the values in survivalGroupingCutoffs ",
                "and ensure that the maximum values that define each group in ",
                "column 2 are equal to, or less than, the minimum value of ",
                "the next group."
            )
        }

        ## Now check for potential gaps in the group - if any are detected,
        ## alert the user with a warning that there may be gaps.
        ## This is a warning because the user may want gaps. For example, they
        ## may want to compare the smallest third vs. the largest third.
        if (
            !all(rowMinEqualToPrevRowMax) |
                min(survivalGroupingCutoffs) != 0 |
                max(survivalGroupingCutoffs) != 1
        ) {
            .warningNoCall(
                "There are gaps in the proportional values which define ",
                "each group in the survivalGroupingCutoffs object and some ",
                "may be omitted from the survival analysis. If this was ",
                "unintended, please check the values in the ",
                "survivalGroupingCutoffs object and ensure there are no gaps ",
                "between the maximum cutoff of the previous group in the ",
                "second column and the minimum cutoff of the next group in ",
                "the first column, the lowest minimum cutoff is 0, and the ",
                "highest maximum cutoff is 1."
            )
        }

        ## Finally, if the user has also set useJenksBreaks to TRUE and/or set
        ## a jenksBreaksGroupCount value, warn them that Jenks breaks won't be
        ## used since a survivalGroupingCutoff object has been specified:
        if (useJenksBreaks || !is.na(jenksBreaksGroupCount)) {
            .warningNoCall(
                "Although useJenksBreaks has been set to TRUE and/or a ",
                "jenksBreaksGroupCount has been set, these values will not be ",
                "utilized as a valid survivalGroupingCutoffs object has been ",
                "supplied to define the groups for survival analysis. If you ",
                "would like to use Jenks natural breaks classification to set ",
                "groups instead, please rerun the analysis and set ",
                "survivalGroupingCutoffs to NA or leave the argument undefined."
            )

            ## Set useJenksBreaks to be FALSE and jenksBreaksGroupCount to be NA
            ## to ensure they aren't used later in the function
            useJenksBreaks <- FALSE
            jenksBreaksGroupCount <- NA
        }
    }

    ## Next if supplied, ensure the jenksBreaksGroupCount is a positive whole
    ## number
    if (!is.na(jenksBreaksGroupCount)) {
        ## First, check that the value is a number to begin with:
        if (!is.numeric(jenksBreaksGroupCount)) {
            .stopNoCall(
                "jenksBreaksGroupCount is not a number. This argument must be ",
                "a positive whole number."
            )
        }

        ## Next, check that the value is positive:
        if (jenksBreaksGroupCount <= 0) {
            .stopNoCall(
                "jenksBreaksGroupCount is not positive. This argument must be ",
                "a positive whole number."
            )
        }

        ## Finally, check that the value is a whole number:
        if (jenksBreaksGroupCount %% 1 != 0) {
            .stopNoCall(
                "jenksBreaksGroupCount is not a whole number. This argument ",
                "must be a positive whole number."
            )
        }
    }

    ## Process the status data of the samples if vitalStatusData is not NA.
    ## Otherwise, look into the supplied colData of the MultiAssayExperiment
    vitalStatusResults <- .importClinicalData(
        userInput = vitalStatusData,
        argumentName = "vitalStatusData",
        clinicalDataColumn = "vital_status",
        returnType = "single",
        TENETMultiAssayExperiment = TENETMultiAssayExperiment
    )

    ## Process the survival time data of the samples if vitalStatusData is
    ## not NA. Otherwise, look into the supplied colData of the
    ## MultiAssayExperiment
    survivalTimeResults <- .importClinicalData(
        userInput = survivalTimeData,
        argumentName = "survivalTimeData",
        clinicalDataColumn = "time",
        returnType = "single",
        TENETMultiAssayExperiment = TENETMultiAssayExperiment
    )

    ## Get gene IDs and names from the MAE, or gene annotation dataset if
    ## provided
    geneIDdf <- .getGeneIDsAndNames(
        TENETMultiAssayExperiment, geneAnnotationDataset
    )

    ## Get the names of the control and case samples in the
    ## methylation data first
    methylationSampleNames <- .getExpOrMetSamplesOfType(
        TENETMultiAssayExperiment,
        "methylation",
        namesOnly = TRUE
    )

    ## Get the methylation values that match with expression values
    ## using the mapping data. This assumes the methylation and expression
    ## values share a clinical data match within the mapping.
    metToExpSampleConversion <- .createMetToExpSampleConversionVector(
        TENETMultiAssayExperiment
    )

    ## Use the conversion vector to get the names of the control and case
    ## samples in the expression data which pair with the
    ## methylationSampleNames
    expressionSampleNames <- metToExpSampleConversion[
        methylationSampleNames
    ]

    ## Start creating a clinical data frame by matching the methylation
    ## and expression sample names with their primary names from the sampleMap
    ## of the TENETMultiAssayExperiment
    clinicalDF <- data.frame(
        "methylationSampleNames" = methylationSampleNames,
        "expressionSampleNames" = expressionSampleNames,
        "primarySampleNames" = TENETMultiAssayExperiment@sampleMap[
            match(
                methylationSampleNames,
                TENETMultiAssayExperiment@sampleMap$colname
            ),
            "primary"
        ],
        "sampleType" = TENETMultiAssayExperiment@sampleMap[
            match(
                methylationSampleNames,
                TENETMultiAssayExperiment@sampleMap$colname
            ),
            "sampleType"
        ],
        stringsAsFactors = FALSE
    )

    ## Add the vitalStatus and time variables that were imported.
    ## Convert vitalStatus to lowercase to ease matching below.
    clinicalDF$vitalStatus <- tolower(vitalStatusResults[
        clinicalDF$primarySampleNames,
        "vital_status"
    ])

    clinicalDF$time <- survivalTimeResults[
        clinicalDF$primarySampleNames,
        "time"
    ]

    ## Format the vitalStatus and time columns, inducing NA values where
    ## samples were not input properly

    ## Format vitalStatus. There are three acceptable formats for the values -
    ## "alive"/"dead", 1/2, and "1"/"2" - which we will harmonize to 1/2
    ## (numeric, for alive/censored and dead, respectively).
    ## Note: 1 (numeric) is equal to 1 (as character) - 1=="1" returns TRUE.
    vitalStatusConversionDF <- data.frame(
        "values" = c("alive", "1", "dead", "2"),
        "return" = c(1, 1, 2, 2),
        stringsAsFactors = FALSE
    )
    rownames(vitalStatusConversionDF) <- vitalStatusConversionDF$values

    clinicalDF$vitalStatus <- ifelse(
        clinicalDF$vitalStatus %in% vitalStatusConversionDF$values,
        vitalStatusConversionDF[
            as.character(clinicalDF$vitalStatus),
            "return"
        ],
        NA
    )

    ## Time should be a numeric variable; other values should be NA.
    ## suppressWarnings() is used here, as it will return a warning
    ## of "NAs introduced by coercion" when non-numeric values are
    ## converted to NA, but that is the intended use of this function.
    clinicalDF$time <- suppressWarnings(as.numeric(clinicalDF$time))

    ## Create an empty list to hold the results from this step 7 function
    resultsList <- list()

    ## Do the analysis for the analysis types selected by the user
    for (hyperHypo in analysisTypes) {
        ## Return results for all genes then TFs for each analysis type
        for (geneOrTF in c("Gene", "TF")) {
            ## Return results for genes then RE DNA methylation sites for all
            ## genes and TFs
            for (genesOrMethSites in c("Genes", "DNAMethylationSites")) {
                ## Finally return statistics then plots for genes and RE DNA
                ## methylation sites
                for (statsOrPlots in c("Stats", "Plots")) {
                    resultsList[[
                        paste0(hyperHypo, "methGplusResults")
                    ]][[
                        paste0("top", geneOrTF, "s")
                    ]][[
                        paste0(
                            "top",
                            genesOrMethSites,
                            "Survival",
                            statsOrPlots
                        )
                    ]] <- .returnSurvivalStatisticsOrGraphs(
                        hyperHypo = hyperHypo,
                        geneIDdf = geneIDdf,
                        clinicalObject = clinicalDF,
                        TENETMultiAssayExperiment = TENETMultiAssayExperiment,
                        topGeneNumber = topGeneNumber,
                        geneOrTF = geneOrTF,
                        genesOrMethSites = genesOrMethSites,
                        statsOrPlots = statsOrPlots,
                        survivalGroupingCutoffs = survivalGroupingCutoffs,
                        useJenksBreaks = useJenksBreaks,
                        jenksBreaksGroupCount = jenksBreaksGroupCount,
                        coreCount = coreCount
                    )
                }
            }
        }
    }

    ## Add the results list to the MultiAssayExperiment object
    TENETMultiAssayExperiment@metadata$step7TopGenesSurvival <- resultsList

    return(TENETMultiAssayExperiment)
}
