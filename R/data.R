#' @title Human transcription factor list
#'
#' @description A character vector of the Ensembl IDs of genes identified
#' as human TFs by Lambert SA et al (PMID: 29425488). Candidate proteins were
#' manually examined by a panel of experts based on available data. Proteins
#' with experimentally demonstrated DNA binding specificity were considered TFs.
#' Other proteins, such as co-factors and RNA binding proteins, were classified
#' as non-TFs. **Citation:** Lambert SA, Jolma A, Campitelli LF, et al. The
#' Human Transcription Factors. Cell. 2018 Feb 8;172(4):650-665. doi:
#' 10.1016/j.cell.2018.01.029. Erratum in: Cell. 2018 Oct 4;175(2):598-599.
#' PMID: 29425488.
#'
#' @usage data("humanTranscriptionFactorList", package = "TENET")
#'
#' @docType data
#'
#' @format A character vector containing 1,639 Ensembl IDs of known human TFs.
#'
#' @examples
#' data("humanTranscriptionFactorList", package = "TENET")
#'
#' @source \url{http://humantfs.ccbr.utoronto.ca/download.php}
"humanTranscriptionFactorList"

#' @title Human transcription factor database
#'
#' @description A data frame with information on genes identified as
#' human TFs by Lambert SA et al (PMID: 29425488). Candidate proteins were
#' manually examined by a panel of experts based on available data. Proteins
#' with experimentally demonstrated DNA binding specificity were considered TFs.
#' Other proteins, such as co-factors and RNA binding proteins, were classified
#' as non-TFs. **Citation:** Lambert SA, Jolma A, Campitelli LF, et al. The
#' Human Transcription Factors. Cell. 2018 Feb 8;172(4):650-665. doi:
#' 10.1016/j.cell.2018.01.029. Erratum in: Cell. 2018 Oct 4;175(2):598-599.
#' PMID: 29425488.
#'
#' @usage data("humanTranscriptionFactorDb", package = "TENET")
#'
#' @format A data frame with 2765 rows and 28 variables.
#' \describe{
#'   \item{\code{Ensembl.ID}}{character Official Ensembl gene ID.}
#'   \item{\code{HGNC.symbol}}{character Official gene name.}
#'   \item{\code{DBD}}{character DNA binding domains contained in protein(s).}
#'   \item{\code{Is.TF.}}{character Is the protein a TF (Yes/No).}
#'   \item{\code{TF.assessment}}{character Assessment of binding activity.}
#'   \item{\code{Binding.mode}}{character Mode of interacting with DNA.}
#'   \item{\code{Motif.status}}{character Current status of motif availability.}
#'   \item{\code{Final.Notes}}{character Final notes, automatically generated.}
#'   \item{\code{Final.Comments}}{character Final comments, manually entered.}
#'   \item{\code{Interpro.ID.s.}}{character Interpro IDs for DBDs.}
#'   \item{\code{EntrezGene.ID}}{character Entrez Gene ID.}
#'   \item{\code{EntrezGene.Description}}{character Entrez Gene Description.}
#'   \item{\code{PDB.ID}}{character Protein Data Bank ID.}
#'   \item{\code{TF.tested.by.HT.SELEX.}}{
#'         character Has the protein been tested for DNA binding in a HT-SELEX
#'         assay in the Taipale lab?
#'   }
#'   \item{\code{TF.tested.by.PBM.}}{
#'         character Has the protein been tested for DNA binding in a PBM assay?
#'   }
#'   \item{\code{Conditional.Binding.Requirements}}{character Notes on
#'     requirements for binding.
#'   }
#'   \item{\code{Original.Comments}}{character Original comments provided by
#'     the primary reviewer of the protein.
#'   }
#'   \item{\code{Vaquerizas.2009.classification}}{character Classification
#'     provided by the Vaquerizas 2009 paper.
#'   }
#'   \item{\code{CisBP.considers.it.a.TF.}}{character Is the protein available
#'     in the CisBP database (build 1.02)?
#'   }
#'   \item{\code{TFCat.classification}}{character Does the TFCat web site
#'     classify the protein as a TF?
#'   }
#'   \item{\code{Is.a.GO.TF.}}{character Does GO (Gene Ontology) classify the
#'     protein as a TF?
#'   }
#'   \item{\code{Initial.assessment}}{character Initial assessment provided by
#'     curators.
#'   }
#'   \item{\code{Curator.1}}{character Name of curator 1.}
#'   \item{\code{Curator.2}}{character Name of curator 2.}
#'   \item{\code{TFclass.considers.it.a.TF.}}{character Does TFclass consider
#'     the protein to be a TF?
#'   }
#'   \item{\code{Go.Evidence}}{character Evidence from GO supporting this
#'     protein being a TF.
#'   }
#'   \item{\code{Pfam.Domains..By.ENSP.ID.}}{character List of Pfam Domains
#'     contained in the protein.
#'   }
#'   \item{\code{Is.C2H2.ZF.KRAB..}}{logical Is the protein a KRAB-containing
#'     Cys2-His2 zinc finger (C2H2-ZF) protein? **Note:** This description is a
#'     guess; Lambert et al did not provide a description for this field.
#'   }
#' }
#'
#' @examples
#' data("humanTranscriptionFactorDb", package = "TENET")
#'
#' @source \url{http://humantfs.ccbr.utoronto.ca/download.php}
"humanTranscriptionFactorDb"
