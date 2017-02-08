#' @importFrom GenomeInfoDb seqlevels
#' @importFrom GenomicAlignments seqlevelsInUse
#' @importFrom GenomicRanges match

filterChromosomes <- function(granges,
                              regex_remove=c("^chrUn_", "_random$", "_alt$")) {
    ## Remove non-standard chromosomes.
    chrs <- GenomeInfoDb::seqlevels(granges)
    pattern <- paste0("(", paste0(regex_remove, collapse="|"), ")")
    chrsBad <- grepl(pattern, chrs)
    chrsGood <- chrs[! chrsBad]
    granges <- granges[GenomeInfoDb::seqnames(granges) %in% chrsGood]
    seqlevels(granges) <- seqlevelsInUse(granges)
    granges
}

annotationsConsensus <- function(genome) {
    ## Read genome from local MySQL database and return GRanges object with
    ## consensus annotations.
    conn <- RMySQL::dbConnect(RMySQL::MySQL(), dbname=genome)
    ## Suppress warnings about unsigned integers being imported as numeric.
    knownGenes <- suppressWarnings(DBI::dbReadTable(conn, "refFlat"))
    invisible(DBI::dbDisconnect(conn))
    ## Convert data.frame into GRanges transcripts object.
    gr <- makeGRangesFromDataFrame(knownGenes,
                                   start.field="txStart",
                                   end.field="txEnd")
    gr$gene_id <- knownGenes$name
    ## Remove non-standard chromosomes, otherwise groHMM will throw warnings
    ## about them.
    gr <- filterChromosomes(gr)
    ## Suppress useless informational messages.
    suppressMessages(makeConsensusAnnotations(gr))
}

calibrate <- function(reads, LtProbB, UTS, consensus) {
    # Computes the HMM and returns comparison of HMM with known annotations.
    capture.output(
        hmm <- groHMM::detectTranscripts(reads=reads,
                                         LtProbB=LtProbB,
                                         UTS=UTS,
                                         mc.cores=1),
        file="/dev/null")
    comparison <- groHMM::evaluateHMMInAnnotations(hmm$transcripts, consensus)
    comparison$eval
}

annotationsExpressed <- function(features, reads) {
    ## Return read counts per gene
    featuresLimited <- limitToXkb(features)
    count <- countOverlaps(featuresLimited, reads)
    features <- features[count > 0,]
}
