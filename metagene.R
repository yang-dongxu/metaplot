library(magrittr)
library(ggplot2)
#' prepare annotation features from GTF file
#'
#' @param gtf_file The gtf file, gzip file is accepted.
#' @param features_name The feature name in gtf file,
#' default for ensembl database gtf file
#' features_name= c("five_prime_utr","CDS","three_prime_utr"), if you use
#' ucsc database gtf file, you should change it to c("5UTR","CDS","3UTR").
#'
#' @return GRanges object
#' @export
prepare_features <- function(gtf_file = NULL,
                             features_name = c("five_prime_utr","CDS","three_prime_utr")){
  # assign names for features
  names(features_name) <- c("UTR5","CDS","UTR3")

  # load gtf file
  gtf <- rtracklayer::import.gff(gtf_file,format = "gtf") %>%
    data.frame() %>%
    dplyr::select(seqnames,start,end,width,strand,type,gene_id,gene_name,transcript_id)

  # filter longest transcript
  longest_trans <- gtf %>%
    dplyr::filter(type %in% "exon") %>%
    dplyr::group_by(gene_id,transcript_id) %>%
    dplyr::summarise(trans_len = sum(width)) %>%
    dplyr::slice_head(n = 1)

  # filter features
  features <- gtf %>% dplyr::filter(transcript_id %in% longest_trans$transcript_id) %>%
    dplyr::filter(type %in% features_name)

  # add new types
  # x = 1
  lapply(seq_along(features_name),function(x){
    tmp <- features %>% dplyr::filter(type %in% features_name[x]) %>%
      dplyr::mutate(type = names(features_name)[x])
  }) %>% do.call("rbind",.) -> features

  # get feature total length
  features_len <- features %>% dplyr::group_by(transcript_id,type) %>%
    dplyr::summarise(f_len = sum(width))

  # add transcript length
  features_len <- features_len %>%
    dplyr::left_join(y = longest_trans,by = "transcript_id")

  features <- features %>%
    dplyr::left_join(y = features_len,by = c("transcript_id","type"),
                     relationship = "many-to-many") %>%
    dplyr::group_by(transcript_id) %>%
    dplyr::arrange(seqnames,start,end) %>%
    dplyr::relocate(f_len,.after = type)

  # caculate for positive and negtive strand
  pos_starnd <- features %>% dplyr::filter(strand == "+") %>%
    dplyr::group_by(transcript_id,type) %>%
    dplyr::mutate(tx_len = cumsum(width),.after = "width")

  neg_starnd <- features %>% dplyr::filter(strand == "-") %>%
    dplyr::arrange(seqnames,-start,-end) %>%
    dplyr::group_by(transcript_id,type) %>%
    dplyr::mutate(tx_len = cumsum(width),.after = "width")

  # combine pos and neg
  mer <- rbind(pos_starnd,neg_starnd) %>% GenomicRanges::GRanges()

  return(mer)
}


#' calculate relative position for each features
#'
#' @param bed_file bed format of peaks files.
#' @param features_anno features_anno object from prepare_features functions.
#' @param scale_region whether do scale for each feature length, default FALSE.
#' @param cut_ratio whether use custom ratio to plot feature regions,
#' for eaxample: c(0.1,0.8,0.1), default NULL.
#'
#' @return data frame
#' @export
calculate_relative_position <- function(bed_file = NULL,
                                        features_anno = NULL,
                                        scale_region = FALSE,
                                        cut_ratio = NULL){
  # load peaks data
  result <- try(rtracklayer::import.bed(bed_file),silent = T)

  # whether is bed file
  if(inherits(result, "try-error")){
    peaks <- utils::read.delim(bed_file,header = F)
    colnames(peaks)[1:6] <- c("seqnames","start","end","name","score","strand")
    peaks <- GenomicRanges::GRanges(peaks)
  }else{
    peaks <- rtracklayer::import.bed(bed_file)
  }

  # overlap
  ov <- IRanges::findOverlaps(query = peaks,subject = features_anno)

  lo <- cbind(as.data.frame(peaks[S4Vectors::queryHits(ov)]),
              as.data.frame(features_anno[S4Vectors::subjectHits(ov)]))

  # make unique names
  names(lo) <- make.names(names(lo),unique = T)

  # caculate repative positions
  lo <- lo %>%
    dplyr::mutate(p_mid = as.integer((start + end)/2),.after = end) %>%
    dplyr::filter(p_mid >= start.1 & p_mid <= end.1) %>%
    dplyr::mutate(rel_pos = dplyr::case_when(strand.1 == "+" ~ (p_mid - end.1 + tx_len)/f_len,
                                             strand.1 == "-" ~ (start.1 - p_mid + tx_len)/f_len)) %>%
    dplyr::mutate(rel_pos = dplyr::case_when(type == "CDS" ~ rel_pos + 1,
                                             type %in% c("UTR3") ~ rel_pos + 2,
                                             .default = rel_pos))

  # ==============================================================================
  # calculate scale factor for 5UTR and 3UTR
  # ==============================================================================
  # whether scale features to its length
  if(scale_region == TRUE){
    df_lo <- lo %>% dplyr::select(type,f_len) %>% unique() %>%
      dplyr::group_by(type) %>%
      dplyr::summarise(median_len = stats::median(f_len))

    utr5.sf <- df_lo[which(df_lo$type == "UTR5"),]$median_len / df_lo[which(df_lo$type == "CDS"),]$median_len
    utr3.sf <- df_lo[which(df_lo$type == "UTR3"),]$median_len / df_lo[which(df_lo$type == "CDS"),]$median_len

    # whether use custom defined ratios
    if(!is.null(cut_ratio)){
      utr5.sf <- cut_ratio[1]/cut_ratio[2]
      utr3.sf <- cut_ratio[3]/cut_ratio[2]
    }

    plot_df <- lo %>% dplyr::select(type,rel_pos) %>%
      dplyr::mutate(rel_pos = dplyr::case_when(type == "UTR5" ~ scales::rescale(rel_pos,to = c(1 - utr5.sf,1),from = c(0,1)),
                                               type == "UTR3" ~ scales::rescale(rel_pos,to = c(2,2 + utr3.sf),from = c(2,3)),
                                               .default = rel_pos))

    attr(plot_df,"utr5.sf") <- utr5.sf
    attr(plot_df,"utr3.sf") <- utr3.sf
  }else{
    plot_df <- lo %>% dplyr::select(type,rel_pos)
  }

  return(plot_df)
}



#' metaPlot
#'
#' @param bed_file bed format of peaks files.
#' @param group_names vector of group names.
#' @param features_anno features_anno object from prepare_features functions.
#' @param scale_region whether do scale for each feature length, default FALSE.
#' @param cut_ratio whether use custom ratio to plot feature regions,
#' for eaxample: c(0.1,0.8,0.1).
#' @param facet whether dispaly plot for each group, default FALSE.
#'
#' @return A list with data and plot
#'
#' @import ggplot2
#'
#' @export
metaPlot <- function(bed_file = NULL,
                     group_names = NULL,
                     features_anno = NULL,
                     scale_region = FALSE,
                     cut_ratio = NULL,
                     facet = FALSE){
  # ==============================================================================
  # overlap peaks and anno
  # ==============================================================================
  if(is.null(group_names)){
    group_names <- bed_file
  }

  lapply(seq_along(bed_file),function(x){
    df <- calculate_relative_position(bed_file = bed_file[[x]],
                                      features_anno = features_anno,
                                      scale_region = scale_region,
                                      cut_ratio = cut_ratio)

    # add group name
    df$group_names <- group_names[x]

    return(df)
  }) %>% do.call("rbind",.) -> plot_df

  # ==============================================================================
  # plot
  # ==============================================================================
  utr5.sf <- attr(plot_df,"utr5.sf")
  utr3.sf <- attr(plot_df,"utr3.sf")

  if(scale_region == TRUE){
    label.x <- c((1 - utr5.sf + 1)/2,1.5,(2 + utr3.sf + 2)/2)
  }else{
    label.x <- c(0.5,1.5,2.5)
  }

  # check samples number
  if(length(bed_file) > 1){
    density_layer <- geom_density(aes(x = rel_pos,fill = group_names),color = NA,alpha = 0.5)
  }else{
    density_layer <- geom_density(aes(x = rel_pos),fill = "#CC0033",color = NA)
  }

  pdensity <-
    ggplot(plot_df) +
    density_layer +
    geom_vline(xintercept = c(1,2),lty = "dashed") +
    theme_classic(base_size = 12) +
    scale_y_continuous(expand = expansion(mult = c(0,0.1))) +
    scale_x_continuous(breaks = label.x,
                       labels = c("5'UTR","CDS","3'UTR")) +
    theme(axis.line = element_line(arrow = arrow(length = unit(0.25,"cm"),
                                                 type = "closed")),
          axis.text = element_text(colour = "black"),
          axis.text.x = element_text(face = "bold")) +
    ylab("Peaks density") + xlab("")

  if(facet == TRUE){
    plot <- pdensity +
      facet_wrap(~group_names,scales = "free")
  }else{
    plot <- pdensity
  }

  # output
  res <- list(data = plot_df,
              plot = plot)

  return(res)
}
