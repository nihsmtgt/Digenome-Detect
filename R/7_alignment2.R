
seedlength = 12

data=read_csv(paste0(output, ".csv")) 

sgRNAseq_seed = sgRNAseq[(length(sgRNAseq) - seedlength + 1):(length(sgRNAseq))]
PAMambiguous = nchar(PAM) - sum(str_count(as.character(PAM), DNA_BASES))

#alignment parameter
submat = pwalign::nucleotideSubstitutionMatrix()
submat[submat == 0] = -1
submat[submat > 0] = 0
gapOpening = 0
gapExtension = 1.1
Alignment = function(pattern, subject) {
  pairwiseAlignment(
    pattern,
    subject,
    type = "global-local",
    gapOpening = gapOpening,
    gapExtension = gapExtension,
    substitutionMatrix = submat
  )
}

#main
data2 = foreach(
  i = 1:nrow(data),
  .packages = c("Biostrings", "tidyverse", "crisprScore"),
  .combine = dplyr::bind_rows
) %dopar% {
  seq= as.character(data[i, "seq"]) |> DNAString()
  revseq = reverseComplement(seq)
  subseq = substr(seq, 10, 40)
  revsubseq = substr(revseq, 10, 40)
  align = Alignment(sgRNAseq, subseq)
  alignr = Alignment(sgRNAseq, revsubseq)
  if (score(align) < score(alignr)) {
    align = alignr
    seq = revseq
    subseq = revsubseq
  }
  temp = tibble(
      mis_2 = nmismatch(align),
      ins_2 = nindel(align)@insertion[, "WidthSum"],
      del_2 = nindel(align)@deletion[, "WidthSum"]
    )
  gap = 0
  if (nchar(align) - nchar(sgRNAseq) < 0) {
    gap = nchar(sgRNAseq) - nchar(align)
  }
  
  PAMseq = as.character(
    subseq(
      seq,
      start = align@subject@range@start + align@subject@range@width + 9,
      end = align@subject@range@start + align@subject@range@width +
        nchar(PAM) - 1 + 9
    )
  )
  spacerseq = toString(
    subseq(
      seq,
      start = if_else(
        align@subject@range@start + align@subject@range@width - 20 +9 > 0,
        align@subject@range@start + align@subject@range@width - 20+9,
        1
      ),
      end = align@subject@range@start + align@subject@range@width -
        1+9
    )
  )
  temp = mutate(
    temp,
    gap2 = gap,
    distance_2 = mis_2 + del_2 + ins_2 + gap,
    PAM_2 = PAMseq,
    spacer = spacerseq,
    PAMmatch_2 = if_else(countPattern(PAM, DNAString(PAM_2), fixed = FALSE) == 0, FALSE, TRUE)
  )
  if (temp$PAMmatch_2 == F) {
    alignPAM = Alignment(DNAString(paste0(sgRNAseq, PAM)), subseq)
    if (nchar(alignPAM) - nchar(paste0(sgRNAseq, PAM)) < 0) {
      gap = nchar(paste0(sgRNAseq, PAM)) - nchar(alignPAM)
    }else{
      gap=0
    }
    if (nedit(align) >= nedit(alignPAM) - PAMambiguous + gap) {
      PAM_2seq = as.character(
        subseq(
          seq,
          start = alignPAM@subject@range@start + alignPAM@subject@range@width - nchar(PAM)+9,
          end = alignPAM@subject@range@start + alignPAM@subject@range@width -
            1+9
        )
      )
      spacerseq = toString(
        subseq(
          seq,
          start = if_else(
            alignPAM@subject@range@start + alignPAM@subject@range@width - 20 - nchar(PAM) +9> 0,
            alignPAM@subject@range@start + alignPAM@subject@range@width - 20 - nchar(PAM) +9,
            1
          ),
          end = alignPAM@subject@range@start + alignPAM@subject@range@width -
            1 - nchar(PAM)+9
        )
      )
      temp = mutate(
        temp,
        mis_2 = nmismatch(alignPAM) - PAMambiguous,
        ins_2 = nindel(alignPAM)@insertion[, "WidthSum"],
        del_2 = nindel(alignPAM)@deletion[, "WidthSum"],
        distance_2 = mis_2 + del_2 + ins_2 + gap,
        PAM_2 = PAM_2seq,
        spacer = spacerseq,
        PAMmatch_2 = if_else(
          countPattern(PAM, DNAString(PAM_2), fixed = FALSE) == 0,
          FALSE,
          TRUE
        )
      )
    }
  }
  
  salign = Alignment(sgRNAseq_seed, align@subject@unaligned[[1]][align@subject@range@start:(align@subject@range@start +
                                                                                              align@subject@range@width - 1)])
  temp = mutate(
    temp,
    seedscore = score(salign),
    seeddist = nedit(salign),
    seedmis = nmismatch(salign),
    seedins = nindel(salign)@insertion[, "WidthSum"],
    seeddel = nindel(salign)@deletion[, "WidthSum"]
  )
  if (str_length(temp$spacer) == 20 & str_length(temp$PAM_2) == 3) {
    temp = mutate(temp, CFD = as.double(
      getCFDScores(
        str_trunc(as.character(sgRNAseq), 20, side = "left", ellipsis=""), #"TCTGCTGGCGAACCACAACA"
        as.character(spacer),
        as.character(PAM_2)
      )[3]
    ))
  }  else{
    temp = mutate(temp, CFD = -1)
  }
  temp
}

result=cbind(data, data2) 

if(doGGGenome==T){
  check=filter(result, distance!=-1 & distance!=distance_2)
  check2=filter(result, distance!=-1 & PAMmatch!=PAMmatch_2)
  result=mutate(result, PAMmatch_2 = if_else(distance_2 == distance, PAMmatch, PAMmatch_2))
}

write_csv(result,
          paste0("./", output, "_align2.csv"))

