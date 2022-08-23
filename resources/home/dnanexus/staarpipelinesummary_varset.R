args <- commandArgs(TRUE)
### mandatory
outfile <- args[1]
test.type <- args[2] # Single, Gene_Centric_Coding, Gene_Centric_Noncoding, ncRNA, Sliding_Window
infile.prefix <- args[3]
### optional
known.varlist.rsID <- args[4]
known.varlist.4columns <- args[5]
nullobj.file <- args[6]
agds.file.name <- args[7]
annotation.name.catalog.file <- args[8]
max.maf <- as.numeric(args[9])
use.stepwise.selection <- args[10]
min.maf.adj <- as.numeric(args[11])
QC_label <- args[12]
variant_type <- args[13]
geno_missing_imputation <- args[14]
Annotation_dir <- args[15]
Use_annotation_weights <- args[16]
Annotation_name <- args[17]

test.type.vals <- c("Single", "Gene_Centric_Coding", "Gene_Centric_Noncoding", "ncRNA", "Sliding_Window")
if(!test.type %in% test.type.vals) stop("Error: test.type must be Single, Gene_Centric_Coding, Gene_Centric_Noncoding, ncRNA, or Sliding_Window")
if(Use_annotation_weights == "YES") {
  Use_annotation_weights = TRUE
} else {
  Use_annotation_weights = FALSE
}
Annotation_name <- unlist(strsplit(Annotation_name,","))

cat("Output file prefix:", outfile, "\n")
cat("Type of test:", test.type, "\n")
cat("Prefix of input results:", infile.prefix, "\n")
cat("Known variant list (rsID):", known.varlist.rsID, "\n")
cat("Known variant list (4 columns):", known.varlist.4columns, "\n")
cat("Null model object:", nullobj.file, "\n")
cat("Name of AGDS files:", agds.file.name, "\n")
cat("Name and the corresponding channel name in the AGDS file:", annotation.name.catalog.file, "\n")
cat("Maximum minor allele frequency to be included for variant-set test:", max.maf, "\n")
cat("Whether to perform stepwise selection for variants in known_varlist_rsID and known_varlist_4columns:", use.stepwise.selection, "\n")
cat("Minimum minor allele frequency for a variant to be adjusted for in conditional analysis:", min.maf.adj, "\n")
cat("Channel name of the QC label in the AGDS file:", QC_label, "\n")
cat("Variants included in the analysis:", variant_type, "\n")
cat("Method of handling missing genotypes:", geno_missing_imputation, "\n")
cat("Channel name of the annotations in the AGDS file:", Annotation_dir, "\n")
cat("Use annotations as weights or not:", Use_annotation_weights, "\n")
cat("Annotations used in STAAR:", Annotation_name, "\n")

suppressMessages(library(gdsfmt))
suppressMessages(library(SeqArray))
suppressMessages(library(SeqVarTools))
suppressMessages(library(STAAR))
suppressMessages(library(SCANG))
suppressMessages(library(STAARpipeline))
suppressMessages(library(STAARpipelineSummary))
suppressMessages(library(dplyr))

if (nullobj.file != "NO_NULL_OBJ") {
  ## Null Model
  nullobj <- get(load(nullobj.file))
  if(class(nullobj) == "GENESIS.nullMixedModel") {
    nullmod$sample.id <- row.names(nullmod$model.matrix)
    nullobj <- genesis2staar_nullmodel(nullmod)
    rm(nullmod); gc()
  }
  
  if(known.varlist.rsID != "NO_KNOWN_VARLIST_RSID") {
    known.varlist.rsID <- read.csv(known.varlist.rsID, as.is=T)
    
    known.varlist.rsID_info <- NULL
    
    if ("CHR" %in% colnames(known.varlist.rsID)) {
      for(chr in unique(known.varlist.rsID[,"CHR"]))
      {
        known.varlist.rsID_chr <- known.varlist.rsID[known.varlist.rsID[,1]==chr,]
        
        if (length(strsplit(agds.file.name,"chr")[[1]]) == 2) {
          gds.path <- paste0(strsplit(agds.file.name,"chr")[[1]][1],"chr",chr,strsplit(agds.file.name,"chr")[[1]][2],".gds")
        } else {
          gds.path <- paste0(strsplit(agds.file.name,"chr")[[1]][1],"chr",chr,".gds")
        }
        genofile <- seqOpen(gds.path)
        
        variant.id <- seqGetData(genofile, "variant.id")
        if(annotation.name.catalog.file == "NO_ANNOTATION_NAME_CATALOG_FILE") stop("Error: Annotation name catalog file cannot be missing when the known variant list is in rsID format")
        Annotation_name_catalog <- read.csv(annotation.name.catalog.file, as.is=T)
        rs_num <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="rs_num")]))
        
        rs_num_in <- rs_num%in%known.varlist.rsID_chr[,"SNPS"]
        if (sum(rs_num_in) > 0) {
          variant.id.in <- variant.id[rs_num_in]
          
          rm(rs_num)
          gc()
          
          rm(variant.id)
          gc()
          
          seqSetFilter(genofile,variant.id=variant.id.in)
          
          ### Genotype of Significant LOCI
          position <- as.numeric(seqGetData(genofile, "position"))
          REF <- as.character(seqGetData(genofile, "$ref"))
          ALT <- as.character(seqGetData(genofile, "$alt"))
          
          known.varlist.rsID_info_chr <- data.frame(CHR=rep(chr,length(position)),POS=position,REF=REF,ALT=ALT)
          
          known.varlist.rsID_info <- rbind(known.varlist.rsID_info,known.varlist.rsID_info_chr)
        }
        
        seqClose(genofile)
      }
    } else {
      for(chr in 1:22)
      {
        if (length(strsplit(agds.file.name,"chr")[[1]]) == 2) {
          gds.path <- paste0(strsplit(agds.file.name,"chr")[[1]][1],"chr",chr,strsplit(agds.file.name,"chr")[[1]][2],".gds")
        } else {
          gds.path <- paste0(strsplit(agds.file.name,"chr")[[1]][1],"chr",chr,".gds")
        }
        genofile <- seqOpen(gds.path)
        
        variant.id <- seqGetData(genofile, "variant.id")
        if(annotation.name.catalog.file == "NO_ANNOTATION_NAME_CATALOG_FILE") stop("Error: Annotation name catalog file cannot be missing when the known variant list is in rsID format")
        Annotation_name_catalog <- read.csv(annotation.name.catalog.file, as.is=T)
        rs_num <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="rs_num")]))
        
        rs_num_in <- rs_num%in%known.varlist.rsID[,"SNPS"]
        if (sum(rs_num_in) > 0) {
          variant.id.in <- variant.id[rs_num_in]
          
          rm(rs_num)
          gc()
          
          rm(variant.id)
          gc()
          
          seqSetFilter(genofile,variant.id=variant.id.in)
          
          ### Genotype of Significant LOCI
          position <- as.numeric(seqGetData(genofile, "position"))
          REF <- as.character(seqGetData(genofile, "$ref"))
          ALT <- as.character(seqGetData(genofile, "$alt"))
          
          known.varlist.rsID_info_chr <- data.frame(CHR=rep(chr,length(position)),POS=position,REF=REF,ALT=ALT)
          
          known.varlist.rsID_info <- rbind(known.varlist.rsID_info,known.varlist.rsID_info_chr)
        }
        
        seqClose(genofile)
      }
    }
    
  } else {
    known.varlist.rsID_info <- NULL
  }
  
  if(known.varlist.4columns != "NO_KNOWN_VARLIST_4COLUMNS") {
    known.varlist.4columns <- read.csv(known.varlist.4columns, as.is=T)
    known.varlist.4columns_info <- known.varlist.4columns[,c("CHR","POS","REF","ALT")]
  } else {
    known.varlist.4columns_info <- NULL
  }
  
  known.varlist_info <- rbind(known.varlist.rsID_info,known.varlist.4columns_info)
  
  if (use.stepwise.selection == "YES") {
    cat("Performing stepwise selection for variants in known_varlist_rsID and known_varlist_4columns...\n")
    
    known_loci_genome <- data.frame(CHR=character(),POS=character(),REF=character(),ALT=character(),stringsAsFactors = FALSE)
    
    for(chr in 1:22)
    {
      cat("Chromosome:", chr, "\n")
      
      if (length(strsplit(agds.file.name,"chr")[[1]]) == 2) {
        gds.path <- paste0(strsplit(agds.file.name,"chr")[[1]][1],"chr",chr,strsplit(agds.file.name,"chr")[[1]][2],".gds")
      } else {
        gds.path <- paste0(strsplit(agds.file.name,"chr")[[1]][1],"chr",chr,".gds")
      }
      genofile <- seqOpen(gds.path)
      
      known_loci <- LD_pruning(chr=chr,genofile=genofile,obj_nullmodel=nullobj,variants_list=known.varlist_info,maf_cutoff = min.maf.adj,
                               QC_label=QC_label,geno_missing_imputation=geno_missing_imputation)
      
      known_loci_genome <- rbind(known_loci_genome,known_loci)
      
      seqClose(genofile)
    }
  } else {
    known_loci_genome <- known.varlist_info
  }
  if (!is.null(known_loci_genome)) {
    known_loci_genome <- known_loci_genome[order(known_loci_genome$CHR, known_loci_genome$POS),]
  }
  write.csv(known_loci_genome, file = paste0(outfile, "_known_loci_genome.csv"), row.names = FALSE)
}


if(test.type == "Single") {
	cat("Summarizing single variant test results, the following arguments will be ignored:\n")
  cat("\tMaximum minor allele frequency to be included for variant-set test:", max.maf, "\n")
  rm(max.maf); gc()
  
  results_single_genome <- NULL
  for (chr in 1:22){
    if (file.exists(paste0(infile.prefix, chr, ".Rdata"))){
      cat("Chromosome:", chr, "\n")
      results <- get(load(paste0(infile.prefix, chr, ".Rdata")))
      results_single_genome <- c(results_single_genome, results)
    }
  }
  results_single_genome <- results_single_genome[sapply(results_single_genome, function(x) class(x) != "try-error")]
  results_single_genome <- do.call("rbind",results_single_genome)
  results_single_genome <- results_single_genome[order(results_single_genome$CHR, results_single_genome$POS),]
  results_single_genome_sig <- results_single_genome[results_single_genome$pvalue < 5E-8, , drop = FALSE]
  save(results_single_genome, file = paste0(outfile, "_results_single_genome.Rdata"))
  rm(results_single_genome)
  gc()
  
  if(nullobj.file == "NO_NULL_OBJ") {
    cat("No null model object, only summarize unconditional analysis results...\n")
    write.csv(results_single_genome_sig, file = paste0(outfile, "_results_single_genome_sig.csv"), row.names = FALSE)
  } else {
    cat("Performing conditional analysis for single variant test results...\n")
    results_single_genome_sig_cond <- NULL
    for(chr in 1:22)
    {
      if(sum(results_single_genome_sig$CHR==chr)>=1)
      {
        cat("Chromosome:", chr, "\n")
        results_single_genome_sig_chr <- results_single_genome_sig[results_single_genome_sig$CHR==chr,]
        
        if (length(strsplit(agds.file.name,"chr")[[1]]) == 2) {
          gds.path <- paste0(strsplit(agds.file.name,"chr")[[1]][1],"chr",chr,strsplit(agds.file.name,"chr")[[1]][2],".gds")
        } else {
          gds.path <- paste0(strsplit(agds.file.name,"chr")[[1]][1],"chr",chr,".gds")
        }
        genofile <- seqOpen(gds.path)
        
        results_single_genome_sig_cond_chr <- Individual_Analysis_cond(chr=chr,individual_results=results_single_genome_sig_chr,genofile,obj_nullmodel=nullobj,known_loci=known_loci_genome,
                                                                       QC_label=QC_label,geno_missing_imputation=geno_missing_imputation)
        
        results_single_genome_sig_cond <- rbind(results_single_genome_sig_cond,results_single_genome_sig_cond_chr)
        
        seqClose(genofile)
      }
    }
    
    write.csv(results_single_genome_sig_cond, file = paste0(outfile, "_results_single_genome_sig.csv"), row.names = FALSE)
  }
}


if(test.type == "Gene_Centric_Coding") {
  cat("Summarizing gene-centric test results of coding functional categories...\n")
  
  results_coding_genome <- NULL
  for (chr in 1:22){
    if (file.exists(paste0(infile.prefix, chr, ".Rdata"))){
      cat("Chromosome:", chr, "\n")
      results <- get(load(paste0(infile.prefix, chr, ".Rdata")))
      results_coding_genome <- c(results_coding_genome, results)
    }
  }
  
  results_plof_genome <- NULL
  results_plof_ds_genome <- NULL
  results_missense_genome <- NULL
  results_disruptive_missense_genome <- NULL
  results_synonymous_genome <- NULL
  
  for(kk in 1:length(results_coding_genome))
  {
    for (k in 1:5){
      results <- results_coding_genome[[kk]][k][[1]]
      if(is.null(results)==FALSE)
      {
        ### plof
        if(results[3]=="plof")
        {
          results_plof_genome <- rbind(results_plof_genome,results)
        }
        ### plof_ds
        if(results[3]=="plof_ds")
        {
          results_plof_ds_genome <- rbind(results_plof_ds_genome,results)
        }
        ### missense
        if(results[3]=="missense")
        {
          results_missense_genome <- rbind(results_missense_genome,results)
        }
        ### disruptive_missense
        if(results[3]=="disruptive_missense")
        {
          results_disruptive_missense_genome <- rbind(results_disruptive_missense_genome,results)
        }
        ### synonymous
        if(results[3]=="synonymous")
        {
          results_synonymous_genome <- rbind(results_synonymous_genome,results)
        }
      }
    }
  }
  
  ### plof
  save(results_plof_genome,file = paste0(outfile, "_results_plof_genome.Rdata"))
  ### plof_ds
  save(results_plof_ds_genome,file = paste0(outfile,"_results_plof_ds_genome.Rdata"))
  ### missense
  save(results_missense_genome,file = paste0(outfile,"_results_missense_genome.Rdata"))
  ### disruptive missense
  save(results_disruptive_missense_genome,file = paste0(outfile,"_results_disruptive_missense_genome.Rdata"))
  ### synonymous
  save(results_synonymous_genome,file = paste0(outfile,"_results_synonymous_genome.Rdata"))
  
  ### plof
  plof_sig <- results_plof_genome[results_plof_genome[,"STAAR-O"] < 2.5E-6, , drop = FALSE]
  ### plof_ds
  plof_ds_sig <- results_plof_ds_genome[results_plof_ds_genome[,"STAAR-O"] < 2.5E-6, , drop = FALSE]
  ### missense
  missense_sig <- results_missense_genome[results_missense_genome[,"STAAR-O"] < 2.5E-6, , drop = FALSE]
  ### disruptive_missense
  disruptive_missense_sig <- results_disruptive_missense_genome[results_disruptive_missense_genome[,"STAAR-O"] < 2.5E-6, , drop = FALSE]
  ### synonymous
  synonymous_sig <- results_synonymous_genome[results_synonymous_genome[,"STAAR-O"] < 2.5E-6, , drop = FALSE]
  ### coding
  coding_sig <- rbind(plof_sig[,c("Gene name","Chr","Category","#SNV","SKAT(1,25)","Burden(1,1)","ACAT-V(1,25)","STAAR-O")],
                      missense_sig[,c("Gene name","Chr","Category","#SNV","SKAT(1,25)","Burden(1,1)","ACAT-V(1,25)","STAAR-O")]) 
  coding_sig <- rbind(coding_sig,synonymous_sig[,c("Gene name","Chr","Category","#SNV","SKAT(1,25)","Burden(1,1)","ACAT-V(1,25)","STAAR-O")])
  coding_sig <- rbind(coding_sig,plof_ds_sig[,c("Gene name","Chr","Category","#SNV","SKAT(1,25)","Burden(1,1)","ACAT-V(1,25)","STAAR-O")])
  coding_sig <- rbind(coding_sig,disruptive_missense_sig[,c("Gene name","Chr","Category","#SNV","SKAT(1,25)","Burden(1,1)","ACAT-V(1,25)","STAAR-O")])
  
  write.csv(plof_sig, file = paste0(outfile, "_results_plof_genome_sig.csv"), row.names = FALSE)
  write.csv(plof_ds_sig, file = paste0(outfile, "_results_plof_ds_genome_sig.csv"), row.names = FALSE)
  write.csv(missense_sig, file = paste0(outfile, "_results_missense_genome_sig.csv"), row.names = FALSE)
  write.csv(disruptive_missense_sig, file = paste0(outfile, "_results_disruptive_missense_genome_sig.csv"), row.names = FALSE)
  write.csv(synonymous_sig, file = paste0(outfile, "_results_synonymous_genome_sig.csv"), row.names = FALSE)
  write.csv(coding_sig, file = paste0(outfile, "_results_coding_genome_sig.csv"), row.names = FALSE)
  
  if(nullobj.file == "NO_NULL_OBJ") {
    cat("No null model object, only summarize unconditional analysis results...\n")
  } else {
    if(annotation.name.catalog.file == "NO_ANNOTATION_NAME_CATALOG_FILE") stop("Error: Annotation name catalog file cannot be missing when test.type is Gene_Centric_Coding")
    Annotation_name_catalog <- read.csv(annotation.name.catalog.file, as.is=T)
    cat("Performing conditional analysis for gene-centric test results of coding functional categories...\n")
    ### plof
    plof_sig_cond <- c()
    if(length(plof_sig)!=0)
    {
      if((class(plof_sig)!="matrix")&(class(plof_sig)!="data.frame"))
      {
        plof_sig <- matrix(plof_sig,nrow=1)
      }
      
      for(k in 1:dim(plof_sig)[1])
      {
        chr <- plof_sig[k,2]
        gene_name <- plof_sig[k,1]
        category <- as.character(plof_sig[k,3])
        
        if (length(strsplit(agds.file.name,"chr")[[1]]) == 2) {
          gds.path <- paste0(strsplit(agds.file.name,"chr")[[1]][1],"chr",chr,strsplit(agds.file.name,"chr")[[1]][2],".gds")
        } else {
          gds.path <- paste0(strsplit(agds.file.name,"chr")[[1]][1],"chr",chr,".gds")
        }
        genofile <- seqOpen(gds.path)
        
        res_cond <- Gene_Centric_Coding_cond(chr=chr,gene_name=gene_name,category=category,genofile=genofile,obj_nullmodel=nullobj,known_loci=known_loci_genome,rare_maf_cutoff=max.maf,
                                             QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                                             Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                                             Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name)
        plof_sig_cond <- rbind(plof_sig_cond,res_cond)
        
        seqClose(genofile)
      }
    }
    
    ### plof_ds
    plof_ds_sig_cond <- c()
    if(length(plof_ds_sig)!=0)
    {
      if((class(plof_ds_sig)!="matrix")&(class(plof_ds_sig)!="data.frame"))
      {
        plof_ds_sig <- matrix(plof_ds_sig,nrow=1)
      }
      
      for(k in 1:dim(plof_ds_sig)[1])
      {
        chr <- plof_ds_sig[k,2]
        gene_name <- plof_ds_sig[k,1]
        category <- as.character(plof_ds_sig[k,3])
        
        if (length(strsplit(agds.file.name,"chr")[[1]]) == 2) {
          gds.path <- paste0(strsplit(agds.file.name,"chr")[[1]][1],"chr",chr,strsplit(agds.file.name,"chr")[[1]][2],".gds")
        } else {
          gds.path <- paste0(strsplit(agds.file.name,"chr")[[1]][1],"chr",chr,".gds")
        }
        genofile <- seqOpen(gds.path)
        
        res_cond <- Gene_Centric_Coding_cond(chr=chr,gene_name=gene_name,category=category,genofile=genofile,obj_nullmodel=nullobj,known_loci=known_loci_genome,rare_maf_cutoff=max.maf,
                                             QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                                             Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                                             Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name)
        plof_ds_sig_cond <- rbind(plof_ds_sig_cond,res_cond)
        
        seqClose(genofile)
      }
    }
    
    ### missense
    missense_sig_cond <- c()
    if(length(missense_sig)!=0)
    {
      if((class(missense_sig)!="matrix")&(class(missense_sig)!="data.frame"))
      {
        missense_sig <- matrix(missense_sig,nrow=1)
      }
      
      for(k in 1:dim(missense_sig)[1])
      {
        chr <- missense_sig[k,2]
        gene_name <- missense_sig[k,1]
        category <- as.character(missense_sig[k,3])
        
        if (length(strsplit(agds.file.name,"chr")[[1]]) == 2) {
          gds.path <- paste0(strsplit(agds.file.name,"chr")[[1]][1],"chr",chr,strsplit(agds.file.name,"chr")[[1]][2],".gds")
        } else {
          gds.path <- paste0(strsplit(agds.file.name,"chr")[[1]][1],"chr",chr,".gds")
        }
        genofile <- seqOpen(gds.path)
        
        res_cond <- Gene_Centric_Coding_cond(chr=chr,gene_name=gene_name,category=category,genofile=genofile,obj_nullmodel=nullobj,known_loci=known_loci_genome,rare_maf_cutoff=max.maf,
                                             QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                                             Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                                             Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name)
        missense_sig_cond <- rbind(missense_sig_cond,res_cond)
        
        seqClose(genofile)
      }
    }
    
    ### disruptive_missense
    disruptive_missense_sig_cond <- c()
    if(length(disruptive_missense_sig)!=0)
    {
      if((class(disruptive_missense_sig)!="matrix")&(class(disruptive_missense_sig)!="data.frame"))
      {
        disruptive_missense_sig <- matrix(disruptive_missense_sig,nrow=1)
      }
      
      for(k in 1:dim(disruptive_missense_sig)[1])
      {
        chr <- disruptive_missense_sig[k,2]
        gene_name <- disruptive_missense_sig[k,1]
        category <- as.character(disruptive_missense_sig[k,3])
        
        if (length(strsplit(agds.file.name,"chr")[[1]]) == 2) {
          gds.path <- paste0(strsplit(agds.file.name,"chr")[[1]][1],"chr",chr,strsplit(agds.file.name,"chr")[[1]][2],".gds")
        } else {
          gds.path <- paste0(strsplit(agds.file.name,"chr")[[1]][1],"chr",chr,".gds")
        }
        genofile <- seqOpen(gds.path)
        
        res_cond <- Gene_Centric_Coding_cond(chr=chr,gene_name=gene_name,category=category,genofile=genofile,obj_nullmodel=nullobj,known_loci=known_loci_genome,rare_maf_cutoff=max.maf,
                                             QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                                             Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                                             Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name)
        disruptive_missense_sig_cond <- rbind(disruptive_missense_sig_cond,res_cond)
        
        seqClose(genofile)
      }
    }
    
    ### synonymous
    synonymous_sig_cond <- c()
    if(length(synonymous_sig)!=0)
    {
      if((class(synonymous_sig)!="matrix")&(class(synonymous_sig)!="data.frame"))
      {
        synonymous_sig <- matrix(synonymous_sig,nrow=1)
      }
      
      for(k in 1:dim(synonymous_sig)[1])
      {
        chr <- synonymous_sig[k,2]
        gene_name <- synonymous_sig[k,1]
        category <- as.character(synonymous_sig[k,3])
        
        if (length(strsplit(agds.file.name,"chr")[[1]]) == 2) {
          gds.path <- paste0(strsplit(agds.file.name,"chr")[[1]][1],"chr",chr,strsplit(agds.file.name,"chr")[[1]][2],".gds")
        } else {
          gds.path <- paste0(strsplit(agds.file.name,"chr")[[1]][1],"chr",chr,".gds")
        }
        genofile <- seqOpen(gds.path)
        
        res_cond <- Gene_Centric_Coding_cond(chr=chr,gene_name=gene_name,category=category,genofile=genofile,obj_nullmodel=nullobj,known_loci=known_loci_genome,rare_maf_cutoff=max.maf,
                                             QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                                             Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                                             Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name)
        synonymous_sig_cond <- rbind(synonymous_sig_cond,res_cond)
        
        seqClose(genofile)
      }
    }
    
    ### coding
    coding_sig_cond <- rbind(plof_sig_cond[,c("Gene name","Chr","Category","#SNV","SKAT(1,25)","Burden(1,1)","ACAT-V(1,25)","STAAR-O")],
                             missense_sig_cond[,c("Gene name","Chr","Category","#SNV","SKAT(1,25)","Burden(1,1)","ACAT-V(1,25)","STAAR-O")]) 
    coding_sig_cond <- rbind(coding_sig_cond,synonymous_sig_cond[,c("Gene name","Chr","Category","#SNV","SKAT(1,25)","Burden(1,1)","ACAT-V(1,25)","STAAR-O")])
    coding_sig_cond <- rbind(coding_sig_cond,plof_ds_sig_cond[,c("Gene name","Chr","Category","#SNV","SKAT(1,25)","Burden(1,1)","ACAT-V(1,25)","STAAR-O")])
    coding_sig_cond <- rbind(coding_sig_cond,disruptive_missense_sig_cond[,c("Gene name","Chr","Category","#SNV","SKAT(1,25)","Burden(1,1)","ACAT-V(1,25)","STAAR-O")])
    
    write.csv(plof_sig_cond, file = paste0(outfile, "_results_plof_genome_sig_cond.csv"), row.names = FALSE)
    write.csv(plof_ds_sig_cond, file = paste0(outfile, "_results_plof_ds_genome_sig_cond.csv"), row.names = FALSE)
    write.csv(missense_sig_cond, file = paste0(outfile, "_results_missense_genome_sig_cond.csv"), row.names = FALSE)
    write.csv(disruptive_missense_sig_cond, file = paste0(outfile, "_results_disruptive_missense_genome_sig_cond.csv"), row.names = FALSE)
    write.csv(synonymous_sig_cond, file = paste0(outfile, "_results_synonymous_genome_sig_cond.csv"), row.names = FALSE)
    write.csv(coding_sig_cond, file = paste0(outfile, "_results_coding_genome_sig_cond.csv"), row.names = FALSE)
  }
}


if(test.type == "Gene_Centric_Noncoding") {
  cat("Summarizing gene-centric test results of noncoding functional categories...\n")
  
  results_noncoding_genome <- NULL
  for (chr in 1:22){
    if (file.exists(paste0(infile.prefix, chr, ".Rdata"))){
      cat("Chromosome:", chr, "\n")
      results <- get(load(paste0(infile.prefix, chr, ".Rdata")))
      results_noncoding_genome <- c(results_noncoding_genome, results)
    }
  }
  
  results_UTR_genome <- NULL
  results_upstream_genome <- NULL
  results_downstream_genome <- NULL
  results_promoter_CAGE_genome <- NULL
  results_promoter_DHS_genome <- NULL
  results_enhancer_CAGE_genome <- NULL
  results_enhancer_DHS_genome <- NULL
  
  for(kk in 1:length(results_noncoding_genome))
  {
    for (k in 1:7){
      results <- results_noncoding_genome[[kk]][k][[1]]
      if(is.null(results)==FALSE)
      {
        ### UTR
        if(results[3]=="UTR")
        {
          results_UTR_genome <- rbind(results_UTR_genome,results)
        }
        ### upstream
        if(results[3]=="upstream")
        {
          results_upstream_genome <- rbind(results_upstream_genome,results)
        }
        ### downstream
        if(results[3]=="downstream")
        {
          results_downstream_genome <- rbind(results_downstream_genome,results)
        }
        ### promoter_CAGE
        if(results[3]=="promoter_CAGE")
        {
          results_promoter_CAGE_genome <- rbind(results_promoter_CAGE_genome,results)
        }
        ### promoter_DHS
        if(results[3]=="promoter_DHS")
        {
          results_promoter_DHS_genome <- rbind(results_promoter_DHS_genome,results)
        }
        ### enhancer_CAGE
        if(results[3]=="enhancer_CAGE")
        {
          results_enhancer_CAGE_genome <- rbind(results_enhancer_CAGE_genome,results)
        }
        ### enhancer_DHS
        if(results[3]=="enhancer_DHS")
        {
          results_enhancer_DHS_genome <- rbind(results_enhancer_DHS_genome,results)
        }
      }
    }
  }
  
  ### UTR
  save(results_UTR_genome,file = paste0(outfile, "_results_UTR_genome.Rdata"))
  ### upstream
  save(results_upstream_genome,file = paste0(outfile,"_results_upstream_genome.Rdata"))
  ### downstream
  save(results_downstream_genome,file = paste0(outfile,"_results_downstream_genome.Rdata"))
  ### promoter_CAGE
  save(results_promoter_CAGE_genome,file = paste0(outfile,"_results_promoter_CAGE_genome.Rdata"))
  ### promoter_DHS
  save(results_promoter_DHS_genome,file = paste0(outfile,"_results_promoter_DHS_genome.Rdata"))
  ### enhancer_CAGE
  save(results_enhancer_CAGE_genome,file = paste0(outfile,"_results_enhancer_CAGE_genome.Rdata"))
  ### enhancer_DHS
  save(results_enhancer_DHS_genome,file = paste0(outfile,"_results_enhancer_DHS_genome.Rdata"))
  
  ### UTR
  UTR_sig <- results_UTR_genome[results_UTR_genome[,"STAAR-O"] < 2.5E-6, , drop = FALSE]
  ### upstream
  upstream_sig <- results_upstream_genome[results_upstream_genome[,"STAAR-O"] < 2.5E-6, , drop = FALSE]
  ### downstream
  downstream_sig <- results_downstream_genome[results_downstream_genome[,"STAAR-O"] < 2.5E-6, , drop = FALSE]
  ### promoter_CAGE
  promoter_CAGE_sig <- results_promoter_CAGE_genome[results_promoter_CAGE_genome[,"STAAR-O"] < 2.5E-6, , drop = FALSE]
  ### promoter_DHS
  promoter_DHS_sig <- results_promoter_DHS_genome[results_promoter_DHS_genome[,"STAAR-O"] < 2.5E-6, , drop = FALSE]
  ### enhancer_CAGE
  enhancer_CAGE_sig <- results_enhancer_CAGE_genome[results_enhancer_CAGE_genome[,"STAAR-O"] < 2.5E-6, , drop = FALSE]
  ### enhancer_DHS
  enhancer_DHS_sig <- results_enhancer_DHS_genome[results_enhancer_DHS_genome[,"STAAR-O"] < 2.5E-6, , drop = FALSE]
  ### noncoding
  noncoding_sig <- rbind(UTR_sig[,c("Gene name","Chr","Category","#SNV","SKAT(1,25)","Burden(1,1)","ACAT-V(1,25)","STAAR-O")],
                         upstream_sig[,c("Gene name","Chr","Category","#SNV","SKAT(1,25)","Burden(1,1)","ACAT-V(1,25)","STAAR-O")]) 
  noncoding_sig <- rbind(noncoding_sig,downstream_sig[,c("Gene name","Chr","Category","#SNV","SKAT(1,25)","Burden(1,1)","ACAT-V(1,25)","STAAR-O")])
  noncoding_sig <- rbind(noncoding_sig,promoter_CAGE_sig[,c("Gene name","Chr","Category","#SNV","SKAT(1,25)","Burden(1,1)","ACAT-V(1,25)","STAAR-O")])
  noncoding_sig <- rbind(noncoding_sig,promoter_DHS_sig[,c("Gene name","Chr","Category","#SNV","SKAT(1,25)","Burden(1,1)","ACAT-V(1,25)","STAAR-O")])
  noncoding_sig <- rbind(noncoding_sig,enhancer_CAGE_sig[,c("Gene name","Chr","Category","#SNV","SKAT(1,25)","Burden(1,1)","ACAT-V(1,25)","STAAR-O")])
  noncoding_sig <- rbind(noncoding_sig,enhancer_DHS_sig[,c("Gene name","Chr","Category","#SNV","SKAT(1,25)","Burden(1,1)","ACAT-V(1,25)","STAAR-O")])
  
  write.csv(UTR_sig, file = paste0(outfile, "_results_UTR_genome_sig.csv"), row.names = FALSE)
  write.csv(upstream_sig, file = paste0(outfile, "_results_upstream_genome_sig.csv"), row.names = FALSE)
  write.csv(downstream_sig, file = paste0(outfile, "_results_downstream_genome_sig.csv"), row.names = FALSE)
  write.csv(promoter_CAGE_sig, file = paste0(outfile, "_results_promoter_CAGE_genome_sig.csv"), row.names = FALSE)
  write.csv(promoter_DHS_sig, file = paste0(outfile, "_results_promoter_DHS_genome_sig.csv"), row.names = FALSE)
  write.csv(enhancer_CAGE_sig, file = paste0(outfile, "_results_enhancer_CAGE_genome_sig.csv"), row.names = FALSE)
  write.csv(enhancer_DHS_sig, file = paste0(outfile, "_results_enhancer_DHS_genome_sig.csv"), row.names = FALSE)
  write.csv(noncoding_sig, file = paste0(outfile, "_results_noncoding_genome_sig.csv"), row.names = FALSE)
  
  if(nullobj.file == "NO_NULL_OBJ") {
    cat("No null model object, only summarize unconditional analysis results...\n")
  } else {
    if(annotation.name.catalog.file == "NO_ANNOTATION_NAME_CATALOG_FILE") stop("Error: Annotation name catalog file cannot be missing when test.type is Gene_Centric_Noncoding")
    Annotation_name_catalog <- read.csv(annotation.name.catalog.file, as.is=T)
    cat("Performing conditional analysis for gene-centric test results of noncoding functional categories...\n")
    ### UTR
    UTR_sig_cond <- c()
    if(length(UTR_sig)!=0)
    {
      
      if((class(UTR_sig)!="matrix")&(class(UTR_sig)!="data.frame"))
      {
        UTR_sig <- matrix(UTR_sig,nrow=1)
      }
      
      for(k in 1:dim(UTR_sig)[1])
      {
        chr <- UTR_sig[k,2]
        gene_name <- UTR_sig[k,1]
        category <- as.character(UTR_sig[k,3])
        
        if (length(strsplit(agds.file.name,"chr")[[1]]) == 2) {
          gds.path <- paste0(strsplit(agds.file.name,"chr")[[1]][1],"chr",chr,strsplit(agds.file.name,"chr")[[1]][2],".gds")
        } else {
          gds.path <- paste0(strsplit(agds.file.name,"chr")[[1]][1],"chr",chr,".gds")
        }
        genofile <- seqOpen(gds.path)
        
        res_cond <- Gene_Centric_Noncoding_cond(chr=chr,gene_name=gene_name,category=category,genofile=genofile,obj_nullmodel=nullobj,known_loci=known_loci_genome,rare_maf_cutoff=max.maf,
                                                QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                                                Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                                                Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name)
        UTR_sig_cond <- rbind(UTR_sig_cond,res_cond)
        
        seqClose(genofile)
      }
    }
    
    ### upstream
    upstream_sig_cond <- c()
    if(length(upstream_sig)!=0)
    {
      
      if((class(upstream_sig)!="matrix")&(class(upstream_sig)!="data.frame"))
      {
        upstream_sig <- matrix(upstream_sig,nrow=1)
      }
      
      for(k in 1:dim(upstream_sig)[1])
      {
        chr <- upstream_sig[k,2]
        gene_name <- upstream_sig[k,1]
        category <- as.character(upstream_sig[k,3])
        
        if (length(strsplit(agds.file.name,"chr")[[1]]) == 2) {
          gds.path <- paste0(strsplit(agds.file.name,"chr")[[1]][1],"chr",chr,strsplit(agds.file.name,"chr")[[1]][2],".gds")
        } else {
          gds.path <- paste0(strsplit(agds.file.name,"chr")[[1]][1],"chr",chr,".gds")
        }
        genofile <- seqOpen(gds.path)
        
        res_cond <- Gene_Centric_Noncoding_cond(chr=chr,gene_name=gene_name,category=category,genofile=genofile,obj_nullmodel=nullobj,known_loci=known_loci_genome,rare_maf_cutoff=max.maf,
                                                QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                                                Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                                                Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name)
        upstream_sig_cond <- rbind(upstream_sig_cond,res_cond)
        
        seqClose(genofile)
      }
    }
    
    ### downstream
    downstream_sig_cond <- c()
    if(length(downstream_sig)!=0)
    {
      if((class(downstream_sig)!="matrix")&(class(downstream_sig)!="data.frame"))
      {
        downstream_sig <- matrix(downstream_sig,nrow=1)
      }
      
      for(k in 1:dim(downstream_sig)[1])
      {
        chr <- downstream_sig[k,2]
        gene_name <- downstream_sig[k,1]
        category <- as.character(downstream_sig[k,3])
        
        if (length(strsplit(agds.file.name,"chr")[[1]]) == 2) {
          gds.path <- paste0(strsplit(agds.file.name,"chr")[[1]][1],"chr",chr,strsplit(agds.file.name,"chr")[[1]][2],".gds")
        } else {
          gds.path <- paste0(strsplit(agds.file.name,"chr")[[1]][1],"chr",chr,".gds")
        }
        genofile <- seqOpen(gds.path)
        
        res_cond <- Gene_Centric_Noncoding_cond(chr=chr,gene_name=gene_name,category=category,genofile=genofile,obj_nullmodel=nullobj,known_loci=known_loci_genome,rare_maf_cutoff=max.maf,
                                                QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                                                Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                                                Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name)
        downstream_sig_cond <- rbind(downstream_sig_cond,res_cond)
        
        seqClose(genofile)
      }
    }
    
    ### promoter_CAGE
    promoter_CAGE_sig_cond <- c()
    if(length(promoter_CAGE_sig)!=0)
    {
      if((class(promoter_CAGE_sig)!="matrix")&(class(promoter_CAGE_sig)!="data.frame"))
      {
        promoter_CAGE_sig <- matrix(promoter_CAGE_sig,nrow=1)
      }
      
      for(k in 1:dim(promoter_CAGE_sig)[1])
      {
        chr <- promoter_CAGE_sig[k,2]
        gene_name <- promoter_CAGE_sig[k,1]
        category <- as.character(promoter_CAGE_sig[k,3])
        
        if (length(strsplit(agds.file.name,"chr")[[1]]) == 2) {
          gds.path <- paste0(strsplit(agds.file.name,"chr")[[1]][1],"chr",chr,strsplit(agds.file.name,"chr")[[1]][2],".gds")
        } else {
          gds.path <- paste0(strsplit(agds.file.name,"chr")[[1]][1],"chr",chr,".gds")
        }
        genofile <- seqOpen(gds.path)
        
        res_cond <- Gene_Centric_Noncoding_cond(chr=chr,gene_name=gene_name,category=category,genofile=genofile,obj_nullmodel=nullobj,known_loci=known_loci_genome,rare_maf_cutoff=max.maf,
                                                QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                                                Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                                                Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name)
        promoter_CAGE_sig_cond <- rbind(promoter_CAGE_sig_cond,res_cond)
        
        seqClose(genofile)
      }
    }
    
    ### promoter_DHS
    promoter_DHS_sig_cond <- c()
    if(length(promoter_DHS_sig)!=0)
    {
      if((class(promoter_DHS_sig)!="matrix")&(class(promoter_DHS_sig)!="data.frame"))
      {
        promoter_DHS_sig <- matrix(promoter_DHS_sig,nrow=1)
      }
      
      
      for(k in 1:dim(promoter_DHS_sig)[1])
      {
        chr <- promoter_DHS_sig[k,2]
        gene_name <- promoter_DHS_sig[k,1]
        category <- as.character(promoter_DHS_sig[k,3])
        
        if (length(strsplit(agds.file.name,"chr")[[1]]) == 2) {
          gds.path <- paste0(strsplit(agds.file.name,"chr")[[1]][1],"chr",chr,strsplit(agds.file.name,"chr")[[1]][2],".gds")
        } else {
          gds.path <- paste0(strsplit(agds.file.name,"chr")[[1]][1],"chr",chr,".gds")
        }
        genofile <- seqOpen(gds.path)
        
        res_cond <- Gene_Centric_Noncoding_cond(chr=chr,gene_name=gene_name,category=category,genofile=genofile,obj_nullmodel=nullobj,known_loci=known_loci_genome,rare_maf_cutoff=max.maf,
                                                QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                                                Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                                                Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name)
        promoter_DHS_sig_cond <- rbind(promoter_DHS_sig_cond,res_cond)
        
        seqClose(genofile)
      }
    }
    
    ### enhancer_CAGE
    enhancer_CAGE_sig_cond <- c()
    if(length(enhancer_CAGE_sig)!=0)
    {
      if((class(enhancer_CAGE_sig)!="matrix")&(class(enhancer_CAGE_sig)!="data.frame"))
      {
        enhancer_CAGE_sig <- matrix(enhancer_CAGE_sig,nrow=1)
      }
      
      for(k in 1:dim(enhancer_CAGE_sig)[1])
      {
        chr <- enhancer_CAGE_sig[k,2]
        gene_name <- enhancer_CAGE_sig[k,1]
        category <- as.character(enhancer_CAGE_sig[k,3])
        
        if (length(strsplit(agds.file.name,"chr")[[1]]) == 2) {
          gds.path <- paste0(strsplit(agds.file.name,"chr")[[1]][1],"chr",chr,strsplit(agds.file.name,"chr")[[1]][2],".gds")
        } else {
          gds.path <- paste0(strsplit(agds.file.name,"chr")[[1]][1],"chr",chr,".gds")
        }
        genofile <- seqOpen(gds.path)
        
        res_cond <- Gene_Centric_Noncoding_cond(chr=chr,gene_name=gene_name,category=category,genofile=genofile,obj_nullmodel=nullobj,known_loci=known_loci_genome,rare_maf_cutoff=max.maf,
                                                QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                                                Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                                                Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name)
        enhancer_CAGE_sig_cond <- rbind(enhancer_CAGE_sig_cond,res_cond)
        
        seqClose(genofile)
      }
    }
    
    ### enhancer_DHS
    enhancer_DHS_sig_cond <- c()
    if(length(enhancer_DHS_sig)!=0)
    {
      if((class(enhancer_DHS_sig)!="matrix")&(class(enhancer_DHS_sig)!="data.frame"))
      {
        enhancer_DHS_sig <- matrix(enhancer_DHS_sig,nrow=1)
      }
      
      for(k in 1:dim(enhancer_DHS_sig)[1])
      {
        chr <- enhancer_DHS_sig[k,2]
        gene_name <- enhancer_DHS_sig[k,1]
        category <- as.character(enhancer_DHS_sig[k,3])
        
        if (length(strsplit(agds.file.name,"chr")[[1]]) == 2) {
          gds.path <- paste0(strsplit(agds.file.name,"chr")[[1]][1],"chr",chr,strsplit(agds.file.name,"chr")[[1]][2],".gds")
        } else {
          gds.path <- paste0(strsplit(agds.file.name,"chr")[[1]][1],"chr",chr,".gds")
        }
        genofile <- seqOpen(gds.path)
        
        res_cond <- Gene_Centric_Noncoding_cond(chr=chr,gene_name=gene_name,category=category,genofile=genofile,obj_nullmodel=nullobj,known_loci=known_loci_genome,rare_maf_cutoff=max.maf,
                                                QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                                                Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                                                Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name)
        enhancer_DHS_sig_cond <- rbind(enhancer_DHS_sig_cond,res_cond)
        
        seqClose(genofile)
      }
    }
    
    ### noncoding
    noncoding_sig_cond <- rbind(UTR_sig_cond[,c("Gene name","Chr","Category","#SNV","SKAT(1,25)","Burden(1,1)","ACAT-V(1,25)","STAAR-O")],
                                upstream_sig_cond[,c("Gene name","Chr","Category","#SNV","SKAT(1,25)","Burden(1,1)","ACAT-V(1,25)","STAAR-O")]) 
    noncoding_sig_cond <- rbind(noncoding_sig_cond,downstream_sig_cond[,c("Gene name","Chr","Category","#SNV","SKAT(1,25)","Burden(1,1)","ACAT-V(1,25)","STAAR-O")])
    noncoding_sig_cond <- rbind(noncoding_sig_cond,promoter_CAGE_sig_cond[,c("Gene name","Chr","Category","#SNV","SKAT(1,25)","Burden(1,1)","ACAT-V(1,25)","STAAR-O")])
    noncoding_sig_cond <- rbind(noncoding_sig_cond,promoter_DHS_sig_cond[,c("Gene name","Chr","Category","#SNV","SKAT(1,25)","Burden(1,1)","ACAT-V(1,25)","STAAR-O")])
    noncoding_sig_cond <- rbind(noncoding_sig_cond,enhancer_CAGE_sig_cond[,c("Gene name","Chr","Category","#SNV","SKAT(1,25)","Burden(1,1)","ACAT-V(1,25)","STAAR-O")])
    noncoding_sig_cond <- rbind(noncoding_sig_cond,enhancer_DHS_sig_cond[,c("Gene name","Chr","Category","#SNV","SKAT(1,25)","Burden(1,1)","ACAT-V(1,25)","STAAR-O")])
    
    write.csv(UTR_sig_cond, file = paste0(outfile, "_results_UTR_genome_sig_cond.csv"), row.names = FALSE)
    write.csv(upstream_sig_cond, file = paste0(outfile, "_results_upstream_genome_sig_cond.csv"), row.names = FALSE)
    write.csv(downstream_sig_cond, file = paste0(outfile, "_results_downstream_genome_sig_cond.csv"), row.names = FALSE)
    write.csv(promoter_CAGE_sig_cond, file = paste0(outfile, "_results_promoter_CAGE_genome_sig_cond.csv"), row.names = FALSE)
    write.csv(promoter_DHS_sig_cond, file = paste0(outfile, "_results_promoter_DHS_genome_sig_cond.csv"), row.names = FALSE)
    write.csv(enhancer_CAGE_sig_cond, file = paste0(outfile, "_results_enhancer_CAGE_genome_sig_cond.csv"), row.names = FALSE)
    write.csv(enhancer_DHS_sig_cond, file = paste0(outfile, "_results_enhancer_DHS_genome_sig_cond.csv"), row.names = FALSE)
    write.csv(noncoding_sig_cond, file = paste0(outfile, "_results_noncoding_genome_sig_cond.csv"), row.names = FALSE)
  }
}


if(test.type == "ncRNA") {
  cat("Summarizing ncRNA test results...\n")
  
  results_ncRNA_genome <- NULL
  for (chr in 1:22){
    if (file.exists(paste0(infile.prefix, chr, ".Rdata"))){
      cat("Chromosome:", chr, "\n")
      results <- get(load(paste0(infile.prefix, chr, ".Rdata")))
      results_ncRNA_genome <- c(results_ncRNA_genome, results)
    }
  }
  results_ncRNA_genome <- do.call("rbind",results_ncRNA_genome)
  ncRNA_sig <- results_ncRNA_genome[results_ncRNA_genome[,"STAAR-O"] < 2.5E-6, , drop = FALSE]
  save(results_ncRNA_genome, file = paste0(outfile, "_results_ncRNA_genome.Rdata"))
  write.csv(ncRNA_sig, file = paste0(outfile, "_results_ncRNA_genome_sig.csv"), row.names = FALSE)
  rm(results_ncRNA_genome)
  gc()
  
  if(nullobj.file == "NO_NULL_OBJ") {
    cat("No null model object, only summarize unconditional analysis results...\n")
  } else {
    if(annotation.name.catalog.file == "NO_ANNOTATION_NAME_CATALOG_FILE") stop("Error: Annotation name catalog file cannot be missing when test.type is ncRNA")
    Annotation_name_catalog <- read.csv(annotation.name.catalog.file, as.is=T)
    cat("Performing conditional analysis for ncRNA test results...\n")
    ### ncRNA
    ncRNA_sig_cond <- c()
    if(length(ncRNA_sig)!=0)
    {
      if((class(ncRNA_sig)!="matrix")&(class(ncRNA_sig)!="data.frame"))
      {
        ncRNA_sig <- matrix(ncRNA_sig,nrow=1)
      }
      
      for(k in 1:dim(ncRNA_sig)[1])
      {
        chr <- ncRNA_sig[k,2]
        gene_name <- ncRNA_sig[k,1]
        category <- as.character(ncRNA_sig[k,3])
        
        if (length(strsplit(agds.file.name,"chr")[[1]]) == 2) {
          gds.path <- paste0(strsplit(agds.file.name,"chr")[[1]][1],"chr",chr,strsplit(agds.file.name,"chr")[[1]][2],".gds")
        } else {
          gds.path <- paste0(strsplit(agds.file.name,"chr")[[1]][1],"chr",chr,".gds")
        }
        genofile <- seqOpen(gds.path)
        
        res_cond <- ncRNA_cond(chr=chr,gene_name=gene_name,genofile=genofile,obj_nullmodel=nullobj,known_loci=known_loci_genome,rare_maf_cutoff=max.maf,
                               QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                               Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                               Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name)
        ncRNA_sig_cond <- rbind(ncRNA_sig_cond,res_cond)
        
        seqClose(genofile)
      }
    }
    
    write.csv(ncRNA_sig_cond, file = paste0(outfile, "_results_ncRNA_genome_sig_cond.csv"), row.names = FALSE)
  }
}


if(test.type == "Sliding_Window") {
  cat("Summarizing sliding window test results...\n")

  results_sliding_window_genome <- NULL
  for (chr in 1:22){
    if (file.exists(paste0(infile.prefix, chr, ".Rdata"))){
      cat("Chromosome:", chr, "\n")
      results <- get(load(paste0(infile.prefix, chr, ".Rdata")))
      results_sliding_window_genome <- c(results_sliding_window_genome, results)
    }
  }
  results_sliding_window_genome <- results_sliding_window_genome[sapply(results_sliding_window_genome, function(x) class(x) != "try-error")]
  results_sliding_window_genome <- do.call("rbind",results_sliding_window_genome)
  sliding_window_sig <- results_sliding_window_genome[results_sliding_window_genome[,"STAAR-O"] < 0.05/dim(results_sliding_window_genome)[1], , drop = FALSE]
  save(results_sliding_window_genome, file = paste0(outfile, "_results_sliding_window_genome.Rdata"))
  write.csv(sliding_window_sig, file = paste0(outfile, "_results_sliding_window_genome_sig.csv"), row.names = FALSE)
  rm(results_sliding_window_genome)
  gc()
  
  if(nullobj.file == "NO_NULL_OBJ") {
    cat("No null model object, only summarize unconditional analysis results...\n")
  } else {
    if(annotation.name.catalog.file == "NO_ANNOTATION_NAME_CATALOG_FILE" && Use_annotation_weights) stop("Error: Annotation name catalog file cannot be missing when test.type is Sliding_Window and Use_annotation_weights is YES")
    Annotation_name_catalog <- read.csv(annotation.name.catalog.file, as.is=T)
    cat("Performing conditional analysis for sliding window test results...\n")
    sliding_window_sig_cond <- NULL
    if(length(sliding_window_sig)!=0)
    {
      for(kk in 1:dim(sliding_window_sig)[1])
      {
        chr <- as.numeric(sliding_window_sig[kk,1])
        start_loc <- as.numeric(sliding_window_sig[kk,2])
        end_loc <- as.numeric(sliding_window_sig[kk,3])
        
        if (length(strsplit(agds.file.name,"chr")[[1]]) == 2) {
          gds.path <- paste0(strsplit(agds.file.name,"chr")[[1]][1],"chr",chr,strsplit(agds.file.name,"chr")[[1]][2],".gds")
        } else {
          gds.path <- paste0(strsplit(agds.file.name,"chr")[[1]][1],"chr",chr,".gds")
        }
        genofile <- seqOpen(gds.path)
        
        res_cond <- Sliding_Window_cond(chr=chr,genofile=genofile,obj_nullmodel=nullobj,start_loc=start_loc,end_loc=end_loc,known_loci=known_loci_genome,rare_maf_cutoff=max.maf,
                                        QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                                        Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                                        Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name)
        sliding_window_sig_cond <- rbind(sliding_window_sig_cond,res_cond)
        
        seqClose(genofile)
      }
    }
    
    write.csv(sliding_window_sig_cond, file = paste0(outfile, "_results_sliding_window_genome_sig_cond.csv"), row.names = FALSE)
  }
}

