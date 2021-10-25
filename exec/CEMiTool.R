#!/usr/bin/env Rscript

"CEMiTool - Co-Expression Modules identification Tool
Modified by Brett Chapman 23rd September 2021

Usage: cemitool.R EXPRSFILE  --output=<dir> [--sample-annot=<annot> --samples-column=<samplecol> --class-column=<classcol> --no-filter (--filter-pval=<p>|--ngenes=<ngenes>) --apply_vst --eps --network-type=<nettype> --tom-type=<tomtype> --interactions=<inter> --pathways=<gmt> --ora-pvalue=<p> --cor-method=<cor> --cor-function=<corfunc> --no-merge --rank-method=<rank> --min-module-size=<min> --diss-thresh=<thresh> --center-func=<fun> --directed --top_hubs=<N> --top_hubs_interact=<N> --verbose]

Input:
  EXPRSFILE                         a normalized expression file .tsv format

Options:
  -h --help                          show this help message
  --version                          show program version
  -s <annot> --sample-annot=<annot>  sample annotation, must have a column with sample names and class
  --samples-column=<samplecol>       the column name containing sample names in template file [default: SampleName]
  --class-column=<classcol>          the column name containing classes in template file [default: Class]
  -i <int> --interactions=<int>      gene interaction file, must have two columns. Interaction file produced from the top 100 hub genes is used as default
  -p <gmt> --pathways=<gmt>          GMT file name (Gene Matrix Transposed format)
  --ora-pvalue=<p>                   p-value cutoff to be used on over representation analysis [default: 0.05]
  --no-filter                        does not filter the expression data.frame
  --filter-pval=<p>                  p-value to be used in the filtering step [default: 0.1]
  --apply_vst                        apply Variance Stabilizing Transformation
  --ngenes=<ngenes>                  number of genes remaining after filtering
  --eps=<eps>                        epsilon [default: 0.1]
  -c <cor> --cor-method=<cor>        correlation method (spearman or pearson) [default: pearson]
  --cor-function=<corfunc>           A character string indicating the correlation function to be used. Supported functions are 'cor' and 'bicor' [default: cor]
  --network-type=<nettype>           network type, 'signed' or 'unsigned' [default: unsigned]
  --tom-type=<nettype>               TOM type, 'signed' or 'unsigned' [default: signed]
  --no-merge                         does not merge related modules based on eigengene similarity
  --rank-method=<rank>               rank method [default: mean]
  --min-module-size=<min>            minimum module size [default: 30]
  --diss-thresh=<thresh>             module merging correlation threshold for eigengene similarity [default: 0.8]
  --center-func=<fun>                metric used for centering [default: mean]
  --directed                         the interactions are directed
  --top_hubs=<N>		     The top N hub genes to list per module and to output as an expression matrix [default: 10] 
  --top_hubs_interact=<N>	     The top N hub genes to use for the basis of the interaction network [default: 10]
  -o <dir> --output=<dir>            output directory
  -v --verbose                       display progress messages

Authors:
  Pedro S T Russo - pedro.russo at usp.br
  Gustavo R Ferreira - gustavo.rodrigues.ferreira at usp.br
  Lucas E Cardozo - lucasecardozo at usp.br
  Matheus C Burger - burger at usp.br
  Thiago D C Hirata - thiagodch at gmail.com
  Diogenes S Lima - diogenes.lima at usp.br
  Fernando M Passos - fmarcon at usp.br
  Raul A Carrasco 
  Melissa Lever - melissalever at gmail.com 
  Vinicius Maracaja-Coutinho
  Helder I Nakaya - hnakaya at usp.br

More information:
  www.csbiology.com
  Computational Systems Biology Laboratory
  University of Sao Paulo, Brazil
" -> doc

if (!interactive()) {
    # Get and check arguments.
    suppressMessages(library("docopt"))
    arg <- docopt(doc, version="1.17.0\n", strict=TRUE)
    arg <- arg[!sapply(arg, is.null)][-(1:2)]  # filter missing, 'help' and 'version'
    clean <- function(s) gsub('-', '_', gsub('^-+', '', tolower(s)))
    names(arg) <- clean(names(arg))
    parameters <- arg

    print(parameters)
    
    library(parallel)

    library(ggplot2)

    library(dplyr)

    ## RUN
    library("CEMiTool")

    # parameters
    p <- list()

    # verbosity
    p$verbose <- parameters[["verbose"]]

    # expression file
    if(p$verbose){
        message("Reading expression file ...")
    }
    p$expr <- data.table::fread(parameters[["exprsfile"]], data.table=FALSE)
    
    # remove the column containing gene symbols
    if(p$verbose) {
        message("Setting row names ...")
    }
    rownames(p$expr) <- p$expr[,1]
    p$expr[,1] <- NULL

    # verify if all columns are numeric
    if(!all(sapply(p$expr, is.numeric))){
        stop("Please make sure that your expression file have only numeric values.")
    }

    # sample annotation file
    if("sample_annot" %in% names(parameters)){
        if(p$verbose){
            message("Reading sample annotation file...")
        }
        p$annot <- data.table::fread(parameters[["sample_annot"]], data.table=FALSE)

        # sample name column in sample annotation file
        p$sample_name_column <- parameters[["samples_column"]]

        # class column in sample annotation file
        p$class_column <- parameters[["class_column"]] 
    }

    # Should filter ?
    p$filter <- !parameters[["no_filter"]]

    # filter p-value or number of genes after filtering
    if("ngenes" %in% names(parameters)) {
        p$n_genes <- as.numeric(parameters[["ngenes"]])
    } else {
        p$filter_pval <- as.numeric(parameters[["filter_pval"]])
    }

    p$apply_vst <- parameters[["apply_vst"]]

    # gmt list
    if("pathways" %in% names(parameters)){
        if(p$verbose){
            message("Reading GMT file ...")
        }
        p$gmt <- read_gmt(parameters[["pathways"]])
    }

    # interactions file
    hub_interact <- FALSE
    if("interactions" %in% names(parameters)){
        if(p$verbose){
            message("Reading interactions file ...")
        }
        p$interactions <- data.table::fread(parameters[["interactions"]], data.table=FALSE)
    } else {
    	hub_interact <- TRUE
    }

    # correlation method
    p$cor_method <- parameters[["cor_method"]]

    # correlation function
    p$cor_function <- parameters[["cor_function"]]

    # epsilon
    p$eps <- parameters[["eps"]]
      
    # network type
    p$network_type <- parameters[["network_type"]]
      
    # TOM type
    p$tom_type <- parameters[["tom_type"]]
      
    # should similar modules be merged ?
    p$merge_similar <- !parameters[["no_merge"]]

    # p-value cutoff to be used on over representation analysis
    p$ora_pval <- as.numeric(parameters[["ora_pvalue"]])

    # minimum size of a module
    p$min_ngen <- as.numeric(parameters[["min_module_size"]])

    # dissimilarity threshold
    p$diss_thresh <- as.numeric(parameters[["diss_thresh"]])

    p$rank_method <- parameters[["rank_method"]]

    p$center_func <- parameters[["center_func"]]
    # should create figures ?
    p$plot <- TRUE

    # Is the interactions file directed ?
    p$directed <- parameters[["directed"]]

    #The top N hub genes to output
    if("top_hubs" %in% names(parameters)){
        if(p$verbose){
            message("Reading top hubs ...")
        }
        top_hubs <- as.numeric(parameters[["top_hubs"]])
    } else {
    	top_hubs <- as.numeric(10)
    }

    #The top N hub interaction genes to output
    if("top_hubs_interact" %in% names(parameters)){
        if(p$verbose){
            message("Reading top hubs interaction ...")
        }
        top_hubs_interact <- as.numeric(parameters[["top_hubs_interact"]])
    } else {
        top_hubs_interact <- as.numeric(10)
    }

    # Increase gene label overlaps to infinity
    options(ggrepel.max.overlaps = Inf)

    # CEMiTool
    if(p$verbose){
        message("Running CEMiTool ...")
    }

    cbindPad <- function(...){
    args <- list(...)
    n <- sapply(args,nrow)
    mx <- max(n)
    pad <- function(x, mx){
    	if (nrow(x) < mx){
        	nms <- colnames(x)
        	padTemp <- matrix(NA, mx - nrow(x), ncol(x))
        	colnames(padTemp) <- nms
        	if (ncol(x)==0) {
        		return(padTemp)
        	} else {
        		return(rbind(x,padTemp))
          	}
    	}
    	else{
       		return(x)
    	}
    }		
    rs <- lapply(args,pad,mx)
    return(do.call(cbind,rs))
    }

    # If true run a first and then second pass of cemitool using the top N hub genes to construct an interaction network
    if(hub_interact == TRUE){
	    cem <- do.call(cemitool, p)
	    cem <- find_modules(cem)

	    # Get top N hub genes to use as the interaction network
	    network <- get_hubs(cem,top_hubs_interact)

	    for(i in 1:length(network)){
            	if(i == 1){
			hub_network <- na.omit(as.data.frame(names(network[[i]])))
			hub_network <- hub_network %>% slice_head(n = top_hubs_interact)
            	} else {
			xx <- na.omit(as.data.frame(names(network[[i]])))
			xx <- xx %>% slice_head(n = top_hubs_interact)
    			hub_network <- cbindPad(hub_network, xx)
            	}
	    }
	    
	    colnames(hub_network) <- names(network)

	    interactions <- list()
	    for(i in 1:length(hub_network)){
	    	y <- hub_network[[i]]
	    	for(n in 1:length(y)){
		    	yy <- y[[n]]
			if(!is.na(yy)){
		    		for(nn in 1:length(y)){
			    		xx <- y[[nn]]
			    		if(!is.na(xx)){
						if(xx != yy){
				    			if(length(interactions) == 0){
					    			interactions[['Gene1']] <- na.omit(as.data.frame(yy))
					    			interactions[['Gene2']] <- na.omit(as.data.frame(xx))
				    			} else {
					    			interactions[['Gene1']] <- rbind(interactions$Gene1, na.omit(as.data.frame(yy)))
								interactions[['Gene2']] <- rbind(interactions$Gene2, na.omit(as.data.frame(xx)))
							}
						}
					}
				}
			}
		}
	    }
	    
	    hub_interactions <- do.call(cbind, interactions)
	    colnames(hub_interactions) <- c('Gene1', 'Gene2')

	    p$interactions <- hub_interactions

	    cem <- do.call(cemitool, p)

	    #Get top N hub genes
	    cem <- find_modules(cem)
	    x <- get_hubs(cem,top_hubs)

	    for(i in 1:length(x)){
		    if(i == 1){
			hubs <- na.omit(as.data.frame(names(x[[i]])))
			hubs <- hubs %>% slice_head(n = top_hubs)
		    } else {
			xx <- na.omit(as.data.frame(names(x[[i]])))
			xx <- xx %>% slice_head(n = top_hubs)
			hubs <- cbindPad(hubs, xx)
		    }
	    }
	    colnames(hubs) <- names(x)
	    top_N_hub_list <- na.omit(unlist(hubs))

	    # extract out expression data, transform and slice out the top hub genes
	    expr_transformed <- filter_genes(cem@expression, apply_vst = cem@parameters$apply_vst)
	    selected_expr_transformed <- expr_transformed[select_genes(expr_transformed, filter_pval = cem@parameters$filter_pval), ]
	    top_N_expr <- selected_expr_transformed[top_N_hub_list, ]
	    top_N_expr <- cbind(gene_name = rownames(top_N_expr),top_N_expr)
            rownames(top_N_expr) <- NULL
	    expr_transformed <- cbind(gene_name = rownames(expr_transformed),expr_transformed)
	    rownames(expr_transformed) <- NULL
	    selected_expr_transformed <- cbind(gene_name = rownames(selected_expr_transformed),selected_expr_transformed)
            rownames(selected_expr_transformed) <- NULL

	    # Write out transformed expression data, filtered data and top N hub gene expression data to CSV
	    dir.create(path = paste(parameters[["output"]], "Tables", sep = '/'), recursive = TRUE)
	    write.table(expr_transformed, paste(parameters[["output"]], "Tables", "transformed_expression.tsv", sep = '/'), row.names = FALSE, sep="\t")
	    write.table(selected_expr_transformed, paste(parameters[["output"]], "Tables", paste("transformed_filtered-",cem@parameters$filter_pval,"_expression.tsv", sep = ''), sep = '/'), row.names = FALSE, sep="\t")
	    write.table(top_N_expr, paste(parameters[["output"]], "Tables", paste("top_",top_hubs,"_transformed_filtered-",cem@parameters$filter_pval,"_expression.tsv", sep = ''), sep = '/'), row.names = FALSE, sep="\t")
    } else {
	    cem <- do.call(cemitool, p)
            cem <- find_modules(cem)
	    x <- get_hubs(cem,top_hubs)
	    
	    for(i in 1:length(x)){
		    if(i == 1){
                    	hubs <- na.omit(as.data.frame(names(x[[i]])))
			hubs <- hubs %>% slice_head(n = top_hubs)
		    } else {
			xx <- na.omit(as.data.frame(names(x[[i]])))
			xx <- xx %>% slice_head(n = top_hubs)
			hubs <- cbindPad(hubs, xx)
		    }
	    }
	    colnames(hubs) <- names(x)
	    top_N_hub_list <- na.omit(unlist(hubs))

	    # extract out expression data, transform and slice out the top hub genes
            expr_transformed <- filter_genes(cem@expression, apply_vst = cem@parameters$apply_vst)
            selected_expr_transformed <- expr_transformed[select_genes(expr_transformed, filter_pval = cem@parameters$filter_pval), ]
            top_N_expr <- selected_expr_transformed[top_N_hub_list, ]
	    top_N_expr <- cbind(gene_name = rownames(top_N_expr),top_N_expr)
	    rownames(top_N_expr) <- NULL
	    expr_transformed <- cbind(gene_name = rownames(expr_transformed),expr_transformed)
            rownames(expr_transformed) <- NULL
            selected_expr_transformed <- cbind(gene_name = rownames(selected_expr_transformed),selected_expr_transformed)
            rownames(selected_expr_transformed) <- NULL

	    # Write out top N hub expression data to CSV
            dir.create(path = paste(parameters[["output"]], "Tables", sep = '/'), recursive = TRUE)
	    write.table(expr_transformed, paste(parameters[["output"]], "Tables", "transformed_expression.tsv", sep = '/'), row.names = FALSE, sep="\t")
            write.table(selected_expr_transformed, paste(parameters[["output"]], "Tables", paste("transformed_filtered-",cem@parameters$filter_pval,"_expression.tsv", sep = ''), sep = '/'), row.names = FALSE, sep="\t")
            write.table(top_N_expr, paste(parameters[["output"]], "Tables", paste("top_",top_hubs,"_transformed_filtered-",cem@parameters$filter_pval,"_expression.tsv", sep = ''), sep = '/'), row.names = FALSE, sep="\t")
    }

    # Write out top N hub genes to CSV
    dir.create(path = paste(parameters[["output"]], "Tables", sep = '/'), recursive = TRUE)
    write.table(hubs, paste(parameters[["output"]], "Tables", paste("top_",top_hubs,"_hub_genes_per_module.tsv", sep = ''), sep = '/'), row.names = FALSE, sep="\t")

    # Save files
    if(p$verbose){
        message("Writing CEMiTool results...")
    }

    # Write out interactions file if based on top hub genes
    if(hub_interact == TRUE){
    	write.table(p$interactions, paste(parameters[["output"]], "Tables", "hub_gene_interactions.tsv", sep = '/'), row.names = FALSE, quote = FALSE, sep="\t")
    }

    # Write out original gene set gmt file
    write.table(data.table::fread(parameters[["pathways"]], data.table=FALSE), paste(parameters[["output"]], "Tables", "gene_sets.gmt", sep = '/'), row.names = FALSE, col.names = FALSE, quote = FALSE, sep="\t")

    # Generate reports
    generate_report(cem, directory=(paste(parameters[["output"]], "Report", sep = '/')), force=TRUE)
   
    # Write analysis results
    write_files(cem, directory=(paste(parameters[["output"]], "Tables", sep = '/')), force=TRUE)

    # save all plots
    save_plots(cem, "all", directory=(paste(parameters[["output"]], "Plots", sep = '/')), force=TRUE)
}
