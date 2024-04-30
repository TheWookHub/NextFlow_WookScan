#!/usr/bin/env Rscript

## Libraries required to run

suppressWarnings(
    suppressPackageStartupMessages(
        {
            library(optparse)
            library(tidyverse)
            library(data.table)
            library(igraph)
            library(doParallel)
            library(foreach)
        }   
    )
)


############################
# Global Variables for use #
############################

optList = c(
    "case_path","threshold","dict_path","total_path",
    "pairwise_path","blast_path","out_path","out_name",
    "avarda_names","show_pb","cores"
)

NUMCORE = 4

######################################################
# Defining arguments to be passed into AVARDA script #
######################################################

option_list = list(
    # Will probably use this option later. But First make the params work
    # make_option(c("-p", "--param_file"), action="store", default=NULL, type='character',
    #             help="parameter file to configure SKAT-O run. Must be in comma separated format."),  
    make_option(
        c("--case_path"),
        action = "store",
        help = paste(
            "Path to file of data to be analyzed, rows are peptides columns are samples,",
            "values are generally binary values (indicating hits)."
        )
    ),
    make_option(
        c("--threshold"),
        action="store",
        help = "Threshold value for input; set to 1 if using binary matrix."
    ),
    make_option(
        c("--dict_path"),
        action="store",
        help = "A csv file indicating all peptide-peptide alignments."
    ),
    make_option(
        c("--total_path"),
        action="store",
        help = paste(
            "Path to a pre-made table (.csv) listing each viruses abundance in the viral",
            "peptide library."
        )
    ),
    make_option(
        c("--pairwise_path"),
        action="store",
        help = paste(
            "Path to a pre-made table (.csv) listing each viruses relative representation",
            "compared with other viruses in your viral peptide library."
        )
    ),
    make_option(
        c("--blast_path"),
        action="store",
        help = "Path to pre-made table (.csv) of peptide-virus alignments storing bitscores."
    ),
    make_option(
        c("--out_path"),
        action="store",
        help = "Path to where you want to write the results to."
    ),
    make_option(
        c("--out_name"),
        action="store",
        help = "File prefix to be added on to file names."
    ),
    make_option(
        c("--avarda_names"),
        action="store",
        help = paste(
            "Table (.csv) with 2 columns. First column stores user peptide id (matching with file in case_path).",
            "Second column are the peptide id from the library."
        )
    ),
    make_option(
        c("--show_pb"),
        action = "store_true",
        default = FALSE,
        help = "Flag to show printing progress bar."
    ),
    make_option(
        c("--cores"),
        action = "store",
        default = 4,
        help = "Number of cores to use. [Default = 4]",
        type = "integer"
    )    
)

help_msg = paste(
    "AVARDA Tool - implemented by Monaco et al.",
    "and modified by Preston L for nextflow usage."
)

opt = parse_args(
    OptionParser(
        option_list = option_list,
        description = help_msg
    )
)

######################
# Define function(s) #
######################

helpMsg = function(input_param){
    print_help(
        OptionParser(
            option_list=option_list,
            description = help_msg
        )
    )
  message("-------------------------------------------")

  missing = setdiff(optList, names(input_param))
  message(paste("The following required parameters has not defined:",paste(missing, collapse = ', ')))
  message("Exit.")
}

AVARDA = function(case_path,thresh,dict_path,total_path,pairwise_path,blast_path,out_path,out_name,avarda_names,showPB){    
    # read in peptide-peptide dictionary
    dict =  data.frame(fread(dict_path,data.table = FALSE))
    # read in null probability table 
    total = data.frame(fread(total_path,data.table = FALSE))
    # read in null probability table
    pairwise = data.frame(fread(pairwise_path,data.table = FALSE),row.names=1)
    # read in virus-peptide alignment filtered table
    blast = data.frame(fread(blast_path,data.table = FALSE),row.names=1)
    # this is the case of interest to be analyzed - changeable
    case = data.frame(fread(case_path, data.table = FALSE,check.names = FALSE, header = TRUE))
    avarda_names =  data.frame(fread(avarda_names, data.table = FALSE,check.names = FALSE))
    # merge case and avarda_names together so that we can replace u_pep_id (WOOKSCAN_000_ATCT.. etc) 
    # with pep_id (1,2,3 etc). This ends up with pep_id a virus to pep_id table
    case[,1] = case[,1:2] %>% as.data.frame() %>%  full_join(avarda_names) %>% select(pep_id)
    
    # --- Filter function --- #
    # independence filter that takes a dictionary (edges) and 
    # a set of pep_ids (vertex) and tells the minimal number of unique epitopes.
    # Reminder dict came from a blastp pep2pep alignment including only eval
    # <100. Indicating that these peptides are overlapping / similar to each other
    # MODIFICATION: 
    #   -> Make it return peptide ids too.    
    filter  = function(edge,vertex){
        # this is a vertex list of peptide ids for a specific virus in one
        # patient
        nodes = unlist(vertex)
        # grab edge-pairs where the 'from' match with the vertex list 
        links_filtered = subset(edge,unlist(edge[,1]) %in% nodes)
        # grab edge-pairs where the 'to' match with the vertex list from links_filtered
        # the aim is to find the peptides ids from nodes that link to each other
        # ie - they're overlapping or similar.
        links_filtered = subset(links_filtered,links_filtered[,2] %in% nodes)
        # if we do find peptides that link to each other..
        if(dim(links_filtered)[1]!=0){
            # create an undirected network
            # imagine a network with N nodes and only a subset of
            # the nodes have edges
            net = as.undirected(
                graph_from_data_frame(
                    d = links_filtered,
                    vertices = nodes, 
                    directed = F
                ) 
            )
            
            # if it creates disconnected islands of network, retrieve them
            # via "decomposing"
            # x = decompose.graph(net) # deprecated
            x = decompose(net)
            
            # get island networks that have node count < 30
            x_1 = x[sapply(x,vcount)<30]
            
            # get the largest (maximum) independence set size
            # x_1_sum  = sum(unlist(lapply(x_1,independence.number))) # deprecated
            x_1_sum  = sum(unlist(lapply(x_1,ivs_size)))            
            
            # here we retrieve the ids of the largest independent vertex set names
            # (i.e. the peptide ids). We take the ids of the very first set (if
            # there are multiple largest sets) because it doesn't really matter
            # which set we take. We could take a random set everytime. Since at the 
            # largest sets, you're basically swapping out one for another, meaning the
            # one being swapped out cannot be included together with the one swapped in
            # due to them being linked to each other (and therefore are overlapping / 
            # similar peptides).
            x_1_ids = unlist(
                lapply(
                    x_1,
                    function(i){
                        all_largest_sets = largest_ivs(i)
                        return(all_largest_sets[[1]]$name)
                    }
                )                
            )
            if(x_1_sum != length(x_1_ids)){
                message("Something weird has happened to getting independent vertex set ids and size.")
            }
            # get island networks that have a node count >= 30
            # usually means loads of peptides are overlapping each other
            # need some simplifications
            x_2 = x[sapply(x,vcount)>=30]
            temp = c()
            x_2_ids = c()
            #x_2 = x
            if(length(x_2) >0){
                # for each island network in the  >= 30 nodes set
                for(R in 1:length(x_2)){                    
                    x_2_r = x_2[[R]]
                    # Continuously checking if the most connected node has degree > 5
                    while(max(degree(x_2_r)>5)){            
                        # "initially remove the most interconnected peptides iteratively 
                        # until the most interconnected peptide(s) has 5 alignments" 
                        # - From the paper
                        toss = degree(x_2_r)==max(degree(x_2_r))
                        x_2_r = delete_vertices(x_2_r, V(x_2_r)[toss][1])
                    }
                    # once again, break it down to see islands of disconnected network
                    x_l = decompose(x_2_r)
                    # for each 
                    # temp[R] = sum(unlist(lapply(x_l,independence.number))) # deprecated
                    temp[R] = sum(unlist(lapply(x_l,ivs_size)))
                    x_2_ids[R] = list(
                        unlist(
                            lapply(
                                x_l,
                                function(i){
                                    all_largest_sets = largest_ivs(i)
                                    return(all_largest_sets[[1]]$name)
                                }
                            )
                        )
                    )
                }
            }
            total_unique_peptides = sum(x_1_sum)+sum(temp)
            all_unique_peptide_ids = list(c(x_1_ids,unlist(x_2_ids)))
            if(length(all_unique_peptide_ids) < 1){
                all_unique_peptide_ids = list(c(""))
            }
            # return(sum(x_1_sum)+sum(temp))
            return(c(total_unique_peptides, all_unique_peptide_ids))
        }
        # if we do not find ANY peptides that link to each other..
        # means peptide (node) is their own network island. 
        # and each is a maximum independent vertex set. Hence,
        # returning the length of the nodes and the node id (peptide id)
        # is good
        if(dim(links_filtered)[1]==0){
            return(c(length(nodes),list(c(""))))
        }
    }

    # --- Binom function --- #
    # v_i       ->  all peptides aligned to virus_i (had bitscore >0)
    # v_xr      ->  all peptides that are cross reactive (had bitscore >0 but <80 to virus_i)
    # v_i_j     ->
    # N_rank    ->  a list of all peptides that are enriched in a patient. And these peptides
    # null_prob ->  the probability (in total_probability) of getting virus_i in the given viral pool
    binom_test  = function(v_i,v_xr,v_i_j,N_rank,null_prob){
         # all virus i aligning minus the xr
        v_total = v_i[!v_i %in% v_xr]
        # all virus i evidence minus any shared with virus j
        v_total = v_total[!v_total %in% v_i_j]
        # this calculates N_rank
        unique_peptides = filter(dict,v_total)
        v_total_f = unique_peptides[[1]] # count
        v_total_ids = unique_peptides[[2]] # ids
        #if(length(N_rank)!=length(unlist(N_rank_2))){
        N_rank = N_rank[!N_rank %in% v_xr]
        # this is by default zero for total binom calculation
        N_rank = N_rank[!N_rank %in% (v_i %in% v_i_j)]
        unique_n_ranks = filter(dict,N_rank)
        N_rank_f = unique_n_ranks[[1]]
        N_rank_ids = unique_n_ranks[[2]]
        #}
        if(N_rank_f == 0){
            return(NULL)
        }
        # x gives the p-value
        x = binom.test(v_total_f,N_rank_f,unlist(null_prob),"greater")[[3]]
        output = list(x,v_total_f,N_rank_f,v_total_ids)
        return(output)
    }
  
    # --- Total calc function --- #
    # sub-setting cases and also ranks peptides
    total_calc  = function(case,column,thresh,total,blast){ 
        # Processing of the blast virus matrix on a case by case 
        total_probs = total
        # retrieve one column of data (1 patient) and make peptide ids (1,2,3 etc..) as rownames
        enriched = data.frame(case[,c(1,column+1)],row.names = 1) 
        # subset on enrichment matrix peptides that are greater than user defined threshold
        enriched = subset(enriched,enriched >=thresh)
        
        # previous code looks for an "enrich" with length > 0. enrich is defined as thresh in the super
        # function. Since Thresh is a numeric,the length will always bit > 0. Hence it should
        # be a typo. The if case should refer to 'enriched' dataframe rather than 'enrich' 
        # the threshold value.
        # if(length(enrich>0)){
        # As long as dataframe has 1 column
        if(length(enriched>0)){
            # if data frame has more than 1600 rows
            # NOTE: unsure why original authors put 1600 as a factor for limiting 
            # number of peptides?
            if(dim(enriched)[1]>1600){
                # reverse the order and keep dataframe structure because the dataframe
                # at this point is only 1 column                
                # enriched = enriched[order(-enriched),1,drop=FALSE] # ORIGINAL
                enriched = enriched[order(-enriched[,1]),1,drop=FALSE] # MODIFIED
                # Then take the top 1600 as enriched.
                enriched = enriched[1:1600,1,drop = FALSE] 
            }
            # if data frame has more than 1 row
            if(dim(enriched)[1] > 0){
                # here we take enriched peptides and then check against which viruses
                # have been picked up via the pep_id
                blast_subset = subset(blast,row.names(blast) %in% row.names(enriched))
                # then we check if the pep_ids have >80 bitscore to that virus. If a virus
                # has more than 2 counts where bitscore >80, then we take the column 
                # (each column is a virus, rows depict peptides, value is bitscore)
                blast_subset=blast_subset[, colSums(ifelse(blast_subset>80, 1, 0)) > 2]
                # Sometimes it might not be a dataframe....? E.g. vector,then it fails.
                if(is.data.frame(blast_subset)==TRUE){
                    # there are at least 1 virus detected
                    if(dim(blast_subset)[2]>0){
                        fullmatrix_sorted = blast_subset
                        v_i_j = NULL
                        order = c()
                        virus = c()
                        if(showPB){
                            pb <- txtProgressBar(min = 0, max =dim(fullmatrix_sorted)[2], style = 3)
                        }                        
                        # For each virus (column) we go through it.
                        for(R in 1:dim(fullmatrix_sorted)[2]){
                            # get pep_ids that align crossreactively to virus_i                            
                            virus_i_xr = rownames(fullmatrix_sorted)[fullmatrix_sorted[,R]>0 & fullmatrix_sorted[,R]<80]
                            # get all peptides align to virus_i
                            v_i = rownames(fullmatrix_sorted)[fullmatrix_sorted[,R]>0]
                            # call the null probability for virus_i using total_pro
                            probability = total_probs[grep(paste0(colnames(fullmatrix_sorted[R]),"$",collapse = ""),unlist(total_probs[,1])),2]
                            # run binom function that takes into account cross-reactivity and shared peptides (which is 0 for this step)
                            x = binom_test(v_i,virus_i_xr,v_i_j,row.names(blast_subset),probability)
                            # reports the p-value for a viruses likelihood of infection
                            order[R] = x[1] 
                            # get which virus is being compared in interation
                            virus[R] = colnames(fullmatrix_sorted[R]) 
                            if(showPB){
                                setTxtProgressBar(pb, R)
                            }
                            
                        }
                        #reorder the subset virus-peptide blast matrix by likelihood of infection
                        results = fullmatrix_sorted[,order(unlist(order))]
                        return(results)
                    }
                    return(NULL)
                }
                return(NULL)
            }
        }
        return(NULL)
    }
    # --- pairwise calc function --- #
    pairwise_calc  = function(total,pairwise,rank1){
        ## just initializing some variables
        unique_probs =pairwise
        total_probs = total
        fullmatrix_sorted_ij = rank1
        final_matrix = data.frame(matrix(0,nrow = dim(fullmatrix_sorted_ij)[1],ncol =dim(fullmatrix_sorted_ij)[2])) #initialize the final matrix of peptide-virus alignments post reassignment
        N_rank = rownames(fullmatrix_sorted_ij)
        z = N_rank
        output = data.frame(matrix(ncol = 11, nrow =  dim(fullmatrix_sorted_ij)[2])) #initialize matrix for final data with significance data
        sim_tag = data.frame(matrix(ncol = 2)) #empty vector to fill with virus pairs that are sim-tagged
        x1 = 1
        ##
        if(showPB){
            pb <- txtProgressBar(min = 0, max = dim(fullmatrix_sorted_ij)[2], style = 3)
        }        
        for(R1 in 1: dim(fullmatrix_sorted_ij)[2]){ # test
            # for(R1 in 1: 13){ # test            
            virus_i = colnames(fullmatrix_sorted_ij)[R1] # get name of virus_i
            virus_i_hits = as.data.frame(fullmatrix_sorted_ij[R1]) #get vector of virus_i peptide alignments
            v_i = subset(virus_i_hits,virus_i_hits>0) # get index of all peptide alignments >0
            v_i_xr = subset(virus_i_hits,virus_i_hits>0 & virus_i_hits<80) # get index of all peptide alignments <80 and >0 (the xr)
            if(sum(virus_i_hits>=80)>=3){
                R2 = 1
                while(R2 <= dim(fullmatrix_sorted_ij)[2]){
                    if(R1 < R2){ #skip iterations of virus_i comparing to previously evaluated viruses
                        binary_i_j = ifelse(fullmatrix_sorted_ij>0, 1, 0)
                        shared = binary_i_j[,R1]*binary_i_j[,R2] # vector multiplication to get the peptides with alignments to virus_i and virus_j
                        virus_j_hits = as.data.frame(fullmatrix_sorted_ij[R2]) #get vector of virus_j alignments to enriched peptides
                        virus_j = colnames(fullmatrix_sorted_ij)[R2] #get name of virus_j
                        if(sum(shared) != 0){ #if both viruses are completely unique no reassignment can occur so skip this
                            probability_i_j = unique_probs[grep(paste0(virus_j,"$",collapse = ""),row.names(unique_probs)),grep(paste0(virus_i,"$",collapse = ""),row.names(unique_probs))] #probabiltiy that a random peptide will align exclusively to virus_i relative to virus_j
                            probability_j_i = unique_probs[grep(paste0(virus_i,"$",collapse = ""),row.names(unique_probs)),grep(paste0(virus_j,"$",collapse = ""),row.names(unique_probs))]
                            v_j = subset(virus_j_hits,virus_j_hits>0) #get all v_j alignments
                            v_j_xr = subset(virus_j_hits,virus_j_hits>0 & virus_j_hits<80) # get v_j xr alignments
                            p_i_j =  binom_test(rownames(v_i),rownames(v_i_xr),rownames(v_j),N_rank,probability_i_j) #likelihood that virus_i has an infection that is unique relative to virus_j
                            p_j_i =  binom_test(rownames(v_j),rownames(v_j_xr),rownames(v_i),N_rank,probability_j_i) #likelihood that virus_j has an infection that is unique relative to virus_i
                            
                            if(p_i_j[1] <= .05 & p_j_i[1] >.05 & p_i_j[2] >=3 ){ #if virus_j loses (and virus_i has at least 3 aligning peptideds
                                fullmatrix_sorted_ij[(shared==1),R2] = 0 #virus_j loses all alignments that are shared with virus_i (even if they are xr in v_i and evidence in v_j)                                
                            }
                            if(p_i_j[1] > .05 & p_j_i[1] <=.05 & p_j_i[2] >=3 ){ #if virus_i loses
                                fullmatrix_sorted_ij[(shared==1),R1] = 0                                
                            }
                            if(p_i_j[1] > .05 & p_j_i[1] > .05 | p_j_i[2] <3 & p_i_j[2] < 3 ){ #if both viruses fail to have enough unique evidence (but both have three peptides)
                                sim_tag[x1,] = cbind(virus_i,virus_j) #note that the pair is indistinguishable
                                #print(x1)
                                x1 = x1+1
                            }
                        }
                    }
                    R2= R2+1
                }                
                ###at this point all virus_j have been compared
                order = c() #initialize
                virus = c() #initialize
                if(R1<dim(fullmatrix_sorted_ij)[2]){
                    z =N_rank
                    virus_i_prob = binom_test(rownames(v_i),rownames(v_i_xr),NULL,z,total_probs[grep(paste0(colnames(fullmatrix_sorted_ij[R1]),"$",collapse = ""),unlist(total_probs[,1])),2])                
                    if(virus_i_prob[1]<=.05 & is.null(dim(fullmatrix_sorted_ij[which(fullmatrix_sorted_ij[R1]>0),(R1+1):dim(fullmatrix_sorted_ij)[2]]))==FALSE){ ##if virus_i was significant then remove all peptides that aligned exclusively to virus_i from future calculations
                        N_remove = names(which(apply(fullmatrix_sorted_ij[which(fullmatrix_sorted_ij[R1]>0),(R1+1):dim(fullmatrix_sorted_ij)[2]],1,max)==0))
                        z =N_rank[!N_rank %in% N_remove] 
                    }                
                    ### Need to rerank remaining viruses after the v_i_j comparisons
                    for(R in (R1+1):dim(fullmatrix_sorted_ij)[2]){ # go iteartively through viruses
                        virus_i_xr2 = rownames(fullmatrix_sorted_ij)[fullmatrix_sorted_ij[,R]>0 & fullmatrix_sorted_ij[,R]<80] # get all peptides total with alignments to virus x
                        v_i2 = rownames(fullmatrix_sorted_ij)[fullmatrix_sorted_ij[,R]>0] # get peptides that are evidence for virus x
                        probability = total_probs[grep(paste0(colnames(fullmatrix_sorted_ij[R]),"$",collapse = ""),unlist(total_probs[,1])),2] # call the null probability for virus x
                        
                        if(length(v_i2[!v_i2 %in% virus_i_xr2])[1] >0){ #if there is at least 1 evidence peptide for virus_x reorder
                            order[R-R1] = binom.test(length(v_i2[!v_i2 %in% virus_i_xr2])[1],length(z),as.numeric(probability),"greater")[[3]]
                        }
                        if(length(v_i2[!v_i2 %in% virus_i_xr2])[1] ==0){ #otherwise just skip
                            order[R-R1] = 1
                        }
                    }                
                    fullmatrix_sorted_ij = cbind(fullmatrix_sorted_ij[1:R1],fullmatrix_sorted_ij[order(order)+R1])
                    # print(dim(fullmatrix_sorted_ij))
                }
            }
            ## at this point all viruses are evaluated so just need final ranking step
            ## the below data is meta data we output when investigating data
            final = fullmatrix_sorted_ij[R1]
            final_i = subset(final,final>0)
            final_xr = subset(final,final>0 & final<80)
            probability_i = total_probs[grep(paste0(colnames(fullmatrix_sorted_ij[R1]),"$",collapse = ""),unlist(total_probs[,1])),2] # call the null probability for virus x
            final_rank = binom_test(rownames(final_i),rownames(final_xr),0,N_rank,probability_i)
            final_i_e = subset(final,final>=80)
            output[R1,] = c(
                virus_i, #1 - virus name
                final_rank[1], #2 - p-value
                paste(rownames(final_i_e),collapse = "|"),#3 - evidence peptide id
                paste(rownames(final_xr),collapse = "|"),#4 - cross reactive peptide id
                length(N_rank),#5 N_rank - Num of peptides considered
                length(rownames(final_i_e)), #6 number of evidence peptides
                length(rownames(final_xr)), #7 number of cross reactive peptides
                final_rank[2],#8 - number of Filtered Evidence peptides
                paste(final_rank[[4]],collapse = "|"),#9 - number of Filtered Evidence peptides
                final_rank[3],#10 - Filtered N-rank number
                probability_i #11 -Null Probability
            )
            N_rank = z    ## set the N_rank for the next virus_i to be reduced by those assigned to the previous virus_i
            if(showPB){
                setTxtProgressBar(pb, R1)
            }
            
        }
        if(showPB){
            close(pb)
        }        
        index = which(output[,2]!=1)
        a123 = cbind(output,1)
        # a123[index,11] = p.adjust(output[index,2],"BH") # bh column
        a123[index,12] = p.adjust(output[index,2],"BH") #
        last = list(a123,sim_tag,fullmatrix_sorted_ij)
        return(last)
    }  
  
    # --- Simtag function --- #
    simtag  = function(last){
        results = as.data.frame(last[1]) #take the virus of two 
        results = cbind(results,0) # 
        table = as.data.frame(last[2]) # The indistinguisables
        names = subset(results[,1],results[,2]<.05)
        table = subset(table,table[,1] %in% names)
        table = subset(table,table[,2] %in% names)
        
        if(dim(table)[1]!=0){
            # create a network of indistinguisables
            net <- as.undirected(graph_from_data_frame(table, directed=F))
            # find the max cliques in the network (min size = 2)
            max =  max_cliques(net,min = 2)
            # for each clique in the network..
            for(R in 1:length(max)){
                # retrieve the vertex information and match the names of 
                # the vertices with the results table. These will be given
                # and indistinguishability tag id which is the clique id
                # from max_clique function for this patient
                names = induced_subgraph(net,max[[R]])
                sim = as.data.frame(vertex_attr(names))
                index = match(as.character(unlist(sim)),results[,1])
                results[index,13] = paste(results[index,13],R,sep = "|")
            }
        }
        return(results)
    }  
    
    # enrich = thresh # later implement so this can be changed??
    registerDoParallel(NUMCORE)    
    # plate = list()
    # cycle through each patient column by column (goal is so be serialized)
    zeta = foreach(R = 1:(dim(case)[2]-1),.combine=rbind) %dopar%{
        # run the subsetting step given input of the total null probs, sample column, 
        # enrichment threshold and hte Virus blast matrix.
        # rank <- total_calc(case,R,enrich,total,blast)
        # print(paste("R:",R,'of',dim(case)[2]-1))
        rank <- total_calc(case,R,thresh,total,blast)
        if(is.null(rank) == FALSE){
            # take subset matrix from above and do all reassignments with pairwise and total null probabilities
            sorted_table=pairwise_calc(total,pairwise, rank)
            # add indistinguishable tags
            sorted_table_2 = simtag(sorted_table) 
            name = colnames(case[R+1])
            output = as.data.frame(sorted_table_2)
            output[,13] =  gsub("^0\\||^0", '', output[,13])
            colnames(output) = c(
                "Virus", #2
                "P-value",#3
                "Evidence_Peptide_Ids",#4
                "XR_Peptides_Ids",
                "N-rank #",
                "Evidence_Peptide_Count",
                "XR_Peptide",
                "Filtered_Evidence_Count", #9
                "Filtered_Evidence_Peptides_Ids",
                "Filtered_N-rank #",#11
                "Null_Probability",#12
                "BH_P-value",#13
                "Indistinguishable_Groups" #14
            )
            fwrite(output,file = paste0(out_path,name,".csv"))
            pool = cbind(name,output)
            return(pool)
        }
    }
    fwrite(as.data.frame(zeta[zeta[,9]>=3 & zeta[,13]<=.05,]),file = paste0(out_path,out_name,"AVARDA_compiled_full_output",".csv"))
    empty_virus_1 = colnames(blast)[which(colnames(blast)%in%zeta$Virus == FALSE)]
    empty_virus = data.frame(matrix(NA,ncol = length(unique(zeta$name))+1,nrow = length(empty_virus_1)))
    empty_virus[,1] = empty_virus_1
    
    asdf = zeta %>% select(name, Virus,`Evidence_Peptide_Ids`) %>% spread(name,`Evidence_Peptide_Ids`,fill = 0)
    fwrite(
        rbindlist(
            list(asdf,empty_virus),
            use.names = FALSE
        ),
        file = paste0(out_path,out_name,"AVARDA_evidence_pep",".csv")
    )
    
    asdf = zeta %>% select(name, Virus,`Evidence_Peptide_Count`) %>% spread(name,`Evidence_Peptide_Count`,fill = 0)
    fwrite(
        rbindlist(
            list(asdf,empty_virus),
            use.names = FALSE
        ),
        file = paste0(out_path,out_name,"AVARDA_unfiltered_evidence_number",".csv")
    )
    
    asdf = zeta %>% select(name, Virus,`Filtered_Evidence_Count`) %>% spread(name,`Filtered_Evidence_Count`,fill = 0)
    fwrite(
        rbindlist(
            list(asdf,empty_virus),
            use.names = FALSE
        ),
        file = paste0(out_path,out_name,"AVARDA_breadth",".csv")
    )
    
    asdf = zeta %>% select(name, Virus,`BH_P-value`) %>% spread(name,`BH_P-value`,fill = 1)
    fwrite(
        rbindlist(
            list(asdf,empty_virus),
            use.names = FALSE
        ),
        file = paste0(out_path,out_name,"AVARDA_post_p_value_BH",".csv")
    )
    
    asdf = zeta %>% select(name, Virus,`P-value`) %>% spread(name,`P-value`,fill = 1)
    fwrite(
        rbindlist(
            list(asdf,empty_virus),
            use.names = FALSE
        ),
        file = paste0(out_path,out_name,"AVARDA_post_p_value",".csv")
    )
}

#################################
# Prepare data from user params #
#################################

if(length(opt) < 10){  
    helpMsg(opt)
}else{
    if(dir.exists(opt$out_path)){
        message(paste(opt$out_path, "exists! Using existing directory.."))
    }else{
        message(paste(opt$out_path,"not found! Creating new directory.."))
        dir.create(opt$out_path)
    }
    # just make sure that this results in a path otherwise
    # you'll get some funky names
    if(length(grep('\\/$', opt$out_path)) < 1){
        fixed_outpath = paste(opt$out_path, "/",sep = "")
    }else{
        fixed_outpath = opt$out_path
    }
    if(length(grep('_$', opt$out_name)) < 1){
        fixed_outname = paste(opt$out_name, "_",sep = "")
    }else{
        fixed_outname = opt$out_name
    }

    MAXCORE = detectCores()
    if(opt$cores != NUMCORE){
        if(opt$cores > MAXCORE){
            NUMCORE = MAXCORE
        }else{
            NUMCORE = opt$cores
        }
    }

    AVARDA(
        opt$case_path,
        as.numeric(opt$threshold),
        opt$dict_path,
        opt$total_path,
        opt$pairwise_path,
        opt$blast_path,
        fixed_outpath,
        fixed_outname,
        opt$avarda_names,
        opt$show_pb
    )
}

  