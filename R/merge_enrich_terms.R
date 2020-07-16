#' @title Merge enriched GO terms.
#' @description combine results from GO enrichment tests obtained with \pkg{topGO} package,
#' for a given ontology (MF, BP, or CC).
#' @importFrom data.table data.table rbindlist := .SD
#' @importFrom biomaRt useDataset getBM useEnsembl
#' @importFrom topGO termStat scoresInTerm
#' @importFrom AnnotationDbi select keys
#' @importFrom GO.db GO.db
#' @family GO_terms
#' @param Input a list containing named elements. Each element must contain the name of \code{\link[topGO]{topGOdata-class}}
#' object created by \code{\link{create_topGOdata}} method and the associated  \code{\link[topGO]{topGOresult-class}} object(s).
#' @param envir objects environment (default to .GlobalEnv).
#' @details This method extracts for each result of GO enrichment test (\code{\link[topGO]{topGOresult-class}} object) and
#' corresponding GO annotations (\code{\link[topGO]{topGOdata-class}} object):
#' informations about GO term (identifiant, name, and description),
#' gene frequency (number of significant genes / Annotated genes), pvalue, -log10(pvalue), significant genes
#' identifiants (GeneID, or Ensembl ID, or uniprot accession), and gene symbols.
#' At the last, this method builds a merged data.table of enriched GO terms (p<0.01)
#' at least once and provides all mentionned columns.
#' @return an \code{\link{enrich_GO_terms-class}} object.
#' @references
#' Alexa A and Rahnenfuhrer J (2016). topGO: Enrichment Analysis for Gene Ontology. R package version 2.28.0.
#'
#' Matt Dowle and Arun Srinivasan (2017). data.table: Extension of data.frame. R package version 1.10.4. https://CRAN.R-project.org/package=data.table
#'
#' Herve Pages, Marc Carlson, Seth Falcon and Nianhua Li (2017). AnnotationDbi: Annotation Database Interface. R package version 1.38.0.
#' @include enrich_GO_terms.R
#' @examples
#' # load genes identifiants (GeneID,ENS...) universe/background (Expressed genes)
#' background_L<-scan(
#'     system.file(
#'         "extdata/data/input",
#'         "background_L.txt",
#'         package = "ViSEAGO"
#'     ),
#'     quiet=TRUE,
#'     what=""
#' )
#'
#' # load Differentialy Expressed (DE) gene identifiants from files
#' PregnantvslactateDE<-scan(
#'     system.file(
#'         "extdata/data/input",
#'         "pregnantvslactateDE.txt",
#'         package = "ViSEAGO"
#'     ),
#'     quiet=TRUE,
#'     what=""
#' )
#'
#' VirginvslactateDE<-scan(
#'     system.file(
#'         "extdata/data/input",
#'         "virginvslactateDE.txt",
#'         package = "ViSEAGO"
#'     ),
#'     quiet=TRUE,
#'     what=""
#' )
#'
#' VirginvspregnantDE<-scan(
#'     system.file(
#'         "extdata/data/input",
#'         "virginvspregnantDE.txt",
#'         package="ViSEAGO"
#'     ),
#'     quiet=TRUE,
#'     what=""
#' )
#' \dontrun{
#' # create topGOdata for BP for each list of DE genes
#' BP_Pregnantvslactate<-ViSEAGO::create_topGOdata(
#'     geneSel=PregnantvslactateDE,
#'     allGenes=background_L,
#'     gene2GO=myGENE2GO,
#'     ont="BP",
#'     nodeSize=5
#' )
#'
#' BP_Virginvslactate<-ViSEAGO::create_topGOdata(
#'     geneSel=VirginvslactateDE,
#'     allGenes=background_L,
#'     gene2GO=myGENE2GO,
#'     ont="BP",
#'     nodeSize=5
#' )
#'
#' BP_Virginvspregnant<-ViSEAGO::create_topGOdata(
#'     geneSel=VirginvspregnantDE,
#'     allGenes=background_L,
#'     gene2GO=myGENE2GO,
#'     ont="BP",
#'     nodeSize=5
#' )
#'
#' # perform TopGO tests
#' elim_BP_Pregnantvslactate<-topGO::runTest(
#'     BP_L_pregnantvslactate,
#'     algorithm ="elim",
#'     statistic = "fisher"
#' )
#'
#' elim_BP_Virginvslactate<-topGO::runTest(
#'     BP_L_virginvslactate,
#'     algorithm ="elim",
#'     statistic = "fisher"
#' )
#'
#' elim_BP_Virginvspregnant<-topGO::runTest(
#'     BP_L_virginvspregnant,
#'     algorithm ="elim",
#'     statistic = "fisher"
#' )
#'
#' # merge topGO results
#' BP_sResults<-ViSEAGO::merge_enrich_terms(
#'     Input=list(
#'         Pregnantvslactate=c("BP_Pregnantvslactate","elim_BP_Pregnantvslactate"),
#'         Virginvslactate=c("BP_Virginvslactate","elim_BP_Virginvslactate"),
#'         Virginvspregnant=c("BP_Virginvspregnant","elim_BP_Virginvspregnant")
#'     )
#' )
#' }
#' @exportMethod merge_enrich_terms
#' @name merge_enrich_terms
#' @rdname merge_enrich_terms-methods
#' @exportMethod merge_enrich_terms
setGeneric(
    name="merge_enrich_terms",
    def=function(
        Input,
        envir=.GlobalEnv
    ){
        standardGeneric("merge_enrich_terms")
    }
)

#' @rdname merge_enrich_terms-methods
#' @aliases merge_enrich_terms
setMethod(
    "merge_enrich_terms",
    signature(
        Input="list"
    ),
    definition=function(Input,envir){

        ## check ontology

        # same ontology ?
        check.onto=unique(
            unlist(
                lapply(seq_along(Input),function(x){

                    # extract quering objects names
                    x=Input[[x]]

                    # check existence
                    values<-ls(envir=envir)

                    # check if available
                    available<-x%in%values

                    # stop if not found
                    if(!all(available)){
                        stop(paste("objects not found:",paste(x[!available],collapse=", "),sep="\n"))
                    }

                    # get objects
                    x=mget(x,envir=envir)

                    # objects type
                    obj.type=vapply(x,class,"")

                    # extract ontology type
                    vapply(seq_along(x),function(y){

                        # extract ontology slot
                        if(obj.type[y]=="topGOdata"){

                            # for topGOdata
                            slot(x[[y]],"ontology")

                        }else{

                            # for topGOresult
                            sub("^.+\nOntology: ","",slot(x[[y]],"description"))
                        }
                    },"")
                })
            )
        )

        # check ontology
        if(length(check.onto)>1){

            # stop if more than one topGOdata by list
            stop("Only one ontology supported")
        }

        ## topGO summary informations

        # build topGO summary
        topGO_summary=lapply(seq_along(Input),function(x){

            # extract  quering objects names
            x=Input[[x]]

            # keep names
            x_names=x

            # get objects
             x=mget(x,envir=envir)

            # objects type
            obj.type=vapply(x,class,"")

            # extract topGO objects summary
            topGO<-lapply(seq_along(x),function(y){

                # extract ontology slot
                if(obj.type[y]=="topGOdata"){

                    # for topGOdata
                     list(

                        # description
                        description=slot(x[[y]],"description"),

                        # availables genes
                        available_genes=length(slot(x[[y]],"allGenes")),

                        # availables genes significant
                        available_genes_significant=table(slot(x[[y]],"allScores"))[2],

                        # feasibles genes
                        feasible_genes=table(slot(x[[y]],"feasible"))[2],

                        # feasibles genes significant
                        feasible_genes_significant=table(slot(x[[y]],"allScores")==1 & slot(x[[y]],"feasible")==TRUE)[2],

                        # nodes with at least  x genes
                        genes_nodeSize=slot(x[[y]],"nodeSize"),

                        # number of nodes
                        nodes_number=length(slot(slot(x[[y]],"graph"),"nodes")),

                        # number of edges
                        edges_number=length(slot(slot(slot(x[[y]],"graph"),"edgeData"),"data"))
                    )

                }else{

                    # for topGO result
                    list(

                        # description
                        description=slot(x[[y]],"description"),

                        # test name
                        test_name=sub(": ","p<",slot(x[[y]],"testName")),

                        # algorithm name
                        algorithm_name=slot(x[[y]],"algorithm"),

                        # scored GOs
                        GO_scored=length(slot(x[[y]],"score")),

                        # significant GOs according cutOff
                        GO_significant=table(slot(x[[y]],"score")<as.numeric(sub("^.+[[:space:]]","",slot(x[[y]],"testName"))))[2],

                        # feasibles genes
                        feasible_genes=slot(x[[y]],"geneData")[1],

                        # feasibles genes significant
                        feasible_genes_significant=slot(x[[y]],"geneData")[2],

                        # nodes with at least  x genes
                        genes_nodeSize=slot(x[[y]],"geneData")[3],

                        # nodes with at least  x genes
                        Nontrivial_nodes=slot(x[[y]],"geneData")[4]
                    )
                }
            })

            # extract topGO objects summary
            names(topGO)<-x_names

            # return topGO
            topGO
        })

        # add names to topGO summary
        names(topGO_summary)<-names(Input)

        # find enrich GOs in a least one comparison
        GOs<-lapply(seq_along(Input),function(x){

            # objects type
            Data=mget(Input[[x]],envir=envir)

            ## checking step

            # objects type
            obj.type=vapply(Data,class,"")

            # objects type
            if(sum(obj.type%in%"topGOdata")>1){

                # stop if more than one godata by list
                stop("Only one topGOdata object is supported by list")
            }

            # objects type
            if(sum(obj.type%in%"topGOresult")>1){

                # stop if more than one godata by list
                stop("Only one topGOresult object is supported by list")
            }

            ## find and extract pvalues

            # find topGOresult
            pos=which(obj.type=="topGOresult")

            # extract scores
            pvalues<-topGO::score(Data[[pos]])
            
            # extract names of enrich terms
            as.vector(
                names(
                    pvalues[pvalues<as.numeric(sub("^.+[[:space:]]","",slot(Data[[pos]],"testName")))]
                )
            )
        })

        # remove redondancy and convert to vector
        GOs<-as.vector(unique(unlist(GOs)))

        ## check genes background

        # extract genes background
        allgenes<-lapply(seq_along(Input),function(x){

            # objects type
            Data=mget(Input[[x]],envir=envir)

            # objects type
            obj.type=vapply(Data,class,"")

            # load GOdata
            pos=which(obj.type=="topGOdata")

            slot(Data[[pos]],"allGenes")
        })

        # check if same gene background
        if(length(Input)>1){
            same_genes_background=all(
                vapply(2:length(allgenes),function(x){
                    identical(sort(allgenes[[1]]),sort(allgenes[[x]]))
                },TRUE)
            )
        }else{
            same_genes_background=TRUE
        }

        # stop if no enrich GO terms
        if(length(GOs)==0){
            stop("No enrich GO terms available in at least one condition")
        }

        # initialize input
        input=list()

        # combine results
        allResults<-lapply(seq_along(Input),function(x){

            ## extract Data

            # objects type
            Data=mget(Input[[x]],envir=envir)

            # objects type
            obj.type=vapply(Data,class,"")

            # load GOdata
            pos=which(obj.type=="topGOdata")

            # load GOdata
            GOdata=Data[[pos]]

            # tested algorithm
            algorithm<-Data[-pos]

            ## extract some statistics from initial GOdata object (before enrichment test)

            # get the GOdatatermStat(GOdata)
            Stats<-termStat(
                GOdata,
                whichGO = GOs
            )

            # convert to data.table
            Stats<-data.table(
                GO.ID=row.names(Stats),
                genes_frequency=paste(
                    round(Stats$Significant/Stats$Annotated*100,digits=2),
                    "% (",Stats$Significant,"/",Stats$Annotated,")",
                    sep=""
                )
            )

            ## extract genes identifiants

            # extract counts genes by term according GeneList
            genes<-scoresInTerm(
                GOdata,GOs,
                use.names = TRUE
            )

            # extract  genes  Ids according GeneList
            genes=lapply(names(genes),function(x){

                # extract significant terms
                val=attr(genes[[x]][genes[[x]]==2],"names")

                # build a table
                data.table(
                    GO.ID=x,
                    Significant_genes=if(length(val)>0){val}else{NA}
                )
            })

            # convert to data.table
            genes<-rbindlist(genes)

            ## add Genes symbols

            # get db
            db=strsplit(slot(GOdata,"description")," ")[[1]]

            # if db match to bioconductor
            if(db[1]=="Bioconductor"){

                # load GeneID and symbols
                annot<-data.table(
                    select(
                        get(db[2]),
                        keys=keys(get(db[2])),
                        columns =c("ENTREZID","SYMBOL")
                    )
                )

                # if found symbols
                if(nrow(annot[!is.na(annot$ENTREZID)])>0){

                    # load GeneID and symbols
                    genes<-merge(
                        genes,
                        annot,
                        by.x="Significant_genes",
                        by.y="ENTREZID",
                        all.x=TRUE
                    )

                }else{

                    # add empty name columns
                    genes[,Name:=NA]
                }
            }

            # if db  match to enstrezGene
            if(db[1]=="EntrezGene"){

                # function for generate data packets defined by the "by" argument
                pos=function(Data,by=""){
                    if(by>length(Data)){
                        by=length(Data)
                    }
                    if(length(Data)<=by){
                        data.table(
                            start=1,
                            end=length(Data)
                        )
                    }else{
                        data.table(
                            start=seq(1,length(Data),by=by),
                            end=unique(c(seq(by,length(Data),by=by),length(Data)))
                        )
                    }
                }

                # pattern.extract
                pattern.extract=function(query,m){
                    query=lapply(seq_along(query),function(i){

                        if(length(na.omit(m[[i]][1]))>0){
                            a=attr(m[[i]],"capture.start")
                            t(substring(query[i],a,a+attr(m[[i]],"capture.length")-1))
                        }else{
                            NA
                        }
                    })
                    query=do.call("rbind",query)
                    query[query%in%""]<-NA
                    query[query%in%"\t"]<-NA
                    query
                }

                # esummary
                esummary<-function(...){

                    # submitted Data
                    Data<-na.omit(unique(unique(...)))

                    # batch size
                    wpos=pos(Data,by=400)

                    # core address
                    core="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?version=2.0&db="

                    # Data submission an retrieve
                    query=lapply(seq_len(nrow(wpos)),function(x){

                        # build block of id
                        Dat<-paste(
                            Data[wpos[x,start]:wpos[x,end]],
                            collapse=","
                        )

                        # submit block
                        query <- paste(core,"gene","&id=",Dat,sep = "")

                        # read and return reults
                        readLines(query)
                    })
                    query=paste(
                        unlist(query),
                        collapse=""
                    )

                    # parse results
                    query<-substring(
                        query,
                        unlist(gregexpr("<DocumentSummary ",query)),
                        unlist(gregexpr("</DocumentSummary>",query))
                    )

                    # extraction pattern
                    pattern="<DocumentSummary uid=\"(?<Id>[[:digit:]]*)\">\t<Name>(?<Name>.*)</Name>"

                    # locate elements
                    m1=gregexpr(
                        paste(pattern,collapse=""),
                        query,
                        perl=TRUE
                    )

                    # extract elements
                    res=data.table(
                        pattern.extract(query,m1)
                    )

                    # add an header
                    names(res)=attr(m1[[1]],"capture.name")

                    # return
                    return(res)
                }

                # get esummary from NCBI
                annot<-esummary(genes$Significant_genes)

                # if found symbols
                if(nrow(annot[!is.na(annot$Id)])>0){

                    # load GeneID and symbols
                    genes<-merge(
                        genes,
                        annot,
                        by.x="Significant_genes",
                        by.y="Id",
                        all.x=TRUE
                    )

                }else{

                    # add empty name columns
                    genes[,Name:=NA]
                }
            }

            # if db  match to Ensembl
            if(db[1]=="Ensembl"){

                # connect to Ensembl
                mart<-useEnsembl(
                    "ensembl",
                    host=db[2],
                    version=db[6]
                )

                # connect to ensembl specified dataset
                myspecies<-useDataset(
                    db[2],
                    mart
                )

                # load Ensembl genes, with GO annotations
                annot<-data.table(
                    getBM(
                        attributes =c("ensembl_gene_id","external_gene_name"),
                        value=TRUE,
                        mart =myspecies
                    )
                )

                # if found symbols
                if(nrow(annot[!"",on="external_gene_name"])>0){

                    # load GeneID and symbols
                    genes<-merge(
                        genes,
                        annot,
                        by.x="Significant_genes",
                        by.y="ensembl_gene_id",
                        all.x=TRUE
                    )

                }else{

                    # add empty name columns
                    genes[,Name:=NA]
                }
            }

            # if db  match to Custom
            if(db[1]=="Custom"){

                # load GeneID and symbols
                annot<-unique(
                    fread(
                        db[3],
                        select=c("gene_id","gene_symbol")
                    )
                )

                # merge gene_symbol
                genes<-merge(
                    genes,
                    annot,
                    by.x="Significant_genes",
                    by.y="gene_id",
                    all.x=TRUE
                )
            }

            # load GeneID and symbols
            names(genes)[3]<-"Significant_genes_symbol"

            # reorder columns
            genes<-genes[,c("GO.ID","Significant_genes","Significant_genes_symbol"),with=FALSE]

            # collapse results
            genes<-genes[,lapply(.SD,function(x){paste(x,collapse=";")}),.SDcols=2:3,by="GO.ID"]

            # replace all blank cells by NA
            genes[genes==""]<-NA

            ## extract pvalue according the algorithm results

            # extract all pvalues from topGOresult
            pvalue<-topGO::score(algorithm[[1]])
            
            # select pvalues from topGOresult
            pvalue<-pvalue[GOs]
            
            # extract pvalue in data.table
            pvalues<-data.table(
                GO.ID=names(pvalue),
                pvalue=as.numeric(format(pvalue,scientific = T)),
                `-log10_pvalue`=round(-log10(pvalue),digits=2)
            )

            # algorithm
            algo=slot(algorithm[[1]],"algorithm")

            # return input params
            assign(
                "input",
                c(input,list(algo)),
                inherits=TRUE
            )

            ## combine results

            # all results in list
            Results<-list(
                data.table(GO.ID=GOs),
                Stats,
                pvalues,
                genes
            )

            # merge all
            Results<-Reduce(
                function(...){
                    merge(
                        ...,
                        by ="GO.ID",
                        sort=FALSE,
                        all=TRUE
                    )
                },
                Results
            )

            # remove NA in GO.Id column
            Results<-Results[!is.na(Results$GO.ID)]

            # extract pvalue threshold
            p<-as.numeric(sub("^.+[[:space:]]","",slot(algorithm[[1]],"testName")))

            # Remove gene ID and symbol if GO term not significant
            Results[Results$pvalue>=p,`:=`(Significant_genes=NA,Significant_genes_symbol=NA)]

            if(!is.null(names(Input))){

                # add GOdata name in the header
                names(Results)[-1]<-paste(
                    names(Input)[x],
                    names(Results)[-1],
                    sep="."
                )
            }

            # return Results
            Results
        })

        # number of Godata
        nb=length(allResults)

        # merge results if more than GOdata input
        if(nb>1){

            # merge all results
            allResults=Reduce(
                function(...){
                    merge(...,
                          by ="GO.ID",
                          sort=FALSE,
                          all=TRUE
                    )
                },
                allResults
            )

        }else{

            # if only one GOdata unlist table
            allResults<-allResults[[1]]
        }

        # add Input names
        names(input)<-names(Input)

        # extract GO terms
        GO<-data.table(
            select(
                GO.db,
                keys=allResults$GO.ID,
                columns = c("GOID","TERM","DEFINITION")
            )
        )

        # add GO term description and definition to sResults
        allResults<-data.table(
            GO,
            allResults[,"GO.ID":=NULL]
        )

        # rename the first 3 columns
        names(allResults)[seq_len(3)]<-c("GO.ID","term","definition")

        # significant results in at least one condition
        new(
            "enrich_GO_terms",
            same_genes_background=same_genes_background,
            input=input,
            ont=check.onto,
            topGO=topGO_summary,
            data=allResults
        )
    }
)
