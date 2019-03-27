#' @title Merge enriched GO terms.
#' @description combine results from GO enrichment tests obtained with \pkg{topGO} package,
#' for a given ontology (MF, BP, or CC).
#' @importFrom data.table data.table rbindlist
#' @importFrom topGO scoresInTerm
#' @importFrom methods setGeneric setMethod slot is
#' @importFrom biomaRt useDataset getBM
#' @importFrom topGO termStat
#' @importFrom AnnotationDbi select
#' @family GO_terms
#' @param Input a list containing named elements. Each element must contain the name of \code{\link[topGO]{topGOdata-class}}
#' object created by \code{\link{create_topGOdata}} method and the associated  \code{\link[topGO]{topGOresult-class}} object(s).
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
#' \dontrun{
#'  ###################
#'  # load genes identifiants (GeneID,ENS...) universe/background (Expressed genes)
#'  background_L<-base::scan(
#'   base::system.file("extdata/data/input","background_L.txt",package = "ViSEAGO"),
#'   quiet=T,what=""
#'  )
#'
#'  ###################
#'  # load Differentialy Expressed (DE) gene identifiants from files
#'  L_pregnantvslactateDE<-base::scan(
#'   base::system.file("extdata/data/input","L_pregnantvslactateDE.txt",package = "ViSEAGO"),
#'   quiet=T,what=""
#'  )
#'
#'  L_virginvslactateDE<-base::scan(
#'   base::system.file("extdata/data/input","L_virginvslactateDE.txt",package = "ViSEAGO"),
#'   quiet=T,what=""
#'  )
#'
#'  L_virginvspregnantDE<-base::scan(
#'   base::system.file("extdata/data/input","L_virginvspregnantDE.txt",package = "ViSEAGO"),
#'   quiet=T,what=""
#'  )
#'
#' ###################
#' # create topGOdata for BP for each list of DE genes
#' BP_L_pregnantvslactate<-ViSEAGO::create_topGOdata(
#'  geneSel=L_pregnantvslactateDE,
#'  allGenes=background_L,
#'  gene2GO=myGENE2GO,
#'  ont="BP",
#'  nodeSize=5
#' )
#'
#' BP_L_virginvslactate<-ViSEAGO::create_topGOdata(
#'  geneSel=L_virginvslactateDE,
#'  allGenes=background_L,
#'  gene2GO=myGENE2GO,
#'  ont="BP",
#'  nodeSize=5
#' )
#'
#' BP_L_virginvspregnant<-ViSEAGO::create_topGOdata(
#'  geneSel=L_virginvspregnantDE,
#'  allGenes=background_L,
#'  gene2GO=myGENE2GO,
#'  ont="BP",
#'  nodeSize=5
#' )
#'
#' ###################
#' # perform TopGO tests
#' elim_BP_L_pregnantvslactate<-topGO::runTest(
#'  BP_L_pregnantvslactate,
#'  algorithm ="elim",
#'  statistic = "fisher"
#' )
#'
#' elim_BP_L_virginvslactate<-topGO::runTest(
#'  BP_L_virginvslactate,
#'  algorithm ="elim",
#'  statistic = "fisher"
#' )
#'
#' elim_BP_L_virginvspregnant<-topGO::runTest(
#'  BP_L_virginvspregnant,
#'  algorithm ="elim",
#'  statistic = "fisher"
#' )
#'
#' ###################
#' # merge topGO results
#' BP_sResults<-ViSEAGO::merge_enrich_terms(
#'  Input=base::list(
#'   L_pregnantvslactate=base::c("BP_L_pregnantvslactate","elim_BP_L_pregnantvslactate"),
#'   L_virginvslactate=base::c("BP_L_virginvslactate","elim_BP_L_virginvslactate"),
#'   L_virginvspregnant=base::c("BP_L_virginvspregnant","elim_BP_L_virginvspregnant")
#'  )
#' )
#' }
#' @export
setGeneric(name="merge_enrich_terms", def=function(Input){standardGeneric("merge_enrich_terms")})

setMethod("merge_enrich_terms",definition=function(Input){

  ###################
  # check ontology
  ###################

    ###################
    # same ontology ?
    check.onto=base::unique(base::c(base::sapply(base::seq_along(Input),function(x){

      ###################
      # extract  quering objects names
      x=Input[[x]]

      ###################
      # check existance
      values<-ls(envir = .GlobalEnv)

      ###################
      # check existance
      available<-x%in%values

      ###################
      # check existance
      if(!base::all(available))base::stop(base::paste("objects not found:",paste(x[!available],collapse=", "),sep="\n"))

      ###################
      # get objects
      x=base::mget(x,envir=.GlobalEnv)

      ###################
      # objects type
      obj.type=base::sapply(x,class)

      ###################
      # extract ontoloy type
      base::sapply(base::seq_along(x),function(y){

        ###################
        # extract ontology slot
        if(obj.type[y]=="topGOdata"){

          ###################
          # for topGOdata
          methods::slot(x[[y]],"ontology")

        }else{

          ###################
          # for topGOresult
          base::sub("^.+\nOntology: ","",methods::slot(x[[y]],"description"))
        }
      })
    })))

    ###################
    # check ontology
    if(base::length(check.onto)>1){

      ###################
      # stop if more than one godata by list
      base::stop("Only one ontology supported")
    }

  ###################
  # topGO summary informations
  ###################

    ###################
    # buil topGO summary
    topGO_summary=base::lapply(base::seq_along(Input),function(x){

      ###################
      # extract  quering objects names
      x=Input[[x]]

      ###################
      # keep names
      x_names=x

      ###################
      # get objects
      x=base::mget(x,envir=.GlobalEnv)

      ###################
      # objects type
      obj.type=base::sapply(x,class)

      ###################
      # extract topGO objects summary
      topGO<-base::lapply(base::seq_along(x),function(y){

        ###################
        # extract ontology slot
        if(obj.type[y]=="topGOdata"){

          ###################
          # for topGOdata
          base::list(

            ###################
            # description
            description=methods::slot(x[[y]],"description"),

            ###################
            # availables genes
            available_genes=base::length(methods::slot(x[[y]],"allGenes")),

            ###################
            # availables genes significant
            available_genes_significant=base::table(methods::slot(x[[y]],"allScores"))[2],

            ###################
            # feasibles genes
            feasible_genes=base::table(methods::slot(x[[y]],"feasible"))[2],

            ###################
            # feasibles genes significant
            feasible_genes_significant=base::table(methods::slot(x[[y]],"allScores")==1 & methods::slot(x[[y]],"feasible")==T)[2],

            ###################
            # nodes with at least  x genes
            genes_nodeSize=methods::slot(x[[y]],"nodeSize"),

            ###################
            # number of nodes
            nodes_number=base::length(methods::slot(methods::slot(x[[y]],"graph"),"nodes")),

            ###################
            # number of edges
            edges_number=base::length(methods::slot(methods::slot(methods::slot(x[[y]],"graph"),"edgeData"),"data"))
          )

        }else{

          ###################
          # for topGO result
          base::list(

            ###################
            # description
            description=methods::slot(x[[y]],"description"),

            ###################
            # test name
            test_name=base::sub(": ","p<",methods::slot(x[[y]],"testName")),

            ###################
            # algorithm name
            algorithm_name=methods::slot(x[[y]],"algorithm"),

            ###################
            # scored GOs
            GO_scored=base::length(methods::slot(x[[y]],"score")),

            ###################
            # significant GOs
            GO_significant=base::table(methods::slot(x[[y]],"score")<0.01)[2],

            ###################
            # feasibles genes
            feasible_genes=methods::slot(x[[y]],"geneData")[1],

            ###################
            # feasibles genes significant
            feasible_genes_significant=methods::slot(x[[y]],"geneData")[2],

            ##################
            # nodes with at least  x genes
            genes_nodeSize=methods::slot(x[[y]],"geneData")[3],

            ##################
            # nodes with at least  x genes
            Nontrivial_nodes=methods::slot(x[[y]],"geneData")[4]
        )
      }

      })

      ###################
      # extract topGO objects summary
      base::names(topGO)<-x_names

      ###################
      # return topGO
      topGO
    })

    ###################
    # add names to topGO summary
    base::names(topGO_summary)<-base::names(Input)

  ###################
  # find enrich GOs in a least one comparison
  GOs<-base::lapply(base::seq_along(Input),function(x){

    ###################
    # objects type
    Data=base::mget(Input[[x]],envir=.GlobalEnv)

    ###################
    # checking step
    ###################

      ###################
      # objects type
      obj.type=base::sapply(Data,class)

      ###################
      # objects type
      if(base::sum(obj.type%in%"topGOdata")>1){

        ###################
        # stop if more than one godata by list
        base::stop("Only one topGOdata object is supported by list")
      }

      ###################
      # objects type
      if(sum(obj.type%in%"topGOresult")<1){

        ###################
        # stop if more than one godata by list
        base::stop("At least one topGOresult object needed by list")
      }

    ###################
    # find and extract pvalues
    ###################

      ###################
      # load GOdata
      pos=base::which(obj.type=="topGOresult")

      ###################
      # tested algorithm
      algorithms<-Data[pos]

      ###################
      # extract significvant pvalues results
      base::unlist(lapply(algorithms,function(y){

        ##################
        # extract scores
        pvalues<-topGO::score(y)

        ##################
        # extract names of enrich terms
        base::as.vector(names(pvalues[pvalues<0.01]))
      }))
    })

  ###################
  # remove redondancy and convert to vector
  GOs<-base::as.vector(base::unique(base::unlist(GOs)))

  ###################
  # check genes background
  ###################

    ###################
    # extract genes background
    allgenes<-base::lapply(base::seq_along(Input),function(x){

      ###################
      # objects type
      Data=base::mget(Input[[x]],envir=.GlobalEnv)

      ###################
      # objects type
      obj.type=base::sapply(Data,class)

      ###################
      # load GOdata
      pos=base::which(obj.type=="topGOdata")

      ###################
      methods::slot(Data[[pos]],"allGenes")
    })

    ###################
    # check if same gene background
    same_genes_background=base::all(
      base::sapply(2:length(allgenes),function(x){
        base::identical(allgenes[[1]],allgenes[[x]])
      })
    )

  ###################
  # stop if no enrich GO terms
  if(base::length(GOs)==0){base::stop("No enrich GO terms available in at least one condition")}

  ###################
  # initialyse input
  input=base::list()

  ###################
  # combine results
  allResults<-base::lapply(base::seq_along(Input),function(x){

    ###################
    # extract Data
    ###################

      ###################
      # objects type
      Data=base::mget(Input[[x]],envir=.GlobalEnv)

      ###################
      # objects type
      obj.type=base::sapply(Data,class)

      ###################
      # load GOdata
      pos=base::which(obj.type=="topGOdata")

      ###################
      # load GOdata
      GOdata=Data[[pos]]

      ###################
      # tested algorithm
      algorithms<-Data[-pos]

    ###################
    # extract some statistics from initial GOdata object (before enrichment test)
    ###################

      ###################
      # get the GOdatatermStat(GOdata)
      Stats<-topGO::termStat(GOdata, whichGO = GOs)

      ###################
      # convert to data.table
      Stats<-data.table::data.table(
        GO.ID=base::row.names(Stats),
        genes_frequency=base::paste(
        base::round(Stats$Significant/Stats$Annotated*100,digits=2),
        "% (",Stats$Significant,"/",Stats$Annotated,")",sep="")
      )

    ###################
    # extract genes identifiants
    ###################

      ##################
      # extract counts genes by term according GeneList
      genes<-topGO::scoresInTerm(GOdata,GOs,use.names = TRUE)

      ##################
      # extract  genes  Ids according GeneList
      genes=base::lapply(base::names(genes),function(x){

        ##################
        # extract significant terms
        val=base::attr(genes[[x]][genes[[x]]==2],"names")

        ##################
        # build a table
        data.table::data.table(
          GO.ID=x,
          Significant_genes=if(base::length(val)>0){val}else{NA}
        )
      })

      ##################
      # convert to data.table
      genes<-data.table::rbindlist(genes)

    ###################
    # add Genes symbols
    ###################

      ##################
       # get db
      db=base::strsplit(methods::slot(GOdata,"description")," ")[[1]]

      ##################
      # if db  match to bioconductor
      if(db[1]=="Bioconductor"){

        ###################
        # load GeneID and symbols
        annot<-data.table::data.table(AnnotationDbi::select(base::get(db[2]),
        keys=AnnotationDbi::keys(base::get(db[2])),columns =c("ENTREZID","SYMBOL")))

        ##################
        # if found symbols
        if(nrow(annot[!is.na(ENTREZID)])>0){

          # ###################
          # load GeneID and symbols
          genes<-merge(genes,annot,by.x="Significant_genes",by.y="ENTREZID",all.x=T)

        }else{

          ###################
          # add empty name columns
          genes[,Name:=NA]
        }
      }

      ##################
      # if db  match to enstrezGene
      if(db[1]=="EntrezGene"){

        ###################
        # function for generate data packets defined by the "by" argument
        pos=function(Data,by=""){
        if(by>base::length(Data)){
          by=base::length(Data)
        }
        if(base::length(Data)<=by){
          data.table::data.table(start=1,end=base::length(Data))
        }else{
          data.table::data.table(start=base::seq(1,length(Data),by=by),
          end=base::unique(base::c(base::seq(by,base::length(Data),by=by),base::length(Data))))
        }
      }

        ###################
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

        ##################
        # esummary
        esummary<-function(...){

          ###################
          # submitted Data
          Data<-stats::na.omit(base::unique(base::unique(...)))

          ###################
          # batch size
          wpos=pos(Data,by=500)

          ###################
          # core address
          core="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?version=2.0&db="

          ###################
          # Data submission an retrieve
          query=base::lapply(1:base::nrow(wpos),function(x){

            ###################
            # build block of id
            Dat<-base::paste(Data[wpos[x,start]:wpos[x,end]],collapse=",")

            ###################
            # submit block
            query <- base::paste(core,"gene","&id=",Dat,sep = "")

            ###################
            # read and return reults
            base::readLines(query)
          })
          query=base::paste(base::unlist(query),collapse="")

          ###################
          # parse results
          query<-base::substring(
            query,base::unlist(base::gregexpr("<DocumentSummary ",query)),
            base::unlist(base::gregexpr("</DocumentSummary>",query)))

          ###################
          # extraction pattern
          pattern="<DocumentSummary uid=\"(?<Id>[[:digit:]]*)\">\t<Name>(?<Name>.*)</Name>"

          ###################
          # locate elements
          m1=base::gregexpr(base::paste(pattern,collapse=""),query,perl=T)

          ###################
          # extract elements
          res=data.table::data.table(pattern.extract(query,m1))

          ###################
          # add an header
          base::names(res)=attr(m1[[1]],"capture.name")

          ###################
          # return
          res
        }

        ##################
        # get esummary from NCBI
        annot<-esummary(genes$Significant_genes)

        ##################
        # if found symbols
        if(nrow(annot[!is.na(Id)])>0){

          ###################
          # load GeneID and symbols
          genes<-merge(genes,annot,by.x="Significant_genes",by.y="Id",all.x=T)

        }else{

          ###################
          # add empty name columns
          genes[,Name:=NA]
        }
      }

      ##################
      # if db  match to bioconductor
      if(db[1]=="Ensembl"){

        ###################
        # connect to Ensembl
        mart<-biomaRt::useEnsembl("ensembl",host=db[2],version=db[6])

        ###################
        # connect to ensembl specified dataset
        myspecies<-biomaRt::useDataset(db[2],mart)

        ###################
        # load Ensembl genes, with GO annotations
        annot<-data.table::data.table(
          biomaRt::getBM(
            attributes =base::c("ensembl_gene_id","external_gene_name"),
            value=TRUE,
            mart =myspecies
          )
        )

        ##################
        # if found symbols
        if(base::nrow(annot[external_gene_name!=""])>0){

          ###################
          # load GeneID and symbols
          genes<-merge(genes,annot,by.x="Significant_genes",by.y="ensembl_gene_id",all.x=T)

        }else{

          ###################
          # add empty name columns
          genes[,Name:=NA]
        }
      }

      ###################
      # load GeneID and symbols
      names(genes)[3]<-"Significant_genes_symbol"

      ###################
      # reorder columns
      genes<-genes[,.(GO.ID,Significant_genes,Significant_genes_symbol)]

      ##################
      # collapse results
      genes<-genes[,
        base::lapply(.SD,function(x){base::paste(x,collapse=";")}),
        .SDcols=2:3,by=GO.ID]

      ##################
      # replace all blank cells by NA
      genes[genes==""]<-NA

    ###################
    # extract pvalue according the algorithm results
    ###################

    ###################
    # extract pvalues results
    pvalues<-base::lapply(base::seq_along(algorithms),function(y){

      ##################
      # extract all pvalues from topGOresult
      pvalue<-topGO::score(algorithms[[y]])

      ##################
      # select pvalues from topGOresult
      pvalue<-pvalue[GOs]

      ###################
      # extract pvalue in data.table
      pvalue<-data.table::data.table(
        GO.ID=base::names(pvalue),
        pvalue=base::as.numeric(base::format(pvalue,scientific = T)),
        base::round(-base::log10(pvalue),digits=2)
      )

      ###################
      # ordering pvalue by go term
      base::names(pvalue)[ncol(pvalue)]="-log10_pvalue"

      if(y>1){

        ###################
        # remove GO.ID
        pvalue[,GO.ID:=NULL]
      }

      ###################
      # return
      pvalue
    })

    ##################
    # convert to data.table
    pvalues=data.table::data.table(base::do.call("cbind",pvalues))

    ##################
    # algoritms
    algorithms=base::sapply(algorithms,function(x)x@algorithm)

    ##################
    # if use of different algorithms
    if(base::length(algorithms)>1){

      ##################
      # add names
      base::names(pvalues)[-1]<-base::paste(base::rep(algorithms,each=2),base::names(pvalues)[-1],sep=".")
    }

    ##################
    # return input params
    input<<-base::c(input,list(algorithms))

    ###################
    # combine results
    ###################

    ##################
    # all results in list
    Results<-base::list(data.table::data.table(GO.ID=GOs),Stats,pvalues,genes)

    ##################
    # merge all
    Results<-base::Reduce(function(...) merge(..., by ="GO.ID", sort=F,all = T),Results)

    ##################
    # remove NA in GO.Id column
    Results<-Results[!base::is.na(GO.ID)]

    ##################
    # Remove gene ID and symbol if GO term not significant
    Results[pvalue>=0.01,
      `:=`(
        Significant_genes=NA,
        Significant_genes_symbol=NA
      )
    ]

    if(!base::is.null(base::names(Input))){

      ##################
      # add GOdata name
      if(x==1){

        ##################
        # add GOdata name in the header
        base::names(Results)[-1]<-base::paste(base::names(Input)[x],
        base::names(Results)[-1],sep=".")

      }else{

        ##################
        # add GOdata name in the header
        base::names(Results)[-1]<-base::paste(base::names(Input)[x],
        base::names(Results)[-1],sep=".")
      }
    }

    ##################
    # return Results
    Results
  })

  ###################
  # number of Godata
  nb=base::length(allResults)

  ###################
  # merge results if more than GOdata input
  if(nb>1){

    ###################
    # merge all results
    allResults=base::Reduce(function(...) merge(..., by ="GO.ID", sort=F,all = T), allResults)

  }else{

    ##################
    # if only one GOdata unlist table
    allResults<-allResults[[1]]
  }

  ##################
  # add Input names
  names(input)<-names(Input)

  ###################
  # extract GO terms
  GO<-data.table::data.table(AnnotationDbi::select(GO.db::GO.db,
    keys=allResults$GO.ID,
    columns = c("GOID","TERM","DEFINITION"))
  )

  ##################
  # add GO term description and definition to sResults
  allResults<-data.table::data.table(GO,allResults[,GO.ID:=NULL])

  ##################
  # rename the fisrt 3 columns
  names(allResults)[1:3]<-base::c("GO.ID","term","definition")

  ##################
  # significant results in at least one condition
  methods::new("enrich_GO_terms",
               same_genes_background=same_genes_background,
               input=input,
               ont=check.onto,
               topGO=topGO_summary,
               data=allResults
  )
})
