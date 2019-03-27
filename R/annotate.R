#' @title Retrieve GO annotations for a specie from genomic ressource database.
#' @description This method retrieves and stores GO annotations for the
#' organism of interest from one of genomic ressource database
#' (Bioconductor, EntrezGene, Ensembl, Uniprot).
#' @importFrom methods setGeneric setMethod new is slot
#' @importFrom AnnotationDbi select keys
#' @importFrom biomaRt useDataset getBM
#' @importFrom data.table data.table
#' @importFrom R.utils gunzip
#' @family genomic_ressource
#' @family GO_terms
#' @param id identifiant corresponding to the organism of interest.
#' This id name is referenced in the first column of the database
#' used (see \code{\link{available_organisms}}).
#' @param object a required \code{\link{genomic_ressource-class}} object created by
#' \code{\link{Bioconductor2GO}}, \code{\link{EntrezGene2GO}},
#' \code{\link{Ensembl2GO}}, or \code{\link{Uniprot2GO}} methods.
#' @param ortholog \code{logical} (default to FALSE). Only available for
#' vertebrates organisms and for object created by \code{\link{EntrezGene2GO}} method (see Details).
#' @details This method uses a \code{\link{genomic_ressource-class}} object to retrieve
#' \href{http://www.geneontology.org/page/ontology-documentation}{GO} annotations for the organism of interest.
#' The stored annotations are structured in 3 slots corresponding to the 3 GO categories: MF (Molecular Function),
#' BP (Biological Process), and CC (Cellular Component). Each slot contains GO terms with
#' associated \href{http://www.geneontology.org/page/guide-go-evidence-codes}{evidence code}.
#'
#' The \code{\link{genomic_ressource-class}} object is created by one of the four available methods:
#' \code{\link{Bioconductor2GO}}, \code{\link{EntrezGene2GO}},
#' \code{\link{Ensembl2GO}}, or \code{\link{Uniprot2GO}}.
#'
#' In the case of vertebrates, setting \code{ortholog} argument to \code{TRUE} is required if you need to add GO terms with experimental
#' \href{http://geneontology.org/page/guide-go-evidence-codes}{evidence codes} from orthologs genes
#' when using \code{\link{EntrezGene2GO}} method. To display organisms supported by NCBI EntrezGene orthologs pipeline,
#' set the arguments \code{id=NULL} and \code{ortholog=TRUE}.
#' This approch is highly similar to the strategy developed by Uniprot-GOA consortium for the Electronic Annotation Method using
#' \href{http://www.ebi.ac.uk/GOA/compara_go_annotations}{Ensembl Compara}.
#' @return \code{annotate} produces an object of \code{\link{gene2GO-class}} required by \code{\link{build_GO_SS}} method.
#' @references
#' Durinck S, Spellman P, Birney E and Huber W (2009). Mapping identifiers for the integration of genomic datasets with the R/Bioconductor
#' package biomaRt. Nature Protocols, 4, pp. 1184-1191.
#'
#' Durinck S, Moreau Y, Kasprzyk A, Davis S, De Moor B, Brazma A and Huber W (2005).
#' BioMart and Bioconductor: a powerful link between biological databases and microarray data analysis. Bioinformatics, 21, pp. 3439-3440.
#'
#' Fong, JH, Murphy, TD, Pruitt, KD (2013). Comparison of RefSeq protein-coding regions in human and vertebrate genomes. BMC Genomics, 14:654.
#'
#' Henrik Bengtsson (2016). R.utils: Various Programming Utilities. R package version 2.5.0. https://CRAN.R-project.org/package=R.utils.
#'
#' Herve Pages, Marc Carlson, Seth Falcon and Nianhua Li (2017). AnnotationDbi: Annotation Database Interface. R package version 1.38.0.
#'
#' Matt Dowle and Arun Srinivasan (2017). data.table: Extension of data.frame. R package version 1.10.4. https://CRAN.R-project.org/package=data.table.
#' @include genomic_ressource.R
#' @examples
#'
#' ###################
#' # load Mus musculus (mouse) GO annotations
#' ###################
#'
#'  ###################
#'  # from Bioconductor
#'  Bioconductor<-ViSEAGO::Bioconductor2GO()
#'  myGENE2GO<-ViSEAGO::annotate("org.Mm.eg.db",Bioconductor)
#'
#'  ###################
#'  # from EntrezGene
#'  EntrezGene<-ViSEAGO::EntrezGene2GO()
#'  myGENE2GO<-ViSEAGO::annotate("10090",EntrezGene)
#'
#'  ###################
#'  # from EntrezGene
#'  Ensembl<-ViSEAGO::Ensembl2GO()
#'  myGENE2GO<-ViSEAGO::annotate("mmusculus_gene_ensembl",Ensembl)
#'
#'  ###################
#'  # from Uniprot
#'  Uniprot<-ViSEAGO::Uniprot2GO()
#'  myGENE2GO<-ViSEAGO::annotate("mouse",Uniprot)
#'
#' ###################
#' # specific options for EntrezGene database
#' ###################
#'
#'  ##################
#'   # Chicken GO annotations without adding orthologs
#'  EntrezGene<-ViSEAGO::EntrezGene2GO()
#'  myGENE2GO<-ViSEAGO::annotate("9031",EntrezGene)
#'
#'  ##################
#'  # Chicken GO annotation with the add of orthologs GO annotations
#'  EntrezGene<-ViSEAGO::EntrezGene2GO()
#'  myGENE2GO<-ViSEAGO::annotate("9031",EntrezGene, ortholog=TRUE)
#'
#'  ##################
#'  # display organisms supported by NCBI EntrezGene orthologs pipeline
#'  EntrezGene<-ViSEAGO::EntrezGene2GO()
#'  ViSEAGO::annotate(NULL,EntrezGene, ortholog=TRUE)
#'
#'  ##################
#'  # Coturnix japonica GO annotations with adding orthologs
#'  EntrezGene<-ViSEAGO::EntrezGene2GO()
#'  myGENE2GO<-ViSEAGO::annotate("93934",EntrezGene, ortholog=TRUE)
#' @export
setGeneric(name="annotate",def=function(id,object,ortholog=FALSE){standardGeneric("annotate")})
#' @importFrom methods setMethod
setMethod("annotate",definition=function(id,object,ortholog){

  ###################
  # check object
  ###################
  if(!methods::is(object,"genomic_ressource"))base::stop("object must be a genomic_ressource class from ViSEAGO::Bioconductor2GO(), ViSEAGO::EntrezGene2GO(), or ViSEAGO::Ensembl2GO()")
  if(methods::slot(object,"db")!="EntrezGene" & ortholog==TRUE)base::stop("ortholog option is only available for genomic_ressource class object from  ViSEAGO::EntrezGene2GO()")
  if(!base::is.null(id) && !base::is.character(id))base::stop("id must be a character value or NULL")

  ###################
  # Annotate
  ###################

    ###################
    # EntrezGene
    if(methods::slot(object,"db")=="EntrezGene"){

      ###################
      # for orthologs
      if(ortholog==TRUE){

        ###################
        # select the target species and ortholog
        gene_group<-ViSEAGO::EntrezGene_orthologs()

        ###################
        # extract Organisms informations
        taxon=ViSEAGO::taxonomy(base::unique(base::c(gene_group$tax_id, gene_group$Other_tax_id)))

        ###################
        # select the target species and ortholog
        if(base::is.null(id)){

          ###################
          # ordering by Scientific name
          data.table::setorder(taxon,ScientificName)

          ###################
          # stop scrip execution
          base::cat("Set id argument from ViSEAGO::annotate(id,EntrezGene,ortholog=TRUE) with an  available taxid (see below), and retry.",
          "Available EntrezGene species with orthologs_groups:\n\n",sep="\n")

          ###################
          # print taxon
          print(taxon,nrows=base::nrow(taxon))

          ###################
          # silencing stop
          opt <- base::options(show.error.messages=FALSE)
          base::on.exit(base::options(opt))
          stop()

        }else{

          ###################
          # check match id
          id=base::match.arg(id,taxon$taxid)

          ####################
          # select the target species and ortholog
          gene_group<-gene_group[relationship=="Ortholog" & (tax_id==id | Other_tax_id==id)]

          ###################
          # select the target species and ortholog from the left part of the table
          group1<-gene_group[tax_id==id,.(GeneID,Other_GeneID)]

          ###################
          # select the target species and ortholog from the right part of the table
          group2<-gene_group[Other_tax_id==id,.(Other_GeneID,GeneID)]

          ###################
          # renames(group2)
          base::names(group2)<-base::names(group1)

          ###################
          # bind and assign to gene_group (GeneID  from target species in first column)
          gene_group<-base::rbind(group1,group2)

          ###################
          # Extract otholog species annotation from data slot
          ortho<-methods::slot(object,"data")[evidence%in%c("EXP","IDA","IPI","IMP","IGI", "IEP")]

          ###################
          # merge experimental GO anotation from orthologs
          ortho<-merge(gene_group,ortho,by.x="Other_GeneID",by.y="gene_id")

          ###################
          # remove ortholog GeneID
          ortho[,Other_GeneID:=NULL]

          ###################
          # rename GeneID to gene_id
          base::names(ortho)[1]<-"gene_id"

          ###################
          # assign target species id to taxid and replac exprerimental evidence by IEA (computationnal)
          ortho[,`:=`(taxid=base::as.numeric(id),evidence="IEA")]

          ###################
          # Extract species annotation from data slot
          annot<-methods::slot(object,"data")[taxid==base::as.numeric(id)]

          ###################
          # add orthologs annotation to species annotation
          annot<-base::unique(base::rbind(annot,ortho))
        }

      }else{

        ###################
        # Extract species annotation from data slot
        annot<-methods::slot(object,"data")[taxid==base::as.numeric(id)]
      }

      ###################
      # Extract species annotation from data slot
      annot[category=="Function",category:="MF"]
      annot[category=="Process",category:="BP"]
      annot[category=="Component",category:="CC"]

      ###################
      # GO database stamp
      stamp=methods::slot(object,"stamp")
    }

    ###################
    # BioConductor
    if(methods::slot(object,"db")=="Bioconductor"){

      ###################
      # check id
      id=base::match.arg(
        id,
        methods::slot(ViSEAGO::Bioconductor2GO(),"organisms")$Package
      )

      ###################
      # install package if needed
      if(!id%in%utils::installed.packages()[,"Package"]){

        ###################
        # bioclite source
        BiocManager::install(id)
      }

      ###################
      # load GO annotations
      base::require(id,character.only=TRUE)

      ###################
      # load GO annotations
      annot<-data.table::data.table(AnnotationDbi::select(base::get(id),
      keys=AnnotationDbi::keys(base::get(id)),columns =c("ENTREZID","GO","ONTOLOGY")))

      ###################
      # rename columns
      base::names(annot)<-c("gene_id","GOID","evidence","category")

      ###################
      # GO database stamp
      stamp=base::eval(base::call(base::sub("\\.db","_dbInfo",id)))[15,2]
    }

    ###################
    # Ensembl
    if(methods::slot(object,"db")=="Ensembl"){

      ###################
      # connect to ensembl specified dataset
      myspecies<-biomaRt::useDataset(id,methods::slot(object,"mart")[[1]])

      ###################
      # with go name according biomart version
      go_filter<- data.table::data.table(biomaRt::listFilters(myspecies))

      ###################
      # load Ensembl genes, with GO annotations
      annot<-data.table::data.table(
        biomaRt::getBM(
          attributes =c("ensembl_gene_id","go_id","go_linkage_type","namespace_1003"),
          filters=go_filter[base::grep("with GO ID",ignore.case = T,go_filter$description),name],
          value=TRUE,
          mart =myspecies
        )
      )

      ###################
      # rename columns
      base::colnames(annot)<-c("gene_id","GOID","evidence","category")

      ###################
      # Extract species annotation from data slot
      annot[category=="molecular_function",category:="MF"]
      annot[category=="biological_process",category:="BP"]
      annot[category=="cellular_component",category:="CC"]

      ###################
      # GO database stamp
      stamp=methods::slot(object,"stamp")
    }

    ###################
    # Uniprot
    if(methods::slot(object,"db")=="Uniprot-GOA"){

      ###################
      # load the file
      utils::download.file(base::paste('ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/',
      base::toupper(id),'/goa_',id,'.gaf.gz',sep=""),destfile = "./data/input/annot.gz",quiet=TRUE)

      ###################
      # unzip
      R.utils::gunzip("./data/input/annot.gz")

      ###################
      # remove the zipped file
      base::unlink("./data/input/annot.gz")

      #################
      # read file
      annot<-base::unique(data.table::fread("./data/input/annot",quote="!",
      select=c(2,5,7,9),col.names=c("gene_id","GOID","evidence","category")))

      ###################
      # Extract species annotation from data slot
      annot[category=="F",category:="MF"]
      annot[category=="P",category:="BP"]
      annot[category=="C",category:="CC"]
    }

  ###################
  # convert to list
  ###################

    ###################
    # convert to list
    Data<-base::lapply(c("MF","BP","CC"),function(x){

      ###################
      # for a category
      Data<-annot[category==x]

      ###################
      # for each
      Dat<-base::lapply(base::unique(Data$gene_id),function(y){

        ###################
        # extract GOID by gene
        values<-Data[gene_id==y,GOID]

        ###################
        # add GO evidence
        base::names(values)<-Data[gene_id==y,evidence]

        ###################
        # return values
        values
      })
      names(Dat)<-base::unique(Data$gene_id)

      ###################
      # return values
      Dat
    })

  ###################
  # create GENE2GO object
  ###################

  methods::new("gene2GO",
    db=methods::slot(object,"db"),
    stamp=stamp,
    organism=id,
    MF=Data[[1]],
    BP=Data[[2]],
    CC=Data[[3]]
  )
})
