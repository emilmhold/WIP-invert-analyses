#ABMI Ambiguous allocation and supplementary scripts
#March 2024
#Ermias T. Azeria
#distribute parents to children based on look up

#The script uses species taxonomic look up table to assign ambiguous or unresolved taxon (parents) to best resolved taxonomy (children) for anlysis
#If analysis grouping is "Distribute to Children"#tmp.lookup$Analysis_Grouping=='Distribute.to.Children'
#Find the corresponding rank, and associated taxonomic name. tmp.lookup$Taxon
#Subset all the Analysis Group names in that taxonomic name

#Two main algorithms/functions: aptc_abmi and dpac_abmi
#aptc_abmi:  At any site having children, assign parent to the most numerous child (cf APTC.sg: Assign Parents To Children-SG )
#dpac_abmi: At any site having children, distribute the count of parent proportionally among children (cf DPAC.sg:Distribute Parents Among Children-SG: )
##If there is no child at a site, assign parents using regional data (regional level treatment) as defined below
#Arguments:
#data: main raw data with a mix of resolved and unresolved to be processed
#aux.data: auxiliary data with 'resolved' children data to be used for resolving taxonomy. For ABMI this refers to UMOS (unique mature organisms search)
#taxalookup: wetland invertebrate taxonomy table with information of parents that need be resolved and taxonomic level information of potential children
#The "Analysis_Grouping" column contains info whether the corrosponding taxonomic levelunder "Taxon" is to be resolved (value="Distribute.to.Children")

#method.region: Two options. "abundant"=allocate to the most abundant; "frequent"= allocate to the most widespread (cf Chao1 and Chao2, respctively)
#reg.prop: logical; if TRUE allocate proportional to argument specified in method.region, in proportion to abundance or frequency
#strata.primary and strata.secondary: when algorithm looks for children in the region, it searches for the abundant or frequent sequentially
#first in the primary spatial strata (eg.,natural subregions), then secondary (e.g., natural regions), and finally across all data
#Ensure strata information matches order of species data input
#Additional functionality
#Add weight (MCN): important for step 2 assignment based on abundance; now data is prepossessed outside the functions

#Additional scripts for preprocessing data
#When there is auxiliary data to be used, the row and columns should match to the main data.
#This is handled by the function "fn_matchdata"

#Main output: A list with the following tables
# "Out.main" = This contains the taxonomy resolved main data. Note this does not include the auxiliary data (i.e., UMOS).
# "Out.aux" = This is the taxonomy resolved data corresponding to the auxiliary data (i.e., UMOS). There were instances where the UMOS also had ambiguous parents.
# "Out.combined" = This data combines the resolved data for both main and auxiliary data tables.
# "Raw.main" = This is the raw main input data with both children and ambiguous parents.
# "Raw.aux" = This is the raw UMOS data.
# "Raw.child"= Children data in the raw data (original data).
# "New.child" = Children data due to assignment of ambiguous parent allocation to children. Note that the sum of raw and new child data gives the Out.main data.
# "Raw.child.aux" = Children data in the auxiliary (UMOS) original data.
# "New.child.aux" = Children data due to assignment of ambiguous parent to children fro auxiliary data (UMOS)
# "AP.info"= Number of children considered for each ambiguous parent.
# "Nochild.data" = Parents with no children. This keeps track of cases when there was no children assignment in the lookup table. This should sum to zero. If there is any parent with no children, then add/edit the lookup table.

#Additional diversity metric: "fn.commdiv" = computes abundance, richness, and diversity metrics

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
aptc_abmi <- function(data=data, aux.data=NULL, taxalookup=taxalookup, strata.primary=NULL, strata.secondary=NULL,method.region="abundant", reg.prop=FALSE, umos=NULL, mcn=NULL){
  require(data.table)

  if (is.null(aux.data)) {aux.raw= NULL} else {aux.raw=aux.data }

  if (!is.null(aux.data)) {
    All.taxa <- unique(c(colnames(data) , colnames(aux.data)))
    All.sites <- unique(c(rownames(data) , rownames(aux.data)))

    #Add missing taxa in main data but in auxiliary, or viceversa to respective data sets
    #with 0 abundance, to simplify downstream analysis
    #Handled now by "fn_matchdata"

    add.y <- setdiff(All.taxa, colnames(data))
    if(length(add.y )>0){
      y.df <- setNames(data.frame( matrix(0,ncol =length( add.y), nrow = nrow(data))), add.y )
      data <- data.frame(data, y.df)
          }
    #add sites if missing
    add.s <- setdiff(All.sites, rownames(data))
    if(length(add.s )>0){
      s.df <- setNames(data.frame( matrix(0, nrow = length( add.s),ncol =ncol(data))), colnames(data))
      rownames( s.df) <- add.s
      data <- rbind(data, s.df)
    }

    #do same for aux
    add.x <- setdiff(All.taxa, colnames(aux.data))
    if(length(add.x )>0){
      x.df <- setNames(data.frame( matrix(0,ncol =length( add.x), nrow = nrow(aux.data))), add.x )
      aux.data <- data.frame(aux.data, x.df)

    }
    add.sx <- setdiff(All.sites, rownames(aux.data))
    if(length(add.sx )>0){
      sx.df <-  setNames(data.frame( matrix(0, nrow = length( add.sx),ncol =ncol(data))), colnames(data))
      rownames( sx.df) <- add.sx
      aux.data <- rbind(aux.data, sx.df)

    }
    #Match sites and species order
    data <- data [ , All.taxa]
    data <- data [ All.sites,  ]
    aux.data <- aux.data[ , All.taxa]
    aux.data<- aux.data [ All.sites,  ]
  } else {
    aux.data <-  data
    aux.data [] <- 0
    All.taxa <- unique(c(colnames(data) , colnames(aux.data)))
     }

#
  All.taxa <- unique(c(colnames(data) , colnames(aux.data)))
  mm <-   match( All.taxa ,  taxalookup$Taxon_Analysis)
  if (any(is.na(mm))) { stop("Taxon_Analysis not found in taxa lookup")}

  tmp.lookup <- taxalookup[mm,  , drop=F]

  #IF UMOS data include in the look up

  ambiguous.parent <-    tmp.lookup$Taxon_Analysis[tmp.lookup$Analysis_Grouping=='Distribute.to.Children']
  ambiguous.rank <-    tmp.lookup$Rank[tmp.lookup$Analysis_Grouping=='Distribute.to.Children']#
  children.keep <- setdiff(All.taxa ,ambiguous.parent )

  #get children belonging to ambiguous parent, iterate among ambiguous group
  #Y.new: empty matrix to collate parent-to-child assignment

  Y.new <-  as.data.frame( matrix(0, nrow=nrow(data), ncol=length(children.keep)))
  rownames(Y.new) <- rownames(data)
  colnames (Y.new) <- children.keep

  #Ditto for original data
  X.orginal <- X.new <- Y.orginal <- Y.new
  comm.child <- intersect(colnames(data),children.keep )
  Y.orginal[ ,  comm.child ] <- data[ , comm.child]
  X.orginal[ ,  comm.child ] <- aux.data[ , comm.child]
  #keep track of no child parents
  Y.nochild <- data[ ,ambiguous.parent, drop=F ]
  Y.nochild [] <- 0
  #keep info for ambiguous parents
  AP.info <- NULL

  for (a in 1:length(ambiguous.parent)){
    ap.nam <-  ambiguous.parent[a]
    ap.rnk <- ambiguous.rank[a]
    ap.tax <- tmp.lookup[which(tmp.lookup$Taxon_Analysis==ap.nam), ap.rnk ]
    child.nam <- tmp.lookup$Taxon_Analysis[which(tmp.lookup[ ,ap.rnk]==ap.tax)]
    child.nam <- setdiff(child.nam ,ambiguous.parent )
    #subset data to taxonomic gp
    #IF no child or single child, return corresponding data
    ap.info <- c(ap.nam,length(child.nam ) )
    if (length(child.nam )<=1){

      if(length(child.nam )==1){ Y.new[ ,child.nam  ] <- Y.new[ ,child.nam   ] + data[ ,ap.nam ]
      X.new[ ,child.nam   ] <- X.new[ ,child.nam   ] + aux.data[ ,ap.nam ]
      }
      if(length(child.nam )==0){ Y.nochild[ ,ap.nam ] <-  data[ ,ap.nam , drop=F]}

    } else{

      Y <- Y.orginal[ ,child.nam ]+ X.orginal[ ,child.nam ]

      if(is.null(strata.primary)){strata.primary <- rep("All", nrow(Y))}

      spp.abu.reg  <- setDT(Y )[, lapply(.SD, sum), by=.( strata.primary), .SDcols=c(child.nam )]
      fn.nonz <- function(x){length (which(x>0))}
      spp.frq.reg  <- setDT(Y )[, lapply(.SD, fn.nonz) , by=.( strata.primary), .SDcols=c(child.nam )]
      setDF(Y)
      setDF(spp.frq.reg)
      rownames(spp.frq.reg) <- spp.frq.reg$strata.primary
      spp.frq.reg$strata.primary <- NULL
      setDF(spp.abu.reg)
      rownames(spp.abu.reg) <- spp.abu.reg$strata.primary
      spp.abu.reg$strata.primary <- NULL

      #If there is no child in strata , use overall frequent/abundant species
      y.abu <- colSums (Y)
      y.frq <- colSums (sign(Y))
      zz <- which(rowSums(spp.frq.reg)==0)
      strata.prim_0 <- names(zz )
      if(length(zz)>=0){
        spp.frq.reg [zz,] <- y.frq
        spp.abu.reg [zz, ] <- y.abu
      }
      #For proportional abu/frq
      spp.frq.reg.p <-  spp.frq.reg/ifelse(rowSums(spp.frq.reg)==0, 1,rowSums(spp.frq.reg))
      spp.abu.reg.p <-  spp.abu.reg/ifelse(rowSums(spp.abu.reg)==0, 1,rowSums(spp.abu.reg))

      #Which species has maximum abundant or frequent regional (primary strata)
      FrqSpp <- spp.frq.reg == do.call(pmax, spp.frq.reg)
      spp.frq.reg[FrqSpp] <- 1
      spp.frq.reg[!FrqSpp] <- 0
      spp.frq.reg <- spp.frq.reg/ifelse(rowSums(spp.frq.reg)==0, 1,rowSums(spp.frq.reg))#If species have a tie, share
      #for abundance
      AbuSpp <- spp.abu.reg == do.call(pmax, spp.abu.reg)
      spp.abu.reg[AbuSpp] <- 1
      spp.abu.reg[!AbuSpp] <- 0
      spp.abu.reg <- spp.abu.reg/ifelse(rowSums(spp.abu.reg)==0, 1,rowSums(spp.abu.reg))#If species have a tie, share

      #expand data to match species dimension
      spp.abu.reg <- spp.abu.reg[strata.primary,  ]
      spp.frq.reg <- spp.frq.reg[strata.primary,  ]

      spp.abu.reg.p <- spp.abu.reg.p[strata.primary,  ]
      spp.frq.reg.p <- spp.frq.reg.p[strata.primary,  ]

      ############
      #If there is secondary strata, calculate corresponding abundance and frequency
      if(!is.null(strata.secondary)){
      spp.abu.sec  <- setDT(Y )[, lapply(.SD, sum), by=.( strata.secondary), .SDcols=c(child.nam )]
      fn.nonz <- function(x){length (which(x>0))}#to obtain species incidence
      spp.frq.sec  <- setDT(Y )[, lapply(.SD, fn.nonz) , by=.( strata.secondary), .SDcols=c(child.nam )]
      setDF(Y)
      setDF(spp.frq.sec)
      rownames(spp.frq.sec) <- spp.frq.sec$strata.secondary
      spp.frq.sec$strata.secondary <- NULL
      setDF(spp.abu.sec)
      rownames(spp.abu.sec) <- spp.abu.sec$strata.secondary
      spp.abu.sec$strata.secondary <- NULL

   #If there is no child in secondary strata (abundance is zero in strata)
      # use overall frequent/abundant species
      y.abu <- colSums (Y)
      y.frq <- colSums (sign(Y))

        zz <- which(rowSums(spp.frq.sec)==0)
        strata.sec_0 <- names(zz )
        if(length(zz)>=0){
          spp.frq.sec [zz,] <- y.frq
          spp.abu.sec [zz, ] <- y.abu
        }
        #EA: in proportion to regional abundance or frequency
        spp.frq.sec.p <-  spp.frq.sec/ifelse(rowSums(spp.frq.sec)==0, 1,rowSums(spp.frq.sec))
        spp.abu.sec.p <-  spp.abu.sec/ifelse(rowSums(spp.abu.sec)==0, 1,rowSums(spp.abu.sec))


        #Which species has maximum abundant or frequent regional (primary strata)
        FrqSpp <- spp.frq.sec == do.call(pmax, spp.frq.sec)
        spp.frq.sec[FrqSpp] <- 1
        spp.frq.sec[!FrqSpp] <- 0
        spp.frq.sec <- spp.frq.sec/ifelse(rowSums(spp.frq.sec)==0, 1,rowSums(spp.frq.sec))#If species have a tie, share
        #for abundance
        AbuSpp <- spp.abu.sec == do.call(pmax, spp.abu.sec)
        spp.abu.sec[AbuSpp] <- 1
        spp.abu.sec[!AbuSpp] <- 0
        spp.abu.sec<- spp.abu.sec/ifelse(rowSums(spp.abu.sec)==0, 1,rowSums(spp.abu.sec))

        spp.abu.sec <- spp.abu.sec[strata.secondary,  ]
        spp.frq.sec <- spp.frq.sec[strata.secondary,  ]

        spp.abu.sec.p <- spp.abu.sec.p[strata.secondary,  ]
        spp.frq.sec.p <- spp.frq.sec.p[strata.secondary,  ]

        id.prim_0 <- which(strata.primary%in%strata.prim_0)

        spp.abu.reg[id.prim_0,  ]  <- spp.abu.sec[id.prim_0,  ]
        spp.frq.reg [id.prim_0,  ] <- spp.frq.sec[id.prim_0,  ]

        spp.abu.reg.p[id.prim_0,  ]  <- spp.abu.sec.p[id.prim_0,  ]
        spp.frq.reg.p [id.prim_0,  ] <- spp.frq.sec.p[id.prim_0,  ]

      }
      #~~~~~~~~~~~~~~~~~
      Y.max <- Y
      S.AbuSpp  <-  Y.max == do.call(pmax,  Y.max)
      Y.max[S.AbuSpp] <- 1
      Y.max[!S.AbuSpp] <- 0
      child.ties <- which (rowSums(Y.max)>1 & rowSums(Y) >0)

      nochild.site<- which (rowSums(Y)==0)

            if(method.region=="abundant") {
        reg.fill <- spp.abu.reg
        reg.fill.p <- spp.abu.reg.p
      }
      if(method.region=="frequent") {
        reg.fill <- spp.frq.reg
        reg.fill.p <- spp.frq.reg.p
      }
         if(reg.prop==TRUE){reg.fill <-reg.fill.p }else{reg.fill <-reg.fill}


      Y.max[nochild.site, ] <- reg.fill[nochild.site, ]
      Y.max.update <- Y.max
      #~~~
      if(length(child.ties)>0){
      Y.max.update [child.ties,] <- Y.max[child.ties,] * reg.fill[  child.ties, ]
      }

      yzz <- which (rowSums(Y.max.update)==0)

      if(length(yzz)>0){Y.max.update[yzz, ] <- Y.max[yzz, ]  }

      Y.max <-  Y.max.update /rowSums( Y.max.update)
      SiteAbu.spp <- Y.max * data[ ,ap.nam ]
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #Fix fraction abundance

      #~~~~~~~~~~

      fn_roundSum <- function(P.frac){
        N <- sum(P.frac)
        p.int <- round(P.frac- 0.005)
        fracs <- P.frac -  p.int
        #
        if(N - sum(p.int)>=1){
          topfrac <- tail(sort(fracs),round(sum(fracs)))
          add.id <- which( fracs%in%topfrac)[1:round(sum(fracs))]
          p.int[add.id] <- p.int[add.id] +1
        }
        p.int
      }


  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SiteAbu.spp <- t(apply( SiteAbu.spp, 1, fn_roundSum  ))
      Y.new[ ,colnames(SiteAbu.spp)] <-   Y.new[ ,colnames(SiteAbu.spp)]+SiteAbu.spp

       X.Abu.spp <- Y.max * aux.data[ ,ap.nam ]
       X.Abu.spp<- t(apply(  X.Abu.spp, 1, fn_roundSum  ))

    }
    AP.info <- rbind( AP.info, ap.info)
  }
  AP.info <- as.data.frame( AP.info )
  colnames(AP.info) <- c("Parent", "Number of children")
  rownames(AP.info) <- AP.info$Parent
  Y.out <- Y.new + Y.orginal
  Y.out  <- Y.out[, sort(colnames(Y.out))]
  X.out <- X.new + X.orginal
  X.out  <- X.out[, sort(colnames(X.out))]
  XY.out <- Y.out + X.out
  XY.out  <- XY.out[, sort(colnames(XY.out))]

  #add no child data to main output
  ncParents <- which(colSums(Y.nochild)>0)
  if(length(ncParents)>0){
    Y.out  <- data.frame( Y.out , Y.nochild[ , ncParents, drop=F])
    X.out  <- data.frame( X.out , Y.nochild[ , ncParents, drop=F])
    XY.out  <- data.frame( XY.out , Y.nochild[ , ncParents, drop=F])
  }


  if(is.null (aux.raw)) {X.out <- XY.out <-aux.data <- X.orginal <- X.new <- aux.raw}

  out <- list (Out.main=Y.out,Out.aux=X.out, Out.combined=XY.out , Raw.main=data, Raw.aux=aux.data,  Raw.child=Y.orginal, New.child= Y.new, Raw.child.aux=X.orginal, New.child.aux= X.new, AP.info = AP.info , Nochild.data= Y.nochild)
  out[sapply( out, is.null)] <- NULL
  out
}



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#dpac ABMI
#Assign parent proportional to child abundance in site
#where there is no child, assign to  more frequent or abundant
#EA consider also in proportion to respective regional abundance or frequency
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dpac_abmi <- function(data=data, aux.data=NULL, taxalookup, strata.primary=NULL, strata.secondary=NULL, method.region="frequent", reg.prop=FALSE, mcn=NULL ){
  require(data.table)

  if (is.null(aux.data)) {aux.raw= NULL} else {aux.raw=aux.data }#for output
  #if auxilary data is added, setup number of species and sites to match for original and filled data
  #Handled now by "fn_matchdata"

  if (!is.null(aux.data)) {
    All.taxa <- unique(c(colnames(data) , colnames(aux.data)))
    All.sites <- unique(c(rownames(data) , rownames(aux.data)))
    add.y <- setdiff(All.taxa, colnames(data))
    if(length(add.y )>0){
      y.df <- setNames(data.frame( matrix(0,ncol =length( add.y), nrow = nrow(data))), add.y )
      data <- data.frame(data, y.df)

    }
    #add sites if missing
    add.s <- setdiff(All.sites, rownames(data))
    if(length(add.s )>0){
      s.df <- setNames(data.frame( matrix(0, nrow = length( add.s),ncol =ncol(data))), colnames(data))
      rownames( s.df) <- add.s
      data <- rbind(data, s.df)

    }

    #do same for aux
    add.x <- setdiff(All.taxa, colnames(aux.data))
    if(length(add.x )>0){
      x.df <- setNames(data.frame( matrix(0,ncol =length( add.x), nrow = nrow(aux.data))), add.x )
      aux.data <- data.frame(aux.data, x.df)

    }
    add.sx <- setdiff(All.sites, rownames(aux.data))
    if(length(add.sx )>0){
      sx.df <-  setNames(data.frame( matrix(0, nrow = length( add.sx),ncol =ncol(data))), colnames(data))
      rownames( sx.df) <- add.sx
      aux.data <- rbind(aux.data, sx.df)

    }
    #Match sites and species order
    data <- data [ , All.taxa]
    data <- data [ All.sites,  ]
    aux.data <- aux.data[ , All.taxa]
    aux.data<- aux.data [ All.sites,  ]
  } else {
    aux.data <-  data
    aux.data [] <- 0
    All.taxa <- unique(c(colnames(data) , colnames(aux.data)))
  }

  #
  #
  All.taxa <- unique(c(colnames(data) , colnames(aux.data)))
  mm <-   match( All.taxa ,  taxalookup$Taxon_Analysis)
  if (any(is.na(mm))) { stop("Taxon_Analysis not found in taxa lookup")}

  tmp.lookup <- taxalookup[mm,  , drop=F]

  ambiguous.parent <-    tmp.lookup$Taxon_Analysis[tmp.lookup$Analysis_Grouping=='Distribute.to.Children']
  ambiguous.rank <-    tmp.lookup$Rank[tmp.lookup$Analysis_Grouping=='Distribute.to.Children']#
  children.keep <- setdiff(All.taxa ,ambiguous.parent )

  #Set up empty matrix to store results

  Y.new <-  as.data.frame( matrix(0, nrow=nrow(data), ncol=length(children.keep)))
  rownames(Y.new) <- rownames(data)
  colnames (Y.new) <- children.keep
  #add same dimension original data

  X.orginal <- X.new <- Y.orginal <- Y.new
  comm.child <- intersect(colnames(data),children.keep )
  Y.orginal[ ,  comm.child ] <- data[ , comm.child]
  X.orginal[ ,  comm.child ] <- aux.data[ , comm.child]
  #keep track of no child parents
  Y.nochild <- data[ ,ambiguous.parent, drop=F ]
  Y.nochild [] <- 0

  AP.info <- NULL #keep info for ambiguous parents
  #iterate among ambiguous group
  for (a in 1:length(ambiguous.parent)){
    ap.nam <-  ambiguous.parent[a]
    ap.rnk <- ambiguous.rank[a]
    ap.tax <- tmp.lookup[which(tmp.lookup$Taxon_Analysis==ap.nam), ap.rnk ]
    child.nam <- tmp.lookup$Taxon_Analysis[which(tmp.lookup[ ,ap.rnk]==ap.tax)]
    child.nam <- setdiff(child.nam ,ambiguous.parent )
    ap.info <- c(ap.nam,length(child.nam ) )
    if (length(child.nam )<=1){
      if(length(child.nam )==1){ Y.new[ ,child.nam   ] <- Y.new[ ,child.nam   ] + data[ ,ap.nam ]
      X.new[ ,child.nam   ] <- X.new[ ,child.nam   ] + aux.data[ ,ap.nam ]
                   }
      if(length(child.nam )==0){ Y.nochild[ ,ap.nam ] <-  data[ ,ap.nam , drop=F]}


         } else {

       Y <- Y.orginal[ ,child.nam ]+ X.orginal[ ,child.nam ]

      if(is.null(strata.primary)){strata.primary <- rep("All", nrow(Y))}

      spp.abu.reg  <- setDT(Y )[, lapply(.SD, sum), by=.( strata.primary), .SDcols=c(child.nam )]
      fn.nonz <- function(x){length (which(x>0))}
      spp.frq.reg  <- setDT(Y )[, lapply(.SD, fn.nonz) , by=.( strata.primary), .SDcols=c(child.nam )]
      setDF(Y)
      setDF(spp.frq.reg)
      rownames(spp.frq.reg) <- spp.frq.reg$strata.primary
      spp.frq.reg$strata.primary <- NULL
      setDF(spp.abu.reg)
      rownames(spp.abu.reg) <- spp.abu.reg$strata.primary
      spp.abu.reg$strata.primary <- NULL

      #overall abundant and frequent
      FrqSpp.reg <- colnames(spp.frq.reg)[apply(spp.frq.reg,1,which.max)]
      AbuSpp.reg <- colnames(spp.abu.reg)[apply(spp.abu.reg,1,which.max)]

      #If child abundance is zero in strata , use overall frequent/abundant species

      y.abu <- colSums (Y)
      y.frq <- colSums (sign(Y))
      zz <- which(rowSums(spp.frq.reg)==0)
      strata.prim_0 <- names(zz ) #keep track of primary strata with no children, to fill with secondary strata info when available
      if(length(zz)>=0){
        spp.frq.reg [zz,] <- y.frq
        spp.abu.reg [zz, ] <- y.abu
      }
      #Proportion to regional abundance or frequency
      spp.frq.reg.p <-  spp.frq.reg/ifelse(rowSums(spp.frq.reg)==0, 1,rowSums(spp.frq.reg))
      spp.abu.reg.p <-  spp.abu.reg/ifelse(rowSums(spp.abu.reg)==0, 1,rowSums(spp.abu.reg))


      FrqSpp <- spp.frq.reg == do.call(pmax, spp.frq.reg) #get child with max abundance or frquency
      spp.frq.reg[FrqSpp] <- 1
      spp.frq.reg[!FrqSpp] <- 0
      spp.frq.reg <- spp.frq.reg/rowSums(spp.frq.reg) #If children have a tie, share

      AbuSpp <- spp.abu.reg == do.call(pmax, spp.abu.reg)
      spp.abu.reg[AbuSpp] <- 1
      spp.abu.reg[!AbuSpp] <- 0
      spp.abu.reg<- spp.abu.reg/rowSums(spp.abu.reg)
      #expand data to match species dimension
      spp.abu.reg <- spp.abu.reg[strata.primary,  ]
      spp.frq.reg <- spp.frq.reg[strata.primary,  ]

      spp.abu.reg.p <- spp.abu.reg.p[strata.primary,  ]
      spp.frq.reg.p <- spp.frq.reg.p[strata.primary,  ]

      ############
      #If there is secondary strata, calculate corresponding abundance and frequency
      if(!is.null(strata.secondary)){
        spp.abu.sec  <- setDT(Y )[, lapply(.SD, sum), by=.( strata.secondary), .SDcols=c(child.nam )]
        fn.nonz <- function(x){length (which(x>0))}#to obtain species incidence
        spp.frq.sec  <- setDT(Y )[, lapply(.SD, fn.nonz) , by=.( strata.secondary), .SDcols=c(child.nam )]
        setDF(Y)
        setDF(spp.frq.sec)
        rownames(spp.frq.sec) <- spp.frq.sec$strata.secondary
        spp.frq.sec$strata.secondary <- NULL
        setDF(spp.abu.sec)
        rownames(spp.abu.sec) <- spp.abu.sec$strata.secondary
        spp.abu.sec$strata.secondary <- NULL

        zz <- which(rowSums(spp.frq.sec)==0)
        strata.sec_0 <- names(zz ) #keep track of primary strata with no children, to fill with secondary strata info when available
        if(length(zz)>=0){
          spp.frq.sec [zz,] <- y.frq
          spp.abu.sec [zz, ] <- y.abu
        }
        #EA: in proportion to regional abundance or frequency
        spp.frq.sec.p <-  spp.frq.sec/ifelse(rowSums(spp.frq.sec)==0, 1,rowSums(spp.frq.sec))
        spp.abu.sec.p <-  spp.abu.sec/ifelse(rowSums(spp.abu.sec)==0, 1,rowSums(spp.abu.sec))


        FrqSpp <- spp.frq.sec == do.call(pmax, spp.frq.sec)
        spp.frq.sec[FrqSpp] <- 1
        spp.frq.sec[!FrqSpp] <- 0
        spp.frq.sec <- spp.frq.sec/rowSums(spp.frq.sec)#If species have a tie, share
        #for abundance
        AbuSpp <- spp.abu.sec == do.call(pmax, spp.abu.sec)
        spp.abu.sec[AbuSpp] <- 1
        spp.abu.sec[!AbuSpp] <- 0
        spp.abu.sec<- spp.abu.sec/rowSums(spp.abu.sec)

        #expand data to match species dimension
        spp.abu.sec <- spp.abu.sec[strata.secondary,  ]
        spp.frq.sec <- spp.frq.sec[strata.secondary,  ]

        spp.abu.sec.p <- spp.abu.sec.p[strata.secondary,  ]
        spp.frq.sec.p <- spp.frq.sec.p[strata.secondary,  ]

        id.prim_0 <- which(strata.primary%in%strata.prim_0)

        spp.abu.reg[id.prim_0,  ]  <- spp.abu.sec[id.prim_0,  ]
        spp.frq.reg [id.prim_0,  ] <- spp.frq.sec[id.prim_0,  ]

        spp.abu.reg.p[id.prim_0,  ]  <- spp.abu.sec.p[id.prim_0,  ]
        spp.frq.reg.p [id.prim_0,  ] <- spp.frq.sec.p[id.prim_0,  ]

      }

      #If no children in strata, use study-wide

      if(method.region=="abundant") {
        reg.fill <- spp.abu.reg
        reg.fill.p <- spp.abu.reg.p
        }
      if(method.region=="frequent") {
        reg.fill <- spp.frq.reg
        reg.fill.p <- spp.frq.reg.p
        }
      #Proportional use?
       if(reg.prop==TRUE){reg.fill <-reg.fill.p }else{reg.fill <-reg.fill}
       SiteProp.spp <- Y/ifelse(rowSums(Y)==0, 1,rowSums(Y))

      nochild.site<- which (rowSums(Y)==0)
      SiteProp.spp[nochild.site, ] <- reg.fill[nochild.site, ]
      SiteAbu.spp <- SiteProp.spp * data[ ,ap.nam ]
      SiteAbu.spp <- as.data.frame(SiteAbu.spp )
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #Fix fractional abundance data
      #~~~~~~~~~~
      fn_roundSum <- function(P.frac){
        N <- sum(P.frac)
        p.int <- round(P.frac- 0.005)
        fracs <- P.frac -  p.int
        #
        if(N - sum(p.int)>=1){
          topfrac <- tail(sort(fracs),round(sum(fracs)))
          add.id <- which( fracs%in%topfrac)[1:round(sum(fracs))]
          p.int[add.id] <- p.int[add.id] +1
        }
        p.int
      }
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SiteAbu.spp <- t(apply( SiteAbu.spp, 1, fn_roundSum  ))

      #Add child to data

      Y.new[ ,colnames(SiteAbu.spp)] <-   Y.new[ ,colnames(SiteAbu.spp)]+SiteAbu.spp

      #Resolve data for aux

      X.Abu.spp <- SiteProp.spp  * aux.data[ ,ap.nam ]
      X.Abu.spp<- t(apply(  X.Abu.spp, 1, fn_roundSum  ))

    }
    AP.info <- rbind( AP.info, ap.info)
  }
  AP.info <- as.data.frame( AP.info )
  colnames(AP.info) <- c("Parent", "Number of children")
  rownames(AP.info) <- AP.info$Parent
  Y.out <- Y.new + Y.orginal
  Y.out  <- Y.out[, sort(colnames(Y.out))]
  X.out <- X.new + X.orginal
  X.out  <- X.out[, sort(colnames(X.out))]
  XY.out <- Y.out + X.out
  XY.out  <- XY.out[, sort(colnames(XY.out))]

  #add no child data to main output
  ncParents <- which(colSums(Y.nochild)>0)
  if(length(ncParents)>0){
    Y.out  <- data.frame( Y.out , Y.nochild[ , ncParents, drop=F])
    X.out  <- data.frame( X.out , Y.nochild[ , ncParents, drop=F])
    XY.out  <- data.frame( XY.out , Y.nochild[ , ncParents, drop=F])
      }

  if(is.null (aux.raw)) {X.out <- XY.out <-aux.data <- X.orginal <- X.new <- aux.raw}

  out <- list (Out.main=Y.out,Out.aux=X.out, Out.combined=XY.out , Raw.main=data, Raw.aux=aux.data,  Raw.child=Y.orginal, New.child= Y.new, Raw.child.aux=X.orginal, New.child.aux= X.new, AP.info = AP.info , Nochild.data= Y.nochild)
  out[sapply( out, is.null)] <- NULL
  out
  }

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#~~~~~~~
#Function to match main and auxiliary data

fn_matchdata <- function(data=data, aux.data=aux.data){

    All.taxa <- unique(c(colnames(data) , colnames(aux.data)))
    All.sites <- unique(c(rownames(data) , rownames(aux.data)))

    #Add missing taxa in main data but in auxiliary, or viceversa, to repective data sets
    #with 0 abundance, to simplify downstream analysis

    add.y <- setdiff(All.taxa, colnames(data))
    if(length(add.y )>0){
      y.df <- setNames(data.frame( matrix(0,ncol =length( add.y), nrow = nrow(data))), add.y )
      data <- data.frame(data, y.df)

    }
    #add sites if missing
    add.s <- setdiff(All.sites, rownames(data))
    if(length(add.s )>0){
      s.df <- setNames(data.frame( matrix(0, nrow = length( add.s),ncol =ncol(data))), colnames(data))
      rownames( s.df) <- add.s
      data <- rbind(data, s.df)

    }

    #do same for aux
    add.x <- setdiff(All.taxa, colnames(aux.data))
    if(length(add.x )>0){
      x.df <- setNames(data.frame( matrix(0,ncol =length( add.x), nrow = nrow(aux.data))), add.x )
      aux.data <- data.frame(aux.data, x.df)

    }
    add.sx <- setdiff(All.sites, rownames(aux.data))
    if(length(add.sx )>0){
      sx.df <-  setNames(data.frame( matrix(0, nrow = length( add.sx),ncol =ncol(data))), colnames(data))
      rownames( sx.df) <- add.sx
      aux.data <- rbind(aux.data, sx.df)

    }
    data <- data [ , All.taxa]
    data <- data [ All.sites,  ]
    aux.data <- aux.data[ , All.taxa]
    aux.data<- aux.data [ All.sites,  ]
    out <- list ( data= data, aux.data=aux.data)
    out
}

#~~~~~~~~~~
#CALCULATE RICHNESS AND DIVERSITY FOR THE ABOVE FINAL DATA

#Most frequent/abundant in region, and based on APTC or DPAC
#For main data and with auxiliary (i.e., UMOS added)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Function to compute richness and diversity measures
#Consider adding  functionality to add coarse data (non-targetd) data
#Align rows between input and coarse data
fn.commdiv <- function (data.in=data.in, coarse.data=NULL){
  require(vegan)
  #check if umos were included
  tmpdf <- data.in$Out.combined
  if(!is.null(tmpdf)){
    Total.abu <- rowSums( tmpdf)
    Taxa.rich <- rowSums(sign(tmpdf))
    ShanonDiv <- diversity(tmpdf)
    ShannonEq <- ShanonDiv/ log(Taxa.rich)
    SimpsonInv <- diversity(tmpdf, index="invsimpson")
    SimpsonInv  [is.infinite(SimpsonInv )] <- 0
    #Ditto for main data only
    tmpdf <- data.in$Out.main
    Total.abu.m <- rowSums( tmpdf)
    Taxa.rich.m <- rowSums(sign(tmpdf))
    ShanonDiv.m <- diversity(tmpdf)
    ShannonEq.m <- ShanonDiv.m/ log(Taxa.rich.m)
    SimpsonInv.m <- diversity(tmpdf, index="invsimpson")
    SimpsonInv.m [is.infinite(SimpsonInv.m )] <- 0
    out<- data.frame(Total.abu,Taxa.rich, ShanonDiv, ShannonEq, SimpsonInv, Total.abu.m ,Taxa.rich.m, ShanonDiv.m, ShannonEq.m, SimpsonInv.m)
  }else{
    tmpdf <- data.in$Out.main
    Total.abu.m <- rowSums( tmpdf)
    Taxa.rich.m <- rowSums(sign(tmpdf))
    ShanonDiv.m <- diversity(tmpdf)
    ShannonEq.m <- ShanonDiv/ log(Taxa.rich)
    SimpsonInv.m <- diversity(tmpdf, index="invsimpson")
    SimpsonInv.m [is.infinite(SimpsonInv.m )] <- 0 #for inverse simpson
    out<- data.frame(Total.abu.m ,Taxa.rich.m, ShanonDiv.m, ShannonEq.m, SimpsonInv.m)
  }
  out[is.na(out)] <- 0
  out <- round(out,3)

  return(out)

}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



#~~~~~~~~~~~~~~~
#function to fix child  site abundance that are high or low than parent due to rounding

fn.roundfix <- function(ff, csa){
  nonz <- which(csa >0)
  if(length(nonz )==0){nonz <- 1:length(csa)
  ss <- sample(nonz,abs(ff), replace=F )
  ffs <- ff/length(ss)
  ffs <- rep(ffs, length(ss))
  csa[ss] <- csa[ss]-ffs
  } else{
    ffs<- round((csa[nonz] /sum(csa[nonz]))*ff)
    csa[nonz ] <- csa[nonz ]-ffs
  }

  csa
}
##~~~~~~~~~~~~~
#use above function to fix values
# if (length(fix.sites) > 0){
#   for (f in fix.sites) {
#     SiteAbu.spp[f, ] <- fn.roundfix (fix.child[f],  as.numeric(SiteAbu.spp[f, ,drop=T] ))
#   }
# }


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
