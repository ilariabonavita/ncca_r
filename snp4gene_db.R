library(DBI)


# the argument is a txt file with this format:
# RSID SNPID position genes
# genes can be 'NA', a string with the name of one gene, a string with the names of more genes separeted by '_'

init.map <- function(rds = '/Users/ilaria_bonavita/statgen/dyslexia/00_data/rs.pos.3.txt'){
  return(read.table(rds,h=T))
}


setup.db <- function(dbname = 'map'){
  ## create empty DB
  mydb <- dbConnect(RSQLite::SQLite(), paste0(dbname,".sqlite"))
  
  ## write pure SNP info (first three columns) to table 'snps'
  dbWriteTable(mydb, 'snps', map[,1:3])
  
  ## create empty 2-col table 'sgmap'
  dbExecute(mydb, 'CREATE TABLE sgmap (snp char16, gene char16);')
  
  ## fill that
  for (i in 1:nrow(map)){
    gene.list <- strsplit(as.character(map[i, 'genes']), '_')[[1]]
    #print(gene.list)
    for (gene in gene.list){
      sql.cmd <- sprintf("INSERT INTO sgmap VALUES ('%s', '%s')", map[i, 'RSID'], gene)
      # print(sql.cmd)
      dbExecute(mydb, sql.cmd)
    }
  }
  return(mydb)
}

snps4gene <- function(db, gene){
  dbGetQuery(db, sprintf("SELECT * FROM sgmap WHERE gene = '%s'", gene))
}


#---- Usage example
#
#map <- init.map()
#mydb <- setup.db(dbname = 'map')
#snp.gene <- snps4gene(mydb,gene)

# Query examples
#gene.list <- dbGetQuery(mapdb, 'SELECT gene FROM sgmap GROUP BY gene HAVING count(gene) > 3 ORDER BY count(gene) DESC')
#gene.count <- dbGetQuery(mydb,"SELECT gene, count(gene) AS cnt FROM sgmap GROUP BY gene HAVING cnt > 3 ORDER BY cnt DESC")
