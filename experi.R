library(rattle)
library(AUC)
library(pROC)

experi <- function(form, ds, dsname, target, modeller, details="", 
                   n=100, control=NULL,
                   keep=FALSE, # Keep the last model built.
                   prob="prob",
                   class="class",
                   log="experi.log")
{
  suppressPackageStartupMessages(require(pROC))
  
  user <- Sys.getenv("LOGNAME")
  node <- Sys.info()[["nodename"]]
  
  wsrpart.model <- modeller=="wsrpart"
  
  numclass <- length(levels(ds[,target]))
  
  start.time <- proc.time()
  
  seeds <- cors <- strs <- aucs <- accs <- NULL
  for (i in seq_len(n))
  {
    loop.time <- proc.time()
    
    seeds  <- c(seeds, seed <- sample(1:1000000, 1))
    set.seed(seed)
    
    train  <- sample(nrow(ds), 0.7*nrow(ds))
    test   <- setdiff(1:nrow(ds), train)
    actual <- ds[test, target]
    
    args   <- append(list(form, data=ds[train,]), control)
    model  <- do.call(modeller, args)
    
    if (numclass==2)
    {
      if (modeller %in% c("ctree", "cforest"))
        pr <- do.call("rbind", predict(model, newdata=ds[test,], type=prob))[,2]
      else
        pr <- predict(model, newdata=ds[test,], type=prob)[,2]
      # For small samples we can get all the same class...
      # Should really ensure we sample both classes
      if (length(unique(actual)) == 1)
        aucs <- c(0.0, aucs)
      else
        aucs <- c(auc(actual, pr), aucs)
    }
    if ("ada" %in% class(model)) class <- "vector"
    if (modeller %in% c("ctree", "cforest")) class <- "response"
    
    #compute cirrelation and strength for mrpart
    if (wsrpart.model)
    {
      cors <- c(correlation(model, ds, form), cors)
      strs <- c(strength(model, ds, form), strs)
    }
    
    cl <- predict(model, newdata=ds[test,], type=class)
    accs <- c(sum(cl==actual, na.rm=TRUE)/length(actual), accs)
    
    if (! is.null(log))
    {
      require(lubridate)
      if (! file.exists(log))
        write.table(data.frame(user=NA, node=NA, ts=NA, ds=NA, model=NA,
                               acc=NA, auc=NA, cor=NA, str=NA,
                               user=NA, elapsed=NA),
                    file=log, sep=",", row.names=FALSE)
      write.table(data.frame(user=user, node=node, ts=now(),
                             ds=dsname, model=modeller,
                             acc=round(accs[1], 4),
                             auc=round(aucs[1], 4),
                             cor=ifelse(wsrpart.model, round(cors[1], 4), NA),
                             str=ifelse(wsrpart.model, round(strs[1], 4), NA),
                             user=round((proc.time()-loop.time)['user.self'], 2),
                             elapsed=round((proc.time()-loop.time)['elapsed'],2)),
                  file=log, sep=",", append=TRUE, col.names=FALSE, row.names=FALSE)
    }
  }
  
  result <- data.frame(modeller=paste0(modeller, "_", details),
                       auc=ifelse(numclass==2, mean(aucs), NA),
                       auc.sd=ifelse(numclass==2, sd(aucs), NA),
                       cor=ifelse(wsrpart.model, mean(cors), NA),
                       cor.sd=ifelse(wsrpart.model , sd(cors), NA),
                       str=ifelse(wsrpart.model, mean(strs), NA),
                       str.sd=ifelse(wsrpart.model , sd(strs), NA),
                       acc=mean(accs), acc.sd=sd(accs), n=n,
                       user=(proc.time()-start.time)['user.self'],
                       elapsed=(proc.time()-start.time)['elapsed'])
  if (wsrpart.model)
    if (numclass==2)
      result[-1]   <- round(result[-1], 2)
  else
    result[-c(1:3)]   <- round(result[-c(1:3)], 2)
  else
    if (numclass==2)
      result[-c(1,4:7)]   <- round(result[-c(1,4:7)], 2)
  else
    result[-c(1:7)]   <- round(result[-c(1:7)], 2)
  
  row.names(result) <- NULL
  
  if (keep)
  {
    if (numclass==2) 
    {
      attr(result, "pr") <- pr
      attr(result, "test") <- test
    }
    attr(result, "model") <- model
  }
  
  return(result)
}