options(scipen=999)
#install.packages("/.../isr_1.0.tar.gz", repos=NULL, type="source")
library(isr)

# lines 63 and 65 need to be updated depending on what link function you are using,
# and for what period you want the slopes
############################################################################################################################################################
############################################################################################################################################################

# updating functions from the ISR package 

summary.of.output<-function(list.of.models, dataset, model.type, resp.variable, expl.variables.without.segm.variable, segm.variable, model.parameters, angle.param=0)
{
  #reconstruct the best model
  final<-list.of.models[[length(list.of.models)]]$best.model
  temp<-fit.model.with.breakpoints(dataset=dataset, model.type=model.type, resp.variable=resp.variable,	
                                   expl.variables.without.segm.variable=expl.variables.without.segm.variable, segm.variable=segm.variable, 
                                   breakpoints=final$transf.breakpoints, model.parameters=model.parameters, angle.param=angle.param)
  best.model<-temp$model
  best.dataset<-temp$dataset
  
  temp.colnames<-c("breakpoint", "left.CI", "right.CI", "when.noticed")
  summary.dataset<-as.data.frame(t(data.frame(rep(NA, length(temp.colnames)))))
  colnames(summary.dataset)<-temp.colnames
  #original angle
  temp.breakpoint<-min(dataset[,segm.variable])
  summary.dataset<-rbind(summary.dataset, c(temp.breakpoint, rep(NA, length(temp.colnames)-1)))
  
  #if at least one breakpoint was found, calculate confidence intervals and detection times
  if (length(list.of.models)>1) {
    for (i in 1:(length(list.of.models)-1))
    {	
      m<-list.of.models[[i]]
      #breakpoint
      number.of.breakpoints<-length(m$best.model$transf.breakpoints)
      temp.breakpoint<-m$best.model$transf.breakpoints[number.of.breakpoints]
      #confidence intervals
      all.breakpoints<-temp.breakpoint
      for (m2 in m$good.breakpoint.models)
        all.breakpoints<-append(all.breakpoints, (m2$transf.breakpoints[number.of.breakpoints]))
      temp.left.CI<-min(all.breakpoints)
      temp.right.CI<-max(all.breakpoints)
      #when noticed
      temp.when.noticed<-m$end.point
      #collate the data together
      summary.dataset<-rbind(summary.dataset, c(temp.breakpoint, temp.left.CI, temp.right.CI, temp.when.noticed))
    }
  }
  summary.dataset<-summary.dataset[2:nrow(summary.dataset),]
  for (c in c("breakpoint", "left.CI", "right.CI", "when.noticed"))
    summary.dataset[,c]<-as.POSIXct(summary.dataset[,c], origin="1970-01-01")
  #95% CI for the segm.variable coefficients
  best.coeffs<-as.data.frame(temp$summary.of.output[(nrow(temp$summary.of.output)-length(final$transf.breakpoints)):nrow(temp$summary.of.output), ])
  if (identical(final$transf.breakpoints, character(0)))
    best.coeffs<-as.data.frame(t(best.coeffs))
  best.coeffs$left.CI<-best.coeffs[,1]-1.96*best.coeffs[,2]
  best.coeffs$right.CI<-best.coeffs[,1]+1.96*best.coeffs[,2]
  summary.dataset$angle<-best.coeffs[,1]
  summary.dataset$angle.left.CI<-best.coeffs[,5]
  summary.dataset$angle.right.CI<-best.coeffs[,6]
  summary.dataset$pvalue<-best.coeffs[,4]
  ##calculate angle per year exp(angle*scaleby)
  scaleby<-3600*24*365.25
  for (c in c("angle", "angle.left.CI", "angle.right.CI"))
    summary.dataset[[c]]<-exp(summary.dataset[[c]]*scaleby)
  return(list(best.model=best.model, dataset=best.dataset, summary.dataset=summary.dataset))
}	

plot.fit<-function(dataset, model, summary.dataset, 
                   original.resp.variable, resp.variable.label, 
                   factor.expl.variables, cts.expl.variables, offset.variable="",
                   segm.variable, segm.variable.label, lambda, ylim.max)
{	
  #replace all continuous variables with their mean value
  if (!(identical(cts.expl.variables, "")))
    for (c in cts.expl.variables)
      dataset[,c]<-mean(dataset[,c], na.rm=TRUE)
    
    #if there are factor variables, generate all their possible combinations
    if (!identical(factor.expl.variables, "") & !identical(factor.expl.variables, character()))
    { 
      factors<-list()
      for (f in factor.expl.variables)
        factors[[f]]<-levels(dataset[[f]])
      factors.df<-expand.grid(factors)
      factors.df.long<-factors.df[rep(seq(nrow(factors.df)), nrow(dataset)),]
      if (length(factor.expl.variables)==1)
      {
        factors.df.long<-as.data.frame(factors.df.long)
        colnames(factors.df.long)<-factor.expl.variables
      }
      
      #create several copies of the dataset: one for each possible combination
      dataset<-dataset[rep(seq(nrow(dataset)), nrow(factors.df)),]
      
      for (f in factor.expl.variables)
        dataset[[f]]<-factors.df.long[[f]]
      
      #record all possible combinations of factor levels
      dataset$factors<-eval(parse(text=paste("paste(", paste("dataset[['", factor.expl.variables, "']]", sep="", collapse=","), ")", sep="")))
      
    }
    #add predict column
    dataset$predict<-predict(model, dataset, type="response")
    dataset<-untransf.resp.variable(dataset=dataset, var="predict", lambda=lambda)$output.dataset	
    
    dataset$yvar<-as.numeric(as.character(dataset[,original.resp.variable]))	
    dataset$xvar<-dataset[,segm.variable]
    
    #scale by offset.variable if necessary
    if (!identical(offset.variable, ""))
    {
      dataset$untransf.predict<-dataset$untransf.predict/dataset[[offset.variable]]
      dataset$yvar<-dataset$yvar/dataset[[offset.variable]]
    }
    
    #plot
    p<-ggplot(dataset, aes(x=xvar))
    #historic data
    p<-p+geom_point(aes(y=yvar,alpha = 1/1000))
    p<-p+scale_alpha(guide = 'none')
    p<-p+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line.x = element_line(color="black", size = 1),
               axis.line.y = element_line(color="black", size = 1),aspect.ratio=4/4)
    p<-p+theme(plot.margin=unit(c(1,1,1,1),"cm"))
    p<-p+ylim(0,ylim.max)
    #predicted values
    #colour by factors if needed
    if (!identical(factor.expl.variables, "") & !identical(factor.expl.variables, character()))
      p<-p+geom_line(aes(y=untransf.predict, colour=factors, group=factors), size=0.5) else
        #        p<-p+geom_line(aes(y=untransf.predict), size=0.5, colour="blue")
        p<-p+geom_line(aes(y=untransf.predict), size=1, colour="black")
    #change x-axis
    p<-p+xlab(segm.variable.label)
    p<-p+theme(text=element_text(family="Arial"),
               axis.text.x=element_text(size=10,colour="black"),
               axis.title.x=element_text(size=10,vjust=-1),
               axis.text.y=element_text(size=10,colour="black"),
               axis.title.y=element_text(size=10),
               legend.text=element_text(size=10))
   #change x-axis
    p<-p+xlab(segm.variable.label)
   #change y-axis
    p<-p+ylab(resp.variable.label)	
    #remove ticks
#    p<-p+theme(axis.ticks=element_blank())
    #if at least one breakpoint was found, plot breakpoints and their confidence intervals and detection times
    if (nrow(summary.dataset)>1)
    {
      #remove first line as it corresponds to the change at the start
      summary.dataset2<-summary.dataset[2:nrow(summary.dataset),]
      #y-coordinate for breakpoints
      summary.dataset2$ycoord<-min(dataset$yvar)
      #breakpoints
      p<-p+geom_point(data=summary.dataset2, aes(x=breakpoint, y=0), colour="blue") 
      #confidence intervals for breakpoints
      p<-p+geom_errorbarh(data=summary.dataset2, aes(xmin=left.CI, xmax=right.CI, x=breakpoint, y=0), colour="blue", height=0)
      #breakpoints detection times
      #p<-p+geom_point(data=summary.dataset2, aes(x=when.noticed, y=0), colour="blue", shape=3)		
    }
    return(list(p=p))
}


master.isr<-function(dataset, 
                     resp.variable, resp.variable.label,  
                     step=FALSE, t.boxcox, 
                     factor.expl.variables, cts.expl.variables, offset.variable="", segm.variable, segm.variable.label, 
                     model.type, model.parameters, angle.param=0,
                     start.point, end.point,  
                     time.unit, interval.length, min.distance,
                     criterion, criterion.difference, ylim.max)
{
  #transform resp.variable using boxcox or identity
  transformation<- if (t.boxcox) "boxcox" else "identity"
  temp<-transf.resp.variable(dataset=dataset, var=resp.variable, transformation=transformation)
  new.dataset<-temp$output.dataset
  new.resp.variable<-temp$transf.var
  lambda<-temp$lambda
  
  #generate explanatory variables (excluding the segmented variable)
  expl.variables.without.segm.variable<-character()
  expl.variables.without.segm.variable<-append(expl.variables.without.segm.variable, factor.expl.variables)
  
  # add step if needed
  if (!identical(step, FALSE))
  { temp<-add.step(dataset=new.dataset, var=resp.variable, time.var=segm.variable, step=step) 
  new.dataset<-temp$output.dataset
  expl.variables.without.segm.variable<-append(expl.variables.without.segm.variable, temp$step.var)
  }
  expl.variables.without.segm.variable<-append(expl.variables.without.segm.variable, cts.expl.variables)
  
  #add offset if needed
  if ((!identical(offset.variable, "")) & (!identical(offset.variable, character())))
  {	log.offset.variable<-paste("offset(log(", offset.variable, "))", sep="")
  expl.variables.without.segm.variable<-append(expl.variables.without.segm.variable, log.offset.variable)
  }
  
  #remove all ""
  expl.variables.without.segm.variable<-expl.variables.without.segm.variable[which(expl.variables.without.segm.variable!='')]
  
  #perform iterative sequential regression
  isr.models<-sequential.grid.search.save.models(dataset=new.dataset, model.type=model.type, resp.variable=new.resp.variable, 
                                                 expl.variables.without.segm.variable=expl.variables.without.segm.variable, 
                                                 segm.variable=segm.variable, 
                                                 model.parameters=model.parameters, angle.param=angle.param,
                                                 start.point=start.point, end.point=end.point, 
                                                 min.distance=min.distance, 
                                                 interval.length=interval.length, time.unit=time.unit,
                                                 criterion=criterion, criterion.difference=criterion.difference)
  
  #summary of the output
  output<-summary.of.output(dataset=new.dataset, list.of.models=isr.models$models, 
                            model.type=model.type, resp.variable=new.resp.variable, 
                            expl.variables.without.segm.variable=expl.variables.without.segm.variable, segm.variable=segm.variable, model.parameters=model.parameters,
                            angle.param=angle.param)
  
  #graph
  require(ggplot2)
  plot.fit<-plot.fit(dataset=output$dataset, model=output$best.model, summary.dataset=output$summary.dataset, 
                     original.resp.variable=resp.variable, resp.variable.label=resp.variable.label, lambda=lambda,
                     factor.expl.variables=factor.expl.variables, cts.expl.variables=cts.expl.variables, 
                     offset.variable=offset.variable, segm.variable=segm.variable,
                     segm.variable.label=segm.variable.label, ylim.max=ylim.max)	
  
  return(list(list.of.models=isr.models$models, output=output, plot.fit=plot.fit))
}

monthstart <- function(x) 
{
  x <- as.POSIXlt(x)
  x$mday <- 1
  as.Date(x)
}

format_pvalues <- function(x)
{
  if(x < 0.001)
    return (paste("p<",0.001,sep=""))
  else
    if(x < 0.0095)
      return (paste("p=",round(x,digits=3),sep=""))
    else
      return (paste("p=",round(x,digits=2),sep=""))
}

test <- c(0.3898675,0.04987,0.007876,0.000032432)
for(i in 1:length(test))
  print(format_pvalues(test[i]))

format_CI <- function(x,y)
  return (paste("95% CI=(",round(x,digits=2),",",round(y,digits=2),")",sep=""))
format_CI(0.1,0.2)

format_IRR <- function(x)
  return(paste("IRR=",round(x,digits=2),sep=""))
format_IRR(0.5)

format_CM <- function(x)
  return(paste("CM=",round(x,digits=2),sep=""))
format_CM(0.557)
