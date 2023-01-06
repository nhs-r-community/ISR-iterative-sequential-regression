bc.transform<-function(var)
{	
	require(MASS)
	bc.transf<-boxcox(var~1, lambda=seq(-5, 5, 1/10))
	#lambda is the value that maximises the log-likelihood
	lambda<-bc.transf$x[which.max(bc.transf$y)]
	#if lambda=0, transformation is log(x)
	#if lambda is non-zero, transformation applied is (x^lambda-1)/lambda
	transf.var<-if (lambda==0) log(var) else (var^lambda-1)/lambda
	return(list(lambda=lambda, transf.var=transf.var))
}

bc.untransform<-function(transf.var, lambda)
{
	#if lambda=0, transformation is exp(x)
	#if lambda is non-zero, transformation applied is (lambda*y+1)^(1/lambda)
	untransf.var<-if (lambda==0) exp(transf.var) else (lambda*transf.var+1)^(1/lambda)
	return(list(untransf.var=untransf.var))
}

transf.resp.variable<-function(dataset, var, transformation)
{ 
	output.dataset<-dataset
	transf.var<-paste("transf", var, sep=".")
	if (transformation=="boxcox")
	{
		temp<-bc.transform(var=dataset[[var]])
		output.dataset[,transf.var]<-temp$transf.var
		lambda<-temp$lambda
	} else {
		output.dataset[,transf.var]<-dataset[,var]
		lambda<-"na"
	}
	return(list(output.dataset=output.dataset, lambda=lambda, transf.var=transf.var))
}

untransf.resp.variable<-function(dataset, var, lambda)
{ 
	output.dataset<-dataset
	untransf.var<-paste("untransf", var, sep=".")
	output.dataset[,untransf.var]<-if (lambda=="na") dataset[,var] else bc.untransform(transf.var=dataset[,var], lambda=lambda)$untransf.var
	return(list(output.dataset=output.dataset, untransf.var=untransf.var))
}

add.step<-function(dataset, var, time.var, step)
{ 
	output.dataset<-dataset
	step.var<-paste(var, "step", sep=".")
	#the step variable takes value 0 before the step change and 1 afterwards
	output.dataset[,step.var]<- ifelse (output.dataset[, time.var]<step, 0, 1)
	return(list(step.var=step.var, output.dataset=output.dataset))
}

next.set.to.test.change.last.coordinate.only<-function(current.set, number.of.intervals, min.distance)
{ 
	n<-length(current.set)
	#if can't move any further up, return "max permitted"
	if (current.set[n]+1+min.distance>number.of.intervals) next.set<-"max permitted" else
	{
		#otherwise increase the last coordinate by 1
		next.set<-current.set
		next.set[n]<-next.set[n]+1
	}
	return(list(next.set=next.set))
}

add.extra.coordinate<-function(current.set, number.of.intervals, min.distance)
{
	n<-length(current.set)
	#read off last coordinate of the current vector
	current.last.coordinate<- if (length(current.set)==0) 0 else current.set[n]
	#the next coordinate must be at least min.distance away
	next.last.coordinate<-current.last.coordinate+min.distance
	#check whether we are still able to move up; if so, adjoin the newly found last coordinate to the current set
	if (next.last.coordinate+min.distance>number.of.intervals) next.set<-"max permitted" else
		next.set<-append(current.set, next.last.coordinate)
	return(list(next.set=next.set))
}

transform.breakpoint<-function(start.point, breakpoint, interval.length, time.unit)
{ 
	if (identical(breakpoint, numeric(0))) transf.breakpoint<-character(0) else
	{ 
		distance<-breakpoint*interval.length
		transf.breakpoint<-lapply(distance, function(x) seq(from=as.Date(start.point), length.out=2, by=paste(x, time.unit))[2])
		transf.breakpoint<-as.POSIXct(as.Date(unlist(transf.breakpoint), origin="1970-01-01"))
	}
	return(list(transf.breakpoint=transf.breakpoint))
}

break.segm.variable<-function(angle.param=0, dataset, segm.variable, breakpoints)
{ 	
	output.dataset<-dataset
	#initialise the first variable, initially equal to the time variable itself
	output.dataset[,paste("char", segm.variable, 0, sep=".")]<- output.dataset[,segm.variable]
	split.variables<-paste("char", segm.variable, 0, sep=".")
	if (angle.param==0) #so that we are interested in changes between angles
	{ if (length(breakpoints)>0)
			for (j in 1:length(breakpoints))
			{ 
				#0-th (initial) time variable remains unchanged
				#j-th time variable equals to 0 before the j-th breakpoint, and difference between current time and j-th breakpoint otherwise 
				output.dataset[,paste("char", segm.variable, j, sep=".")]<-
						ifelse(output.dataset[,segm.variable]<breakpoints[j],
								difftime(dataset[,segm.variable], dataset[,segm.variable]),
								difftime(dataset[,segm.variable], breakpoints[j], units="secs"))
				
				split.variables<-append(split.variables, paste("char", segm.variable, j, sep="."))
			}
	} else { #so that we are interested in new angles
		if(length(breakpoints)>=1)  
		{	
			#0-th (initial) time variable equal to current time before the 1st breakpoint, and 1st breakpoint afterwards 
			output.dataset[,paste("char", segm.variable, 0, sep=".")]<- pmin(output.dataset[,segm.variable], breakpoints[1])
			for (j in 1:length(breakpoints))
			{ 
				lower.bound<-breakpoints[j]
				upper.bound<-if (j<length(breakpoints)) breakpoints[j+1] else max(dataset[,segm.variable])
				#j-th time variable equals to 0 before j-th breakpoint, (j+1)st breakpoint after the (j+1)st breakpoint, and difference between current time and j-th breakpoint in between 
				output.dataset[,paste("char", segm.variable, j, sep=".")]<-
						difftime(pmin(pmax(output.dataset[,segm.variable], lower.bound), upper.bound), lower.bound)
				split.variables<-append(split.variables, paste("char", segm.variable, j, sep="."))
			}
		}
	}	
	return(list(output.dataset=output.dataset, split.variables=split.variables))
}

fit.model<-function(dataset, model.type, resp.variable, expl.variables, model.parameters)
{ 
	#reconstruct the formula for modelling
	formula.to.model<-paste(paste(resp.variable, collapse='+'), "~", paste(expl.variables, collapse='+'))
	if (model.parameters!='' && !is.null(model.parameters)) formula.to.model<-paste(formula.to.model, ",", model.parameters)
	formula.to.model<-paste(model.type, "(", formula.to.model, ",", "data=dataset", ")")
	model<-eval(parse(text=formula.to.model)) 
	return(list(coefficients=summary(model, se="iid")$coefficients[,1], significance=summary(model, se="iid")$coefficients[,4], 
					formula.to.model=formula.to.model, model=model))
}

fit.model.with.breakpoints<-
		function(dataset, model.type, resp.variable, expl.variables.without.segm.variable, segm.variable, breakpoints, model.parameters, angle.param=0)
{	
	#split the time variable into several, according to the breakpoints
	#NB as before, angle.param=0 if we are interested in change between the new angle and the old angle, and 1 if we are interested in the new angle  
	temp<-break.segm.variable(angle.param=angle.param, dataset=dataset, segm.variable=segm.variable, breakpoints=breakpoints)
	dataset<-temp$output.dataset
	extra.expl.variables<-temp$split.variable
	#add the new time variables to the list of response variables
	expl.variables<-append(expl.variables.without.segm.variable, extra.expl.variables)  
	expl.variables<-expl.variables[which(expl.variables!="")]
	#fit the model, taking the new time variables into account
	fit<-fit.model(dataset=dataset, model.type=model.type, resp.variable=resp.variable, expl.variables=expl.variables, model.parameters=model.parameters)
	
	#reproduce model coefficients: estimates, p-values etc.
	summary.of.output<-summary(fit$model, se="iid")$coefficients
	#time variables coefficients
	#estimates
	segm.var.coeffs<-c()
	for (v in extra.expl.variables)
		segm.var.coeffs<-append(segm.var.coeffs, fit$coefficients[v])
	#p-values
	segm.var.significance<-c()
	for (v in extra.expl.variables)
		segm.var.significance<-append(segm.var.significance, fit$significance[v])
		
	return(list(dataset=dataset, coefficients=fit$coefficients, significance=fit$significance, 
			segm.var.coeffs=segm.var.coeffs, segm.var.significance=segm.var.significance, formula.to.model=fit$formula.to.model, model=fit$model,
			summary.of.output=summary.of.output))
}

find.smallest<-function(list.of.models, criterion)
{ 
	value<-list.of.models[[1]][[criterion]]
	index<-1
	for (i in 1:length(list.of.models))
		if (list.of.models[[i]][[criterion]]<value)
		{ index<-i
			value<-list.of.models[[i]][[criterion]]
		}
	return(list(value=value, index=index))
} 

breakpoint.models.fixed.interval<-function(dataset, model.type, 
		resp.variable, expl.variables.without.segm.variable, segm.variable, 
		model.parameters, angle.param=0, fixed.breakpoints,
		start.point, end.point, min.distance, interval.length, time.unit)
{	
	number.of.intervals<-length(seq(from=start.point, to=end.point, by=paste(interval.length, time.unit)))-1
	#generate all possible breakpoints that include the fixed breakpoints
	#start with a "straight line" - ie only the fixed breakpoints, and no change in trend afterwards
	list.of.breakpoints<-list(straight.line.breakpoints=fixed.breakpoints)
	#change in trend assumed, so need to add an extra breakpoint
	next.breakpoints<-add.extra.coordinate(current.set=fixed.breakpoints, number.of.intervals=number.of.intervals, 
			min.distance=min.distance)$next.set
	#add all possible sets of breakpoints constaining the fixed breakpoints plus an extra breakpoint between the last fixed breakpoint & current end point 
	while(!identical(next.breakpoints,"max permitted"))
	{
		list.of.breakpoints<-append(list.of.breakpoints, list(next.breakpoints))
		next.breakpoints<-
				next.set.to.test.change.last.coordinate.only(current.set=next.breakpoints, 
						number.of.intervals=number.of.intervals, min.distance=min.distance)$next.set
	}		
	#for each set of breakpoints, fit a model
	list.of.models<-list()	
	for (b in list.of.breakpoints)
	{ 		
		#transform breakpoints into the POSIXct format
		transf.breakpoints<-
			transform.breakpoint(start.point=start.point, breakpoint=b, 
			interval.length=interval.length, time.unit=time.unit)$transf.breakpoint
		#fit a model corresponding to the current breakpoints
	  	temp.fit<-fit.model.with.breakpoints(dataset=dataset, model.type=model.type, resp.variable=resp.variable, 
			expl.variables.without.segm.variable=expl.variables.without.segm.variable, segm.variable=segm.variable, 
			breakpoints=transf.breakpoints, model.parameters=model.parameters, angle.param=angle.param)
		#for each model its breakpoints, time var coefficients, AIC and BIC are recorded
	 	list.of.models[[paste("breakpoints=", b, collapse=',')]]<-list(number.of.breakpoints=length(b), 
			breakpoints=b, transf.breakpoints=transf.breakpoints, 
			segm.var.coeffs=temp.fit$segm.var.coeffs, segm.var.significance=temp.fit$segm.var.significance,
			aic=AIC(temp.fit$model), bic=AIC(temp.fit$model, k=log(nrow(dataset))), 
			formula.to.model=temp.fit$formula.to.model)
	}
	#save the "straight line" model
	#will then be compared with all "change of trend" models
	straight.line.model<-list.of.models[[paste("breakpoints=",fixed.breakpoints, collapse=',')]]
	
	return(list(list.of.breakpoints=list.of.breakpoints,  
			list.of.models=list.of.models, straight.line.model=straight.line.model))
}

search.for.new.breakpoint.fixed.interval<-function(dataset, model.type, 
		resp.variable, expl.variables.without.segm.variable, segm.variable, 
		model.parameters, angle.param=0, fixed.breakpoints,
		start.point, end.point, min.distance,
		interval.length, time.unit,
		criterion, criterion.difference)
{
	#fit the "straight line" and the "change in trend" models
	alpha<-breakpoint.models.fixed.interval(dataset=dataset, model.type=model.type,
			resp.variable=resp.variable, 
			expl.variables.without.segm.variable=expl.variables.without.segm.variable, 
			segm.variable=segm.variable, 
			model.parameters=model.parameters, angle.param=angle.param, fixed.breakpoints=fixed.breakpoints,
			start.point=start.point, end.point=end.point, 
			min.distance=min.distance, 
			interval.length=interval.length, time.unit=time.unit)
	#find the best model, ie the one with the smallest AIC/BIC
	best.model.index<-find.smallest(alpha$list.of.models, criterion)$index
	best.model<-alpha$list.of.models[[best.model.index]]
	#check whether the model is statistically better than the straight line one
	#"statistically better" means that its BIC is more than criterion.difference away from the BIC of the straight line model	
	#if so, we conclude that there is a change in trend on the interval 
	found.new.breakpoint<-FALSE
	new.fixed.breakpoints<-fixed.breakpoints
	#record all models that are statistically equivalent to the best model
	#will only be relevant (non-empty) if a change in trend is found
	good.breakpoint.models<-list()
	if (best.model[[criterion]]<alpha$straight.line.model[[criterion]]-criterion.difference) #so that the best model is statistically better than the straight line one
	{	
		#flag that a new breakpoint had been found
		found.new.breakpoint<- TRUE
		#fix its time
		new.fixed.breakpoints<-best.model$breakpoints 
		#record all models that are statistically equivalent to the best model
		for (m in alpha$list.of.models)
			if (m[[criterion]]-best.model[[criterion]]<criterion.difference)
				good.breakpoint.models<-append(good.breakpoint.models, list(m))
	}
	return(list(best.model=best.model, 
			good.breakpoint.models=good.breakpoint.models,
			straight.line.model=alpha$straight.line.model,
			found.new.breakpoint=found.new.breakpoint, new.fixed.breakpoints=new.fixed.breakpoints))
}

sequential.search.for.new.breakpoint.whole.interval<-function(dataset, model.type, 
		resp.variable, expl.variables.without.segm.variable, segm.variable, 
		model.parameters, angle.param=0, fixed.breakpoints, last.end.point,
		start.point, end.point, min.distance, 
		interval.length, time.unit,
		criterion, criterion.difference)
{
	
	last.fixed.breakpoint<-if (identical(fixed.breakpoints, numeric())) 0 else fixed.breakpoints[length(fixed.breakpoints)]
	found.new.breakpoint<-FALSE
	#determine next endpoint; needs to be at least 2*min.distance away from the last fixed point
	temp.end.point<-max(last.end.point, last.fixed.breakpoint+2*min.distance)	
	#transform the endpoint into the POSIXct format 
	transf.temp.end.point<-
			transform.breakpoint(start.point=as.POSIXct(start.point), breakpoint=temp.end.point, 
			interval.length=interval.length,
			time.unit=time.unit)$transf.breakpoint

	#only consider the data between the start.point and the temp.end.point	
	temp.dataset<-eval(parse(text=paste("subset(dataset, ", segm.variable, " <=transf.temp.end.point)")))
	good.breakpoint.models<-list()
	#sequential search
	#on each interval compare the "straight line" model and the "change in trend" models
	while (transf.temp.end.point<=end.point)
	{		
		alpha<-search.for.new.breakpoint.fixed.interval(dataset=temp.dataset, 
			model.type=model.type, 
			resp.variable=resp.variable, 
			expl.variables.without.segm.variable=expl.variables.without.segm.variable, 
			segm.variable=segm.variable, 
			model.parameters=model.parameters, angle.param=angle.param, fixed.breakpoints=fixed.breakpoints,
			start.point=start.point, end.point=transf.temp.end.point, min.distance=min.distance, 
			interval.length=interval.length,
			time.unit=time.unit,
			criterion=criterion, criterion.difference=criterion.difference)
		if (alpha$found.new.breakpoint==TRUE)
		#stop if found new breakpoint 
		{
			found.new.breakpoint<-TRUE
			new.best.model<-alpha$best.model
			new.best.model.dataset<-temp.dataset
			#record the time when the change was detected
			new.end.point<-temp.end.point
			#add the new breakpoint to the list of fixed breakpoint
			new.fixed.breakpoints<-alpha$new.fixed.breakpoints
			#record the models that are statistically equivalent to the best one
			good.breakpoint.models<-append(good.breakpoint.models, alpha$good.breakpoint.models)
			break
		}
		#if there is no evidence for a change in trend at this interval (ie no new breakpoint found), move forward and repeat
		temp.end.point<-temp.end.point + 1
		transf.temp.end.point<-
				transform.breakpoint(start.point=start.point, breakpoint=temp.end.point, 
				interval.length=interval.length, time.unit=time.unit)$transf.breakpoint
		temp.dataset<-eval(parse(text=paste("subset(dataset, ", segm.variable, " <transf.temp.end.point)")))
	}
	#if reached the end without finding any extra breakpoints, there is no evidence of change of trend apart from the ones that have already been found
	#so the final model is the model with the breakpoints that were fixed initially and straight line afterwards
	if (found.new.breakpoint==FALSE)
	{
		new.end.point<-temp.end.point
		transf.breakpoints<-
				transform.breakpoint(start.point=start.point, breakpoint=fixed.breakpoints, 
					interval.length=interval.length, time.unit=time.unit)$transf.breakpoint
		new.best.model.dataset<-dataset
		temp.fit<-fit.model.with.breakpoints(dataset=new.best.model.dataset, 
				model.type=model.type, resp.variable=resp.variable, 
				expl.variables.without.segm.variable=expl.variables.without.segm.variable, segm.variable=segm.variable, 
				breakpoints=transf.breakpoints, model.parameters=model.parameters, angle.param=angle.param)
		new.best.model<-list(number.of.breakpoints=length(fixed.breakpoints), 
				breakpoints=fixed.breakpoints, transf.breakpoints=transf.breakpoints, 
				segm.var.coeffs=temp.fit$segm.var.coeffs, segm.var.significance=temp.fit$segm.var.significance,
				aic=AIC(temp.fit$model), bic=AIC(temp.fit$model, k=log(nrow(dataset))), 
				formula.to.model=temp.fit$formula.to.model)	
		good.breakpoint.models<-list(new.best.model)
		new.fixed.breakpoints<-fixed.breakpoints
	}		
	return(list(found.new.breakpoint=found.new.breakpoint, new.end.point=new.end.point,
			new.fixed.breakpoints=new.fixed.breakpoints,
			new.best.model.dataset=new.best.model.dataset,
			new.best.model=new.best.model, good.breakpoint.models=good.breakpoint.models))
}

sequential.grid.search.save.models<-function(dataset, model.type, 
		resp.variable, expl.variables.without.segm.variable, segm.variable, 
		model.parameters, angle.param=0, start.point, end.point, min.distance, 
		interval.length, time.unit,
		criterion, criterion.difference)
{
	models<-list()
	fixed.breakpoints<-numeric()
	last.end.point<-1
	found.new.breakpoint<-TRUE
	while(found.new.breakpoint)
	{
		alpha<-sequential.search.for.new.breakpoint.whole.interval(dataset=dataset,
			model.type=model.type, resp.variable=resp.variable, 
			expl.variables.without.segm.variable=expl.variables.without.segm.variable, 
			segm.variable=segm.variable, 
			model.parameters=model.parameters, angle.param=angle.param, 
			fixed.breakpoints=fixed.breakpoints, 
			last.end.point=last.end.point,
			start.point=start.point, end.point=end.point, 
			min.distance=min.distance, 
			interval.length=interval.length, time.unit=time.unit,
			criterion=criterion, criterion.difference=criterion.difference)
	    found.new.breakpoint<-alpha$found.new.breakpoint
		fixed.breakpoints<-alpha$new.fixed.breakpoints
		last.end.point<-alpha$new.end.point
		models<-append(models, list(list(best.model=alpha$new.best.model, good.breakpoint.models=alpha$good.breakpoint.models,
						start.point=start.point, 
						end.point=transform.breakpoint(start.point=start.point, 
								breakpoint=alpha$new.end.point, 
								interval.length=interval.length, 
								time.unit=time.unit)$transf.breakpoint)))		
	}
	return(list(models=models))
}

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
	##calculate angle per year
	#scaleby<-3600*24*365
	#for (c in c("angle", "angle.left.CI", "angle.right.CI"))
	#	summary.dataset[[c]]<-summary.dataset[[c]]*scaleby
	
	return(list(best.model=best.model, dataset=best.dataset, summary.dataset=summary.dataset))
		
}	
	
plot.fit<-function(dataset, model, summary.dataset, 
		original.resp.variable, resp.variable.label, 
		factor.expl.variables, cts.expl.variables, offset.variable="",
	      segm.variable, segm.variable.label, lambda)
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
	p<-p+geom_point(aes(y=yvar))
	#predicted values
	#colour by factors if needed
	if (!identical(factor.expl.variables, "") & !identical(factor.expl.variables, character()))
		p<-p+geom_line(aes(y=untransf.predict, colour=factors, group=factors), size=0.5) else
			p<-p+geom_line(aes(y=untransf.predict), size=0.5, colour="blue")
	#change x-axis
	p<-p+xlab(segm.variable.label)
	p<-p+opts(axis.text.x=theme_text(angle=45))
	#change y-axis
	p<-p+ylab(resp.variable.label)	
	#remove ticks
	p<-p+opts(axis.ticks=theme_blank())
	#if at least one breakpoint was found, plot breakpoints and their confidence intervals and detection times
	if (nrow(summary.dataset)>1)
	{
		#remove first line as it corresponds to the change at the start
		summary.dataset2<-summary.dataset[2:nrow(summary.dataset),]
		#y-coordinate for breakpoints
		summary.dataset2$ycoord<-min(dataset$yvar)
		#breakpoints
		p<-p+geom_point(data=summary.dataset2, aes(x=breakpoint, y=ycoord), colour="blue") 
		#confidence intervals for breakpoints
		p<-p+geom_errorbarh(data=summary.dataset2, aes(xmin=left.CI, xmax=right.CI, x=breakpoint, y=ycoord), colour="blue", height=0)
		#breakpoints detection times
		p<-p+geom_point(data=summary.dataset2, aes(x=when.noticed, y=ycoord), colour="blue", shape=3)		
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
	criterion, criterion.difference)
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
			segm.variable.label=segm.variable.label)	
		
	return(list(list.of.models=isr.models$models, output=output, plot.fit=plot.fit))
}