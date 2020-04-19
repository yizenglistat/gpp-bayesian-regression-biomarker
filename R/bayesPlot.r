bayesPlot <- function(path='./output/', plot_beta=FALSE, plot_sigma=FALSE,plot_ise=FALSE, which_beta=1, model=1, sizes=sizes,
	true_params=list(sigma_true=sigma_true,model_list=model_list),
	plot_config = list(line_col=c("#28B463",'#A93226'), band_col='gray', band_opacity=0.3)){
	
	debug <- FALSE
	if(debug){
		path='./output/';
		plot_beta=FALSE;
		plot_sigma=FALSE;
		which_beta=1;
		model=1;
		sizes=c(2,4);
		true_params=list(sigma_true=sigma_true,model_list=model_list);
		plot_config = list(line_col=c("#28B463",'#A93226'), band_col='gray', band_opacity=0.3)
	}

	sigma_true 	<- true_params$sigma_true
	beta_list 	<- model_list[[paste0('model',model)]]
	nbeta 		<- length(beta_list)

	if(substr(path,nchar(path),nchar(path))!='/') path <- paste0(path,'/')

	# require pacakges
	require(lattice)
	require(latticeExtra)
	require(data.table)
	graphics.off()

	# read/load data automatically for beta_knots_summ and theta_summ
	readFiles <- function(path,model,size,homo=FALSE,beta_list){
		path <- paste0(path,'model',model,'/size',size,'/')
		if(homo) path <- paste0(path,'homo/')
		files <- list.files(path)
		
		beta_knots_names 	<- c()
		beta_knots_total 	<- 0
		beta_knots_count 	<- 0

		theta_names 		<- c()
		theta_total 		<- 0
		theta_count 		<- 0
		theta_median 		<- c()
		beta_ise 		<- c()

		for(filename in files){
			if(grepl('beta_knots_summ',filename)){
				name <- unlist(strsplit(filename,'.csv'))
				filepath <- paste0(path, filename)
				assign(name, as.matrix(fread(filepath)))
				beta_knots_total <- beta_knots_total + get(name)
				t_knots <- get(name)[,1]
				beta_true <- betaTrue(t_knots,beta_list,center_intc=FALSE)
				beta_ise <- rbind(beta_ise,apply((get(name)[,4:5]-beta_true)^2,2,mean))
				beta_knots_names <- c(beta_knots_names, name)
				beta_knots_count <- beta_knots_count + 1
			}
			if(grepl('theta_summ',filename)){
				name <- unlist(strsplit(filename,'.csv'))
				filepath <- paste0(path, filename)
				assign(name, as.matrix(fread(filepath)))
				theta_total <- theta_total + get(name)
				theta_median<- cbind(theta_median,get(name)[1,c(2,4,5)])
				theta_names <- c(theta_names, name)
				theta_count <- theta_count + 1
			}
		}

		beta_knots_avg 		<- beta_knots_total / beta_knots_count
		theta_avg 			<- theta_total / theta_count
		if(homo){
			colnames(theta_avg) 	<- paste0(colnames(theta_avg),'_homo')
			theta_avg 				<- theta_avg[,order(ncol(theta_avg):1)]
		}
		return(list(beta_knots_avg=beta_knots_avg,beta_ise=beta_ise,theta_avg=theta_avg,theta_median=theta_median))
	}

	createThetaView <- function(theta_avg,theta_homo_avg){

		theta_true 			<- c(sigma_true,	
								rep(NA,2*nbeta))

		theta_res <- data.frame(theta_homo_avg,true=theta_true, theta_avg,
								row.names=c('sigma',
											paste0('tau',1:nbeta),
											paste0('phi',1:nbeta)))
		return(theta_res[,order(ncol(theta_res):1)])
	}
	
	createThetaTable <- function(theta_med,theta_med_homo,type){
		nkeep <- ncol(theta_med)
		sigma_res <- data.frame(true=sigma_true,
								theta = theta_med[1,],
								lower = theta_med[2,],
								upper = theta_med[3,],
								type=type,label=rep('Random',each=nkeep),varname='sigma')
		colnames(sigma_res) <- c('true','theta','lower','upper','type','label','varname')
		sigma_homo_res <- data.frame(true=sigma_true,
								theta = theta_med_homo[1,],
								lower = theta_med_homo[2,],
								upper = theta_med_homo[3,],
								type=type,label=rep('Homogeneous',each=nkeep),varname='sigma')
		colnames(sigma_homo_res) <- c('true','theta','lower','upper','type','label','varname')
		sigma_res <- rbind(sigma_res,sigma_homo_res)
		sigma_res$type <- as.factor(sigma_res$type)
		sigma_res$label <- as.factor(sigma_res$label)

		theta_res <- sigma_res

		return(theta_res)
	}

	createBetaTable <- function(beta_knots_avg,beta_knots_homo_avg,type,which_beta,beta_list=beta_list){
		# beta varying function visualization
		t_knots 			<- beta_knots_avg[,'t_knots']
		nknots 				<- length(t_knots)
		beta_true 			<- betaTrue(t_knots, beta_list, center_intc=FALSE)
		beta_index 			<- seq(2,nbeta*5,nbeta) + which_beta
		beta_res <- data.frame(t_knots=t_knots,
							beta_true=beta_true[,which_beta+1],
							beta_knots_avg[,beta_index],
							type=type,
							label=rep('Random',each=nknots))
		colnames(beta_res) <- c('t','true','mean','median','sd','lower','upper','type','label')

		beta_homo_res <- data.frame(t_knots=t_knots,
							beta_true=beta_true[,which_beta+1],
							beta_knots_homo_avg[,beta_index],
							type=type,
							label=rep('Homogeneous',each=nknots))
		colnames(beta_homo_res) <- c('t','true','mean','median','sd','lower','upper','type','label')
		
		beta_res <- rbind(beta_res, beta_homo_res)
		beta_res$type <- as.factor(beta_res$type)
		beta_res$label <- as.factor(beta_res$label)

		return(beta_res)
	}

	createBetaISE <- function(beta_ise,beta_ise_homo, type){
		nbeta <- ncol(beta_ise)
		nkeep <- nrow(beta_ise)
		beta_ise_res <- data.frame(ise=as.vector(beta_ise),
							varname=rep(paste0('beta',(1:nbeta)-1),each=nkeep),
							type=type,
							label='Random')
		colnames(beta_ise_res) <- c('ise','varname','type','label')
		beta_ise_homo_res <- data.frame(ise_homo=as.vector(beta_ise_homo),
							varname=rep(paste0('beta',(1:nbeta)-1),each=nkeep),
							type=type,
							label='Homogeneous')
		colnames(beta_ise_homo_res) <- c('ise','varname','type','label')
		beta_ise_res <- rbind(beta_ise_res, beta_ise_homo_res)
		beta_ise_res$type <- as.factor(beta_ise_res$type)
		beta_ise_res$label <- as.factor(beta_ise_res$label)
		beta_ise_res$varname <- factor(beta_ise_res$varname, levels=unique(paste0('beta',(1:nbeta)-1)))
		return(beta_ise_res)
	}

	theta_vars <- c()
	theta_homo_vars <- c()
	
	theta_med_vars <- c()
	theta_med_homo_vars <- c()

	beta_knots_vars <- c()
	beta_knots_homo_vars <- c()

	beta_ise_vars <- c()
	beta_ise_homo_vars <- c()

	theta_res<- c()
	beta_res <- c()
	beta_ise_res <- c()

	for(size in sizes){
			
		theta_vars <- c(theta_vars,paste0('theta_size',size))
		theta_homo_vars <- c(theta_homo_vars,paste0('theta_homo_size',size))

		theta_med_vars <- c(theta_med_vars,paste0('theta_med_size',size))
		theta_med_homo_vars <- c(theta_med_homo_vars,paste0('theta_med_homo_size',size))
		
		beta_knots_vars <- c(beta_knots_vars,paste0('beta_knots_size',size))
		beta_knots_homo_vars <- c(beta_knots_homo_vars,paste0('beta_knots_homo_size',size))
		
		beta_ise_vars <- c(beta_ise_vars,paste0('beta_ise_size',size))
		beta_ise_homo_vars <- c(beta_ise_homo_vars,paste0('beta_ise_homo_size',size))

		out_random <- readFiles(path,model=model,size=size,homo=FALSE,beta_list=beta_list)
		out_homo <- readFiles(path,model=model,size=size,homo=TRUE,beta_list=beta_list)
		assign(tail(theta_vars,1),out_random$theta_avg)
		assign(tail(theta_homo_vars,1),out_homo$theta_avg)

		assign(tail(theta_med_vars,1),out_random$theta_median)
		assign(tail(theta_med_homo_vars,1),out_homo$theta_median)

		assign(tail(beta_knots_vars,1),out_random$beta_knots_avg)
		assign(tail(beta_knots_homo_vars,1),out_homo$beta_knots_avg)
		
		assign(tail(beta_ise_vars,1),out_random$beta_ise)
		assign(tail(beta_ise_homo_vars,1),out_homo$beta_ise)

		theta_tab <- createThetaView(get(tail(theta_vars,1)),get(tail(theta_homo_vars,1)))
		savepath <- paste0(path,'model',model,'/')
		write.csv(theta_tab, paste0(savepath,'theta_size',size,'.csv'),row.names=FALSE)
		theta_temp_res <- createThetaTable(get(tail(theta_med_vars,1)),get(tail(theta_med_homo_vars,1)),
											type=size)

		theta_res <- rbind(theta_res, theta_temp_res)

		beta_temp_res <- createBetaTable(beta_knots_avg=get(tail(beta_knots_vars,1)),beta_knots_homo_avg=get(tail(beta_knots_homo_vars,1)),
										type=paste0('Pool Size = ',size),which_beta=which_beta,beta_list=beta_list)
		beta_res 	<- rbind(beta_res,beta_temp_res)

		beta_temp_ise_res <- createBetaISE(beta_ise=get(tail(beta_ise_vars,1)),beta_ise_homo=get(tail(beta_ise_homo_vars,1)), type=size)
		beta_ise_res <- rbind(beta_ise_res,beta_temp_ise_res)

	}

	# confidence interval/prediction interval band
	mybands <- function(x, y, upper, lower, col, subscripts, ..., font, fontface){
		upper <- upper[subscripts]
		lower <- lower[subscripts]
		panel.polygon(c(x, rev(x)), c(upper, rev(lower)), col = band_col, alpha=band_opacity, border = FALSE,...)
	}

	# panel setting
	mypanel <- function(x, y, ...){
		panel.superpose(x, y, panel.groups=mybands, type='l',...)
		panel.xyplot(x, y, type='l', cex=0.6, lty=c(1,2,1,1),lwd=c(3,3,0.15,0.15), col=c(line_col,'black','black'),...)
	}
	# strip setting
	mystrip <- function(which.given, ..., factor.levels) {
		levs <- if (which.given == 1) factor.levels
		else   c('Random','Homogeneous')#c('IT','MPT','DT','AT')
		strip.default(which.given, ..., factor.levels = levs)
	}
	# configure xyplot setting
	line_col <- plot_config$line_col
	band_col <- plot_config$band_col
	band_opacity <- plot_config$band_opacity

	bwpanel <- function(...){
		panel.bwplot(...)
		#panel.points(x=c(1,0.5), pch="x", cex=2,col='yellow')
		panel.points(1,0.5, pch=16, cex=1.4,col='green',alpha=0.5)
		panel.points(2,0.5, pch=16, cex=1.4,col='green',alpha=0.5)
		panel.points(3,0.5, pch=16, cex=1.4,col='green',alpha=0.5)
		panel.points(4,0.5, pch=16, cex=1.4,col='green',alpha=0.5)
		panel.points(5,0.5, pch=16, cex=1.4,col='green',alpha=0.5)
	}

	# generate boxplot
	# box_opacity <- 0.4
	# sigma_res <- theta_res[theta_res['varname']=='sigma',]
	# boxfigure <- bwplot(theta~type,data=sigma_res,
	# 	col="salmon",pch=0,panel=bwpanel,cex=1.4,
	# 	xlab='Pool Size', ylab=bquote(sigma),
	# 	scales=list(y=list(rot=0)),
	# 	key=list(space="top",
	# 			column=2,
	# 			text=list(label=c('Random','Homogeneous'),cex=0.8,col="darkgrey"),
	# 			points=list(col=c("salmon","dodgerblue"),cex=1.4,pch=c(0,5),alpha=box_opacity),
	# 			rectangles=list(col=c("salmon","dodgerblue"),rectangles=c("salmon","dodgerblue"),alpha=box_opacity,border=c(FALSE,TRUE),size=3)),
	# 			par.settings=list(box.rectangle=list(col="salmon",fill="salmon",alpha=box_opacity+0.2,border=FALSE),
	# 								box.umbrella=list(col="salmon",alpha=.3),
	# 								plot.symbol=list(col="salmon",alpha=box_opacity+0.2)))

	# boxfigure <- boxfigure + as.layer(bwplot(theta~type,data=sigma_res,
	# 								col="dodgerblue",pch=5,panel=bwpanel,cex=1.4,
	# 								par.settings=list(box.rectangle=list(col="dodgerblue",fill="dodgerblue",alpha=box_opacity),
	# 													box.umbrella=list(col="dodgerblue",alpha=0.3),
	# 													plot.symbol=list(col="dodgerblue",alpha=box_opacity))))
	# if(plot_sigma){
	# 	print(boxfigure)
	# }
	
	bwstrip <- function(which.given, ..., factor.levels) {
		levs <- if (which.given == 1) factor.levels
		else parse(text=paste0('beta[',1:nbeta-1,']'))
		strip.default(which.given, ..., factor.levels = levs)
	}

	ise_figure <- bwplot(
	  ise~type | label*varname, 
	  data=beta_ise_res,
	  layout=c(length(sizes),2),
	  xlab='Pool Size',
	  ylab=bquote(.('ISE of') ~ beta),
	  box.ratio=1,
	  #panel = function(x, box.ratio=box.ratio,box.width = box.ratio / (1 + box.ratio),...) {
	    #level <- dimnames(trellis.last.object())[["varname"]][packet.number()]
	    #panel.abline(h = abline_df$true[abline_df$varname==level],lty=2,col='green',lwd=2);
	    #panel.bwplot(x, box.ratio=box.ratio,...);
	  #},
	  #par.settings=list(plot.symbol=list(alpha=0.4,cex=0.5),fontsize = list(text = 8, points = 4)),
	  #key=list(space='top',column=1,text=list('True',cex=0.8),lines=list(lty=2,col='green',lwd=2)),
	  par.settings =list(strip.background=list(col=c('#E5E7E9','#9CC3D5FF'))),
	  strip=bwstrip
	)
	
	sigma_figure <- bwplot(
	  theta~type | label, 
	  data=theta_res,
	  #layout=,
	  xlab='Pool Size',
	  ylab=bquote(sigma),
	  box.ratio=1,
	  panel = function(x, box.ratio=box.ratio,box.width = box.ratio / (1 + box.ratio),...) {
	    #level <- dimnames(trellis.last.object())[["varname"]][packet.number()]
	    #panel.abline(h = abline_df$true[abline_df$varname==level],lty=2,col='green',lwd=2);
	    panel.abline(h=sigma_true,lty=2,col='green',lwd=2)
	    panel.bwplot(x, box.ratio=box.ratio,...);
	  },
	  #par.settings=list(plot.symbol=list(alpha=0.4,cex=0.5),fontsize = list(text = 8, points = 4)),
	  #key=list(space='top',column=1,text=list('True',cex=0.8),lines=list(lty=2,col='green',lwd=2)),
	  #strip=bwstrip
	  par.settings =list(strip.background=list(col=c('#E5E7E9','#9CC3D5FF'))),
	)
	
	# ci_panel <- function(x, y, upper, lower, fill, col, subscripts, ..., font, fontface){
	#   upper <- upper[subscripts]
	#   lower <- lower[subscripts]
	#   panel.polygon(c(x, rev(x)), c(upper, rev(lower)),
	#                 col =c('darkgray'),alpha=0.4,border = FALSE,...)
	# }

	# xyplot(theta~type|label,groups=label,
	#        data=theta_res,
	#        xlab='Pool Size',
	#        ylab=bquote(theta),
	#        lower = theta_res$lower,
	#        upper = theta_res$upper,
	#        #layout=c(5,2),
	#        ylim=c(0.4,0.6),
	#        panel = function(x, y, ...){
	#          #level <- dimnames(trellis.last.object())[["varname"]][panel.number()]
	#          #panel.abline(h = abline_df$true[abline_df$varname==level],lty=2,col='green',lwd=2)
	#          panel.superpose(x, y, 
	#                          panel.groups = ci_panel, 
	#                          type='l',...)
	#          panel.xyplot(x, y, type='b',cex=0.7, lty=1,...)
	#        },
	#        par.settings =list(strip.background=list(col=c('#E5E7E9','#9CC3D5FF'))),
	#        #key=list(space='top',column=1,text=list('True',cex=0.8),lines=list(lty=2,col='green',lwd=1.2)),
	# )



	# generate figure for beta
	beta_figure <- xyplot(true+median+lower+upper~t |type*label,
			 		data=beta_res,
			 		#data = replace(beta_res, "label", structure(beta_res$label, levels = )),
			 		auto.key = list(space='top',
			 					border=FALSE, 
			 					points=F, 
			 					lines=T, 
			 					#rectangles=T,
			 					columns = 3, 
			 					text=c('True','Median Estimate','95% Credible Interval Region')),
			   		lower=beta_res$lower,upper=beta_res$upper,
			   		panel = mypanel,
			   		par.settings =list(superpose.line = list(col=c(line_col,'#E5E7E9'), lwd = c(3,3,15), lty=c(1,2,1)),
			   							#superpose.polygon = list(col=c(NA,NA,'#D5D8DC'),border=c(0,0,0)),
			   							#strip.background=list(col=c('#9CC3D5FF','#3498DB'))),
			   							#strip.background=list(col=c('#3498DB','#9CC3D5FF'))),
			   							#strip.background=list(col=c('#F1F4FFFF','#A2A2A1FF'))),
			   							strip.background=list(col=c('#E5E7E9','#9CC3D5FF'))),
			 		xlab='age',
			 		ylab=bquote(beta[.(which_beta)]),
			 		#scales=list(y=list(relation="free")),
			 		#ylim=c(-0.5,1),
			 		strip=mystrip,
			 		layout=c(length(sizes),2))
	
	if(plot_beta){
		trellis.device(png,width=540,height=527,file=paste0(path,'model',model,'/beta_res.png'))
		print(beta_figure)
		dev.off()
		print(beta_figure)
		outpath <- paste0(path,'model',model,'/beta',which_beta,'_res.csv') 
		write.csv(beta_res, outpath)
	}
	if(plot_sigma){
		trellis.device(png,width=540,height=527,file=paste0(path,'model',model,'/sigma_res.png'))
		print(sigma_figure)
		dev.off()
		print(sigma_figure)
		outpath <- paste0(path,'model',model,'/sigma_res.csv') 
		write.csv(theta_res, outpath)
	}
	if(plot_ise){
		trellis.device(png,width=540,height=527,file=paste0(path,'model',model,'/ise_res.png'))
		print(ise_figure)
		dev.off()
		print(ise_figure)
		outpath <- paste0(path,'model',model,'/ise_res.csv') 
		write.csv(beta_ise_res, outpath)
	}
	
}