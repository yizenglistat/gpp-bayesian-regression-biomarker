bayesPlot <- function(path='./output/', plot_figure=FALSE, which_beta=1, models=models, sizes=sizes,
	true_params=list(sigma_true=sigma_true,alpha_true=alpha_true,beta_list=beta_list),
	plot_config = list(line_col=c("#28B463",'#A93226'), band_col='gray', band_opacity=0.3)){
	
	sigma_true 	<- true_params$sigma_true
	alpha_true 	<- true_params$alpha_true
	beta_list 	<- data_lst$beta_list
	nbeta 		<- length(beta_list)
	nalpha 		<- length(alpha_true)

	if(substr(path,nchar(path),nchar(path))!='/') path <- paste0(path,'/')

	# require pacakges
	require(lattice)
	require(data.table)
	graphics.off()

	# read/load data automatically for beta_knots_summ and theta_summ
	readFiles <- function(path,model,size,homo=FALSE){
		path <- paste0(path,'model',model,'/size',size,'/')
		if(homo) path <- paste0(path,'homo/')
		files <- list.files(path)
		
		beta_knots_names 	<- c()
		beta_knots_total 	<- 0
		beta_knots_count 	<- 0

		theta_names 		<- c()
		theta_total 		<- 0
		theta_count 		<- 0

		for(filename in files){
			if(grepl('beta_knots_summ',filename)){
				name <- unlist(strsplit(filename,'.csv'))
				filepath <- paste0(path, filename)
				assign(name, as.matrix(fread(filepath)))
				beta_knots_total <- beta_knots_total + get(name)
				beta_knots_names <- c(beta_knots_names, name)
				beta_knots_count <- beta_knots_count + 1
			}
			if(grepl('theta_summ',filename)){
				name <- unlist(strsplit(filename,'.csv'))
				filepath <- paste0(path, filename)
				assign(name, as.matrix(fread(filepath)))
				theta_total <- theta_total + get(name)
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
		return(list(beta_knots_avg=beta_knots_avg,theta_avg=theta_avg))
	}

	createThetaTable <- function(theta_avg,theta_homo_avg){

		theta_true 			<- c(sigma_true,	
								alpha_true,
								rep(NA,2*nbeta))

		theta_res <- data.frame(theta_homo_avg,true=theta_true, theta_avg,
								row.names=c('sigma',
											paste0('alpha',1:nalpha),
											paste0('tau',1:nbeta),
											paste0('phi',1:nbeta)))
		return(theta_res)
	}

	createBetaTable <- function(beta_knots_avg,beta_knots_homo_avg,type,which_beta,beta_list=beta_list){
		# beta varying function visualization
		t_knots 			<- beta_knots_avg[,'t_knots']
		nknots 				<- length(t_knots)
		beta_true 			<- betaTrue(t_knots, beta_list, center_intc=TRUE)
		beta_index 			<- seq(2,nbeta*5,nbeta) + which_beta
		beta_res <- data.frame(t_knots=t_knots,
							beta_true=beta_true[,which_beta+1],
							beta_knots_avg[,beta_index],
							type=type,
							label=rep('Heterogeneous',each=nknots))
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

	theta_vars <- c()
	theta_homo_vars <- c()
	
	beta_knots_vars <- c()
	beta_knots_homo_vars <- c()

	beta_res <- c()

	for(model in models){

		for(size in sizes){
			
			theta_vars <- c(theta_vars,paste0('theta_size',size))
			theta_homo_vars <- c(theta_homo_vars,paste0('theta_size',size))
			
			beta_knots_vars <- c(beta_knots_vars,paste0('beta_knots_size',size))
			beta_knots_homo_vars <- c(beta_knots_homo_vars,paste0('beta_knots_size',size))
			
			assign(tail(theta_vars,1),readFiles(path,model=model,size=size,homo=FALSE)$theta_avg)
			assign(tail(theta_homo_vars,1),readFiles(path,model=model,size=size,homo=TRUE)$theta_avg)
			
			assign(tail(beta_knots_vars,1),readFiles(path,model=model,size=size,homo=FALSE)$beta_knots_avg)
			assign(tail(beta_knots_homo_vars,1),readFiles(path,model=model,size=size,homo=TRUE)$beta_knots_avg)
			
			theta_res <- createThetaTable(get(tail(theta_vars,1)),get(tail(theta_homo_vars,1)))
			savepath <- paste0(path,'model',model)
			write.csv(theta_res, paste0(savepath,'theta_size',size,'.csv'))
			
			beta_temp_res <- createBetaTable(beta_knots_avg=get(tail(beta_knots_vars,1)),beta_knots_homo_avg=get(tail(beta_knots_homo_vars,1)),
											type=paste0('Pool Size = ',size),which_beta=which_beta,beta_list=beta_list)
			beta_res 	<- rbind(beta_res,beta_temp_res)

		}
	
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
		else   c('Heterogeneous','Homogeneous')#c('IT','MPT','DT','AT')
		strip.default(which.given, ..., factor.levels = levs)
	}
	# configure xyplot setting
	line_col <- plot_config$line_col
	band_col <- plot_config$band_col
	band_opacity <- plot_config$band_opacity

	# generate figure 
	figure <- xyplot(true+median+lower+upper~t |type*label,
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
			 		layout=c(5,2))
	
	if(plot_figure){
		trellis.device(png,width=540,height=527,file=paste0(path,'model',model,'/beta_res.png'))
		print(figure)
		dev.off()
		print(figure)
	}
	outpath <- paste0(path,'model',model,'/beta',which_beta,'_res.csv') 
	write.csv(beta_res, outpath)
}