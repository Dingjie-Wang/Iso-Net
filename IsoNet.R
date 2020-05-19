# Aim: Construction of co-expression networks based on exon-level data -----------
# Author: Dingjie Wang, Xiufen Zou and Kin Fai Au-----------------
# Date: May 13, 2020--------------------
# Contact: Dingjie.Wang@osumc.edu-------


# libraries to be loaded----------------
# Please install packages: CompQuadForm and survey
library(CompQuadForm);
library(survey); 

#--------useful function----------------
abstract_exp = function(name1,name2,exp){
	N = dim(exp)[2];
	X1 = exp[which(exp[,1]==name1),3:N];
	X2 = exp[which(exp[,1]==name2),3:N];
	return(list(t(X1),t(X2)));
}

#------------MINet function-------------
MINet = function(gene1,gene2){
	N = dim(gene1)[1];
	p = dim(gene1)[2];
	q = dim(gene2)[2];
	Cor = cor(gene1,gene2);
	Cor1 = cor(gene1);
	Cor2 = cor(gene2);
	lambda1 = eigen(Cor1,only.values=T)$values;
	lambda2 = eigen(Cor2,only.values=T)$values;
	Lamda = rep(lambda1,each=q)*lambda2;
	Index = which(abs(Lamda)<1e-08);
	if(length(Index)>0){
		Lamda[which(abs(Lamda)<1e-08)] = 0;
	}
	
      
      MICor = -0.5*log(1-Cor^2);
	MIFnorm = (N-2.5)*sum(MICor^2);

	results = davies(q=MIFnorm,lambda=Lamda,lim=1000000,acc=1e-9);
	if(results$ifault==0 | results$ifault==2){
		if(results$Qq>=0){
			p.value_MIFnorm = results$Qq;
		}else{
			p.value_MIFnorm = pchisqsum(x=MIFnorm, df=rep(1,length(Lamda)), a=Lamda, lower.tail = FALSE,method='saddlepoint');
		}
	}else{
		p.value_MIFnorm = pchisqsum(x=MIFnorm, df=rep(1,length(Lamda)), a=Lamda, lower.tail = FALSE,method='saddlepoint');
	}
	return(data.frame(MIFnorm=MIFnorm,p.value_MIFnorm=p.value_MIFnorm));
}

#--------------RVNet function---------------------
RVNet = function(gene1,gene2){
     RVCor = sum(diag(gene1%*%t(gene1)%*%gene2%*%t(gene2)))/sqrt(sum(diag(gene1%*%t(gene1)%*%gene1%*%t(gene1)))*sum(diag(gene2%*%t(gene2)%*%gene2%*%t(gene2))));
     return(data.frame(RVcorrelation=RVCor));
}


#---------------------------------------

#Outputfile: this is the name of the output file with ".txt" extension.
#Inputfile: this is the name of the input file with ".txt" extension.

PerformIsoNet = function(Inputfile,Outputfile,sample_size){

	Exp	= read.table(Inputfile,header=T,sep="\t");
	Exp[,c(1,2)] = as.character(as.matrix(Exp[,c(1,2)]));
	N = dim(Exp)[2]-2;
	Var = apply(Exp[,3:dim(Exp)[2]],1,var);
	Exp = Exp[which(Var>0),]
	gene  = unique(Exp[,1]);
	result = data.frame();
	
	for(i in 1:(length(gene)-1)){
		for(j in (i+1):length(gene)){
			X1 = abstract_exp(gene[i],gene[j],Exp)[[1]];
			X2 = abstract_exp(gene[i],gene[j],Exp)[[2]];
                  if (sample_size<=20){
                       result = rbind(result,data.frame(gene1=gene[i],gene2=gene[j],RVNet(X1,X2))); 
                  } else{
                       result = rbind(result,data.frame(gene1=gene[i],gene2=gene[j],MINet(X1,X2)));
                  }		
		}
	}

	write.table(result,file=Outputfile,sep="\t",quote=F,row.names=F,col.names=T);
}
