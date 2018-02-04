# Model-based imputaion  - last save 23/06/16
# Ref:  "A statistical framework for protein quantitation in bottom-up MS-based
#        proteomics. Karpievitch Y, Stanley J, Taverner T, Huang J, Adkins JN,
#        Ansong C, Heffron F, Metz TO, Qian WJ, Yoon H, Smith RD, Dabney AR.
#        Bioinformatics 2009
#
# Written by Yuliya Karpievitch, Tom Taverner, and Shelley Herbrich
# for TAMU, PPNNL and community

#   source('MBfilter.r')
#   source('MBImpute.r')
#   fnameLum = 'D:/yuliya/pnnl_cod/my_functions/data/proteins2_lumican_anthithr.txt'
#   datasetLum = read.table(fnameLum, header=TRUE, sep="\t")
#   trLum = c(1,1,1,1,1,1,1,1,1,1, 2,2,2,2,2,2,2,2,2,2)
#   resMBfilter = MBfilter(datasetLum[,-(1:2)], trLum, datasetLum[,1:2], pr_ppos=2, my.pi=.05)
#
#   resMBimpute = MBimpute(resMBfilter$y_filtered,trLum,resMBfilter$ft_prot.info,pr_ppos=2)



#' Model-Based Imputation of missing values
#'
#' Impute missing values based on information from multiple peptides within a protein
#' Expects the data to be filtered to contain at least one observation per treatment
#' group.
#'
#' @param mm number of peptides x number of samples matrix of intensities
#' @param treatment vector indicating the treatment group of each sample eg
#'        as.factor(c('CG','CG','CG', 'mCG','mCG','mCG')) or c(1,1,1,1,2,2,2,2)
#' @param pr_ppos column index for protein ID in prot.info
#' @param my.pi PI value, estimate of the proportion of peptides missign completely at
#'              random, as compared to censored at lower abundance levels default values
#'              of 0.05 is usually reasoanble for missing completely at random values
#'              in proteomics data
#' @param compute_pi TRUE/FALSE (default=FALSE) estimate Pi is set to TRUE, otherwise
#'              use the provided value. We consider Pi=0.05 a reasonable estimate for
#'              onservations missing completely at random in proteomics experiments.
#'              Thus values is set to NOT estimate Pi by default. Note: spline smoothing
#'              can sometimes produce values of Pi outside the range of possible values.
#'
#' @return A structure with multiple components
#' \describe{
#'   \item{y_imputed}{number of peptides x m matrix of peptides with no missing data}
#'   \item{imp_prot.info}{imputed protein info, 2+ columns: peptide ID, protein IDs, etc
#'                  Dimentions should be the same as passed in}
#'}
#' @examples
#' data(mm_peptides)
#' head(mm_peptides)
#' intsCols = 8:13 # different from parameter names as R uses outer name spaces if variable is undefined
#' metaCols = 1:7 # reusing this variable
#' m_logInts = make_intencities(mm_peptides, intsCols)  # will reuse the name
#' m_prot.info = make_meta(mm_peptides, metaCols)
#' m_logInts = convert_log2(m_logInts)
#' grps = as.factor(c('CG','CG','CG', 'mCG','mCG','mCG')) # 3 samples for CG and 3 for mCG
#' mm_m_ints_eig1 = eig_norm1(m=m_logInts,treatment=grps,prot.info=m_prot.info)
#' mm_m_ints_eig1$h.c # check the number of bias trends detected
#' mm_m_ints_norm = eig_norm2(rv=mm_m_ints_eig1)
#' mm_prot.info = mm_m_ints_norm$normalized[,1:7]
#' mm_norm_m =  mm_m_ints_norm$normalized[,8:13]
#' imp_mm = MBimpute(mm_norm_m, grps, prot.info=mm_prot.info, pr_ppos=2, my.pi=0.05,
#'                   compute_pi=FALSE, sseed=131)
#' @export
MBimpute = function(mm, treatment, prot.info, pr_ppos=2, my.pi=0.05, compute_pi=FALSE, sseed=123457){

  set.seed(sseed) # must be rst for reproducibility
  # calculate PI or use one passed in as parameter:
  if (compute_pi){
    print('Estimating Pi')
  	my.pi = eigen_pi(mm, toplot=T)
  }

  # (From EigenMS) check if treatment is a 'factor' vs data.frame', i.e. single vs multiple factors
  if(class(treatment) == "factor") { # TRUE if one factor
    n.treatment = 1 # number of factors in the model
    n.u.treatment = length(unique(treatment))[1] # number of levels within the only factor
  } else { # data.frame
    n.treatment = dim(treatment)[2]
    n.u.treatment = dim(unique(treatment))[1] # all possible treatment combinations
  }

  # Match to protein
  all.proteins = unique(prot.info[,pr_ppos])
  y_imputed = NULL
  imp_prot.info = NULL
  prot.info.tmp = NULL
  cat("Imputing...\n")

  for (kk in 1:length(all.proteins)){
    #if(kk == 191) browser()
    prot = all.proteins[kk]
    pmid.matches = prot.info[prot.info[,pr_ppos]==prot,1]

    curr_prot.info = prot.info[prot.info[,pr_ppos]==prot,]
    idx.prot = which(prot.info[,1] %in% pmid.matches)  # do not obligate
                              # to be rownames, as those have to be Unique
    y_raw = mm[idx.prot,,drop=F]
    cat(paste("Protein: ", prot, ": ", dim(y_raw)[1], " peptides (", kk, "/", length(all.proteins), ")", sep="" ))
    y_info = prot.info[idx.prot,,drop=F]

    if (nrow(y_raw) == 0) next  # yuliya: this should not happen here (NO observations)

    # peptides and proteins of poor quality are removed prior to analysis
    # "poor quality" was defined (Karpievitch et al. 2009) as having little "information"
    # abt grp differences compared to other peptides estimate data parameters.
    # To speed up execution we allow EigenMS to eliminat all peptides with fewer than
    # one observation per treatment group and impte the rest.
    n.peptide = nrow(y_raw)
    yy = as.vector(t(y_raw))
    nn = length(yy)
    peptide =rep(1:n.peptide, each=dim(data.frame(treatment))[1])

    # filter out min.missing
    n.present = array(NA, c(n.peptide, n.u.treatment))
    colnames(n.present) = unique(treatment)

    for(jj in 1:n.u.treatment) {
       n.present[,jj] = rowSums(!is.na(y_raw[,treatment==unique(treatment)[jj],drop=F]))
    }
    # remove peptides with completely missing group(s)
    present.min = apply(n.present, 1, min)
    ii = present.min > 0
    y_raw = y_raw[ii,,drop=F] # reassign Y_raw to a submatrix of 1+ observations in each group

    # keep track of pepIDs and prIDs here...
    if (nrow(y_raw) == 0) next
    c.guess = min(yy, na.rm=T)
	  peptide = rep(1:n.peptide, each=dim(data.frame(treatment))[1])

    # make column names for the n.present matrix
    tmp = unique(treatment)
    nrow_tmp = dim(tmp)[1] # yuliya: this does nto work if a factor!!
    # R does not remove the variable col_names1 from name space outside of the if/else..
    if(class(tmp)== 'factor') {
     col_names1 = tmp
    } else {
      col_names1 = vector('character', nrow_tmp)
      bob = data.frame(lapply(tmp, as.character), stringsAsFactors=FALSE)
      for(ii in 1:nrow_tmp) {
       col_names1[ii] = paste(bob[ii,1], bob[ii,2], sep='_')
      }
    }

    # calculate pooled variance for each protein
    grp = array(NA, c(1, n.u.treatment))
    for (jj in 1:n.u.treatment){
      grp[jj] = sum(n.present[, jj])
    }
    mpos = which.max(grp)
    pep_var = 0 # yuliya: was not declared in impute only
    go = protein_var(y_raw, treatment) # yuliya: function, see below
    overall_var = go$overall_var

    num = 0
    den = 0
    # calculate pooled variance for each peptide
    # if only 1 onservation in a dx group assign the overall variance
    for(ii in 1:n.peptide)
    {
      # which treatment has more obsertations? - use one of the 2 groups
      pep = yy[peptide==ii]
      most_obs = which(max(n.present[ii,]) == n.present[ii,])
      lala = length(most_obs)
      if(lala > 1) { # more than 1 gorup has max # of observations, pick one at random
        most_obs = most_obs[sample(lala,1)]
      }
      tmp = data.frame(most_obs)
      treat_to_use = rownames(tmp)
      y.i = pep[treatment == treat_to_use]
      p2 = var(y.i, na.rm=TRUE)
      if (is.na(p2)) { p2 = 0 }
      present = length(y.i)
      num = num + (p2 * (present-1))
      den = den + (present-1)
    }
    pep_var = num/den

    if(is.na(pep_var)) {
      pep_var = overall_var
      cat(idx.prot)
      cat('\t')
    } # only occurs if we have 1 pep per group
    if (pep_var == 0) { pep_var = overall_var }
    if (pep_var == 0) { pep_var = .0001 } # if overall var is 0...
    peptides.missing = rowSums(is.na(y_raw))

    f.treatment = factor(rep(treatment, n.peptide))
    f.peptide = factor(peptide)

    # estimate rough model parameters
    # create model matrix for each protein
    # remove any missing values from consideration
    iii = (1:nn)[is.na(yy)] # positions of the NA's
    if (n.peptide != 1){
      X  = model.matrix(~f.peptide + f.treatment, contrasts = list(f.treatment="contr.sum", f.peptide="contr.sum") )
    } else {
      X = model.matrix(~f.treatment, contrasts=list(f.treatment="contr.sum"))
    }
    if(length(iii) > 0){
      y.c = yy[-iii] # remove NA value
      X.c = X[-iii,] # remove rows for NA values from the model matrix
    } else {
      y.c = yy
      X.c = X
    }

    # calculate initial beta values and residuals
    beta = drop(solve(t(X.c) %*% X.c) %*% t(X.c) %*% y.c)

    # compute initial delta's
    peptides.missing[peptides.missing==0] = 0.9
    delta.y = as.numeric(1/sqrt(pep_var*peptides.missing))
    dd = delta.y[as.numeric(peptide)]

    c_hat = apply(y_raw,1,min,na.rm=TRUE)
    c_h = c_hat[as.numeric(peptide)]

    if(n.peptide==1){
      y.predict = model.matrix(~f.treatment, contrasts=list(f.treatment="contr.sum"))%*% beta
    } else {
      y.predict = model.matrix(~f.peptide + f.treatment,contrasts = list(f.treatment="contr.sum", f.peptide="contr.sum"))%*% beta
    }

    zeta = dd*(c_h - y.predict)
    PHI = pnorm(zeta, 0, 1)
    # prob.cen = pnorm(zeta, 0, 1)/(my.pi + (1-my.pi) * pnorm(zeta, 0, 1))
    prob.cen = PHI / ( (my.pi + (1-my.pi) * PHI) )

    choose.cen = runif(nn) < prob.cen
    set.cen = is.na(yy) & choose.cen
    set.mar = is.na(yy) &! choose.cen
    kappa = my.pi + (1 - my.pi) * dnorm(zeta,0, 1)

    # Imputation: Replace missing values with random numbers drawn from the estimated
    # likelihood model
    sigma = 1/dd
    y.impute = t(y_raw)
    if(sum(set.cen) > 0) # censored
      mus = y.predict[set.cen]
      ss = sigma[set.cen]
      cutoff = c_h[set.cen] # rep(c.guess, nn)[set.cen]
      y.impute[set.cen] = rnorm.trunc(sum(set.cen), mus, ss, hi=cutoff)

    if(sum(set.mar) > 0) # randomly missing
      y.impute[set.mar] = rnorm(nn, y.predict, sigma)[set.mar]

    y.impute.return = t(y.impute)
    imp_prot.info = rbind(imp_prot.info,curr_prot.info)
    y_imputed = rbind(y_imputed, y.impute.return)
  } # end for each protein

  colnames(y_imputed) = colnames(mm)
  cat("Done imputing.\n")
  return(list(y_imputed=y_imputed,imp_prot.info=imp_prot.info,pi=as.matrix(my.pi)))
}
############# end imputation ###################




############## function pi #############
#' Compute PI - proportion of observations missing completely at random
#'
#' @param m matrix of abundances, numsmaples x numpeptides
#' @param toplot TRUE/FALSE plot mean vs protportion missing curve and PI
#' @return pi estimate of the proportion of observations missing completely at random
#'
#' Contributed by Shelley Herbrich & Tom Taverner for Karpievitch et al. 2009
#' @export
eigen_pi = function(m, toplot=T)
{
  # (1) compute 1) ave of the present values from each petide
  #             2) number of missing and present values for each peptide

  #remove completely missing rows
  m = m[rowSums(m, na.rm=T)!=0,]

  pepmean = apply(m, 1, mean, na.rm=T)
  propmiss = rowSums(is.na(m))/ncol(m)

  smooth_span = (0.4)
  fit = lowess(pepmean, propmiss, f=smooth_span)
  PI = fit$y[fit$x==max(pepmean)]

  count = 1
  while (PI<=0){
    smooth_span = smooth_span-.1
    fit = lowess(pepmean, propmiss, f=smooth_span)
    PI = fit$y[fit$x==max(pepmean)]
    count = count + 1
    if (count > 500) break
  }

  if (toplot){
  st = paste("PI: ", PI)
  plot(pepmean, propmiss, xlab="x", ylab="y", cex=0.5) #plot data point
  lines(fit)
  title("Lowess Regression", sub = st,
      cex.main = 2,   font.main= 3, col.main= "purple",
      cex.sub = 1, font.sub = 3, col.sub = "red")
  }
  return (pi=PI)
}


######################################################
protein_var = function(Y_raw, treatment){
# estimates coefficients for all peptides from a single protein
# a portion of what get_coeffs does , without information computations
# that is nto needed in imputation
#
# Input:
#   Y_raw: m peptides by n samples arrays matrix of expression data
#          from a given protein
#   treatment: treatment groups
#
# Output:
#   overall_var:  need to check...

  n.peptide = nrow(Y_raw)
  y = as.vector(t(Y_raw))
  n = length(y)
  n.treatment = length(treatment)
  n.u.treatment = length(unique(treatment))
  peptide =rep(rep(1:n.peptide, each=n.treatment))

  n.present = array(NA, c(n.peptide, n.u.treatment))
  colnames(n.present) = unique(treatment)
  for(i in 1:n.peptide) for(j in 1:n.u.treatment) n.present[i,j] = sum(!is.na(y [peptide==i & treatment==unique(treatment)[j]]))

  peptides.missing = rowSums(is.na(Y_raw))

  f.treatment = factor(rep(treatment, n.peptide))
  f.peptide = factor(peptide)

  # estimate rough model parameters
  # create model matrix for each protein and
  # remove any peptides with missing values
  ii = (1:n)[is.na(y)]
  if (n.peptide != 1){
    X  = model.matrix(~f.peptide + f.treatment, contrasts = list(f.treatment="contr.sum", f.peptide="contr.sum") )
  } else {
    X = model.matrix(~f.treatment, contrasts=list(f.treatment="contr.sum"))
  }
  if(length(ii) > 0){
    y.c = y[-ii]
    X.c = X[-ii,]
  } else {
    y.c = y
    X.c = X
  }

  # calculate initial beta values and residuals
  beta = drop(solve(t(X.c) %*% X.c) %*% t(X.c) %*% y.c)
  Y_hat = X.c %*% beta
  Y_temp = Y_raw
  # Y_temp = as.numeric(t(Y_temp))
  Y_temp = as.vector(t(Y_temp)) # yuliya, same as in filter now
  Y_temp[!is.na(Y_temp)] = Y_hat
  Y_temp = matrix(Y_temp, nrow = n.peptide, byrow = T)
  Y_hat = Y_temp
  ee = Y_raw - Y_hat

  effects = X.c %*% beta
  resid = y.c - effects
  overall_var = var(resid)
  return(list(overall_var=det(overall_var)))
}
################################################
my.Psi = function(x, my.pi){
# calculates Psi
exp(log(1-my.pi)  + dnorm(x, 0, 1, log=T) - log(my.pi + (1 - my.pi) * pnorm(x, 0, 1) ))
}
# end my.Psi

my.Psi.dash = function(x, my.pi){
# calculates the derivative of Psi
-my.Psi(x, my.pi) * (x + my.Psi(x, my.pi))
}
# end my.Psi.dash

phi = function(x){dnorm(x)}

rnorm.trunc = function (n, mu, sigma, lo=-Inf, hi=Inf){
# Calculates truncated noraml
  p.lo = pnorm (lo, mu, sigma)
  p.hi = pnorm (hi, mu, sigma)
  u = runif (n, p.lo, p.hi)
  return (qnorm (u, mu, sigma))
}
# End rnorm.trunc



################################################
# dialog for DanteR
################################################
MBimpute.dialog = list( title='Model-based Imputation',
  m.dataframeItem=NULL, label='Data to impute',
  signal = list("default", "get.dataset.factors", "treatment", user.data=list(include.primary=FALSE)),
  signal = list("default", "get.dataset.row.metadata.fields", "protein_group", user.data=list(include.primary=FALSE)),
  treatment.choiceItem=NULL, label='Factors',
  protein_group.choiceItem = NULL, label = "Field Containing Proteins",
  compute_pi.trueFalseItem=FALSE, label='Manually Estimate PI',
   tooltip = "If this is checked, the user can set the estimated random missingness probability",
  signal = c("default", "toggle.sensitive", "my.pi"),
    my.pi.numericItem="0.05", label='PI', indent=10
)




get_pi = function (choice){
		if (choice=="Yes"){
		my.pi = my.pi
	}else {
		my.pi = run.dialog(enter_pi)$retval
	}
  return(my.pi)
}

enter_pi = function(new_val){
	my.pi = new_val
	return(my.pi)
}
enter_pi.dialog=list(title='New PI Value', new_val.numericItem=NULL, label='Enter new PI value:')


