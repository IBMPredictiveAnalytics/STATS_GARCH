#/***********************************************************************
# * Licensed Materials - Property of IBM 
# *
# * IBM SPSS Products: Statistics Common
# *
# * (C) Copyright IBM Corp. 2014
# *
# * US Government Users Restricted Rights - Use, duplication or disclosure
# * restricted by GSA ADP Schedule Contract with IBM Corp. 
# ************************************************************************/

# author__ = "SPSS, JKP"
# version__ = "1.0.0"

# History
# 25-Jun-2014 original version


helptext='Estimate a GARCH model

Note: For mathematical details of the models and methods
used here, see
Introduction to the rugarch package.pdf.
This will be found under the R installation directory
in library/rugarch/doc.

STATS GARCH 
VARIABLE=varname MODELSOURCE=ESTIMATE* or FILE or WORKSPACE
MODELFILE="file specification"
START=number number FREQUENCY=frequency if start has two numbers

/FIT 
HOLDOUT=number of holdout cases
SOLVER=nlminb or solnp or lbfgs or gosolnp or nloptr or hybrid
NRESTARTS=number of restarts PARALLEL=YES or NO*
PACKAGE=snowfall or multicore
NCORES=number of cores to use
NSIM=number
SEED=random number seed.

/SPECIFICATION
VARMODEL=SGARCH or FGARCH or EGARCH or GJRGARCH or APARCH
  or IGARCH or CSGARCH
VARGARCHORDER=arorder maorder
SUBMODEL=GARCH or TGARCH or AVGARCH or NGARCH or NAGARCH or APARCH
or GJRGARCH or ALLGARCH
VARREGRESSORS=list of variables
VARTARGETING = YES or NO*
MEANORDER= arorder maorder
INCLUDEMEAN=YES* or NO
MEANARCHPOW=STDDEV* or VAR
MEANARCHM=YES or NO  MEANARFIMA=YES or NO  MEANARCHEX=YES or NO
MEANREGRESSORS=list of variables
DISTRIBUTION = NORM* or SNORM or STDT or SSTD or GED or SGED or NORMINVG
  or GHYP or JOHNSONSU

/DISPLAY 
PLOTS=ALL or chosen from this list:
  SERIES2SD, SERIESWVARLIMITS, CONDSD, OBSACF, SQACF, ABSACF, CROSSCORR,
  RESIDDENSITY, RESQQ, RESACF, RES2ACF, NEWSIMPACT
FORECASTPLOTS=ALL or chosen from this list:
  TSUNCOND, TSROLLING, SIGMAUNCOND, SIGMAROLLING

/SAVE
RFDATASET=new dataset name  FDATASET=new dataset name
ID=variable name
FHORIZON=number  FROLL=number
WORKSPACEACTION=RETAIN or CLEAR
WORKSPACEOUTFILE="file specification"

STATS GARCH /HELP prints this help and does nothing else.

Example:  
STATS GARCH VARIABLE=y MODELSOURCE=ESTIMATE
/FIT  HOLDOUT=10
/SPECIFICATION VARMODEL=FGARCH VARGARCHORDER=2 0 SUBMODEL=TGARCH
VARTARGETING=NO MEANORDER=1 1 INCLUDEMEAN=YES DISTRIBUTION=STDT
/SAVE WORKSPACEACTION=RETAIN.

The only required specification is VARIABLE.

* indicates default.

MODELSOURCE
This command works in two modes.  If MODELSOURCE=ESTIMATE, the model is
estimated.  If it is FILE or WORKSPACE, forecasts are produced from a
previously estimated model.  With FILE, the saved model workspace is
indicated with MODELFILE.  With WORKSPACE, the workspace must have
been created in the current session and preserved with 
WORKSPACEACTION=RETAIN.

For saving a dataset, an ID variable can be specified or
START and FREQUENCY parameters can be specified for inclusion of
generated date values in the dataset.  
START has the form START=number or START=number number
depending on the time structure.  With a subperiod, such as years and
months, the second parameter indicates when within the major period
the data start, numbering from 1.  FREQUENCY is the number of minor
periods within the major period, 12 in this example.

FIT.
HOLDOUT specifies the number of cases at the end of the estimation data
to be held out from estimation for use in forecast testing.  
The default is 0.
SEED sets the random number seed and can be used to ensure that
the results can be reproduced.
SOLVER specifies which solver to use.

The remaining parameters specify details of the numerical 
solvers and apply to solver GOSOLNP.
NRESTARTS specifies the required number of solver restarts.
The default is 1. 
PARALLEL specifies parallel processing using
NCORES number of CPUs.
PKG specifies the multiprocessing package.
NSIM specifies the number of simulated parameter vectors 
to generate for NRESTARTS restarts.
The HYBRID solver tries several different solvers until 
convergence is achieved.

SPECIFICATION.
This subcommand handles the details of the model.
VARMODEL specifies the variance model.  The model types are
SGARCH: standard GARCH
IGARCH: integrated GARCH
EGARCH: exponential GARCH
GJRGARCH: asymmetric positive and negative shocks
APARCH: asymmetric power ARCH
FGARCH: family or omnibus GARCH
csGARCH: standard GARCH decomposed into permanent and 
transitory components.

If VARMODEL=FGARCH, a submodel can be specified in SUBMODEL with the
choices as listed above.
VARGARCHORDER specifies the ar and ma orders of the variance model.
VARREGRESSORS specifies additional variables for the variance model.
VARTARGETING specifies variance targeting for the
conditional variance intercept.

MEANORDER specifies the ar and ma orders for the mean model.
INCLUDEMEAN specifies whether or not to include the mean.
MEANARCHPOW specifies standard deviation or variance for
the ARCH in the mean regression.
MEANARFIMA specifies fractional differencing.
MEANARCHM specifies ARCH volatility inclusion in the 
mean regression.
MEANREGRESSORS lists the variables to be included in the
mean regression.
MEANARCHEX specifies whether to multiply the external
regressors by the conditional standard deviation.
DISTRIBUTION specifies the conditional distribution
to use for innovations.  The choices are
NORM: normal
SNORM: skew-normal
STDT: Student t
SSTD: skew Student t
GED: generalized error
SGED for skew generalized error
NORMINVG: normal inverse gaussian
GHYP: generalized hyperbolic
JOHNSONSU: Johnsons SU

DISPLAY.
This subcommand specifies which estimation-time and forecast-time
plots to produce.  Choose as many as desired.  The
names are reasonably descriptive, but actually seeing
the plots should resolve any uncertainty.

For forecast-time plots (FORECASTPLOTS) some plots will not
appear unless FROLL > 0.

SAVE.
This command specifies the datasets and files to be created.
RFDATASET specifies a dataset of residuals and fitted values.
FDATASET specifies a dataset of forecast (out of sample) values.
Dataset names must not already be in use.
ID specifies an ID variable to be included in the RFDATASET.
For forecasts, the cases are labelled as "T+1", "T+2", ...
FHORIZON specifies how far into the future to forecast.  The
default is 0.
FROLL specifies rolling, i.e., one step ahead forecasts
based on the prior cases.  The default is 0, which means
static FHORIZON forecasts.  FROLL must be <= HOLDOUT.

WORKSPACEACTION specifies whether the workspace containing
the model is retained or cleared.  If retained, it can be
used for forecasting in the current session without loading 
the model from a file.

WORKSPACEOUTFILE specifies a file where the model is saved
for future use in forecasting with MODELSOURCE=FILE.

This extension command is implemented via the R rugarch
package by Alexios Ghalanos.

'
varmodellist = list(sgarch="sGARCH", fgarch="fGARCH", egarch="eGARCH", gjrgarch="gjrGARCH", 
    aparch="apARCH", igarch="iGARCH", csgarch="csGARCH")
submodellist=list(garch="GARCH", tgarch="TGARCH", avgarch="AVGARCH", ngarch="NGARCH", 
    nagarch="NAGARCH", aparch="APARCH", gjrgarch="GJRGARCH", allgarch="ALLGARCH")
distlist = list(norm="norm", snorm="snorm", stdt="std", sstd="sstd", 
        ged="ged", sged="sged", norminvg="nig", ghyp="ghyp", johnsonsu="jsu")



### MAIN ROUTINE ###
doGARCH = function(variable=NULL, start=NULL, frequency=1, modelsource="estimate", modelfile=NULL,
    holdout=0, 
    solver="solnp", nrestarts=1, parallel=FALSE, package="snowfall", ncores=1,
    nsim=1, seed=NULL,
    varmodel="sgarch", vargarchorder=c(1,1),
    submodel=NULL, varregressors=NULL, vartargeting=FALSE, 
    meanorder=c(1,1), includemean=TRUE,
    meanarchpow="stddev", meanarfima=FALSE, meanarchex=FALSE, meanarchm=FALSE,
    meanregressors=NULL, distribution="norm", plots=list(), forecastplots=NULL, displayignore=FALSE,
    rdataset=NULL, fdataset=NULL, id=NULL, fhorizon=5, froll=0,
    workspaceaction="clear", workspaceoutfile=NULL) {
    
    # Estimate or forecast a GARCH model
    
    setuplocalization("STATS_GARCH")
    
    # A warnings proc name is associated with the regular output
    # (and the same omsid), because warnings/errors may appear in
    # a separate procedure block following the regular output
    procname=gtxt("GARCH")
    warningsprocname = gtxt("GARCH: Warnings")
    omsid="STATSGARCH"
    warns = Warn(procname=warningsprocname,omsid=omsid)

    tryCatch(library(rugarch), error=function(e){
        warns$warn(gtxtf("The R %s package is required but could not be loaded.", "rugarch"),dostop=TRUE)
        }
    )

    # error checking
    if (!is.null(fdataset) || !is.null(rdataset)) {
        alldatasets = spssdata.GetDataSetList()
        if ("*" %in% alldatasets) {
            warns$warn(gtxt("The active dataset must have a name in order to use this procedure"),
                dostop=TRUE)
        }
        if (any(c(fdataset, rdataset) %in% alldatasets)) {
            warns$warn(gtxt("The output dataset names must not already be in use"),
                dostop=TRUE)
        }
    }
    if (modelsource != "estimate" && is.null(fdataset)) {
        warns$warn(gtxt("Forecasts were requested but no forecast dataset name was specified"), dostop=TRUE)
    }
    if (modelsource == "file" && is.null(modelfile)) {
        warns$warn(gtxt("A file was specified as the model source, but no name was given"), dostop=TRUE)
    }
    if (modelsource == "workspace" && !is.null(modelfile)) {
        warns$warn(gtxt("The model source was given as workspace, but an input model file was specified"),
            dostop=TRUE)
    }
    if (!is.null(submodel) && varmodel != "fgarch") {
        warns$warn(gtxt("A submodel was specified fo a model that does not support submodels"),
            dostop=TRUE)
    }
    if (is.null(variable) && modelsource == "estimate") {
        warns$warn(gtxt("Model estimation was specified, but no variable was specified"),
            dostop=TRUE)
    }
    if (length(vargarchorder) != 2 || length(meanorder) != 2) {
        warns$warn(gtxt("The ARMA order for the variance and mean must specify two values"),
            dostop=TRUE)
    }
    if (modelsource != "estimate" && !is.null(rdataset)) {
        warns$warn(gtxt("The dataset of residuals and fitted values cannot be requested except when estimating"),
            dostop=TRUE)
    }


    varmodelx = varmodellist[[varmodel]]
    if (!is.null(submodel)) {
        submodelx = submodellist[[submodel]]
    } else {
        submodelx = NULL
    }
    distributionx = distlist[[distribution]]
    meanarchpowx = ifelse(meanarchpow=="stddev", 1, 2)
    thismodel = modelsource
    thismodelfile = modelfile
    # allargsnow is the variable set right now
    # allargsest will be the variable set at estimation time
    allargsnow = as.list(environment())

    if (modelsource != "estimate") {
        if (thismodel == "modelfile") {
            load(thismodelfile)
        }
        if (!exists("resgarch") || !("uGARCHfit" %in% class(resgarch)) ) {
            warns$warn(gtxt("The workspace does not contain a GARCH model"),
            dostop=TRUE)
        }
        allvars = c(varregressors, meanregressors)
    } else {
        allvars = c(variable, varregressors, meanregressors)
        allargsest = as.list(environment())  # current args
    }
    if (!is.null(allvars)) { # forecasting with no external variables implies no data
        dta = tryCatch(spssdata.GetDataFromSPSS(allvars, row.label=id, missingValueToNA=TRUE,
            factorMode="levels"),
            error=function(e) {warns$warn(e$message, dostop=TRUE)}
        )
    
        if (!is.null(start)) {
            if (length(start) > 2) {
                warns$warn(gtxt("The start time must have one or two values"),
                dostop=TRUE)
            }
            dta = ts(dta, start=start, frequency=frequency)
            dta = data.frame(dta, row.names=time(dta))
        }
    } else {
        dta = NULL
    }
    if (modelsource == "estimate") {
        if (nrow(dta) < 100) {
            warns$warn(gtxt("Dataset too small.  At least 100 cases are required"),
            dostop=TRUE)
        }
        if (solver != "gosolnp") {
            solveropts=NULL
        } else {
            solveropts = list(n.restarts=nrestarts, parallel=parallel, pkg=package,
                cores=ncores, n.sim=nsim, rseed=seed)
        }
        vmodel = list(model=varmodelx, garchOrder=vargarchorder,
                      submodel=submodelx,
                      variance.targeting = vartargeting)
        if (!is.null(varregressors)) {
            vmodel[["external.regressors"]] = as.matrix(dta[,c(unlist(varregressors))])
        }
        mmodel = list(armaOrder=meanorder, include.mean=includemean,
                      archm=meanarchm, archpow=meanarchpowx, arfima=meanarfima,
                      external.regressors=meanregressors)
        if (!is.null(meanregressors)) {
            mmodel[["external.regressors"]] = as.matrix(dta[,c(unlist(meanregressors))])
            if (meanarchex) {
                mmodel[["archex"]] = length(meanregressors)
            }
        }
        
        spec = ugarchspec(
            variance.model = vmodel,
            mean.model= mmodel,
            distribution.model=distributionx)

        # TODO: scale, fit control

        resgarch = tryCatch(ugarchfit(spec=spec, data=dta, out.sample=holdout,
            solver=solver, solver.control=solveropts), 
            error=function(e) {warns$warn(e$message, dostop=TRUE)})

        if (convergence(resgarch) != 0) {
            warns$warn(paste(capture.output(resgarch), collapse="\n"), dostop=TRUE)
        }
        allargsest$modeldate = date()
    }

    if (!is.null(fdataset)) {
        if (froll > allargsest$holdout) {
            warns$warn(gtxt("The roll value must be less than or equal to the holdout value"),
                dostop=TRUE)
        }
        fres = tryCatch(ugarchforecast(resgarch, n.ahead=fhorizon, n.roll=froll),
            error=function(e) {warns$warns(e$message, dostop=TRUE)}
        )
    } else {
        fres = NULL
    }
    
    displayresults(allargsest, allargsnow, resgarch, fres, warns)

    if (!is.null(fdataset) || !is.null(rdataset)) {
        savepred(allargsest, allargsnow, resgarch, fres, dta, warns)
    }
    if (!is.null(workspaceoutfile)) {
        if (modelsource != "estimate") {
            warns$warn(gtxt("Ignoring workspace save request because command is not estimating the model"),
                       dostop=FALSE)
        } else {
            rm(thismodel, thismodelfile)
            save(resgarch, allargsest, allargsnow, file=workspaceoutfile)
        }
    }
    if (workspaceaction == "clear") {
        rm(list=ls())
        rm(list=ls(envir=.GlobalEnv), envir=.GlobalEnv)
    } else {
        if (modelsource == "estimate") {
            assign("resgarch", resgarch, envir=.GlobalEnv)
            assign("allargsest", allargsest, envir=.GlobalEnv)

        }
    }
}

plotlist=list("series2sd", "serieswvarlimits", "condsd",
              "obsacf", "sqacf", "absacf", "crosscorr", "residdensity", "resqq",
              "resacf", "res2acf", "newsimpact", "all")
forecastplotlist = list("tsuncond", "tsrolling", "sigmauncond", "sigmarolling", "all")

displayresults = function(allargsest, allargsnow, res, fres, warns) {
    # display results
    # allargsest is the list of input parameters at estimation time
    # allargsnow is the list of arguments now
    # res is the estimated model
    # fres is the forecast object (may be NULL)


    
    StartProcedure(allargsnow[["procname"]], allargsnow[["omsid"]])
    
    lbls = c(
        gtxt("Variable"),
        gtxt("Model Source"),
        gtxt("GARCH Model"),
        gtxt("Variance External Regressors"),
        gtxt("Mean Model Arma Order"),
        gtxt("Submodel"),
        gtxt("Mean External Regressors"),
        gtxt("Distribution"),
        gtxt("Variance Targeting"),
        gtxt("Holdout Cases"),
        gtxt("Solver"),
        gtxt("Convergence"),
        gtxt("Log Likelihood"),
        gtxt("AIC"),
        gtxt("BIC"),
        gtxt("Shibata"),
        gtxt("Hannan-Quinn"),
        gtxt("Forecast Dataset"),
        gtxt("Residual Dataset"),
        gtxt("Workspace Action"),
        gtxt("Saved Workspace"),
        gtxt("Random Number Seed"),
        gtxt("Estimation Date")
    )
    rfit = res@fit
    rmod = res@model
    vmodel = rmod$modeldesc$vmodel
    values = c(
            allargsest$variable,
            allargsnow$modelsource,
            paste(vmodel, "(", rmod$modelinc[8],
                ",", rmod$modelin[9], ")", sep=""),  #GARCH model
            ifelse(is.null(allargsest$varregressors), gtxt("--NA--"), 
                paste(allargsest$varregressors, collapse=" ")),
            paste(allargsest$meanorder, collapse=", "),
            ifelse(vmodel == "fGARCH", rmod$modeldesc$vsubmodel,
                gtxt("--NA--")),
            ifelse(is.null(allargsest$meanregressors), gtxt("--NA--"), 
                paste(allargsest$meanregressors, collapse=" ")),
            allargsest$distribution,   #7
            ifelse(allargsest$vartargeting, gtxt("Yes"), gtxt("No")),
            allargsest$holdout,
            allargsest$solver,
            ifelse(convergence(res) ==0, gtxt("Yes"), gtxt("NO")),
            round(rfit$LLH, 6),  # log likelihood
            round(infocriteria(res)[c(1,2,3,4)], 6),
            ifelse(is.null(allargsnow$fdataset), gtxt("--NA--"), allargsnow$fdataset),
            ifelse(is.null(allargsnow$rdataset), gtxt("--NA--"), allargsnow$rdataset),
            allargsnow$workspaceaction,
            ifelse(is.null(allargsest$workspaceoutfile), gtxt("--NA--"), allargsest$workspaceoutfile),
            ifelse(is.null(allargsest$seed), gtxt("--NA--"), allargsest$seed),
            allargsest$modeldate
    )

    names(values) = lbls
    summarydf = data.frame(cbind(values))
    colnames(summarydf) = gtxt("Values")

    spsspivottable.Display(summarydf, title=gtxt("GARCH Summary"), 
                           templateName="GARCHSUMMARY",
                           caption=gtxt("Results computed by R rugarch package"),
                           isSplit=FALSE
    )

    if (allargsnow$modelsource == "estimate") {
        # notation table
        lbls = c(
            "ar1, ar2,...",
            "ma1, ma2, ...",
            "mu",
            "archm",
            "arfima",
            "mxreg1, mxreg2,...",
            "alpha1, alpha2,...",
            "beta1, beta2, ...",
            "omega",
            "vxreg1, vxreg2,...",
            "gamma1, gamm2,..."
        )
        values = c(
            gtxt("AR"),
            gtxt("MA"),
            gtxt("Mean"),
            gtxt("Archm"),
            gtxt("Arfima"),
            gtxt("External Regressor"),
            gtxt("ARCH(q)"),
            gtxt("GARCH(p)"),
            gtxt("Variance Intercept"),
            gtxt("External Regressor"),
            gtxt("Leverage")
        )
        names(values) = lbls
        summarydf = data.frame(cbind(values))
        colnames(summarydf) = gtxt("Parameter Definition")
        
        spsspivottable.Display(summarydf, title=gtxt("Output Notation"),
           rowdim=gtxt("Term"),
           hiderowdimtitle=FALSE,
           templateName="GARCHNOTATION",
           isSplit=FALSE
        )
        # coefficients
        df = data.frame(rfit$matcoef, data.frame(rfit$robust.matcoef[,c(2,3,4)]))
        names(df) = c(gtxt("Estimate"), gtxt("Std. Error"), gtxt("t"), gtxt("Sig."),
            gtxt("Robust Std. Error"), gtxt("Robust t"), gtxt("Robust Sig."))
        spsspivottable.Display(df, title=gtxt("Coefficients"),
            templateName="GARCHCOEF")
    
        # goodness-of-fit tests
        df = data.frame(gof(res, c(20,30,40,50, 100)))
        names(df) = c(gtxt("Bin Size"), gtxt("Statistic"), gtxt("Sig."))
        spsspivottable.Display(df, gtxt("Goodness of Fit Tests"),
            templateName="GARCHGOF")
    
        #Hansen-Nyblom stability test
        nyb = nyblom(res)
        df = data.frame(rbind(nyb[[1]], nyb$JointStat))
        row.names(df)[nrow(df)]= gtxt("Joint")
        names(df) = c(gtxt("Statistic"))
    
        # make significance table into caption
        caption1 = paste(gtxt("Individual Critical Values:"), 
            paste(names(nyb$IndividualCritical), nyb$IndividualCritical, collapse=" "))
        caption2 = paste(gtxt("Joint Critical Values:"), 
            paste(names(nyb$JointCritical),nyb$JointCritical, collapse=" "))
        caption = paste(caption1, caption2, sep="\n")
        spsspivottable.Display(df, title=gtxt("Hansen-Nyblom Stability Test"),
            templateName="GARCHNYBLOM",
            caption=caption
        )
    
        # signbias test of Engle and Ng (remove significance ** column)
        df = data.frame(signbias(res))[-3]
        if (length(rownames(df)) == 4) {  # expose names for translation iff there are 4
            rownames(df) = c(
                gtxt("Sign Bias"),
                gtxt("Negative Sign Bias"),
                gtxt("Positive Sign Bias"),
                gtxt("Joint Effect")
            )
        }
        names(df) = c(gtxt("t"), gtxt("Sig."))
        spsspivottable.Display(df, title=gtxt("Engle-Ng Sign Bias Test"),
            templateName="GARCHSIGNBIAS"
        )
        
            unc = uncvariance(res)
            uncm = uncmean(res)
            persist = persistence(res)
            half = halflife(res)
        df = data.frame(rbind(unc, uncm, persist, half))
        names(df) = gtxt("Statistic")
        row.names(df) = c(gtxt("Long Run Unconditional Variance"),
            gtxt("Unconditional Mean of Conditional Mean Equation"),
            gtxt("Model Persistence"),
            gtxt("Halflife of Fit Variance")
        )
        spsspivottable.Display(df, title=gtxt("Other Statistics"),
            templateName="GARCHMISC")

        # plots by number
        if ("all" %in% allargsnow$plots) {
            pnums = 1:12
        } else {
            pnums = match(allargsnow$plots, plotlist)  # can't be NA's in result
        }

        if (length(pnums) > 0) {
            for (p in pnums) {
                plot(res, which=p)
            }
        }
    }
    # forecast plots by number
    if (!is.null(fres)) {
        if ("all" %in% allargsnow$forecastplots) {
            pnums = 1:4
        } else {
            pnums = match(allargsnow$forecastplots, forecastplotlist)
        }
        # some plots will fail if roll < 5
        if (length(pnums) > 0) {
            for (p in pnums) {
                if (p %in% c(2,4) && allargsnow$froll < 5) {
                    next
                }
                tryCatch(plot(fres, which=p),
                    error = function(e) {print(e$message)}
                )
            }
        }
    }
    spsspkg.EndProcedure()
}

savepred = function(allargsest, allargsnow, res, fres, dta, warns) {
    # save residuals and fitted values
    # res is the fitted objec
    # fres is the forecasted object

    # residuals dataset structure will have an id, residuals, and fitted values columns
    # forecast dataset will have and id, predicted means, and 
    # predicted sigmas (nroll+1 columns)

    if (!is.null(allargsnow$rdataset)) {   # can't be requested at forecast time
        residdf = data.frame(residuals(res), fitted(res))
        # allow for holdout in estimation
        row.names(residdf) = row.names(dta)[1:nrow(residdf)]

        # row names will always be nominal strings
        # guessing that we need length tripling for Unicode expansion
        rnlength = max(nchar(row.names(residdf)), na.rm=TRUE) * 3
        idlabel = ""
        if (!is.null(allargsnow$id)) {
            idlabel = allargsnow$id
        } else if (!is.null(allargsnow$start)) {
            idlabel = gtxt("Time")
        }
        dict = list()
        dict[[1]] = c("ID", idlabel, rnlength, paste("A", rnlength, sep=""), "nominal")
        dict[[2]] = c("Residuals", gtxt("Residuals"), 0, "F8.2", "scale")
        dict[[3]] = c("FittedValues", gtxt("Fitted Values"), 0, "F8.2", "scale")
        dict = spssdictionary.CreateSPSSDictionary(dict)
        spssdictionary.SetDictionaryToSPSS(allargsnow$rdataset, dict)
        tryCatch(spssdata.SetDataToSPSS(allargsnow$rdataset, data.frame(row.names(residdf), residdf)),
            error=function(e) {warns$warn(e$message, dostop=TRUE)}
        )
    }
    # forecasting based on already estimated model
    if (!is.null(fres)) {
        fdf = data.frame(fitted(fres), sigma(fres))
        rnlength = max(nchar(row.names(fdf)), na.rm=TRUE) * 3
        dict = list()
        dict[[1]] = c("Time", gtxt("Time"), rnlength, paste("A", rnlength, sep=""), "nominal")
        flen2 = length(fdf)/2
        for (v in seq(1, flen2)) {
            dict[[v+1]] = c(paste("Forecast", (v-1), sep=""), 
                paste(gtxt("Forecast, roll="), (v-1), sep=""), 0, "F8.2", "scale")
            dict[[flen2 +v+1]] = c(paste("Sigma", (v-1), sep=""), 
                paste(gtxt("Sigma, roll="), (v-1), sep=""), 0, "F8.2", "scale")
        }
        dict = spssdictionary.CreateSPSSDictionary(dict)
        spssdictionary.SetDictionaryToSPSS(allargsnow$fdataset, dict)
        tryCatch(spssdata.SetDataToSPSS(allargsnow$fdataset, data.frame(row.names(fdf), fdf)),
                 error=function(e) {warns$warn(e$message, dostop=TRUE)}
        )
    }
    
    spssdictionary.EndDataStep()
}

Warn = function(procname, omsid) {
    # constructor (sort of) for message management
    lcl = list(
        procname=procname,
        omsid=omsid,
        msglist = list(),  # accumulate messages
        msgnum = 0
    )
    # This line is the key to this approach
    lcl = list2env(lcl) # makes this list into an environment

    lcl$warn = function(msg=NULL, dostop=FALSE, inproc=FALSE) {
        # Accumulate messages and, if dostop or no message, display all
        # messages and end procedure state
        # If dostop, issue a stop.

        if (!is.null(msg)) { # accumulate message
            assign("msgnum", lcl$msgnum + 1, envir=lcl)
            # There seems to be no way to update an object, only replace it
            m = lcl$msglist
            m[[lcl$msgnum]] = msg
            assign("msglist", m, envir=lcl)
        } 

        if (is.null(msg) || dostop) {
            lcl$display(inproc)  # display messages and end procedure state
            if (dostop) {
                stop(gtxt("End of procedure"), call.=FALSE)  # may result in dangling error text
            }
        }
    }
    
    lcl$display = function(inproc=FALSE) {
        # display any accumulated messages as a warnings table or as prints
        # and end procedure state, if any

        if (lcl$msgnum == 0) {   # nothing to display
            if (inproc) {
                spsspkg.EndProcedure()
            }
        } else {
            if (!inproc) {
                procok =tryCatch({
                    StartProcedure(lcl$procname, lcl$omsid)
                    TRUE
                    },
                    error = function(e) {
                        FALSE
                    }
                )
            }
            if (procok) {  # build and display a Warnings table if we can
                table = spss.BasePivotTable("Warnings ","Warnings") # do not translate this
                rowdim = BasePivotTable.Append(table,Dimension.Place.row, 
                    gtxt("Message Number"), hideName = FALSE,hideLabels = FALSE)

                for (i in 1:lcl$msgnum) {
                    rowcategory = spss.CellText.String(as.character(i))
                    BasePivotTable.SetCategories(table,rowdim,rowcategory)
                    BasePivotTable.SetCellValue(table,rowcategory, 
                        spss.CellText.String(lcl$msglist[[i]]))
                }
                spsspkg.EndProcedure()   # implies display
            } else { # can't produce a table
                for (i in 1:lcl$msgnum) {
                    print(lcl$msglist[[i]])
                }
            }
        }
    }
    return(lcl)
}

# localization initialization
setuplocalization = function(domain) {
    # find and bind translation file names
    # domain is the root name of the extension command .R file, e.g., "SPSSINC_BREUSCH_PAGAN"
    # This would be bound to root location/SPSSINC_BREUSCH_PAGAN/lang

    fpath = Find(file.exists, file.path(.libPaths(), paste(domain, ".R", sep="")))
    bindtextdomain(domain, file.path(dirname(fpath), domain, "lang"))
} 
# override for api to account for extra parameter in V19 and beyond
StartProcedure <- function(procname, omsid) {
    if (substr(spsspkg.GetSPSSVersion(),1, 2) >= 19) {
        spsspkg.StartProcedure(procname, omsid)
    }
    else {
        spsspkg.StartProcedure(omsid)
    }
}

gtxt <- function(...) {
    return(gettext(...,domain="STATS_GARCH"))
}

gtxtf <- function(...) {
    return(gettextf(...,domain="STATS_GARCH"))
}


Run = function(args) {
    #Execute the STATS GARCH command

    cmdname = args[[1]]
    args = args[[2]]

    oobj = spsspkg.Syntax(list(
        spsspkg.Template("VARIABLE", subc="", ktype="existingvarlist", var="variable"),
        spsspkg.Template("MODELSOURCE", subc="", ktype="str", var="modelsource",
            vallist=list("estimate", "file", "workspace")),
        spsspkg.Template("MODELFILE", subc="", ktype="literal", var="modelfile"),
        spsspkg.Template("START", subc="", ktype="int", var="start", islist=TRUE),
        spsspkg.Template("FREQUENCY", subc="", ktype="int", var="frequency"),
        
        spsspkg.Template("HOLDOUT", subc="FIT", ktype="int", var="holdout",
            vallist=list(0)),
        spsspkg.Template("SOLVER", subc="FIT", ktype="str", var="solver",
            vallist=list("nlminb", "solnp", "lbfgs", "gosolnp", "nloptr",
            "hybrid")),
        spsspkg.Template("NRESTARTS", subc="FIT", ktype="int", var="nrestarts",
            vallist=list(1)),
        spsspkg.Template("PARALLEL", subc="FIT", ktype="bool", var="parallel"),
        spsspkg.Template("PACKAGE", subc="FIT", ktype="str", var="package",
            vallist=list("snowfall", "multicore")),
        spsspkg.Template("NCORES", subc="FIT", ktype="int", var="ncores",
            vallist=list(1)),
        spsspkg.Template("NSIM", subc="FIT", ktype="int", var="nsim",
            vallist=list(1)),
        spsspkg.Template("SEED", subc="FIT", ktype="float", var="seed"),
        
        spsspkg.Template("VARMODEL", subc="SPECIFICATION", ktype="str", 
            var="varmodel", islist=FALSE,
            vallist=list("sgarch", "fgarch", "egarch", "gjrgarch", "aparch", "igarch", "csgarch")),
        spsspkg.Template("VARGARCHORDER", subc="SPECIFICATION", ktype="int", islist=TRUE),
        spsspkg.Template("SUBMODEL", subc="SPECIFICATION", ktype="str", var="submodel",
            vallist=list("garch", "tgarch", "avgarch", "ngarch", "nagarch", 
            "aparch", "gjrgarch", "allgarch")),
        spsspkg.Template("VARREGRESSORS", subc="SPECIFICATION", ktype="existingvarlist",
            var="varregressors", islist=TRUE),
        spsspkg.Template("VARTARGETING", subc="SPECIFICATION", ktype="bool",
            var="vartargeting"),
        spsspkg.Template("MEANORDER", subc="SPECIFICATION", ktype="int", var="meanorder",
            islist=TRUE, vallist=list(0)),
        spsspkg.Template("INCLUDEMEAN", subc="SPECIFICATION", ktype="bool",var="includemean"),
        spsspkg.Template("MEANARCHPOW", subc="SPECIFICATION", ktype="str", var="meanarchpow",
            vallist=list("stddev", "var")),
        spsspkg.Template("MEANARCHM", subc="SPECIFICATION", ktype="bool", var="meanarchm"),
        spsspkg.Template("MEANARFIMA", subc="SPECIFICATION", ktype="bool", var="meanarfima"),
        spsspkg.Template("MEANARCHEX", subc="SPECIFICATION", ktype="bool",var="meanarchex"),
        spsspkg.Template("MEANREGRESSORS", subc="SPECIFICATION", ktype="existingvarlist",
            var="meanregressors"),
        spsspkg.Template("DISTRIBUTION", subc="SPECIFICATION", ktype="str", var="distribution",
            vallist=list("norm", "snorm", "stdt", "sstd", "ged", "sged", "norminvg",
            "ghyp", "johnsonsu")),
        
        spsspkg.Template("PLOTS", subc="DISPLAY", ktype="str", var="plots", islist=TRUE,
            vallist=list("all", "series2sd", "serieswvarlimits", "condsd",
                "obsacf", "sqacf", "absacf", "crosscorr", "residdensity", "resqq",
                "resacf", "res2acf", "newsimpact")),
        spsspkg.Template("FORECASTPLOTS", "DISPLAY", ktype="str", var="forecastplots", islist=TRUE,
            vallist=list("tsuncond", "tsrolling", "sigmauncond", "sigmarolling", "all")),
        spsspkg.Template("IGNORE", subc="DISPLAY", ktype="bool", var="displayignore"),
        
        spsspkg.Template("RFDATASET", subc="SAVE", ktype="varname", var="rdataset"),
        spsspkg.Template("FDATASET", subc="SAVE", ktype="varname", var="fdataset"),
        spsspkg.Template("ID", subc="SAVE", ktype="existingvarlist", var="id"),
        spsspkg.Template("FHORIZON", subc="SAVE", ktype="int", var="fhorizon",
            vallist=list(1)),
        spsspkg.Template("FROLL", subc="SAVE", ktype="int", var="froll",
            vallist=list(0)),
        spsspkg.Template("WORKSPACEACTION", subc="SAVE", ktype="str", var="workspaceaction",
            vallist=list("retain", "clear")),
        spsspkg.Template("WORKSPACEOUTFILE", subc="SAVE", ktype="literal", 
            var="workspaceoutfile")
    ))

    # A HELP subcommand overrides all else
    if ("HELP" %in% attr(args,"names")) {
        #writeLines(helptext)
        helper(cmdname)
    }
    else {
        res <- spsspkg.processcmd(oobj, args, "doGARCH")
    }
}

helper = function(cmdname) {
    # find the html help file and display in the default browser
    # cmdname may have blanks that need to be converted to _ to match the file
    
    fn = gsub(" ", "_", cmdname, fixed=TRUE)
    thefile = Find(file.exists, file.path(.libPaths(), fn, "markdown.html"))
    if (is.null(thefile)) {
        print("Help file not found")
    } else {
        browseURL(paste("file://", thefile, sep=""))
    }
}
if (exists("spsspkg.helper")) {
assign("helper", spsspkg.helper)
}