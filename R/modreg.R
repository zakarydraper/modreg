modreg <- function(dv, idv, mod1, mod2 = NULL, multMOD1 = 1, multMOD2 = 1, ci = 95, xlab = NULL, ylab = NULL) {

	# one continuous moderator.

	if (is.null(mod2) && !is.factor(mod1)) {

		datam <- data.frame(idv, mod1, idv*mod1, dv)

		# n, mean, sd, & correlation matrix

		# vcv <- var(datam)

		npr <- ncol(datam) - 1
		n <- nrow(datam)
		sd <- apply(datam, 2, sd)
		mn <- apply(datam, 2, mean)
		lm1 <- lm(dv ~ idv * mod1)
		lm2 <- lm(dv ~ idv + mod1)
		aov <- anova(lm1, lm2)
		cr <- cor(datam)

		beta <- solve(cr[1:npr, 1:npr]) %*% cr[1:npr, npr + 1]
		b <- lm1$coef[2:length(lm1$coef)]

		# Overall regression coefficients.
		a <- lm1$coef[1]
		r2chXn <- summary(lm1)$r.squared - summary(lm2)$r.squared
		fsquared <- (summary(lm1)$r.squared - summary(lm2)$r.squared) / (1 - summary(lm1)$r.squared)
		F <- aov$F[2]
		dferror <- summary(lm1)$df[2]
		pF <- summary(lm1)$coef[npr+1,4]

		# simple slopes info.
		modlo <- mean(mod1) - (sd(mod1) * multMOD1)
		modmd <- mean(mod1)
		modhi <- mean(mod1) + (sd(mod1) * multMOD1)

		# intercepts.
		aslopes <- rbind((b[2] * modlo + a), (b[2] * modmd + a), (b[2] * modhi + a))

		# raw b.
		slopes <- rbind((b[1] + (b[3] * modlo)), (b[1] + (b[3] * modmd)), (b[1] + (b[3] * modhi)))

		# t-values.
		mse <- (n/(n - npr)) * (sd[(npr + 1)]^2) * (1 - summary(lm1)$r.squared)
		Sb <- as.numeric(mse) * solve((diag(sd[1:npr]) %*% cr[1:npr, 1:npr] %*% diag(sd[1:npr])) * (n - 1))

		SEslopes <- rbind((sqrt(cbind(1, 0, modlo) %*% Sb %*% t(cbind(1, 0, modlo)))),
		 	 (sqrt(cbind(1, 0, modmd) %*% Sb %*% t(cbind(1, 0, modmd)))),
		 	 (sqrt(cbind(1, 0, modhi) %*% Sb %*% t(cbind(1, 0, modhi)))))
		tslopes <- slopes/SEslopes

		# dfs
		df <- c(rep(summary(lm1)$df[2],3))

		# p-values.
		pslopes <- (1 - pt(abs(tslopes), summary(lm1)$df[2])) * 2

		# standardized beta weights.
		zslopes <- slopes * (sd[1]/sd[4])

		# standard error.
		zSE <- SEslopes * (sd[1]/sd[4])

		# Confidence Intervals for standardized slopes.
		ci2 <- 1 - ((1 - (ci/100))/2)
		tabledT <- qt(ci2,summary(lm1)$df[2])
		confidlo <- zslopes - tabledT * zSE
		confidhi <- zslopes + tabledT * zSE

		# Printing results.

		grp.labels <- c("Mod=low", "Mod=med", "Mod=high")

		coeffs <- cbind(r2chXn, F, 1, summary(lm1)$df[2], fsquared, pF)
		colnames(coeffs) <- c("Rsq. ch.", "F", "df num.", "df denom.", "fsquared", "Sig. F")

		betaw <- cbind(b, beta)
		rownames(betaw) <- c("idv", "mod", "Xn")
		colnames(betaw) <- c("raw b", "std. beta")

		slopecoeffs <- cbind(aslopes, slopes, tslopes, df, pslopes)
		rownames(slopecoeffs) <- grp.labels
		colnames(slopecoeffs) <- c("a", "raw b", "t test", "df", "Sig. t")

		z.s.slopes <- cbind(zslopes, zSE, confidlo, confidhi)
		rownames(z.s.slopes) <- grp.labels
		colnames(z.s.slopes) <- c("std. beta", "SE", paste0(ci,"% Low"), paste0(ci, "% High"))

		s.slope.dv <- (b[1]/b[3]) * (-1)

		s.reglines <- (b[2]/b[3]) * (-1)

		plot(idv, dv, xlab = xlab, ylab = ylab)
		abline(aslopes[1, 1], slopes[1, 1])
		abline(aslopes[2, 1], slopes[2, 1])
		abline(aslopes[3, 1], slopes[3, 1])

		out <- list(lm = coeffs,
		            beta = betaw,
		            slopes = slopecoeffs,
		            std_slopes = z.s.slopes,
		            slope_dv = s.slope.dv,
		            reg_lines = s.reglines)

		return(out)

	}  # End is.null(mod2) condition.

  # Any categorical variables or two moderator variables.

  if (!is.null(mod2) || "factor" %in% sapply(data.frame(idv, dv, mod1), class)) {

    stop("You have either included a factor variable or two moderator variables. This function currently only supports models with one continuous outcome, one continuous predictor, and one continuous moderator variable. Support for two moderators and categorical predictor and moderator variables will be added in a future release.")

    }  # End categorical/two moderator conditional.

}  # End function.
