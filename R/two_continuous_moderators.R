simple_slope_2 <- function(idv, dv, mod1, mod2 = NULL, multMOD1 = 1, multMOD2 = 1, multiIDV = 2, ci = 95, xlab = NULL, ylab = NULL) {



	# two continuous moderators.
	if (!is.null(mod2) && !is.factor(mod1) && !is.factor(mod2)) {
		datam <- data.frame(idv, mod1, mod2,
							idv*mod1, idv*mod2, mod1*mod2,
							idv*mod1*mod2,
							dv)

		n <- nrow(datam)
		sd <- apply(datam, 2, sd)
		mn <- apply(datam, 2, mean)
		cr <- cor(datam)
		npr <- ncol(datam) - 1
		nvar <- ncol(datam)

		lm1 <- lm(dv ~ idv * mod1 * mod2)
		lm2 <- lm(dv ~ idv + mod1 + mod2)

		aov <- anova(lm1, lm2)

		beta <- solve(cr[1:npr, 1:npr]) %*% cr[1:npr, nvar]
		b <- lm1$coef[2:nvar]

		# Overall regression coefficients.
		a <- lm1$coeff[1]
		r2all <- t(beta) %*% cr[1:npr, nvar]

		# for the X1 interaction, with main effects, X2, & X3 in the equation.

		if (!is.factor(mod1) && !is.factor(mod2)) {
			j1 <- 3
			j2 <- 4
			j3 <- 5
			j4 <- 6
		} else if (is.factor(mod1) && !is.factor(mod2)) {
			j1 <- 3 + ndc1 - 1		# main effects
			j2 <- j1 + ndc1			# main effects, idv-md1
			j3 <- j1 + ndc1 + 1		# main effects, idv-md1, idv-md2
			j4 <- j3 + ndc1			# main effects, idv-md1, idv-md2, md1-md2
		} else if (is.factor(mod1) && is.factor(mod2)) {
			nmain <- 3 + (ndc1 - 1) + (ndc2 - 1)
		}


		r2 <- t( (solve(cr[1:(npr - 1), 1:(npr - 1)]) %*% cr[1:(npr - 1), nvar]) ) %*% cr[1:(npr - 1), nvar]
		cr1 <- rbind( cbind(cr[1:j1, 1:j1], cr[1:j1, j3:j4]), cbind(cr[j3:j4, 1:j1], cr[j3:j4, j3:j4]) )
		r21 <- t( solve(cr1) %*% c(cr[1:j1,nvar],cr[j3:j4,nvar]) ) %*% c(cr[1:j1,nvar],cr[j3:j4,nvar])
		r2chX1 <- r2 - r21
		fsquare1 <- (r2 - r21) / (1 - r2)
		df6 <- summary(lm1)$df[2]
		F1 <- ( (r2 - r21) / 1 ) / ( (1 - r2) / (df6) )
		pF1 <- 1 - pf(F1, 1, df6)

		# for the X2 interaction, with main effects, X1, & X3 in the equation.
		cr2 <- rbind( cbind(cr[1:j2, 1:j2], cr[1:j2, j4]), c(cr[j4, 1:j2], cr[j4,j4]) )
		r22 <- t( solve(cr2) %*% c(cr[1:j2, nvar], cr[j4, nvar]) ) %*% c(cr[1:j2, nvar], cr[j4, nvar])
		r2chX2 <- r2 - r22
		fsquare2 <- (r2chX2) / (1 - r2)
		F2 <- ( (r2chX2) / 1 ) / ( (1 - r2) / df6 )
		pF2 <- 1 - pf(F2, 1, df6)

		# for the X3 interaction, with main effects, X1, & X2 in the equation.
		r23 <- t( solve(cr[1:j3, 1:j3]) %*% cr[1:j3, nvar] ) %*% cr[1:j3, nvar]
		r2chX3 <- r2 - r23
		fsquare3 <- r2chX3 / (1 - r2)
		F3 <- ( (r2chX3) / 1 ) / ( (1 - r2) / df6 )
		pF3 <- 1 - pf(F3, 1, df6)

		# for X4, the 3-way interaction.
		r2main <- t(solve(cr[1:j4, 1:j4]) %*% cr[1:j4, nvar]) %*% cr[1:j4, nvar]
		r2chXn <- r2all - r2main
		fsquared <- r2chXn / (1 - r2all)
		dferror <- n - 7 - 1
		F <- ( r2chXn / 1 ) / ( (1 - r2all) / dferror )
		pF <- 1 - pf(F, 1, dferror)

		if (!is.factor(mod1) && !is.factor(mod2)) {
			mod1lo <- mean(mod1) - (sd(mod1) * multMOD1)
			mod1md <- mean(mod1)
			mod1hi <- mean(mod1) + (sd(mod1) * multMOD1)

			mod2lo <- mean(mod2) - (sd(mod2) * multMOD2)
			mod2hi <- mean(mod2) + (sd(mod2) * multMOD2)
			slopes <- rbind(
			 (b[1] + b[4] * mod1lo + b[5] * mod2lo + b[7] * mod1lo * mod2lo),
			 (b[1] + b[4] * mod1lo + b[5] * mod2hi + b[7] * mod1lo * mod2hi),
			 (b[1] + b[4] * mod1hi + b[5] * mod2lo + b[7] * mod1hi * mod2lo),
			 (b[1] + b[4] * mod1hi + b[5] * mod2hi + b[7] * mod1hi * mod2hi))
			ints <- rbind(
			 b[2] * mod1lo + b[3] * mod2lo + b[6] * mod1lo * mod2lo + a,
			 b[2] * mod1lo + b[3] * mod2hi + b[6] * mod1lo * mod2hi + a,
			 b[2] * mod1hi + b[3] * mod2lo + b[6] * mod1hi * mod2lo + a,
			 b[2] * mod1hi + b[3] * mod2hi + b[6] * mod1hi * mod2hi + a)
		} else if (is.factor(mod1) && !is.factor(mod2)) {
			mod2lo <- mean(mod2) - (sd(mod2) * multMOD2)
			mod2hi <- mean(mod2) + (sd(mod2) * multMOD2)

			slopes <- b[1] + b[5 + (ndc1 - 1)*2]*mod2lo
			slopes <- c(slopes, (b[1] + b[5 + (ndc1 - 1)*2]*mod2hi))
			for (i in 1:(ndc1)) {
				slopes <- c(slopes,
				 (b[1] + b[4 + (ndc1 - 1) + (i - 1)] + b[5 + (ndc1 - 1)*2]*mod2lo + b[7 + (ndc1 - 1)*3 + (i - 1)]*mod2lo)
				)
			}
			for (i in 1:(ndc1)) {
				slopes <- c(slopes,
				 (b[1] + b[4 + (ndc1 - 1) + (i - 1)] + b[5 + (ndc1 - 1)*2]*mod2hi + b[7 + (ndc1 - 1)*3 + (i - 1)]*mod2hi)
				)
			}
		} else if (is.factor(mod1) && is.factor(mod2)) {
			print("This feature is not currently implemented.")
		}

		mse <- (n / (n - 7)) * sd[nvar]^2 * (1 - r2all)
		Sb <- as.numeric(mse) * solve((diag(sd[1:npr]) %*% cr[1:npr, 1:npr] %*% diag(sd[1:npr])) * (n - 1))
		# STOP

		if (!is.factor(mod1) && !is.factor(mod2)) {
			SEslopes <- rbind(
			 sqrt (cbind(1, 0, 0, mod1lo, mod2lo, 0, mod1lo%*%mod2lo) %*% Sb %*%
			 t(cbind(1, 0, 0, mod1lo, mod2lo, 0, mod1lo%*%mod2lo))),
			 sqrt (cbind(1, 0, 0, mod1lo, mod2hi, 0, mod1lo%*%mod2hi) %*% Sb %*%
			 t(cbind(1, 0, 0, mod1lo, mod2hi, 0, mod1lo%*%mod2hi))),
			 sqrt (cbind(1, 0, 0, mod1hi, mod2lo, 0, mod1hi%*%mod2lo) %*% Sb %*%
			 t(cbind(1, 0, 0, mod1hi, mod2lo, 0, mod1hi%*%mod2lo))),
			 sqrt (cbind(1, 0, 0, mod1hi, mod2hi, 0, mod1hi%*%mod2hi) %*% Sb %*%
			 t(cbind(1, 0, 0, mod1hi, mod2hi, 0, mod1hi%*%mod2hi))))
		} else if (is.factor(mod1) && !is.factor(mod2)) {
			mod1lo <- 0
			mod1hi <- 1

			m <- matrix(0, ndc1*2 + 1, nvar) # ?????
			rownames(m) <- c("idv:mod1g1_g2", "idv:mod1g1_g3", "idv:mod1g1_g4", "idv:mod2", "3way1", "3way2", "3way3")
			m[,1] <- 1
			m[1,] <- c(1, rep(0, 10), mod1lo, mod2lo, 0, mod1lo*mod2lo)

		}
		# dfs
		dfs <- n - npr - 1
		df <- rbind(dfs, dfs, dfs, dfs)

		tslopes <- slopes/SEslopes


		# p-values.
		pslopes <- (1 - pt(abs(tslopes), dfs)) * 2

		# standardized beta weights.
		zslopes <- slopes * (sd[1]/sd[nvar])

		# standard error.
		zSE <- SEslopes * (sd[1]/sd[nvar])

		# Confidence Intervals for standardized slopes.
		ci2 <- 1 - ((1 - (ci/100))/2)
		tabledT <- qt(ci2,dfs)
		confidlo <- zslopes - tabledT * zSE
		confidhi <- zslopes + tabledT * zSE

		# printing results.

		coeffs2way <- cbind(rbind(r2chX1, r2chX2, r2chX3), rbind(F1, F2, F3), rbind(1,1,1), rbind(df6,df6,df6), rbind(fsquare1, fsquare2, fsquare3), rbind(pF1, pF2, pF3))
		rownames(coeffs2way) <- c("idv-md1", "idv-md2", "md2-md3")
		colnames(coeffs2way) <- c("Rsq. ch.", "F", "df num.", "df denom.", "fsquared", "Sig. F")
		cat("\n\nCoefficients for the 2-way Interactions:\n\n")
		print(coeffs2way)

		coeffs3way <- cbind(r2chXn, F, 1, dferror, fsquared, pF)
		colnames(coeffs3way) <- c("Rsq. ch.", "F", "df num.", "df denom.", "fsquared", "Sig. F")
		cat("\n\nCoefficients for the 3-way Interaction:\n\n")
		print(coeffs3way)

		betaw <- cbind(b, beta)
		rownames(betaw) <- c("idv", "md1", "md2", "idv-md1", "idv-md2", "md2-md3", "3-way")
		colnames(betaw) <- c("raw b", "std.beta")
		cat("\n\nBeta weights for the full interaction:\n\n")
		print(betaw)

		cat("\n\nThe intercept is:\n\n", a)

		# simple slopes for when Moderator 1 is low.
		mod1low <- cbind(ints[1:2, ], slopes[1:2, ], tslopes[1:2, ], df[1:2, ], pslopes[1:2, ])
		rownames(mod1low) <- c("md2=low", "md2=high")
		colnames(mod1low) <- c("a", "raw b", "t test", "df", "Sig. T")
		cat("\n\nSimple Slope Coefficients at 2 levels of Moderator 2 when moderator 1 is low:\n\n")
		print(mod1low)

		ssci <- cbind(zslopes[1:2], zSE[1:2], confidlo[1:2], confidhi[1:2])
		rownames(ssci) <- c("md2=low", "md2=hi")
		colnames(ssci) <- c("std.beta", "SE", paste0(ci,"% Low"), paste0(ci,"% Hi"))
		cat("\n\nStandardized Simple Slopes & Confidence Intervals for the above:\n\n")
		print(ssci)

		ssint <- as.numeric((b[3] + (b[6] * mod1lo)) / (b[5] + (b[7] * mod1lo)) * (-1))
		cat("\n\nWhen Mod1=low, the simple slopes at Mod2=high and Mod2=low intersect at IDV =\n\n")
		print(ssint)

		ssflat <- as.numeric((b[1] + b[4] * mod1lo) / (b[5] + b[7] * mod1lo) * (-1))
		cat("\n\nWhen Mod1=low the simple slope for the DV on the IDV is zero (flat) at Mod 2 =\n\n")
		print(ssflat)

		# simple slopes for when Moderator 1 is high.
		mod1high <- cbind(ints[3:4, ], slopes[3:4, ], tslopes[3:4, ], df[3:4, ], pslopes[3:4, ])
		rownames(mod1high) <- c("md2=low", "md2=high")
		colnames(mod1high) <- c("a", "raw b", "t test", "df", "Sig. T")
		cat("\n\nSimple Slope Coefficients at 2 levels of Moderator 2 when Moderator 1 is high:\n\n")
		print(mod1high)

		ssci2 <- cbind(zslopes[1:2], zSE[1:2], confidlo[1:2], confidhi[1:2])
		rownames(ssci2) <- c("md2=low", "md2=hi")
		colnames(ssci2) <- c("std.beta", "SE", paste0(ci, "% Low"), paste0(ci, "% Hi"))
		cat("\n\nStandardized Simple Slopes & Confidence Intervals for the above:\n\n")
		print(ssci2)

		ssint2 <- as.numeric((b[3] + (b[6] * mod1hi)) / (b[5] + (b[7] * mod1hi)) * (-1))
		cat("\n\nWhen Mod1=high, the simple slopes at Mod2=high and Mod2=low intersect at IDV =\n\n")
		print(ssint2)

		ssflat2 <- as.numeric((b[1] + b[4] * mod1hi) / (b[5] + b[7] * mod1hi) * (-1))
		cat("\n\nWhen Mod1=high the simple slope for the DV on the IDV is zero (flat) at Mod 2 =\n\n")
		print(ssflat2)

		plot(idv, dv, xlab = xlab, ylab = ylab)
		abline

	}  # End !is.null(mod2) conditional.

}  # End function.
