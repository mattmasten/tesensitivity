capture log close
log using "vignette", smcl replace
//_1
sysuse lalonde1986, clear
//_2
local outcome "re78"
local treatment "treat"
local controls "married age black hispanic education re74 re75 re74pos re75pos"
//_3
teffects ipw (`outcome') (`treatment' `controls') if sample1
//_4
teffects ipw (`outcome') (`treatment' `controls') if sample3
//_5
tesensitivity cpi (`outcome' `controls') (`treatment' `controls') if sample1, ate
//_6
tesensitivity cpiplot
//_7q
qui graph export ate1.png, width(500) replace
//_8
tesensitivity cscale
//_9
tesensitivity cscale education, density 
//_10q
qui graph export edu_density.png, width(500) replace
//_11
tesensitivity cpi (`outcome' `controls') (`treatment' `controls') if sample1, ate cref
//_12
tesensitivity cpiplot, creflines
//_13q
qui graph export ate_reflines.png, width(500) replace
//_14
estat summarize
//_15
qui tesensitivity cpi (`outcome' `controls') (`treatment' `controls') if sample1, ate
estimates store ate1
estimates replay ate1
    
//_16
tesensitivity cpi (`outcome' `controls') (`treatment' `controls') if sample3, ate
estimates store ate3
//_17
tesensitivity cpitable ate1 ate3
//_18
tesensitivity cpiplot ate1 ate3
//_19q
qui graph export ate_compare.png, width(500) replace
//_20
qui tesensitivity cpi (`outcome' `controls') (`treatment' `controls') if sample1, atet
qui estimates store atet1
qui tesensitivity cpi (`outcome' `controls') (`treatment' `controls') if sample3, atet
qui estimates store atet3
tesensitivity cpiplot atet1 atet3
//_21
qui graph export ATET_compare.png, width(500) replace
//_22
qui tesensitivity cpi (`outcome' `controls') (`treatment' `controls') if sample1, qte
qui estimates store qte1_50
tesensitivity cpitable ate1 qte1_50
tesensitivity cpiplot ate1 qte1_50
//_23q
qui graph export ate_qte_compare.png, width(500) replace
//_24
matrix cov = (1)
matrix colnames cov = married
qui tesensitivity cpi (`outcome' `controls') (`treatment' `controls') ///
    if sample1, cate median cov(cov)
qui estimates store cate1_50

matrix cov[1,1] = 0
matrix qcov = (.9)
matrix colnames qcov = re74
qui tesensitivity cpi (`outcome' `controls') (`treatment' `controls') ///
    if sample1, cate median qcov(qcov) cov(cov)
qui estimates store cate1_90 

tesensitivity cpiplot cate1_50 cate1_90 
//_25q
qui graph export cate_compare.png, width(500) replace
//_26
estimates restore cate1_50
matrix list e(covsupp)
    
estimates restore cate1_90
matrix list e(covsupp)
//_27
matrix cov = (1)
matrix colnames cov = married
qui tesensitivity cpi (`outcome' `controls') (`treatment' `controls') ///
    if sample1, cqte median cov(cov)
qui estimates store cqte1_50

matrix cov[1,1] = 0
matrix qcov = (.9)
matrix colnames qcov = re74
qui tesensitivity cpi (`outcome' `controls') (`treatment' `controls') ///
    if sample1, cqte median qcov(qcov) cov(cov)
qui estimates store cqte1_90 

tesensitivity cpiplot cqte1_50 cqte1_90 
//_28q
qui graph export cqte_compare.png, width(500) replace
//_29
gen re78pos = re78 > 0
tesensitivity cpi (re78pos `controls') (`treatment' `controls') if sample1, ate
//_30
tesensitivity cpi (re78pos `controls') (`treatment' `controls') if sample1, ate verbose
//_31
qui tesensitivity cpiplot cqte1_50 cqte1_90, ///
                  xtitle(c-dependence) ytitle(conditional quantile treatment effect) ///
                  graphregion(color(%8) margin(vsmall)) ///
                  title(CQTE at 90th and 50th percentiles of income)
//_32q
qui graph export cqte_compare_reformat1.png, width(500) replace
//_33
qui tesensitivity cpiplot cqte1_50 cqte1_90, ///
                  xtitle(c-dependence) ytitle(conditional quantile treatment effect) ///
                  graphregion(color(%8) margin(vsmall)) ///
                  title(CQTE at 20th and 50th percentiles of income) ///
                  boundcolors(navy ltblue) boundpatterns(solid) boundoptions(lwidth(vthick))
//_34q
qui graph export cqte_compare_reformat2.png, width(500) replace
//_35
qui tesensitivity cpiplot cqte1_50 cqte1_90, ///
                  xtitle(c-dependence) ytitle(conditional quantile treatment effect) ///
                  graphregion(color(%8) margin(vsmall)) ///
                  title(CQTE at 20th and 50th percentiles of income) ///
                  boundcolors(navy ltblue) boundpatterns(solid) boundoptions(lwidth(vthick)) ///
                  breakdownoptions(lcolor(orange))
//_36q
qui graph export cqte_compare_reformat3.png, width(500) replace
//_37
qui tesensitivity cpiplot cqte1_50 cqte1_90, ///
                  xtitle(c-dependence) ytitle(conditional quantile treatment effect) ///
                  graphregion(color(%8) margin(vsmall)) ///
                  title(CQTE at 20th and 50th percentiles of income) ///
                  boundcolors(navy ltblue) boundpatterns(solid) boundoptions(lwidth(vthick)) ///
                  breakdownoptions(lcolor(orange)) nolegend
//_38q
qui graph export cqte_compare_reformat4.png, width(500) replace
//_39
qui tesensitivity cpiplot cqte1_50 cqte1_90, ///
                  xtitle(c-dependence) ytitle(conditional quantile treatment effect) ///
                  graphregion(color(%8) margin(vsmall)) ///
                  title(CQTE at 20th and 50th percentiles of income) ///
                  boundcolors(navy ltblue) boundpatterns(solid) boundoptions(lwidth(vthick)) ///
                  breakdownoptions(lcolor(orange)) legoptions(region(fcolor(%8)))
//_40q
qui graph export cqte_compare_reformat5.png, width(500) replace
//_^
log close
