* find a way to store the estimates from the HRS, ARIC, and meta-analyzed and then bootstrap the whole thing.
* Call with
* do prelimrct6 500
args N
drop _all
set obs `N'
generate HRSb0=.
generate HRSb1=.
generate ARICb0=.
generate ARICb1=.
generate pooledb0=.
generate pooledb1=.

local i 0

while `i'<`N' {

local i=`i'+1
preserve
clear
*set seed 128
* gen hrs
set obs 8888
gen id=_n
gen data="hrs"
gen diabetic=id<1838
gen male=runiform()<.46 if diabetic==1
replace male=runiform()<.39 if diabetic==0
*gen hba1c=7+(1.3)*rnormal() if diabetic==1
* beta mean=a/(a+b)
* beta variance=(a*b)/(a+b)^2(a+b+1)
gen hba1c=5.5+9*rbeta(1.5,5) if diabetic==1
*hist hba1c if diabetic==1
replace hba1c=5.5+(.42)*rnormal() if diabetic==0
replace hba1c=6.5 if hba1c>6.5 & diabetic==0
replace hba1c=4.5 if hba1c<4.5
replace hba1c=18 if hba1c>18
sum
gen age=67.4+8.8*rnormal() if diabetic==1
replace age=66.9+9.6*rnormal() if diabetic==0
gen age_c=(age-70)/10
gen db=hba1c>6.5
gen hba1c_c=hba1c-6.5
gen hba1c_spline=0 if hba1c<6.5
replace hba1c_spline=hba1c-6.5 if hba1c>=6.5
sum if diabetic==1
sum if diabetic==0
*hist hba1c if diabetic==1

gen ageslope=.1*rnormal()-.494
gen hba1cslope=.027
gen agehba1cslope=.09*rnormal()-.052
gen hba1csplineslope=.05*rnormal()+-.044
gen agehba1csplineslope=.1*rnormal()+.066
gen covs=.0001*rnormal() if diabetic==0
replace covs=.0001*rnormal() if diabetic==1
gen intercept=rnormal()*.5+hba1cslope*hba1c_c+hba1csplineslope*hba1c_spline+ageslope*age_c /*
  */+ agehba1cslope*age_c*hba1c_c + agehba1csplineslope*age_c*hba1c_spline

gen numobs=int(rnormal()+4.5)
replace numobs=4 if numobs>4
replace numobs=1 if numobs==0
tab numobs
sum numobs
sum 

* make longitudinal
expand numobs
sort id
by id: gen wave=_n-1
gen agenow=age_c+.2*wave
gen age_hba1c=agenow*hba1c_c
gen age_hba1c_spline=agenow*hba1c_spline
gen time=wave*.2
gen time_hba1c=time*hba1c_c
gen time_hba1cspline=time*hba1c_spline
gen timedb=time*db
gen mem_fix=intercept+time*ageslope+time_hba1c*agehba1cslope /*
  */ +time_hba1cspline*agehba1csplineslope
gen mem=mem_fix+rnormal()*.5

gen agenow2=agenow/2
gen age2_hba1c=agenow2*hba1c_c
gen age2_covs=agenow2*covs

mixed mem covs age2_covs agenow2 hba1c_c age2_hba1c  || id: if diabetic==0
matrix b0=e(b)
scalar HRSb0=b0[1,5]
dis HRSb0

mixed mem covs age2_covs agenow2 hba1c_c age2_hba1c  || id: if diabetic==1
matrix b1=e(b)
scalar HRSb1=b1[1,5]
dis HRSb1

save HRS , replace

clear
*set seed 897
set obs 12854
gen id=_n
gen data="aric"
gen diabetic=id<1293

*using hrs distbn of male
gen male=runiform()<.46 if diabetic==1
replace male=runiform()<.39 if diabetic==0
*gen hba1c=7+(1.3)*rnormal() if diabetic==1
* beta mean=a/(a+b)
* beta variance=(a*b)/(a+b)^2(a+b+1)
gen hba1c=5+15*rbeta(1.5,5) if diabetic==1
hist hba1c if diabetic==1
*replace hba1c=5.7+(.42)*rnormal() if diabetic==0
replace hba1c=4.5+15*rbeta(1.,5) if diabetic==0
*replace hba1c=6.5 if hba1c>6.5 & diabetic==0
replace hba1c=4.5 if hba1c<4.5
replace hba1c=18 if hba1c>18
sum
gen age=57.8+20*(runiform()-.64) if diabetic==1
replace age=56.5+20*(runiform()-.575) if diabetic==0
gen age_c=(age-70)/10
gen db=hba1c>6.5
gen hba1c_c=hba1c-6.5
gen hba1c_spline=0 if hba1c<6.5
replace hba1c_spline=hba1c-6.5 if hba1c>=6.5
sum if diabetic==1
sum if diabetic==0
hist hba1c if diabetic==1

gen timesp1slope=0*rnormal()+.0327247 if diabetic==0
replace timesp1slope=-.0030159 if diabetic==1

gen timesp2slope=0*rnormal()-.060707 if diabetic==0
replace timesp2slope=-.0998912 if diabetic==1

gen hba1cslope=.01*rnormal()-.0507143 if diabetic==1
replace hba1cslope=.01*rnormal()+.2151937 if diabetic==0

gen age1hba1cslope=.01*rnormal()-.0085427   if diabetic==0
replace age1hba1cslope=.01*rnormal()+.0055946 if diabetic==1

gen age2hba1cslope=.01*rnormal()-.0020503  if diabetic==0
replace age2hba1cslope=.01*rnormal()+.0004765 if diabetic==1

gen covs=hba1c_c+.25*rnormal() if diabetic==0
replace covs=hba1c_c*.7+rnormal() if diabetic==1

*ignoring age effects on intercept in aric
gen intercept=rnormal()+hba1cslope*hba1c_c 

gen numobs=int(rnormal()+2.9) if diabetic==0
replace numobs=int(rnormal()+2.3) if diabetic==1
replace numobs=3 if numobs>3
replace numobs=1 if numobs<1

* make longitudinal
expand numobs
sort id
by id: gen wave=_n-1
gen agenow=age_c if wave==0
replace agenow=age_c+.6 if wave==1
replace agenow=age_c+2.0 if wave==2

gen timesp1=0 if wave==0
replace timesp1=6 if wave==1 | wave==2
gen timesp2=0 if wave==0 | wave==1
replace timesp2=14 if wave==2

gen age_hba1c=agenow*hba1c_c
gen age_hba1c_spline=agenow*hba1c_spline
gen time=agenow-age_c
gen time_hba1c=time*hba1c_c
gen time_hba1cspline=time*hba1c_spline
gen timedb=time*db
gen mem_fix=intercept+timesp1*timesp1slope+timesp2*timesp2slope + /*
           */ timesp1*hba1c_c*age1hba1cslope + timesp2*hba1c_c*age2hba1cslope
gen mem=mem_fix+rnormal()*.5

gen agenow2=agenow/2
gen age2_hba1c=agenow2*hba1c_c
gen age2_covs=agenow2*covs
mixed mem covs age2_covs agenow2 hba1c_c age2_hba1c  || id:  if diabetic==0
matrix b0=e(b)
scalar ARICb0=b0[1,5]

mixed mem covs age2_covs agenow2 hba1c_c age2_hba1c  || id:  if diabetic==1
matrix b1=e(b)
scalar ARICb1=b1[1,5]

replace id=id+8889
save ARIC , replace

append using HRS
gen HRS=(data=="hrs")
gen HRS_agenow2=HRS*agenow2
gen HRS_hba1c=HRS*hba1c_c
gen HRS_covs=HRS*covs
gen HRS_age2_covs=HRS*age2_covs
tab data
sum 
mixed mem covs age2_covs HRS_covs HRS_age2_covs HRS HRS_agenow2 HRS_hba1c agenow2 hba1c_c age2_hba1c || id: if diabetic==0
matrix b0=e(b)
scalar pooledb0=b0[1,10]
dis pooledb0

mixed mem covs age2_covs HRS_covs HRS_age2_covs HRS HRS_agenow2 HRS_hba1c agenow2 hba1c_c age2_hba1c  || id: if diabetic==1
matrix b1=e(b)
scalar pooledb1=b1[1,10]
dis pooledb1

restore
quietly replace HRSb0=scalar(HRSb0) in `i'
quietly replace HRSb1=scalar(HRSb1) in `i'
quietly replace ARICb0=scalar(ARICb0) in `i'
quietly replace ARICb1=scalar(ARICb1) in `i'
quietly replace pooledb0=scalar(pooledb0) in `i'
quietly replace pooledb1=scalar(pooledb1) in `i'
noisily dis in red `i'
}
*/
