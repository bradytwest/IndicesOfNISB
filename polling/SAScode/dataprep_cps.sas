/* Prep CP8S data	*/
/* 6/28/2022 RA		*/

options nodate nocenter nonumber ls=256 formdlim=" " formchar="|----|+|---+=|-/\<>*";

libname cps ".\data_cps";

data cps.cps_all(keep = weight state
				likely
				male age1 age2 age3 age4 age5 age6 
				hsless somecoll coll postcoll white black hisp other
			where=(likely > . and state in("AK","AZ","FL","MI","MN","NC","PA","WI")));
length state $2.;
set cps.cps_nov2020;
*design variables;
state = fipstate(statefip);
weight = VOSUPPWT; /* voter supplement weight (same as WTFINL) */
*voting status; *variable called "likely" but really it is self-report of whether voted in Nov 2020 election;
if VOTED = 2 then likely = 1; *codes as missing 96=refused/97=DK/98=not reported/not available/99=not in universe;
	else if VOTED = 1 then likely = 0;
*binary gender;
male = (SEX = 1); *no missing data;
*age;
age1 = (18 <= AGE <= 24); *no missing data - ages < 18 will be excluded in analysis;
age2 = (25 <= AGE <= 29);
age3 = (30 <= AGE <= 39);
age4 = (40 <= AGE <= 49);
age5 = (50 <= AGE <= 64);
age6 = (AGE >= 65);
*education; *no missing data;
if EDUC not in(1,999) then do; *codes as missing 1=not in universe or blank/999=missing/unknown;
	hsless = (2 <= EDUC <= 73);
	somecoll = (EDUC in(81,91,92)); *81=some college but no degree/91=associates degree vocational/92=assoc degree academic;
	coll = (EDUC = 111);
	postcoll = (EDUC in(123,124,125)); *123=masters/124=professional/125=doctoral;
end;
*race/ethnicity;
if HISPAN not in(901,902) then not_hispanic = (HISPAN = 0);
if RACE ~= 999 then do; *codes as missing 999=blank;
	white = (RACE = 100 and not_hispanic = 1);
	black = (RACE = 200 and not_hispanic = 1);
	hisp = (not_hispanic = 0);
	other = (300 <= RACE <= 830 and not_hispanic = 1);
end;
proc sort;
by state;
run;
*Save as CSV;
proc export data=cps.cps_all outfile="./data_cps/cps_all.csv" DBMS=CSV replace;
run;

/* descriptive stats, weighted and unweighted */
proc freq data=cps.cps_all;
table state;
run;
/********
                                  Cumulative    Cumulative
state    Frequency     Percent     Frequency      Percent
----------------------------------------------------------
AK            767        6.38           767         6.38
AZ           1194        9.94          1961        16.32
FL           2902       24.15          4863        40.47
MI           1596       13.28          6459        53.75
MN            970        8.07          7429        61.82
NC           1502       12.50          8931        74.32
PA           1975       16.44         10906        90.75
WI           1111        9.25         12017       100.00
********/
proc surveyfreq data=cps.cps_all;
by state;
tables likely*(male--other)/ row cl nototal;
ods output crosstabs=stats_unwt;
run;
proc surveyfreq data=cps.cps_all;
by state;
weight weight;
tables likely*(male--other)/ row cl nototal;
ods output crosstabs=stats_wt;
run;
data stats;
length type $15. state $2. variable $20. level 8.;
set stats_unwt(in=inU) stats_wt(in=inW);
if inU then type = "CPS Unweighted";
	else if inW then type = "CPS Weighted";
*pull varible name out from table description;
variable = scan(table,-1);
*collapse level of outcome into one variable;
level = max(of male age1--other);
*convert percentages into decimals;
pct = RowPercent/100;
pctSE = RowStdErr/100;
pctLB = RowLowerCL/100;
pctUB = RowUpperCL/100;
keep type state variable level likely frequency WgtFreq pct:;
format _NUMERIC_;
run;
proc export data=stats outfile="..\stats_cps.csv" DBMS=CSV replace;
run;
