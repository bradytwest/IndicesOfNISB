/* Prep UK Poll data	*/
/* 7/11/2022 RA			*/

options nodate nocenter nonumber ls=256 formdlim=" " formchar="|----|+|---+=|-/\<>*";

libname uk ".\data_UKpolls";

*Read in Stata dataset;
proc import datafile="./data_UKpolls/PooledPolls.dta" out=ukpolls replace;
run;
*Data manipulations to prep for analysis in R;
data uk.ukpolls_all(keep = weight poll company consvote 
					male age1 age2 age3 age4 age5 age6
					vote2010_:
					region_: region5_: region7_:);
set ukpolls(keep=WfullO poll company vote turnout gender agegroup region: vote2010);
*design variables;
weight = WfullO;
*conservative vote indicator;
if vote > . then consvote = (vote = 1);
*binary gender; *no missing values;
male = (gender = 1);
*age; *no missing values;
age1 = (agegroup = 1);
age2 = (agegroup = 2);
age3 = (agegroup = 3);
age4 = (agegroup = 4);
age5 = (agegroup = 5);
age6 = (agegroup = 6);
*party ID based on 2010 vote;
if vote2010 > . then do;
	vote2010_cons = (vote2010 = 1);
	vote2010_labour = (vote2010 = 2);
	vote2010_libdem = (vote2010 = 3);
	vote2010_other = (vote2010 in(4,5,6,7,8));
	vote2010_novote = (vote2010 = 9);
end;
*region;
if region > . then do;
	region_1 = (region = 1); *North East;
	region_2 = (region = 2); *North West;
	region_3 = (region = 3); *Yorkshire and the Humber;
	region_4 = (region = 4); *East Midlands;
	region_5 = (region = 5);* West Midlands;
	region_6 = (region = 6); *East of England;
	region_7 = (region = 7); *London;
	region_8 = (region = 8); *South East;
	region_9 = (region = 9); *South West;
	region_10 = (region = 10); *Wales;
	region_11 = (region = 11); *Scotland;
end;
if region5 > . then do;
	region5_1 = (region5 = 1); *North;
	region5_2 = (region5 = 2); *Midlands;
	region5_3 = (region5 = 3); *South East;
	region5_4 = (region5 = 4); *South West & Wales;
	region5_5 = (region5 = 5); *Scotland;
end;
if region7 > . then do;
	region7_1 = (region7 = 1); *North;
	region7_2 = (region7 = 2); *Midlands;
	region7_3 = (region7 = 3); *Eastern;
	region7_4 = (region7 = 4); *London;
	region7_5 = (region7 = 5); *South;
	region7_6 = (region7 = 6); *Scotland;
	region7_7 = (region7 = 7); *Wales;
end;
run;
*Save as CSV;
proc export data=uk.ukpolls_all outfile="./data_UKpolls/UKpolls_all.csv" DBMS=CSV replace;
run;

/* descriptive stats, weighted and unweighted */
/* using all records with weight > 0 (not subset to likely voters) */
/* includes estimated percentage of votes for conservative party among likely voters */
proc surveyfreq data=uk.ukpolls_all;
by poll;
where weight > 0;
tables consvote--vote2010_novote region:/ cl nototal;
ods output oneway=stats_unwt;
run;
proc surveyfreq data=uk.ukpolls_all;
by poll;
where weight > 0;
weight weight;
tables consvote--vote2010_novote region:/ cl nototal;
ods output oneway=stats_wt;
run;
data stats;
length type $25. poll 8. variable $20. level 8.;
set stats_unwt(in=inU) stats_wt(in=inW);
if inU then type = "UK Polls Unweighted";
	else if inW then type = "UK Polls Weighted";
*pull varible name out from table description;
variable = scan(table,-1);
*collapse level of outcome into one variable;
level = max(of consvote male age1-age6 vote2010_cons vote2010_labour vote2010_libdem vote2010_other vote2010_novote region:);
*convert percentages into decimals;
pct = Percent/100;
pctSE = StdErr/100;
pctLB = LowerCL/100;
pctUB = UpperCL/100;
keep type poll variable level frequency WgtFreq pct:;
format _NUMERIC_;
run;
proc export data=stats outfile="..\stats_ukpolls.csv" DBMS=CSV replace;
run;

*check sample sizes for region variables;
data temp;
set ukpolls(where=(Poll in(13,24,33,43,53,63,73,83,93)));
run;
proc freq data=temp;
table Company*region /*Company*region5*/ Company*region7/list missing;
run;

proc surveyfreq data=ukpolls;
where poll=73;
weight WfullO;
table vote;
run;
