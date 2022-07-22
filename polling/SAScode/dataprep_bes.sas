/* Prep UK BES data	*/
/* 7/11/2022 RA		*/

options nodate nocenter nonumber ls=256 formdlim=" " formchar="|----|+|---+=|-/\<>*";

libname bes ".\data_bes";

*Read in Stata dataset;
proc import datafile="./data_bes/bes_f2f_2015_v4.0.dta" out=bes replace;
run;
*Data manipulations to prep for analysis in R;
data bes.bes(keep = weight
					likely
					male age1 age2 age3 age4 age5 age6
					vote2010_:
					region_: region5_: region7_:);
set bes;
*design variables;
weight = wt_combined_main_capped;
*likely voter indicator; *codes as missing -2=refused/-1=DK;
if b01 in(1,2) then likely = (b01 = 1);
*binary gender; *no missing values;
male = (y09 = 1);
*age; *codes as missing -1 (presumably = refused);
if age > -1 then do;*18-24, 25-34, 35-44, 45-54, 55-64, and 65-;
	age1 = (18 <= age <= 24);
	age2 = (25 <= age <= 34);
	age3 = (35 <= age <= 44);
	age4 = (45 <= age <= 54);
	age5 = (55 <= age <= 64);
	age6 = (age >= 65);
end;
/**party ID; *codes as missing -2=refused/-1=DK;*/
/*if d01 > -1 then do;*/
/*	party_cons = (d01 = 3);*/
/*	party_labour = (d01 = 2);*/
/*	party_libdem = (d01 = 4);*/
/*	party_other = (d01 in(5,6,7,8,9,10));*/
/*	party_novote = (d01 = 1);*/
/*end;*/
*2010 vote; *codes as missing -2=refused/-1=DK/11=note eligible or too young to vote;
if u05 not in(-2,-1,11) then do;
	vote2010_cons = (u05 = 3);
	vote2010_labour = (u05 = 2);
	vote2010_libdem = (u05 = 4);
	vote2010_other = (u05 in(5,6,7,8,9,10));
	vote2010_novote = (u05 = 1);
end;
*region;
region_1 = (GOR = 1); *North East;
region_2 = (GOR = 2); *North West;
region_3 = (GOR = 3); *Yorkshire and the Humber;
region_4 = (GOR = 4); *East Midlands;
region_5 = (GOR = 5);* West Midlands;
region_6 = (GOR = 6); *East of England;
region_7 = (GOR = 7); *London;
region_8 = (GOR = 8); *South East;
region_9 = (GOR = 9); *South West;
region_10 = (GOR = 11); *Wales;  ***NOTE REVERSAL OF WALES AND SCOTLAND IN BES DATA COMPARED TO POLLS;
region_11 = (GOR = 10); *Scotland;
*5-level region;
region5_1 = (GOR in(1,2,3)); *North = North East + North West + Yorkshire and the Humber;
region5_2 = (GOR in(4,5)); *Midlands = East Midlands + West Midlands;
region5_3 = (GOR in(6,7,8)); *South East = East of England + London + South East;
region5_4 = (GOR in(9,10)); *South West & Wales = South West + Wales;
region5_5 = (GOR = 11); *Scotland;
*7-level region;
region7_1 = (GOR in(1,2,3)); *North = North East + North West + Yorkshire and the Humber;
region7_2 = (GOR in(4,5)); *Midlands = East Midlands + West Midlands;
region7_3 = (GOR = 6); *Eastern;
region7_4 = (GOR = 7); *London;
region7_5 = (GOR in(8,9)); *South = South East + South West;
region7_6 = (GOR = 11); *Scotland;
region7_7 = (GOR = 10); *Wales;
run;
*Save as CSV;
proc export data=bes.bes outfile="./data_bes/bes.csv" DBMS=CSV replace;
run;

/* descriptive stats, weighted and unweighted */
/* includes estimated percentage of votes for conservative party among likely voters */
proc surveyfreq data=bes.bes;
tables likely*(male--vote2010_novote region:)/ row cl nototal;
ods output crosstabs=stats_unwt;
run;
proc surveyfreq data=bes.bes;
weight weight;
tables likely*(male--vote2010_novote region:)/ row cl nototal;
ods output crosstabs=stats_wt;
run;
data stats(where=(likely=1));
length type $20. variable $20. level 8.;
set stats_unwt(in=inU) stats_wt(in=inW);
if inU then type = "BES Unweighted";
	else if inW then type = "BES Weighted";
*pull varible name out from table description;
variable = scan(table,-1);
*collapse level of outcome into one variable;
level = max(of male age1-age6 vote2010_cons vote2010_labour vote2010_libdem vote2010_other vote2010_novote region:);
*convert percentages into decimals;
pct = RowPercent/100;
pctSE = RowStdErr/100;
pctLB = RowLowerCL/100;
pctUB = RowUpperCL/100;
keep type variable level likely frequency WgtFreq pct:;
format _NUMERIC_;
run;
proc export data=stats outfile="..\stats_bes.csv" DBMS=CSV replace;
run;
