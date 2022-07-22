/* Prep Roper Poll data	*/
/* 6/28/2022 RA			*/

options nodate nocenter nonumber ls=256 formdlim=" " formchar="|----|+|---+=|-/\<>*";

libname roper ".\data_roper";

*Read in Stata datasets;
proc import datafile="./data_roper/AZ_31118620.dta" out=az replace;
proc import datafile="./data_roper/AK_31118619.dta" out=ak replace;
run;
*Data manipulations to prep for analysis in R;
data roper.roper_all(keep = weight state trumpind 
					likely
					male age1 age2 age3 age4 age5 age6 hsless somecoll coll postcoll 
					white black hisp other 
					demparty_nolean repparty_nolean otherparty_nolean);
set ak(drop=age4) az(drop=age4);
*design variables;
weight = FINAL_WE;
*trump indicator; *codes as missing 7=not voting/9=don't know/refused;
*NOTE: QHR2 = voting preference including people who initially said dk/refused (QPRES) but endorsed leaning towards one candidate (QPRESLN);
if QHR2 = 2 then trumpind = 1;
	else if QHR2 in(1,3,6) then trumpind = 0; *1=Biden/3=Jorgenson/4=Someone else;
*binary gender;
if GENDER = 1 then male = 1;
	else if GENDER = 2 then male = 0; *codes as missing 3=other/refused;
*age;
if AGE > 0 then do; *code as missing 0=refused;
	age1 = (18 <= AGE <= 24);
	age2 = (25 <= AGE <= 29);
	age3 = (30 <= AGE <= 39);
	age4 = (40 <= AGE <= 49);
	age5 = (50 <= AGE <= 64);
	age6 = (AGE >= 65);
end;
*education;
if QEDUC4 ~= 9 then do; *codes as missing 9=refused;
	hsless = (QEDUC4 = 1);
	somecoll = (QEDUC4 = 2);
	coll = (QEDUC4 = 3);
	postcoll = (QEDUC4 = 4);
end;
*race/ethnicity;
if QRACER ~= 9 then do; *codes as missing 9=refused;
	white = (QRACER = 1);
	black = (QRACER = 2);
	hisp = (QRACER = 3);
	other = (QRACER in(4,5,6)); *4=Asian, 5=Other, 6="6" only in AK data, probably AK natives?;
end;
*likely voter;
if QLIKELY ~= 9 then do; *codes as missing 9=don't know/refused;
	if QLIKELY in(1,2,3,6) then likely = 1; *1=Almost certain/2=Very likely/3=Somewhat likely/6=Already voted (6 is only in AK data);
		else if QLIKELY in(4,5) then likely = 0; *4=Not very likely/5=Not at all likely;
end;
*party identification;
*QPARTYID = party ID with no leaner question: 1=democrat/2=republican/3=independnt/4=another party/9=DK/refused;
if QPARTYID in(1,2,3,4) then do;
	demparty_nolean = (QPARTYID = 1); *codes 9=don't know/refused as 0s for all variables;
	repparty_nolean = (QPARTYID = 2);
	otherparty_nolean = (QPARTYID in(3,4));
end;
run;
*Save as CSV;
proc export data=roper.roper_all outfile="./data_roper/roper_all.csv" DBMS=CSV replace;
run;


/* descriptive stats, weighted and unweighted */
/* includes estimated percentage of votes for trump among likely voters */
proc surveyfreq data=roper.roper_all;
by state;
tables likely*(trumpind--other demparty_nolean--otherparty_nolean)/ row cl nototal;
*where trumpind ne . and male ne . and age1 ne . and age2 ne . and age3 ne . and age4 ne . and age5 ne . and age6 ne . 
	and hsless ne . and somecoll ne . and coll ne . and postcoll ne . 
	and white ne . and black ne . and hisp ne . and other ne . 
	and demparty_nolean ne . and repparty_nolean ne . and otherparty_nolean ne . and weight ne .; 
ods output crosstabs=stats_unwt;
run;
proc surveyfreq data=roper.roper_all;
by state;
weight weight;
tables likely*(trumpind--other demparty_nolean--otherparty_nolean)/ row cl nototal;
*where trumpind ne . and male ne . and age1 ne . and age2 ne . and age3 ne . and age4 ne . and age5 ne . and age6 ne . 
	and hsless ne . and somecoll ne . and coll ne . and postcoll ne . 
	and white ne . and black ne . and hisp ne . and other ne . 
	and demparty_nolean ne . and repparty_nolean ne . and otherparty_nolean ne . and weight ne .; 
ods output crosstabs=stats_wt;
run;
data stats(where=(likely=1));
length type $20. state $2. variable $20. level 8.;
set stats_unwt(in=inU) stats_wt(in=inW);
if inU then type = "Roper Unweighted";
	else if inW then type = "Roper Weighted";
*pull varible name out from table description;
variable = scan(table,-1);
*collapse level of outcome into one variable;
level = max(of trumpind male age1-age6 hsless somecoll coll postcoll white black hisp other demparty_nolean repparty_nolean otherparty_nolean);
*convert percentages into decimals;
pct = RowPercent/100;
pctSE = RowStdErr/100;
pctLB = RowLowerCL/100;
pctUB = RowUpperCL/100;
keep type state variable level likely frequency WgtFreq pct:;
format _NUMERIC_;
run;
proc export data=stats outfile="..\stats_roper.csv" DBMS=CSV replace;
run;

/* descriptive stats, weighted and unweighted */
/* using all records (not subset to likely voters) */
proc surveyfreq data=roper.roper_all;
by state;
tables trumpind--other demparty_nolean--otherparty_nolean/ cl nototal;
ods output oneway=stats_unwt;
run;
proc surveyfreq data=roper.roper_all;
by state;
weight weight;
tables trumpind--other demparty_nolean--otherparty_nolean/ cl nototal;
ods output oneway=stats_wt;
run;
data stats;
length type $22. state $2. variable $20. level 8.;
set stats_unwt(in=inU) stats_wt(in=inW);
if inU then type = "Roper (all) Unweighted";
	else if inW then type = "Roper (all) Weighted";
*pull varible name out from table description;
variable = scan(table,-1);
*collapse level of outcome into one variable;
level = max(of trumpind male age1-age6 hsless somecoll coll postcoll white black hisp other demparty_nolean repparty_nolean otherparty_nolean);
*convert percentages into decimals;
pct = Percent/100;
pctSE = StdErr/100;
pctLB = LowerCL/100;
pctUB = UpperCL/100;
keep type state variable level frequency WgtFreq pct:;
format _NUMERIC_;
run;
proc export data=stats outfile="..\stats_roper_allresp.csv" DBMS=CSV replace;
run;
