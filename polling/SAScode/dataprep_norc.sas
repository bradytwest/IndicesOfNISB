/* Prep AP/NORC data	*/
/* 6/29/2022 RA			*/

options nodate nocenter nonumber ls=256 formdlim=" " formchar="|----|+|---+=|-/\<>*";

libname apnorc ".\data_apnorc";

data apnorc.norc_all(keep = weight state 
				likely
				male age1 age2 age3 age4 age5 age6 
				hsless somecoll coll postcoll white black hisp other 
				lib mod con noideo 
				demparty repparty otherparty demparty_nolean repparty_nolean otherparty_nolean
				rename=(p_state=state));
set apnorc.AP_VOTECAST_2020_GENERAL(where=(p_state in("FL","MI","PA","WI","AZ","MN","NC","AK")));
*design variables;
weight = finalvote_state_weight;
state = p_state;
*likely voter; *no missing/DK/ref data for this variable;
likely = (likelyvoter = 1);
*binary gender;
if gender = 1 then male = 1;
	else if gender = 2 then male = 0;
*age;
age = age65;
if age = 99 then age = .;
if age = 1 then age1 = 1; /* 18-24 */
	else if age ne . then age1 = 0;
if age = 2 then age2 = 1; /* 25-29 */
	else if age ne . then age2 = 0;
if age = 3 then age3 = 1; /* 30-39 */
	else if age ne . then age3 = 0;
if age = 4 then age4 = 1; /* 40-49 */
	else if age ne . then age4 = 0;
if age = 5 then age5 = 1; /* 50-64 */
	else if age ne . then age5 = 0;
if age = 6 then age6 = 1; /* 64+ */
	else if age ne . then age6 = 0;
*education;
if educ ~= 99 then do;
	hsless = (educ = 1);
	somecoll = (educ = 2);
	coll = (educ = 3);
	postcoll = (educ = 4);
end;
*race;
if raceth5 not in (88,99) then do;
	white = (raceth5 = 1);
	black = (raceth5 = 2);
	hisp = (raceth5 = 3);
	other = (raceth5 in (4,5));
end;
*political ideology;
*ideo=5-pt ideology scale: 1=very liberal/2=somewhat liberal/3=moderate/4=somewhat conservative/5=conservative/99=refused;
lib = (ideo in (1,2)); *very liberal/liberal;
mod = (ideo = 3);
con = (ideo in(4,5));
noideo = (ideo = 99);
*party ID with leaners;
demparty = (partyfull = 1); *there are no missing/refused, only 3 options of democrat/republican/indep;
repparty = (partyfull = 2);
otherparty = (partyfull = 3);
*party ID without leaners;
if party ~= 99 then do; *codes as missing 99=refused;
	demparty_nolean = (party = 1);
	repparty_nolean = (party = 2);
	otherparty_nolean = (party = 3);
end;
label p_state = " ";
format _NUMERIC_;
proc sort; 
by state;
run;

*Save as CSV;
proc export data=apnorc.norc_all outfile="./data_apnorc/norc_all.csv" DBMS=CSV replace;
run;

/* descriptive stats, weighted and unweighted */
proc freq data=apnorc.norc_all;
where likely = 1;
table state;
run;
/**********
                                  Cumulative    Cumulative
state    Frequency     Percent     Frequency      Percent
----------------------------------------------------------
AK            696        2.53           696         2.53
AZ           3835       13.92          4531        16.45
FL           3895       14.14          8426        30.59
MI           3713       13.48         12139        44.07
MN           3703       13.44         15842        57.51
NC           3862       14.02         19704        71.53
PA           4248       15.42         23952        86.96
WI           3593       13.04         27545       100.00
**********/
proc surveyfreq data=apnorc.norc_all;
by state;
tables likely*(male--otherparty_nolean)/ row cl nototal;
ods output crosstabs=stats_unwt;
run;
proc surveyfreq data=apnorc.norc_all;
by state;
weight weight;
tables likely*(male--otherparty_nolean)/ row cl nototal;
ods output crosstabs=stats_wt;
run;
data stats;
length type $15. state $2. variable $20. level 8.;
set stats_unwt(in=inU) stats_wt(in=inW);
if inU then type = "NORC Unweighted";
	else if inW then type = "NORC Weighted";
*pull varible name out from table description;
variable = scan(table,-1);
*collapse level of outcome into one variable;
level = max(of male age1--otherparty_nolean);
*convert percentages into decimals;
pct = RowPercent/100;
pctSE = RowStdErr/100;
pctLB = RowLowerCL/100;
pctUB = RowUpperCL/100;
keep type state variable level likely frequency WgtFreq pct:;
format _NUMERIC_;
run;
proc export data=stats outfile="..\stats_norc.csv" DBMS=CSV replace;
run;
