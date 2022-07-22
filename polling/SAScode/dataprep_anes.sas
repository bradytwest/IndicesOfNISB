/* Prep ANES data	*/
/* 7/8/2022 RA		*/

options nodate nocenter nonumber ls=256 formdlim=" " formchar="|----|+|---+=|-/\<>*";

libname anes ".\data_anes";

*Read in SPSS dataset;
proc import datafile="./data_anes/anes_timeseries_2020_spss_20220210.sav" out=anes replace;
run;
data anes.anes_pre_all(keep = weight stratum cluster state state_registered
				registered likely
				male age1 age2 age3 age4 age5 age6 
				hsless somecoll coll postcoll white black hisp other 
				lib con mod noideo 
				demparty repparty otherparty demparty_nolean repparty_nolean otherparty_nolean
				where=(state in("AK","AZ","FL","MI","MN","NC","PA","WI")));
length state state_registered $2.;
set anes;
*design variables;
weight = V200010a; /* full sample ANES weight for pre-election in 2020 */
stratum = V200010d;
cluster = V200010c;
statefip = V203000; *state from sample info;
if 1 <= V201014b <= 56 then statefip_registered = V201014b; *state the person is registered in;
	else statefip_registered = statefip; *use sample state if not registered;
state = fipstate(statefip);
state_registered = fipstate(statefip_registered);
*PRE-ELECTION VARIABLES;
*voting registration status;
if V201025x in(2,3,4) then registered = 1; *codes as missing -4=technical error; *2=not registered but plans to / 3=registered but hasn't voted yet / 4=registered and already voted early;
	else if V201025x = 1 then registered = 0; *1=not registered and doesn't plan to;
*likely voting status;
if V201025x = 4 then likely = 1; *already voted early;
	else if V201100 in(1,2,3) then likely = 1; *1=extremely likely/2=very likely/3=moderately likely;
	else if V201100 in(4,5) then likely = 0; *4=slightly likely/5=not likely at all;
*follow politics; *codes as missing -9=refused;
if V201005 in(1,2,3) then followpolitics = 1; *1=always/2=most of the time/3=about half the time;
	else if V201005 in(4,5) then followpolitics = 0; *4=some of the time/5=never;
*binary gender;
if V201600 = 1 then male = 1; *codes as missing -9=refused;
	else if V201600 = 2 then male = 0;
*age;
if V201507x > -9 then do; *code as missing -9=refused;
	age1 = (18 <= V201507x <= 24);
	age2 = (25 <= V201507x <= 29);
	age3 = (30 <= V201507x <= 39);
	age4 = (40 <= V201507x <= 49);
	age5 = (50 <= V201507x <= 64);
	age6 = (V201507x >= 65);
end;
*education;
if V201511x not in(-9,-8,-2) then do; *codes as missing -9=refused/-8=DK/-2=missing/other specify not coded;
	hsless = (V201511x in(1,2));
	somecoll = (V201511x = 3);
	coll = (V201511x = 4);
	postcoll = (V201511x = 5);
end;
*race/ethnicity;
if V201549x not in(-9,-8) then do; *codes as missing -9=refused/-8=DK;
	white = (V201549x = 1);
	black = (V201549x = 2);
	hisp = (V201549x = 3);
	other = (V201549x in(4,5,6)); *4=Asian/NHOPI non-hispanic, 5=Native American/Alaska Native or other race non-hispanic, 6=multiple races non-hispanic;
end;
*ideology;
*V201200=7 pt ideology scale: 1=extremely lib/2=lib/3=slightly lib/4=moderate/5=slightly con/6=con/7=extremely con/99=haven't thought about it much/-9=refused/-8=DK;
if V201200 in(1,2,3,4,5,6,7,99) then do; *codes as missing -9=refused/-8=DK;
	lib = (V201200 in(1,2,3));
	mod = (V201200 = 4);
	con = (V201200 in(5,6,7));
	noideo = (V201200= 99);
end;
*party identification;
*V201231x: 1=strong dem/2=not v strong dem/3=indep dem/4=indep/5=indep-repub/6=not v strong repub/7=strong repub/-9=refused/-8=DK;
if V201231x in(1,2,3,4,5,6,7) then do;
	demparty = (V201231x in(1,2,3));
	repparty = (V201231x in(5,6,7));
	otherparty = (V201231x = 4);
end;
*party identification no leaners;
*V201228: 1=dem/2=repub/3=ind/5=other/0=no pref/-9=ref/-8=DK/-4=tech error;
if V201228 in(1,2,3,5,0) then do;
	demparty_nolean = (V201228 = 1);
	repparty_nolean = (V201228 = 2);
	otherparty_nolean = (V201228 in(3,5,0));
end;
proc sort;
by state state_registered;
run;

*Save as CSV;
proc export data=anes.anes_pre_all outfile="./data_anes/anes_pre_all.csv" DBMS=CSV replace;
run;

/* descriptive stats, weighted and unweighted */
proc freq data=anes.anes_pre_all;
table state;
run;
/********
                                  Cumulative    Cumulative
state    Frequency     Percent     Frequency      Percent
----------------------------------------------------------
AK              9        0.45             9         0.45
AZ            177        8.81           186         9.26
FL            499       24.85           685        34.11
MI            289       14.39           974        48.51
MN            183        9.11          1157        57.62
NC            287       14.29          1444        71.91
PA            359       17.88          1803        89.79
WI            205       10.21          2008       100.00
********/
/* Not enough sample in AK, so exclude */
proc surveyfreq data=anes.anes_pre_all(where=(state ~= "AK"));
by state;
tables likely*(male--otherparty_nolean)/ row cl nototal;
ods output crosstabs=stats_unwt;
run;
proc surveyfreq data=anes.anes_pre_all(where=(state ~= "AK"));
by state;
weight weight;
tables likely*(male--otherparty_nolean)/ row cl nototal;
ods output crosstabs=stats_wt;
run;
data stats;
length type $15. state $2. variable $20. level 8.;
set stats_unwt(in=inU) stats_wt(in=inW);
if inU then type = "ANES Unweighted";
	else if inW then type = "ANES Weighted";
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
proc export data=stats outfile="..\stats_anes_pre.csv" DBMS=CSV replace;
run;
