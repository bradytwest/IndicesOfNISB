options nodate nocenter nonumber ls=256 formdlim=" " formchar="|----|+|---+=|-/\<>*";

libname abc "X:\Brady\Selection Bias R21 Grant\Paper 5 - 2020 Election Outcomes\data_abc";
libname apnorc "X:\Brady\Selection Bias R21 Grant\Paper 5 - 2020 Election Outcomes\data_apnorc";

/****************************************/
/* create ABC polling data 				*/
/****************************************/
data abc.all (keep = state trumpind male age1 age2 age3 age4 age5 age6 hsless somecoll coll postcoll 
					white black hisp other lib noideo demparty likely weight);
set abc.florida(in=inFL) abc.michigan(in=inMI) abc.penn(in=inPA) abc.wisconsin(in=inWI)
   abc.nc(in=inNC) abc.arizona(in=inAZ) abc.minn(in=inMN);
if inFL then state ="FL";
	else if inMI then state = "MI";
	else if inPA then state = "PA";
	else if inWI then state = "WI";
	else if inNC then state = "NC";
	else if inAZ then state = "AZ";
	else if inMN then state = "MN";
if q5aq6net = 1 then trumpind = 1;
	else if q5aq6net >= 2 and q5aq6net <= 6 then trumpind = 0;
if inAZ or inMN then do;
    if q6net = 1 then trumpind = 1;
	else if q6net >= 2 and q6net <= 4 then trumpind = 0;
end; 
if q924net = 1 then male = 1;
	else if q924net = 2 then male = 0;
age = q910;
if age = 99 then age = .;
if 18 <= age <= 24 then age1 = 1;
	else if age ne . then age1 = 0;
if 25 <= age <= 29 then age2 = 1;
	else if age ne . then age2 = 0;
if 30 <= age <= 39 then age3 = 1;
	else if age ne . then age3 = 0;
if 40 <= age <= 49 then age4 = 1;
	else if age ne . then age4 = 0;
if 50 <= age <= 64 then age5 = 1;
	else if age ne . then age5 = 0;
if 65 <= age then age6 = 1;
	else if age ne . then age6 = 0;
if educnew = 1 then hsless = 1;
	else if educnew ne 8 then hsless = 0;
if educnew = 2 then somecoll = 1;
	else if educnew ne 8 then somecoll = 0;
if educnew = 3 then coll = 1;
	else if educnew ne 8 then coll = 0;
if educnew = 4 then postcoll = 1;
	else if educnew ne 8 then postcoll = 0;
if q918 = 1 then white = 1;
	else if q918 not in (8, 9) then white = 0;
if q918 = 2 then black = 1;
	else if q918 not in (8, 9) then black = 0;
if q918 = 3 then whiteh = 1;
	else if q918 not in (8, 9) then whiteh = 0;
if q918 = 4 then blackh = 1;
	else if q918 not in (8, 9) then blackh = 0;
if q918 = 5 then hispNR = 1;
	else if q918 not in (8, 9) then hispNR = 0;
if whiteh = 1 or blackh = 1 or hispNR = 1 then hisp = 1;
	else if whiteh = 0 and blackh = 0 and hispNR = 0 then hisp = 0;
if q918 in (6,7) then other = 1;
	else if q918 not in (8, 9) then other = 0;
if q905 in (1,3) and q2 in (1,2,3,6) then likely = 1;
	else if q905 = 2 or q2 in (4,5) then likely = 0;
if ideo5 in (1,2) then lib = 1;
	else if ideo5 in (3,4,5,6,8) then lib = 0;
if ideo5 in (6,8) then noideo = 1;
	else if ideo5 in (1,2,3,4,5) then noideo = 0;
if partlean = 1 then demparty = 1;
	else if partlean in (2,3,4,8) then demparty = 0;
run;

*export to CSV for R;
proc export data=abc.all outfile="X:\Brady\Selection Bias R21 Grant\Paper 5 - 2020 Election Outcomes\data_abc\abc_all.csv" DBMS=CSV replace;
run;

proc sort data = abc.all;
by state;
run;

proc freq data = abc.all;
   tables state;
run;

/* descriptive stats, weighted and unweighted */
/* includes estimated percentage of votes for trump among likely voters */
proc surveyfreq data=abc.all;
by state;
tables likely*(trumpind--other lib--demparty)/ row cl nototal;
where trumpind ne . and male ne . and age1 ne . and age2 ne . and age3 ne . and age4 ne . and age5 ne . 
   and age6 ne . and hsless ne . and somecoll ne . and coll ne . and postcoll ne . 
   and white ne . and black ne . and hisp ne . and other ne . and lib ne . and noideo ne . 
   and demparty ne . and weight ne .; 
ods output crosstabs=stats_unwt;
run;
proc surveyfreq data=abc.all;
by state;
weight weight;
tables likely*(trumpind--other lib--demparty)/ row cl nototal;
where trumpind ne . and male ne . and age1 ne . and age2 ne . and age3 ne . and age4 ne . and age5 ne . 
   and age6 ne . and hsless ne . and somecoll ne . and coll ne . and postcoll ne . 
   and white ne . and black ne . and hisp ne . and other ne . and lib ne . and noideo ne . 
   and demparty ne . and weight ne .;
ods output crosstabs=stats_wt;
run;
data stats(where=(likely=1));
length type $15. state $2. variable $10. level 8.;
set stats_unwt(in=inU) stats_wt(in=inW);
if inU then type = "ABC Unweighted";
	else if inW then type = "ABC Weighted";
*pull varible name out from table description;
variable = scan(table,-1);
*collapse level of outcome into one variable;
level = max(of trumpind male--demparty);
*convert percentages into decimals;
pct = RowPercent/100;
pctSE = RowStdErr/100;
pctLB = RowLowerCL/100;
pctUB = RowUpperCL/100;
keep type state variable level likely frequency WgtFreq pct:;
format _NUMERIC_;
run;
proc export data=stats outfile="X:\Brady\Selection Bias R21 Grant\Paper 5 - 2020 Election Outcomes\abc_stats.csv" DBMS=CSV replace;
run;
proc datasets nolist;
delete stats:;
quit;
/* evaluate ability of covariates to predict Trump indicator */
*ods html file="X:\Brady\Selection Bias R21 Grant\Paper 5 - 2020 Election Outcomes\abc_logistic_models.xls";
	proc logistic data = abc.all;
	by state;
	model trumpind (event = "1") = male age1 age2 age3 age4 age5 somecoll coll postcoll black hisp other lib noideo demparty / rsq;
	where likely = 1;
	run;
	proc surveylogistic data = abc.all;
	model trumpind (event = "1") = male age1 age2 age3 age4 age5 somecoll coll postcoll black hisp other lib noideo demparty;
	weight weight;
	where likely = 1;
	run;
*ods html close;

/****************************************/
/* create AP/NORC VoteCast data			*/
/****************************************/
options nofmterr;

proc freq data = apnorc.AP_VOTECAST_2020_GENERAL;
tables raceth5;
run;

data apnorc.norc_all(keep = finalvote_state_weight p_state male age1 age2 age3 age4 age5 age6 
				hsless somecoll coll postcoll white black hisp other lib noideo demparty
				rename=(p_state=state));
set apnorc.AP_VOTECAST_2020_GENERAL(where=(likelyvoter=1 and p_state in("FL","MI","PA","WI","AZ","MN","NC")));
if gender = 1 then male = 1;
	else if gender = 2 then male = 0;
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
if educ ~= 99 then do;
	hsless = (educ = 1);
	somecoll = (educ = 2);
	coll = (educ = 3);
	postcoll = (educ = 4);
end;
if raceth5 not in (88,99) then do;
	white = (raceth5 = 1);
	black = (raceth5 = 2);
	hisp = (raceth5 = 3);
	other = (raceth5 in (4,5));
end;
lib = (ideo in (1,2)); *very liberal/liberal;
noideo = (ideo = 99); *refused; 
demparty = (partyfull = 1);
label FINALVOTE_STATE_WEIGHT = " " p_state = " ";
format _NUMERIC_;
proc sort; 
by state;
run;
*export to CSV for R;
proc export data=apnorc.norc_all outfile="X:\Brady\Selection Bias R21 Grant\Paper 5 - 2020 Election Outcomes\data_apnorc\norc_all.csv" DBMS=CSV replace;
run;

/* descriptive stats, weighted and unweighted */
ods trace on;
proc surveyfreq data=apnorc.norc_all;
by state;
tables male--demparty/ cl nototal;
ods output oneway=stats_unwt;
run;
proc surveyfreq data=apnorc.norc_all;
by state;
weight finalvote_state_weight;
tables male--demparty/ cl nototal;
ods output oneway=stats_wt;
run;
data stats;
length type $15. state $2. variable $10. level 8.;
set stats_unwt(in=inU) stats_wt(in=inW);
if inU then type = "NORC Unweighted";
	else if inW then type = "NORC Weighted";
*pull varible name out from table description;
variable = scan(table,-1);
*collapse level of outcome into one variable;
level = max(of male age1--demparty);
*convert percentages into decimals;
pct = Percent/100;
pctSE = StdErr/100;
pctLB = LowerCL/100;
pctUB = UpperCL/100;
keep type state variable level frequency WgtFreq pct:;
format _NUMERIC_;
run;
proc export data=stats outfile="X:\Brady\Selection Bias R21 Grant\Paper 5 - 2020 Election Outcomes\norc_stats.csv" DBMS=CSV replace;
run;
proc datasets nolist;
delete stats:;
quit;
