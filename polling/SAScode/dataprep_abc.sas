/* Prep ABC/WP data	*/
/* 6/28/2022 RA		*/

options nodate nocenter nonumber ls=256 formdlim=" " formchar="|----|+|---+=|-/\<>*";

libname abc ".\data_abc";

*Read in SPSS data;
proc import datafile="./data_abc/1216 W2 AZ SPSS.sav" out=az replace;
proc import datafile="./data_abc/1216 W6 FL SPSS.sav" out=fl replace;
proc import datafile="./data_abc/1216 W5 MI SPSS.sav" out=mi replace;
proc import datafile="./data_abc/1216 W1 MN SPSS.sav" out=mn replace;
proc import datafile="./data_abc/1216 W4 NC SPSS.sav" out=nc replace;
proc import datafile="./data_abc/1216 W6 PA SPSS.sav" out=pa replace;
proc import datafile="./data_abc/1216 W5 WI SPSS.sav" out=wi replace;
run;

*Data manipulations;
data abc.abc_all(keep = weight state trumpind 
					likely
					male age1 age2 age3 age4 age5 age6 
					hsless somecoll coll postcoll white black hisp other 
					lib mod con noideo 
					demparty repparty otherparty demparty_nolean repparty_nolean otherparty_nolean);
set az(in=inAZ) fl(in=inFL) mi(in=inMI) mn(in=inMN) nc(in=inNC) pa(in=inPA) wi(in=inWI);
if inFL then state ="FL";
	else if inMI then state = "MI";
	else if inPA then state = "PA";
	else if inWI then state = "WI";
	else if inNC then state = "NC";
	else if inAZ then state = "AZ";
	else if inMN then state = "MN";
*trump indicator;
*FL, MI, NC, PA, WI = based on composite variable that includes leaners;
if state in("FL", "MI", "NC", "PA", "WI") then do; *codes as missing 7=would not vote/8=DK/no opinion;
	if q5aq6net = 1 then trumpind = 1;
		else if q5aq6net in(2,3,4,5,6) then trumpind = 0; *2=Biden/3=Jorgensen/4=Hawkins/5=Other/6=None of these;
end;
else if state in("AZ", "MN") then do; *codes as missing 5=would not vote/8=DK/no opinion;
    if q6net = 1 then trumpind = 1;
		else if q6net in(2,3,4) then trumpind = 0; *2=Biden/3=Jorgensen/4=Hawkins;
end;
*binary gender; *no missing or refused;
if q924net = 1 then male = 1;
	else if q924net = 2 then male = 0;
*age;
if q910 not in(.,99) then do; *codes as missing 99=refused and system missing;
	age1 = (18 <= q910 <= 24);
	age2 = (25 <= q910 <= 29);
	age3 = (30 <= q910 <= 39);
	age4 = (40 <= q910 <= 49);
	age5 = (50 <= q910 <= 64);
	age6 = (q910 >= 65);
end;
*education;
if educnew not in(.,8) then do; *codes as missing 8=DK/no opinion and system missing;
	hsless = (educnew = 1);
	somecoll = (educnew = 2);
	coll = (educnew = 3);
	postcoll = (educnew = 4);
end;
*race/ethnicity;
if q918 not in(.,8) then do; *codes as missing 8=DK/no opinion and system missing;
	white = (q918 = 1);
	black = (q918 = 2);
	hisp = (q918 in(3,4,5)); *3=white hisp/4=black hisp/5=hisp no race given;
	other = (q918 in(6,7)); *6=Asian, 7=Other;
end;
*likely voter;
*q905=registered voter status: 1=registered/2=No/3=will register by voting day/8=DK/no opinion;
*q2=voting likelihood: 1=absolutely certain to vote/2=will probably vote/3=chances 50-50/4=less than that/5=don't think will vote/6=already voted/8=DK/no opinion;
if q905 in (1,3) and q2 in (1,2,3,6) then likely = 1;
	else if q905 = 2 or q2 in (4,5) then likely = 0;
*ideology;
*ideo5=5-pt ideology scale: 1=very liberal/2=somewhat liberal/3=moderate/4=somewhat conservative/5=conservative/6=don't think in those terms/8=DK/no opinion;
if ideo5 > . then do;
	lib = (ideo5 in(1,2));
	mod = (ideo5 = 3);
	con = (ideo5 in(4,5));
	noideo = (ideo5 in(6,8));
end;
*party identification;
*partlean=party ID with leaners: *1=democrat/2=republican/3=independent/4=other/8=DK/no opinion;
if partlean in(1,2,3,4) then do; *codes as missing 8=DK/no opinion;
	demparty = (partlean = 1);
	repparty = (partlean = 2);
	otherparty = (partlean in(3,4));
end;
*party identification no leaners;
if q901 in(1,2,3,4) then do; *codes as missing 8=DK/no opinion;
	demparty_nolean = (q901 = 1);
	repparty_nolean = (q901 = 2);
	otherparty_nolean = (q901 in(3,4));
end;
format _ALL_;
run;
*export to CSV for R;
proc export data=abc.abc_all outfile=".\data_abc\abc_all.csv" DBMS=CSV replace;
run;

/* descriptive stats, weighted and unweighted */
/* includes estimated percentage of votes for trump among likely voters */
proc surveyfreq data=abc.abc_all;
by state;
tables likely*(trumpind--other lib--otherparty_nolean)/ row cl nototal;
*where trumpind ne . and male ne . and age1 ne . and age2 ne . and age3 ne . and age4 ne . and age5 ne . and age6 ne . 
	and hsless ne . and somecoll ne . and coll ne . and postcoll ne . 
	and white ne . and black ne . and hisp ne . and other ne . 
	and lib ne . and mod ne . and con ne . and and noideo ne . 
	and demparty ne . and repparty ne . and otherparty ne . 
	and demparty_nolean ne . and repparty_nolean ne . and otherparty_nolean ne . 
	and weight ne .; 
ods output crosstabs=stats_unwt;
run;
proc surveyfreq data=abc.abc_all;
by state;
weight weight;
tables likely*(trumpind--other lib--otherparty_nolean)/ row cl nototal;
*where trumpind ne . and male ne . and age1 ne . and age2 ne . and age3 ne . and age4 ne . and age5 ne . and age6 ne . 
	and hsless ne . and somecoll ne . and coll ne . and postcoll ne . 
	and white ne . and black ne . and hisp ne . and other ne . 
	and lib ne . and mod ne . and con ne . and and noideo ne . 
	and demparty ne . and repparty ne . and otherparty ne . 
	and demparty_nolean ne . and repparty_nolean ne . and otherparty_nolean ne . 
	and weight ne .; 
ods output crosstabs=stats_wt;
run;
data stats(where=(likely=1));
length type $15. state $2. variable $20. level 8.;
set stats_unwt(in=inU) stats_wt(in=inW);
if inU then type = "ABC Unweighted";
	else if inW then type = "ABC Weighted";
*pull varible name out from table description;
variable = scan(table,-1);
*collapse level of outcome into one variable;
level = max(of trumpind male--otherparty_nolean);
*convert percentages into decimals;
pct = RowPercent/100;
pctSE = RowStdErr/100;
pctLB = RowLowerCL/100;
pctUB = RowUpperCL/100;
keep type state variable level likely frequency WgtFreq pct:;
format _NUMERIC_;
run;
proc export data=stats outfile="..\stats_abc.csv" DBMS=CSV replace;
run;
