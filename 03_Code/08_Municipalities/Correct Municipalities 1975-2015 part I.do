
** replace wrong codings by BfS in the years 1975-2015

replace municipalityno=96 if name=="Kuhn" & firstname=="Peter" & municipality=="Adlikon/Regensdorf"& year==1987	 
replace municipalityno=115 if name=="Brugger" & firstname=="Kurt" & year==1987	 
replace municipalityno=115 if name=="Romann" & firstname=="Christine" & year==1987	 

* replace checked codings of BfS (see file Quality Check Gemeinde für Florence.xlsx)

replace municipalityno=4024 if  firstname=="Beda"& name=="Humbel" & votes==32405
replace municipalityno=4027 if  firstname=="Otto"& name=="Winter" & votes==4402
replace municipalityno=4027 if  firstname=="Kenny"& name=="Bhatti" & votes==831
replace municipalityno=4100 if  firstname=="Claudia"& name=="Leuenberger" & votes==1017
replace municipalityno=4100 if  firstname=="Claudia"& name=="Leuenberger" & votes==1327
replace municipalityno=4312 if  firstname=="Mierta"& name=="Bremer-Clavuot" & votes==1391
replace municipalityno=4312 if  firstname=="Mario"& name=="Schönenberger" & votes==12737
replace municipalityno=4312 if  firstname=="Werner"& name=="Laube" & votes==34569
replace municipalityno=4312 if  firstname=="Kurt"& name=="Schmid" & votes==6245
replace municipalityno=4312 if  firstname=="Julia"& name=="Suter" & votes==1807
replace municipalityno=4235 if  firstname=="Ralf"& name=="Bucher" & votes==20336
replace municipalityno=4207 if  firstname=="Viktor"& name=="Andermatt" & votes==827
replace municipalityno=4288 if firstname=="Reinhard"& name=="Müller" & votes==30160
replace municipalityno=4182 if  firstname=="Anita"& name=="Lenzin-Tschudi" & votes==13285
replace municipalityno=4182 if  firstname=="Peter"& name=="Bircher" & votes==27210
replace municipalityno=561 if  firstname=="Fred"& name=="Rubi" & votes==105397
replace municipalityno=561 if  firstname=="Konrad"& name=="Hari" & votes==10549
replace municipalityno=630 if  firstname=="Christian"& name=="Daumüller" & votes==572
replace municipalityno=630 if  firstname=="Daniel"& name=="Pulfer" & votes==4689
replace municipalityno=371 if  firstname=="Pia"& name=="Riedwyl-Lörtscher" & votes==7261
replace municipalityno=384 if  firstname=="Marc"& name=="Schwab" & votes==3964
replace municipalityno=763 if  firstname=="Arnold"& name=="Zum Wald" & votes==21923
replace municipalityno=590 if  firstname=="Paul"& name=="Günter" & votes==56493
replace municipalityno=590 if  firstname=="Paul"& name=="Günter" & votes==63946
replace municipalityno=327 if  firstname=="Maria Magdalena"& name=="Moser" & votes==7477
replace municipalityno=956 if  firstname=="Hans"& name=="Grunder" & votes==104042
replace municipalityno=929 if  firstname=="Peter"& name=="Dütschler" & votes==23796
replace municipalityno=929 if  firstname=="Susanne"& name=="Gerber" & votes==982
replace municipalityno=929 if  firstname=="Katharina"& name=="Berger" & votes==13554
replace municipalityno=868 if firstname=="Thomas"& name=="Feuz" & votes==1434
replace municipalityno=870 if  firstname=="Heinz"& name=="Landolf" & votes==20841
replace municipalityno=870 if  firstname=="Linus"& name=="Zimmermann" & votes==64257
replace municipalityno=870 if  firstname=="Hans-Ueli"& name=="Aebi" & votes==21656
replace municipalityno=723 if  firstname=="Paul-Emile"& name=="Marti" & votes==19784
replace municipalityno=723 if  firstname=="Martin"& name=="Lehmann" & votes==7148
replace municipalityno=723 if  firstname=="Jacques"& name=="Hirt" & votes==5843
replace municipalityno=723 if  firstname=="Roland"& name=="Maître" & votes==5635
replace municipalityno=723 if  firstname=="Jean-Jacques"& name=="Kellerhals" & votes==5033
replace municipalityno=726 if  firstname=="Marc"& name=="Früh" & votes==8169
replace municipalityno=698 if firstname=="René"& name=="Schaller" & votes==17363
replace municipalityno=934 if  firstname=="Patrick"& name=="Minder" & votes==2423
replace municipalityno=934 if  firstname=="Thomas"& name=="Heuberger" & votes==20797
replace municipalityno=766 if  firstname=="Andreas"& name=="Gafner" & votes==9439
replace municipalityno=567 if  firstname=="Annemarie"& name=="Kempf Schluchter" & votes==53733
replace municipalityno=567 if  firstname=="Annemarie"& name=="Kempf Schluchter" & votes==30055
replace municipalityno=567 if  firstname=="Hans-Ulrich"& name=="Trachsel" & votes==9491
replace municipalityno=567 if  firstname=="Markus"& name=="Grossen" & votes==7221
replace municipalityno=668 if firstname=="Markus"& name=="Horst" & votes==5823
replace municipalityno=668 if firstname=="Dori"& name=="Schaer-Born" & votes==11955
replace municipalityno=551 if  firstname=="Hugo"& name=="Gerber" & votes==2015
replace municipalityno=551 if  firstname=="Elmar"& name=="Jergen" & votes==569
replace municipalityno=957 if  firstname=="Christian"& name=="Waber" & votes==4790
replace municipalityno=957 if  firstname=="Christian"& name=="Waber" & votes==11145
replace municipalityno=995 if  firstname=="Regine"& name=="Bohny" & votes==68577
replace municipalityno=947 if  firstname=="Heidi"& name=="Iseli-Klossner" & votes==26617
replace municipalityno=947 if  firstname=="Ulrich"& name=="Krummenacher" & votes==4969
replace municipalityno=2846 if  firstname=="Edith"& name=="Stauber" & votes==932
replace municipalityno=2326 if  firstname=="Maryline"& name=="Vial" & votes==2117
replace municipalityno=2281 if  firstname=="Philippe"& name=="Chautems" & votes==3597
replace municipalityno=2281 if  firstname=="Philippe"& name=="Chautems" & votes==3008
replace municipalityno=6643 if  firstname=="Jean-Claude"& name=="Vaudroz" & votes==10422
replace municipalityno=6621 if  firstname=="Werner"& name=="Marti" & votes==264
replace municipalityno=6622 if  firstname=="Jean-Jacques"& name=="Griessen" & votes==874
replace municipalityno=6623 if  firstname=="Denise"& name=="Kessler-Nicolet" & votes==16156
replace municipalityno=6623 if  firstname=="Hélène"& name=="Braun-Roth" & votes==16960
replace municipalityno=6621 if  firstname=="Albert"& name=="Anor" & votes==798
replace municipalityno=6621 if  firstname=="Georges"& name=="Meylan" & votes==633
replace municipalityno=6621 if  firstname=="Antonio"& name=="Hodgers" & votes==14187
replace municipalityno=1604 if  firstname=="Fritz"& name=="Hösli" & votes==4089
replace municipalityno=3955 if  firstname=="Andreas"& name=="Thöny" & votes==989
replace municipalityno=3955 if  firstname=="Thomas"& name=="Bigliel" & votes==223
replace municipalityno=3955 if  firstname=="Andreas"& name=="Thöny" & votes==5030
replace municipalityno=3955 if  firstname=="Ernst"& name=="Nigg" & votes==6449
replace municipalityno=3955 if  firstname=="Agnes"& name=="Brandenburger-Caderas" & votes==1897
replace municipalityno=3962 if  firstname=="Jann-Andrea"& name=="Thöny" & votes==528
replace municipalityno=3669 if  firstname=="Simeon"& name=="Bühler" & votes==5726
replace municipalityno=1005 if  firstname=="Manfred"& name=="Aregger" & votes==35859
replace municipalityno=1066 if  firstname=="Markus"& name=="Belser" & votes==608
replace municipalityno=1009 if  firstname=="Manuel"& name=="Schneider" & votes==1557
replace municipalityno=1148 if firstname=="Hilmar"& name=="Gernet" & votes==34666
replace municipalityno=6405 if  firstname=="Laurent"& name=="Debrot" & votes==4468
replace municipalityno=3442 if  firstname=="Emil"& name=="Keller" & votes==10217
replace municipalityno=3442 if  firstname=="Otto"& name=="Köppel" & votes==4875
replace municipalityno=3442 if  firstname=="Susanne"& name=="Vincenz-Stauffacher" & votes==12640
replace municipalityno=3231 if  firstname=="Manfred"& name=="Messmer" & votes==3422
replace municipalityno=3392 if  firstname=="Daniela"& name=="Dürrenmatt" & votes==1015
replace municipalityno=3421 if  firstname=="Hans"& name=="Ruckstuhl" & votes==47001
replace municipalityno=3377 if  firstname=="Peter"& name=="Ledergerber" & votes==9061
replace municipalityno=2458 if firstname=="Anita"& name=="Hug" & votes==15726
replace municipalityno=1322 if  firstname=="Alfred"& name=="Böni" & votes==14206
replace municipalityno=1322 if  firstname=="Dominik"& name=="Zehnder" & votes==4976
replace municipalityno=1322 if  firstname=="Ilias"& name=="Läber" & votes==1414
replace municipalityno=1322 if  firstname=="Albert"& name=="Meile" & votes==1087
replace municipalityno=1322 if  firstname=="Hedy"& name=="Jager-Stählin" & votes==3655
replace municipalityno=1372 if  firstname=="Karl"& name=="Weber" & votes==10726
replace municipalityno=1372 if  firstname=="Josef-Mariä"& name=="Steiner" & votes==1443
replace municipalityno=1372 if  firstname=="Xaver"& name=="Schuler" & votes==16304
replace municipalityno=1372 if  firstname=="Karl"& name=="Weber" & votes==10527
replace municipalityno=4891 if  firstname=="Josef"& name=="Kressibucher" & votes==475
replace municipalityno=4841 if  firstname=="Willy"& name=="Schmidhauser" & votes==3871
replace municipalityno=4836 if  firstname=="Willy"& name=="Schmidhauser" & votes==5505
replace municipalityno=4836 if  firstname=="Willy J."& name=="Schmidhauser" & votes==5572
replace municipalityno=4836 if  firstname=="Ruedi"& name=="Buzek" & votes==6087
replace municipalityno=4508 if firstname=="Werner"& name=="Messmer" & votes==9118
replace municipalityno=4687 if  firstname=="Margrit"& name=="Beck-Föhn" & votes==4162
replace municipalityno=4533 if  firstname=="Andreas"& name=="Aeberhardt" & votes==4359
replace municipalityno=5625 if  firstname=="Nicolas"& name=="Daïna" & votes==31167
replace municipalityno=5886 if  firstname=="Robert"& name=="Gurtner" & votes==959
replace municipalityno=5905 if  firstname=="Alice"& name=="Glauser" & votes==37643
replace municipalityno=5582 if  firstname=="Michèle"& name=="Gay Vallotton" & votes==28976
replace municipalityno=5715 if  firstname=="Anna"& name=="Blanchoud" & votes==1518
replace municipalityno=5721 if  firstname=="Patrick"& name=="Vallat" & votes==4700
replace municipalityno=5721 if  firstname=="David"& name=="Mayer" & votes==1702
replace municipalityno=5474 if  firstname=="Estella"& name=="Liniger" & votes==11577
replace municipalityno=5886 if  firstname=="Catherine"& name=="Buchet Buillard" & votes==29727
replace municipalityno=5611 if  firstname=="Philippe"& name=="Bornand" & votes==2944
replace municipalityno=6136 if  firstname=="Pascal"& name=="Couchepin" & votes==28318
replace municipalityno=6136 if  firstname=="Vital"& name=="Darbellay" & votes==36730
replace municipalityno=6137 if  firstname=="Christophe"& name=="Darbellay" & votes==46197
replace municipalityno=6137 if  firstname=="Christophe"& name=="Mariéthoz" & votes==615
replace municipalityno=6137 if  firstname=="Christophe"& name=="Darbellay" & votes==40241
replace municipalityno=6152 if  firstname=="Julie"& name=="Delaloye" & votes==858
replace municipalityno=96 if  firstname=="Heinz"& name=="Ryser" & votes==11859
replace municipalityno=96 if  firstname=="Hans Rudolf"& name=="Metz" & votes==14163
replace municipalityno=96 if  firstname=="Tanja"& name=="Potzmader" & votes==711
replace municipalityno=96 if  firstname=="Matthias"& name=="Rubi" & votes==1276
replace municipalityno=96 if  firstname=="Matthias"& name=="Rubi" & votes==902
replace municipalityno=21 if  firstname=="Werner"& name=="Gmür" & votes==573
replace municipalityno=131 if  firstname=="Hanspeter"& name=="Clesle" & votes==8506
replace municipalityno=131 if  firstname=="Davide"& name=="Loss" & votes==62814
replace municipalityno=131 if  firstname=="Lukas"& name=="Gloor" & votes==1024
replace municipalityno=131 if  firstname=="Clemens"& name=="Ruckstuhl" & votes==15912
replace municipalityno=115 if  firstname=="Marcel"& name=="Lenggenhager" & votes==20228
replace municipalityno=242 if  firstname=="Reinhard"& name=="Walther" & votes==701
replace municipalityno=192 if  firstname=="Konrad"& name=="Basler" & votes==44516
replace municipalityno=192 if  firstname=="Konrad"& name=="Basler" & votes==57545
replace municipalityno=192 if  firstname=="Hans Jörg"& name=="Fischer" & votes==5700
replace municipalityno=192 if  firstname=="Hans Jörg"& name=="Fischer" & votes==9511
replace municipalityno=192 if  firstname=="Hans Jörg"& name=="Fischer" & votes==3731
replace municipalityno=231 if  firstname=="Paul"& name=="Wüthrich" & votes==1617
replace municipalityno=231 if  firstname=="Kurt"& name=="Rutz" & votes==23766
replace municipalityno=231 if  firstname=="Marcus"& name=="Nägeli" & votes==5465
replace municipalityno=120 if  firstname=="Alfred"& name=="Schaffner" & votes==203
replace municipalityno=92 if  firstname=="Melanie (Melä)"& name=="Amstad" & votes==213
replace municipalityno=231 if  firstname=="Hans"& name=="Baumann" & votes==566
replace municipalityno=118 if firstname=="Ursula"& name=="Hänni-Hauser" & votes==2480
replace municipalityno=118 if firstname=="Peter"& name=="Honegger" & votes==10011
replace municipalityno=118 if firstname=="Alexandre-Pierre"& name=="Frick" & votes==654
replace municipalityno=118 if firstname=="Rolf"& name=="Strasser" & votes==698
replace municipalityno=118 if firstname=="Anton"& name=="Melliger" & votes==13647
replace municipalityno=118 if firstname=="Daniel"& name=="Hänni" & votes==630
replace municipalityno=118 if firstname=="Stephan"& name=="Berndt" & votes==609
replace municipalityno=42 if  firstname=="Carmen"& name=="Schwager" & votes==1227
replace municipalityno=200 if  firstname=="Rolf Erich"& name=="Peter" & votes==2565
replace municipalityno=200 if  firstname=="Andrea"& name=="Barbadimos" & votes==1179
replace municipalityno=112 if  firstname=="Erich"& name=="Vontobel" & votes==7990
replace municipalityno=112 if  firstname=="Diego"& name=="Gehrig" & votes==1752
replace municipalityno=112 if  firstname=="Heinz"& name=="Lüscher" & votes==8224
replace municipalityno=112 if  firstname=="Monika"& name=="Stahl" & votes==2014
replace municipalityno=261 if  firstname=="Andreas"& name=="Hugi" & votes==7397
replace municipalityno=261 if  firstname=="Peter"& name=="Püntener" & votes==6643
