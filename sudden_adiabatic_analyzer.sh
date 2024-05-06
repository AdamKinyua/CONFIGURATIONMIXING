#! /bin/bash
echo 'Mo of ethylene at time t as the point charge moves with time - reference mo'
grep -A26 'Initial HOMO-7 at' output.txt | grep ' 1 ' > ihomo-7;
grep -A26 'Initial HOMO-6 at' output.txt | grep ' 2 ' > ihomo-6;
grep -A26 'Initial HOMO-5 at' output.txt | grep ' 3 ' > ihomo-5;
grep -A26 'Initial HOMO-4 at' output.txt | grep ' 4 ' > ihomo-4;
grep -A26 'Initial HOMO-3 at' output.txt | grep ' 5 ' > ihomo-3;
grep -A26 'Initial HOMO-2 at' output.txt | grep ' 6 ' > ihomo-2;
grep -A26 'Initial HOMO-1 at' output.txt | grep ' 7 ' > ihomo-1;
grep -A26 'Initial HOMO at' output.txt | grep ' 8 ' > ihomo;
grep -A26 'Initial LUMO at' output.txt | grep ' 9 ' > ilumo;
grep -A26 'Initial LUMO+1 at' output.txt | grep ' 10 ' > ilumo+1;
grep -A26 'Initial LUMO+2 at' output.txt | grep ' 11 ' > ilumo+2;

echo 'Configuration mixing between mo-i and mo-i at time t - tracks how each mo changes with time';
grep -A26 'HOMO-7 vector at' output.txt | grep ' 1 ' > homo-7-homo-7;
grep -A26 'HOMO-6 vector at' output.txt | grep ' 2 ' > homo-6-homo-6;
grep -A26 'HOMO-5 vector at' output.txt | grep ' 3 ' > homo-5-homo-5;
grep -A26 'HOMO-4 vector at' output.txt | grep ' 4 ' > homo-4-homo-4;
grep -A26 'HOMO-3 vector at' output.txt | grep ' 5 ' > homo-3-homo-3;
grep -A26 'HOMO-2 vector at' output.txt | grep ' 6 ' > homo-2-homo-2;
grep -A26 'HOMO-1 vector at' output.txt | grep ' 7 ' > homo-1-homo-1;
grep -A26 'HOMO vector at' output.txt | grep ' 8 ' > homo0-homo;
grep -A26 'LUMO vector at' output.txt | grep ' 9 ' > lumo-lumo;
grep -A26 'LUMO+1 vector at' output.txt | grep ' 10 ' > lumo+1-lumo+1;
grep -A26 'LUMO+2 vector at' output.txt | grep ' 11 ' > lumo+2-lumo+2;

echo 'Configuration mixing between homo and mo-i at time t';
grep -A26 'HOMO-7 vector at' output.txt | grep ' 8 ' > homohomo7;
grep -A26 'HOMO-6 vector at' output.txt | grep ' 8 ' > homohomo6;
grep -A26 'HOMO-5 vector at' output.txt | grep ' 8 ' > homohomo5;
grep -A26 'HOMO-4 vector at' output.txt | grep ' 8 ' > homohomo4;
grep -A26 'HOMO-3 vector at' output.txt | grep ' 8 ' > homohomo3;
grep -A26 'HOMO-2 vector at' output.txt | grep ' 8 ' > homohomo2;
grep -A26 'HOMO-1 vector at' output.txt | grep ' 8 ' > homohomo1;
grep -A26 'HOMO vector at' output.txt | grep ' 8 ' > homohomo;
grep -A26 'LUMO vector at' output.txt | grep ' 8 ' > homolumo;
grep -A26 'LUMO+1 vector at' output.txt | grep ' 8 ' > homolumo1;
grep -A26 'LUMO+2 vector at' output.txt | grep ' 8 ' > homolumo2;

echo 'done!'
