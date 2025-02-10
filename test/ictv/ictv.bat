set EXE="../../src/x64/Release/clusty.exe"
set PY=C:\Python\Python312\python.exe


%EXE% --algo single --objects-file ictv.list --similarity --min ani 0.95  ictv.ani ictv.single.95.csv 
%EXE% --algo single --objects-file ictv.list --similarity --min ani 0.95  ictv.ani ictv.single.95.reps.csv --out-representatives
%EXE% --algo single --objects-file ictv.list --similarity --min ani 0.70  ictv.ani ictv.single.70.csv
%EXE% --algo single --objects-file ictv.list --similarity --min ani 0.70  ictv.ani ictv.single.70.reps.csv --out-representatives
%PY% ../cmp.py ictv.single.95.csv ictv.single.95.python.csv
%PY% ../cmp.py ictv.single.95.reps.csv ictv.single.95.python.csv
%PY% ../cmp.py ictv.single.70.csv ictv.single.70.python.csv
%PY% ../cmp.py ictv.single.70.reps.csv ictv.single.70.python.csv


%EXE% --algo complete --objects-file ictv.list --similarity --min ani 0.95  ictv.ani ictv.complete.95.csv 
%EXE% --algo complete --objects-file ictv.list --similarity --min ani 0.95  ictv.ani ictv.complete.95.reps.csv --out-representatives
%EXE% --algo complete --objects-file ictv.list --similarity --min ani 0.70  ictv.ani ictv.complete.70.csv
%EXE% --algo complete --objects-file ictv.list --similarity --min ani 0.70  ictv.ani ictv.complete.70.reps.csv --out-representatives
%PY% ../cmp.py ictv.complete.95.csv ictv.complete.95.python.csv
%PY% ../cmp.py ictv.complete.95.reps.csv ictv.complete.95.python.csv
%PY% ../cmp.py ictv.complete.70.csv ictv.complete.70.python.csv
%PY% ../cmp.py ictv.complete.70.reps.csv ictv.complete.70.python.csv


%EXE% --algo complete --objects-file ictv.list --similarity --min ani 0.95 --numeric-ids ictv.num ictv.num.complete.95.csv 
%EXE% --algo complete --objects-file ictv.list --similarity --min ani 0.95 --numeric-ids ictv.num ictv.num.complete.95.reps.csv --out-representatives
%EXE% --algo complete --objects-file ictv.list --similarity --min ani 0.70 --numeric-ids ictv.num ictv.num.complete.70.csv
%EXE% --algo complete --objects-file ictv.list --similarity --min ani 0.70 --numeric-ids  ictv.num ictv.num.complete.70.reps.csv --out-representatives
%PY% ../cmp.py ictv.num.complete.95.csv ictv.complete.95.python.csv
%PY% ../cmp.py ictv.num.complete.95.reps.csv ictv.complete.95.python.csv
%PY% ../cmp.py ictv.num.complete.70.csv ictv.complete.70.python.csv
%PY% ../cmp.py ictv.num.complete.70.reps.csv ictv.complete.70.python.csv


%EXE% --algo set-cover --objects-file ictv.list --similarity --min ani 0.70 ictv.ani ictv.set-cover.70.csv
%EXE% --algo cd-hit --objects-file ictv.list --similarity --min ani 0.70 ictv.ani ictv.cd-hit.70.csv
%EXE% --algo uclust --objects-file ictv.list --similarity --min ani 0.70 ictv.ani ictv.uclust.70.csv
%EXE% --algo leiden --objects-file ictv.list --similarity --min ani 0.70 ictv.ani ictv.leiden.70.csv

%EXE% --algo set-cover --objects-file ictv.list --similarity --min ani 0.70 ictv.ani ictv.set-cover.70.reps.csv --out-representatives
%EXE% --algo cd-hit --objects-file ictv.list --similarity --min ani 0.70 ictv.ani ictv.cd-hit.70.reps.csv --out-representatives
%EXE% --algo uclust --objects-file ictv.list --similarity --min ani 0.70 ictv.ani ictv.uclust.70.reps.csv --out-representatives
%EXE% --algo leiden --objects-file ictv.list --similarity --min ani 0.70 ictv.ani ictv.leiden.70.reps.csv --out-representatives

