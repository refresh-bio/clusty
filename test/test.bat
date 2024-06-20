set EXE="../src/x64/Release/rapid-cluster.exe"
set PY=C:\Python38\python.exe



#%EXE% --algo slink --objects-file vir61.list --similarity --min ani 0.95 --id-cols id2 id1 --distance-col ani vir61.ani vir61.slink.95 --out-representatives
#%EXE% --algo slink --objects-file vir61.list --similarity --min ani 0.70 --id-cols id2 id1 --distance-col ani vir61.ani vir61.slink.70 --out-csv
#%PY% cmp.py vir61.slink.95 vir61.single.95.csv
#%PY% cmp.py vir61.slink.70 vir61.single.70.csv


#%EXE% --algo single --objects-file vir61.list --similarity --min ani 0.95 --id-cols id2 id1 --distance-col ani vir61.ani vir61.single.95 --out-representatives
#%EXE% --algo single --objects-file vir61.list --similarity --min ani 0.70 --id-cols id2 id1 --distance-col ani vir61.ani vir61.single.70
#%PY% cmp.py vir61.single.95 vir61.single.95.csv
#%PY% cmp.py vir61.single.70 vir61.single.70.csv

#%EXE% --algo complete --objects-file vir61.list --similarity --min ani 0.95 --id-cols id2 id1 --distance-col ani vir61.ani vir61.complete.95
#%EXE% --algo complete --objects-file vir61.list --similarity --min ani 0.70 --id-cols id2 id1 --distance-col ani vir61.ani vir61.complete.70
#%PY% cmp.py vir61.complete.95 vir61.complete.95.csv
#%PY% cmp.py vir61.complete.70 vir61.complete.70.csv

#%EXE% --algo slink --objects-file ictv.list --similarity --min ani 0.95  ictv.ani ictv.slink.95.csv 
#%EXE% --algo slink --objects-file ictv.list --similarity --min ani 0.95  ictv.ani ictv.slink.95.reps.csv --out-representatives
#%EXE% --algo slink --objects-file ictv.list --similarity --min ani 0.70  ictv.ani ictv.slink.70.csv
#%EXE% --algo slink --objects-file ictv.list --similarity --min ani 0.70  ictv.ani ictv.slink.70.reps.csv --out-representatives
#%PY% cmp.py ictv.slink.95.csv ictv.single.95.python.csv
#%PY% cmp.py ictv.slink.70.reps.csv ictv.single.70.python.csv

%EXE% --algo single --objects-file ictv.list --similarity --min ani 0.95  ictv.ani ictv.single.95.csv 
%EXE% --algo single --objects-file ictv.list --similarity --min ani 0.95  ictv.ani ictv.single.95.reps.csv --out-representatives
%EXE% --algo single --objects-file ictv.list --similarity --min ani 0.70  ictv.ani ictv.single.70.csv
%EXE% --algo single --objects-file ictv.list --similarity --min ani 0.70  ictv.ani ictv.single.70.reps.csv --out-representatives
%PY% cmp.py ictv.single.95.csv ictv.single.95.python.csv
%PY% cmp.py ictv.single.70.reps.csv ictv.single.70.python.csv


%EXE% --algo complete --objects-file ictv.list --similarity --min ani 0.95  ictv.ani ictv.complete.95.csv 
%EXE% --algo complete --objects-file ictv.list --similarity --min ani 0.95  ictv.ani ictv.complete.95.reps.csv --out-representatives
%EXE% --algo complete --objects-file ictv.list --similarity --min ani 0.70  ictv.ani ictv.complete.70.csv
%EXE% --algo complete --objects-file ictv.list --similarity --min ani 0.70  ictv.ani ictv.complete.70.reps.csv --out-representatives
%PY% cmp.py ictv.complete.95.csv ictv.complete.95.python.csv
%PY% cmp.py ictv.complete.70.reps.csv ictv.complete.70.python.csv


%EXE% --algo set-cover --objects-file ictv.list --similarity --min ani 0.70 ictv.ani ictv.set-cover.70.csv
%EXE% --algo cd-hit --objects-file ictv.list --similarity --min ani 0.70 ictv.ani ictv.cd-hit.70.csv
%EXE% --algo uclust --objects-file ictv.list --similarity --min ani 0.70 ictv.ani ictv.uclust.70.csv
%EXE% --algo leiden --objects-file ictv.list --similarity --min ani 0.70 ictv.ani ictv.leiden.70.csv

%EXE% --algo set-cover --objects-file ictv.list --similarity --min ani 0.70 ictv.ani ictv.set-cover.70.reps.csv --out-representatives
%EXE% --algo cd-hit --objects-file ictv.list --similarity --min ani 0.70 ictv.ani ictv.cd-hit.70.reps.csv --out-representatives
%EXE% --algo uclust --objects-file ictv.list --similarity --min ani 0.70 ictv.ani ictv.uclust.70.reps.csv --out-representatives
%EXE% --algo leiden --objects-file ictv.list --similarity --min ani 0.70 ictv.ani ictv.leiden.70.reps.csv --out-representatives

