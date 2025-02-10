set EXE="../../src/x64/Release/clusty.exe"
#set EXE="../../../clusty/src/x64/Release/clusty.exe"
set PY=C:\Python\Python312\python.exe

%EXE% --algo single --objects-file vir61.list --similarity --min ani 0.95 --id-cols id2 id1 --distance-col ani vir61.ani vir61.single.95 --out-csv
%EXE% --algo single --objects-file vir61.list --similarity --min ani 0.95 --id-cols id2 id1 --distance-col ani vir61.ani vir61.single.95.reps --out-representatives
%EXE% --algo single --objects-file vir61.list --similarity --min ani 0.70 --id-cols id2 id1 --distance-col ani vir61.ani vir61.single.70
%EXE% --algo single --objects-file vir61.list --similarity --min ani 0.70 --id-cols id2 id1 --distance-col ani vir61.ani vir61.single.70.reps --out-representatives
%PY% ../cmp.py vir61.single.95.python.csv vir61.single.95
%PY% ../cmp.py vir61.single.95.python.csv vir61.single.95.reps
%PY% ../cmp.py vir61.single.70.python.csv vir61.single.70
%PY% ../cmp.py vir61.single.70.python.csv vir61.single.70.reps

%EXE% --algo complete --objects-file vir61.list --similarity --min ani 0.95 --id-cols id2 id1 --distance-col ani vir61.ani vir61.complete.95
%EXE% --algo complete --objects-file vir61.list --similarity --min ani 0.95 --id-cols id2 id1 --distance-col ani vir61.ani vir61.complete.95.reps --out-representatives
%EXE% --algo complete --objects-file vir61.list --similarity --min ani 0.70 --id-cols id2 id1 --distance-col ani vir61.ani vir61.complete.70
%EXE% --algo complete --objects-file vir61.list --similarity --min ani 0.70 --id-cols id2 id1 --distance-col ani vir61.ani vir61.complete.70.reps --out-representatives
%PY% ../cmp.py vir61.complete.95.python.csv vir61.complete.95
%PY% ../cmp.py vir61.complete.95.python.csv vir61.complete.95.reps
%PY% ../cmp.py vir61.complete.70.python.csv vir61.complete.70
%PY% ../cmp.py vir61.complete.70.python.csv vir61.complete.70.reps
