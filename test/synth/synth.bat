set EXE="D:/Projekty/clusty-dev/src/x64/Release/clusty.exe"

%EXE% --id-cols id1 id2 --distance-col ani --similarity --min ani 0.70 --numeric-ids synth.ani numeric.clusty
%EXE% --id-cols id1 id2 --distance-col ani --similarity --min ani 0.70 --numeric-ids synth.ani numeric.objs.clusty --objects-file synth.ids 

%EXE% --id-cols id1 id2 --distance-col ani --similarity --min ani 0.70 --numeric-ids synth.ani numeric.reps.clusty --out-representatives 
%EXE% --id-cols id1 id2 --distance-col ani --similarity --min ani 0.70 --numeric-ids synth.ani numeric.objs.reps.clusty --objects-file synth.ids --out-representatives 

%EXE% --id-cols name1 name2 --distance-col ani --similarity --min ani 0.70 synth.ani named.clusty
%EXE% --id-cols name1 name2 --distance-col ani --similarity --min ani 0.70 synth.ani named.objs.clusty --objects-file synth.ids 

%EXE% --id-cols name1 name2 --distance-col ani --similarity --min ani 0.70 synth.ani named.reps.clusty --out-representatives 
%EXE% --id-cols name1 name2 --distance-col ani --similarity --min ani 0.70 synth.ani named.objs.reps.clusty --objects-file synth.ids --out-representatives 