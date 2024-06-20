cd libs

rmdir /S /Q winflexbison
mkdir winflexbison
cd winflexbison
curl -L --ssl-no-revoke https://github.com/refresh-bio-dependencies/winflexbison/releases/download/v2.5.25/win_flex_bison-2.5.25.zip --output win_flex_bison-2.5.25.zip
tar -xf win_flex_bison-2.5.25.zip
set PATH=%cd%;%PATH%

cd ../igraph
rmdir /S /Q build
mkdir build 
cd build

cmake ..
cmake --build . --config Debug
cmake --build . --config Release
cd ../../..